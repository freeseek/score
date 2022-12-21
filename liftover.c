/* The MIT License

   Copyright (C) 2022 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/kseq.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "bcftools.h"
#include "regidx.h"
#include "htslib/khash.h" // required to reset the contigs dictionary and table
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)

#define LIFTOVER_VERSION "2022-12-21"

#define FLIP_TAG "FLIP"
#define SWAP_TAG "SWAP"
#define DROP_TAGS "INFO/AC,FMT/AC"
#define REVERSE_TAGS "INFO/AF:1,FMT/AF:1,FMT/DS:2,FMT/AP1:1,FMT/AP2:1"
#define FLIP_TAGS "FMT/EZ,FMT/ES,FMT/ED"
#define GENOTYPE_TAGS "INFO/ALLELE_A,INFO/ALLELE_B"

static const char revnt[128] = {
    '\x00', '\x01', '\x02', '\x03', '\x04', '\x05', '\x06', '\a',   '\b',   '\t',   '\n',   '\v',   '\f',
    '\r',   '\x0e', '\x0f', '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17', '\x18', '\x19',
    '\x1a', '\e',   '\x1c', '\x1d', '\x1e', '\x1f', ' ',    '!',    '\"',   '#',    '$',    '%',    '&',
    '\'',   '(',    ')',    '*',    '+',    ',',    '-',    '.',    '/',    '0',    '1',    '2',    '3',
    '4',    '5',    '6',    '7',    '8',    '9',    ':',    ';',    '<',    '=',    '>',    '\?',   '@',
    'T',    'B',    'G',    'D',    'E',    'F',    'C',    'H',    'I',    'J',    'K',    'L',    'M',
    'N',    'O',    'P',    'Q',    'R',    'S',    'A',    'U',    'V',    'W',    'X',    'Y',    'Z',
    '[',    '\\',   ']',    '^',    '_',    '`',    't',    'b',    'g',    'd',    'e',    'f',    'c',
    'h',    'i',    'j',    'k',    'l',    'm',    'n',    'o',    'p',    'q',    'r',    's',    'a',
    'u',    'v',    'w',    'x',    'y',    'z',    '{',    '|',    '}',    '~',    '\x7f'};

typedef struct {
    int tStart;
    int qStart;
    int size;
    int chain_ind;
    int tStart_gap; // size of gap to reach previous contiguous block
    int tEnd_gap;   // size of gap to reach next contiguous block
} block_t;

typedef struct {
    uint64_t score;
    int t_rid;
    int tSize;
    int tStart;
    int tEnd;
    int q_rid;
    int qSize;
    int qStrand;
    int qStart;
    int qEnd;
    int id;
    int block_ind;
    int n_blocks;
} chain_t;

typedef struct {
    int gt_id;
    const char *flip_tag;
    const char *swap_tag;

    int n_drop_tags;
    int *drop_tags_int_id;
    int *drop_tags_coltype;

    int n_reverse_tags;
    int *reverse_tags_int_id;
    int *reverse_tags_coltype;
    float *reverse_tags_sum;

    int n_flip_tags;
    int *flip_tags_int_id;
    int *flip_tags_coltype;

    int n_genotype_tags;
    int *genotype_tags_int_id;
    int *genotype_tags_coltype;

    int n_GR_tags;
    int *GR_tags_int_id;
    int *GR_tags_coltype;
} tags_t;

typedef struct {
    bcf_hdr_t *in_hdr;
    bcf_hdr_t *out_hdr;
    faidx_t *src_fai;
    faidx_t *dst_fai;
    int n_chains;
    chain_t *chains;
    block_t *blocks;
    regidx_t *idx;
    regitr_t *itr;
    htsFile *reject_fh;

    int indel_win;
    int aln_win;
    int in_mt_rid;
    int out_mt_rid;
    int lift_mt;
    int no_left_align;
    int write_src;
    tags_t tags;

    int warning_symbolic;
    int warning_indel;
    int ntotal;
    int nswapped;
    int nref_added;
    int nrejected;

    kstring_t tmp_kstr;
    kstring_t *tmp_als;
    int m_tmp_als;
    void *tmp_arr;
    int m_tmp_arr;
} args_t;

args_t *args;

static void tags_destroy(tags_t *tags) {
    free(tags->drop_tags_int_id);
    free(tags->drop_tags_coltype);
    free(tags->reverse_tags_int_id);
    free(tags->reverse_tags_coltype);
    free(tags->reverse_tags_sum);
    free(tags->flip_tags_int_id);
    free(tags->flip_tags_coltype);
    free(tags->genotype_tags_int_id);
    free(tags->genotype_tags_coltype);
    free(tags->GR_tags_int_id);
    free(tags->GR_tags_coltype);
}

/****************************************
 * CHAIN FUNCTIONS                      *
 ****************************************/

static inline int bcf_hdr_name2id_flexible(const bcf_hdr_t *hdr, char *chr) {
    if (!chr) return -1;
    char buf[] = {'c', 'h', 'r', '\0', '\0', '\0'};
    int rid = bcf_hdr_name2id(hdr, chr);
    if (rid >= 0) return rid;
    if (strncmp(chr, "chr", 3) == 0) rid = bcf_hdr_name2id(hdr, chr + 3);
    if (rid >= 0) return rid;
    strncpy(buf + 3, chr, 2);
    rid = bcf_hdr_name2id(hdr, buf);
    if (rid >= 0) return rid;
    if (strcmp(chr, "23") == 0 || strcmp(chr, "25") == 0 || strcmp(chr, "XY") == 0 || strcmp(chr, "XX") == 0
        || strcmp(chr, "PAR1") == 0 || strcmp(chr, "PAR2") == 0) {
        rid = bcf_hdr_name2id(hdr, "X");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrX");
    } else if (strcmp(chr, "24") == 0) {
        rid = bcf_hdr_name2id(hdr, "Y");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrY");
    } else if (strcmp(chr, "26") == 0 || strcmp(chr, "MT") == 0 || strcmp(chr, "chrM") == 0) {
        rid = bcf_hdr_name2id(hdr, "MT");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrM");
    }
    return rid;
}

// load the chain file (see http://genome.ucsc.edu/goldenPath/help/chain.html)
static int read_chains(htsFile *fp, const bcf_hdr_t *in_hdr, const bcf_hdr_t *out_hdr, int max_snp_gap,
                       int max_indel_gap, chain_t **chains, block_t **blocks) {
    int n_chains = 0;
    int n_blocks = 0;
    int m_chains = 0;
    int m_blocks = 0;
    char *tmp = NULL;
    kstring_t str = {0, 0, NULL};
    int moff = 0, *off = NULL;
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        hts_expand(chain_t, n_chains + 1, m_chains, *chains);
        chain_t *chain = &(*chains)[n_chains++];
        int ncols = ksplit_core(str.s, ' ', &moff, &off);
        if (ncols != 13) error("Wrong number of columns in the chain file: %s\n", fp->fn);

        // read Header Lines
        assert(strcmp(&str.s[off[0]], "chain") == 0);
        chain->score = (uint64_t)strtoll(&str.s[off[1]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[1]], fp->fn);
        chain->t_rid = bcf_hdr_name2id_flexible(in_hdr, &str.s[off[2]]);
        chain->tSize = strtol(&str.s[off[3]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[3]], fp->fn);
        if (chain->t_rid >= 0) {
            uint64_t len = in_hdr->id[BCF_DT_CTG][chain->t_rid].val->info[0];
            if (len != 0 && chain->tSize != len)
                fprintf(stderr,
                        "Warning: source contig %s has length %" PRId64 " in the VCF and length %s in the chain file\n",
                        &str.s[off[2]], len, &str.s[off[3]]);
        }
        assert(str.s[off[4]] == '+');
        chain->tStart = strtol(&str.s[off[5]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[5]], fp->fn);
        chain->tEnd = strtol(&str.s[off[6]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[6]], fp->fn);
        chain->q_rid = bcf_hdr_name2id_flexible(out_hdr, &str.s[off[7]]);
        chain->qSize = strtol(&str.s[off[8]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[8]], fp->fn);
        if (chain->q_rid >= 0) {
            uint64_t len = out_hdr->id[BCF_DT_CTG][chain->q_rid].val->info[0];
            if (len != 0 && chain->qSize != len)
                fprintf(stderr,
                        "Warning: query contig %s has length %" PRId64 " in the VCF and length %s in the chain file\n",
                        &str.s[off[7]], len, &str.s[off[8]]);
        }
        assert(str.s[off[9]] == '+' || str.s[off[9]] == '-');
        chain->qStrand = str.s[off[9]] == '-';
        chain->qStart = strtol(&str.s[off[10]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[10]], fp->fn);
        chain->qEnd = strtol(&str.s[off[11]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[11]], fp->fn);
        chain->id = strtol(&str.s[off[12]], &tmp, 0);
        if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[12]], fp->fn);

        // read Data Lines
        int tStart = 0;
        int qStart = 0;
        int dt = -1;
        int dq = -1;
        chain->block_ind = n_blocks;
        chain->n_blocks = 0;
        block_t *block;
        int merge_in_progress = 0;
        while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
            int ncols = ksplit_core(str.s, '\t', &moff, &off);
            if (ncols != 1 && ncols != 3) error("Wrong number of columns in the chain file: %s\n", fp->fn);
            int size = strtol(&str.s[off[0]], &tmp, 0);
            if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[0]], fp->fn);

            if (merge_in_progress) {
                block->size += size;
            } else {
                hts_expand(block_t, n_blocks + 1, m_blocks, *blocks);
                block = &(*blocks)[n_blocks++];
                block->tStart = tStart;
                block->qStart = qStart;
                block->size = size;
                block->chain_ind = n_chains - 1;
                block->tStart_gap = dt <= max_indel_gap ? dt : -1;
                block->tEnd_gap = -1;
                chain->n_blocks++;
            }

            if (ncols == 1) {
                tStart += size;
                assert(chain->tStart + tStart == chain->tEnd);
                qStart += size;
                assert(chain->qStart + qStart == chain->qEnd);
                block->tEnd_gap = -1;
            } else {
                dt = strtol(&str.s[off[1]], &tmp, 0);
                if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[1]], fp->fn);
                dq = strtol(&str.s[off[2]], &tmp, 0);
                if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[2]], fp->fn);
                tStart += size + dt;
                qStart += size + dq;
                merge_in_progress = dt == dq && dt <= max_snp_gap;
                if (merge_in_progress)
                    block->size += dt;
                else
                    block->tEnd_gap = dt <= max_indel_gap ? dt : -1;
            }
        }
    }

    free(off);
    free(str.s);
    return n_chains;
}

static void write_chains(FILE *stream, const bcf_hdr_t *in_hdr, const bcf_hdr_t *out_hdr, const chain_t *chains,
                         int n_chains, block_t *blocks) {
    for (int i = 0; i < n_chains; i++) {
        const chain_t *chain = &chains[i];
        const char *tName = bcf_hdr_id2name(in_hdr, chain->t_rid);
        const char *qName = bcf_hdr_id2name(out_hdr, chain->q_rid);
        fprintf(stream, "chain %" PRIhts_pos " %s %d + %d %d %s %d %c %d %d %d\n", chain->score, tName, chain->tSize,
                chain->tStart, chain->tEnd, qName, chain->qSize, chain->qStrand ? '-' : '+', chain->qStart, chain->qEnd,
                chain->id);
        for (int j = 0; j < chain->n_blocks; j++) {
            const block_t *block = &blocks[chain->block_ind + j];
            fprintf(stream, "%s\t%d\t%d\t%d\t%" PRIhts_pos "\t%c\t%s\t%d\t%d\n", tName, chain->tStart + block->tStart,
                    chain->tStart + block->tStart + block->size, chain->id, chain->score, chain->qStrand ? '-' : '+',
                    qName,
                    chain->qStrand ? chain->qSize - chain->qStart - block->qStart - block->size
                                   : chain->qStart + block->qStart,
                    chain->qStrand ? chain->qSize - chain->qStart - block->qStart
                                   : chain->qStart + block->qStart + block->size);
        }
    }
}

static regidx_t *regidx_init_chains(const bcf_hdr_t *in_hdr, const chain_t *chains, int n_chains, block_t *blocks) {
    regidx_t *idx = regidx_init(NULL, NULL, NULL, sizeof(uint64_t), NULL);
    for (int i = 0; i < n_chains; i++) {
        const chain_t *chain = &chains[i];
        if (chain->t_rid < 0 || chain->q_rid < 0) continue;
        const char *name = bcf_hdr_id2name(in_hdr, chain->t_rid);
        int len = strlen(name);
        for (int j = 0; j < chain->n_blocks; j++) {
            int block_ind = chain->block_ind + j;
            const block_t *block = &blocks[block_ind];
            regidx_push(idx, (char *)name, (char *)name + len, chain->tStart + block->tStart + 1,
                        chain->tStart + block->tStart + block->size, (void *)&block_ind);
        }
    }
    return idx;
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Lift over a VCF from one genome build to another.\n"; }

const char *usage(void) {
    return "\n"
           "About: Lift over a VCF from one genome build to another. "
           "(version " LIFTOVER_VERSION
           " https://github.com/freeseek/score)\n"
           "Usage: bcftools +liftover [General Options] -- [Plugin Options]\n"
           "Options:\n"
           "   run \"bcftools plugin\" for a list of common options\n"
           "\n"
           "Plugin options:\n"
           "   -s, --src-fasta-ref <file>      source reference sequence in fasta format\n"
           "   -f, --fasta-ref <file>          destination reference sequence in fasta format\n"
           "       --set-cache-size <int>      select fasta cache size in bytes\n"
           "   -c, --chain <file>              UCSC liftOver chain file\n"
           "       --max-snp-gap <int>         maximum distance to merge contiguous blocks separated by same distance "
           "[1]\n"
           "       --max-indel-gap <int>       maximum distance between contiguous blocks to pad alleles [20]\n"
           "       --indel-win <int>           maximum distance between two edges of an indel to accept liftover "
           "[250]\n"
           "       --lift-mt                   force liftover of MT/chrMT [automatically determined from contig "
           "lengths]\n"
           "       --print-blocks <file>       output contiguous blocks used for the liftOver\n"
           "       --no-left-align             do not attempt to left align indels after liftover\n"
           "       --reject <file>             output variants that cannot be lifted over\n"
           "   -O, --reject-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "       --write-src                 write the source contig/position/alleles for lifted variants\n"
           "\n"
           "Options for how to update INFO/FORMAT records:\n"
           "       --flip-tag <string>         INFO annotation flag to record whether alleles are flipped [FLIP]\n"
           "       --swap-tag <string>         INFO annotation to record when alleles are swapped [SWAP]\n"
           "       --tags-to-drop <list>       INFO and FORMAT tags to drop when alleles are swapped [" DROP_TAGS
           "]\n"
           "       --tags-to-reverse <list>    INFO and FORMAT tags to be reversed when alleles are swapped (must be "
           "Number=A,Type=Float)\n"
           "                                   [" REVERSE_TAGS
           "]\n"
           "       --tags-to-flip <list>       INFO and FORMAT tags that have the sign flipped when alleles are "
           "swapped (must be Number=A)\n"
           "                                   [" FLIP_TAGS
           "]\n"
           "       --tags-genotype <list>      INFO and FORMAT tags with genotype integers like FORMAT/GT (must be "
           "Type=Integer)\n"
           "                                   [" GENOTYPE_TAGS
           "]\n"
           "\n"
           "Examples:\n"
           "      bcftools +liftover -Ob -o output.hg38.bcf input.hg19.bcf -- \\\n"
           "        -s human_g1k_v37.fasta -f Homo_sapiens_assembly38.fasta -c hg19ToHg38.over.chain.gz\n"
           "      bcftools +liftover -Oz -o chm13v2.0_dbSNPv155.vcf.gz GRCh38_dbSNPv155.vcf.gz -- \\\n"
           "        -s Homo_sapiens_assembly38.fasta -f chm13v2.0.fa -c hg38-chm13v2.over.chain.gz\n"
           "\n"
           "To obtain UCSC liftOver chain files:\n"
           "      wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz\n"
           "      wget http://hgdownload.cse.ucsc.edu/goldenPath/hs1/liftOver/hg38-chm13v2.over.chain.gz\n"
           "\n";
}

static FILE *get_file_handle(const char *str) {
    FILE *ret;
    if (strcmp(str, "-") == 0)
        ret = stdout;
    else {
        ret = fopen(str, "w");
        if (!ret) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

static inline int bcf_is_symbolic(const bcf1_t *rec) {
    for (int i = 0; i < rec->n_allele; i++)
        if (rec->d.allele[i][0] == '<' || rec->d.allele[i][0] == '*') return 1;
    return 0;
}

static int read_tags(const bcf_hdr_t *hdr, const char *tags, char data_sep, int length_flag, int field_type_flag,
                     int **tags_int_id, int **tags_coltype, float **tags_sum) {
    char *s = strdup(tags);
    int moff = 0, *off = NULL;
    int n = ksplit_core(s, ',', &moff, &off);
    *tags_int_id = (int *)malloc(sizeof(int) * n);
    *tags_coltype = (int *)malloc(sizeof(int) * n);
    if (data_sep) *tags_sum = (float *)malloc(sizeof(float) * n);
    int n_tags = 0;
    for (int i = 0; i < n; i++) {
        char *ss = &s[off[i]];
        char *ptr;
        if (data_sep) {
            ptr = strchr(ss, data_sep);
            if (!ptr) error("The tag \"%s\" does not have information defined by the \"%c\" character\n", ss, data_sep);
            *ptr = '\0';
        }
        int coltype = -1;
        if (!strncasecmp("INFO/", ss, 5)) {
            coltype = BCF_HL_INFO;
            ss += 5;
        } else if (!strncasecmp("INF/", ss, 4)) {
            coltype = BCF_HL_INFO;
            ss += 4;
        } else if (!strncasecmp("FORMAT/", ss, 7)) {
            coltype = BCF_HL_FMT;
            ss += 7;
        } else if (!strncasecmp("FMT/", ss, 4)) {
            coltype = BCF_HL_FMT;
            ss += 4;
        }
        if (coltype == -1) error("The tag \"%s\" is not specified as either an INFO or a FORMAT tag\n", ss);
        int int_id = bcf_hdr_id2int(hdr, BCF_DT_ID, ss);
        if (!bcf_hdr_idinfo_exists(hdr, coltype, int_id)) continue;
        int length = bcf_hdr_id2length(hdr, coltype, int_id);
        if (!(1 << length & length_flag))
            error("The %s tag \"%s\" is not the correct AGR tag\n", coltype == BCF_HL_INFO ? "INFO" : "FORMAT", ss);
        int field_type = bcf_hdr_id2type(hdr, coltype, int_id);
        if (!(1 << field_type & field_type_flag))
            error("The %s tag \"%s\" is not the correct field type tag\n", coltype == BCF_HL_INFO ? "INFO" : "FORMAT",
                  ss);
        (*tags_int_id)[n_tags] = int_id;
        (*tags_coltype)[n_tags] = coltype;
        if (data_sep) {
            char *tmp;
            (*tags_sum)[n_tags] = strtof(++ptr, &tmp);
            if (*tmp) error("Could not parse float: %s\n", ptr);
        }
        n_tags++;
    }
    free(off);
    free(s);
    return n_tags;
}

// find all tags with a given length
static int find_tags(const bcf_hdr_t *hdr, int length_flag, int field_type_flag, int **tags_int_id,
                     int **tags_coltype) {
    const int coltypes[2] = {BCF_HL_INFO, BCF_HL_FMT};
    int n_tags = 0;
    int m_tags_int_id = 0;
    int m_tags_coltype = 0;
    *tags_int_id = NULL;
    *tags_coltype = NULL;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < hdr->n[BCF_DT_ID]; j++) {
            if (!hdr->id[BCF_DT_ID][j].val || !bcf_hdr_idinfo_exists(hdr, coltypes[i], j)) continue;
            int length = bcf_hdr_id2length(hdr, coltypes[i], j);
            int field_type = bcf_hdr_id2type(hdr, coltypes[i], j);
            if (1 << length & length_flag && 1 << field_type & field_type_flag) {
                hts_expand(int, n_tags + 1, m_tags_int_id, *tags_int_id);
                hts_expand(int, n_tags + 1, m_tags_coltype, *tags_coltype);
                (*tags_int_id)[n_tags] = j;
                (*tags_coltype)[n_tags] = coltypes[i];
                n_tags++;
            }
        }
    }
    return n_tags;
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out) {
    args = (args_t *)calloc(1, sizeof(args_t));
    args->in_hdr = in;
    args->out_hdr = out;
    args->aln_win = 100;
    args->indel_win = 250;
    args->tags.flip_tag = FLIP_TAG;
    args->tags.swap_tag = SWAP_TAG;
    int cache_size = 0;
    int output_type = FT_VCF;
    int clevel = -1;
    int max_snp_gap = 1;    // maximum distance between two contiguous blocks to allow merging
    int max_indel_gap = 20; // maximum distance between two contiguous blocks to allow padding of alleles
    char *tmp = NULL;
    const char *src_ref_fname = NULL;
    const char *dst_ref_fname = NULL;
    const char *chain_fname = NULL;
    const char *blocks_fname = NULL;
    const char *reject_fname = NULL;
    char *drop_tags = DROP_TAGS;
    char *reverse_tags = REVERSE_TAGS;
    char *flip_tags = FLIP_TAGS;
    char *genotype_tags = GENOTYPE_TAGS;

    static struct option loptions[] = {{"src-fasta-ref", required_argument, NULL, 's'},
                                       {"fasta-ref", required_argument, NULL, 'f'},
                                       {"set-cache-size", required_argument, NULL, 1},
                                       {"chain", required_argument, NULL, 'c'},
                                       {"max-snp-gap", required_argument, NULL, 2},
                                       {"max-indel-gap", required_argument, NULL, 3},
                                       {"indel-win", required_argument, NULL, 4},
                                       {"lift-mt", no_argument, NULL, 5},
                                       {"print-blocks", required_argument, NULL, 6},
                                       {"no-left-align", no_argument, NULL, 7},
                                       {"reject", required_argument, NULL, 8},
                                       {"reject-type", required_argument, NULL, 'O'},
                                       {"write-src", no_argument, NULL, 9},
                                       {"flip-tag", required_argument, NULL, 10},
                                       {"swap-tag", required_argument, NULL, 11},
                                       {"tags-to-drop", required_argument, NULL, 12},
                                       {"tags-to-reverse", required_argument, NULL, 13},
                                       {"tags-to-flip", required_argument, NULL, 14},
                                       {"tags-genotype", required_argument, NULL, 15},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?s:f:c:O:", loptions, NULL)) >= 0) {
        switch (c) {
        case 's':
            src_ref_fname = optarg;
            break;
        case 'f':
            dst_ref_fname = optarg;
            break;
        case 1:
            cache_size = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --set-cache-size %s\n", optarg);
            break;
        case 'c':
            chain_fname = optarg;
            break;
        case 2:
            max_snp_gap = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --max-snp-gap %s\n", optarg);
            break;
        case 3:
            max_indel_gap = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --max-indel-gap %s\n", optarg);
            break;
        case 4:
            args->indel_win = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --indel-win %s\n", optarg);
            break;
        case 5:
            args->lift_mt = 1;
            break;
        case 6:
            blocks_fname = optarg;
            break;
        case 7:
            args->no_left_align = 1;
            break;
        case 8:
            reject_fname = optarg;
            break;
        case 'O':
            switch (optarg[0]) {
            case 'b':
                output_type = FT_BCF_GZ;
                break;
            case 'u':
                output_type = FT_BCF;
                break;
            case 'z':
                output_type = FT_VCF_GZ;
                break;
            case 'v':
                output_type = FT_VCF;
                break;
            default: {
                clevel = strtol(optarg, &tmp, 10);
                if (*tmp || clevel < 0 || clevel > 9) error("The output type \"%s\" not recognised\n", optarg);
            }
            }
            if (optarg[1]) {
                clevel = strtol(optarg + 1, &tmp, 10);
                if (*tmp || clevel < 0 || clevel > 9)
                    error("Could not parse argument: --compression-level %s\n", optarg + 1);
            }
            break;
        case 9:
            args->write_src = 1;
            break;
        case 10:
            args->tags.flip_tag = optarg;
            break;
        case 11:
            args->tags.swap_tag = optarg;
            break;
        case 12:
            drop_tags = optarg;
            break;
        case 13:
            reverse_tags = optarg;
            break;
        case 14:
            flip_tags = optarg;
            break;
        case 15:
            genotype_tags = optarg;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage());
            break;
        }
    }

    if (!in || !out) error("Expected input VCF\n%s", usage());

    // load target reference file
    if (src_ref_fname) {
        args->src_fai = fai_load(src_ref_fname);
        if (!args->src_fai) error("Could not load the reference %s\n", src_ref_fname);
        if (cache_size) fai_set_cache_size(args->src_fai, cache_size);
    }

    // does not perform liftover ... only performs the indel extension
    if (!dst_ref_fname || !chain_fname) return 0;

    // load query reference file
    args->dst_fai = fai_load(dst_ref_fname);
    if (!args->dst_fai) error("Could not load the reference %s\n", dst_ref_fname);
    if (cache_size) fai_set_cache_size(args->dst_fai, cache_size);

    // reset contig table for the output header
    bcf_hdr_remove(out, BCF_HL_CTG, NULL);
    // required to reset the contigs dictionary and table
    kh_clear(vdict, out->dict[BCF_DT_CTG]);
    for (int i = 0; i < out->n[BCF_DT_CTG]; i++) free((void *)out->id[BCF_DT_CTG][i].key);
    out->n[BCF_DT_CTG] = 0;

    int n = faidx_nseq(args->dst_fai);
    for (int i = 0; i < n; i++) {
        const char *seq = faidx_iseq(args->dst_fai, i);
        int len = faidx_seq_len(args->dst_fai, seq);
        bcf_hdr_printf(out, "##contig=<ID=%s,length=%d>", seq, len);
    }
    if (bcf_hdr_printf(
            out, "##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Whether alleles flipped strand during liftover\">",
            args->tags.flip_tag)
        < 0)
        error_errno("Failed to add \"%s\" INFO header", args->tags.flip_tag);
    if (bcf_hdr_printf(out,
                       "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Which alternate allele became the reference "
                       "during liftover (-1 for new reference)\">",
                       args->tags.swap_tag)
        < 0)
        error_errno("Failed to add \"%s\" INFO header", args->tags.swap_tag);
    if (args->write_src) {
        if (bcf_hdr_printf(out,
                           "##INFO=<ID=SRC_CHROM,Number=1,Type=String,Description=\"The name of the source contig of "
                           "the variant prior to liftover\">")
            < 0)
            error_errno("Failed to add \"SRC_CHROM\" INFO header");
        if (bcf_hdr_printf(out,
                           "##INFO=<ID=SRC_POS,Number=1,Type=Integer,Description=\"The position of the variant on the "
                           "source contig prior to liftover\">")
            < 0)
            error_errno("Failed to add \"SRC_POS\" INFO header");
        if (bcf_hdr_printf(out,
                           "##INFO=<ID=SRC_REF_ALT,Number=.,Type=String,Description=\"A list of the original alleles "
                           "of the variant prior to liftover\">")
            < 0)
            error_errno("Failed to add \"SRC_REF_ALT\" INFO header");
    }
    if (bcf_hdr_sync(out) < 0) error_errno("Failed to update header");

    // load chain file
    if (max_snp_gap > 20)
        fprintf(
            stderr,
            "Warning: merging contiguous blocks farther than 20 bp apart (--max-snp-gap %d used) is not recommended\n",
            max_snp_gap);
    if (max_indel_gap > 20)
        fprintf(stderr,
                "Warning: padding alleles across contiguous blocks farther than 20 bp apart (--max-indel-gap %d used) "
                "is not recommended\n",
                max_indel_gap);

    htsFile *fp = hts_open(chain_fname, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fp->fn, strerror(errno));
    args->n_chains = read_chains(fp, in, out, max_snp_gap, max_indel_gap, &args->chains, &args->blocks);
    if (hts_close(fp) < 0) error("Close failed: %s\n", fp->fn);

    if (blocks_fname) {
        FILE *blocks_file = get_file_handle(blocks_fname);
        write_chains(blocks_file, in, out, args->chains, args->n_chains, args->blocks);
        if (blocks_file != stdout && blocks_file != stderr) fclose(blocks_file);
    }

    args->idx = regidx_init_chains(in, args->chains, args->n_chains, args->blocks);
    args->itr = regitr_init(args->idx);

    tags_t *tags = &args->tags;
    tags->gt_id = bcf_hdr_id2int(in, BCF_DT_ID, "GT");
    tags->n_drop_tags = read_tags(in, drop_tags, '\0',
                                  1 << BCF_VL_FIXED | 1 << BCF_VL_VAR | 1 << BCF_VL_A | 1 << BCF_VL_G | 1 << BCF_VL_R,
                                  1 << BCF_HT_FLAG | 1 << BCF_HT_INT | 1 << BCF_HT_REAL | 1 << BCF_HT_STR,
                                  &tags->drop_tags_int_id, &tags->drop_tags_coltype, NULL);
    tags->n_reverse_tags = read_tags(in, reverse_tags, ':', 1 << BCF_VL_A, 1 << BCF_HT_REAL, &tags->reverse_tags_int_id,
                                     &tags->reverse_tags_coltype, &tags->reverse_tags_sum);
    tags->n_flip_tags =
        read_tags(in, flip_tags, '\0', 1 << BCF_VL_A, 1 << BCF_HT_INT | 1 << BCF_HT_REAL | 1 << BCF_HT_STR,
                  &tags->flip_tags_int_id, &tags->flip_tags_coltype, NULL);
    tags->n_genotype_tags = read_tags(
        in, genotype_tags, '\0', 1 << BCF_VL_FIXED | 1 << BCF_VL_VAR | 1 << BCF_VL_A | 1 << BCF_VL_G | 1 << BCF_VL_R,
        1 << BCF_HT_INT, &tags->genotype_tags_int_id, &tags->genotype_tags_coltype, NULL);
    tags->n_GR_tags = find_tags(in, 1 << BCF_VL_G | 1 << BCF_VL_R, 1 << BCF_HT_INT | 1 << BCF_HT_REAL | 1 << BCF_HT_STR,
                                &tags->GR_tags_int_id, &tags->GR_tags_coltype);

    // mitochrondria is handled in a special way as GRCh37 has a different mitochondria than hg19
    args->in_mt_rid = bcf_hdr_name2id_flexible(in, "MT");
    args->out_mt_rid = bcf_hdr_name2id_flexible(out, "chrM");
    if (!args->lift_mt && args->in_mt_rid >= 0 && args->out_mt_rid >= 0) {
        int in_mt_length = in->id[BCF_DT_CTG][args->in_mt_rid].val->info[0];
        int out_mt_length = out->id[BCF_DT_CTG][args->out_mt_rid].val->info[0];
        if (in_mt_length == 0 || in_mt_length != out_mt_length) args->lift_mt = 1;
    }

    if (reject_fname) {
        char wmode[8];
        set_wmode(wmode, output_type, (char *)reject_fname, clevel);
        args->reject_fh = hts_open(reject_fname, hts_bcf_wmode(output_type));
        if (args->reject_fh == NULL || bcf_hdr_write(args->reject_fh, args->in_hdr) < 0)
            error("Error: cannot write to \"%s\": %s\n", args->reject_fh->fn, strerror(errno));
    }

    return 0;
}

/****************************************
 * ALLELE EXTENSION FRAMEWORK           *
 ****************************************/

// class to realign alleles using the reference
typedef struct {
    const bcf_hdr_t *hdr;
    bcf1_t *rec;
    kstring_t *als;
    faidx_t *fai;
    int aln_win;
    char *ref;
    hts_pos_t beg;
    hts_pos_t end;
} bcf1_realign_t;

// Petr Danecek's code from bcftools/vcfnorm.c
static inline void seq_to_upper(char *seq, int len) {
    int i;
    for (i = 0; i < len; i++) seq[i] = nt_to_upper(seq[i]);
}

// Petr Danecek's code from bcftools/vcfnorm.c
static inline int replace_iupac_codes(char *seq, int nseq) {
    // Replace ambiguity codes with N for now, it awaits to be seen what the VCF spec codifies in the end
    int i, n = 0;
    for (i = 0; i < nseq; i++) {
        char c = toupper(seq[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            seq[i] = 'N';
            n++;
        }
    }
    return n;
}

static void fetch_seq64(bcf1_realign_t *this) {
    free(this->ref);
    hts_pos_t len;
    this->ref = faidx_fetch_seq64(this->fai, (char *)bcf_seqname(this->hdr, this->rec), this->beg, this->end, &len);
    if (!this->ref || len == 1)
        error("faidx_fetch_seq failed at %s:%" PRIhts_pos "\n", bcf_seqname(this->hdr, this->rec), this->beg + 1);
    seq_to_upper(this->ref, len);
    replace_iupac_codes(this->ref, len);
}

static void shift_left(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    if (rec->n_allele == 1) return;
    while (rec->pos > 0) {
        char last_ref = this->als[0].s[this->als[0].l - 1];
        for (int i = 1; i < rec->n_allele; i++)
            if (this->als[i].s[this->als[i].l - 1] != last_ref) return;
        if (this->beg == 0 || rec->pos - 1 < this->beg) {
            this->beg = rec->pos - this->aln_win;
            fetch_seq64(this);
        }
        char first_ref = this->ref[rec->pos - this->beg - 1];
        for (int i = 0; i < rec->n_allele; i++) {
            memmove(this->als[i].s + 1, this->als[i].s, this->als[i].l);
            this->als[i].s[0] = first_ref;
        }
        rec->pos--;
    }
}

static void trim_left(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    if (this->als[0].l <= 1) return;
    int i, j;
    for (i = 0; i < this->als[0].l - 1; i++) {
        char first_ref = this->als[0].s[i];
        for (j = 1; j < rec->n_allele; j++)
            if (i >= this->als[j].l - 1 || this->als[j].s[i] != first_ref) break;
        if (j < rec->n_allele) break;
    }
    for (j = 0; j < rec->n_allele; j++) {
        this->als[j].l -= i;
        memmove(this->als[j].s, this->als[j].s + i, this->als[j].l);
    }
    rec->pos += i;
}

static void trim_right(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    if (this->als[0].l <= 1) return;
    int i, j;
    for (i = 0; i < this->als[0].l - 1; i++) {
        char last_ref = this->als[0].s[this->als[0].l - 1 - i];
        for (j = 1; j < rec->n_allele; j++)
            if (i >= this->als[j].l - 1 || this->als[j].s[this->als[j].l - 1 - i] != last_ref) break;
        if (j < rec->n_allele) break;
    }
    for (j = 0; j < rec->n_allele; j++) this->als[j].l -= i;
}

static void pad_left(bcf1_realign_t *this, int npad) {
    bcf1_t *rec = this->rec;
    if (this->beg == 0 || rec->pos - npad < this->beg) {
        this->beg = rec->pos - npad;
        fetch_seq64(this);
    }
    const char *ptr = &this->ref[rec->pos - this->beg - npad];
    for (int i = 0; i < rec->n_allele; i++) {
        ks_resize(&this->als[i], this->als[i].l + npad);
        memmove(this->als[i].s + npad, this->als[i].s, this->als[i].l);
        memcpy(this->als[i].s, ptr, npad);
        this->als[i].l += npad;
    }
    rec->pos -= npad;
}

static void pad_right(bcf1_realign_t *this, int npad) {
    bcf1_t *rec = this->rec;
    if (this->end < rec->pos + this->als[0].l + npad - 1) {
        this->end = rec->pos + this->als[0].l + npad - 1;
        fetch_seq64(this);
    }
    const char *ptr = &this->ref[rec->pos - this->beg + this->als[0].l];
    for (int i = 0; i < rec->n_allele; i++) kputsn(ptr, npad, &this->als[i]);
}

static void pad_from_right(bcf1_realign_t *this, const char *s_ptr, int d) {
    bcf1_t *rec = this->rec;
    int npad = 0;
    const char *l_ptr = &this->ref[rec->pos - this->beg + this->als[0].l];
    while (1) {
        if (rec->pos + this->als[0].l + npad > this->end) {
            // extract more sequence from the reference
            this->end += this->aln_win;
            fetch_seq64(this);
            l_ptr = &this->ref[rec->pos - this->beg + this->als[0].l + npad];
            if (npad >= d) s_ptr = &this->ref[rec->pos - this->beg + this->als[0].l + npad - d];
        }
        if (npad == d) s_ptr = &this->ref[rec->pos - this->beg + this->als[0].l];
        npad++;
        if (*l_ptr != *s_ptr) break;
        l_ptr++;
        s_ptr++;
    }
    const char *ptr = &this->ref[rec->pos - this->beg + this->als[0].l];
    for (int i = 0; i < rec->n_allele; i++) kputsn(ptr, npad, &this->als[i]);
}

static void pad_from_left(bcf1_realign_t *this, const char *s_ptr, int d) {
    bcf1_t *rec = this->rec;
    int npad = 0;
    const char *l_ptr = &this->ref[rec->pos - this->beg - 1];
    while (1) {
        if (rec->pos - npad <= this->beg) {
            // extract more sequence from the reference
            this->beg -= this->aln_win;
            fetch_seq64(this);
            l_ptr = &this->ref[rec->pos - this->beg - npad - 1];
            if (npad >= d) s_ptr = &this->ref[rec->pos - this->beg - npad + d - 1];
        }
        if (npad == d) s_ptr = &this->ref[rec->pos - this->beg - 1];
        npad++;
        if (*l_ptr != *s_ptr) break;
        l_ptr--;
        s_ptr--;
    }
    for (int i = 0; i < rec->n_allele; i++) {
        ks_resize(&this->als[i], this->als[i].l + npad);
        memmove(this->als[i].s + npad, this->als[i].s, this->als[i].l);
        memcpy(this->als[i].s, l_ptr, npad);
        this->als[i].l += npad;
    }
    rec->pos -= npad;
}

static void initialize_alleles(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    for (int i = 0; i < rec->n_allele; i++) {
        this->als[i].l = 0;
        kputsn_(rec->d.allele[i], strlen(rec->d.allele[i]), &this->als[i]);
    }

    this->beg = rec->pos - 1;
    this->end = rec->pos + this->als[0].l;
    fetch_seq64(this);

    // make sure the reference allele matches the reference
    if (strncasecmp(&this->ref[rec->pos - this->beg], this->als[0].s, this->als[0].l) != 0)
        error("Error: the reference allele %s does not match the reference at %s:%" PRIhts_pos "\n", this->als[0].s,
              bcf_seqname(this->hdr, rec), this->beg + 1);
}

// updated the alleles in the bcf1_t structure
static void update_alleles(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    kstring_t *str = &this->als[0];
    for (int i = 1; i < rec->n_allele; i++) {
        kputc_(',', str);
        kputsn_(this->als[i].s, this->als[i].l, str);
    }
    kputc_('\0', str);
    str->l--;
    bcf_update_alleles_str(this->hdr, rec, str->s);
    free(this->ref);
}

// inspired by Petr Danecek's code from realign() in bcftools/vcfnorm.c
static void extend_alleles(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    initialize_alleles(this);

    // find allele pairs that need to be extended
    for (int i = 0; i < rec->n_allele; i++) {
        for (int j = i + 1; j < rec->n_allele; j++) {
            if (this->als[i].s[0] != this->als[j].s[0]) pad_left(this, 1);
            if (this->als[i].s[this->als[i].l - 1] != this->als[j].s[this->als[j].l - 1]) pad_right(this, 1);
            if (this->als[i].l == this->als[j].l) continue;
            kstring_t *l_str = this->als[i].l > this->als[j].l ? &this->als[i] : &this->als[j];
            kstring_t *s_str = this->als[i].l > this->als[j].l ? &this->als[j] : &this->als[i];
            if (strncasecmp(l_str->s, s_str->s, s_str->l) == 0)
                pad_from_right(this, l_str->s + s_str->l, l_str->l - s_str->l);
            if (strncasecmp(l_str->s + l_str->l - s_str->l, s_str->s, s_str->l) == 0) {
                pad_from_left(this, l_str->s + l_str->l - s_str->l - 1, l_str->l - s_str->l);
            }
        }
    }

    // updated the alleles in the bcf1_t structure
    update_alleles(this);
}

/****************************************
 * PROCESS RECORDS                      *
 ****************************************/

// check whether alleles in the record require extension to avoid ambiguity
static int is_extension_needed(bcf1_t *rec, void **tmp_arr, int *m_tmp_arr) {
    hts_expand(char, sizeof(int) * rec->n_allele, *m_tmp_arr, *tmp_arr);
    int *len = (int *)(*tmp_arr);
    for (int i = 0; i < rec->n_allele; i++) len[i] = strlen(rec->d.allele[i]);

    for (int i = 0; i < rec->n_allele; i++) {
        for (int j = i + 1; j < rec->n_allele; j++) {
            if (rec->d.allele[i][0] != rec->d.allele[j][0]) return 1;
            if (rec->d.allele[i][len[i] - 1] != rec->d.allele[j][len[j] - 1]) return 1;
            if (len[i] == len[j]) continue;
            int ret;
            if (len[i] > len[j]) {
                ret = (strncasecmp(rec->d.allele[i], rec->d.allele[j], len[j]) == 0)
                      || (strncasecmp(rec->d.allele[i] + len[i] - len[j], rec->d.allele[j], len[j]) == 0);
            } else {
                ret = (strncasecmp(rec->d.allele[i], rec->d.allele[j], len[i]) == 0)
                      || (strncasecmp(rec->d.allele[i], rec->d.allele[j] + len[j] - len[i], len[i]) == 0);
            }
            if (ret) return ret;
        }
    }
    return 0;
}

// Tan, A., Abecasis, G. R., Kang, H. M., Unified representation of genetic variants.
// Bioinformatics, (2015) https://doi.org/10.1093/bioinformatics/btv112
static int is_left_aligned(bcf1_t *rec, void **tmp_arr, int *m_tmp_arr) {
    if (rec->n_allele == 1) return strlen(rec->d.allele[0]) <= 1;

    hts_expand(char, sizeof(int) * rec->n_allele, *m_tmp_arr, *tmp_arr);
    int *len = (int *)(*tmp_arr);
    int i, min_len_is_one = 0;
    for (i = 0; i < rec->n_allele; i++) {
        len[i] = strlen(rec->d.allele[i]);
        if (len[i] == 1) min_len_is_one = 1;
    }

    // 1. The alleles end with at least two different nucleotides
    char last_ref = rec->d.allele[0][len[0] - 1];
    for (i = 0; i < rec->n_allele; i++)
        if (rec->d.allele[i][len[i] - 1] != last_ref) break;
    if (i == rec->n_allele) return 0;

    // 2. The alleles start with at least two different nucleotides, or the shortest allele has length 1
    if (!min_len_is_one) {
        char first_ref = rec->d.allele[0][0];
        for (i = 0; i < rec->n_allele; i++)
            if (rec->d.allele[i][0] != first_ref) break;
        if (i == rec->n_allele) return 0;
    }

    return 1;
}

// t_pos has to be provided in 1-based coordinates
// if successful it also returns the number of bases until the next contiguous block
static int liftover_snp(regidx_t *idx, regitr_t *itr, const char *chr, hts_pos_t t_pos, int t_strand, int *q_rid,
                        hts_pos_t *q_pos, int *q_strand) {
    int ret = -1;
    if (regidx_overlap(idx, chr, (uint32_t)t_pos, (uint32_t)t_pos, itr)) {
        for (int i = 0; regitr_overlap(itr); i++) {
            if (i > 0)
                fprintf(stderr, "Warning: more than one contiguous block overlaps with position %s:%" PRIhts_pos "\n",
                        chr, t_pos);
            int block_ind = regitr_payload(itr, int);
            block_t *block = &args->blocks[block_ind];
            chain_t *chain = &args->chains[block->chain_ind];

            *q_rid = chain->q_rid;
            *q_strand = chain->qStrand;
            int block_pos = t_pos - itr->beg;
            if (*q_strand) // - strand
                *q_pos = (hts_pos_t)(chain->qSize - chain->qStart - block->qStart - block_pos - 1);
            else // + strand
                *q_pos = (hts_pos_t)(chain->qStart + block->qStart + block_pos);
            if (t_strand) // pad sequence to the left
                ret = block->tStart_gap < 0 ? 0 : t_pos - itr->beg + block->tStart_gap + 1;
            else // pad sequence to the right
                ret = block->tEnd_gap < 0 ? 0 : itr->end - t_pos + block->tEnd_gap + 1;
        }
    }
    return ret;
}

static int liftover_indel(regidx_t *idx, regitr_t *itr, const char *chr, hts_pos_t t_beg, int *len, int indel_win,
                          int *rid, hts_pos_t *pos, int *strand, int *npad) {
    int left_rid = -1, right_rid = -1;
    hts_pos_t left_pos = -1, right_pos = -1;
    int left_strand = 0, right_strand = 0;

    int ret_left = liftover_snp(idx, itr, chr, t_beg, 0, &left_rid, &left_pos, &left_strand);
    int ret_right = liftover_snp(idx, itr, chr, t_beg + *len - 1, 1, &right_rid, &right_pos, &right_strand);
    if (ret_left < 0 && ret_right < 0) return -1; // both edges of indel failed to lift over

    // strategy to pad alleles to get to the very first base on the next contiguous block
    if (ret_left < 0 && ret_right > 0 && ret_right < indel_win) { // pad sequence to the left
        *npad = -(ret_right - *len + 1);
        //    assert(liftover_snp(idx, itr, chr, t_beg + *npad + 1, 0, &left_rid, &left_pos, &left_strand) < 0);
        ret_left = liftover_snp(idx, itr, chr, t_beg + *npad, 0, &left_rid, &left_pos, &left_strand);
    } else if (ret_left > 0 && ret_left < indel_win && ret_right < 0) { // pad sequence to the right
        *npad = ret_left - *len + 1;
        //    assert(liftover_snp(idx, itr, chr, t_beg + ret_left - 1, 1, &right_rid, &right_pos, &right_strand) < 0);
        ret_right = liftover_snp(idx, itr, chr, t_beg + ret_left, 1, &right_rid, &right_pos, &right_strand);
    } else if (ret_left < 0 || ret_right < 0) { // only one edge lifted over
        return -1;
    }

    if (left_rid != right_rid) return -1;
    if (abs(right_pos - left_pos) > indel_win) return -1;
    if (left_strand != right_strand) return -1;

    *rid = left_strand == 0 ? left_rid : right_rid;
    *pos = left_strand == 0 ? left_pos : right_pos;
    *strand = left_strand;
    *len = (left_strand ? left_pos - right_pos : right_pos - left_pos) + 1;
    return 0;
}

// flip all alleles, except those that are symbolic
static void flip_alleles(bcf1_t *rec, const bcf_hdr_t *hdr, kstring_t *str) {
    str->l = 0;
    for (int i = 0; i < rec->n_allele; i++) {
        char *ptr, *allele = rec->d.allele[i];
        int len = strlen(allele);
        if (allele[0] == '<' || allele[0] == '*') {
            kputsn_(allele, len, str);
        } else {
            ks_resize(str, str->l + len);
            for (ptr = allele + len - 1; ptr >= allele; ptr--) str->s[str->l++] = revnt[*ptr & 0x7F];
        }
        kputc_(',', str);
    }
    str->l--;
    str->s[str->l] = '\0';
    bcf_update_alleles_str(hdr, rec, str->s);
}

static int find_reference(bcf1_t *rec, const bcf_hdr_t *hdr, const faidx_t *fai, int ref_len, kstring_t *str) {
    const char *chr = bcf_seqname(hdr, rec);
    hts_pos_t len;
    char *ref = faidx_fetch_seq64(fai, chr, rec->pos, rec->pos + ref_len - 1, &len);
    if (!ref) error("faidx_fetch_seq failed at %s:%" PRIhts_pos "\n", bcf_seqname(hdr, rec), rec->pos + 1);
    seq_to_upper(ref, len);
    replace_iupac_codes(ref, len);

    int swap = -1;
    for (int i = 0; i < rec->n_allele; i++) {
        int len = strlen(rec->d.allele[i]);
        if (strncasecmp(ref, rec->d.allele[i], len) == 0) {
            if (swap >= 0) {
                fprintf(stderr,
                        "Warning: as option --src-fasta-ref is missing it is impossible to infer which allele is the "
                        "reference allele at position %s:%" PRIhts_pos "\n",
                        bcf_seqname(hdr, rec), rec->pos + 1);
                break;
            } else {
                swap = i;
            }
        }
    }

    // update alleles if necessary
    if (swap) {
        if (swap == -1) {
            str->l = 0;
            kputsn_(ref, ref_len, str);
            kputc_(',', str);
            for (int i = 0; i < rec->n_allele; i++) {
                char *allele = rec->d.allele[i];
                int len = strlen(allele);
                kputsn_(allele, len, str);
                kputc_(',', str);
            }
            str->l--;
            str->s[str->l] = '\0';
            bcf_update_alleles_str(hdr, rec, str->s);
        } else {
            char *tmp = rec->d.allele[0];
            rec->d.allele[0] = rec->d.allele[swap];
            rec->d.allele[swap] = tmp;
            bcf_update_alleles(args->out_hdr, rec, (const char **)rec->d.allele, rec->n_allele);
        }
    }

    free(ref);
    return swap;
}

// this function updates FORMAT/GT records if the number or order of alleles has changed
static void update_genotypes(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int swap, void **tmp_arr, int *m_tmp_arr) {
    if (swap < 0) {
        int ngt = bcf_get_genotypes(args->out_hdr, rec, &args->tmp_arr, &args->m_tmp_arr);
        if (ngt <= 0) return;
        int *gts = (int *)(args->tmp_arr);
        for (int i = 0; i < ngt; i++)
            if (!bcf_gt_is_missing(gts[i]) && !(gts[i] == bcf_int32_vector_end)) gts[i] += 2;
        bcf_update_genotypes(hdr, rec, gts, ngt);
    } else {
        bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
        if (!fmt) return;
        int n = rec->n_sample * fmt->n;
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (int i = 0; i < n; i++) {                                                                                  \
            if (bcf_gt_allele(p[i]) == 0)                                                                              \
                p[i] += 2 * swap;                                                                                      \
            else if (bcf_gt_allele(p[i]) == swap)                                                                      \
                p[i] -= 2 * swap;                                                                                      \
        }                                                                                                              \
    }
        switch (fmt->type) {
        case BCF_BT_INT8:
            BRANCH(int8_t);
            break;
        case BCF_BT_INT16:
            BRANCH(int16_t);
            break;
        case BCF_BT_INT32:
            BRANCH(int32_t);
            break;
        default:
            error("Unexpected type %d\n", fmt->type);
        }
#undef BRANCH
    }
}

static void update_genotype_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int swap, void **tmp_arr,
                                   int *m_tmp_arr) {
    int n, type, len;
    void *p;
    const char *key = bcf_hdr_int2id(hdr, BCF_DT_ID, int_id);
    int field_type = bcf_hdr_id2type(hdr, coltype, int_id);

    if (coltype == BCF_HL_INFO) { // INFO fields
        n = 1;
        bcf_info_t *info = bcf_get_info_id(rec, int_id);
        if (!info) return;
        len = info->len;
        type = info->type;
        p = info->vptr;
    } else { // FORMAT fields
        n = rec->n_sample;
        bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
        if (!n || !fmt) return;
        len = fmt->n;
        type = fmt->type;
        p = fmt->p;
    }
    assert(field_type == BCF_HT_INT);

    hts_expand(char, sizeof(int32_t) * n * len, *m_tmp_arr, *tmp_arr);

    if (swap < 0) {
#define BRANCH(type_t, is_missing, is_vector_end)                                                                      \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        int32_t *dst = (int32_t *)(*tmp_arr);                                                                          \
        for (int i = 0; i < n * len; i++) {                                                                            \
            if (is_missing)                                                                                            \
                *dst = bcf_int32_missing;                                                                              \
            else if (is_vector_end)                                                                                    \
                *dst = bcf_int32_vector_end;                                                                           \
            else                                                                                                       \
                *dst = *src + 1;                                                                                       \
            src++;                                                                                                     \
            dst++;                                                                                                     \
        }                                                                                                              \
    }
        switch (type) {
        case BCF_BT_INT8:
            BRANCH(int8_t, *src == bcf_int8_missing, *src == bcf_int8_vector_end);
            break;
        case BCF_BT_INT16:
            BRANCH(int16_t, *src == bcf_int16_missing, *src == bcf_int16_vector_end);
            break;
        case BCF_BT_INT32:
            BRANCH(int32_t, *src == bcf_int32_missing, *src == bcf_int32_vector_end);
            break;
        default:
            fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
            exit(1);
        }
#undef BRANCH
        if (coltype == BCF_HL_INFO) // INFO fields
            bcf_update_info_int32(hdr, rec, key, *tmp_arr, len);
        else // FORMAT fields
            bcf_update_format_int32(hdr, rec, key, *tmp_arr, len);
    } else {
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (int i = 0; i < n * len; i++) {                                                                            \
            if (ptr[i] == 0)                                                                                           \
                ptr[i] = swap;                                                                                         \
            else if (ptr[i] == swap)                                                                                   \
                ptr[i] = 0;                                                                                            \
        }                                                                                                              \
        if (coltype == BCF_HL_INFO && len == 1) {                                                                      \
            bcf_info_t *info = bcf_get_info_id(rec, int_id);                                                           \
            info->v1.i = ((type_t *)info->vptr)[0];                                                                    \
        }                                                                                                              \
    }
        switch (type) {
        case BCF_BT_INT8:
            BRANCH(int8_t);
            break;
        case BCF_BT_INT16:
            BRANCH(int16_t);
            break;
        case BCF_BT_INT32:
            BRANCH(int32_t);
            break;
        default:
            error("Unexpected type %d\n", type);
        }
#undef BRANCH
    }
}

// this function updates INFO/FORMAT Number=G/Number=R records if the number or order of alleles has changed
// it does not currently support Number=G for Type=Character/Type=String
static void update_AGR_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int swap, void **tmp_arr,
                              int *m_tmp_arr) {
    int n, type, len;
    void *p;
    int n_als = rec->n_allele;
    int length = bcf_hdr_id2length(hdr, coltype, int_id);
    const char *key = bcf_hdr_int2id(hdr, BCF_DT_ID, int_id);
    int field_type = bcf_hdr_id2type(hdr, coltype, int_id);
    char *s;

    if (coltype == BCF_HL_INFO) { // INFO fields
        n = 1;
        bcf_info_t *info = bcf_get_info_id(rec, int_id);
        if (!info) return;
        len = info->len;
        type = info->type;
        p = info->vptr;
    } else { // FORMAT fields
        n = rec->n_sample;
        bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
        if (!n || !fmt) return;
        len = fmt->n;
        type = fmt->type;
        p = fmt->p;
    }

    int size = type == BCF_BT_CHAR ? n * (len + 2) + 1 : sizeof(int) * n * (len + n_als);
    hts_expand(char, size, *m_tmp_arr, *tmp_arr);

    if (swap < 0) {               // add reference allele entry
        if (length == BCF_VL_G) { // define how to add an element to a VL_G array
            if (type != BCF_BT_CHAR) assert(len == (n_als - 1) * n_als / 2);
#define BRANCH(type_t, set_missing, set_vector_end, is_missing, is_vector_end, out_type_t)                             \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        out_type_t *dst = (out_type_t *)(*tmp_arr);                                                                    \
        for (int i = 0; i < n; i++) {                                                                                  \
            for (int j = 0; j < n_als; j++) {                                                                          \
                set_missing;                                                                                           \
                dst++;                                                                                                 \
                for (int k = 0; k < j; k++) {                                                                          \
                    if (is_missing)                                                                                    \
                        set_missing;                                                                                   \
                    else if (is_vector_end)                                                                            \
                        set_vector_end;                                                                                \
                    else                                                                                               \
                        *dst = *src;                                                                                   \
                    src++;                                                                                             \
                    dst++;                                                                                             \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
    }
            switch (type) {
            case BCF_BT_INT8:
                BRANCH(int8_t, *dst = bcf_int32_missing, *dst = bcf_int32_vector_end, *src == bcf_int8_missing,
                       *src == bcf_int8_vector_end, int);
                break;
            case BCF_BT_INT16:
                BRANCH(int16_t, *dst = bcf_int32_missing, *dst = bcf_int32_vector_end, *src == bcf_int16_missing,
                       *src == bcf_int16_vector_end, int);
                break;
            case BCF_BT_INT32:
                BRANCH(int32_t, *dst = bcf_int32_missing, *dst = bcf_int32_vector_end, *src == bcf_int32_missing,
                       *src == bcf_int32_vector_end, int);
                break;
            case BCF_BT_FLOAT:
                BRANCH(float, bcf_float_set_missing(*dst), bcf_float_set_vector_end(*dst), bcf_float_is_missing(*src),
                       bcf_float_is_vector_end(*src), float);
                break;
            case BCF_BT_CHAR:
                error(
                    "[E::%s] Lifting over of Number=G strings (in your case %s/%s) is not supported yet, sorry!\n"
                    "Note that using FORMAT strings is not a good idea in general - it is slow to parse and does not "
                    "compress\n"
                    "well, it is better to use integer codes instead. If you don't really need it, use\n"
                    "`bcftools annotate -x` to remove the annotation before lifting over.\n",
                    __func__, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
                break;
            default:
                fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
                exit(1);
            }
#undef BRANCH
            len = n * n_als * (n_als + 1) / 2;
        } else if (length == BCF_VL_A || length == BCF_VL_R) { // define how to add an element to a VL_A or VL_R array
            if (type != BCF_BT_CHAR) assert(len == n_als - 1 - (length == BCF_VL_A));
#define BRANCH(type_t, set_missing, set_vector_end, is_missing, is_vector_end, out_type_t)                             \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        out_type_t *dst = (out_type_t *)(*tmp_arr);                                                                    \
        for (int i = 0; i < n; i++) {                                                                                  \
            set_missing;                                                                                               \
            dst++;                                                                                                     \
            for (int j = 0; j < len; j++) {                                                                            \
                if (is_missing)                                                                                        \
                    set_missing;                                                                                       \
                else if (is_vector_end)                                                                                \
                    set_vector_end;                                                                                    \
                else                                                                                                   \
                    *dst = *src;                                                                                       \
                src++;                                                                                                 \
                dst++;                                                                                                 \
            }                                                                                                          \
        }                                                                                                              \
    }
            switch (type) {
            case BCF_BT_INT8:
                BRANCH(int8_t, *dst = bcf_int32_missing, *dst = bcf_int32_vector_end, *src == bcf_int8_missing,
                       *src == bcf_int8_vector_end, int);
                break;
            case BCF_BT_INT16:
                BRANCH(int16_t, *dst = bcf_int32_missing, *dst = bcf_int32_vector_end, *src == bcf_int16_missing,
                       *src == bcf_int16_vector_end, int);
                break;
            case BCF_BT_INT32:
                BRANCH(int32_t, *dst = bcf_int32_missing, *dst = bcf_int32_vector_end, *src == bcf_int32_missing,
                       *src == bcf_int32_vector_end, int);
                break;
            case BCF_BT_FLOAT:
                BRANCH(float, bcf_float_set_missing(*dst), bcf_float_set_vector_end(*dst), bcf_float_is_missing(*src),
                       bcf_float_is_vector_end(*src), float);
                break;
            case BCF_BT_CHAR:
                s = (char *)(*tmp_arr);
                for (int i = 0; i < n; i++) {
                    *s++ = '.';
                    *s++ = ',';
                    memcpy(s, p, len);
                    s += len;
                    p += len;
                }
                if (n == 1) *s = '\0'; // required by bcf_update_info() as it runs strlen()
                break;
            default:
                fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
                exit(1);
            }
#undef BRANCH
            len = n * (type == BCF_BT_CHAR ? len + 2 : n_als - (length == BCF_VL_A));
        } else {
            error("Error: only possible to handle BCF_VL_G and BCF_VL_R fields\n");
        }

        if (coltype == BCF_HL_INFO) // INFO fields
            bcf_update_info(hdr, rec, key, *tmp_arr, len, field_type);
        else // FORMAT fields
            bcf_update_format(hdr, rec, key, *tmp_arr, len, field_type);
    } else {                      // swap reference allele entry
        if (length == BCF_VL_G) { // define how to swap two elements in a VL_G array
            if (type != BCF_BT_CHAR) assert(len == n_als * (n_als + 1) / 2);
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (int i = 0; i < n; i++) {                                                                                  \
            type_t tmp = ptr[0];                                                                                       \
            ptr[0] = ptr[swap * (swap + 3) / 2];                                                                       \
            ptr[swap * (swap + 3) / 2] = tmp;                                                                          \
            for (int j = 1; j < swap; j++) {                                                                           \
                tmp = ptr[j * (j + 1) / 2];                                                                            \
                ptr[j * (j + 1) / 2] = ptr[swap * (swap + 1) / 2 + j];                                                 \
                ptr[swap * (swap + 1) / 2 + j] = tmp;                                                                  \
            }                                                                                                          \
            for (int j = swap + 1; j < n_als; j++) {                                                                   \
                tmp = ptr[j * (j + 1) / 2 + swap];                                                                     \
                ptr[j * (j + 1) / 2 + swap] = ptr[j * (j + 1) / 2];                                                    \
                ptr[j * (j + 1) / 2] = tmp;                                                                            \
            }                                                                                                          \
            ptr += n_als * (n_als + 1) / 2;                                                                            \
        }                                                                                                              \
    }
            switch (type) {
            case BCF_BT_INT8:
                BRANCH(int8_t);
                break;
            case BCF_BT_INT16:
                BRANCH(int16_t);
                break;
            case BCF_BT_INT32:
                BRANCH(int32_t);
                break;
            case BCF_BT_FLOAT:
                BRANCH(float);
                break;
            case BCF_BT_CHAR:
                error(
                    "[E::%s] Lifting over of Number=G strings (in your case %s/%s) is not supported yet, sorry!\n"
                    "Note that using FORMAT strings is not a good idea in general - it is slow to parse and does not "
                    "compress\n"
                    "well, it is better to use integer codes instead. If you don't really need it, use\n"
                    "`bcftools annotate -x` to remove the annotation before lifting over.\n",
                    __func__, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
                break;
            default:
                fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
                exit(1);
            }
#undef BRANCH
        } else if (length == BCF_VL_R) { // define how to swap two elements in a VL_R array
            if (type != BCF_BT_CHAR) assert(len == n_als);
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (int i = 0; i < n; i++) {                                                                                  \
            type_t tmp = ptr[0];                                                                                       \
            ptr[0] = ptr[swap];                                                                                        \
            ptr[swap] = tmp;                                                                                           \
            ptr += len;                                                                                                \
        }                                                                                                              \
    }
            switch (type) {
            case BCF_BT_INT8:
                BRANCH(int8_t);
                break;
            case BCF_BT_INT16:
                BRANCH(int16_t);
                break;
            case BCF_BT_INT32:
                BRANCH(int32_t);
                break;
            case BCF_BT_FLOAT:
                BRANCH(float);
                break;
            case BCF_BT_CHAR:
                s = (char *)(*tmp_arr);
                for (int i = 0; i < n; i++) {
                    memcpy(s, p, len);
                    char *ptr1 = strchr(s, ',');
                    if (!ptr1) error("Error: string %s contains less then %d elements", s, n_als);
                    char *ptr2 = ptr1;
                    for (int j = 1; ptr2 && j < swap; j++) ptr2 = strchr(ptr2 + 1, ',');
                    if (!ptr2) error("Error: string %s contains less then %d elements", s, n_als);
                    char *ptr3 = swap == n_als - 1 ? s + len : strchr(ptr2 + 1, ',');
                    if (!ptr3) error("Error: string %s contains less then %d elements", s, n_als);
                    memcpy(p, ptr2 + 1, ptr3 - ptr2 - 1);
                    p += ptr3 - ptr2 - 1;
                    memcpy(p, ptr1, ptr2 - ptr1 + 1);
                    p += ptr2 - ptr1 + 1;
                    memcpy(p, s, ptr1 - s);
                    p += len + ptr1 - ptr3;
                    s += len;
                }
                break;
            default:
                fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
                exit(1);
            }
#undef BRANCH
        }
    }
}

// this function updates INFO/FORMAT Number=A Type=Float records if the number or order of alleles has changed
static void reverse_A_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, float sum, int swap,
                             void **tmp_arr, int *m_tmp_arr) {
    int n, type, len;
    float *p;
    int n_als = rec->n_allele;
    int length = bcf_hdr_id2length(hdr, coltype, int_id);   // it should be BCF_VL_A
    int field_type = bcf_hdr_id2type(hdr, coltype, int_id); // it should be BCF_HT_REAL
    const char *key = bcf_hdr_int2id(hdr, BCF_DT_ID, int_id);

    if (coltype == BCF_HL_INFO) { // INFO fields
        n = 1;
        bcf_info_t *info = bcf_get_info_id(rec, int_id);
        if (!info) return;
        len = info->len;
        type = info->type;
        p = (float *)info->vptr;
    } else { // FORMAT fields
        n = rec->n_sample;
        bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
        if (!n || !fmt) return;
        len = fmt->n;
        type = fmt->type;
        p = (float *)fmt->p;
    }
    assert(length == BCF_VL_A && field_type == BCF_HT_REAL && type == BCF_BT_FLOAT);

    int size = sizeof(float) * n * (n_als - 1);
    hts_expand(char, size, *m_tmp_arr, *tmp_arr);

    if (swap < 0) {
        assert(len == n_als - 2);
        float *dst = (float *)(*tmp_arr);
        for (int i = 0; i < n; i++) {
            dst[0] = sum;
            for (int j = 0; j < len; j++) {
                dst[0] -= p[j];
                dst[j + 1] = p[j];
            }
            p += len;
            dst += n_als - 1;
        }
        len = n * (n_als - 1);

        if (coltype == BCF_HL_INFO) // INFO fields
            bcf_update_info_float(hdr, rec, key, *tmp_arr, len);
        else // FORMAT fields
            bcf_update_format_float(hdr, rec, key, *tmp_arr, len);
    } else {
        assert(len == n_als - 1);
        for (int i = 0; i < n; i++) {
            float ref = sum;
            for (int j = 0; j < len; j++) ref -= p[j];
            p[swap - 1] = ref;
            p += len;
        }
        if (coltype == BCF_HL_INFO && len == 1) {
            bcf_info_t *info = bcf_get_info_id(rec, int_id);
            info->v1.f = ((float *)info->vptr)[0];
        }
    }
}

// this function updates INFO/FORMAT records if the number or order of alleles has changed
static void flip_A_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int swap) {
    int n, type, len;
    float *p;
    int n_als = rec->n_allele;
    int length = bcf_hdr_id2length(hdr, coltype, int_id);
    char *s;

    if (coltype == BCF_HL_INFO) { // INFO fields
        n = 1;
        bcf_info_t *info = bcf_get_info_id(rec, int_id);
        if (!info) return;
        len = info->len;
        type = info->type;
        p = (float *)info->vptr;
    } else { // FORMAT fields
        n = rec->n_sample;
        bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
        if (!n || !fmt) return;
        len = fmt->n;
        type = fmt->type;
        p = (float *)fmt->p;
    }
    assert(swap > 0 && length == BCF_VL_A);

    if (type != BCF_BT_CHAR) assert(len == n_als - 1);
#define BRANCH(type_t, easy_access)                                                                                    \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (int i = 0; i < n; i++) {                                                                                  \
            ptr[swap - 1] = -ptr[swap - 1];                                                                            \
            ptr += len;                                                                                                \
        }                                                                                                              \
        if (coltype == BCF_HL_INFO && len == 1) {                                                                      \
            bcf_info_t *info = bcf_get_info_id(rec, int_id);                                                           \
            easy_access = ((type_t *)info->vptr)[0];                                                                   \
        }                                                                                                              \
    }
    switch (type) {
    case BCF_BT_INT8:
        BRANCH(int8_t, info->v1.i);
        break;
    case BCF_BT_INT16:
        BRANCH(int16_t, info->v1.i);
        break;
    case BCF_BT_INT32:
        BRANCH(int32_t, info->v1.i);
        break;
    case BCF_BT_FLOAT:
        BRANCH(float, info->v1.f);
        break;
    case BCF_BT_CHAR:
        s = (char *)p;
        for (int i = 0; i < n; i++) {
            char *ptr = s;
            if (swap > 1) {
                ptr = strchr(ptr, ',');
                for (int j = 0; ptr && j < swap - 2; j++) ptr = strchr(ptr + 1, ',');
                ptr++;
            }
            while (ptr < s + len && *ptr != ',') {
                if (*ptr == '+')
                    *ptr = '-';
                else if (*ptr == '-')
                    *ptr = '+';
                ptr++;
            }
            s += len;
        }
        break;
    default:
        fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
        exit(1);
    }
}

bcf1_t *process(bcf1_t *rec) {
    int is_snp = bcf_is_snp(rec);
    int is_symbolic = bcf_is_symbolic(rec);
    if (is_symbolic && !args->warning_symbolic) {
        fprintf(stderr, "Warning: input VCF includes symbolic alleles that might not properly lift over\n");
        args->warning_symbolic = 1;
    }

    const char *src_chr = bcf_hdr_id2name(args->in_hdr, rec->rid);
    hts_pos_t src_pos = rec->pos + 1;
    if (args->reject_fh || args->write_src) {
        args->tmp_kstr.l = 0;
        for (int i = 0; i < rec->n_allele; i++) {
            char *allele = rec->d.allele[i];
            int len = strlen(allele);
            kputsn_(allele, len, &args->tmp_kstr);
            kputc_(',', &args->tmp_kstr);
        }
        args->tmp_kstr.l--;
        args->tmp_kstr.s[args->tmp_kstr.l] = '\0';
    }

    if (!is_snp && !is_symbolic) {
        if (args->src_fai) {
            if (is_extension_needed(rec, &args->tmp_arr, &args->m_tmp_arr)) {
                hts_expand0(kstring_t, rec->n_allele, args->m_tmp_als, args->tmp_als);
                bcf1_realign_t this = {args->in_hdr, rec, args->tmp_als, args->src_fai, args->aln_win, NULL, 0, 0};
                extend_alleles(&this);
            }
        } else {
            if (!args->warning_indel) {
                fprintf(stderr,
                        "Warning: input VCF includes indels but option --src-fasta-ref is missing which is not "
                        "recommended\n");
                args->warning_indel = 1;
            }
        }
    }

    args->ntotal++;
    if (!args->idx) return rec;

    if (!args->lift_mt && rec->rid == args->in_mt_rid) {
        rec->rid = args->out_mt_rid;
        return rec;
    }

    // lift over record coordinates
    int ret, rid, strand, npad = 0;
    hts_pos_t pos;
    int len = strlen(rec->d.allele[0]);
    if (is_snp || is_symbolic)
        ret = liftover_snp(args->idx, args->itr, src_chr, rec->pos + 1, 0, &rid, &pos, &strand);
    else
        ret = liftover_indel(args->idx, args->itr, src_chr, rec->pos + 1, &len, args->indel_win, &rid, &pos, &strand,
                             &npad);
    if (ret < 0 || (npad && !args->src_fai)) {
        args->nrejected++;
        if (args->reject_fh) {
            // restore original position and alleles
            if (!is_snp && !is_symbolic) {
                rec->pos = src_pos - 1;
                bcf_update_alleles_str(args->in_hdr, rec, args->tmp_kstr.s);
            }
            if (bcf_write(args->reject_fh, args->in_hdr, rec) < 0) error("Error: Unable to write to output VCF file\n");
        }
        return NULL;
    }

    // pad the sequence to reach coordinates shared by both references (this happens rarely)
    if (npad && args->src_fai) {
        hts_expand0(kstring_t, rec->n_allele, args->m_tmp_als, args->tmp_als);
        bcf1_realign_t this = {args->in_hdr, rec, args->tmp_als, args->src_fai, args->aln_win, NULL, 0, 0};
        initialize_alleles(&this);
        if (npad < 0)
            pad_left(&this, -npad);
        else
            pad_right(&this, npad);
        update_alleles(&this);
    }

    if (args->write_src) {
        bcf_update_info_string(args->out_hdr, rec, "SRC_CHROM", src_chr);
        bcf_update_info_int32(args->out_hdr, rec, "SRC_POS", &src_pos, 1);
        bcf_update_info_string(args->out_hdr, rec, "SRC_REF_ALT", args->tmp_kstr.s);
    }

    // flip alleles
    if (strand) {
        bcf_update_info_flag(args->out_hdr, rec, args->tags.flip_tag, NULL, 1);
        flip_alleles(rec, args->out_hdr, &args->tmp_kstr);
    }

    rec->rid = rid;
    rec->pos = pos;
    if (is_symbolic) return rec;

    int swap = find_reference(rec, args->out_hdr, args->dst_fai, len, &args->tmp_kstr);

    // left align indels
    if (!is_snp && !args->no_left_align && !is_left_aligned(rec, &args->tmp_arr, &args->m_tmp_arr)) {
        hts_expand0(kstring_t, rec->n_allele, args->m_tmp_als, args->tmp_als);
        bcf1_realign_t this = {args->out_hdr, rec, args->tmp_als, args->dst_fai, args->aln_win, NULL, 0, 0};
        initialize_alleles(&this);
        trim_right(&this);
        trim_left(&this);
        shift_left(&this);
        update_alleles(&this);
    }

    // address records that have the reference allele added or swapped
    if (swap) {
        if (swap < 0)
            args->nref_added++;
        else
            args->nswapped++;
        bcf_update_info_int32(args->out_hdr, rec, args->tags.swap_tag, &swap, 1);
        tags_t *tags = &args->tags;
        // address tags to be dropped
        for (int i = 0; i < tags->n_drop_tags; i++) {
            int int_id = tags->drop_tags_int_id[i];
            int coltype = tags->drop_tags_coltype[i];
            int field_type = bcf_hdr_id2type(args->out_hdr, coltype, int_id);
            const char *key = bcf_hdr_int2id(args->out_hdr, BCF_DT_ID, int_id);
            if (coltype == BCF_HL_INFO)
                bcf_update_info(args->out_hdr, rec, key, NULL, 0, field_type);
            else
                bcf_update_format(args->out_hdr, rec, key, NULL, 0, field_type);
        }
        // address genotypes
        update_genotypes(args->out_hdr, rec, tags->gt_id, swap, &args->tmp_arr, &args->m_tmp_arr);
        // address GR tags
        for (int i = 0; i < tags->n_GR_tags; i++) {
            int int_id = tags->GR_tags_int_id[i];
            int coltype = tags->GR_tags_coltype[i];
            update_AGR_record(args->out_hdr, rec, int_id, coltype, swap, &args->tmp_arr, &args->m_tmp_arr);
        }
        // address A tags to be reversed
        for (int i = 0; i < tags->n_reverse_tags; i++) {
            int int_id = tags->reverse_tags_int_id[i];
            int coltype = tags->reverse_tags_coltype[i];
            float sum = tags->reverse_tags_sum[i];
            reverse_A_record(args->out_hdr, rec, int_id, coltype, sum, swap, &args->tmp_arr, &args->m_tmp_arr);
        }
        // address A tags to be flipped
        for (int i = 0; i < tags->n_flip_tags; i++) {
            int int_id = tags->flip_tags_int_id[i];
            int coltype = tags->flip_tags_coltype[i];
            if (swap < 0)
                update_AGR_record(args->out_hdr, rec, int_id, coltype, swap, &args->tmp_arr, &args->m_tmp_arr);
            else
                flip_A_record(args->out_hdr, rec, int_id, coltype, swap);
        }
        // address tags with allele indexes
        for (int i = 0; i < tags->n_genotype_tags; i++) {
            int int_id = tags->genotype_tags_int_id[i];
            int coltype = tags->genotype_tags_coltype[i];
            update_genotype_record(args->out_hdr, rec, int_id, coltype, swap, &args->tmp_arr, &args->m_tmp_arr);
        }
    }

    return rec;
}

void destroy(void) {
    fprintf(stderr, "Lines   total/swapped/reference added/rejected:\t%d/%d/%d/%d\n", args->ntotal, args->nswapped,
            args->nref_added, args->nrejected);
    free(args->tmp_kstr.s);
    for (int i = 0; i < args->m_tmp_als; i++) free(args->tmp_als[i].s);
    free(args->tmp_als);
    free(args->tmp_arr);
    if (args->reject_fh && hts_close(args->reject_fh) < 0) error("Close failed: %s\n", args->reject_fh->fn);
    tags_destroy(&args->tags);
    if (args->idx) regidx_destroy(args->idx);
    if (args->itr) regitr_destroy(args->itr);
    free(args->chains);
    free(args->blocks);
    fai_destroy(args->src_fai);
    fai_destroy(args->dst_fai);
    free(args);
}
