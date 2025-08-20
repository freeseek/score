/* The MIT License

   Copyright (C) 2022-2025 Giulio Genovese

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
#include <htslib/khash.h> // required to reset the contigs dictionary and table
#include "bcftools.h"
#include "regidx.h" // cannot use htslib/regdix.h see http://github.com/samtools/htslib/pull/761
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)

#define LIFTOVER_VERSION "2025-08-20"

#define FLIP_TAG "FLIP"
#define SWAP_TAG "SWAP"
#define DROP_TAGS "."
#define AC_TAGS "INFO/AC,FMT/AC"
#define AF_TAGS "INFO/AF,FMT/AF,FMT/AP1,FMT/AP2"
#define DS_TAGS "FMT/DS"
#define GT_TAGS "INFO/ALLELE_A,INFO/ALLELE_B"
#define ES_TAGS "FMT/EZ,FMT/ES,FMT/ED"

#define NO_RULE 0
#define RULE_DROP 1
#define RULE_AC 2
#define RULE_AF 3
#define RULE_DS 4
#define RULE_GT 5
#define RULE_ES 6
#define RULE_AGR 7
const char *rules_str[] = {"DROP", "AC", "AF", "DS", "GT", "FLIP", "AGR"};

static inline char rev_nt(char iupac) {
    static const char iupac_complement[128] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F,
        0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, '-',  0x2E, '/',
        0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E, 0x3F,
        0x40, 'T',  'V',  'G',  'H',  0x45, 0x46, 'C',  'D',  0x49, 0x4A, 'M',  0x4C, 'K',  'N',  0x4F,
        0x50, 0x51, 'Y',  'S',  'A',  0x55, 'B',  'W',  0x58, 'R',  0x5A, ']',  0x5C, '[',  0x5E, 0x5F,
        0x60, 't',  'v',  'g',  'h',  0x65, 0x66, 'c',  'd',  0x69, 0x6A, 'm',  0x6C, 'k',  'n',  0x6F,
        0x70, 0x71, 'y',  's',  'a',  0x75, 'b',  'w',  0x78, 'r',  0x7A, 0x7B, 0x7C, 0x7D, 0x7E, 0x7F,
    };
    return iupac_complement[(int)(iupac & 0x7F)];
}

static inline void reverse_complement(char *str) {
    int i, len = strlen(str);
    for (i = 0; i < len / 2; i++) {
        char tmp = str[i];
        str[i] = rev_nt(str[len - i - 1]);
        str[len - i - 1] = rev_nt(tmp);
    }
    if (len % 2 == 1) str[len / 2] = rev_nt(str[len / 2]);
}

typedef struct {
    int tStart;
    int qStart;
    int size;
    int chain_ind;  // index of chain the block belongs to
    int tStart_gap; // size of gap to reach previous contiguous block in target space
    int tEnd_gap;   // size of gap to reach next contiguous block in target space
    int qStart_gap; // size of gap to reach previous contiguous block in query space
    int qEnd_gap;   // size of gap to reach next contiguous block in query space
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
    int block_ind; // index of first block of the chain
    int n_blocks;  // number of blocks in the chain
} chain_t;

// return previous block from the same chain (NULL if it is the first block)
static inline const block_t *prev_block(const block_t *block, const block_t *blocks, const chain_t *chains) {
    int ind = (block - blocks) - chains[block->chain_ind].block_ind;
    return ind == 0 ? NULL : block - 1;
}

// return next block from the same chain (NULL if it is the last block)
static inline const block_t *next_block(const block_t *block, const block_t *blocks, const chain_t *chains) {
    int ind = (block - blocks) - chains[block->chain_ind].block_ind;
    return ind == chains[block->chain_ind].n_blocks - 1 ? NULL : block + 1;
}

typedef struct {
    int int_id;  // VCF header int_id
    int coltype; // whether BCF_HL_INFO or BCF_HL_FMT
    int rule;    // which rule should apply (DROP, AC, AF, DS, GT, FLIP, or AGR)
} tag_t;

typedef struct {
    bcf_hdr_t *in_hdr;
    bcf_hdr_t *out_hdr;
    faidx_t *src_fai;
    faidx_t *dst_fai;
    int n_ctgs;
    int n_chains;
    chain_t *chains;
    block_t *blocks;
    regidx_t *idx;
    regitr_t *itr;
    htsFile *reject_fh;

    int max_indel_inc;
    int aln_win;
    int in_mt_rid;
    int out_mt_rid;
    int lift_mt;
    int no_left_align;
    int write_src;
    int write_fail;
    int write_nw;
    int reject_filter;
    int lift_end;

    int info_end_id;
    int info_an_id;
    int fmt_gt_id;
    int fmt_an_id;
    int *af_arr;
    int *ploidy_arr;
    const char *flip_tag;
    const char *swap_tag;
    int n_tags;
    tag_t *tags;

    int warning_symbolic;
    int warning_indel;
    int ntotal;
    int nswapped;
    int nref_added;
    int nrejected;

    kstring_t tmp_kstr;
    kstring_t tmp_pad;
    kstring_t *tmp_als;
    int m_tmp_als;
    int8_t *tmp_arr;
    int m_tmp_arr;
    int32_t *int32_arr;
    int m_int32_arr;
} args_t;

args_t *args;

/****************************************
 * CHAIN FUNCTIONS                      *
 ****************************************/

static inline int bcf_hdr_name2id_flexible(const bcf_hdr_t *hdr, char *chr) {
    if (!chr || strlen(chr) < 1) return -1;
    int rid = bcf_hdr_name2id(hdr, chr);
    if (rid >= 0) return rid;
    if (strncmp(chr, "chr", 3) == 0) rid = bcf_hdr_name2id(hdr, chr + 3);
    if (rid >= 0) return rid;
    char buf[] = {'c', 'h', 'r', chr[0], chr[1], '\0'};
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
                       chain_t **chains, block_t **blocks) {
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
        int ncols = ksplit_core(str.s, 0, &moff, &off);
        if (ncols != 13) error("Wrong number of columns in the chain file: %s\n", fp->fn);

        // read Header Lines
        if (strcmp(&str.s[off[0]], "chain") != 0)
            error("Chain line should start with word \"chain\" but \"%s\" found in the chain file: %s\n",
                  &str.s[off[0]], fp->fn);
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
        if (str.s[off[4]] != '+')
            error("Chain line fifth column should be \"+\" but \"%s\" found in the chain file: %s\n", &str.s[off[4]],
                  fp->fn);
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
        if (str.s[off[9]] != '+' && str.s[off[9]] != '-')
            error("Chain line tenth column should be \"+\" or \"-\" but \"%s\" found in the chain file: %s\n",
                  &str.s[off[9]], fp->fn);
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
            int ncols = ksplit_core(str.s, 0, &moff, &off);
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
                block->tStart_gap = dt;
                block->tEnd_gap = -1;
                block->qStart_gap = dq;
                block->qEnd_gap = -1;
                chain->n_blocks++;
            }

            if (ncols == 1) {
                tStart += size;
                if (chain->tStart + tStart != chain->tEnd)
                    error("Chain malformed as target interval %s:%d-%d not fully covered in the chain file: %s\n",
                          bcf_hdr_id2name(in_hdr, chain->t_rid), chain->tStart, chain->tEnd, fp->fn);
                qStart += size;
                if (chain->qStart + qStart != chain->qEnd)
                    error("Chain malformed as query interval %s:%d-%d not fully covered in the chain file: %s\n",
                          bcf_hdr_id2name(in_hdr, chain->q_rid), chain->qStart, chain->qEnd, fp->fn);
            } else {
                dt = strtol(&str.s[off[1]], &tmp, 0);
                if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[1]], fp->fn);
                dq = strtol(&str.s[off[2]], &tmp, 0);
                if (*tmp) error("Could not parse integer %s in the chain file: %s\n", &str.s[off[2]], fp->fn);
                tStart += size + dt;
                qStart += size + dq;
                merge_in_progress = dt == dq && dt <= max_snp_gap;
                if (merge_in_progress) {
                    block->size += dt;
                } else {
                    block->tEnd_gap = dt;
                    block->qEnd_gap = dq;
                }
            }
        }
    }

    free(off);
    free(str.s);
    return n_chains;
}

static void write_chains(FILE *stream, const bcf_hdr_t *in_hdr, const bcf_hdr_t *out_hdr, const chain_t *chains,
                         int n_chains, block_t *blocks) {
    int i, j;
    for (i = 0; i < n_chains; i++) {
        const chain_t *chain = &chains[i];
        const char *tName = bcf_hdr_id2name(in_hdr, chain->t_rid);
        const char *qName = bcf_hdr_id2name(out_hdr, chain->q_rid);
        fprintf(stream, "chain %" PRIhts_pos " %s %d + %d %d %s %d %c %d %d %d\n", chain->score, tName, chain->tSize,
                chain->tStart, chain->tEnd, qName, chain->qSize, chain->qStrand ? '-' : '+', chain->qStart, chain->qEnd,
                chain->id);
        for (j = 0; j < chain->n_blocks; j++) {
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

KHASH_MAP_INIT_INT(32, char)
static regidx_t *regidx_init_chains(const bcf_hdr_t *in_hdr, const chain_t *chains, int n_chains, block_t *blocks) {
    khash_t(32) *h = kh_init(32);
    regidx_t *idx = regidx_init(NULL, NULL, NULL, sizeof(uint64_t), NULL);
    int i, j;
    for (i = 0; i < n_chains; i++) {
        const chain_t *chain = &chains[i];
        if (chain->t_rid < 0 || chain->q_rid < 0) continue;

        // check whether the chain has already been added as Ensembl duplicates chains in the chain file
        if (chain->id) {
            khiter_t k = kh_get(32, h, chain->id);
            int ret, is_missing = (k == kh_end(h));
            if (!is_missing) continue;
            kh_put(32, h, chain->id, &ret);
        }

        const char *name = bcf_hdr_id2name(in_hdr, chain->t_rid);
        int len = strlen(name);
        for (j = 0; j < chain->n_blocks; j++) {
            int block_ind = chain->block_ind + j;
            const block_t *block = &blocks[block_ind];
            regidx_push(idx, (char *)name, (char *)name + len, chain->tStart + block->tStart + 1,
                        chain->tStart + block->tStart + block->size, (void *)&block_ind);
        }
    }
    kh_destroy(32, h);
    return idx;
}

/****************************************
 * TAGS FUNCTIONS                       *
 ****************************************/

// assign AGR tags rules
static void find_AGR_tags(const bcf_hdr_t *hdr, int *info_rules, int *fmt_rules) {
    const int coltypes[2] = {BCF_HL_INFO, BCF_HL_FMT};
    int i, j;
    for (i = 0; i < 2; i++) {
        int *rules = i == 0 ? info_rules : fmt_rules;
        for (j = 0; j < hdr->n[BCF_DT_ID]; j++) {
            // http://github.com/samtools/htslib/issues/1538
            if (!hdr->id[BCF_DT_ID][j].val || bcf_hdr_id2coltype(hdr, coltypes[i], j) == 0xf) continue;
            int length = bcf_hdr_id2length(hdr, coltypes[i], j);
            if (length == BCF_VL_A || length == BCF_VL_G || length == BCF_VL_R) rules[j] = RULE_AGR;
        }
    }
}

// assign special tags rules
static void assign_tags(const bcf_hdr_t *hdr, const char *tags, int length_flag, int field_type_flag, int rule,
                        int *info_rules, int *fmt_rules) {
    if (strcmp(tags, ".") == 0) return;
    char *s = strdup(tags);
    int i, moff = 0, *off = NULL;
    int n = ksplit_core(s, ',', &moff, &off);
    for (i = 0; i < n; i++) {
        char *ss = &s[off[i]];
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
        int int_id = bcf_hdr_id2int(hdr, BCF_DT_ID, ss);
        if (coltype == -1) {
            if (bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, int_id)) error("Error: did you mean INFO/%s?\n", ss);
            if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, int_id)) error("Error: did you mean FORMAT/%s?\n", ss);
            error("No matching tag  %s\n", ss);
        }
        if (!bcf_hdr_idinfo_exists(hdr, coltype, int_id)) continue;
        unsigned int length = bcf_hdr_id2length(hdr, coltype, int_id);
        if (!(1 << length & length_flag)) {
            if (length) {
                char bcf_vl = length == BCF_VL_A ? 'A' : length == BCF_VL_G ? 'G' : length == BCF_VL_R ? 'R' : '.';
                error("The %s tag \"%s\" is a Number=%c tag which is not the expected length for this tag\n",
                      coltype == BCF_HL_INFO ? "INFO" : "FORMAT", ss, bcf_vl);
            } else {
                int number = bcf_hdr_id2number(hdr, coltype, int_id);
                error("The %s tag \"%s\" is a Number=%d tag which is not the expected length for this tag\n",
                      coltype == BCF_HL_INFO ? "INFO" : "FORMAT", ss, number);
            }
        }
        int field_type = bcf_hdr_id2type(hdr, coltype, int_id);
        if (!(1 << field_type & field_type_flag)) {
            const char *bcf_ht = field_type == BCF_HT_FLAG   ? "Flag"
                                 : field_type == BCF_HT_INT  ? "Integer"
                                 : field_type == BCF_HT_REAL ? "Float"
                                 : field_type == BCF_HT_STR  ? "String"
                                                             : ".";
            error("The %s tag \"%s\" is a Type=%s tag which is not the expected type for this tag\n",
                  coltype == BCF_HL_INFO ? "INFO" : "FORMAT", ss, bcf_ht);
        }
        int *rules = coltype == BCF_HL_INFO ? info_rules : fmt_rules;
        rules[int_id] = rule;
    }
    free(off);
    free(s);
}

// make list of tags rules
static int compress_tags(int *info_rules, int *fmt_rules, int n, tag_t **tags) {
    const int coltypes[2] = {BCF_HL_INFO, BCF_HL_FMT};
    int i, j, n_tags = 0, m_tags = 0;
    *tags = NULL;
    for (i = 0; i < 2; i++) {
        int *rules = i == 0 ? info_rules : fmt_rules;
        for (j = 0; j < n; j++) {
            if (!rules[j]) continue;
            hts_expand(tag_t, n_tags + 1, m_tags, *tags);
            tag_t *tag = &(*tags)[n_tags++];
            tag->int_id = j;
            tag->coltype = coltypes[i];
            tag->rule = rules[j];
        }
    }
    return n_tags;
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Lift over a VCF from one genome build to another.\n"; }

const char *usage(void) {
    return "\n"
           "About: Lift over a VCF from one genome build to another. "
           "(version " LIFTOVER_VERSION
           " http://github.com/freeseek/score)\n"
           "[ Genovese, G., et al. BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants\n"
           "across genome assemblies. Bioinformatics 40, Issue 2 (2024) http://doi.org/10.1093/bioinformatics/btae038 "
           "]\n"
           "\n"
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
           "       --max-indel-inc <int>       maximum distance used to increase the size an indel during liftover "
           "[250]\n"
           "       --lift-mt                   force liftover of MT/chrMT [automatically determined from contig "
           "lengths]\n"
           "       --print-blocks <file>       output contiguous blocks used for the liftOver\n"
           "       --no-left-align             do not attempt to left align indels after liftover\n"
           "       --reject <file>             output variants that cannot be lifted over\n"
           "   -O, --reject-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "       --write-src                 write the source contig/position/alleles for lifted variants\n"
           "       --write-fail                write whether the 5' and 3' anchors have failed to lift\n"
           "       --write-nw                  write the Needleman-Wunsch alignments when required\n"
           "       --write-reject              write the reason variants cannot be lifted over\n"
           "\n"
           "Options for how to update INFO/FORMAT records:\n"
           "       --fix-tags                  fix Number type for INFO/AC, INFO/AF, FORMAT/GP, and FORMAT/DS tags\n"
           "       --lift-end                  lift the position of the INFO/END tag instead of recomputing it\n"
           "       --flip-tag <string>         INFO annotation flag to record whether alleles are flipped [FLIP]\n"
           "       --swap-tag <string>         INFO annotation to record when alleles are swapped [SWAP]\n"
           "       --drop-tags <list>          tags to drop when alleles are swapped [" DROP_TAGS
           "]\n"
           "       --ac-tags <list>            AC-like tags (must be Number=A,Type=Integer/Float) [" AC_TAGS
           "]\n"
           "       --af-tags <list>            AF-like tags (must be Number=A,Type=Float) [" AF_TAGS
           "]\n"
           "       --ds-tags <list>            DS-like tags (must be Number=A,Type=Float) [" DS_TAGS
           "]\n"
           "       --gt-tags <list>            tags with integers like FORMAT/GT (must be Type=Integer) [" GT_TAGS
           "]\n"
           "       --es-tags <list>            GWAS-VCF tags (must be Number=A) [" ES_TAGS
           "]\n"
           "\n"
           "Examples:\n"
           "      bcftools +liftover -Ou input.hg19.bcf -- -s hg19.fa -f hg38.fa \\\n"
           "        -c hg19ToHg38.over.chain.gz | bcftools sort -Ob -o output.hg38.bcf -W\n"
           "      bcftools +liftover -Ou GRCh38_dbSNPv156.vcf.gz -- -s hg38.fa -f chm13v2.0.fa \\\n"
           "        -c hg38ToHs1.over.chain.gz | bcftools sort -Oz -o chm13v2.0_dbSNPv156.vcf.gz -W=tbi\n"
           "\n"
           "To obtain liftover chain files:\n"
           "      wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz\n"
           "      wget http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz\n"
           "      wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz\n"
           "\n";
}

static inline FILE *get_file_handle(const char *str) {
    FILE *ret;
    if (strcmp(str, "-") == 0)
        ret = stdout;
    else {
        ret = fopen(str, "w");
        if (!ret) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out) {
    args = (args_t *)calloc(1, sizeof(args_t));
    args->in_hdr = in;
    args->out_hdr = out;
    args->max_indel_inc = 250;
    args->aln_win = 100;
    args->flip_tag = FLIP_TAG;
    args->swap_tag = SWAP_TAG;
    int i, j, k;
    int cache_size = 0;
    int output_type = FT_VCF;
    int clevel = -1;
    int max_snp_gap = 1; // maximum distance between two contiguous blocks to allow merging
    int fix_tags = 0;
    char *tmp = NULL;
    const char *src_ref_fname = NULL;
    const char *dst_ref_fname = NULL;
    const char *chain_fname = NULL;
    const char *blocks_fname = NULL;
    const char *reject_fname = NULL;
    char *drop_tags = DROP_TAGS;
    char *ac_tags = AC_TAGS;
    char *af_tags = AF_TAGS;
    char *ds_tags = DS_TAGS;
    char *gt_tags = GT_TAGS;
    char *es_tags = ES_TAGS;

    static struct option loptions[] = {{"src-fasta-ref", required_argument, NULL, 's'},
                                       {"fasta-ref", required_argument, NULL, 'f'},
                                       {"set-cache-size", required_argument, NULL, 1},
                                       {"chain", required_argument, NULL, 'c'},
                                       {"max-snp-gap", required_argument, NULL, 2},
                                       {"max-indel-inc", required_argument, NULL, 3},
                                       {"lift-mt", no_argument, NULL, 4},
                                       {"print-blocks", required_argument, NULL, 5},
                                       {"no-left-align", no_argument, NULL, 6},
                                       {"reject", required_argument, NULL, 7},
                                       {"reject-type", required_argument, NULL, 'O'},
                                       {"write-src", no_argument, NULL, 8},
                                       {"write-fail", no_argument, NULL, 9},
                                       {"write-nw", no_argument, NULL, 10},
                                       {"write-reject", no_argument, NULL, 11},
                                       {"fix-tags", no_argument, NULL, 12},
                                       {"lift-end", no_argument, NULL, 13},
                                       {"flip-tag", required_argument, NULL, 14},
                                       {"swap-tag", required_argument, NULL, 15},
                                       {"drop-tags", required_argument, NULL, 16},
                                       {"ac-tags", required_argument, NULL, 17},
                                       {"af-tags", required_argument, NULL, 18},
                                       {"ds-tags", required_argument, NULL, 19},
                                       {"gt-tags", required_argument, NULL, 20},
                                       {"es-tags", required_argument, NULL, 21},
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
            args->max_indel_inc = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --max-indel-inc %s\n", optarg);
            break;
        case 4:
            args->lift_mt = 1;
            break;
        case 5:
            blocks_fname = optarg;
            break;
        case 6:
            args->no_left_align = 1;
            break;
        case 7:
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
        case 8:
            args->write_src = 1;
            break;
        case 9:
            args->write_fail = 1;
            break;
        case 10:
            args->write_nw = 1;
            break;
        case 11:
            args->reject_filter = 1;
            break;
        case 12:
            fix_tags = 1;
            break;
        case 13:
            args->lift_end = 1;
            break;
        case 14:
            args->flip_tag = optarg;
            break;
        case 15:
            args->swap_tag = optarg;
            break;
        case 16:
            drop_tags = optarg;
            break;
        case 17:
            ac_tags = optarg;
            break;
        case 18:
            af_tags = optarg;
            break;
        case 19:
            ds_tags = optarg;
            break;
        case 20:
            gt_tags = optarg;
            break;
        case 21:
            es_tags = optarg;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage());
            break;
        }
    }

    if (!in || !out) error("Expected input VCF\n%s", usage());
    args->n_ctgs = in->n[BCF_DT_CTG];

    // load target reference file
    if (src_ref_fname) {
        args->src_fai = fai_load(src_ref_fname);
        if (!args->src_fai) error("Could not load the reference %s\n", src_ref_fname);
        if (cache_size) fai_set_cache_size(args->src_fai, cache_size);
    } else {
        if (!dst_ref_fname || !chain_fname)
            error("At least one of --src-fasta-ref or --fasta-ref and --chain options required\n");
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
    for (i = 0; i < out->n[BCF_DT_CTG]; i++) free((void *)out->id[BCF_DT_CTG][i].key);
    out->n[BCF_DT_CTG] = 0;

    int n = faidx_nseq(args->dst_fai);
    for (i = 0; i < n; i++) {
        const char *seq = faidx_iseq(args->dst_fai, i);
        int len = faidx_seq_len(args->dst_fai, seq);
        bcf_hdr_printf(out, "##contig=<ID=%s,length=%d>", seq, len);
    }

    // fix INFO/AC, INFO/AF, FORMAT/GP, and FORMAT/DS tags from the Michigan imputation server
    if (fix_tags) {
        bcf_hdr_t *hdrs[] = {in, out};
        const char *tags[] = {"AC", "AF", "GP", "DS"};
        int coltypes[] = {BCF_HL_INFO, BCF_HL_INFO, BCF_HL_FMT, BCF_HL_FMT};
        int numbers[] = {1, 1, 3, 1};
        int lengths[] = {BCF_VL_A, BCF_VL_A, BCF_VL_G, BCF_VL_A};
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 4; j++) {
                int int_id = bcf_hdr_id2int(hdrs[i], BCF_DT_ID, tags[j]);
                if (int_id < 0) continue;
                bcf_idinfo_t *val = (bcf_idinfo_t *)hdrs[i]->id[BCF_DT_ID][int_id].val;
                if ((val->info[coltypes[j]] & 0xf) == 0xf) continue;
                if ((val->info[coltypes[j]] >> 8 & 0xf) == 0 && val->info[coltypes[j]] >> 12 == numbers[j]) {
                    val->info[coltypes[j]] &= 0xff;
                    val->info[coltypes[j]] += lengths[j] << 8;
                    bcf_hrec_t *hrec = val->hrec[coltypes[j]];
                    for (k = 0; k < hrec->nkeys; k++) {
                        if (!strcmp(hrec->keys[k], "Number")) {
                            if (!strcmp(hrec->vals[k], "1"))
                                hrec->vals[k][0] = 'A';
                            else if (!strcmp(hrec->vals[k], "3"))
                                hrec->vals[k][0] = 'G';
                        }
                    }
                }
            }
        }
    }

    if (bcf_hdr_printf(
            out, "##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Whether alleles flipped strand during liftover\">",
            args->flip_tag)
        < 0)
        error_errno("Failed to add \"%s\" INFO header", args->flip_tag);
    if (bcf_hdr_printf(out,
                       "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Which alternate allele became the reference "
                       "during liftover (-1 for new reference)\">",
                       args->swap_tag)
        < 0)
        error_errno("Failed to add \"%s\" INFO header", args->swap_tag);
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
    if (args->write_fail) {
        if (bcf_hdr_printf(
                out,
                "##INFO=<ID=FAIL5,Number=0,Type=Flag,Description=\"Whether the original 5' anchor failed to lift\">")
            < 0)
            error_errno("Failed to add \"FAIL5\" INFO header");
        if (bcf_hdr_printf(
                out,
                "##INFO=<ID=FAIL3,Number=0,Type=Flag,Description=\"Whether the original 3' anchor failed to lift\">")
            < 0)
            error_errno("Failed to add \"FAIL3\" INFO header");
    }
    if (args->write_nw) {
        if (bcf_hdr_printf(out, "##INFO=<ID=NW,Number=.,Type=String,Description=\"Needleman-Wunsch alignments\">") < 0)
            error_errno("Failed to add \"NW\" INFO header");
    }
    if (src_ref_fname)
        bcf_hdr_printf(out, "##liftover_target_reference=%s",
                       strrchr(src_ref_fname, '/') ? strrchr(src_ref_fname, '/') + 1 : src_ref_fname);
    if (dst_ref_fname)
        bcf_hdr_printf(out, "##liftover_query_reference=%s",
                       strrchr(dst_ref_fname, '/') ? strrchr(dst_ref_fname, '/') + 1 : dst_ref_fname);
    if (chain_fname)
        bcf_hdr_printf(out, "##liftover_chain_file=%s",
                       strrchr(chain_fname, '/') ? strrchr(chain_fname, '/') + 1 : chain_fname);
    if (bcf_hdr_sync(out) < 0) error_errno("Failed to update header");

    // load chain file
    if (max_snp_gap > 20)
        fprintf(
            stderr,
            "Warning: merging contiguous blocks farther than 20 bp apart (--max-snp-gap %d used) is not recommended\n",
            max_snp_gap);
    if (args->max_indel_inc > 250)
        fprintf(stderr,
                "Warning: dealing with indels with edges farther apart than 250 bp apart (--max-indel-inc %d used) is "
                "not recommended\nThe complexity of the Needleman-Wunsch algorithm used to realign indel sequences "
                "near chain gaps grows quadratically in this number",
                args->max_indel_inc);

    htsFile *fp = hts_open(chain_fname, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fp->fn, strerror(errno));
    args->n_chains = read_chains(fp, in, out, max_snp_gap, &args->chains, &args->blocks);
    if (hts_close(fp) < 0) error("Close failed: %s\n", fp->fn);

    if (blocks_fname) {
        FILE *blocks_file = get_file_handle(blocks_fname);
        write_chains(blocks_file, in, out, args->chains, args->n_chains, args->blocks);
        if (blocks_file != stdout && blocks_file != stderr) fclose(blocks_file);
    }

    args->idx = regidx_init_chains(in, args->chains, args->n_chains, args->blocks);
    args->itr = regitr_init(args->idx);

    args->info_end_id = bcf_hdr_id2int(in, BCF_DT_ID, "END");
    if (!bcf_hdr_idinfo_exists(in, BCF_HL_INFO, args->info_end_id)) args->info_end_id = -1;
    args->info_an_id = bcf_hdr_id2int(in, BCF_DT_ID, "AN");
    if (!bcf_hdr_idinfo_exists(in, BCF_HL_INFO, args->info_an_id)) args->info_an_id = -1;
    args->fmt_gt_id = bcf_hdr_id2int(in, BCF_DT_ID, "GT");
    if (!bcf_hdr_idinfo_exists(in, BCF_HL_FMT, args->fmt_gt_id)) args->fmt_gt_id = -1;
    args->fmt_an_id = bcf_hdr_id2int(in, BCF_DT_ID, "AN");
    if (!bcf_hdr_idinfo_exists(in, BCF_HL_FMT, args->fmt_an_id)) args->fmt_an_id = -1;
    args->af_arr = (int *)malloc(sizeof(int) * bcf_hdr_nsamples(in));
    for (i = 0; i < bcf_hdr_nsamples(in); i++) args->af_arr[i] = 1;
    args->ploidy_arr = (int *)malloc(sizeof(int) * bcf_hdr_nsamples(in));

    int *info_rules = (int *)calloc(sizeof(int), in->n[BCF_DT_ID]);
    int *fmt_rules = (int *)calloc(sizeof(int), in->n[BCF_DT_ID]);
    find_AGR_tags(in, info_rules, fmt_rules);
    assign_tags(in, ac_tags, 1 << BCF_VL_A, 1 << BCF_HT_INT | 1 << BCF_HT_REAL, RULE_AC, info_rules, fmt_rules);
    assign_tags(in, af_tags, 1 << BCF_VL_A, 1 << BCF_HT_REAL, RULE_AF, info_rules, fmt_rules);
    assign_tags(in, ds_tags, 1 << BCF_VL_A, 1 << BCF_HT_REAL, RULE_DS, info_rules, fmt_rules);
    assign_tags(in, gt_tags, 1 << BCF_VL_FIXED | 1 << BCF_VL_VAR | 1 << BCF_VL_A | 1 << BCF_VL_G | 1 << BCF_VL_R,
                1 << BCF_HT_INT, RULE_GT, info_rules, fmt_rules);
    assign_tags(in, es_tags, 1 << BCF_VL_A, 1 << BCF_HT_INT | 1 << BCF_HT_REAL | 1 << BCF_HT_STR, RULE_ES, info_rules,
                fmt_rules);
    assign_tags(in, drop_tags, 1 << BCF_VL_FIXED | 1 << BCF_VL_VAR | 1 << BCF_VL_A | 1 << BCF_VL_G | 1 << BCF_VL_R,
                1 << BCF_HT_FLAG | 1 << BCF_HT_INT | 1 << BCF_HT_REAL | 1 << BCF_HT_STR, RULE_DROP, info_rules,
                fmt_rules);
    args->n_tags = compress_tags(info_rules, fmt_rules, in->n[BCF_DT_ID], &args->tags);
    free(info_rules);
    free(fmt_rules);
    for (i = 0; i < args->n_tags; i++) {
        tag_t *tag = &args->tags[i];
        const char *key = bcf_hdr_int2id(in, BCF_DT_ID, tag->int_id);
        fprintf(stderr, "%s/%s is handled by %s rule\n", tag->coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key,
                rules_str[tag->rule - 1]);
    }

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
        if (args->reject_filter) {
            if (bcf_hdr_printf(args->in_hdr,
                               "##FILTER=<ID=MissingContig,Description=\"Contig not defined in the header\">")
                < 0)
                error_errno("Failed to add \"MissingContig\" FILTER header");
            if (bcf_hdr_printf(args->in_hdr, "##FILTER=<ID=UnmappedAnchors,Description=\"Anchors unmapped\">") < 0)
                error_errno("Failed to add \"UnmappedAnchors\" FILTER header");
            if (bcf_hdr_printf(args->in_hdr, "##FILTER=<ID=UnmappedAnchor5,Description=\"5' anchor unmapped\">") < 0)
                error_errno("Failed to add \"UnmappedAnchor5\" FILTER header");
            if (bcf_hdr_printf(args->in_hdr, "##FILTER=<ID=UnmappedAnchor3,Description=\"3' anchor unmapped\">") < 0)
                error_errno("Failed to add \"UnmappedAnchor3\" FILTER header");
            if (bcf_hdr_printf(args->in_hdr,
                               "##FILTER=<ID=MismatchAnchors,Description=\"Anchors mappings inconsistent\">")
                < 0)
                error_errno("Failed to add \"MismatchAnchors\" FILTER header");
            if (bcf_hdr_printf(args->in_hdr, "##FILTER=<ID=ApartAnchors,Description=\"Anchors mapping too far apart\">")
                < 0)
                error_errno("Failed to add \"ApartAnchors\" FILTER header");
            if (!args->src_fai
                && bcf_hdr_printf(
                       args->in_hdr,
                       "##FILTER=<ID=MissingFasta,Description=\"Reference allele sequence could not be extended\">")
                       < 0)
                error_errno("Failed to add \"MissingFasta\" FILTER header");
        }
        if (args->reject_fh == NULL || bcf_hdr_write(args->reject_fh, args->in_hdr) < 0)
            error("Error: cannot write to \"%s\": %s\n", reject_fname, strerror(errno));
    }

    return 0;
}

/****************************************
 * FETCH REFERENCE GENOME SEQUENCE      *
 ****************************************/

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

// fetches sequence in the same way samtools faidx would do using 1-based coordinate system
static inline char *fetch_sequence(const faidx_t *fai, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i) {
    hts_pos_t len;
    char *ref = faidx_fetch_seq64(fai, c_name, p_beg_i - 1, p_end_i - 1, &len);
    if (!ref || len != p_end_i - p_beg_i + 1) {
        free(ref);
        return NULL;
    }
    seq_to_upper(ref, len);
    replace_iupac_codes(ref, len);
    return ref;
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

static void safe_fetch_sequence(bcf1_realign_t *this) {
    int seq_len = faidx_seq_len(this->fai, bcf_seqname(this->hdr, this->rec));
    if (seq_len < 0) error("Unable to fetch sequence for contig %s\n", bcf_seqname(this->hdr, this->rec));
    if (this->beg + 1 < 1) this->beg = 0;
    if (this->end + 1 > seq_len) this->end = seq_len - 1;
    this->ref = fetch_sequence(this->fai, bcf_seqname(this->hdr, this->rec), this->beg + 1, this->end + 1);
    if (this->ref == NULL)
        error("Unable to fetch sequence at %s:%" PRIhts_pos "-%" PRIhts_pos "\n", bcf_seqname(this->hdr, this->rec),
              this->beg + 1, this->end + 1);
}

static void shift_left(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    int i, n_allele = 1;
    for (i = 1; i < rec->n_allele; i++)
        if (this->als[i].s[0] != '*') n_allele++;
    if (n_allele == 1) return;
    while (rec->pos > 0) {
        char last_ref = this->als[0].s[this->als[0].l - 1];
        for (i = 1; i < rec->n_allele; i++) {
            if (this->als[i].s[0] == '*') continue;
            if (this->als[i].s[this->als[i].l - 1] != last_ref) return;
        }
        if (this->beg == 0 || rec->pos - 1 < this->beg) {
            this->beg = rec->pos - this->aln_win;
            free(this->ref);
            safe_fetch_sequence(this);
        }
        char first_ref = this->ref[rec->pos - this->beg - 1];
        for (i = 0; i < rec->n_allele; i++) {
            if (this->als[i].s[0] == '*') continue;
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
        for (j = 1; j < rec->n_allele; j++) {
            if (this->als[j].s[0] == '*') continue;
            if (i >= this->als[j].l - 1 || this->als[j].s[i] != first_ref) break;
        }
        if (j < rec->n_allele) break;
    }
    for (j = 0; j < rec->n_allele; j++) {
        if (this->als[j].s[0] == '*') continue;
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
        for (j = 1; j < rec->n_allele; j++) {
            if (this->als[j].s[0] == '*') continue;
            if (i >= this->als[j].l - 1 || this->als[j].s[this->als[j].l - 1 - i] != last_ref) break;
        }
        if (j < rec->n_allele) break;
    }
    for (j = 0; j < rec->n_allele; j++) {
        if (this->als[j].s[0] == '*') continue;
        this->als[j].l -= i;
    }
}

static void pad_left(bcf1_realign_t *this, int npad) {
    bcf1_t *rec = this->rec;
    if (this->beg == 0 || rec->pos - npad < this->beg) {
        this->beg = rec->pos - npad;
        free(this->ref);
        safe_fetch_sequence(this);
    }
    const char *ptr = &this->ref[rec->pos - this->beg - npad];
    int i;
    for (i = 0; i < rec->n_allele; i++) {
        if (this->als[i].s[0] == '*') continue;
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
        free(this->ref);
        safe_fetch_sequence(this);
    }
    const char *ptr = &this->ref[rec->pos - this->beg + this->als[0].l];
    int i;
    for (i = 0; i < rec->n_allele; i++) {
        if (this->als[i].s[0] == '*') continue;
        kputsn(ptr, npad, &this->als[i]);
    }
}

static void pad_from_right(bcf1_realign_t *this, const char *s_ptr, int d) {
    bcf1_t *rec = this->rec;
    int npad = 0;
    const char *l_ptr = &this->ref[rec->pos - this->beg + this->als[0].l];
    while (1) {
        if (rec->pos + this->als[0].l + npad > this->end) {
            // extract more sequence from the reference
            this->end += this->aln_win;
            free(this->ref);
            safe_fetch_sequence(this);
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
    int i;
    for (i = 0; i < rec->n_allele; i++) {
        if (this->als[i].s[0] == '*') continue;
        kputsn(ptr, npad, &this->als[i]);
    }
}

static void pad_from_left(bcf1_realign_t *this, const char *s_ptr, int d) {
    bcf1_t *rec = this->rec;
    int npad = 0;
    const char *l_ptr = &this->ref[rec->pos - this->beg - 1];
    while (1) {
        if (rec->pos - npad <= this->beg) {
            // extract more sequence from the reference
            this->beg -= this->aln_win;
            free(this->ref);
            safe_fetch_sequence(this);
            l_ptr = &this->ref[rec->pos - this->beg - npad - 1];
            if (npad >= d) s_ptr = &this->ref[rec->pos - this->beg - npad + d - 1];
        }
        if (npad == d) s_ptr = &this->ref[rec->pos - this->beg - 1];
        npad++;
        if (*l_ptr != *s_ptr) break;
        l_ptr--;
        s_ptr--;
    }
    int i;
    for (i = 0; i < rec->n_allele; i++) {
        if (this->als[i].s[0] == '*') continue;
        ks_resize(&this->als[i], this->als[i].l + npad);
        memmove(this->als[i].s + npad, this->als[i].s, this->als[i].l);
        memcpy(this->als[i].s, l_ptr, npad);
        this->als[i].l += npad;
    }
    rec->pos -= npad;
}

static void initialize_alleles(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    int i;
    for (i = 0; i < rec->n_allele; i++) {
        this->als[i].l = 0;
        kputsn_(rec->d.allele[i], strlen(rec->d.allele[i]), &this->als[i]);
    }

    this->beg = rec->pos - 1;
    this->end = rec->pos + this->als[0].l;
    free(this->ref);
    safe_fetch_sequence(this);

    // make sure the reference allele matches the reference
    if (strncasecmp(&this->ref[rec->pos - this->beg], this->als[0].s, this->als[0].l) != 0)
        error("Error: the reference allele %.*s does not match the reference %.*s at %s:%" PRIhts_pos "\n",
              (int)this->als[0].l, this->als[0].s, (int)this->als[0].l, &this->ref[rec->pos - this->beg],
              bcf_seqname(this->hdr, rec), rec->pos + 1);
}

// updated the alleles in the bcf1_t structure
static void update_alleles(bcf1_realign_t *this) {
    bcf1_t *rec = this->rec;
    kstring_t *str = &this->als[0];
    int i;
    for (i = 1; i < rec->n_allele; i++) {
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

    // covers weird ClinVar cases
    if (rec->n_allele == 1 && this->als[0].l == 1) {
        pad_right(this, 1);
        pad_left(this, 1);
    }

    // find allele pairs that need to be extended
    int i, j;
    for (i = 0; i < rec->n_allele; i++) {
        if (this->als[i].s[0] == '*') continue;
        for (j = i + 1; j < rec->n_allele; j++) {
            if (this->als[j].s[0] == '*') continue;
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
 * LIFTOVER FUNCTIONS                   *
 ****************************************/

// check whether alleles in the record require extension to avoid ambiguity
static int is_extension_needed(bcf1_t *rec, int8_t **tmp_arr, int *m_tmp_arr) {
    hts_expand(int8_t, sizeof(int) * rec->n_allele, *m_tmp_arr, *tmp_arr);
    int i, j, *len = (int *)(*tmp_arr);
    for (i = 0; i < rec->n_allele; i++) len[i] = strlen(rec->d.allele[i]);
    if (rec->n_allele == 1 && len[0] == 1) return 1; // covers weird ClinVar cases

    for (i = 0; i < rec->n_allele; i++) {
        if (rec->d.allele[i][0] == '*') continue;
        for (j = i + 1; j < rec->n_allele; j++) {
            if (rec->d.allele[j][0] == '*') continue;
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
// Bioinformatics, (2015) http://doi.org/10.1093/bioinformatics/btv112
static int is_left_aligned(bcf1_t *rec, int8_t **tmp_arr, int *m_tmp_arr) {
    if (rec->n_allele == 1) return strlen(rec->d.allele[0]) <= 1;

    hts_expand(int8_t, sizeof(int) * rec->n_allele, *m_tmp_arr, *tmp_arr);
    int *len = (int *)(*tmp_arr);
    int i, min_len_is_one = 0;
    for (i = 0; i < rec->n_allele; i++) {
        if (rec->d.allele[i][0] == '*') continue;
        len[i] = strlen(rec->d.allele[i]);
        if (len[i] == 1) min_len_is_one = 1;
    }

    // 1. The alleles end with at least two different nucleotides
    char last_ref = rec->d.allele[0][len[0] - 1];
    for (i = 0; i < rec->n_allele; i++) {
        if (rec->d.allele[i][0] == '*') continue;
        if (rec->d.allele[i][len[i] - 1] != last_ref) break;
    }
    if (i == rec->n_allele) return 0;

    // 2. The alleles start with at least two different nucleotides, or the shortest allele has length 1
    if (!min_len_is_one) {
        char first_ref = rec->d.allele[0][0];
        for (i = 0; i < rec->n_allele; i++) {
            if (rec->d.allele[i][0] == '*') continue;
            if (rec->d.allele[i][0] != first_ref) break;
        }
        if (i == rec->n_allele) return 0;
    }

    return 1;
}

// lifts over a single base pair and returns the block index
static int liftover_bp(regidx_t *idx, regitr_t *itr, const char *t_chr, hts_pos_t t_pos, int *q_rid, hts_pos_t *q_pos,
                       int *q_strand) {
    int block_ind = -1;
    if (regidx_overlap(idx, t_chr, (uint32_t)t_pos, (uint32_t)t_pos, itr)) {
        int i;
        for (i = 0; regitr_overlap(itr); i++) {
            if (i > 0)
                fprintf(stderr, "Warning: more than one contiguous block overlaps with position %s:%" PRIhts_pos "\n",
                        t_chr, t_pos);
            block_ind = regitr_payload(itr, int);
            block_t *block = &args->blocks[block_ind];
            assert(block->size == itr->end - itr->beg + 1);
            chain_t *chain = &args->chains[block->chain_ind];

            *q_rid = chain->q_rid;
            *q_strand = chain->qStrand;
            int block_pos = t_pos - itr->beg;
            if (*q_strand) // - strand
                *q_pos = (hts_pos_t)(chain->qSize - chain->qStart - block->qStart - block_pos);
            else // + strand
                *q_pos = (hts_pos_t)(chain->qStart + block->qStart + block_pos + 1);
        }
    }
    return block_ind;
}

// position are 1-based
// returns -1 if neither anchor can be mapped
// returns -2 if the 5' anchor is not mappable
// returns -3 if the 3' anchor is not mappable
// returns -4 if anchors mapped but in an inconsistent way
// returns -5 if anchors mapped too far from each other
static int liftover_indel(regidx_t *idx, regitr_t *itr, const char *src_chr, hts_pos_t src_pos5, hts_pos_t src_pos3,
                          int max_indel_inc, int *dst_rid, hts_pos_t *dst_pos5, hts_pos_t *dst_pos3, int *strand,
                          int *npad) {
    int rid5 = -1, rid3 = -1, strand5 = 0, strand3 = 0;
    hts_pos_t pos5 = -1, pos3 = -1;
    int block_ind5 = liftover_bp(idx, itr, src_chr, src_pos5, &rid5, &pos5, &strand5);
    const block_t *block5 = block_ind5 < 0 ? NULL : &args->blocks[block_ind5];
    int block_ind3 = liftover_bp(idx, itr, src_chr, src_pos3, &rid3, &pos3, &strand3);
    const block_t *block3 = block_ind3 < 0 ? NULL : &args->blocks[block_ind3];
    if (!block5 && !block3) return -1; // both anchors of the indel failed to lift over

    // strategy to pad alleles to get to the very first base on the next contiguous block
    if (!block5) { // pad sequence to the left
        const chain_t *aux_chain = &args->chains[block3->chain_ind];
        int dst_chr_size = aux_chain->qSize;
        *dst_rid = rid3;
        *strand = strand3;

        // identify new 5' anchor
        const block_t *aux_block = prev_block(block3, args->blocks, args->chains); // identify the previous block
        *npad = aux_block ? aux_chain->tStart + aux_block->tStart + aux_block->size - src_pos5 : -max_indel_inc;
        if (*npad < -max_indel_inc || *npad >= 0) {
            *npad = -max_indel_inc; // rare cases where the indel spans the whole previous block
            if (*npad > src_pos5 - src_pos3) return -2;
        }
        if (src_pos5 + *npad < 1) *npad = 1 - src_pos5; // hit left edge on the source chromosome

        // attempt to liftover the new anchor
        block_ind5 = liftover_bp(idx, itr, src_chr, src_pos5 + *npad, &rid5, &pos5, &strand5);
        int dst_npad = *strand ? pos5 - pos3 : pos3 - pos5; // check that the other anchor did not liftover too far
        if (rid5 != rid3 || strand5 != strand3 || dst_npad < 0 || dst_npad > max_indel_inc) {
            if (*strand) {
                // hit right edge on the destination chromosome
                if (pos3 - *npad > dst_chr_size) *npad = pos3 - dst_chr_size;
                pos5 = pos3 - *npad; // dst_chr_size
            } else {
                // hit left edge on the destination chromosome
                if (pos3 - (src_pos3 - src_pos5) + *npad < 1) *npad = 1 - pos3 + (src_pos3 - src_pos5);
                pos5 = pos3 - (src_pos3 - src_pos5) + *npad; // 1
            }
        }
    } else if (!block3) { // pad sequence to the right
        const chain_t *aux_chain = &args->chains[block5->chain_ind];
        int src_chr_size = aux_chain->tSize;
        int dst_chr_size = aux_chain->qSize;
        *dst_rid = rid5;
        *strand = strand5;

        // identify new 3' anchor
        const block_t *aux_block = next_block(block5, args->blocks, args->chains); // identify the next block
        *npad = aux_block ? aux_chain->tStart + aux_block->tStart + 1 - src_pos3 : max_indel_inc;
        if (*npad > max_indel_inc || *npad <= 0) {
            *npad = max_indel_inc; // rare cases where the indel spans the whole next block
            if (*npad < src_pos3 - src_pos5) return -3;
        }
        if (src_pos3 + *npad > src_chr_size) *npad = src_chr_size - src_pos3; // hit right edge on the source chromosome

        // attempt to liftover the new anchor
        block_ind3 = liftover_bp(idx, itr, src_chr, src_pos3 + *npad, &rid3, &pos3, &strand3);
        int dst_npad = *strand ? pos5 - pos3 : pos3 - pos5; // check that the other anchor did not liftover too far
        if (rid5 != rid3 || strand5 != strand3 || dst_npad < 0 || dst_npad > max_indel_inc) {
            if (*strand) {
                // hit left edge on the destination chromosome
                if (pos5 - *npad < 1) *npad = pos5 - 1;
                pos3 = pos5 - *npad; // 1
            } else {
                // hit rigth edge on the destination chromosome
                if (pos5 + (src_pos3 - src_pos5) + *npad > dst_chr_size)
                    *npad = dst_chr_size - pos5 - (src_pos3 - src_pos5);
                pos3 = pos5 + (src_pos3 - src_pos5) + *npad; // dst_chr_size
            }
        }
    } else {
        if (block5->chain_ind != block3->chain_ind) return -4;
        if (rid5 != rid3) return -4;
        *dst_rid = rid5;
        if (strand5 != strand3) return -4;
        *strand = strand5;
        if (abs((int)(pos3 - pos5)) > src_pos3 - src_pos5 + max_indel_inc) return -5;
    }

    *dst_pos5 = *strand == 0 ? pos5 : pos3;
    *dst_pos3 = *strand == 0 ? pos3 : pos5;
    assert(*dst_pos5 <= *dst_pos3);
    return 0;
}

static int find_reference(bcf1_t *rec, const bcf_hdr_t *hdr, const char *ref, int is_snp, kstring_t *str) {
    int i, swap = -1;
    int ref_len = strlen(ref);
    for (i = 0; i < rec->n_allele; i++) {
        if (rec->d.allele[i][0] == '*') continue;
        int len = strlen(rec->d.allele[i]);
        if (len != ref_len) continue;
        if (ref_len == 1 ? (ref[0] == rec->d.allele[i][0]) : !strncasecmp(ref + 1, rec->d.allele[i] + 1, len - 2)) {
            if (swap >= 0) {
                // two matches can only happen if alleles were not maximally extended
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

    // this is a special hack to avoid SNPs becoming complex variants (e.g. rs17554556)
    // when the new assembly is longer due to an insertion right or left of the SNP
    // and the anchors are ambiguous (saves ~30% of SNPs from becoming complex variants)
    if (swap < 0 && is_snp && ref_len > 2 && ref_len != strlen(rec->d.allele[0])) {
        int len, left_match, right_match, n_matches = 0;
        for (i = 0; i < rec->n_allele; i++) {
            if (rec->d.allele[i][0] == '*') continue;
            len = strlen(rec->d.allele[i]);
            if (len < ref_len) {
                left_match = !strncasecmp(ref, rec->d.allele[i], len - 1);
                right_match = !strncasecmp(ref + 1 + ref_len - len, rec->d.allele[i] + 1, len);
                n_matches += left_match + right_match;
                if (n_matches > 1) break; // too many matches
                if (left_match || right_match) swap = i;
            }
        }
        if (n_matches == 1) {
            len = strlen(rec->d.allele[swap]);
            left_match = !strncasecmp(ref, rec->d.allele[swap], len - 1);
            right_match = !strncasecmp(ref + 1 + ref_len - len, rec->d.allele[swap] + 1, len);
            if (left_match) {
                ref_len = len;
            } else { // right match
                ref += ref_len - len;
                rec->pos += ref_len - len;
                ref_len = len;
            }
        } else {
            swap = -1;
        }
    }

    // if needed, update left and right anchors
    // this is because sometimes the anchors can include SNPs
    // that should not be incorporated in the liftover
    if (ref_len > 1) { // don't update if the variant is a SNP
        char ref5 = ref[0];
        char ref3 = ref[ref_len - 1];
        for (i = 0; i < rec->n_allele; i++) {
            if (rec->d.allele[i][0] == '*') continue;
            int len = strlen(rec->d.allele[i]);
            rec->d.allele[i][len - 1] = ref3;
            rec->d.allele[i][0] = ref5;
        }
    }

    // update alleles if necessary
    if (swap) {
        if (swap < 0) {
            str->l = 0;
            kputsn_(ref, ref_len, str);
            kputc_(',', str);
            for (i = 0; i < rec->n_allele; i++) {
                char *allele = rec->d.allele[i];
                kputsn_(allele, strlen(allele), str);
                kputc_(',', str);
            }
            str->l--;
            str->s[str->l] = '\0';
            bcf_update_alleles_str(hdr, rec, str->s);
        } else {
            char *tmp = rec->d.allele[0];
            rec->d.allele[0] = rec->d.allele[swap];
            rec->d.allele[swap] = tmp;
            bcf_update_alleles(hdr, rec, (const char **)rec->d.allele, rec->n_allele);
        }
    }

    return swap;
}

/****************************************
 * NEEDLEMAN–WUNSCH FUNCTIONS           *
 ****************************************/

#define _M 0
#define _D 1
#define _I 2
// this version left-aligns deletions and insertions
#define MAX_IND(v) ((v[_D]) < (v[_I]) ? (v[_I]) < (v[_M]) ? (_M) : (_I) : (v[_D]) < (v[_M]) ? (_M) : (_D))

// performs sequence alignment using the Needleman–Wunsch algorithm using affine gap penalty
// the scoring values for matches and gap penalties are drawn from bwa mem
// return a sequence of values of length at most s_l + t_l defining the sequence match
static int nw(const char *s, size_t s_l, const char *t, size_t t_l, kstring_t *path) {
    if (s_l < 0 || t_l < 0) return -1;
    int a = 1;  // score for a sequence match
    int b = -4; // penalty for a mismatch
    int o = -6; // gap open penalties for deletions and insertions
    int e = -1; // gap extension penalty
    int n = s_l + 1;
    int m = t_l + 1;
    int i, j, cell, ind, score[3];
    int *A[3];
    char *B[3];
    ks_resize(path, n + m - 1);
    for (i = 0; i < 3; i++) {
        A[i] = (int *)malloc(sizeof(int) * n * m);
        B[i] = (char *)malloc(sizeof(char) * n * m);
    }

    // build three scoring matrices corresponding to the three possible states
    A[_M][0] = 0;
    A[_D][0] = INT32_MIN / 2;
    A[_I][0] = INT32_MIN / 2;
    for (j = 1; j < m; j++) {
        A[_M][j] = INT32_MIN / 2;
        A[_D][j] = INT32_MIN / 2;
        A[_I][j] = o + e * j;
        B[_I][j] = _I;
    }
    for (i = 1; i < n; i++) {
        A[_M][i * m] = INT32_MIN / 2;
        A[_D][i * m] = o + e * i;
        B[_D][i * m] = _D;
        A[_I][i * m] = INT32_MIN / 2;
        for (j = 1; j < m; j++) {
            cell = i * m + j;
            score[_M] = A[_M][cell - m - 1];
            score[_D] = A[_D][cell - m - 1];
            score[_I] = A[_I][cell - m - 1];
            ind = MAX_IND(score);
            A[_M][cell] = (s[i - 1] == t[j - 1] || s[i - 1] == 'N' || t[j - 1] == 'N' ? a : b) + score[ind];
            B[_M][cell] = ind;

            score[_M] = o + A[_M][cell - m];
            score[_D] = A[_D][cell - m];
            score[_I] = o + A[_I][cell - m];
            ind = MAX_IND(score);
            A[_D][cell] = e + score[ind];
            B[_D][cell] = ind;

            score[_M] = o + A[_M][cell - 1];
            score[_D] = o + A[_D][cell - 1];
            score[_I] = A[_I][cell - 1];
            ind = MAX_IND(score);
            A[_I][cell] = e + score[ind];
            B[_I][cell] = ind;
        }
    }

    // adjust the values at the end to account for 3' clipping
    cell = m * n - 1;
    score[_M] = A[_M][cell];
    score[_D] = A[_D][cell];
    score[_I] = A[_I][cell];
    ind = MAX_IND(score);
    int ret = score[ind];

    // reconstructs the path by backtracking
    i = n + m - 2;
    path->s[i] = '\0';
    while (cell > 0) {
        path->s[--i] = 1 + ind;
        int shift = (ind == _M || ind == _D ? m : 0) + (ind == _M || ind == _I ? 1 : 0);
        ind = B[ind][cell];
        cell -= shift;
    }
    assert(cell == 0);

    // i here represents the number of matches
    if (i) memmove(path->s, &path->s[i], n + m - 1 - i);
    path->l = n + m - i;

    for (i = 0; i < 3; i++) {
        free(A[i]);
        free(B[i]);
    }

    return ret;
}

static int get_shift(char *path, int npad) {
    int shift = 0;
    if (npad < 0) {
        char *ptr = path;
        while (npad < 0) {
            int ind = (*ptr++) - 1;
            if (ind == _M || ind == _D) shift++;
            if (ind == _M || ind == _I) npad++;
        }
    } else {
        char *ptr = &path[strlen(path) - 1];
        while (npad > 0) {
            int ind = (*ptr--) - 1;
            if (ind == _M || ind == _D) shift--;
            if (ind == _M || ind == _I) npad--;
        }
    }
    return shift;
}

// remove the pad by re-aligning the gap using Needleman-Wunsch
// this does not happen very often but it is necessary to recognize what
// the previous reference interval corresponds to in the new reference
static void clip_pad(bcf1_t *rec, const char *ref, char *pad, int npad, hts_pos_t *pos5, hts_pos_t *pos3,
                     kstring_t *str, int write_nw) {
    int score = INT32_MIN;
    int i, shift = 0;
    kstring_t path = {0, 0, NULL};
    kstring_t nw_s = {0, 0, NULL};

    for (i = 0; i < rec->n_allele; i++) {
        if (rec->d.allele[i][0] == '*') continue;

        // generate source reference sequence
        str->l = 0;
        if (npad < 0) {
            kputsn_(pad, -npad, str);
            kputs(rec->d.allele[i], str);
        } else {
            kputsn_(rec->d.allele[i], strlen(rec->d.allele[i]), str);
            kputsn(pad, npad, str);
        }

        // pairwise align the source and destination sequences excluding the anchors
        int new_score = nw(ref + 1, *pos3 - *pos5 - 1, str->s + 1, str->l - 2, &path);

        if (write_nw) {
            if (i > 0) kputc_(',', &nw_s);
            const char *p = path.s;
            const char *s = ref + 1;
            for (; *p; p++) {
                int ind = *p - 1;
                kputc_(ind == _M || ind == _D ? *s++ : '-', &nw_s);
            }
            kputc_('|', &nw_s);
            p = path.s;
            const char *t = str->s + 1;
            for (; *p; p++) {
                int ind = *p - 1;
                kputc_(ind == _M || ind == _I ? *t++ : '-', &nw_s);
            }
        }

        if (new_score > score) {
            score = new_score;
            shift = get_shift(path.s, npad);
        }
    }

    if (write_nw) {
        str->l = 0;
        kputsn(nw_s.s, nw_s.l, str);
        free(nw_s.s);
    }

    if (npad < 0)
        *(pos5) += shift;
    else
        *(pos3) += shift;
    free(path.s);
}

/****************************************
 * RECORD UPDATING FUNCTIONS            *
 ****************************************/

// this function updates FORMAT/GT records if the number or order of alleles has changed
static void update_genotypes(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int swap, int32_t **int32_arr,
                             int *m_int32_arr) {
    int i;
    if (swap < 0) {
        int ngt = bcf_get_genotypes(args->out_hdr, rec, &args->int32_arr, &args->m_int32_arr);
        if (ngt <= 0) return;
        int *gts = (int *)(args->int32_arr);
        for (i = 0; i < ngt; i++)
            if (!bcf_gt_is_missing(gts[i]) && !(gts[i] == bcf_int32_vector_end)) gts[i] += 2;
        bcf_update_genotypes(hdr, rec, gts, ngt);
    } else {
        bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
        if (!fmt) return;
        int n = rec->n_sample * fmt->n;
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (i = 0; i < n; i++) {                                                                                      \
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

static void update_ploidy(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int *ploidy_arr) {
    int i, j;
    for (i = 0; i < bcf_hdr_nsamples(hdr); i++) ploidy_arr[i] = 2;
    bcf_fmt_t *fmt = bcf_get_fmt_id(rec, int_id);
    if (!fmt) {
        return;
        if (rec->n_sample != bcf_hdr_nsamples(hdr))
            error("Number %d of samples in the VCF record does not match the number %d of samples in the header\n",
                  rec->n_sample, bcf_hdr_nsamples(hdr));
#define BRANCH(type_t, vector_end)                                                                                     \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (i = 0; i < rec->n_sample; i++) {                                                                          \
            ploidy_arr[i] = fmt->n;                                                                                    \
            for (j = 0; j < fmt->n; j++)                                                                               \
                if (p[j] == vector_end) {                                                                              \
                    ploidy_arr[i] = j;                                                                                 \
                    break;                                                                                             \
                }                                                                                                      \
            p += fmt->n;                                                                                               \
        }                                                                                                              \
    }
        switch (fmt->type) {
        case BCF_BT_INT8:
            BRANCH(int8_t, bcf_int8_vector_end);
            break;
        case BCF_BT_INT16:
            BRANCH(int16_t, bcf_int16_vector_end);
            break;
        case BCF_BT_INT32:
            BRANCH(int32_t, bcf_int32_vector_end);
            break;
        default:
            error("Unexpected type %d\n", fmt->type);
        }
#undef BRANCH
    }
    return;
}

static void update_genotype_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int swap,
                                   int8_t **tmp_arr, int *m_tmp_arr) {
    int i, n, type, len;
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

    hts_expand(int8_t, sizeof(int32_t) * n * len, *m_tmp_arr, *tmp_arr);

    if (swap < 0) {
#define BRANCH(type_t, is_missing, is_vector_end)                                                                      \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        int32_t *dst = (int32_t *)(*tmp_arr);                                                                          \
        for (i = 0; i < n * len; i++) {                                                                                \
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
        for (i = 0; i < n * len; i++) {                                                                                \
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
static void update_AGR_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int swap, int8_t **tmp_arr,
                              int *m_tmp_arr) {
    int i, j, k, n, type, len;
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
    hts_expand(int8_t, size, *m_tmp_arr, *tmp_arr);

    if (swap < 0) {               // add reference allele entry
        if (length == BCF_VL_G) { // define how to add an element to a VL_G array
            if (type != BCF_BT_CHAR && len != (n_als - 1) * n_als / 2)
                error(
                    "Number of elements in the VCF record %s should be %d but %d found\nUse --drop-tags %s/%s to skip "
                    "this error\n",
                    key, (n_als - 1) * n_als / 2, len, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);

#define BRANCH(type_t, set_missing, set_vector_end, is_missing, is_vector_end, out_type_t)                             \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        out_type_t *dst = (out_type_t *)(*tmp_arr);                                                                    \
        for (i = 0; i < n; i++) {                                                                                      \
            for (j = 0; j < n_als; j++) {                                                                              \
                set_missing;                                                                                           \
                dst++;                                                                                                 \
                for (k = 0; k < j; k++) {                                                                              \
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
            if (type != BCF_BT_CHAR && len != n_als - 1 - (length == BCF_VL_A))
                error(
                    "Number of elements in the VCF record %s should be %d but %d found\nUse --drop-tags %s/%s to skip "
                    "this error\n",
                    key, n_als - 1 - (length == BCF_VL_A), len, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
#define BRANCH(type_t, set_missing, set_vector_end, is_missing, is_vector_end, out_type_t)                             \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        out_type_t *dst = (out_type_t *)(*tmp_arr);                                                                    \
        for (i = 0; i < n; i++) {                                                                                      \
            set_missing;                                                                                               \
            dst++;                                                                                                     \
            for (j = 0; j < len; j++) {                                                                                \
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
                for (i = 0; i < n; i++) {
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
            if (type != BCF_BT_CHAR && len != n_als * (n_als + 1) / 2)
                error(
                    "Number of elements in the VCF record %s should be %d but %d found\nUse --drop-tags %s/%s to skip "
                    "this error\n",
                    key, n_als * (n_als + 1) / 2, len, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (i = 0; i < n; i++) {                                                                                      \
            type_t tmp = ptr[0];                                                                                       \
            ptr[0] = ptr[swap * (swap + 3) / 2];                                                                       \
            ptr[swap * (swap + 3) / 2] = tmp;                                                                          \
            for (j = 1; j < swap; j++) {                                                                               \
                tmp = ptr[j * (j + 1) / 2];                                                                            \
                ptr[j * (j + 1) / 2] = ptr[swap * (swap + 1) / 2 + j];                                                 \
                ptr[swap * (swap + 1) / 2 + j] = tmp;                                                                  \
            }                                                                                                          \
            for (j = swap + 1; j < n_als; j++) {                                                                       \
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
            if (type != BCF_BT_CHAR && len != n_als)
                error(
                    "Number of elements in the VCF record %s should be %d but %d found\nUse --drop-tags %s/%s to skip "
                    "this error\n",
                    key, n_als, len, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
#define BRANCH(type_t)                                                                                                 \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (i = 0; i < n; i++) {                                                                                      \
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
                for (i = 0; i < n; i++) {
                    memcpy(s, p, len);
                    char *ptr1 = strchr(s, ',');
                    if (!ptr1) error("Error: string %s contains less then %d elements", s, n_als);
                    char *ptr2 = ptr1;
                    for (j = 1; ptr2 && j < swap; j++) ptr2 = strchr(ptr2 + 1, ',');
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

// this function updates INFO/FORMAT Number=A Type=Integer/Float records if the number or order of alleles has changed
static void reverse_A_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int *sum, int swap,
                             int8_t **tmp_arr, int *m_tmp_arr) {
    int i, j, n, type, len;
    void *p;
    int n_als = rec->n_allele;
    int length = bcf_hdr_id2length(hdr, coltype, int_id);   // it should be BCF_VL_A
    int field_type = bcf_hdr_id2type(hdr, coltype, int_id); // it should be BCF_HT_INT or BCF_HT_REAL
    const char *key = bcf_hdr_int2id(hdr, BCF_DT_ID, int_id);

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
    assert(length == BCF_VL_A && type != BCF_BT_CHAR);

    int size = (type == BCF_HT_INT ? sizeof(int32_t) : sizeof(float)) * n * (n_als - 1);
    hts_expand(int8_t, size, *m_tmp_arr, *tmp_arr);

    if (len != n_als - 1 - (swap < 0))
        error(
            "Number of elements in the VCF record %s should be %d but %d found\nUse --drop-tags %s/%s to skip this "
            "error\n",
            key, n_als - 1 - (swap < 0), len, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
#define BRANCH(type_t, set_missing, set_vector_end, is_missing, is_vector_end, out_type_t)                             \
    {                                                                                                                  \
        type_t *src = (type_t *)p;                                                                                     \
        out_type_t *dst = (out_type_t *)(*tmp_arr);                                                                    \
        for (i = 0; i < n; i++) {                                                                                      \
            out_type_t value = (out_type_t)sum[i];                                                                     \
            out_type_t *swap_dst = swap < 0 ? dst++ : &dst[swap - 1];                                                  \
            for (j = 0; j < len; j++) {                                                                                \
                if (is_missing)                                                                                        \
                    set_missing;                                                                                       \
                else if (is_vector_end)                                                                                \
                    set_vector_end;                                                                                    \
                else {                                                                                                 \
                    value -= (out_type_t) * src;                                                                       \
                    *dst = *src;                                                                                       \
                }                                                                                                      \
                src++;                                                                                                 \
                dst++;                                                                                                 \
            }                                                                                                          \
            *swap_dst = value;                                                                                         \
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
    default:
        fprintf(stderr, "TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type);
        exit(1);
    }
#undef BRANCH

    len = n * (n_als - 1);
    if (coltype == BCF_HL_INFO) // INFO fields
        bcf_update_info(hdr, rec, key, *tmp_arr, len, field_type);
    else // FORMAT fields
        bcf_update_format(hdr, rec, key, *tmp_arr, len, field_type);
}

// this function updates INFO/FORMAT records if the number or order of alleles has changed
static void flip_A_record(const bcf_hdr_t *hdr, bcf1_t *rec, int int_id, int coltype, int swap) {
    int i, j, n, type, len;
    float *p;
    int n_als = rec->n_allele;
    int length = bcf_hdr_id2length(hdr, coltype, int_id);
    const char *key = bcf_hdr_int2id(hdr, BCF_DT_ID, int_id);
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

    if (type != BCF_BT_CHAR && len != n_als - 1)
        error(
            "Number of elements in the VCF record %s should be %d but %d found\nUse --drop-tags %s/%s to skip this "
            "error\n",
            key, n_als - 1, len, coltype == BCF_HL_INFO ? "INFO" : "FORMAT", key);
#define BRANCH(type_t, easy_access)                                                                                    \
    {                                                                                                                  \
        type_t *ptr = (type_t *)p;                                                                                     \
        for (i = 0; i < n; i++) {                                                                                      \
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
        for (i = 0; i < n; i++) {
            char *ptr = s;
            if (swap > 1) {
                ptr = strchr(ptr, ',');
                for (j = 0; ptr && j < swap - 2; j++) ptr = strchr(ptr + 1, ',');
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

static inline int bcf_is_symbolic(const bcf1_t *rec) {
    if (rec->d.allele[0][0] == '*') return 1;
    int i;
    for (i = 0; i < rec->n_allele; i++)
        if (rec->d.allele[i][0] == '<') return 1;
    return 0;
}

/****************************************
 * PROCESS RECORDS                      *
 ****************************************/

bcf1_t *process(bcf1_t *rec) {
    int i, is_snp = bcf_is_snp(rec);
    int is_symbolic = bcf_is_symbolic(rec);
    if (is_symbolic && !args->warning_symbolic) {
        fprintf(stderr, "Warning: input VCF includes symbolic alleles that might not properly lift over\n");
        args->warning_symbolic = 1;
    }

    if (rec->errcode == BCF_ERR_CTG_UNDEF)
        fprintf(stderr, "Warning: variants from contig %s cannot be lift over\n", bcf_seqname(args->in_hdr, rec));

    const char *src_chr = bcf_hdr_id2name(args->in_hdr, rec->rid);
    hts_pos_t src_pos = rec->pos + 1;
    if (args->reject_fh || args->write_src) {
        args->tmp_kstr.l = 0;
        for (i = 0; i < rec->n_allele; i++) {
            char *allele = rec->d.allele[i];
            int len = strlen(allele);
            kputsn_(allele, len, &args->tmp_kstr);
            kputc_(',', &args->tmp_kstr);
        }
        args->tmp_kstr.l--;
        args->tmp_kstr.s[args->tmp_kstr.l] = '\0';
    }

    // lift over record coordinates
    int dst_rid = -1;
    hts_pos_t dst_pos5 = -1; // 1-based coordinate system
    hts_pos_t dst_pos3 = -1; // 1-based coordinate system
    int ret = 0, strand, is_difficult_snp = 0, npad = 0;
    if (args->idx) {
        ret = liftover_bp(args->idx, args->itr, src_chr, rec->pos + 1, &dst_rid, &dst_pos5, &strand);
        dst_pos3 = dst_pos5;
        is_difficult_snp = ret < 0;
    }

    if ((!is_snp || is_difficult_snp) && !is_symbolic) {
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

    // if no chain file was provided return the record as is
    if (!args->idx) return rec;
    args->ntotal++;

    if (!args->lift_mt && rec->rid == args->in_mt_rid) {
        rec->rid = args->out_mt_rid;
        return rec;
    }

    hts_pos_t src_pos5 = rec->pos + 1;                        // 1-based coordinate system
    hts_pos_t src_pos3 = rec->pos + strlen(rec->d.allele[0]); // 1-based coordinate system
    if ((!is_snp || is_difficult_snp) && !is_symbolic) {
        ret = liftover_indel(args->idx, args->itr, src_chr, src_pos5, src_pos3, args->max_indel_inc, &dst_rid,
                             &dst_pos5, &dst_pos3, &strand, &npad);
    }
    if (ret < 0 || (npad && !args->src_fai)) {
        args->nrejected++;
        if (args->reject_fh) {
            // restore original position and alleles
            if ((!is_snp || is_difficult_snp) && !is_symbolic) {
                rec->pos = src_pos - 1;
                bcf_update_alleles_str(args->in_hdr, rec, args->tmp_kstr.s);
            }
            // include explanation for rejection
            if (args->reject_filter) {
                if (rec->rid >= args->n_ctgs) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "MissingContig"));
                } else if (ret == -1) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "UnmappedAnchors"));
                } else if (ret == -2) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "UnmappedAnchor5"));
                } else if (ret == -3) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "UnmappedAnchor3"));
                } else if (ret == -4) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "MismatchAnchors"));
                } else if (ret == -5) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "ApartAnchors"));
                } else if (ret == 0) {
                    bcf_add_filter(args->in_hdr, rec, bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "MissingFasta"));
                }
            }
            if (bcf_write(args->reject_fh, args->in_hdr, rec) < 0) error("Error: Unable to write to output VCF file\n");
        }
        return NULL;
    }

    if (args->write_src) {
        bcf_update_info_string(args->out_hdr, rec, "SRC_CHROM", src_chr);
        bcf_update_info_int32(args->out_hdr, rec, "SRC_POS", &src_pos, 1);
        bcf_update_info_string(args->out_hdr, rec, "SRC_REF_ALT", args->tmp_kstr.s);
    }

    if (args->write_fail) {
        if (npad < 0) bcf_update_info_flag(args->out_hdr, rec, "FAIL5", NULL, 1);
        if (npad > 0) bcf_update_info_flag(args->out_hdr, rec, "FAIL3", NULL, 1);
    }

    // collect the required pad sequence in case we might need it
    if (npad && args->src_fai) {
        hts_pos_t npad_pos5 = npad < 0 ? src_pos5 + npad : src_pos3 + 1;
        hts_pos_t npad_pos3 = npad < 0 ? src_pos5 - 1 : src_pos3 + npad;
        char *ref = fetch_sequence(args->src_fai, bcf_seqname(args->in_hdr, rec), npad_pos5, npad_pos3);
        if (ref == NULL)
            error("Unable to fetch sequence from the source reference at %s:%" PRIhts_pos "-%" PRIhts_pos
                  " while processing variant at position %s:%" PRIhts_pos "\n",
                  bcf_seqname(args->in_hdr, rec), npad_pos5, npad_pos3, src_chr, src_pos);
        args->tmp_pad.l = 0;
        kputs(ref, &args->tmp_pad);
        free(ref);
    }

    // flip alleles
    if (strand) {
        bcf_update_info_flag(args->out_hdr, rec, args->flip_tag, NULL, 1);
        for (i = 0; i < rec->n_allele; i++) {
            char *allele = rec->d.allele[i];
            if (allele[0] != '<' && allele[0] != '*') reverse_complement(allele);
        }
        bcf_update_alleles(args->out_hdr, rec, (const char **)rec->d.allele, rec->n_allele);
        if (npad) {
            reverse_complement(args->tmp_pad.s);
            npad = -npad;
        }
    }

    // adjust the destination reference allele if padding was required during liftover
    rec->rid = dst_rid;
    if (npad && dst_pos3 > dst_pos5) {
        char *ref = fetch_sequence(args->dst_fai, bcf_seqname(args->out_hdr, rec), dst_pos5, dst_pos3);
        if (ref == NULL)
            error("Unable to fetch sequence from the destination reference at %s:%" PRIhts_pos "-%" PRIhts_pos
                  " while processing variant at position %s:%" PRIhts_pos "\n",
                  bcf_seqname(args->out_hdr, rec), dst_pos5, dst_pos3, src_chr, src_pos);
        clip_pad(rec, ref, args->tmp_pad.s, npad, &dst_pos5, &dst_pos3, &args->tmp_kstr, args->write_nw);
        if (args->write_nw) bcf_update_info_string(args->out_hdr, rec, "NW", args->tmp_kstr.s);
        free(ref);
    }
    rec->pos = dst_pos5 - 1;

    int swap = 0;
    if (!is_symbolic) {
        char *ref = fetch_sequence(args->dst_fai, bcf_seqname(args->out_hdr, rec), dst_pos5, dst_pos3);
        if (ref == NULL)
            error("Unable to fetch sequence from the destination reference at %s:%" PRIhts_pos "-%" PRIhts_pos
                  " while processing variant at position %s:%" PRIhts_pos "\n",
                  bcf_seqname(args->out_hdr, rec), dst_pos5, dst_pos3, src_chr, src_pos);
        swap = find_reference(rec, args->out_hdr, ref, is_snp, &args->tmp_kstr);
        free(ref);
    }

    // left align indels
    if ((!is_snp || is_difficult_snp) && !is_symbolic && !args->no_left_align
        && !is_left_aligned(rec, &args->tmp_arr, &args->m_tmp_arr)) {
        hts_expand0(kstring_t, rec->n_allele, args->m_tmp_als, args->tmp_als);
        bcf1_realign_t this = {args->out_hdr, rec, args->tmp_als, args->dst_fai, args->aln_win, NULL, 0, 0};
        initialize_alleles(&this);
        trim_right(&this);
        trim_left(&this);
        shift_left(&this);
        update_alleles(&this);
    }

    // update INFO/END field if present
    if (!is_symbolic) rec->rlen = strlen(rec->d.allele[0]);
    if (args->info_end_id >= 0) {
        bcf_info_t *end_info = bcf_get_info_id(rec, args->info_end_id);
        if (end_info) {
            if (args->lift_end) {
                int end_rid, end_strand;
                hts_pos_t end_pos;
                ret = liftover_bp(args->idx, args->itr, src_chr, end_info->v1.i, &end_rid, &end_pos, &end_strand);
                if (ret >= 0 && end_rid == rec->rid && end_strand == strand
                    && ((strand && end_pos <= rec->pos + 1) || (!strand && end_pos >= rec->pos + 1)))
                    end_info->v1.i = end_pos;
                else
                    end_info->v1.i = bcf_int64_missing;
            } else {
                end_info->v1.i = rec->pos + rec->rlen;
            }
        }
    }

    if (!swap) return rec;

    // address records that have the reference allele added or swapped
    if (swap < 0)
        args->nref_added++;
    else
        args->nswapped++;
    bcf_update_info_int32(args->out_hdr, rec, args->swap_tag, &swap, 1);

    // address genotypes
    update_genotypes(args->out_hdr, rec, args->fmt_gt_id, swap, &args->int32_arr, &args->m_int32_arr);
    // estimate ploidy
    int ploidy = 2;
    update_ploidy(args->out_hdr, rec, args->fmt_gt_id, args->ploidy_arr);

    // extract INFO/AN and FORMAT/AN information
    int an = -1;
    if (args->info_an_id >= 0) {
        bcf_info_t *an_info = bcf_get_info_id(rec, args->info_an_id);
        if (an_info) an = (float)an_info->v1.i;
    }
    int *an_arr = NULL;
    if (args->fmt_an_id >= 0) {
        bcf_fmt_t *an_fmt = bcf_get_fmt_id(rec, args->fmt_an_id);
        if (an_fmt) {
            int n_an = bcf_get_format_int32(args->out_hdr, rec, "AN", &args->int32_arr, &args->m_int32_arr);
            if (n_an >= 0) {
                if (n_an != bcf_hdr_nsamples(args->out_hdr))
                    error(
                        "Number %d of FORMAT/AN values in the VCF record does not match the number %d of samples in "
                        "the header\n",
                        n_an, bcf_hdr_nsamples(args->out_hdr));
                an_arr = (int *)args->int32_arr;
            }
        }
    }
    int af = 1;

    for (i = 0; i < args->n_tags; i++) {
        tag_t *tag = &args->tags[i];

        int field_type = bcf_hdr_id2type(args->out_hdr, tag->coltype, tag->int_id);
        const char *key = bcf_hdr_int2id(args->out_hdr, BCF_DT_ID, tag->int_id);
        switch (tag->rule) {
        case RULE_DROP:
            if (tag->coltype == BCF_HL_INFO)
                bcf_update_info(args->out_hdr, rec, key, NULL, 0, field_type);
            else
                bcf_update_format(args->out_hdr, rec, key, NULL, 0, field_type);
            break;
        case RULE_AC:
            if ((tag->coltype == BCF_HL_INFO && an < 0) || (tag->coltype == BCF_HL_FMT && !an_arr))
                update_AGR_record(args->out_hdr, rec, tag->int_id, tag->coltype, swap, &args->tmp_arr,
                                  &args->m_tmp_arr);
            else
                reverse_A_record(args->out_hdr, rec, tag->int_id, tag->coltype,
                                 tag->coltype == BCF_HL_INFO ? &an : an_arr, swap, &args->tmp_arr, &args->m_tmp_arr);
            break;
        case RULE_AF:
            reverse_A_record(args->out_hdr, rec, tag->int_id, tag->coltype, BCF_HL_INFO ? &af : args->af_arr, swap,
                             &args->tmp_arr, &args->m_tmp_arr);
            break;
        case RULE_DS:
            reverse_A_record(args->out_hdr, rec, tag->int_id, tag->coltype,
                             tag->coltype == BCF_HL_INFO ? &ploidy : args->ploidy_arr, swap, &args->tmp_arr,
                             &args->m_tmp_arr);
            break;
        case RULE_GT:
            update_genotype_record(args->out_hdr, rec, tag->int_id, tag->coltype, swap, &args->tmp_arr,
                                   &args->m_tmp_arr);
            break;
        case RULE_ES:
            if (swap < 0)
                update_AGR_record(args->out_hdr, rec, tag->int_id, tag->coltype, swap, &args->tmp_arr,
                                  &args->m_tmp_arr);
            else
                flip_A_record(args->out_hdr, rec, tag->int_id, tag->coltype, swap);
            break;
        case RULE_AGR:
            update_AGR_record(args->out_hdr, rec, tag->int_id, tag->coltype, swap, &args->tmp_arr, &args->m_tmp_arr);
            break;
        default:
            error("Unexpected rule %d\n", tag->rule);
        }
    }

    return rec;
}

void destroy(void) {
    if (args->idx)
        fprintf(stderr, "Lines   total/swapped/reference added/rejected:\t%d/%d/%d/%d\n", args->ntotal, args->nswapped,
                args->nref_added, args->nrejected);
    free(args->af_arr);
    free(args->ploidy_arr);
    free(args->tmp_kstr.s);
    free(args->tmp_pad.s);
    int i;
    for (i = 0; i < args->m_tmp_als; i++) free(args->tmp_als[i].s);
    free(args->tmp_als);
    free(args->tmp_arr);
    free(args->int32_arr);
    if (args->reject_fh && hts_close(args->reject_fh) < 0) error("Close failed: %s\n", args->reject_fh->fn);
    free(args->tags);
    if (args->idx) regidx_destroy(args->idx);
    if (args->itr) regitr_destroy(args->itr);
    free(args->chains);
    free(args->blocks);
    fai_destroy(args->src_fai);
    fai_destroy(args->dst_fai);
    free(args);
}
