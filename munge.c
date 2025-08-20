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
#include <ctype.h>
#include <htslib/kseq.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "bcftools.h"
#include "score.h"

#define MUNGE_VERSION "2025-08-19"

#define IFFY_TAG "IFFY"
#define MISMATCH_TAG "REF_MISMATCH"

// http://github.com/MRCIEU/gwas-vcf-specification
#define NS 0
#define EZ 1
#define SI 2
#define NC 3
#define ES 4
#define SE 5
#define LP 6
#define AF 7
#define AC 8
#define NE 9
#define I2 10
#define CQ 11
#define SIZE 12
static const char *id_str[SIZE + 1] = {"NS", "EZ", "SI", "NC", "ES", "SE", "LP", "AF", "AC", "NE", "I2", "CQ", "ED"};
static const char *desc_str[SIZE + 1] = {
    "Variant-specific number of samples/individuals with called genotypes used to test association with specified "
    "trait",                                                                                 // NS
    "Z-score provided if it was used to derive the ES and SE fields",                        // EZ
    "Accuracy score of association statistics imputation",                                   // SI
    "Variant-specific number of cases used to estimate genetic effect (binary traits only)", // NC
    "Effect size estimate relative to the alternative allele",                               // ES
    "Standard error of effect size estimate",                                                // SE
    "-log10 p-value for effect estimate",                                                    // LP
    "Alternative allele frequency in trait subset",                                          // AF
    "Alternative allele count in the trait subset",                                          // AC
    "Variant-specific effective sample size",                                                // NE
    "Cochran's I^2 statistics",                                                              // I2
    "Cochran's Q -log10 p-value",                                                            // CQ
    "Effect size direction across studies"};                                                 // ED

int tsv_read_allele(tsv_t *tsv, bcf1_t *rec, void *usr) {
    kstring_t *str = (kstring_t *)usr;
    if (tsv->ss[0] == '<' || tsv->ss[0] == '*') {
        kputsn(tsv->ss, tsv->se - tsv->ss, str);
        return 0;
    }

    int symbolic = (tsv->ss[0] == 'D' || tsv->ss[0] == 'I' || tsv->ss[0] == 'd' || tsv->ss[0] == 'i');
    char *ptr;
    for (ptr = tsv->ss; ptr < tsv->se; ptr++)
        if (*ptr == '+') symbolic = 1;
    if (symbolic) {
        str->l = tsv->se - tsv->ss + 2;
        ks_resize(str, str->l + 1);
        str->s[0] = '<';
        memcpy(str->s + 1, tsv->ss, tsv->se - tsv->ss);
        str->s[str->l - 1] = '>';
    } else if (tsv->se - tsv->ss == 1 && ((tsv->ss[0] - '1') & 0xFC) == 0) {
        str->l = 1;
        ks_resize(str, 2);
        str->s[0] = "ACGT"[tsv->ss[0] - '1'];
        str->s[1] = '\0';
    } else {
        str->l = tsv->se - tsv->ss;
        ks_resize(str, str->l + 1);
        int i;
        for (i = 0; i < str->l; i++) str->s[i] = toupper(tsv->ss[i]);
    }
    str->s[str->l] = '\0';
    return 0;
}

int tsv_read_float_and_mult_two(tsv_t *tsv, bcf1_t *rec, void *usr) {
    float *single = (float *)usr;
    char *endptr;
    *single = (float)strtof(tsv->ss, &endptr);
    if (endptr != tsv->se) {
        if (*tsv->ss == '.' || (*tsv->ss == 'N' && *(tsv->ss + 1) == 'A'))
            *single = NAN;
        else
            return -1;
    }
    *single *= 2;
    return 0;
}

int tsv_read_string(tsv_t *tsv, bcf1_t *rec, void *usr) {
    kstring_t *str = (kstring_t *)usr;
    kputsn(tsv->ss, tsv->se - tsv->ss, str);
    return 0;
}

static int (*tsv_setters[])(tsv_t *tsv, bcf1_t *rec, void *usr) = {tsv_setter_id_flexible,         // SNP
                                                                   tsv_setter_pos_flexible,        // BP
                                                                   tsv_setter_chrom_flexible,      // CHR
                                                                   tsv_read_allele,                // A1
                                                                   tsv_read_allele,                // A2
                                                                   tsv_read_float_and_minus_log10, // P
                                                                   tsv_read_float,                 // Z
                                                                   tsv_read_float_and_log,         // OR
                                                                   tsv_read_float,                 // BETA
                                                                   tsv_read_float,                 // N
                                                                   tsv_read_float,                 // N_CAS
                                                                   tsv_read_float,                 // N_CON
                                                                   tsv_read_float,                 // INFO
                                                                   tsv_read_float,                 // FRQ
                                                                   tsv_read_allele,                // A0
                                                                   tsv_read_float,                 // SE
                                                                   tsv_read_float,                 // LP
                                                                   tsv_read_float,                 // AC
                                                                   tsv_read_float,                 // NEFF
                                                                   tsv_read_float_and_mult_two,    // NEFFDIV2
                                                                   tsv_read_float,                 // NET_I2
                                                                   tsv_read_float_and_minus_log10, // HET_P
                                                                   tsv_read_float,                 // HET_LP
                                                                   tsv_read_string};               // DIRE

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Convert summary statistics to GWAS-VCF.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Convert summary statistics to GWAS-VCF. "
           "(version " MUNGE_VERSION
           " http://github.com/freeseek/score)\n"
           "\n"
           "Usage: bcftools +munge [options] <score.gwas.ssf.tsv>\n"
           "Plugin options:\n"
           "   -c, --columns <preset>          column headers from preset "
           "(PLINK/PLINK2/REGENIE/SAIGE/BOLT/METAL/PGS/SSF)\n"
           "   -C, --columns-file <file>       column headers from tab-delimited file\n"
           "   -f, --fasta-ref <file>          reference sequence in fasta format\n"
           "       --fai <file>                reference sequence .fai index\n"
           "       --set-cache-size <int>      select fasta cache size in bytes\n"
           "       --iffy-tag <string>         FILTER annotation tag to record whether reference allele could not be "
           "determined [" IFFY_TAG
           "]\n"
           "       --mismatch-tag <string>     FILTER annotation tag to record whether reference does not match any "
           "allele [" MISMATCH_TAG
           "]\n"
           "   -s, --sample-name <string>      sample name for the phenotype [SAMPLE]\n"
           "       --ns <float>                number of samples\n"
           "       --nc <float>                number of cases\n"
           "       --ne <float>                effective sample size\n"
           "       --no-version                do not append version and command line to the header\n"
           "   -o, --output <file>             write output to a file [no output]\n"
           "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "       --threads <int>             use multithreading with INT worker threads [0]\n"
           "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
           "\n"
           "Examples:\n"
           "      bcftools +munge -c PLINK -f human_g1k_v37.fasta -Ob -o score.bcf score.assoc\n"
           "      bcftools +munge -C colheaders.tsv -f human_g1k_v37.fasta -s SCZ_2022 -Ob -o PGC3_SCZ.bcf "
           "PGC3_SCZ.tsv.gz\n"
           "\n";
}

int run(int argc, char **argv) {
    float ns = 0.0f;
    float nc = 0.0f;
    float ne = 0.0f;
    int i, idx;
    int cache_size = 0;
    int record_cmd_line = 1;
    int write_index = 0;
    int output_type = FT_VCF;
    int clevel = -1;
    int n_threads = 0;
    char *tmp = NULL;
    const char *columns_preset = NULL;
    const char *columns_fname = NULL;
    const char *ref_fname = NULL;
    const char *fai_fname = NULL;
    const char *iffy_tag = IFFY_TAG;
    const char *mismatch_tag = MISMATCH_TAG;
    const char *sample = "SAMPLE";
    const char *output_fname = "-";
    char *index_fname;
    faidx_t *fai;
    htsFile *out_fh = NULL;

    static struct option loptions[] = {{"columns", required_argument, NULL, 'c'},
                                       {"columns-file", required_argument, NULL, 'C'},
                                       {"fasta-ref", required_argument, NULL, 'f'},
                                       {"fai", required_argument, NULL, 1},
                                       {"set-cache-size", required_argument, NULL, 2},
                                       {"iffy-tag", required_argument, NULL, 3},
                                       {"mismatch-tag", required_argument, NULL, 4},
                                       {"sample-name", required_argument, NULL, 's'},
                                       {"ns", required_argument, NULL, 5},
                                       {"nc", required_argument, NULL, 6},
                                       {"ne", required_argument, NULL, 7},
                                       {"no-version", no_argument, NULL, 8},
                                       {"output", required_argument, NULL, 'o'},
                                       {"output-type", required_argument, NULL, 'O'},
                                       {"threads", required_argument, NULL, 9},
                                       {"write-index", optional_argument, NULL, 'W'},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?c:C:f:s:o:O:W::", loptions, NULL)) >= 0) {
        switch (c) {
        case 'c':
            columns_preset = optarg;
            break;
        case 'C':
            columns_fname = optarg;
            break;
        case 'f':
            ref_fname = optarg;
            break;
        case 1:
            fai_fname = optarg;
            break;
        case 2:
            cache_size = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --set-cache-size %s\n", optarg);
            break;
        case 3:
            iffy_tag = optarg;
            break;
        case 4:
            mismatch_tag = optarg;
            break;
        case 's':
            sample = optarg;
            break;
        case 5:
            ns = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --ns %s\n", optarg);
            break;
        case 6:
            nc = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --nc %s\n", optarg);
            break;
        case 7:
            ne = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --ne %s\n", optarg);
            break;
        case 8:
            record_cmd_line = 0;
            break;
        case 'o':
            output_fname = optarg;
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
            n_threads = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --threads %s\n", optarg);
            break;
        case 'W':
            if (!(write_index = write_index_parse(optarg))) error("Unsupported index format '%s'\n", optarg);
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }

    char *input_fname = NULL;
    if (optind == argc) {
        if (!isatty(fileno((FILE *)stdin))) {
            input_fname = "-"; // reading from stdin
        } else {
            error("%s", usage_text());
        }
    } else if (optind + 1 != argc) {
        error("%s", usage_text());
    } else {
        input_fname = argv[optind];
    }

    char wmode[8];
    set_wmode(wmode, output_type, (char *)output_fname, clevel);
    out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
    if (out_fh == NULL) error("Error: cannot write to \"%s\": %s\n", output_fname, strerror(errno));
    if (n_threads) hts_set_threads(out_fh, n_threads);
    if (!ref_fname && !fai_fname) error("Expected the -f or --fai option\n");
    fai = fai_load3(ref_fname ? ref_fname : fai_fname, fai_fname, NULL, FAI_CREATE);
    if (!fai) error("Could not load the reference %s\n", ref_fname);
    if (cache_size) fai_set_cache_size(fai, cache_size);
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    int n = faidx_nseq(fai);
    for (i = 0; i < n; i++) {
        const char *seq = faidx_iseq(fai, i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq, len);
    }
    if (bcf_hdr_printf(hdr, "##FILTER=<ID=%s,Description=\"Reference allele could not be determined\">", iffy_tag) < 0)
        error_errno("Failed to add \"%s\" FILTER header", iffy_tag);
    int iffy_id = bcf_hdr_id2int(hdr, BCF_DT_ID, iffy_tag);
    if (bcf_hdr_printf(hdr, "##FILTER=<ID=%s,Description=\"Reference does not match any allele\">", mismatch_tag) < 0)
        error_errno("Failed to add \"%s\" FILTER header", mismatch_tag);
    int mismatch_id = bcf_hdr_id2int(hdr, BCF_DT_ID, mismatch_tag);

    if ((!columns_preset && !columns_fname) || (columns_preset && columns_fname))
        error("Error: one of --columns or --columns-file should be given, not both\n%s", usage_text());
    int mapping_n = 0;
    mapping_t *mapping =
        columns_preset ? mapping_preset_init(columns_preset, &mapping_n) : mapping_file_init(columns_fname, &mapping_n);
    if (columns_preset && !mapping)
        error("Error: preset not recognized with --columns %s\n%s", columns_preset, usage_text());

    kstring_t str = {0, 0, NULL};
    htsFile *fp = hts_open(input_fname, "r");
    if (fp == NULL) error("Could not open %s: %s\n", input_fname, strerror(errno));
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", input_fname);
    while (str.s[0] == '#' && strncmp(str.s, "#CHR", 4) != 0 && strncmp(str.s, "#ID", 3) != 0)
        hts_getline(fp, KS_SEP_LINE, &str);

    // remove leading # if present
    if (str.s[0] == '#') {
        memmove(str.s, str.s + 1, str.l - 1);
        str.l--;
    }

    // remove _[0-9]*$ from the header to make sure FRQ_U is recognized
    char *ptr = strchr(str.s, '_');
    while (ptr) {
        char *ptr2 = ptr + 1;
        while (isdigit(*ptr2)) ptr2++;
        if (ptr2 > ptr + 1 && (isspace(*ptr2) || *ptr2 == '\0')) {
            memmove(ptr, ptr2, str.l - (ptr2 - str.s));
            str.l -= ptr2 - ptr;
            str.s[str.l] = '\0';
        }
        ptr = strchr(ptr + 1, '_');
    }

    // some formats are tab-delimited, some are comma-separated, and some formats (e.g. PLINK and SBayesR) are not here
    // we make a determination based on the first header row
    char delimiter = strchr(str.s, '\t') ? '\t' : strchr(str.s, ',') ? ',' : '\0';
    kstring_t alleles[2] = {{0, 0, NULL}, {0, 0, NULL}};
    kstring_t esd_str = {0, 0, NULL};
    tsv_t *tsv = tsv_init_delimiter(str.s, delimiter);

    int chr = 0;
    int pos = 0;
    int alt = 0;
    int ref = 0;
    int output[SIZE] = {0};
    int output_nco = 0;
    int output_esd = 0;
    float val_nco;
    float val[SIZE];
    for (i = 0; i < mapping_n; i++) {
        void *usr;
        switch (mapping[i].hdr_num) {
        case HDR_SNP:
        case HDR_CHR:
            usr = (void *)hdr;
            break;
        case HDR_A1:
            usr = (void *)&alleles[1];
            break;
        case HDR_A2:
        case HDR_A0:
            usr = (void *)&alleles[0];
            break;
        case HDR_P:
        case HDR_LP:
            usr = (void *)&val[LP];
            break;
        case HDR_Z:
            usr = (void *)&val[EZ];
            break;
        case HDR_OR:
        case HDR_BETA:
            usr = (void *)&val[ES];
            break;
        case HDR_N:
            usr = (void *)&val[NS];
            break;
        case HDR_N_CAS:
            usr = (void *)&val[NC];
            break;
        case HDR_N_CON:
            usr = (void *)&val_nco;
            break;
        case HDR_INFO:
            usr = (void *)&val[SI];
            break;
        case HDR_FRQ:
            usr = (void *)&val[AF];
            break;
        case HDR_SE:
            usr = (void *)&val[SE];
            break;
        case HDR_AC:
            usr = (void *)&val[AC];
            break;
        case HDR_NEFF:
        case HDR_NEFFDIV2:
            usr = (void *)&val[NE];
            break;
        case HDR_HET_I2:
            usr = (void *)&val[I2];
            break;
        case HDR_HET_P:
        case HDR_HET_LP:
            usr = (void *)&val[CQ];
            break;
        case HDR_DIRE:
            usr = (void *)&esd_str;
            break;
        default:
            usr = NULL;
        }
        int ret = tsv_register(tsv, mapping[i].hdr_str, tsv_setters[mapping[i].hdr_num], usr);
        if (ret < 0) continue;
        switch (mapping[i].hdr_num) {
        case HDR_CHR:
            chr = 1;
            break;
        case HDR_BP:
            pos = 1;
            break;
        case HDR_A1:
            alt = 1;
            break;
        case HDR_A2:
        case HDR_A0:
            ref = 1;
            break;
        case HDR_P:
        case HDR_LP:
            output[LP] = 1;
            break;
        case HDR_Z:
            output[EZ] = 1;
            break;
        case HDR_OR:
        case HDR_BETA:
            output[ES] = 1;
            break;
        case HDR_N:
            output[NS] = 1;
            break;
        case HDR_N_CAS:
            output[NC] = 1;
            break;
        case HDR_N_CON:
            output_nco = 1;
            break;
        case HDR_INFO:
            output[SI] = 1;
            break;
        case HDR_FRQ:
            output[AF] = 1;
            break;
        case HDR_SE:
            output[SE] = 1;
            break;
        case HDR_AC:
            output[AC] = 1;
            break;
        case HDR_NEFF:
        case HDR_NEFFDIV2:
            output[NE] = 1;
            break;
        case HDR_HET_I2:
            output[I2] = 1;
            break;
        case HDR_HET_P:
        case HDR_HET_LP:
            output[CQ] = 1;
            break;
        case HDR_DIRE:
            output_esd = 1;
            break;
        default:
            usr = NULL;
        }
        if (ns) output[NS] = 1;
        if (nc) output[NC] = 1;
        if (output[NC] && output_nco) output[NS] = 1;
    }
    if (!chr) error("Could not find chromosome column in input file\n");
    if (!pos) error("Could not find position column in input file\n");
    if (!ref) error("Could not find reference allele column in input file\n");
    if (!alt) error("Could not find alternate allele column in input file\n");
    if (!output[ES]) fprintf(stderr, "Warning: could not find column to compute beta in input file\n");
    if (!output[SE]) fprintf(stderr, "Warning: could not find standard error column in input file\n");
    if (!output[LP]) fprintf(stderr, "Warning: could not find column to compute -log10 p-value in input file\n");
    for (idx = 0; idx < SIZE; idx++)
        if (output[idx]
            && bcf_hdr_printf(hdr, "##FORMAT=<ID=%s,Number=A,Type=Float,Description=\"%s\">", id_str[idx],
                              desc_str[idx])
                   < 0)
            error_errno("Failed to add \"%s\" FORMAT header", id_str[idx]);
    if (output_esd
        && bcf_hdr_printf(hdr, "##FORMAT=<ID=%s,Number=A,Type=String,Description=\"%s\">", id_str[SIZE], desc_str[SIZE])
               < 0)
        error_errno("Failed to add \"%s\" FORMAT header", id_str[SIZE]);
    if (record_cmd_line) bcf_hdr_append_version(hdr, argc, argv, "bcftools_munge");
    for (idx = 0; idx < SIZE; idx++) {
        if (output[idx]) {
            bcf_hdr_add_sample(hdr, sample);
            break;
        }
    }
    if (bcf_hdr_write(out_fh, hdr) < 0) error("Unable to write to output VCF file\n");
    if (init_index2(out_fh, hdr, output_fname, &index_fname, write_index) < 0)
        error("Error: failed to initialise index for %s\n", output_fname);

    bcf1_t *rec = bcf_init();
    bcf_update_id(NULL, rec, NULL);
    bcf_float_set_missing(rec->qual);
    while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
        if (str.s[0] == '#') continue; // skip comments
        rec->rid = -1;
        rec->pos = -1;
        alleles[0].l = 0;
        alleles[1].l = 0;
        bcf_update_filter(hdr, rec, NULL, 0);
        for (idx = 0; idx < SIZE; idx++) val[idx] = NAN;
        esd_str.l = 0;
        if (ns) val[NS] = ns;
        if (nc) val[NC] = nc;
        if (ne) val[NE] = ne;
        if (tsv_parse_delimiter(tsv, rec, str.s, delimiter) < 0) error("Could not parse line: %s\n", str.s);
        if (rec->rid < 0) {
            fprintf(stderr, "Warning: could not convert record\n%s\n", str.s);
            continue;
        }
        if (output[NC] && output_nco) val[NS] = val[NC] + val_nco;
        int swap = 0;
        if (ref_fname) { // swap ref and alt alleles if necessary
            int len;
            char *ref =
                faidx_fetch_seq(fai, bcf_hdr_id2name(hdr, rec->rid), rec->pos,
                                rec->pos + (alleles[0].l > alleles[1].l ? alleles[0].l : alleles[1].l) - 1, &len);
            if (!ref || len < 1)
                error("faidx_fetch_seq failed at %s:%" PRId64 " (are you using the correct reference genome?)\n",
                      bcf_seqname(hdr, rec), rec->pos + 1);
            int ref_match = strncasecmp(ref, alleles[0].s, alleles[0].l) == 0;
            int alt_match = strncasecmp(ref, alleles[1].s, alleles[1].l) == 0;
            if (!ref_match && alt_match) {
                swap = 1;
                val[EZ] = -val[EZ];
                val[ES] = -val[ES];
                val[AF] = 1.0f - val[AF];
                val[AC] = 2.0f * val[NS] - val[AC];
            } else if (ref_match && alt_match) {
                bcf_update_filter(hdr, rec, &iffy_id, 1);
            } else if (!ref_match && !alt_match) {
                bcf_update_filter(hdr, rec, &mismatch_id, 1);
            }
            free(ref);
        }
        if (swap) {
            kputc(',', &alleles[1]);
            kputs(alleles[0].s, &alleles[1]);
            bcf_update_alleles_str(hdr, rec, alleles[1].s);
            if (output_esd) {
                char *ptr;
                for (ptr = esd_str.s; ptr < esd_str.s + esd_str.l; ptr++) {
                    if (*ptr == '+')
                        *ptr = '-';
                    else if (*ptr == '-')
                        *ptr = '+';
                }
            }
        } else {
            kputc(',', &alleles[0]);
            kputs(alleles[1].s, &alleles[0]);
            bcf_update_alleles_str(hdr, rec, alleles[0].s);
        }
        for (idx = 0; idx < SIZE; idx++) {
            if (output[idx]) {
                if (isnan(val[idx])) bcf_float_set_missing(val[idx]);
                bcf_update_format_float(hdr, rec, id_str[idx], &val[idx], 1);
            }
        }
        if (output_esd) bcf_update_format_char(hdr, rec, id_str[SIZE], esd_str.s, esd_str.l);
        if (bcf_write(out_fh, hdr, rec) < 0) error("Unable to write to output VCF file\n");
    }

    hts_close(fp);
    bcf_destroy(rec);
    tsv_destroy(tsv);
    free(alleles[0].s);
    free(alleles[1].s);
    free(esd_str.s);
    free(str.s);
    if (columns_fname) {
        for (i = 0; i < mapping_n; i++) free(mapping[i].hdr_str);
        free(mapping);
    }
    bcf_hdr_destroy(hdr);
    fai_destroy(fai);
    if (write_index) {
        if (bcf_idx_save(out_fh) < 0) {
            if (hts_close(out_fh) != 0) error("Close failed %s\n", strcmp(output_fname, "-") ? output_fname : "stdout");
            error("Error: cannot write to index %s\n", index_fname);
        }
        free(index_fname);
    }
    if (hts_close(out_fh) < 0) error("Close failed: %s\n", out_fh->fn);
    return 0;
}
