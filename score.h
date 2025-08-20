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

#include "tsv2vcf.h"

// see http://github.com/neurogenomics/MungeSumstats#future-enhancements
#define HDR_SNP 0
#define HDR_BP 1
#define HDR_CHR 2
#define HDR_A1 3
#define HDR_A2 4
#define HDR_P 5
#define HDR_Z 6
#define HDR_OR 7
#define HDR_BETA 8
#define HDR_N 9
#define HDR_N_CAS 10
#define HDR_N_CON 11
#define HDR_INFO 12
#define HDR_FRQ 13
#define HDR_A0 14
#define HDR_SE 15
#define HDR_LP 16
#define HDR_AC 17
#define HDR_NEFF 18
#define HDR_NEFFDIV2 19
#define HDR_HET_I2 20
#define HDR_HET_P 21
#define HDR_HET_LP 22
#define HDR_DIRE 23
#define HDR_SIZE 24
static const char *col_headers[HDR_SIZE] = {"SNP",  "BP", "CHR",   "A1",       "A2",     "P",     "Z",      "OR",
                                            "BETA", "N",  "N_CAS", "N_CON",    "INFO",   "FRQ",   "A*",     "SE",
                                            "LP",   "AC", "NEFF",  "NEFFDIV2", "HET_I2", "HET_P", "HET_LP", "DIRE"};

typedef struct {
    char *hdr_str;
    int hdr_num;
} mapping_t;

// see http://www.cog-genomics.org/plink/1.9/formats#assoc
static mapping_t plink_mapping[] = {{"SNP", HDR_SNP},   {"BP", HDR_BP},   {"CHR", HDR_CHR}, {"A1", HDR_A1},
                                    {"A2", HDR_A2},     {"P", HDR_P},     {"OR", HDR_OR},   {"BETA", HDR_BETA},
                                    {"INFO", HDR_INFO}, {"F_U", HDR_FRQ}, {"FRQ", HDR_FRQ}, {"SE", HDR_SE}};

// see http://www.cog-genomics.org/plink/2.0/formats#glm_logistic
static mapping_t plink2_mapping[] = {
    {"ID", HDR_SNP},      {"POS", HDR_BP},   {"CHROM", HDR_CHR},     {"A1", HDR_A1},     {"AX", HDR_A2},
    {"P", HDR_P},         {"Z_STAT", HDR_Z}, {"OR", HDR_OR},         {"BETA", HDR_BETA}, {"MACH_R2", HDR_INFO},
    {"A1_FREQ", HDR_FRQ}, {"SE", HDR_SE},    {"LOG(OR)_SE", HDR_SE}, {"LOG10_P", HDR_LP}};

// see print_header_output_single() in http://github.com/rgcgithub/regenie/blob/master/src/Step2_Models.cpp
static mapping_t regenie_mapping[] = {
    {"ID", HDR_SNP},     {"GENPOS", HDR_BP},     {"CHROM", HDR_CHR},        {"ALLELE1", HDR_A1},   {"BETA", HDR_BETA},
    {"N", HDR_N},        {"N_CASES", HDR_N_CAS}, {"N_CONTROLS", HDR_N_CON}, {"INFO", HDR_INFO},    {"A1FREQ", HDR_FRQ},
    {"ALLELE0", HDR_A0}, {"SE", HDR_SE},         {"LOG10P", HDR_LP},        {"AC_ALLELE1", HDR_AC}};

// see http://saigegit.github.io/SAIGE-doc/docs/single_step2.html#output-file
// see http://github.com/saigegit/SAIGE/blob/main/src/Main.cpp
static mapping_t saige_mapping[] = {{"SNPID", HDR_SNP},      {"markerID", HDR_SNP}, {"POS", HDR_BP},
                                    {"CHR", HDR_CHR},        {"Allele2", HDR_A1},   {"Allele1", HDR_A2},
                                    {"p.value", HDR_P},      {"BETA", HDR_BETA},    {"N", HDR_N},
                                    {"N_case", HDR_N_CAS},   {"N_ctrl", HDR_N_CON}, {"imputationInfo", HDR_INFO},
                                    {"AF_Allele2", HDR_FRQ}, {"SE", HDR_SE},        {"AC_Allele2", HDR_AC}};

// see http://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html#x1-490009.1
static mapping_t bolt_mapping[] = {{"SNP", HDR_SNP},      {"BP", HDR_BP},      {"CHR", HDR_CHR},
                                   {"ALLELE1", HDR_A1},   {"P_LINREG", HDR_P}, {"P_BOLT_LMM_INF", HDR_P},
                                   {"P_BOLT_LMM", HDR_P}, {"BETA", HDR_BETA},  {"A1FREQ", HDR_FRQ},
                                   {"ALLELE0", HDR_A0},   {"SE", HDR_SE}};

// static mapping_t genesis_mapping[] = {{"rsid", HDR_SNP}, {"pos", HDR_BP}, {"chr", HDR_CHR},
//                                       {"alt", HDR_A1}, {"P", HDR_P}, {"Z", HDR_Z},
//                                       {"Est", HDR_BETA}, {"n.obs", HDR_N}, {"eff.frq", HDR_FRQ},
//                                       {"ref", HDR_A0}, {"Est.SE", HDR_SE}};

// see Analyze() in http://github.com/statgen/METAL/blob/master/metal/Main.cpp
static mapping_t metal_mapping[] = {{"MarkerName", HDR_SNP}, {"Position", HDR_BP},   {"Chromosome", HDR_CHR},
                                    {"Allele1", HDR_A1},     {"Allele2", HDR_A2},    {"P-value", HDR_P},
                                    {"Zscore", HDR_Z},       {"Effect", HDR_BETA},   {"N", HDR_N},
                                    {"Freq1", HDR_FRQ},      {"StdErr", HDR_SE},     {"log(P)", HDR_LP},
                                    {"Weight", HDR_NEFF},    {"HetISq", HDR_HET_I2}, {"HetPVal", HDR_HET_P},
                                    {"logHetP", HDR_HET_LP}, {"Direction", HDR_DIRE}};

// see http://www.pgscatalog.org/downloads/
static const mapping_t pgs_mapping[] = {{"chr_name", HDR_CHR},
                                        {"chr_position", HDR_BP},
                                        {"rsID", HDR_SNP},
                                        {"effect_allele", HDR_A1},
                                        {"other_allele", HDR_A2},
                                        {"OR", HDR_OR},
                                        {"HR", HDR_OR},
                                        {"effect_weight", HDR_BETA},
                                        {"allelefrequency_effect", HDR_FRQ}};

// see http://doi.org/10.1101/2022.07.15.500230
// and http://www.ebi.ac.uk/gwas/docs/summary-statistics-format#format
// there is a mismatch between rsid and rs_id between the two definitions
// TODO support the encouraged ref_allele column
static const mapping_t ssf_mapping[] = {{"chromosome", HDR_CHR},   {"base_pair_location", HDR_BP},
                                        {"variant_id", HDR_SNP},   {"rsid", HDR_SNP},
                                        {"rs_id", HDR_SNP},        {"effect_allele", HDR_A1},
                                        {"other_allele", HDR_A2},  {"p_value", HDR_P},
                                        {"odds_ratio", HDR_OR},    {"hazard_ratio", HDR_OR},
                                        {"beta", HDR_BETA},        {"n", HDR_N},
                                        {"info", HDR_INFO},        {"effect_allele_frequency", HDR_FRQ},
                                        {"standard_error", HDR_SE}};

static inline mapping_t *mapping_preset_init(const char *preset, int *n) {
    if (strcasecmp(preset, "PLINK") == 0) {
        *n = sizeof(plink_mapping) / sizeof(mapping_t);
        return (mapping_t *)plink_mapping;
    } else if (strcasecmp(preset, "PLINK2") == 0) {
        *n = sizeof(plink2_mapping) / sizeof(mapping_t);
        return (mapping_t *)plink2_mapping;
    } else if (strcasecmp(preset, "REGENIE") == 0) {
        *n = sizeof(regenie_mapping) / sizeof(mapping_t);
        return (mapping_t *)regenie_mapping;
    } else if (strcasecmp(preset, "SAIGE") == 0) {
        *n = sizeof(saige_mapping) / sizeof(mapping_t);
        return (mapping_t *)saige_mapping;
    } else if (strcasecmp(preset, "BOLT") == 0) {
        *n = sizeof(bolt_mapping) / sizeof(mapping_t);
        return (mapping_t *)bolt_mapping;
    } else if (strcasecmp(preset, "METAL") == 0) {
        *n = sizeof(metal_mapping) / sizeof(mapping_t);
        return (mapping_t *)metal_mapping;
    } else if (strcasecmp(preset, "PGS") == 0) {
        *n = sizeof(pgs_mapping) / sizeof(mapping_t);
        return (mapping_t *)pgs_mapping;
    } else if (strcasecmp(preset, "SSF") == 0) {
        *n = sizeof(ssf_mapping) / sizeof(mapping_t);
        return (mapping_t *)ssf_mapping;
    } else {
        return NULL;
    }
}

static inline mapping_t *mapping_file_init(const char *fname, int *n) {
    mapping_t *mapping = NULL;
    int mapping_n = 0;
    int mapping_m = 0;
    kstring_t str = {0, 0, NULL};
    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fname, strerror(errno));
    while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
        char *ret = strchr(str.s, '\t');
        if (!ret) error("Could not parse line: %s\n", str.s);
        *ret = '\0';
        int i;
        for (i = 0; i < HDR_SIZE; i++) {
            if (strcmp(ret + 1, col_headers[i]) != 0) continue;
            hts_expand(mapping_t, mapping_n + 1, mapping_m, mapping);
            mapping[mapping_n].hdr_str = strdup(str.s);
            mapping[mapping_n].hdr_num = i;
            mapping_n++;
        }
    }
    free(str.s);
    hts_close(fp);
    *n = mapping_n;
    return mapping;
}

/****************************************
 * TSV FUNCTION                         *
 ****************************************/

// adapted from Petr Danecek's implementation of tsv_init() in
// bcftools/tsv2vcf.c
static inline tsv_t *tsv_init_delimiter(const char *str, char delimiter) {
    tsv_t *tsv = (tsv_t *)calloc(1, sizeof(tsv_t));
    kstring_t tmp = {0, 0, 0};
    const char *ss = str, *se = ss;
    tsv->ncols = 0;
    while (*ss) {
        if (delimiter == '\0')
            while (*se && !isspace(*se)) se++;
        else
            while (*se && *se != delimiter) se++;
        tsv->ncols++;
        tsv->cols = (tsv_col_t *)realloc(tsv->cols, sizeof(tsv_col_t) * tsv->ncols);
        tsv->cols[tsv->ncols - 1].name = NULL;
        tsv->cols[tsv->ncols - 1].setter = NULL;
        tmp.l = 0;
        kputsn(ss, se - ss, &tmp);
        if (strcasecmp("-", tmp.s)) tsv->cols[tsv->ncols - 1].name = strdup(tmp.s);
        if (!*se) break;
        se++;
        if (delimiter == '\0')
            while (*se && isspace(*se)) se++;
        ss = se;
    }
    free(tmp.s);
    return tsv;
}

// adapted from Petr Danecek's implementation of tsv_init() in
// bcftools/tsv2vcf.c
static int tsv_parse_delimiter(tsv_t *tsv, bcf1_t *rec, char *str, char delimiter) {
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while (*tsv->ss && tsv->icol < tsv->ncols) {
        if (delimiter == '\0')
            while (*tsv->se && !isspace(*tsv->se)) tsv->se++;
        else
            while (*tsv->se && *tsv->se != delimiter) tsv->se++;
        if (tsv->cols[tsv->icol].setter) {
            int ret = tsv->cols[tsv->icol].setter(tsv, rec, tsv->cols[tsv->icol].usr);
            if (ret < 0) return -1;
            status++;
        }
        if (*tsv->se) {
            tsv->se++;
            if (delimiter == '\0')
                while (*tsv->se && isspace(*tsv->se)) tsv->se++;
        }
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

static int tsv_setter_id_flexible(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char tmp = *tsv->se;
    *tsv->se = 0;
    if (strcasecmp("NA", tsv->ss))
        bcf_update_id((bcf_hdr_t *)usr, rec, tsv->ss);
    else
        bcf_update_id(NULL, rec, NULL);
    *tsv->se = tmp;
    return 0;
}

static int tsv_setter_pos_flexible(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char *endptr;
    rec->pos = strtol(tsv->ss, &endptr, 10) - 1;
    if (endptr != tsv->se) return -1;
    return 0;
}

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

static int tsv_setter_chrom_flexible(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char tmp = *tsv->se;
    *tsv->se = 0;
    rec->rid = bcf_hdr_name2id_flexible((bcf_hdr_t *)usr, tsv->ss);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_float_and_minus_log10(tsv_t *tsv, bcf1_t *rec, void *usr) {
    float *single = (float *)usr;
    char *ptr, *endptr;
    for (ptr = tsv->ss; *ptr != 'e' && *ptr != 'E' && ptr < tsv->se; ptr++);
    char tmp = *ptr;
    *ptr = '\0';
    *single = strtof(tsv->ss, &endptr);
    if (endptr != ptr) {
        if (*tsv->ss == '.' || (*tsv->ss == 'N' && *(tsv->ss + 1) == 'A'))
            *single = NAN;
        else
            return -1;
    } else {
        *single = 0.0 - log10f(*single);
        if (ptr != tsv->se) {
            float exponent = strtof(ptr + 1, &endptr);
            if (endptr != tsv->se) return -1;
            *single -= exponent;
        }
    }
    *ptr = tmp;
    return 0;
}

int tsv_read_float(tsv_t *tsv, bcf1_t *rec, void *usr) {
    float *single = (float *)usr;
    char *endptr;
    *single = (float)strtof(tsv->ss, &endptr);
    if (endptr != tsv->se) {
        if (*tsv->ss == '.' || (*tsv->ss == 'N' && *(tsv->ss + 1) == 'A'))
            *single = NAN;
        else
            return -1;
    }
    return 0;
}

int tsv_read_float_and_log(tsv_t *tsv, bcf1_t *rec, void *usr) {
    float *single = (float *)usr;
    char *endptr;
    *single = (float)strtof(tsv->ss, &endptr);
    if (endptr != tsv->se) {
        if (*tsv->ss == '.' || (*tsv->ss == 'N' && *(tsv->ss + 1) == 'A'))
            *single = NAN;
        else
            return -1;
    } else {
        *single = logf(*single);
    }
    return 0;
}
