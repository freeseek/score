/* The MIT License

   Copyright (C) 2021-2025 Giulio Genovese

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
#include <dirent.h>
#include <getopt.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include "bcftools.h"
#include "filter.h"
#include "score.h"

#define SCORE_VERSION "2025-08-19"

#define FLT_INCLUDE (1 << 0)
#define FLT_EXCLUDE (1 << 1)
#define Q_SCORE_THR (1 << 2)
#define TSV_MODE (1 << 3)
#define VARIANT_ID_MODE (1 << 4)

// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">
#define SCORE_GT 1
// ##FORMAT=<ID=DS,Number=A,Type=Float,Description="Genotype dosage">
#define SCORE_DS 2 // DS = AP1 + AP2
// ##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate
// Allele Dosage ">
#define SCORE_HDS 3
// ##FORMAT=<ID=AP1,Number=A,Type=Float,Description="ALT allele probability of
// first haplotype">
// ##FORMAT=<ID=AP2,Number=A,Type=Float,Description="ALT allele probability of
// second haplotype">
#define SCORE_AP 4
// ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype
// Probability">
#define SCORE_GP 5
// ##FORMAT=<ID=AS,Number=1,Type=Integer,Description="Allelic shift (1/-1 if the
// alternate allele is over/under represented)">
#define SCORE_AS 6

KHASH_MAP_INIT_INT64(64, int64_t)

static inline int is_missing(float f) { return isnan(f) || bcf_float_is_missing(f) || bcf_float_is_vector_end(f); }

/****************************************
 * HELPER FUNCTIONS                     *
 ****************************************/

static inline char **get_file_list(const char *pathname, int *nfiles) {
    char **filenames = NULL;
    DIR *d = opendir(pathname);
    if (d) { // check if d is a directory
        struct dirent *dir;
        *nfiles = 0;
        int mfiles = 0;
        int p = strlen(pathname);
        while ((dir = readdir(d))) {
            if (strcmp(dir->d_name, ".") == 0 || strcmp(dir->d_name, "..") == 0) continue;
            hts_expand0(char *, *nfiles + 1, mfiles, filenames);
            int q = strlen(dir->d_name);
            filenames[*nfiles] = (char *)malloc((p + q + 2) * sizeof(char));
            memcpy(filenames[*nfiles], pathname, p);
            filenames[*nfiles][p] = '/';
            memcpy(filenames[*nfiles] + p + 1, dir->d_name, q + 1);
            (*nfiles)++;
        }
        closedir(d);
    } else {
        filenames = hts_readlines(pathname, nfiles);
        if (!filenames) error("Failed to read from file %s\n", pathname);
    }
    if (*nfiles == 0) error("No files found in %s\n", pathname);
    return filenames;
}

/****************************************
 * ALLELES IMPLEMENTATION               *
 ****************************************/

typedef struct {
    int n;
    int m;
    char **str;
    void *str2int;
} alleles_t;

alleles_t *alleles;

static void alleles_destroy(alleles_t *alleles) {
    khash_str2int_destroy_free(alleles->str2int);
    free(alleles->str);
    free(alleles);
}

static int tsv_read_allele(tsv_t *tsv, bcf1_t *rec, void *usr) {
    int *idx = (int *)usr;
    if (tsv->se == tsv->ss) {
        *idx = -1;
    } else {
        char tmp = *tsv->se;
        *tsv->se = 0;
        char *s;
        for (s = tsv->ss; s < tsv->se; s++) *s = toupper((unsigned char)*s);
        if (khash_str2int_get(alleles->str2int, tsv->ss, idx) < 0) {
            hts_expand(char *, alleles->n + 1, alleles->m, alleles->str);
            alleles->str[alleles->n] = strdup(tsv->ss);
            khash_str2int_inc(alleles->str2int, alleles->str[alleles->n]);
            *idx = alleles->n;
            alleles->n++;
        }
        *tsv->se = tmp;
    }
    return 0;
}

static int (*tsv_setters[])(tsv_t *tsv, bcf1_t *rec, void *usr) = {tsv_setter_id_flexible,         // SNP
                                                                   tsv_setter_pos_flexible,        // BP
                                                                   tsv_setter_chrom_flexible,      // CHR
                                                                   tsv_read_allele,                // A1
                                                                   NULL,                           // A2
                                                                   tsv_read_float_and_minus_log10, // P
                                                                   NULL,                           // Z
                                                                   tsv_read_float_and_log,         // OR
                                                                   tsv_read_float,                 // BETA
                                                                   NULL,                           // N
                                                                   NULL,                           // N_CAS
                                                                   NULL,                           // N_CON
                                                                   NULL,                           // INFO
                                                                   NULL,                           // FRQ
                                                                   NULL,                           // A0
                                                                   NULL,                           // SE
                                                                   NULL,                           // LP
                                                                   NULL,                           // AC
                                                                   NULL,                           // NEFF
                                                                   NULL,                           // NEFFDIV2
                                                                   NULL,                           // NET_I2
                                                                   NULL,                           // HET_P
                                                                   NULL,                           // HET_LP
                                                                   NULL};                          // DIRE

/****************************************
 * PGS FILE IMPLEMENTATION              *
 ****************************************/

typedef struct {
    int a1_idx;
    float es;
    float lp;
} marker_t;

typedef struct {
    int use_snp;
    int snp_tsv;
    int chr_tsv;
    int bp_tsv;
    int a1_tsv;
    int beta_tsv;
    int or_tsv;
    int p_tsv;
    void *rid_pos2idx;
    void *id2idx;
    marker_t *markers;
    int n_markers;
    int m_markers;
    int all_markers;
} summary_t;

static summary_t *summary_init(const char *fn, bcf_hdr_t *hdr, mapping_t *mapping, int mapping_n, int flags) {
    summary_t *summary = (summary_t *)calloc(1, sizeof(summary_t));

    htsFile *fp = hts_open(fn, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fn, strerror(errno));

    kstring_t str = {0, 0, NULL};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", fn);
    while (str.s[0] == '#' && strncmp(str.s, "#CHR", 4) != 0 && strncmp(str.s, "#ID", 3) != 0)
        hts_getline(fp, KS_SEP_LINE, &str);

    // remove leading # if present
    if (str.s[0] == '#') {
        memmove(str.s, str.s + 1, str.l - 1);
        str.l--;
    }

    // some formats are tab-delimited and some formats (e.g. PLINK and SBayesR)
    // are not here we make a determination based on the first header row
    char delimiter = strchr(str.s, '\t') ? '\t' : '\0';
    marker_t *marker = (marker_t *)calloc(1, sizeof(marker_t));
    tsv_t *tsv = tsv_init_delimiter(str.s, delimiter);

    int i;
    for (i = 0; i < mapping_n; i++) {
        void *usr;
        if (!tsv_setters[mapping[i].hdr_num]) continue;
        switch (mapping[i].hdr_num) {
        case HDR_SNP:
        case HDR_CHR:
            usr = (void *)hdr;
            break;
        case HDR_A1:
            usr = (void *)&marker->a1_idx;
            break;
        case HDR_P:
        case HDR_LP:
            usr = (void *)&marker->lp;
            break;
        case HDR_OR:
        case HDR_BETA:
            usr = (void *)&marker->es;
            break;
        default:
            usr = NULL;
        }
        int ret = tsv_register(tsv, mapping[i].hdr_str, tsv_setters[mapping[i].hdr_num], usr);
        if (ret < 0) continue;
        switch (mapping[i].hdr_num) {
        case HDR_SNP:
            summary->snp_tsv = 1;
        case HDR_CHR:
            summary->chr_tsv = 1;
            break;
        case HDR_BP:
            summary->bp_tsv = 1;
            break;
        case HDR_A1:
            summary->a1_tsv = 1;
            break;
        case HDR_P:
        case HDR_LP:
            summary->p_tsv = 1;
            break;
        case HDR_OR:
        case HDR_BETA:
            summary->beta_tsv = 1;
            break;
        default:
            usr = NULL;
        }
    }

    if ((flags & VARIANT_ID_MODE) && !summary->snp_tsv)
        error(
            "Column for marker name is not provided in file %s but required with "
            "option --snp\n",
            fn);

    if ((flags & VARIANT_ID_MODE) || !summary->chr_tsv || !summary->bp_tsv) summary->use_snp = 1;
    if (summary->use_snp) {
        summary->id2idx = khash_str2int_init();
    } else {
        summary->rid_pos2idx = kh_init(64);
    }
    if (!summary->snp_tsv && summary->use_snp)
        error(
            "Columns for chromosome and position required if column for marker "
            "name is not provided in file %s\n",
            fn);

    if (!summary->a1_tsv) error("Column for effect allele required but missing from file %s\n", fn);

    if (!summary->beta_tsv)
        error(
            "Column for either effect weight or odds ratio required but "
            "missing from file %s\n",
            fn);

    if ((flags & Q_SCORE_THR) && !summary->p_tsv)
        fprintf(stderr,
                "Warning: p-value column needed with --q-score-thr missing from "
                "file %s\n",
                fn);

    int allele_warning = 0, match_warning = 0, duplicate_warning = 0;
    bcf1_t *rec = bcf_init();
    while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
        if (str.s[0] == '#') continue; // skip comments
        summary->all_markers++;
        rec->rid = -1;
        rec->pos = -1;
        if (tsv_parse_delimiter(tsv, rec, str.s, delimiter) < 0) error("Could not parse line: %s\n", str.s);
        if (marker->a1_idx < 0) {
            if (!allele_warning) {
                fprintf(stderr,
                        "Warning: could not recognize the effect allele in "
                        "line:\n%s\n(This warning is printed "
                        "only once.)\n",
                        str.s);
                allele_warning = 1;
            }
            continue;
        }
        if (summary->use_snp) {
            if (!*rec->d.id) {
                if (!match_warning) {
                    fprintf(stderr,
                            "Warning: missing variant id in line:\n%s\n(This warning "
                            "is printed only once.)\n",
                            str.s);
                    match_warning = 1;
                }
                continue;
            }
            char *key = strdup(rec->d.id);
            int size = khash_str2int_size(summary->id2idx);
            int ret = khash_str2int_inc(summary->id2idx, key);
            if (ret < 0) error("Unable to insert key %s in hash table\n", key);
            if (ret == size) {
                hts_expand(marker_t, summary->n_markers + 1, summary->m_markers, summary->markers);
                memcpy((void *)&summary->markers[summary->n_markers], (const void *)marker, sizeof(marker_t));
                summary->n_markers++;
            } else if (!duplicate_warning) {
                fprintf(stderr,
                        "Warning: could not include marker as variant id %s present "
                        "multiple times\n(This warning "
                        "is printed only once.)\n",
                        rec->d.id);
                free(key);
                duplicate_warning = 1;
            }
        } else {
            if (rec->rid < 0 || rec->pos < 0) {
                if (!match_warning) {
                    fprintf(stderr,
                            "Warning: could not recognize chromosome or position in "
                            "line:\n%s\n(This warning is "
                            "printed only once.)\n",
                            str.s);
                    match_warning = 1;
                }
                continue;
            }
            int ret;
            khash_t(64) *hash = (khash_t(64) *)summary->rid_pos2idx;
            khiter_t k = kh_put(64, hash, (((hts_pos_t)rec->rid) << 44) + rec->pos, &ret);
            if (ret < 0)
                error("Unable to insert key chromosome position %s %" PRId64 " in hash table\n",
                      bcf_hdr_id2name(hdr, rec->rid), rec->pos + 1);
            if (ret > 0) {
                kh_val(hash, k) = kh_size(hash) - 1;
                hts_expand(marker_t, summary->n_markers + 1, summary->m_markers, summary->markers);
                memcpy((void *)&summary->markers[summary->n_markers], (const void *)marker, sizeof(marker_t));
                summary->n_markers++;
            } else if (!duplicate_warning) {
                fprintf(stderr,
                        "Warning: could not include marker as chromosome position %s "
                        "%" PRId64 " present multiple times\n(This warning is printed only once.)\n",
                        bcf_hdr_id2name(hdr, rec->rid), rec->pos + 1);
                duplicate_warning = 1;
            }
        }
    }
    bcf_destroy(rec);
    tsv_destroy(tsv);
    free(marker);
    free(str.s);
    hts_close(fp);

    return summary;
}

static void summary_destroy(summary_t *summary) {
    if (summary->rid_pos2idx) kh_destroy(64, summary->rid_pos2idx);
    khash_str2int_destroy_free(summary->id2idx);
    free(summary->markers);
    free(summary);
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Compute polygenic scores from GWAS-VCF summary statistics.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Compute polygenic scores from GWAS-VCF summary statistics. (version " SCORE_VERSION
           " http://github.com/freeseek/score)\n"
           "\n"
           "Usage: bcftools +score [options] <in.vcf.gz> [<score1.gwas.vcf.gz> <score2.gwas.vcf.gz> ...]\n"
           "Plugin options:\n"
           "       --use <tag>               FORMAT tag to use to compute allele dosages: GP, AP, HDS, DS, GT, AS\n"
           "       --summaries <dir|file>    summary statistics files from directory or list from file\n"
           "       --q-score-thr LIST        comma separated list of p-value thresholds\n"
           "       --counts                  include SNP counts in the output table\n"
           "   -o, --output <file.tsv>       write output to a file [standard output]\n"
           "       --sample-header           output header for sample ID column [SAMPLE]\n"
           "   -e, --exclude <expr>          exclude sites for which the expression is true\n"
           "   -f, --apply-filters <list>    require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n"
           "   -i, --include <expr>          select sites for which the expression is true\n"
           "   -r, --regions <region>        restrict to comma-separated list of regions\n"
           "   -R, --regions-file <file>     restrict to regions listed in a file\n"
           "       --regions-overlap 0|1|2   Include if POS in the region (0), record overlaps (1), variant overlaps "
           "(2) [1]\n"
           "   -t, --targets [^]<region>     restrict to comma-separated list of regions. Exclude regions with \"^\" "
           "prefix\n"
           "   -T, --targets-file [^]<file>  restrict to regions listed in a file. Exclude regions with \"^\" prefix\n"
           "       --targets-overlap 0|1|2   Include if POS in the region (0), record overlaps (1), variant overlaps "
           "(2) [0]\n"
           "   -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" "
           "prefix)\n"
           "   -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n"
           "       --force-samples           only warn about unknown subset samples\n"
           "\n"
           "TSV Summary Statistics Options:\n"
           "   -c, --columns <preset>        column headers from preset "
           "(PLINK/PLINK2/REGENIE/SAIGE/BOLT/METAL/PGS/SSF)\n"
           "   -C, --columns-file <file>     column headers from tab-delimited file\n"
           "       --use-variant-id          use variant_id to match variants rather than chromosome and "
           "base_pair_location\n"
           "\n"
           "Examples:\n"
           "   bcftools +score --use DS -o scores.tsv input.bcf -c PLINK score.assoc\n"
           "   bcftools +score --use DS -o scores.tsv input.bcf -C colheaders.tsv "
           "PGC3_SCZ_wave3_public.clumped.v2.tsv.gz\n"
           "   bcftools +score --use GT -o scores.tsv --q-score-thr 1e-8,1e-7,1e-6,1e-5,1e-4,0.001,0.01,0.05 input.bcf "
           "-c GWAS-SSF PGS000001.txt.gz\n"
           "   bcftools +score --use DS -o scores.tsv -i 'INFO>0.8 && AF>0.01 && AF<0.99' input.bcf -c GWAS-SSF "
           "PGS000001.txt.gz PGS000002.txt.gz\n"
           "\n";
}

static double *parse_list(const char *str, int *n) {
    char *endptr;
    char **s = hts_readlist(str, 0, n);
    double *v = (double *)malloc(*n * sizeof(double));
    int i;
    for (i = 0; i < *n; i++) {
        v[i] = strtof(s[i], &endptr);
        if (*endptr) error("Could not parse element: %s\n", s[i]);
        free(s[i]);
    }
    free(s);
    return v;
}

int run(int argc, char **argv) {
    int i, j, k, idx, ap;
    int use_tag = 0;
    int display_cnts = 0;
    int filter_logic = 0;
    int regions_is_file = 0;
    int regions_overlap = 1;
    int targets_is_file = 0;
    int targets_overlap = 0;
    int sample_is_file = 0;
    int force_samples = 0;
    int flags = 0;
    const char *q_score_thr_str = NULL;
    const char *pathname = NULL;
    const char *output_fname = "-";
    const char *sample_header = "SAMPLE";
    const char *filter_str = NULL;
    const char *regions_list = NULL;
    const char *targets_list = NULL;
    const char *sample_names = NULL;
    const char *columns_preset = NULL;
    const char *columns_fname = NULL;

    filter_t *filter = NULL;
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);

    static struct option loptions[] = {{"use", required_argument, NULL, 1},
                                       {"summaries", required_argument, NULL, 2},
                                       {"q-score-thr", required_argument, NULL, 3},
                                       {"counts", no_argument, NULL, 4},
                                       {"output", required_argument, NULL, 'o'},
                                       {"sample-header", required_argument, NULL, 5},
                                       {"exclude", required_argument, NULL, 'e'},
                                       {"apply-filters", required_argument, NULL, 'f'},
                                       {"include", required_argument, NULL, 'i'},
                                       {"regions", required_argument, NULL, 'r'},
                                       {"regions-file", required_argument, NULL, 'R'},
                                       {"regions-overlap", required_argument, NULL, 6},
                                       {"targets", required_argument, NULL, 't'},
                                       {"targets-file", required_argument, NULL, 'T'},
                                       {"targets-overlap", required_argument, NULL, 7},
                                       {"samples", required_argument, NULL, 's'},
                                       {"samples-file", required_argument, NULL, 'S'},
                                       {"force-samples", no_argument, NULL, 8},
                                       {"columns", required_argument, NULL, 'c'},
                                       {"columns-file", required_argument, NULL, 'C'},
                                       {"use-variant-id", no_argument, NULL, 9},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?o:e:f:i:r:R:t:T:s:S:c:C:", loptions, NULL)) >= 0) {
        switch (c) {
        case 1:
            if (!strcasecmp(optarg, "GT"))
                use_tag = SCORE_GT;
            else if (!strcasecmp(optarg, "DS"))
                use_tag = SCORE_DS;
            else if (!strcasecmp(optarg, "HDS"))
                use_tag = SCORE_HDS;
            else if (!strcasecmp(optarg, "AP"))
                use_tag = SCORE_AP;
            else if (!strcasecmp(optarg, "GP"))
                use_tag = SCORE_GP;
            else if (!strcasecmp(optarg, "AS"))
                use_tag = SCORE_AS;
            else
                error(
                    "The argument not recognised, expected --use GT, DS, HDS, AP, "
                    "GP, or AS: %s\n",
                    optarg);
            break;
        case 2:
            pathname = optarg;
            break;
        case 3:
            flags |= Q_SCORE_THR;
            q_score_thr_str = optarg;
            break;
        case 4:
            display_cnts = 1;
            break;
        case 'o':
            output_fname = optarg;
            break;
        case 5:
            sample_header = optarg;
            break;
        case 'e':
            filter_str = optarg;
            filter_logic |= FLT_EXCLUDE;
            break;
        case 'f':
            sr->apply_filters = optarg;
            break;
        case 'i':
            filter_str = optarg;
            filter_logic |= FLT_INCLUDE;
            break;
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 6:
            if (!strcasecmp(optarg, "0"))
                regions_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                regions_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                regions_overlap = 2;
            else
                error("Could not parse: --regions-overlap %s\n", optarg);
            break;
        case 't':
            targets_list = optarg;
            break;
        case 'T':
            targets_list = optarg;
            targets_is_file = 1;
            break;
        case 7:
            if (!strcasecmp(optarg, "0"))
                targets_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                targets_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                targets_overlap = 2;
            else
                error("Could not parse: --targets-overlap %s\n", optarg);
            break;
        case 's':
            sample_names = optarg;
            break;
        case 'S':
            sample_names = optarg;
            sample_is_file = 1;
            break;
        case 8:
            force_samples = 1;
            break;
        case 'c':
            flags |= TSV_MODE;
            columns_preset = optarg;
            break;
        case 'C':
            flags |= TSV_MODE;
            columns_fname = optarg;
            break;
        case 9:
            flags |= VARIANT_ID_MODE;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }

    if ((pathname && optind + 1 != argc) || (!pathname && optind + 2 > argc)) error("%s", usage_text());

    if (columns_preset && columns_fname)
        error("Error: either --columns or --columns-file should be given, not both\n%s", usage_text());

    if (filter_logic == (FLT_EXCLUDE | FLT_INCLUDE)) error("Only one of --include or --exclude can be given.\n");
    if (regions_list) {
        bcf_sr_set_opt(sr, BCF_SR_REGIONS_OVERLAP, regions_overlap);
        if (bcf_sr_set_regions(sr, regions_list, regions_is_file) < 0)
            error("Failed to read the regions: %s\n", regions_list);
    }
    if (targets_list) {
        bcf_sr_set_opt(sr, BCF_SR_TARGETS_OVERLAP, targets_overlap);
        if (bcf_sr_set_targets(sr, targets_list, targets_is_file, 0) < 0)
            error("Failed to read the targets: %s\n", targets_list);
        sr->collapse |= COLLAPSE_BOTH;
    }

    if (!bcf_sr_add_reader(sr, argv[optind]))
        error("Error opening %s: %s\n", argv[optind], bcf_sr_strerror(sr->errnum));

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if (filter_str) filter = filter_init(hdr, filter_str);
    int n_q_score_thr = 1;
    double *q_score_thr = q_score_thr_str ? parse_list(q_score_thr_str, &n_q_score_thr) : NULL;
    if (q_score_thr)
        for (i = 0; i < n_q_score_thr; i++) q_score_thr[i] = -log10(q_score_thr[i]);

    // subset VCF file
    if (sample_names) {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if (ret < 0)
            error("Error parsing the sample list\n");
        else if (ret > 0) {
            if (force_samples)
                fprintf(stderr, "Warn: sample #%d not found in the header... skipping\n", ret);
            else
                error(
                    "Error: sample #%d not found in the header. Use "
                    "\"--force-samples\" to "
                    "ignore this error\n",
                    ret);
        }
        if (bcf_hdr_nsamples(hdr) == 0) error("Error: subsetting has removed all samples\n");
    }
    int n_smpls = bcf_hdr_nsamples(hdr);

    int n_files, n_prs, *prs2vcf = NULL, m_prs2vcf = 0, *prs2idx = NULL, m_prs2idx = 0;
    char **filenames = NULL;
    if (pathname) {
        filenames = get_file_list(pathname, &n_files);
    } else {
        n_files = argc - optind - 1;
        filenames = argv + optind + 1;
    }
    summary_t **summaries = NULL;
    char **prs_names = NULL;
    int m_prs_names = 0;
    if (!(flags & TSV_MODE)) {
        n_prs = 0;
        for (i = 0; i < n_files; i++) {
            if (!bcf_sr_add_reader(sr, filenames[i]))
                error("Error opening %s: %s\n", filenames[i], bcf_sr_strerror(sr->errnum));
            hdr = bcf_sr_get_header(sr, i + 1);
            hts_expand(int, n_prs + bcf_hdr_nsamples(hdr), m_prs2vcf, prs2vcf);
            hts_expand(int, n_prs + bcf_hdr_nsamples(hdr), m_prs2idx, prs2idx);
            hts_expand(char *, n_prs + bcf_hdr_nsamples(hdr), m_prs_names, prs_names);
            for (j = 0; j < bcf_hdr_nsamples(hdr); j++) {
                prs2vcf[n_prs + j] = i + 1;
                prs2idx[n_prs + j] = j;
                prs_names[n_prs + j] = strdup(hdr->samples[j]);
            }
            n_prs += bcf_hdr_nsamples(hdr);
            int es_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ES");
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, es_id))
                error(
                    "VCF summary statistics file %s does not include the ES FORMAT "
                    "field\n",
                    filenames[i]);
            if (q_score_thr) {
                int lp_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "LP");
                if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, lp_id))
                    error(
                        "VCF summary statistics file %s does not include the LP FORMAT "
                        "field\n",
                        filenames[i]);
            }
        }
    } else {
        n_prs = n_files;
        alleles = (alleles_t *)calloc(1, sizeof(alleles_t));
        alleles->str2int = khash_str2int_init();
        summaries = (summary_t **)malloc((n_prs) * sizeof(summary_t *));
        prs_names = (char **)malloc((n_prs) * sizeof(char *));
        hdr = bcf_sr_get_header(sr, 0);

        int mapping_n = 0;
        mapping_t *mapping = columns_preset ? mapping_preset_init(columns_preset, &mapping_n)
                                            : mapping_file_init(columns_fname, &mapping_n);
        if (columns_preset && !mapping)
            error("Error: preset not recognized with --columns %s\n%s", columns_preset, usage_text());
        for (i = 0; i < n_prs; i++) {
            summaries[i] = summary_init(filenames[i], hdr, mapping, mapping_n, flags);
            fprintf(stderr, "Loaded %d out of %d markers from file %s and matching by %s\n", summaries[i]->n_markers,
                    summaries[i]->all_markers, filenames[i],
                    summaries[i]->use_snp ? "marker name" : "chromosome position");
            char *ptr, *ext_str[] = {"gz", "txt", "tsv", "vcf", "bcf"};
            int j = 0;
            while (j < sizeof(ext_str) / sizeof(char *) && (ptr = strrchr(filenames[i], '.')))
                for (j = 0; j < sizeof(ext_str) / sizeof(char *); j++)
                    if (strcmp(ptr + 1, ext_str[j]) == 0) {
                        *ptr = '\0';
                        break;
                    }
            prs_names[i] = strdup(strrchr(filenames[i], '/') ? strrchr(filenames[i], '/') + 1 : filenames[i]);
        }
        if (columns_fname) {
            for (i = 0; i < mapping_n; i++) free(mapping[i].hdr_str);
            free(mapping);
        }
    }

    hdr = bcf_sr_get_header(sr, 0);
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    int ds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS");
    int hds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HDS");
    int ap1_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AP1");
    int ap2_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AP2");
    int gp_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GP");
    int as_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AS");

    if (!use_tag) {
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id)) use_tag = SCORE_GT;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, hds_id)) use_tag = SCORE_HDS;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap1_id) && bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap2_id))
            use_tag = SCORE_AP;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gp_id)) use_tag = SCORE_GP;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ds_id)) use_tag = SCORE_DS;
        if (!use_tag)
            error(
                "VCF file %s does not include any of the GT, DS, HDS, AP1/AP2, or "
                "DS FORMAT fields\n",
                argv[optind]);
    } else {
        switch (use_tag) {
        case SCORE_GT:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id))
                error("VCF file %s does not include the GT FORMAT field\n", argv[optind]);
            break;
        case SCORE_DS:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ds_id))
                error("VCF file %s does not include the DS FORMAT field\n", argv[optind]);
            break;
        case SCORE_HDS: // only for Minimac4 VCFs
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, hds_id))
                error("VCF file %s does not include the HDS FORMAT field\n", argv[optind]);
            break;
        case SCORE_AP:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap1_id) || !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap2_id))
                error(
                    "VCF file %s does not include either the AP1 or the AP2 FORMAT "
                    "fields\n",
                    argv[optind]);
            break;
        case SCORE_GP:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gp_id))
                error("VCF file %s does not include the GP FORMAT field\n", argv[optind]);
            break;
        case SCORE_AS:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, as_id))
                error("VCF file %s does not include the AS FORMAT field\n", argv[optind]);
            break;
        }
    }

    fprintf(stderr, "Using %s to compute polygenic scores\n",
            use_tag == SCORE_GT    ? "genotypes (GT)"
            : use_tag == SCORE_GP  ? "genotype probabilities (GP)"
            : use_tag == SCORE_AP  ? "ALT haplotype probabilities (AP)"
            : use_tag == SCORE_HDS ? "haploid alternate allele dosage (HDS)"
            : use_tag == SCORE_DS  ? "genotype dosages (DS)"
                                   : "allelic shifts (AS)");

    FILE *out_fh = strcmp("-", output_fname) ? fopen(output_fname, "w") : stdout;
    if (!out_fh) error("Error: cannot write to %s\n", output_fname);

    int m_int32 = 0, m_float = 0, n_float, m_aps = 2 * n_smpls;
    int32_t *int32_arr = NULL;
    float *float_arr = NULL;
    char *str = NULL;
    float *aps = (float *)malloc(m_aps * sizeof(float));
    int *n_matched = (int *)calloc(n_prs, sizeof(int));
    float *missing = (float *)malloc(n_smpls * sizeof(float));
    int *idxs = (flags & TSV_MODE) ? (int *)malloc(n_prs * sizeof(int)) : NULL;
    float *scores = (float *)calloc(n_prs * n_q_score_thr * n_smpls, sizeof(float));
    int *cnts = display_cnts ? (int *)calloc(n_prs * n_q_score_thr * n_smpls, sizeof(int)) : NULL;

    while (bcf_sr_next_line(sr)) {
        if (!bcf_sr_has_line(sr, 0)) continue;
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if (filter) {
            int ret = filter_test(filter, line, NULL);
            if ((filter_logic == FLT_INCLUDE && !ret) || ret) continue;
        }

        int skip_line = 1;
        for (i = 0; i < n_prs; i++) {
            if (!(flags & TSV_MODE)) {
                if (!bcf_sr_has_line(sr, prs2vcf[i])) continue;
                skip_line = 0;
            } else {
                idxs[i] = -1;
                if (summaries[i]->use_snp) {
                    if (khash_str2int_get(summaries[i]->id2idx, line->d.id, &idxs[i]) < 0) continue;
                } else {
                    khash_t(64) *hash = (khash_t(64) *)summaries[i]->rid_pos2idx;
                    khiter_t k = kh_get(64, hash, (((hts_pos_t)line->rid) << 44) + line->pos);
                    if (k == kh_end(hash)) continue;
                    idxs[i] = kh_val(hash, k);
                }
                skip_line = 0;
            }
        }
        if (skip_line) continue;

        hdr = bcf_sr_get_header(sr, 0);
        hts_expand(float, line->n_allele *n_smpls, m_aps, aps);
        memset((void *)aps, 0, line->n_allele * n_smpls * sizeof(float));
        memset((void *)missing, 0, n_smpls * sizeof(int));
        int number;
        char *ap_str[] = {"AP1", "AP2"};
        switch (use_tag) {
        case SCORE_GT:
            number = bcf_get_genotypes(hdr, line, &int32_arr, &m_int32);
            if (number <= 0) continue;
            number /= bcf_hdr_nsamples(hdr);
            assert(number == 2);
            for (k = 0; k < n_smpls; k++) {
                int32_t *ptr = int32_arr + (number * k);
                if (bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1])) {
                    missing[k] = 1;
                } else {
                    size_t allele;
                    if (ptr[0] != bcf_int32_vector_end) {
                        allele = bcf_gt_allele(ptr[0]);
                        if (allele < line->n_allele) aps[allele * n_smpls + k]++;
                    }
                    if (ptr[1] != bcf_int32_vector_end) {
                        allele = bcf_gt_allele(ptr[1]);
                        if (allele < line->n_allele) aps[allele * n_smpls + k]++;
                    }
                }
            }
            break;
        case SCORE_DS:
            number = bcf_get_format_float(hdr, line, "DS", &float_arr, &m_float);
            if (number <= 0) continue;
            number /= bcf_hdr_nsamples(hdr);
            assert(number == line->n_allele - 1);
            if (number == 1) { // line->n_allele == 2
                for (k = 0; k < n_smpls; k++) {
                    if (is_missing(float_arr[k])) {
                        missing[k] = 1;
                    } else {
                        aps[k] += 2.0f - float_arr[k];
                        aps[n_smpls + k] += float_arr[k];
                    }
                }
            } else {
                for (k = 0; k < n_smpls; k++) {
                    float *ptr = float_arr + (number * k);
                    aps[k] += 2.0f;
                    for (idx = 0; idx < number; idx++) {
                        if (is_missing(ptr[idx])) {
                            missing[k] = 1;
                        } else {
                            aps[k] -= ptr[idx];
                            aps[(idx + 1) * n_smpls + k] += ptr[idx];
                        }
                    }
                }
            }
            break;
        case SCORE_HDS: // only for Minimac4 VCFs
            number = bcf_get_format_float(hdr, line, "HDS", &float_arr, &m_float);
            if (number <= 0) continue;
            number /= bcf_hdr_nsamples(hdr);
            assert(number == 2 && line->n_allele == 2);
            for (k = 0; k < n_smpls; k++) {
                if (is_missing(float_arr[2 * k]) || is_missing(float_arr[2 * k + 1])) {
                    missing[k] = 1;
                } else {
                    aps[k] += 2.0f - float_arr[2 * k] - float_arr[2 * k + 1];
                    aps[n_smpls + k] += float_arr[2 * k] + float_arr[2 * k + 1];
                }
            }
            break;
        case SCORE_AP:
            for (ap = 0; ap < sizeof(ap_str) / sizeof(char *); ap++) {
                number = bcf_get_format_float(hdr, line, ap_str[ap], &float_arr, &m_float);
                if (number <= 0) continue;
                number /= bcf_hdr_nsamples(hdr);
                assert(number == line->n_allele - 1);
                if (number == 1) { // line->n_allele == 2
                    for (k = 0; k < n_smpls; k++) {
                        if (is_missing(float_arr[k])) {
                            missing[k] = 1;
                        } else {
                            aps[k] += 1.0f - float_arr[k];
                            aps[n_smpls + k] += float_arr[k];
                        }
                    }
                } else {
                    for (k = 0; k < n_smpls; k++) {
                        float *ptr = float_arr + (number * k);
                        aps[k] += 1.0f;
                        for (idx = 0; idx < number; idx++) {
                            if (is_missing(ptr[idx])) {
                                missing[k] = 1;
                            } else {
                                aps[k] -= ptr[idx];
                                aps[(idx + 1) * n_smpls + k] += ptr[idx];
                            }
                        }
                    }
                }
            }
            break;
        case SCORE_GP:
            number = bcf_get_format_float(hdr, line, "GP", &float_arr, &m_float);
            if (number <= 0) continue;
            number /= bcf_hdr_nsamples(hdr);
            assert(number == (line->n_allele) * (line->n_allele + 1) / 2);
            if (number == 3) { // line->n_allele == 2
                for (k = 0; k < n_smpls; k++) {
                    float *ptr = float_arr + (number * k);
                    if (is_missing(ptr[0]) || is_missing(ptr[1]) || is_missing(ptr[2])) {
                        missing[k] = 1;
                    } else {
                        aps[k] += 2.0f * ptr[0] + ptr[1];
                        aps[n_smpls + k] += ptr[1] + 2.0f * ptr[2];
                    }
                }
            } else {
                for (k = 0; k < n_smpls; k++) {
                    float *ptr = float_arr + (number * k);
                    // The Variant Call Format Specification
                    // for P=2 and N=2, the ordering is 00,01,11,02,12,22
                    // for P=2, the index of the genotype “a/b”, where a≤b, is b(b+ 1)/2
                    // +a
                    int a, b;
                    for (b = 0; b < line->n_allele; b++) {
                        for (a = 0; a <= b; a++) {
                            int idx = b * (b + 1) / 2 + a;
                            if (is_missing(ptr[idx])) {
                                missing[k] = 1;
                            } else {
                                aps[a * n_smpls + k] += ptr[idx];
                                aps[b * n_smpls + k] += ptr[idx];
                            }
                        }
                    }
                }
            }
            break;
        case SCORE_AS:
            number = bcf_get_format_int32(hdr, line, "AS", &int32_arr, &m_int32);
            if (number <= 0) continue;
            number /= bcf_hdr_nsamples(hdr);
            assert(number == 1 && line->n_allele == 2);
            for (k = 0; k < n_smpls; k++) {
                if (int32_arr[k] == 0 || int32_arr[k] == bcf_int32_missing) {
                    missing[k] = 1;
                } else {
                    aps[k] -= (float)int32_arr[k];
                    aps[n_smpls + k] += (float)int32_arr[k];
                }
            }
            break;
        }

        float es, lp = 0.0f;
        for (i = 0; i < n_prs; i++) {
            char *a1;
            if (!(flags & TSV_MODE)) {
                if (!bcf_sr_has_line(sr, prs2vcf[i])) continue;
                hdr = bcf_sr_get_header(sr, prs2vcf[i]);
                line = bcf_sr_get_line(sr, prs2vcf[i]);
                a1 = line->d.allele[1];
                n_float = bcf_get_format_float(hdr, line, "ES", &float_arr, &m_float);
                if (n_float <= 0) continue;
                if (n_float != bcf_hdr_nsamples(hdr))
                    error(
                        "VCF file %s has incorrect number of ES fields at position "
                        "%" PRId64 "\n",
                        filenames[prs2vcf[i] - 1], line->pos + 1);
                es = float_arr[prs2idx[i]];
                if (is_missing(es)) continue;
                if (q_score_thr) {
                    n_float = bcf_get_format_float(hdr, line, "LP", &float_arr, &m_float);
                    if (n_float <= 0) continue;
                    if (n_float != bcf_hdr_nsamples(hdr))
                        error(
                            "VCF file %s has incorrect number of LP fields at position "
                            "%" PRId64 "\n",
                            filenames[prs2vcf[i] - 1], line->pos + 1);
                    lp = float_arr[prs2idx[i]];
                    if (is_missing(lp)) continue;
                }
            } else {
                if (idxs[i] < 0) continue;
                marker_t *marker = &summaries[i]->markers[idxs[i]];
                a1 = alleles->str[marker->a1_idx];
                es = marker->es;
                lp = marker->lp;
            }
            // find effect allele
            int idx_allele;
            for (idx_allele = 0; idx_allele < line->n_allele; idx_allele++)
                if (strcmp(a1, line->d.allele[idx_allele]) == 0) break;
            if (idx_allele == line->n_allele) continue;
            n_matched[i]++;
            for (j = 0; j < n_q_score_thr; j++) {
                if (q_score_thr && lp < q_score_thr[j]) continue;
                float *ptr = scores + (i * n_q_score_thr + j) * n_smpls;
                float *ptr2 = aps + idx_allele * n_smpls;
                if (display_cnts) {
                    int *ptr3 = cnts + (i * n_q_score_thr + j) * n_smpls;
                    for (k = 0; k < n_smpls; k++) {
                        if (missing[k]) continue;
                        ptr[k] += es * ptr2[k];
                        ptr3[k]++;
                    }
                } else {
                    for (k = 0; k < n_smpls; k++) {
                        if (missing[k]) continue;
                        ptr[k] += es * ptr2[k];
                    }
                }
            }
        }
    }

    if (!(flags & TSV_MODE)) {
        for (i = 0; i < n_prs; i++)
            fprintf(stderr, "Matched %d markers for summary statistic %s\n", n_matched[i], prs_names[i]);
    } else {
        for (i = 0; i < n_prs; i++)
            fprintf(stderr, "Matched %d of %d loaded markers for summary statistic %s\n", n_matched[i],
                    summaries[i]->n_markers, prs_names[i]);
        alleles_destroy(alleles);
    }

    hdr = bcf_sr_get_header(sr, 0);
    fprintf(out_fh, "%s", sample_header);
    for (i = 0; i < n_prs; i++) {
        for (j = 0; j < n_q_score_thr; j++) {
            fprintf(out_fh, "\t%s", prs_names[i]);
            if (q_score_thr) fprintf(out_fh, "_p%.6g", exp(-M_LN10 * q_score_thr[j]));
            if (display_cnts) {
                fprintf(out_fh, "\t%s_CNT", prs_names[i]);
                if (q_score_thr) fprintf(out_fh, "_p%.6g", exp(-M_LN10 * q_score_thr[j]));
            }
        }
    }
    fprintf(out_fh, "\n");
    for (k = 0; k < n_smpls; k++) {
        fprintf(out_fh, "%s", hdr->samples[k]);
        for (i = 0; i < n_prs; i++)
            for (j = 0; j < n_q_score_thr; j++) {
                fprintf(out_fh, "\t%#.6g", scores[(i * n_q_score_thr + j) * n_smpls + k]);
                if (display_cnts) fprintf(out_fh, "\t%d", cnts[(i * n_q_score_thr + j) * n_smpls + k]);
            }
        fprintf(out_fh, "\n");
    }

    if (filter) filter_destroy(filter);
    if (pathname) {
        for (i = 0; i < n_files; i++) free(filenames[i]);
        free(filenames);
    }
    for (i = 0; i < n_prs; i++) free(prs_names[i]);
    free(prs_names);
    if (out_fh != stdout) fclose(out_fh);
    free(scores);
    free(cnts);
    free(idxs);
    free(aps);
    free(n_matched);
    free(missing);
    free(str);
    free(int32_arr);
    free(float_arr);
    free(q_score_thr);
    free(prs2vcf);
    free(prs2idx);
    if (summaries) {
        for (i = 0; i < n_prs; i++) summary_destroy(summaries[i]);
        free(summaries);
    }
    bcf_sr_destroy(sr);

    return 0;
}
