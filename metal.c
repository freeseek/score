/* The MIT License

   Copyright (C) 2022-2024 Giulio Genovese

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
#include <getopt.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include "bcftools.h"
#include "filter.h"

#define METAL_VERSION "2024-05-05"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

// https://github.com/MRCIEU/gwas-vcf-specification
#define NS 0
#define EZ 1
#define NC 2
#define ES 3
#define SE 4
#define LP 5
#define AF 6
#define AC 7
#define NE 8
#define I2 9
#define CQ 10
#define SIZE 11
static const char *id_str[SIZE + 1] = {"NS", "EZ", "NC", "ES", "SE", "LP", "AF", "AC", "NE", "I2", "CQ", "ED"};
static const char *desc_str[SIZE + 1] = {
    "Variant-specific number of samples/individuals with called genotypes used to test association with specified "
    "trait",                                                                                 // NS
    "Z-score provided if it was used to derive the ES and SE fields",                        // EZ
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

/****************************************
 * FUNCTION TO COMPUTE LOG P-VALUES     *
 ****************************************/

// Cody, W. J. Algorithm 715: SPECFUN–a portable FORTRAN package of special function routines and test drivers. ACM
// Trans. Math. Softw. 19, 22–30 (1993). https://doi.org/10.1145/151271.151273 ANORM function see pnorm_both() in
// https://github.com/wch/r-source/blob/trunk/src/nmath/pnorm.c see logndist() in
// https://github.com/statgen/METAL/blob/master/libsrc/MathStats.cpp
#define M_SQRT_32 5.656854249492380195206754896838    /* sqrt(32) */
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /* 1/sqrt(2pi) */
// this function is equivalent to pnorm(-z, log.p = TRUE) but with lower precision for -37.5193 < z < -0.67448975
static double log_ndist(double z) {
    const double a0 = 2.2352520354606839287E0;
    const double a1 = 1.6102823106855587881E2;
    const double a2 = 1.0676894854603709582E3;
    const double a3 = 1.8154981253343561249E4;
    const double a4 = 6.5682337918207449113E-2;
    const double b0 = 4.720258190468824187E1;
    const double b1 = 9.7609855173777669322E2;
    const double b2 = 1.0260932208618978205E4;
    const double b3 = 4.5507789335026729956E4;
    const double c0 = 3.9894151208813466764E-1;
    const double c1 = 8.8831497943883759412E0;
    const double c2 = 9.3506656132177855979E1;
    const double c3 = 5.9727027639480026226E2;
    const double c4 = 2.4945375852903726711E3;
    const double c5 = 6.8481904505362823326E3;
    const double c6 = 1.1602651437647350124E4;
    const double c7 = 9.8427148383839780218E3;
    const double c8 = 1.0765576773720192317E-8;
    const double d0 = 2.2266688044328115691E1;
    const double d1 = 2.3538790178262499861E2;
    const double d2 = 1.519377599407554805E3;
    const double d3 = 6.485558298266760755E3;
    const double d4 = 1.8615571640885098091E4;
    const double d5 = 3.4900952721145977266E4;
    const double d6 = 3.8912003286093271411E4;
    const double d7 = 1.9685429676859990727E4;
    const double p0 = 2.1589853405795699E-1;
    const double p1 = 1.274011611602473639E-1;
    const double p2 = 2.2235277870649807E-2;
    const double p3 = 1.421619193227893466E-3;
    const double p4 = 2.9112874951168792E-5;
    const double p5 = 2.307344176494017303E-2;
    const double q0 = 1.28426009614491121E0;
    const double q1 = 4.68238212480865118E-1;
    const double q2 = 6.59881378689285515E-2;
    const double q3 = 3.78239633202758244E-3;
    const double q4 = 7.29751555083966205E-5;
    double y, xsq, xnum, xden, temp, del;
    y = fabs(z);
    if (z < -0.67448975) {
        return log(ldexp(kf_erfc(z / M_SQRT2), -1));
    } else if (z <= 0.67448975) {
        xsq = z * z;
        xnum = (((a4 * xsq + a0) * xsq + a1) * xsq + a2) * xsq + a3;
        xden = (((xsq + b0) * xsq + b1) * xsq + b2) * xsq + b3;
        temp = z * xnum / xden;
        return log(0.5 - temp);
    } else if (z <= M_SQRT_32) {
        xnum = (((((((c8 * y + c0) * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7;
        xden = (((((((y + d0) * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7;
        temp = xnum / xden;
    } else {
        xsq = 1.0 / (z * z);
        xnum = ((((p5 * xsq + p0) * xsq + p1) * xsq + p2) * xsq + p3) * xsq;
        xden = ((((xsq + q0) * xsq + q1) * xsq + q2) * xsq + q3) * xsq;
        temp = xsq * (xnum + p4) / (xden + q4);
        temp = (M_1_SQRT_2PI - temp) / y;
    }
    xsq = ldexp(trunc(ldexp(y, 4)), -4);
    del = (y - xsq) * (y + xsq);
    return (-xsq * ldexp(xsq, -1)) - ldexp(del, -1) + log(temp);
}

#define M_2PI 6.283185307179586476925286766559 /* 2*pi */

// Wichura, M. J. Algorithm AS 241: The Percentage Points of the Normal Distribution. Applied Statistics 37, 477 (1988).
// https://doi.org/10.2307/2347330 PPND16 function (algorithm AS241) http://lib.stat.cmu.edu/apstat/241
// see qnorm5() in https://github.com/wch/r-source/blob/trunk/src/nmath/qnorm.c
// see ninv() in https://github.com/statgen/METAL/blob/master/libsrc/MathStats.cpp
// this function is equivalent to qnorm(log_p, log.p = TRUE)
static double inv_log_ndist(double log_p) {
    const double a0 = 3.3871328727963666080E0;
    const double a1 = 1.3314166789178437745E2;
    const double a2 = 1.9715909503065514427E3;
    const double a3 = 1.3731693765509461125E4;
    const double a4 = 4.5921953931549871457E4;
    const double a5 = 6.7265770927008700853E4;
    const double a6 = 3.3430575583588128105E4;
    const double a7 = 2.5090809287301226727E3;
    const double b1 = 4.2313330701600911252E1;
    const double b2 = 6.8718700749205790830E2;
    const double b3 = 5.3941960214247511077E3;
    const double b4 = 2.1213794301586595867E4;
    const double b5 = 3.9307895800092710610E4;
    const double b6 = 2.8729085735721942674E4;
    const double b7 = 5.2264952788528545610E3;
    const double c0 = 1.42343711074968357734E0;
    const double c1 = 4.63033784615654529590E0;
    const double c2 = 5.76949722146069140550E0;
    const double c3 = 3.64784832476320460504E0;
    const double c4 = 1.27045825245236838258E0;
    const double c5 = 2.41780725177450611770E-1;
    const double c6 = 2.27238449892691845833E-2;
    const double c7 = 7.74545014278341407640E-4;
    const double d1 = 2.05319162663775882187E0;
    const double d2 = 1.67638483018380384940E0;
    const double d3 = 6.89767334985100004550E-1;
    const double d4 = 1.48103976427480074590E-1;
    const double d5 = 1.51986665636164571966E-2;
    const double d6 = 5.47593808499534494600E-4;
    const double d7 = 1.05075007164441684324E-9;
    const double e0 = 6.65790464350110377720E0;
    const double e1 = 5.46378491116411436990E0;
    const double e2 = 1.78482653991729133580E0;
    const double e3 = 2.96560571828504891230E-1;
    const double e4 = 2.65321895265761230930E-2;
    const double e5 = 1.24266094738807843860E-3;
    const double e6 = 2.71155556874348757815E-5;
    const double e7 = 2.01033439929228813265E-7;
    const double f1 = 5.99832206555887937690E-1;
    const double f2 = 1.36929880922735805310E-1;
    const double f3 = 1.48753612908506148525E-2;
    const double f4 = 7.86869131145613259100E-4;
    const double f5 = 1.84631831751005468180E-5;
    const double f6 = 1.42151175831644588870E-7;
    const double f7 = 2.04426310338993978564E-15;
    double p, q, r, x;
    p = exp(log_p);
    q = p - 0.5;
    if (fabs(q) <= 0.425) {
        r = 0.180625 - q * q;
        return q * (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3) * r + a2) * r + a1) * r + a0)
               / (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3) * r + b2) * r + b1) * r + 1.0);
    }
    r = q < 0 ? sqrt(-log_p) : sqrt(-log(1.0 - p));
    if (r <= 5.0) { // for p >= 1.389e−11
        r -= 1.6;
        x = (((((((c7 * r + c6) * r + c5) * r + c4) * r + c3) * r + c2) * r + c1) * r + c0)
            / (((((((d7 * r + d6) * r + d5) * r + d4) * r + d3) * r + d2) * r + d1) * r + 1.0);
    } else if (r <= 27) { // for p >= 2.51e-317
        r -= 5.0;
        x = (((((((e7 * r + e6) * r + e5) * r + e4) * r + e3) * r + e2) * r + e1) * r + e0)
            / (((((((f7 * r + f6) * r + f5) * r + f4) * r + f3) * r + f2) * r + f1) * r + 1.0);
    } else if (r < 6.4e8) {               // improvement from Martin Maechler
        double s2 = -2 * log_p;           // = -2*lp = 2s
        double x2 = s2 - log(M_2PI * s2); // = xs_1
        if (r < 36000) {
            x2 = s2 - log(M_2PI * x2) - 2 / (2 + x2);                                  // == xs_2
            if (r < 840) {                                                             // 27 < r < 840
                x2 = s2 - log(M_2PI * x2) + 2 * log1p(-(1 - 1 / (4 + x2)) / (2 + x2)); // == xs_3
                if (r < 109) {                                                         // 27 < r < 109
                    x2 = s2 - log(M_2PI * x2) + 2 * log1p(-(1 - (1 - 5 / (6 + x2)) / (4 + x2)) / (2 + x2)); // == xs_4
                    if (r < 55) { // 27 < r < 55
                        x2 = s2 - log(M_2PI * x2)
                             + 2 * log1p(-(1 - (1 - (5 - 9 / (8 + x2)) / (6 + x2)) / (4 + x2)) / (2 + x2)); // == xs_5
                    }
                }
            }
        }
        x = sqrt(x2);
    } else {
        return r * M_SQRT2;
    }
    return q < 0 ? -x : x;
}

// see pchisq() in https://github.com/wch/r-source/blob/trunk/src/nmath/pchisq.c
// see pgamma() in https://github.com/wch/r-source/blob/trunk/src/nmath/pgamma.c
// see chidist() in https://github.com/statgen/METAL/blob/master/libsrc/MathStats.cpp
// see kf_gammaq() in https://github.com/samtools/htslib/blob/develop/kfunc.c
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290
// regularized lower incomplete gamma function, by series expansion
static double _kf_gammap(double s, double z) {
    double sum, x;
    int k;
    for (k = 1, sum = x = 1.; k < 100; ++k) {
        sum += (x *= z / (s + k));
        if (x / sum < KF_GAMMA_EPS) break;
    }
    return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}
// regularized upper incomplete gamma function, by continued fraction
static double _kf_log_gammaq(double s, double z) {
    int j;
    double C, D, f;
    f = 1. + z - s;
    C = f;
    D = 0.;
    // Modified Lentz's algorithm for computing continued fraction
    // See Numerical Recipes in C, 2nd edition, section 5.2
    for (j = 1; j < 100; ++j) {
        double a = j * (s - j), b = (j << 1) + 1 + z - s, d;
        D = b + a * D;
        if (D < KF_TINY) D = KF_TINY;
        C = b + a / C;
        if (C < KF_TINY) C = KF_TINY;
        D = 1. / D;
        d = C * D;
        f *= d;
        if (fabs(d - 1.) < KF_GAMMA_EPS) break;
    }
    return s * log(z) - z - kf_lgamma(s) - log(f);
}
static double kf_log_gammaq(double s, double z) {
    return z <= 1. || z < s ? log(1.0 - _kf_gammap(s, z)) : _kf_log_gammaq(s, z);
}
static double log_chidist(double x, double df) { return kf_log_gammaq(ldexp(df, -1), ldexp(x, -1)); }

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Run meta-analysis from GWAS-VCF summary statistics.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Run meta-analysis from GWAS-VCF summary statistics. "
           "(version " METAL_VERSION
           " https://github.com/freeseek/score)\n"
           "[ Willer, C. J., Li, Y. & Abecasis, G. R. METAL: fast and efficient meta-analysis of genomewide\n"
           "association scans. Bioinformatics 26, 2190–2191 (2010) http://doi.org/10.1093/bioinformatics/btq340 ]\n"
           "\n"
           "Usage: bcftools +metal [options] <score1.gwas.vcf.gz> "
           "<score2.gwas.vcf.gz> [<score3.gwas.vcf.gz> ...]\n"
           "Plugin options:\n"
           "       --summaries <file>          list of summary statistics VCFs from file\n"
           "   -e, --exclude EXPR              Exclude sites for which the expression is true (see man page for "
           "details)\n"
           "   -i, --include EXPR              Select sites for which the expression is true (see man page for "
           "details)\n"
           "       --szw                       perform meta-analysis based on sample-size weighted scheme\n"
           "                                   rather than inverse-variance weighted scheme\n"
           "       --het                       perform heterogenity analysis\n"
           "       --esd                       output effect size direction across studies\n"
           "       --no-version                do not append version and command line to the header\n"
           "   -o, --output <file>             write output to a file [no output]\n"
           "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "   -r, --regions <region>          restrict to comma-separated list of regions\n"
           "   -R, --regions-file <file>       restrict to regions listed in a file\n"
           "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps "
           "(2) [1]\n"
           "   -t, --targets [^]<region>       restrict to comma-separated list of regions. Exclude regions with \"^\" "
           "prefix\n"
           "   -T, --targets-file [^]<file>    restrict to regions listed in a file. Exclude regions with \"^\" "
           "prefix\n"
           "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps "
           "(2) [0]\n"
           "       --threads <int>             use multithreading with INT worker threads [0]\n"
           "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
           "\n"
           "Examples:\n"
           "      bcftools +metal -Ob -o ukb_mvp.gwas.bcf -i ukb.gwas.bcf mvp.gwas.bcf\n"
           "      bcftools +metal -Ob -o ukb_mvp.gwas.bcf -i 'NS>1000 & AF>0.01 & AF<0.99' ukb.gwas.bcf mvp.gwas.bcf\n"
           "      bcftools +metal -Ob -o ukb_mvp.gwas.bcf -i 'ID=\"rs1234\" || ID=\"rs123456\" || ID=\"rs123\"' "
           "ukb.gwas.bcf mvp.gwas.bcf\n"
           "\n";
}

static inline int filter_test_with_logic(filter_t *filter, bcf1_t *line, uint8_t **smpl_pass, int filter_logic) {
    if (!filter) return 1;
    int i, pass = filter_test(filter, line, (const uint8_t **)smpl_pass);
    if (filter_logic & FLT_EXCLUDE) {
        if (pass) {
            pass = 0;
            if (!(*smpl_pass)) return pass;
            for (i = 0; i < line->n_sample; i++)
                if ((*smpl_pass)[i])
                    (*smpl_pass)[i] = 0;
                else {
                    (*smpl_pass)[i] = 1;
                    pass = 1;
                }
        } else {
            pass = 1;
            if ((*smpl_pass))
                for (i = 0; i < line->n_sample; i++) (*smpl_pass)[i] = 1;
        }
    }
    return pass;
}

int run(int argc, char **argv) {
    int i, j, k, l, rid, idx;
    int filter_logic = 0;
    int szw = 0;
    int het = 0;
    int esd = 0;
    int record_cmd_line = 1;
    int write_index = 0;
    int output_type = FT_VCF;
    int clevel = -1;
    int regions_is_file = 0;
    int regions_overlap = 1;
    int targets_is_file = 0;
    int targets_overlap = 0;
    int n_threads = 0;
    char *tmp = NULL;
    const char *pathname = NULL;
    const char *output_fname = "-";
    char *index_fname;
    const char *regions_list = NULL;
    const char *targets_list = NULL;
    const char *filter_str = NULL;
    int *passes = NULL;
    uint8_t **smpl_passes = NULL;
    filter_t **filters = NULL;
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);
    htsFile *out_fh = NULL;

    static struct option loptions[] = {{"summaries", required_argument, NULL, 1},
                                       {"exclude", required_argument, NULL, 'e'},
                                       {"include", required_argument, NULL, 'i'},
                                       {"szw", no_argument, NULL, 2},
                                       {"het", no_argument, NULL, 3},
                                       {"esd", no_argument, NULL, 4},
                                       {"no-version", no_argument, NULL, 8},
                                       {"output", required_argument, NULL, 'o'},
                                       {"output-type", required_argument, NULL, 'O'},
                                       {"regions", required_argument, NULL, 'r'},
                                       {"regions-file", required_argument, NULL, 'R'},
                                       {"regions-overlap", required_argument, NULL, 5},
                                       {"targets", required_argument, NULL, 't'},
                                       {"targets-file", required_argument, NULL, 'T'},
                                       {"targets-overlap", required_argument, NULL, 6},
                                       {"threads", required_argument, NULL, 9},
                                       {"write-index", optional_argument, NULL, 'W'},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?e:i:o:O:r:R:t:T:W::", loptions, NULL)) >= 0) {
        switch (c) {
        case 1:
            pathname = optarg;
            break;
        case 'e':
            if (filter_str) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
            filter_str = optarg;
            filter_logic |= FLT_EXCLUDE;
            break;
        case 'i':
            if (filter_str) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
            filter_str = optarg;
            filter_logic |= FLT_INCLUDE;
            break;
        case 2:
            szw = 1;
            break;
        case 3:
            het = 1;
            break;
        case 4:
            esd = 1;
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
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 5:
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
        case 6:
            if (!strcasecmp(optarg, "0"))
                targets_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                targets_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                targets_overlap = 2;
            else
                error("Could not parse: --targets-overlap %s\n", optarg);
            break;
        case 9:
            n_threads = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --threads %s\n", optarg);
            break;
        case 'W':
            if (!(write_index = write_index_parse(optarg)))
                error("Unsupported index format '%s'\n", optarg);
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }

    if ((pathname && optind != argc) || (!pathname && optind + 2 > argc)) error("%s", usage_text());

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

    if (n_threads && bcf_sr_set_threads(sr, n_threads) < 0) error("Failed to create threads\n");

    int n_files;
    char **filenames = NULL;
    if (pathname) {
        filenames = hts_readlines(pathname, &n_files);
        if (!filenames) error("Failed to read from file %s\n", pathname);
    } else {
        n_files = argc - optind;
        filenames = argv + optind;
    }
    if (n_files < 2) error("At least 2 summary statistics files required as input while only %d provided\n", n_files);

    for (j = 0; j < n_files; j++)
        if (!bcf_sr_add_reader(sr, filenames[j]))
            error("Error opening %s: %s\n", filenames[j], bcf_sr_strerror(sr->errnum));

    bcf_hdr_t *hdr = NULL;

    if (filter_str) {
        passes = (int *)malloc(sizeof(int) * n_files);
        smpl_passes = (uint8_t **)malloc(sizeof(uint8_t *) * n_files);
        filters = (filter_t **)malloc(sizeof(filter_t *) * n_files);
        for (j = 0; j < n_files; j++) {
            hdr = bcf_sr_get_header(sr, j);
            filters[j] = filter_init(hdr, filter_str);
        }
    }

    hdr = bcf_sr_get_header(sr, 0);
    bcf_hdr_t *out_hdr = bcf_hdr_init("w");
    int *id = (int *)malloc(sizeof(int) * SIZE * n_files);
    int output[SIZE] = {0};
    for (j = 0; j < n_files; j++) {
        hdr = bcf_sr_get_header(sr, j);
        // copy filters information in new header
        for (i = 0; i < hdr->nhrec; i++) {
            bcf_hrec_t *hrec = hdr->hrec[i];
            if (hrec->type == BCF_HL_FLT) {
                // copied from htslib/vcf.c
                int k = bcf_hrec_find_key(hrec, "ID");
                assert(k >= 0); // this should always be true for valid VCFs
                bcf_hrec_t *rec = bcf_hdr_get_hrec(out_hdr, hrec->type, "ID", hrec->vals[k], NULL);
                if (!rec) {
                    int res = bcf_hdr_add_hrec(out_hdr, bcf_hrec_dup(hrec));
                    assert(res == 1);
                }
            }
        }
        // copy contigs information in new header
        for (rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
            const char *seq = hdr->id[BCF_DT_CTG][rid].key;
            if (bcf_hdr_name2id(out_hdr, seq) != -1) continue;
            uint64_t len = hdr->id[BCF_DT_CTG][rid].val->info[0];
            bcf_hdr_printf(out_hdr, len ? "##contig=<ID=%s,length=%" PRIu64 ">" : "##contig=<ID=%s>", seq, len);
        }
        // copy samples information in new header
        for (l = 0; l < bcf_hdr_nsamples(hdr); l++)
            if (bcf_hdr_id2int(out_hdr, BCF_DT_SAMPLE, hdr->samples[l]) < 0)
                bcf_hdr_add_sample(out_hdr, hdr->samples[l]);

        for (idx = 0; idx < SIZE; idx++) {
            id[j * SIZE + idx] = bcf_hdr_id2int(hdr, BCF_DT_ID, id_str[idx]);
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, id[j * SIZE + idx])) id[j * SIZE + idx] = -1;
        }

        if (szw) { // sample-size weighted scheme
            if (id[j * SIZE + NS] < 0 && id[j * SIZE + NE] < 0)
                error("NS or NE FORMAT fields required to compute study weight missing from file %s\n", filenames[j]);
            if (id[j * SIZE + EZ] < 0 && (id[j * SIZE + ES] < 0 || id[j * SIZE + LP] < 0))
                error("EZ or ES and LP FORMAT fields required to compute study Z-score missing from file %s\n",
                      filenames[j]);
        } else { // inverse-variance weighted scheme
            if (id[j * SIZE + SE] < 0)
                error("SE FORMAT field required to compute study weight missing from file %s\n", filenames[j]);
            if (id[j * SIZE + ES] < 0)
                error("ES FORMAT field required to compute study effect size missing from file %s\n", filenames[j]);
        }
    }

    if (bcf_hdr_sync(out_hdr) < 0) error_errno("Failed to update header");
    if (bcf_hdr_nsamples(out_hdr) == 0) error("No summary statistics in the input files");

    // create sample map
    // i is an index 1..n_smpl for the samples in the output VCF
    // j is an index 1..n_files for the input VCFs
    // k is an index 1..i2n[i] for which input VCFs contribute to a given sample in the output VCF
    // l is an index for the samples in the input VCFs
    int n_smpl = bcf_hdr_nsamples(out_hdr);
    int *i2n = (int *)calloc(sizeof(int), n_smpl);
    int **i_k2j = (int **)malloc(sizeof(int *) * n_smpl);
    int **i_k2l = (int **)malloc(sizeof(int *) * n_smpl);
    for (i = 0; i < n_smpl; i++) {
        for (j = 0; j < n_files; j++) {
            hdr = bcf_sr_get_header(sr, j);
            if (bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, out_hdr->samples[i]) != -1) i2n[i]++;
        }
        i_k2j[i] = (int *)malloc(sizeof(int) * i2n[i]);
        i_k2l[i] = (int *)malloc(sizeof(int) * i2n[i]);
        int output_ns = 1;
        int output_nc = 1;
        int output_af = 1;
        int output_ac = 1;
        for (j = 0, k = 0; j < n_files; j++) {
            hdr = bcf_sr_get_header(sr, j);
            int l = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, out_hdr->samples[i]);
            if (l < 0) continue;
            i_k2j[i][k] = j;
            i_k2l[i][k] = l;
            // if one input study misses the FORMAT field there will be no output for that meta-analysis
            if (id[j * SIZE + NS] < 0) output_ns = 0;
            if (id[j * SIZE + NC] < 0) output_nc = 0;
            if (id[j * SIZE + AF] < 0) output_af = 0;
            if (id[j * SIZE + AC] < 0) output_ac = 0;
            k++;
        }
        // if one output meta-analysis can compute a FORMAT field there will be output
        if (output_ns) output[NS] = 1;
        if (output_nc) output[NC] = 1;
        if (output_af) output[AF] = 1;
        if (output_ac) output[AC] = 1;
    }

    // add FORMAT header fields
    if (szw) { // sample-size weighted scheme
        output[EZ] = 1;
        output[NE] = 1;
    } else { // inverse-variance weighted scheme
        output[ES] = 1;
        output[SE] = 1;
        if (output[NC]) output[NE] = 1;
    }
    output[LP] = 1;
    if (het) { // Cochran's Q test
        output[CQ] = 1;
        output[I2] = 1;
    }
    for (idx = 0; idx < SIZE; idx++)
        if (output[idx]
            && bcf_hdr_printf(out_hdr, "##FORMAT=<ID=%s,Number=A,Type=Float,Description=\"%s\">", id_str[idx],
                              desc_str[idx])
                   < 0)
            error_errno("Failed to add \"%s\" FORMAT header", id_str[idx]);
    if (esd
        && bcf_hdr_printf(out_hdr, "##FORMAT=<ID=%s,Number=A,Type=String,Description=\"%s\">", id_str[SIZE],
                          desc_str[SIZE])
               < 0)
        error_errno("Failed to add \"%s\" FORMAT header", id_str[SIZE]);
    if (bcf_hdr_sync(out_hdr) < 0) error_errno("Failed to update header");

    char wmode[8];
    set_wmode(wmode, output_type, output_fname, clevel);
    out_fh = hts_open(output_fname, wmode);
    if (out_fh == NULL) error("Error: cannot write to \"%s\": %s\n", output_fname, strerror(errno));
    if (n_threads) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
    if (record_cmd_line) bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_metal");
    if (bcf_hdr_write(out_fh, out_hdr) < 0) error("Unable to write to output VCF file\n");
    if (init_index2(out_fh, hdr, output_fname, &index_fname, write_index) < 0)
        error("Error: failed to initialise index for %s\n", output_fname);

    // process GWAS-VCF rows
    float *val_arr = (float *)malloc(sizeof(float) * SIZE * n_smpl);
    char *esd_arr = (char *)malloc(sizeof(float) * n_files * n_smpl);
    for (i = 0; i < n_files * n_smpl; i++) esd_arr[i] = bcf_str_vector_end;
    bcf1_t *out_line = bcf_init();
    bcf_float_set_missing(out_line->qual);
    while (bcf_sr_next_line(sr)) {
        if (filters) {
            int pass = 0;
            for (j = 0; j < n_files; j++) {
                if (!bcf_sr_has_line(sr, j)) continue;
                bcf1_t *line = bcf_sr_get_line(sr, j);
                hdr = bcf_sr_get_header(sr, j);
                passes[j] = filter_test_with_logic(filters[j], line, &smpl_passes[j], filter_logic);
                if (passes[j]) pass = 1;
            }
            if (!pass) continue; // skip the line for all input VCFs
        }

        for (i = 0; i < SIZE * n_smpl; i++) bcf_float_set_missing(val_arr[i]);

        for (i = 0; i < n_smpl; i++) {
            double xnum = 0.0;
            double xden = 0.0;
            double ns_sum = 0.0;
            double nc_sum = 0.0;
            double af_sum = 0.0;
            double ac_sum = 0.0;
            double ne_sum = 0.0;
            double cq_sum = 0.0;
            int fill_line = 0;
            int df = -1;
            bcf_update_id(NULL, out_line, NULL);
            bcf_update_filter(out_hdr, out_line, NULL, 0);
            for (k = 0; k < i2n[i]; k++) {
                esd_arr[n_files * i + k] = '?';
                int j = i_k2j[i][k];
                if (!bcf_sr_has_line(sr, j)) continue;
                bcf1_t *line = bcf_sr_get_line(sr, j);
                hdr = bcf_sr_get_header(sr, j);
                if (!fill_line) {
                    const char *name = bcf_hdr_id2name(hdr, line->rid);
                    out_line->rid = bcf_hdr_name2id(out_hdr, name);
                    out_line->pos = line->pos;
                    bcf_update_alleles(out_hdr, out_line, (const char **)line->d.allele, line->n_allele);
                    fill_line = 1;
                }
                int l = i_k2l[i][k];
                if (filters && (!passes[j] || (smpl_passes[j] && !smpl_passes[j][l])))
                    continue; // skip the line for one input VCF

                double val[SIZE]; // set all values to NAN
                for (idx = 0; idx < SIZE; idx++) {
                    bcf_fmt_t *fmt = bcf_get_fmt_id(line, id[j * SIZE + idx]);
                    if (fmt && !bcf_float_is_missing(((float *)fmt->p)[l])
                        && !bcf_float_is_vector_end(((float *)fmt->p)[l]))
                        val[idx] = (double)((float *)fmt->p)[l];
                    else
                        val[idx] = NAN;
                }

                if (esd) {
                    int effect = szw ? EZ : ES;
                    if (!isnan(val[effect]))
                        esd_arr[n_files * i + k] = val[effect] == 0.0 ? '0' : (val[effect] > 0.0 ? '+' : '-');
                }

                if (isnan(val[NE]) && !isnan(val[NS])) {
                    // compute effective sample size for binary traits
                    val[NE] = isnan(val[NC]) ? val[NS] : 4.0 * (val[NS] - val[NC]) * val[NC] / val[NS];
                }

                double wt;
                if (szw) { // sample-size weighted scheme
                    if (isnan(val[NE])) continue;
                    wt = val[NE];
                    if (isnan(val[EZ])) {
                        if (isnan(val[ES]) || isnan(val[LP])) continue;
                        val[EZ] = -inv_log_ndist(-val[LP] * M_LN10 - M_LN2);
                        if (val[ES] < 0) val[EZ] = -val[EZ];
                    }
                    xnum += val[EZ] * sqrt(wt);
                    cq_sum += val[EZ] * val[EZ];
                } else { // inverse-variance weighted scheme
                    if (isnan(val[ES]) || isnan(val[SE])) continue;
                    wt = 1.0 / (val[SE] * val[SE]);
                    xnum += val[ES] * wt;
                    cq_sum += val[ES] * val[ES] * wt;
                }

                xden += wt;
                ns_sum += val[NS];
                nc_sum += val[NC];
                af_sum += val[AF] * wt;
                ac_sum += val[AC];
                ne_sum += val[NE];

                if (line->d.id[0] != '.' || line->d.id[1]) bcf_add_id(NULL, out_line, line->d.id);
                bcf_unpack(line, BCF_UN_FLT);
                for (l = 0; l < line->d.n_flt; l++) {
                    const char *flt = hdr->id[BCF_DT_ID][line->d.flt[l]].key;
                    int flt_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, flt);
                    bcf_add_filter(out_hdr, out_line, flt_id);
                }
                df++;
            }

            if (df < 0) continue;
            if (szw) { // sample-size weighted scheme
                val_arr[n_smpl * EZ + i] = (float)xnum / sqrt(xden);
                val_arr[n_smpl * NE + i] = (float)xden;
            } else { // inverse-variance weighted scheme
                val_arr[n_smpl * ES + i] = (float)(xnum / xden);
                val_arr[n_smpl * SE + i] = (float)sqrt(1.0 / xden);
            }
            val_arr[n_smpl * LP + i] = (float)((-M_LN2 - log_ndist(fabs(xnum / sqrt(xden)))) / M_LN10);
            val_arr[n_smpl * NS + i] = (float)ns_sum;
            val_arr[n_smpl * NC + i] = (float)nc_sum;
            val_arr[n_smpl * AF + i] = (float)(af_sum / xden);
            val_arr[n_smpl * AC + i] = (float)ac_sum;
            val_arr[n_smpl * NE + i] = (float)ne_sum;

            if (het && df) { // Cochran's Q test
                cq_sum -= xnum * xnum / xden;
                val_arr[n_smpl * I2 + i] = (float)(cq_sum < (double)df ? 0.0 : (cq_sum - (double)df) / cq_sum * 100.0);
                val_arr[n_smpl * CQ + i] = (float)(-log_chidist(cq_sum, (double)df) / M_LN10);
            }

            for (idx = 0; idx < SIZE; idx++)
                if (isnan(val_arr[n_smpl * idx + i])) bcf_float_set_missing(val_arr[n_smpl * idx + i]);
        }

        for (idx = 0; idx < SIZE; idx++)
            if (output[idx]) bcf_update_format_float(out_hdr, out_line, id_str[idx], &val_arr[n_smpl * idx], n_smpl);
        if (esd) bcf_update_format_char(out_hdr, out_line, id_str[SIZE], esd_arr, n_files * n_smpl);
        if (bcf_write(out_fh, out_hdr, out_line) < 0) error("Unable to write to output VCF file\n");
    }

    free(val_arr);
    free(esd_arr);
    free(id);
    for (i = 0; i < n_smpl; i++) {
        free(i_k2j[i]);
        free(i_k2l[i]);
    }
    free(i_k2j);
    free(i_k2l);
    free(i2n);
    free(passes);
    free(smpl_passes);
    if (filters) {
        for (j = 0; j < n_files; j++) filter_destroy(filters[j]);
        free(filters);
    }
    if (pathname) {
        for (j = 0; j < n_files; j++) free(filenames[j]);
        free(filenames);
    }
    if (write_index) {
        if (bcf_idx_save(out_fh) < 0) {
            if (hts_close(out_fh) != 0) error("Close failed %s\n", strcmp(output_fname, "-") ? output_fname : "stdout");
            error("Error: cannot write to index %s\n", index_fname);
        }
        free(index_fname);
    }
    if (hts_close(out_fh) < 0) error("Close failed: %s\n", out_fh->fn);
    bcf_sr_destroy(sr);
    bcf_destroy(out_line);
    bcf_hdr_destroy(out_hdr);

    return 0;
}
