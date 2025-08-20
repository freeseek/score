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
#include <getopt.h>
#include <htslib/ksort.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include "bcftools.h"
#include "filter.h"

#define BLUP_VERSION "2025-08-19"

#define AVERAGE_LD_SCORE_DFLT 72.6
#define MEDIAN_CHISQ 0.45493642311957283031 // qchisq(.5,1)

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

// http://github.com/MRCIEU/gwas-vcf-specification
#define NS 0
#define EZ 1
#define NC 2
#define ES 3
#define SE 4
#define LP 5
#define NE 6
#define SIZE 7
static const char *id_str[SIZE] = {"NS", "EZ", "NC", "ES", "SE", "LP", "NE"};

static const char *desc_str[SIZE] = {
    "Variant-specific number of samples/individuals with called genotypes used to test association with specified "
    "trait",                                                                                 // NS
    "Z-score provided if it was used to derive the ES and SE fields",                        // EZ
    "Variant-specific number of cases used to estimate genetic effect (binary traits only)", // NC
    "Effect size estimate relative to the alternative allele",                               // ES
    "Standard error of effect size estimate",                                                // SE
    "-log10 p-value for effect estimate",                                                    // LP
    "Variant-specific effective sample size"};                                               // NE

/****************************************
 * FUNCTION TO COMPUTE Z FROM LOG P     *
 ****************************************/

#define M_2PI 6.283185307179586476925286766559 /* 2*pi */

// Wichura, M. J. Algorithm AS 241: The Percentage Points of the Normal Distribution. Applied Statistics 37, 477 (1988).
// http://doi.org/10.2307/2347330 PPND16 function (algorithm AS241) http://lib.stat.cmu.edu/apstat/241
// see qnorm5() in http://github.com/wch/r-source/blob/trunk/src/nmath/qnorm.c
// see ninv() in http://github.com/statgen/METAL/blob/master/libsrc/MathStats.cpp
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

// this macro from ksort.h defines the function
// double ks_ksmall_double(size_t n, double arr[], size_t kk);
KSORT_INIT_GENERIC(double)

// compute the median of a vector using the ksort library (with iterator)
double get_median(const double *v, int n, int shift) {
    if (n == 0) return NAN;
    double *w = (double *)malloc(n * sizeof(double));
    int i;
    for (i = 0; i < n; i++) w[i] = v[i * shift];
    double ret = ks_ksmall_double((size_t)n, w, (size_t)n / 2);
    if (n % 2 == 0) ret = (ret + w[n / 2 - 1]) * 0.5f;
    free(w);
    return ret;
}

// compute the median of a vector using the ksort library (with iterator)
double get_median2(const double *v, int n, int shift) {
    if (n == 0) return NAN;
    double *w = (double *)malloc(n * sizeof(double));
    int i;
    for (i = 0; i < n; i++) w[i] = v[i * shift] * v[i * shift];
    double ret = ks_ksmall_double((size_t)n, w, (size_t)n / 2);
    if (n % 2 == 0) ret = (ret + w[n / 2 - 1]) * 0.5f;
    free(w);
    return ret;
}

/****************************************
 * BASIC SPARSE MATRIX MANIPULATION     *
 ****************************************/

// A sparse matrix in COOrdinate format
// see http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
typedef struct {
    int i;
    int j;
    double x;
} coo_cell_t;

typedef struct {
    int nrow;
    int nnz;
    int m_d;
    double *d; // elements on the diagonal
    int m;
    coo_cell_t *cell;
} coo_matrix_t;

// Compressed Sparse Row matrix structure
// see http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
typedef struct {
    int nrow;  // number of rows
    double *d; // elements on the diagonal
    int *p;    // ptr to the row starts, length n+1
    int *j;    // column index, length j[nrow]
    double *x; // cell values, length j[nrow]
} csr_matrix_t;

static void coo_clear(coo_matrix_t *coo) {
    coo->nrow = 0;
    coo->nnz = 0;
    memset((void *)coo->d, 0, sizeof(double) * coo->m_d);
}

static void coo_destroy(coo_matrix_t *coo) {
    free(coo->d);
    free(coo->cell);
}

static void csr_destroy(csr_matrix_t *csr) {
    free(csr->d);
    free(csr->p);
    free(csr->j);
    free(csr->x);
}

static inline void append_diag(int i, double x, coo_matrix_t *coo) {
    hts_expand0(double, i + 1, coo->m_d, coo->d);
    coo->d[i] = x;
}

static inline void append_nnz(int i, int j, double x, coo_matrix_t *coo) {
    coo->nnz++;
    hts_expand(coo_cell_t, coo->nnz, coo->m, coo->cell);
    coo_cell_t *cell = &coo->cell[coo->nnz - 1];
    cell->i = i;
    cell->j = j;
    cell->x = x;
}

// see http://github.com/rgl-epfl/cholespy/blob/main/src/cholesky_solver.cpp
static inline void coo_to_csr(const coo_matrix_t *coo, csr_matrix_t *csr) {
    int k;
    csr->nrow = coo->nrow;
    csr->d = csr->nrow ? (double *)calloc(sizeof(double), csr->nrow) : NULL;
    memcpy(csr->d, coo->d, sizeof(double) * (coo->nrow > coo->m_d ? coo->m_d : coo->nrow));
    csr->p = (int *)calloc(sizeof(int), csr->nrow + 1);
    csr->j = (int *)malloc(sizeof(int) * coo->nnz);
    csr->x = (double *)malloc(sizeof(double) * coo->nnz);

    for (k = 0; k < coo->nnz; k++) csr->p[coo->cell[k].i + 1]++;

    csr->p[0] = 0;
    for (k = 0; k < csr->nrow; k++) csr->p[k + 1] += csr->p[k];

    for (k = 0; k < coo->nnz; k++) {
        int row = coo->cell[k].i;
        int dst = csr->p[row];
        csr->j[dst] = coo->cell[k].j;
        csr->x[dst] = coo->cell[k].x;
        csr->p[row]++;
    }

    for (k = csr->nrow; k > 0; k--) csr->p[k] = csr->p[k - 1];
    csr->p[0] = 0;
}

/****************************************
 * MATRIX MULTIPLICATION AND DIVISION   *
 ****************************************/

// return P += S
static inline void add_matrix(const coo_matrix_t *S, coo_matrix_t *P) {
    if (S->nrow != P->nrow) error("Error: Sigma and Precision matrix have different dimensions\n");

    // add diagonal elements
    int k, n = S->nrow > S->m_d ? S->m_d : S->nrow;
    hts_expand0(double, n, P->m_d, P->d);
    for (k = 0; k < n; k++) P->d[k] += S->d[k];

    // add non-diagonal elements
    for (k = 0; k < S->nnz; k++) append_nnz(S->cell[k].i, S->cell[k].j, S->cell[k].x, P);
}

// return A ./ x
static inline void matdiv(coo_matrix_t *A, double x) {
    int k, n = A->nrow > A->m_d ? A->m_d : A->nrow;
    for (k = 0; k < n; k++) A->d[k] /= x;
    for (k = 0; k < A->nnz; k++) A->cell[k].x /= x;
}

// return A * x for square matrices with diagonal elements
static inline void sdmult(const csr_matrix_t *A, const double *x, double *y) {
    int i, j;
    for (i = 0; i < A->nrow; i++) {
        y[i] = A->d[i] * x[i];
        for (j = A->p[i]; j < A->p[i + 1]; j++) y[i] += A->x[j] * x[A->j[j]];
    }
}

// return A * x for rectangular matrices without diagonal elements
static inline void rect_sdmult(const csr_matrix_t *A, const double *x, double *y) {
    int i, j;
    for (i = 0; i < A->nrow; i++) {
        y[i] = 0.0;
        for (j = A->p[i]; j < A->p[i + 1]; j++) y[i] += A->x[j] * x[A->j[j]];
    }
}

// return x' * y
static inline double dot(const double *x, const double *y, int n) {
    double ret = 0.0;
    int i;
    for (i = 0; i < n; i++) ret += x[i] * y[i];
    return ret;
}

// conjugate gradient method with Jacobi preconditioner
// http://en.wikipedia.org/wiki/Conjugate_gradient_method#Example_code_in_MATLAB_/_GNU_Octave
// http://en.wikipedia.org/wiki/Preconditioner#Jacobi_(or_diagonal)_preconditioner
// http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
// alternative approach to http://github.com/awohns/ldgm/blob/main/MATLAB/precisionDivide.m
static int pcg(const csr_matrix_t *A, double *x, double tol, int jacobi) {
    int i, iter, n = A->nrow;
    double rsold, rsnew, tol2 = tol * tol;
    double *p = (double *)malloc(sizeof(double) * n);
    double *r = (double *)malloc(sizeof(double) * n);
    double *Ap = (double *)malloc(sizeof(double) * n);
    double *z = jacobi ? (double *)malloc(sizeof(double) * n) : r;

    sdmult(A, x, Ap);
    for (i = 0; i < n; i++) r[i] = x[i] - Ap[i];
    if (jacobi)
        for (i = 0; i < n; i++) z[i] = r[i] / A->d[i]; // Jacobi preconditioning
    for (i = 0; i < n; i++) p[i] = z[i];
    rsold = dot(r, z, n);
    for (iter = 0; iter < n; iter++) {
        sdmult(A, p, Ap);
        double alpha = rsold / dot(p, Ap, n);
        for (i = 0; i < n; i++) x[i] += alpha * p[i];
        for (i = 0; i < n; i++) r[i] -= alpha * Ap[i];
        if (jacobi)
            for (i = 0; i < n; i++) z[i] = r[i] / A->d[i]; // Jacobi preconditioning
        rsnew = dot(r, z, n);
        if (rsnew < tol2) break;
        double beta = rsnew / rsold;
        for (i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
        rsold = rsnew;
    }

    free(p);
    free(r);
    free(Ap);
    if (jacobi) free(z);
    return iter;
}

// see http://github.com/awohns/ldgm/blob/main/MATLAB/precisionMultiply.m
// computes x = (P/P00)y where P = [P00, P01; P10, P11] and P/P00 is the Schur complement
// http://en.wikipedia.org/wiki/Schur_complement
static int precision_multiply(csr_matrix_t schur[][2], const double *y1, double tol, int jacobi, double *x1) {
    int n0 = schur[0][0].nrow;
    int n1 = schur[1][1].nrow;
    double *tmp = (double *)malloc(sizeof(double) * (n0 > n1 ? n0 : n1));
    rect_sdmult(&schur[0][1], y1, tmp);
    int i, n_iter = pcg(&schur[0][0], tmp, tol, jacobi);
    rect_sdmult(&schur[1][0], tmp, x1);
    sdmult(&schur[1][1], y1, tmp);
    for (i = 0; i < n1; i++) x1[i] = tmp[i] - x1[i];
    free(tmp);
    return n_iter;
}

/****************************************
 * LDGM-VCF ROUTINES                    *
 ****************************************/

// LDGM-VCF fields
typedef struct {
    int ld_node;
    int n_line; // map to the location of the corresponding GWAS-VCF line, if found
    double ez_deriv;
    double sqrt_het; // sqrt(2pq)
    double sd;       // SE(beta) = 1 / (2pq) / N_eff
    double ne;       // effective sample size
} row_t;

// GWAS-VCF fields
typedef struct {
    int ld_node;
    int aa;
    bcf1_t *line;
} line_t;

typedef struct {
    int imap; // map to the summary statistic sample number in the GWAS-VCF file
    const char *seqname;
    int ld_block;
    double neff;
    double mean_neff;
    int row_ptr;
    int n_missing;
    int n_iter; // number of iterations required to perform precisionMultiply()

    row_t *rows;
    int m_rows;
    int *node2row; // map to the location of the LDGM-VCF fields, if found
    int n_node2row;
    int m_node2row;

    coo_matrix_t coo;
    coo_matrix_t coo_schur[2][2];
    int m_schur_imap;
    int *schur_imap;

    // statistics relevant across multiple LD blocks
    double all_trace_inf;
    int all_n_missing;
    int all_n_non_missing;
} ld_block_t;

static inline void ld_block_clear(ld_block_t *block) {
    block->n_missing = 0;
    block->n_node2row = 0;
    memset((void *)block->node2row, 0, sizeof(int) * block->m_node2row);
    coo_clear(&block->coo);
    int k;
    for (k = 0; k < 4; k++) coo_clear(&block->coo_schur[k / 2][k % 2]);
}

static inline void ld_block_destroy(ld_block_t *block) {
    free(block->rows);
    free(block->node2row);
    coo_destroy(&block->coo);
    int k;
    for (k = 0; k < 4; k++) coo_destroy(&block->coo_schur[k / 2][k % 2]);
    free(block->schur_imap);
}

// adapted from Petr Danecek's implementation of remove_format() in bcftools/vcfannotate.c
static void bcf_remove_format(bcf1_t *line) {
    // remove all FORMAT fields
    if (!(line->unpacked & BCF_UN_FMT)) bcf_unpack(line, BCF_UN_FMT);

    int i;
    for (i = 0; i < line->n_fmt; i++) {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        if (fmt->p_free) {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        line->d.indiv_dirty = 1;
        fmt->p = NULL;
    }
}

// verify the GWAS-VCF file has sufficient information to compute Z-score using one of the following strategies:
// use EZ if available
// use ES/SE if available
// use -qnorm(10^-LP/2) * sign(ES) if available
// verify the GWAS-VCF file has sufficient information to compute the effective population size:
// use --sample-sizes if available
// use NE if available
// use NS/NC if available
// use NS if available
static inline void check_gwas(bcf_hdr_t *hdr, int n_sample_sizes) {
    int idx, id[SIZE];
    for (idx = 0; idx < SIZE; idx++) {
        id[idx] = bcf_hdr_id2int(hdr, BCF_DT_ID, id_str[idx]);
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, id[idx])) id[idx] = -1;
    }
    if (id[EZ] < 0 && (id[ES] < 0 || (id[SE] < 0 && id[LP] < 0)))
        error(
            "Error: Either the FORMAT field EZ, or ES and SE, or ES and LP must be defined in the header of the "
            "GWAS-VCF summary statistics file\n");
    if (!n_sample_sizes && id[NE] < 0 && id[NS] < 0)
        error(
            "Error: Either the FORMAT field NE or NS must be defined in the header of the GWAS-VCF summary statistics "
            "file\nor else effective sample sizes must be input with the --sample-sizes option\n");
}

// verify a LDGM-VCF file header is compliant
static inline void check_ldgm(bcf_hdr_t *hdr) {
    static const char *info[] = {"AA", "AF", "LD_block", "LD_node", "LD_diagonal", "LD_neighbors", "LD_weights"};
    int i;
    for (i = 0; i < sizeof(info) / sizeof(char *); i++)
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, info[i])))
            error("Error: The INFO field %s is not defined in the header of the LDGM-VCF precision matrix file\n",
                  info[i]);
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

static int read_ld_block(bcf_srs_t *sr, ld_block_t *blocks, int n_pops, double alpha_param, line_t **lines,
                         int *n_lines, int *m_lines, filter_t *filter, int filter_logic) {
    int pop, idx, i, k;
    int *int_arr = (int *)calloc(sizeof(int), 1);
    int n_int_arr, m_int_arr = 1;
    float *float_arr = NULL;
    int n_float_arr, m_float_arr = 0;

    bcf1_t *line = NULL;
    bcf_hdr_t *hdr = NULL;
    double *ez = (double *)malloc(sizeof(double) * n_pops);
    double *lp = (double *)malloc(sizeof(double) * n_pops);
    double *ne = (double *)malloc(sizeof(double) * n_pops);

    int aa, ld_block, ld_node;
    float ld_diagonal;
    double af;
    int block_started = 1;
    int block_ended = 0;
    int curr_ld_block = -1;
    int ret = bcf_sr_has_line(sr, 0);
    for (pop = 0; pop < n_pops; pop++) {
        ret += bcf_sr_has_line(sr, 1 + pop);
        ld_block_clear(&blocks[pop]);
    }

    do {
        // populate GWAS-VCF data in temporary structures if data available
        int pass = 0;
        uint8_t *smpl_pass = NULL;
        if (bcf_sr_has_line(sr, 0)) {
            line = bcf_sr_get_line(sr, 0);
            hdr = bcf_sr_get_header(sr, 0);
            pass = filter_test_with_logic(filter, line, &smpl_pass, filter_logic);
        }

        if (pass) {
            // check if the VCF line has enough information to compute the Z-score
            bcf_fmt_t *fmt[SIZE];
            for (idx = 0; idx < SIZE; idx++) fmt[idx] = bcf_get_fmt(hdr, line, id_str[idx]);

            for (pop = 0; pop < n_pops; pop++) {
                int pop_ind = blocks[pop].imap;
                int ind_pass = !smpl_pass || smpl_pass[pop_ind];
                double val[SIZE];
                for (idx = 0; idx < SIZE; idx++)
                    if (ind_pass && fmt[idx] && !bcf_float_is_missing(((float *)fmt[idx]->p)[pop_ind])
                        && !bcf_float_is_vector_end(((float *)fmt[idx]->p)[pop_ind]))
                        val[idx] = (double)((float *)fmt[idx]->p)[pop_ind];
                    else
                        val[idx] = NAN;
                if (isnan(val[EZ]) && !isnan(val[ES])) {
                    if (!isnan(val[SE])) {
                        val[EZ] = val[ES] / val[SE];
                    } else if (!isnan(val[LP])) {
                        val[EZ] = -inv_log_ndist(-val[LP] * M_LN10 - M_LN2);
                        if (val[ES] < 0) val[EZ] = -val[EZ];
                    }
                }
                if (blocks[pop].neff) { // force effective sample size regardless of what found in the GWAS-VCF
                    val[NE] = blocks[pop].neff;
                } else if (isnan(val[NE]) && !isnan(val[NS])) {
                    // compute effective sample size for binary traits
                    val[NE] = isnan(val[NC]) ? val[NS] : 4.0 * (val[NS] - val[NC]) * val[NC] / val[NS];
                }
                ez[pop] = val[EZ];
                lp[pop] = val[LP];
                ne[pop] = val[NE];
            }
        } else {
            for (pop = 0; pop < n_pops; pop++) {
                ez[pop] = NAN;
                lp[pop] = NAN;
                ne[pop] = NAN;
            }
        }

        // load LDGM-VCF data
        int save_line = 0;
        for (pop = 0; pop < n_pops; pop++) {
            if (!bcf_sr_has_line(sr, 1 + pop)) continue; // no line in the LDGM-VCF file
            line = bcf_sr_get_line(sr, 1 + pop);
            hdr = bcf_sr_get_header(sr, 1 + pop);

            // load data for the VCF record while verifying compliancy with the LDGM-VCF specification
            n_int_arr = bcf_get_info_int32(hdr, line, "AA", &int_arr, &m_int_arr);
            if (n_int_arr != 1 || (int_arr[0] != 0 && int_arr[0] != 1))
                error("Error: AA INFO field from file %s is nonconformal at %s:%" PRId64 "\n",
                      (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);
            aa = int_arr[0];
            n_float_arr = bcf_get_info_float(hdr, line, "AF", &float_arr, &m_float_arr);
            if (n_float_arr != 1)
                error("Error: AF INFO field from file %s is nonconformal at %s:%" PRId64 "\n",
                      (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);
            af = (double)float_arr[0];
            n_int_arr = bcf_get_info_int32(hdr, line, "LD_block", &int_arr, &m_int_arr);
            if (n_int_arr != 1)
                error("Error: LD_block INFO field from file %s is nonconformal at %s:%" PRId64 "\n",
                      (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);
            ld_block = int_arr[0];
            n_int_arr = bcf_get_info_int32(hdr, line, "LD_node", &int_arr, &m_int_arr);
            if (n_int_arr != 1 || int_arr[0] < 0)
                error("Error: LD_node INFO field from file %s is nonconformal at %s:%" PRId64 "\n",
                      (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);
            ld_node = int_arr[0];
            n_float_arr = bcf_get_info_float(hdr, line, "LD_diagonal", &float_arr, &m_float_arr);
            if (n_float_arr != 1 || float_arr[0] < 1.0f)
                error("Error: LD_diagonal INFO field from file %s is nonconformal at %s:%" PRId64 "\n",
                      (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);
            ld_diagonal = float_arr[0];
            n_int_arr = bcf_get_info_int32(hdr, line, "LD_neighbors", &int_arr, &m_int_arr);
            for (i = 0; i < n_int_arr; i++)
                if (int_arr[i] <= ld_node)
                    error("Error: LD_neighbors INFO field from file %s is nonconformal at %s:%" PRId64 "\n",
                          (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);
            n_float_arr = bcf_get_info_float(hdr, line, "LD_weights", &float_arr, &m_float_arr);
            //      this currently happens, though really it should not
            //      for (i=0; i<n_float_arr; i++)
            //        if (float_arr[i] == 0.0f)
            //          error("Error: LD_weights INFO field is nonconformal at %s:%"PRId64"\n", bcf_seqname(hdr, line),
            //          (int64_t)line->pos+1);
            if (n_int_arr != n_float_arr)
                error("Error: arrays LD_neighbors and LD_weights from file %s have different lengths at %s:%" PRId64
                      "\n",
                      (bcf_sr_get_reader(sr, 1 + pop))->fname, bcf_seqname(hdr, line), (int64_t)line->pos + 1);

            if (block_started && blocks[pop].ld_block == ld_block)
                continue; // skip first line if the reader is still stuck on the previous LDGM
            // skip line if the reader reached the next LDGM and stop reading more lines
            if (!block_started && curr_ld_block != ld_block) {
                block_ended = 1;
                continue;
            }
            curr_ld_block = ld_block;

            ld_block_t *block = &blocks[pop];
            if (ld_node >= block->n_node2row) block->n_node2row = ld_node + 1;
            hts_expand0(int, block->n_node2row, block->m_node2row, block->node2row);

            // if not already loaded, add LDGM-VCF entry
            row_t *row;
            if (block->node2row[ld_node] == 0) {
                hts_expand(row_t, block->coo.nrow + 1, block->m_rows, block->rows);
                block->node2row[ld_node] = block->coo.nrow + 1;
                row = &block->rows[block->coo.nrow];

                row->ld_node = ld_node;
                // flip the allele frequency if the alternate allele is the ancestral allele
                row->sqrt_het = sqrt(2.0 * (double)af * (1.0 - (double)af));
                row->sd = alpha_param == 0.0 ? row->sqrt_het : pow(row->sqrt_het, alpha_param + 1.0);

                if (bcf_sr_has_line(sr, 0) && !isnan(ez[pop])) {
                    row->n_line = *n_lines + 1;
                    save_line = 1;
                    // flip the Z-score if the alternate allele is the ancestral allele
                    row->ez_deriv = aa ? -ez[pop] : ez[pop];
                    row->ne = ne[pop];
                } else {
                    row->n_line = 0;
                    row->ez_deriv = NAN;
                    row->ne = NAN;
                }

                append_diag(block->coo.nrow, (double)ld_diagonal, &block->coo);
                for (i = 0; i < n_int_arr; i++) {
                    if (float_arr[i] == 0.0f)
                        continue; // there should not be need for this check once they fix the LDGM precision matrices
                    append_nnz(ld_node, int_arr[i], (double)float_arr[i], &block->coo);
                    append_nnz(int_arr[i], ld_node, (double)float_arr[i], &block->coo);
                }
                block->coo.nrow++;
            } else {
                row = &block->rows[block->node2row[ld_node] - 1];
                if (!row->n_line && bcf_sr_has_line(sr, 0) && !isnan(ez[pop])) {
                    row->n_line = *n_lines + 1;
                    save_line = 1;
                    // flip the Z-score if the alternate allele is the ancestral allele
                    row->ez_deriv = aa ? -ez[pop] : ez[pop];
                    row->ne = ne[pop];
                }
            }
        }
        if (ret > bcf_sr_has_line(sr, 0)) {
            block_started = 0;
            for (pop = 0; pop < n_pops; pop++) {
                blocks[pop].seqname = bcf_hdr_id2name(hdr, line->rid);
                blocks[pop].ld_block = curr_ld_block;
            }
        }

        // load GWAS-VCF data if it is references in any of the LDGM-VCF structures
        if (save_line) {
            hts_expand0(line_t, *n_lines + 1, *m_lines, *lines);
            line_t *curr_line = &(*lines)[*n_lines];
            curr_line->ld_node = ld_node;
            curr_line->aa = aa;
            line = bcf_sr_get_line(sr, 0);
            curr_line->line = bcf_dup(line);
            // remove all format fields from the line
            bcf_remove_format(curr_line->line);
            (*n_lines)++;
        }
    } while (!block_ended && (ret = bcf_sr_next_line(sr)));

    for (pop = 0; pop < n_pops; pop++) {
        ld_block_t *block = &blocks[pop];
        block->row_ptr = pop == 0 ? 0 : blocks[pop - 1].row_ptr + blocks[pop - 1].coo.nrow;

        // compute number of missing rows from the summary statistics and average sample size
        block->mean_neff = 0.0;
        int l = 0;
        for (k = 0; k < block->coo.nrow; k++) {
            row_t *row = &block->rows[k];
            if (row->n_line == 0)
                block->n_missing++;
            else if (!isnan(row->ne)) {
                block->mean_neff += row->ne;
                l++;
            }
        }
        block->mean_neff /= l;
        block->all_n_missing += block->n_missing;
        block->all_n_non_missing += block->coo.nrow - block->n_missing;

        // compress all LDGM edges so that they refer to rows rather than LD nodes
        for (k = 0; k < block->coo.nnz; k++) {
            coo_cell_t *cell = &block->coo.cell[k];
            if (cell->j >= block->n_node2row || block->node2row[cell->j] == 0)
                error(
                    "Error: LDGM node %d in LD_block %d from file %s lists neighbor node %d which was not observed in "
                    "the LDGM\n",
                    cell->i, block->ld_block, (bcf_sr_get_reader(sr, 1 + pop))->fname, cell->j);
            cell->i = block->node2row[cell->i] - 1;
            cell->j = block->node2row[cell->j] - 1;
        }
    }

    free(int_arr);
    free(float_arr);
    free(ez);
    free(lp);
    free(ne);
    return ret;
}

static void schur_split(ld_block_t *block, csr_matrix_t schur[][2]) {
    int k, n0 = 0, n1 = 0;
    const coo_matrix_t *coo = &block->coo;
    hts_expand(int, coo->nrow, block->m_schur_imap, block->schur_imap);

    // computes the sizes of P00 and P11 and add diagonal elements
    for (k = 0; k < coo->nrow; k++) {
        block->schur_imap[k] = block->rows[k].n_line ? n1++ : n0++;
        if (k >= coo->m_d) break;
        int kk = block->rows[k].n_line > 0;
        append_diag(block->schur_imap[k], coo->d[k], &block->coo_schur[kk][kk]);
    }

    // add non-diagonal elements to P00, P01, P10, and P11
    for (k = 0; k < coo->nnz; k++) {
        int i = coo->cell[k].i;
        int j = coo->cell[k].j;
        int ii = block->rows[i].n_line > 0;
        int jj = block->rows[j].n_line > 0;
        append_nnz(block->schur_imap[i], block->schur_imap[j], coo->cell[k].x, &block->coo_schur[ii][jj]);
    }

    // convert P00, P01, P10, and P11 from COO to CSR representation
    for (k = 0; k < 4; k++) {
        block->coo_schur[k / 2][k % 2].nrow = k < 2 ? n0 : n1;
        coo_to_csr(&block->coo_schur[k / 2][k % 2], &schur[k / 2][k % 2]);
    }
}

static void make_sigma(ld_block_t *blocks, int n_pops, double beta_cov, double cross_corr, coo_matrix_t *out) {
    coo_clear(out);
    out->nrow = blocks[n_pops - 1].row_ptr + blocks[n_pops - 1].coo.nrow;

    // compute the largest node number
    int pop, k, ld_node, pop1, pop2, n_ld_nodes = 0;
    for (pop = 0; pop < n_pops; pop++) {
        ld_block_t *block = &blocks[pop];
        if (block->n_node2row > n_ld_nodes) n_ld_nodes = block->n_node2row;
    }

    // fill diagonal elements by looping across all rows
    for (pop = 0; pop < n_pops; pop++) {
        ld_block_t *block = &blocks[pop];
        for (k = 0; k < block->coo.nrow; k++) {
            row_t *row = &block->rows[k];
            if (row->n_line == 0) continue;
            append_diag(block->row_ptr + k, beta_cov * row->sd * row->sd, out);
        }
    }

    // fill non-diagonal elements by looping across all nodes
    for (ld_node = 0; ld_node < n_ld_nodes; ld_node++) {
        for (pop1 = 0; pop1 < n_pops; pop1++) {
            ld_block_t *block1 = &blocks[pop1];
            // SNP missing from the summary statistics for population pop1
            if (ld_node >= block1->n_node2row || block1->node2row[ld_node] == 0) continue;
            row_t *row1 = &block1->rows[block1->node2row[ld_node] - 1];
            if (row1->n_line == 0) continue;
            for (pop2 = pop1 + 1; pop2 < n_pops; pop2++) {
                ld_block_t *block2 = &blocks[pop2];
                // SNP missing from the summary statistics for population pop2
                if (ld_node >= block2->n_node2row || block2->node2row[ld_node] == 0) continue;
                row_t *row2 = &block2->rows[block2->node2row[ld_node] - 1];
                if (row2->n_line == 0) continue;
                int row = block1->row_ptr + block1->node2row[ld_node] - 1;
                int col = block2->row_ptr + block2->node2row[ld_node] - 1;
                append_nnz(row, col, cross_corr * beta_cov * row1->sd * row2->sd, out);
                append_nnz(col, row, cross_corr * beta_cov * row1->sd * row2->sd, out);
            }
        }
    }
}

static void concatenate(const ld_block_t *blocks, int n_pops, coo_matrix_t *out) {
    coo_clear(out);
    out->nrow = blocks[n_pops - 1].row_ptr + blocks[n_pops - 1].coo.nrow;

    hts_expand0(double, out->nrow, out->m_d, out->d);
    out->nnz = 0;
    int pop, k;
    for (pop = 0; pop < n_pops; pop++) {
        const ld_block_t *block = &blocks[pop];
        const coo_matrix_t *coo = &block->coo;
        memcpy(&out->d[block->row_ptr], coo->d, sizeof(double) * (coo->nrow > coo->m_d ? coo->m_d : coo->nrow));
        for (k = 0; k < coo->nnz; k++)
            append_nnz(block->row_ptr + coo->cell[k].i, block->row_ptr + coo->cell[k].j, coo->cell[k].x, out);
    }
}

static void write_ld_block(htsFile *fh, bcf_hdr_t *hdr, line_t *lines, int n_lines, ld_block_t *blocks, int n_pops) {
    float *es_arr = (float *)malloc(sizeof(float) * n_pops);

    int k, pop;
    for (k = 0; k < n_lines; k++) {
        for (pop = 0; pop < n_pops; pop++) {
            ld_block_t *block = &blocks[pop];
            // check there is a matching row in the LDGM matrix and that the row maps to this line to avoid duplicating
            // loadings
            if (lines[k].ld_node >= block->n_node2row || block->node2row[lines[k].ld_node] == 0
                || block->rows[block->node2row[lines[k].ld_node] - 1].n_line - 1 != k
                || isnan(block->rows[block->node2row[lines[k].ld_node] - 1].ez_deriv)) {
                bcf_float_set_missing(es_arr[pop]);
            } else {
                row_t *row = &block->rows[block->node2row[lines[k].ld_node] - 1];
                double ez = lines[k].aa ? -row->ez_deriv : row->ez_deriv;
                es_arr[pop] = (float)ez;
            }
        }
        bcf_update_format_float(hdr, lines[k].line, id_str[ES], es_arr, n_pops);
        if (bcf_write(fh, hdr, lines[k].line) < 0) error("Error: Unable to write to output VCF file\n");
        bcf_destroy(lines[k].line);
    }

    free(es_arr);
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Compute best linear unbiased predictor from GWAS-VCF summary statistics.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Compute best linear unbiased predictor from GWAS-VCF summary statistics.\n"
           "(version " BLUP_VERSION
           " http://github.com/freeseek/score)\n"
           "[ Nowbandegani, P. S., Wohns, A. W., et al. Extremely sparse models of linkage disequilibrium in "
           "ancestrally\n"
           "diverse association studies. Nat Genet, 55, 1494–1502, (2023) http://doi.org/10.1038/s41588-023-01487-8 ]\n"
           "\n"
           "Usage: bcftools +blup [options] <score.gwas.vcf.gz> [<ldgm.vcf.gz> <ldgm2.vcf.gz> ...]\n"
           "Plugin options:\n"
           "   -v, --verbose                   verbose output\n"
           "       --ldgm-vcfs <list>          List of LDGM-VCF files to use\n"
           "       --ldgm-vcfs-file <file>     File of list of LDGM-VCF files to use\n"
           "   -e, --exclude EXPR              Exclude sites for which the expression is true (see man page for "
           "details)\n"
           "   -i, --include EXPR              Select sites for which the expression is true (see man page for "
           "details)\n"
           "       --no-version                do not append version and command line to the header\n"
           "   -o, --output <file>             write output to a file [no output]\n"
           "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "   -l, --log <file>                write log to file [standard error]\n"
           "   -r, --regions <region>          restrict to comma-separated list of regions\n"
           "   -R, --regions-file <file>       restrict to regions listed in a file\n"
           "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps "
           "(2) [1]\n"
           "   -s, --samples <list>            List of summary statitics to include\n"
           "   -S, --samples-file <file>       File of list of summary statistics to include\n"
           "   -t, --targets [^]<region>       restrict to comma-separated list of regions. Exclude regions with \"^\" "
           "prefix\n"
           "   -T, --targets-file [^]<file>    restrict to regions listed in a file. Exclude regions with \"^\" "
           "prefix\n"
           "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps "
           "(2) [0]\n"
           "       --threads <int>             use multithreading with INT worker threads [0]\n"
           "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
           "\n"
           "Model options:\n"
           "       --stats-only                only compute suggested summary options for a given alpha parameter\n"
           "       --average-ld-score <float>  average LD score per marker [72.6]\n"
           "   -a, --alpha-param <float>       alpha parameter [-0.5]\n"
           "   -b, --beta-cov <float>          frequency-dependent architecture parameter [1e-7]\n"
           "   -x, --cross-corr <float>        cross ancestry correlation parameter [0.9]\n"
           "       --sample-sizes <list>       List of sample sizes for each input summary statistic [estimated from "
           "NS/NC/NE fields]\n"
           "\n"
           "Linear algebra options:\n"
           "       --tolerance <float>         Tolerance threshold for the conjugate gradient [1e-6]\n"
           "       --no-jacobi                 Do not use Jacobi preconditioning when solving linear systems with "
           "conjugate gradient\n"
           "\n"
           "Examples:\n"
           "      bcftools +blup -Ob -o ukb.blup.bcf -b 2e-7 ukb.gwas.bcf 1kg_ldgm.EUR.bcf\n"
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

int run(int argc, char **argv) {
    int i, pop, k, l;
    double average_ld_score = AVERAGE_LD_SCORE_DFLT;
    int verbose = 0;
    int n_sample_sizes = 0;
    char **sample_sizes = NULL;
    int n_files = 0;
    char **filenames = NULL;
    int filter_logic = 0;
    int record_cmd_line = 1;
    int write_index = 0;
    int output_type = FT_VCF;
    int clevel = -1;
    int regions_is_file = 0;
    int regions_overlap = 1;
    int sample_is_file = 0;
    int targets_is_file = 0;
    int targets_overlap = 0;
    int n_threads = 0;
    int stats_only = 0;
    double alpha_param = -0.5;
    double beta_cov = NAN;
    double cross_corr = 0.9;
    double tol = 1e-6;
    int jacobi = 1;
    char *tmp = NULL;
    const char *output_fname = "-";
    char *index_fname;
    const char *regions_list = NULL;
    const char *sample_list = NULL;
    const char *targets_list = NULL;
    const char *filter_str = NULL;
    filter_t *filter = NULL;
    FILE *log_file = stderr;
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);

    static struct option loptions[] = {{"verbose", no_argument, NULL, 'v'},
                                       {"ldgm-vcfs", required_argument, NULL, 1},
                                       {"ldgm-vcfs-file", required_argument, NULL, 2},
                                       {"exclude", required_argument, NULL, 'e'},
                                       {"include", required_argument, NULL, 'i'},
                                       {"no-version", no_argument, NULL, 8},
                                       {"output", required_argument, NULL, 'o'},
                                       {"output-type", required_argument, NULL, 'O'},
                                       {"log", required_argument, NULL, 'l'},
                                       {"regions", required_argument, NULL, 'r'},
                                       {"regions-file", required_argument, NULL, 'R'},
                                       {"regions-overlap", required_argument, NULL, 3},
                                       {"samples", required_argument, NULL, 's'},
                                       {"samples-file", required_argument, NULL, 'S'},
                                       {"targets", required_argument, NULL, 't'},
                                       {"targets-file", required_argument, NULL, 'T'},
                                       {"targets-overlap", required_argument, NULL, 4},
                                       {"threads", required_argument, NULL, 9},
                                       {"write-index", optional_argument, NULL, 'W'},
                                       {"stats-only", no_argument, NULL, 5},
                                       {"average-ld-score", required_argument, NULL, 6},
                                       {"alpha-param", required_argument, NULL, 'a'},
                                       {"beta-cov", required_argument, NULL, 'b'},
                                       {"cross-corr", required_argument, NULL, 'x'},
                                       {"sample-sizes", required_argument, NULL, 7},
                                       {"tolerance", required_argument, NULL, 10},
                                       {"no-jacobi", no_argument, NULL, 11},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?ve:i:o:O:l:r:R:s:S:t:T:W::a:b:x:", loptions, NULL)) >= 0) {
        switch (c) {
        case 'v':
            verbose++;
            break;
        case 1:
            if (filenames)
                error(
                    "Error: only one --ldgm-vcfs or --ldgm-vcfs-file inputs can be given, and they cannot be "
                    "combined\n");
            filenames = hts_readlist(optarg, 0, &n_files);
            break;
        case 2:
            if (filenames)
                error(
                    "Error: only one --ldgm-vcfs or --ldgm-vcfs-file inputs can be given, and they cannot be "
                    "combined\n");
            filenames = hts_readlist(optarg, 1, &n_files);
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
        case 'l':
            log_file = get_file_handle(optarg);
            break;
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 3:
            if (!strcasecmp(optarg, "0"))
                regions_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                regions_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                regions_overlap = 2;
            else
                error("Could not parse: --regions-overlap %s\n", optarg);
            break;
        case 's':
            sample_list = optarg;
            break;
        case 'S':
            sample_list = optarg;
            sample_is_file = 1;
            break;
        case 't':
            targets_list = optarg;
            break;
        case 'T':
            targets_list = optarg;
            targets_is_file = 1;
            break;
        case 4:
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
            n_threads = (int)strtol(optarg, &tmp, 10);
            if (*tmp) error("Could not parse: --threads %s\n", optarg);
            break;
        case 'W':
            if (!(write_index = write_index_parse(optarg))) error("Unsupported index format '%s'\n", optarg);
            break;
        case 5:
            stats_only = 1;
            break;
        case 6:
            average_ld_score = strtod(optarg, &tmp);
            if (*tmp) error("Could not parse: --average-ld-score %s\n", optarg);
            break;
        case 'a':
            alpha_param = strtod(optarg, &tmp);
            if (*tmp) error("Could not parse: --alpha-param %s\n", optarg);
            if (alpha_param < -1 || alpha_param > 0) error("The alpha parameter must be a number between -1 and 0\n");
            break;
        case 'b':
            beta_cov = strtod(optarg, &tmp);
            if (*tmp) error("Could not parse: --beta-cov %s\n", optarg);
            break;
        case 'x':
            cross_corr = strtod(optarg, &tmp);
            if (*tmp) error("Could not parse: --cross-corr %s\n", optarg);
            break;
        case 7:
            sample_sizes = hts_readlist(optarg, 0, &n_sample_sizes);
            break;
        case 10:
            tol = strtod(optarg, &tmp);
            if (*tmp) error("Could not parse: --tolerance %s\n", optarg);
            break;
        case 11:
            jacobi = 0;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }

    if (optind == argc) // no inputs provided
        error("%s", usage_text());

    if (filenames && optind + 1 != argc) // LDGM-VCFs input in two different ways
        error("Cannot input LDGM-VCFs in more than one way\n");

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

    if (n_files == 0) {
        n_files = argc - optind - 1;
        if (n_files == 0) error("No LDGM-VCF file was input\n");
        filenames = argv + optind + 1;
    }

    if (n_sample_sizes && n_sample_sizes != n_files)
        error("List of sample sizes does not match the number of LDGM-VCF files\n");

    if (!bcf_sr_add_reader(sr, argv[optind]))
        error("Error opening GWAS-VCF file %s: %s\n", argv[optind], bcf_sr_strerror(sr->errnum));
    check_gwas(bcf_sr_get_header(sr, 0), n_sample_sizes);

    for (i = 0; i < n_files; i++) {
        if (!bcf_sr_add_reader(sr, filenames[i]))
            error("Error opening LDGM-VCF file %s: %s\n", filenames[i], bcf_sr_strerror(sr->errnum));
        check_ldgm(bcf_sr_get_header(sr, 1 + i));
    }

    if (filenames != argv + optind + 1) {
        for (pop = 0; pop < n_files; pop++) free(filenames[pop]);
        free(filenames);
    }

    if (!stats_only && isnan(beta_cov)) {
        fprintf(log_file, "Warning: --beta-cov option should be specified rather than using default value\n");
        beta_cov = 1e-7;
        fprintf(log_file, "Warning: use option --stats-only to first identify value for option --beta-cov\n");
    }

    char wmode[8];
    set_wmode(wmode, output_type, output_fname, clevel);
    htsFile *out_fh = NULL;
    if (!stats_only) {
        out_fh = hts_open(output_fname, wmode);
        if (out_fh == NULL) error("Error: cannot write to \"%s\": %s\n", output_fname, strerror(errno));
        if (n_threads) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
    }

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if (bcf_hdr_nsamples(hdr) < n_files)
        error_errno("GWAS-VCF header file has only %d samples while %d required\n", bcf_hdr_nsamples(hdr), n_files);
    if (filter_str) filter = filter_init(hdr, filter_str);

    // subset input GWAS-VCF file to required summary statistics only
    char **samples;
    int *imap = (int *)malloc(n_files * sizeof(int));
    if (sample_list) {
        int n_samples;
        samples = hts_readlist(sample_list, sample_is_file, &n_samples);
        if (!samples) error("Could not read the list: \"%s\"\n", sample_list);
        if (n_samples != n_files)
            error("List of summary statistics has only %d samples while %d required\n", n_samples, n_files);
    } else {
        samples = hdr->samples;
    }
    for (pop = 0; pop < n_files; pop++) {
        imap[pop] = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, samples[pop]);
        if (imap[pop] < 0)
            error("Summary statistic %s not found in the GWAS-VCF file %s\n", samples[pop],
                  (bcf_sr_get_reader(sr, 0))->fname);
    }
    if (sample_list) {
        for (i = 0; i < n_files; i++) free(samples[i]);
        free(samples);
    }

    bcf_hdr_t *out_hdr = NULL;
    if (!stats_only) {
        out_hdr = bcf_hdr_subset(hdr, 0, 0, 0);
        bcf_hdr_remove(out_hdr, BCF_HL_FMT, NULL);
        kstring_t str = {0, 0, NULL};
        for (i = 0; i < n_files; i++) {
            str.l = 0;
            ksprintf(&str, n_files > 1 ? "%s_blupx_a%g_b%.2g" : "%s_blup_a%g_b%.2g", hdr->samples[imap[i]],
                     0.0 - alpha_param, beta_cov);
            bcf_hdr_add_sample(out_hdr, str.s);
        }
        free(str.s);

        if (bcf_hdr_sync(out_hdr) < 0)
            error_errno("[%s] Failed to update header", __func__); // updates the number of samples
        if (bcf_hdr_printf(out_hdr, "##FORMAT=<ID=%s,Number=A,Type=Float,Description=\"%s\">", id_str[ES], desc_str[ES])
            < 0)
            error_errno("[%s] Failed to add \"%s\" FORMAT header", id_str[ES], __func__);
        if (record_cmd_line) bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_blup");
        if (bcf_hdr_write(out_fh, out_hdr) < 0) error("Unable to write to output VCF file\n");
        if (init_index2(out_fh, hdr, output_fname, &index_fname, write_index) < 0)
            error("Error: failed to initialise index for %s\n", output_fname);
    }

    double *alpha_hat_1 = NULL;
    int m_alpha_hat_1 = 0;
    double *beta_hat_1 = NULL;
    int m_beta_hat_1 = 0;
    double *beta_hat = NULL;
    int m_beta_hat = 0;
    double *beta_blup = NULL;
    int m_beta_blup = 0;

    fprintf(log_file, "BLUP " BLUP_VERSION " http://github.com/freeseek/score BCFTOOLS %s HTSLIB %s\n",
            bcftools_version(), hts_version());

    if (!stats_only) {
        fprintf(log_file, "=== PARAMETERS ===\n");
        fprintf(log_file, "alpha: %.4g\n", alpha_param);
        fprintf(log_file, "betaCov: %.4g\n", beta_cov);
    }

    // allocate structures needed across ancestries
    ld_block_t *blocks = (ld_block_t *)calloc(sizeof(ld_block_t), n_files);
    for (pop = 0; pop < n_files; pop++) {
        blocks[pop].imap = imap[pop];
        blocks[pop].ld_block = -1;
    }
    free(imap);
    line_t *lines = NULL;
    int n_lines = 0;
    int m_lines = 0;

    if (sample_sizes) {
        for (pop = 0; pop < n_files; pop++) {
            blocks[pop].neff = strtod(sample_sizes[pop], &tmp);
            if (*tmp) error("Could not parse element: %s\n", sample_sizes[pop]);
            free(sample_sizes[pop]);
        }
        free(sample_sizes);
    }

    coo_matrix_t coo_S = {0};
    coo_matrix_t coo_P = {0};
    csr_matrix_t schur_P[2][2];
    csr_matrix_t S, S_P;

    // compute statistics across LD blocks
    int n_blocks = 0;
    double *medians_alpha_hat2 = NULL;
    int m_medians_alpha_hat2 = 0;
    double *means_neff = NULL;
    int m_means_neff = 0;

    if (!stats_only && verbose) fprintf(log_file, "=== LD_BLOCKS ===\n");

    // see http://github.com/awohns/ldgm/blob/main/MATLAB/BLUPxldgm.m
    int ret;
    do {
        n_lines = 0;
        ret = read_ld_block(sr, blocks, n_files, alpha_param, &lines, &n_lines, &m_lines, filter, filter_logic);
        if (stats_only || !verbose) fprintf(log_file, "\33[2K\r%s ld_block=%d", blocks[0].seqname, blocks[0].ld_block);

        int nrow = blocks[n_files - 1].row_ptr + blocks[n_files - 1].coo.nrow;
        hts_expand(double, nrow, m_alpha_hat_1, alpha_hat_1);
        hts_expand(double, nrow, m_beta_hat_1, beta_hat_1);
        hts_expand(double, nrow, m_beta_hat, beta_hat);
        hts_expand(double, nrow, m_beta_blup, beta_blup);
        hts_expand(double, (n_blocks + 1) * n_files, m_medians_alpha_hat2, medians_alpha_hat2);
        hts_expand(double, (n_blocks + 1) * n_files, m_means_neff, means_neff);

        for (pop = 0; pop < n_files; pop++) {
            ld_block_t *block = &blocks[pop];

            // create Schur complement
            schur_split(block, schur_P);

            // import data vector from LD block structure
            int l = 0;
            for (k = 0; k < block->coo.nrow; k++) {
                row_t *row = &block->rows[k];
                if (row->n_line == 0) continue;
                block->all_trace_inf += row->sd * row->sd;
                alpha_hat_1[block->row_ptr + l] = row->ez_deriv / sqrt(row->ne);
                l++;
            }
            medians_alpha_hat2[n_blocks * n_files + pop] = get_median2(&alpha_hat_1[block->row_ptr], l, 1);
            means_neff[n_blocks * n_files + pop] = block->mean_neff;
            if (stats_only) {
                for (k = 0; k < 4; k++) csr_destroy(&schur_P[k / 2][k % 2]);
                continue;
            }

            // run precisionMultiply()
            block->n_iter =
                precision_multiply(schur_P, &alpha_hat_1[block->row_ptr], tol, jacobi, &beta_hat_1[block->row_ptr]);
            for (k = 0; k < 4; k++) csr_destroy(&schur_P[k / 2][k % 2]);

            // divide precision matrix by n
            matdiv(&block->coo, block->mean_neff);

            if (verbose)
                fprintf(log_file, "%s neff=%.0f nnz=%d rows=%d missing=%d cg_multiply=%d\n", hdr->samples[block->imap],
                        block->mean_neff, block->coo.nnz, block->coo.nrow, block->n_missing, block->n_iter);
        }
        if (stats_only) {
            for (k = 0; k < n_lines; k++) bcf_destroy(lines[k].line);
            n_blocks++;
            continue;
        }

        // concatenate data vector
        for (pop = 0; pop < n_files; pop++) {
            ld_block_t *block = &blocks[pop];
            for (k = 0, l = 0; k < block->coo.nrow; k++) {
                row_t *row = &block->rows[k];
                if (row->n_line == 0)
                    beta_hat[block->row_ptr + k] = 0.0;
                else {
                    beta_hat[block->row_ptr + k] = beta_hat_1[block->row_ptr + l];
                    l++;
                }
            }
        }

        // create concatenated sigma
        make_sigma(blocks, n_files, beta_cov, cross_corr, &coo_S);
        coo_to_csr(&coo_S, &S);

        // create concatenated sigma + precision matrix
        concatenate(blocks, n_files, &coo_P);
        add_matrix(&coo_S, &coo_P);
        coo_to_csr(&coo_P, &S_P);

        // run precisionDivide()
        int n_iter = pcg(&S_P, beta_hat, tol, jacobi);
        csr_destroy(&S_P);

        // multiply by sigma
        sdmult(&S, beta_hat, beta_blup);
        csr_destroy(&S);

        // export data vector into LD block structure
        for (pop = 0; pop < n_files; pop++) {
            ld_block_t *block = &blocks[pop];
            for (k = 0; k < block->coo.nrow; k++) {
                row_t *row = &block->rows[k];
                if (row->n_line == 0) continue;
                row->ez_deriv = beta_blup[block->row_ptr + k] / row->sqrt_het;
            }
        }

        if (verbose)
            fprintf(log_file, "%s ld_block=%d rows=%d cg_divide=%d\n", blocks[0].seqname, blocks[0].ld_block,
                    coo_P.nrow, n_iter);

        // write BLUP GWAS-VCF files
        write_ld_block(out_fh, out_hdr, lines, n_lines, blocks, n_files);
        n_blocks++;
    } while (ret);

    fprintf(log_file, "\33[2K\r=== SUMMARY ===\n");

    // print trace(Sigma) estimate
    for (pop = 0; pop < n_files; pop++) {
        ld_block_t *block = &blocks[pop];
        double median_alpha_hat2 = get_median(&medians_alpha_hat2[pop], n_blocks, n_files);
        double sample_size = get_median(&means_neff[pop], n_blocks, n_files);
        double lambda_GC = sample_size * median_alpha_hat2 / MEDIAN_CHISQ;
        double proportion_non_missing =
            (double)block->all_n_non_missing / (double)(block->all_n_non_missing + block->all_n_missing);
        double sigmasq_inf = (lambda_GC - 1.0) / (sample_size * average_ld_score * proportion_non_missing);
        if (sigmasq_inf < 0.0) sigmasq_inf = 0.0;
        const char *fname = (bcf_sr_get_reader(sr, 1 + pop))->fname;
        fprintf(log_file, "%s %s non_missing=%d missing=%d lambda_GC=%.4f sigmasqInf=%.4g\n", hdr->samples[block->imap],
                strrchr(fname, '/') ? strrchr(fname, '/') + 1 : fname, block->all_n_non_missing, block->all_n_missing,
                lambda_GC, sigmasq_inf);
        if (!stats_only) fprintf(log_file, "Tr(S_inf)=%g\n", block->all_trace_inf * beta_cov);
        double meanPerSNPh2 = block->all_trace_inf / (double)block->all_n_non_missing;
        fprintf(log_file, "Advised options: --alpha-param %.4f --beta-cov %.4g\n", alpha_param,
                sigmasq_inf / meanPerSNPh2);
    }

    free(medians_alpha_hat2);
    free(means_neff);
    free(alpha_hat_1);
    free(beta_hat_1);
    free(beta_hat);
    free(beta_blup);
    for (pop = 0; pop < n_files; pop++) ld_block_destroy(&blocks[pop]);
    free(blocks);
    free(lines);
    coo_destroy(&coo_S);
    coo_destroy(&coo_P);

    if (filter) filter_destroy(filter);
    if (!stats_only) {
        if (write_index) {
            if (bcf_idx_save(out_fh) < 0) {
                if (hts_close(out_fh) != 0)
                    error("Close failed %s\n", strcmp(output_fname, "-") ? output_fname : "stdout");
                error("Error: cannot write to index %s\n", index_fname);
            }
            free(index_fname);
        }
        if (hts_close(out_fh) < 0) error("Close failed: %s\n", out_fh->fn);
        bcf_hdr_destroy(out_hdr);
    }
    bcf_sr_destroy(sr);
    if (log_file != stdout && log_file != stderr) fclose(log_file);

    return 0;
}
