/**
 * @file    radar_detection.c
 * @brief   Radar detection theory: CFAR, Pd/Pfa, Marcum Q, Albersheim equation
 */
#include "radar_detection.h"
#include "radar_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "radar_detection: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

int cfar_config_init(cfar_config_t *cfg, cfar_type_t type,
                     size_t guard_n, size_t ref_n, double pfa)
{
    if (!cfg || pfa <= 0.0 || pfa >= 1.0 || ref_n == 0) return -1;
    memset(cfg, 0, sizeof(*cfg));
    cfg->type = type;
    cfg->num_guard_cells = guard_n;
    cfg->num_ref_cells = ref_n;
    cfg->pfa_desired = pfa;
    cfg->use_log_detector = 0;
    size_t N = 2 * ref_n;
    if (pfa < 1e-15) pfa = 1e-15;
    cfg->threshold_factor = (double)N * (pow(pfa, -1.0/(double)N) - 1.0);
    cfg->os_rank = (size_t)(0.75 * (double)N);
    if (cfg->os_rank < 1) cfg->os_rank = 1;
    if (cfg->os_rank > N) cfg->os_rank = N;
    return 0;
}

int cfar_detect(const double *rp, size_t nbins,
                const cfar_config_t *cfg, int *dets)
{
    if (!rp || !cfg || !dets || nbins == 0) return -1;
    size_t g = cfg->num_guard_cells, nr = cfg->num_ref_cells;
    size_t w = g + nr;
    double alpha = cfg->threshold_factor;

    for (size_t i = 0; i < nbins; i++) {
        dets[i] = 0;
        int hl = (i >= w), hr = (i + w < nbins);
        if (!hl && !hr) continue;

        double nest = 0.0; size_t nu = 0;

        if (cfg->type == CFAR_GO || cfg->type == CFAR_SO) {
            double lm = 0.0, rm = 0.0; size_t nl = 0, nrr = 0;
            if (hl) {
                for (size_t j = i - g - nr; j + g < i; j++)
                { lm += rp[j]; nl++; }
                if (nl > 0) lm /= (double)nl;
            }
            if (hr) {
                for (size_t j = i + g + 1;
                     j <= i + g + nr && j < nbins; j++)
                { rm += rp[j]; nrr++; }
                if (nrr > 0) rm /= (double)nrr;
            }
            if (cfg->type == CFAR_GO) {
                if (hl && hr) nest = (lm > rm) ? lm : rm;
                else if (hl) nest = lm; else nest = rm;
            } else {
                if (hl && hr) nest = (lm < rm) ? lm : rm;
                else if (hl) nest = lm; else nest = rm;
            }
        } else {
            if (hl) {
                for (size_t j = i - g - nr; j + g < i; j++)
                { nest += rp[j]; nu++; }
            }
            if (hr) {
                for (size_t j = i + g + 1;
                     j <= i + g + nr && j < nbins; j++)
                { nest += rp[j]; nu++; }
            }
            if (nu > 0) nest /= (double)nu;
        }

        double thresh = alpha * nest;
        if (thresh > 0.0 && rp[i] > thresh) dets[i] = 1;
    }
    return 0;
}

int cfar_os_detect(const double *rp, size_t nbins,
                   const cfar_config_t *cfg, int *dets)
{
    if (!rp || !cfg || !dets || nbins == 0) return -1;
    size_t g = cfg->num_guard_cells, nr = cfg->num_ref_cells;
    size_t total_ref = 2 * nr;
    size_t rank = cfg->os_rank;
    double alpha_os = cfg->threshold_factor;
    {
        double pfa = cfg->pfa_desired;
        double lo = 0.0, hi = 100.0;
        for (int iter = 0; iter < 50; iter++) {
            double mid = (lo + hi) / 2.0;
            double prod = 1.0;
            for (size_t j = rank; j <= total_ref; j++)
                prod *= (double)j / ((double)j + mid);
            if (prod > pfa) lo = mid; else hi = mid;
        }
        alpha_os = (lo + hi) / 2.0;
    }

    for (size_t i = 0; i < nbins; i++) {
        dets[i] = 0;
        size_t w = g + nr;
        int hl = (i >= w), hr = (i + w < nbins);
        if (!hl || !hr) continue;

        double *ref = safe_malloc(total_ref * sizeof(double));
        size_t nc = 0;
        for (size_t j = i - g - nr; j + g < i; j++)
            ref[nc++] = rp[j];
        for (size_t j = i + g + 1;
             j <= i + g + nr && j < nbins; j++)
            ref[nc++] = rp[j];

        radar_dsort(ref, nc);
        size_t kidx = rank;
        if (kidx >= nc) kidx = nc - 1;
        double nest = ref[kidx];
        free(ref);

        double thresh = alpha_os * nest;
        if (thresh > 0.0 && rp[i] > thresh) dets[i] = 1;
    }
    return 0;
}

double radar_detection_threshold(double pfa, size_t n_pulses)
{
    if (pfa <= 0.0 || pfa >= 1.0 || n_pulses == 0) return 0.0;
    double N = (double)n_pulses;
    double target = 1.0 - pfa;
    double lo = 0.0, hi = 100.0 + 20.0 * N;
    for (int iter = 0; iter < 60; iter++) {
        double mid = (lo + hi) / 2.0;
        double val = radar_gamma_inc(N, mid / 2.0);
        if (val < target) lo = mid; else hi = mid;
    }
    return (lo + hi) / 2.0;
}

double radar_marcum_q(int m, double a, double b)
{
    if (m < 1 || a < 0.0 || b < 0.0) return 0.0;
    if (b == 0.0) return 1.0;
    if (a == 0.0) {
        double sum = 1.0, term = 1.0;
        for (int k = 1; k < m; k++) {
            term *= (b * b) / (2.0 * (double)k);
            sum += term;
        }
        return sum * exp(-b * b / 2.0);
    }
    double s = a * b;
    double q1;
    if (s < 30.0) {
        double term = 1.0, sum = 1.0;
        for (int k = 1; k < 100; k++) {
            term *= s / (double)k;
            sum += term;
            if (term / sum < 1e-12) break;
        }
        double i0 = sum * exp(-s);
        sum = 0.0; term = 1.0;
        for (int k = 0; k < 100; k++) {
            term *= (s / 2.0) / (double)(k + 1);
            sum += term;
            if (term / (sum + 1e-30) < 1e-12) break;
        }
        double i1 = (s / 2.0) * sum * exp(-s);
        double t = b / a;
        q1 = exp(-(a*a + b*b) / 2.0) *
             (i0 + t*i1 + t*t*(i0 + (2.0/s)*i1)/2.0);
        if (q1 < 0.0) q1 = 0.0;
        if (q1 > 1.0) q1 = 1.0;
    } else {
        double d = b - a;
        double sigma2 = (a + b) / (2.0 * s);
        q1 = 0.5 * erfc(d / sqrt(2.0 * sigma2));
    }
    if (m == 1) return q1;
    double q = q1;
    for (int k = 2; k <= m; k++) {
        double ratio = pow(b / a, (double)(k - 1));
        double exp_term = exp(-(a*a + b*b) / 2.0);
        double ik = (s < 20.0)
            ? pow(s/2.0, (double)(k-1)) * exp(-s)
            : exp(s) / sqrt(2.0 * M_PI * s);
        q += ratio * exp_term * ik;
        if (q > 1.0) q = 1.0;
    }
    return q;
}

double radar_pd_marcum(double snr_linear, double pfa, size_t n_pulses)
{
    if (snr_linear <= 0.0 || pfa <= 0.0 || pfa >= 1.0 || n_pulses == 0) return 0.0;
    double threshold = radar_detection_threshold(pfa, n_pulses);
    double a = sqrt(2.0 * (double)n_pulses * snr_linear);
    double b = sqrt(threshold);
    return radar_marcum_q((int)n_pulses, a, b);
}

double radar_pd_swerling1(double snr_linear, double pfa)
{
    if (snr_linear <= 0.0 || pfa <= 0.0 || pfa >= 1.0) return 0.0;
    return pow(pfa, 1.0 / (1.0 + snr_linear));
}

double radar_pd_swerling2(double snr_linear, double pfa, size_t n_pulses)
{
    if (snr_linear <= 0.0 || pfa <= 0.0 || pfa >= 1.0 || n_pulses == 0) return 0.0;
    double threshold = radar_detection_threshold(pfa, n_pulses);
    double x = threshold / (2.0 * (1.0 + snr_linear));
    double ginc = radar_gamma_inc((double)n_pulses, x);
    return 1.0 - ginc;
}

double radar_albersheim_snr(double pd, double pfa, size_t n_pulses)
{
    if (pd <= 0.0 || pd >= 1.0 || pfa <= 0.0 || pfa >= 1.0 || n_pulses == 0) return 0.0;
    double N = (double)n_pulses;
    double A = log(0.62 / pfa);
    double B = log(pd / (1.0 - pd));
    double term1 = -5.0 * log10(N);
    double term2 = 6.2 + 4.54 / sqrt(N + 0.44);
    double term3 = log10(A + 0.12 * A * B + 1.7 * B);
    return term1 + term2 * term3;
}

int radar_roc_curve(double snr_linear, double *pfa_arr, double *pd_arr,
                    size_t n_points, size_t n_pulses, int swerling_case)
{
    if (!pfa_arr || !pd_arr || n_points == 0 || snr_linear <= 0.0 || n_pulses == 0)
        return -1;
    double pfa_min = 1e-6, pfa_max = 0.5;
    for (size_t i = 0; i < n_points; i++) {
        double t = (double)i / (double)(n_points - 1);
        pfa_arr[i] = pfa_min * pow(pfa_max / pfa_min, t);
    }
    for (size_t i = 0; i < n_points; i++) {
        switch (swerling_case) {
        case 0: pd_arr[i] = radar_pd_marcum(snr_linear, pfa_arr[i], n_pulses); break;
        case 1: pd_arr[i] = radar_pd_swerling1(snr_linear, pfa_arr[i]); break;
        case 2: pd_arr[i] = radar_pd_swerling2(snr_linear, pfa_arr[i], n_pulses); break;
        default: pd_arr[i] = radar_pd_marcum(snr_linear, pfa_arr[i], n_pulses); break;
        }
    }
    return 0;
}

double radar_gamma_inc(double a, double x)
{
    if (a <= 0.0 || x < 0.0) return 0.0;
    if (x == 0.0) return 0.0;
    if (x > 700.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 0; n < 100; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * 1e-12) break;
        }
        return sum * exp(-x + a * log(x) - lgamma(a));
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / 1e-30;
        double d = 1.0 / b;
        double h = d;
        for (int n = 1; n < 100; n++) {
            double an = -(double)n * ((double)n - a);
            b += 2.0;
            d = an * d + b;
            if (fabs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (fabs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (fabs(del - 1.0) < 1e-12) break;
        }
        return 1.0 - exp(-x + a * log(x) - lgamma(a)) * h;
    }
}

static int dbl_cmp(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

void radar_dsort(double *arr, size_t n)
{
    if (arr && n > 1) qsort(arr, n, sizeof(double), dbl_cmp);
}

/* L7: Shnidman's equation for required SNR (extended Albersheim) */
double radar_shnidman_snr(double pd, double pfa, size_t n_pulses, int swerling)
{
    if (pd <= 0.0 || pd >= 1.0 || pfa <= 0.0 || pfa >= 1.0 || n_pulses == 0)
        return 0.0;
    double snr_albersheim = radar_albersheim_snr(pd, pfa, n_pulses);
    /* Shnidman correction factors for Swerling cases */
    double C_sw[5] = {0.0, 0.5, 1.0, -0.5, -1.0};
    double C_n[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    int idx = (swerling >= 0 && swerling <= 4) ? swerling : 0;
    return snr_albersheim + C_sw[idx] + C_n[idx] / sqrt((double)n_pulses);
}

/* L7: Compute optimal detection threshold via Neyman-Pearson */
double radar_optimal_threshold(double pfa, double noise_power)
{
    if (pfa <= 0.0 || pfa >= 1.0 || noise_power <= 0.0) return 0.0;
    /* For square-law detector in AWGN: T = -noise_power * ln(Pfa) */
    return -noise_power * log(pfa);
}

/* L7: Binary integration (M-of-N detection) */
double radar_m_of_n_pd(double single_pd, size_t m, size_t n)
{
    if (m > n || n == 0) return 0.0;
    /* P(M of N) = sum_{k=m}^{N} C(N,k) * Pd^k * (1-Pd)^{N-k} */
    double total = 0.0;
    for (size_t k = m; k <= n; k++) {
        /* Compute binomial coefficient C(n,k) */
        double binom = 1.0;
        for (size_t j = 1; j <= k; j++)
            binom *= (double)(n - k + j) / (double)j;
        total += binom * pow(single_pd, (double)k) *
                 pow(1.0 - single_pd, (double)(n - k));
    }
    return total;
}

/* L7: Compute CFAR loss in dB */
double radar_cfar_loss_db(size_t n_ref_cells, double pfa)
{
    if (n_ref_cells == 0 || pfa <= 0.0 || pfa >= 1.0) return 0.0;
    double N = (double)n_ref_cells;
    /* CFAR loss = 1 + (1/N)*(Pfa^{-1/N} - 1)/(ln(Pfa^{-1/N})) */
    double pfa_inv_N = pow(pfa, -1.0 / N);
    double loss_lin = 1.0 + (pfa_inv_N - 1.0) / (N * log(pfa_inv_N));
    return 10.0 * log10(loss_lin);
}

/* L7: Detection range at specific Pd/Pfa */
double radar_detection_range(const radar_params_t *params,
                              double rcs_m2, double pd, double pfa,
                              size_t n_pulses, int swerling)
{
    if (!params) return 0.0;
    double snr_req_db = radar_shnidman_snr(pd, pfa, n_pulses, swerling);
    double snr_req = db2lin(snr_req_db);
    return radar_max_range(params, rcs_m2, snr_req);
}

/* L7: Integration loss due to non-ideal detector */
double radar_detector_loss_db(const char *detector_type)
{
    /* Approximate losses relative to ideal matched filter */
    if (!detector_type) return 0.0;
    if (strcmp(detector_type, "square-law") == 0) return 0.0;
    if (strcmp(detector_type, "linear") == 0)     return 0.5;
    if (strcmp(detector_type, "log") == 0)        return 1.0;
    return 0.0;
}

/* L8: Adaptive detection threshold using local statistics */
void radar_adaptive_threshold(const double *range_profile, size_t n_bins,
                               size_t window_size, double pfa,
                               double *thresholds)
{
    if (!range_profile || !thresholds || window_size < 2) return;
    for (size_t i = 0; i < n_bins; i++) {
        double sum = 0.0, sumsq = 0.0;
        size_t count = 0;
        size_t start = (i >= window_size / 2) ? i - window_size / 2 : 0;
        size_t end = (i + window_size / 2 < n_bins) ?
                      i + window_size / 2 : n_bins - 1;
        for (size_t j = start; j <= end; j++) {
            sum += range_profile[j];
            sumsq += range_profile[j] * range_profile[j];
            count++;
        }
        double mean = (count > 0) ? sum / (double)count : 0.0;
        double var = (count > 1) ? (sumsq - sum * mean) / (double)(count - 1) : 0.0;
        if (var < 0.0) var = 0.0;
        thresholds[i] = mean + sqrt(var) * sqrt(-2.0 * log(pfa));
    }
}
