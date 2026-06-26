/**
 * @file    ir_detection.c
 * @brief   Infrared Target Detection - CFAR, Thresholding, Matched Filter
 *
 * Implements fundamental IR target detection algorithms:
 * CA-CFAR (Cell-Averaging Constant False Alarm Rate),
 * background subtraction, Otsu thresholding, matched filter,
 * and Signal-to-Clutter Ratio (SCR) estimation.
 *
 * References:
 *   Richards, Scheer & Holm (2010) "Principles of Modern Radar", Ch.14
 *   Kay, S.M. (1998) "Fundamentals of Statistical Signal Processing:
 *     Detection Theory", Prentice Hall
 *   Otsu, N. (1979) "A Threshold Selection Method from Gray-Level
 *     Histograms", IEEE Trans. SMC, 9(1), 62-66
 */

#include "ir_core.h"
#include "ir_image.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ========================================================================
 * CA-CFAR (Cell-Averaging Constant False Alarm Rate) Detector
 *
 * Used extensively in IRST (Infrared Search and Track) and missile warning
 * systems. Adaptively sets detection threshold based on local background.
 *
 * Threshold = mu_bg * alpha  where alpha = N * (P_fa^(-1/N) - 1)
 * ======================================================================== */

int ir_cfar_detect_1d(const double *signal, int length,
                       int test_idx, int guard_cells,
                       int reference_cells, double p_fa,
                       double *threshold_out) {
    if (!signal || length <= 0 || test_idx < 0 || test_idx >= length)
        return -1;
    if (guard_cells < 0 || reference_cells <= 0 || p_fa <= 0.0 || p_fa >= 1.0)
        return -1;

    double sum = 0.0;
    int count = 0;

    /* Left reference window */
    for (int i = test_idx - guard_cells - reference_cells;
         i < test_idx - guard_cells; i++) {
        if (i >= 0 && i < length) { sum += signal[i]; count++; }
    }
    /* Right reference window */
    for (int i = test_idx + guard_cells + 1;
         i <= test_idx + guard_cells + reference_cells; i++) {
        if (i >= 0 && i < length) { sum += signal[i]; count++; }
    }

    if (count == 0) return 0;

    double mu_bg = sum / count;
    double alpha = count * (pow(p_fa, -1.0 / count) - 1.0);
    double threshold = mu_bg * alpha;

    if (threshold_out) *threshold_out = threshold;
    return (signal[test_idx] > threshold) ? 1 : 0;
}

/* ========================================================================
 * Background Subtraction
 * ======================================================================== */

int ir_background_subtract(const double *frame, const double *background,
                            int length, double *diff_out) {
    if (!frame || !background || !diff_out || length <= 0) return -1;
    for (int i = 0; i < length; i++)
        diff_out[i] = fabs(frame[i] - background[i]);
    return 0;
}

int ir_background_update(const double *frame, double *background,
                          int length, double alpha) {
    if (!frame || !background || length <= 0) return -1;
    if (alpha < 0.0 || alpha > 1.0) return -1;
    for (int i = 0; i < length; i++)
        background[i] = alpha * frame[i] + (1.0 - alpha) * background[i];
    return 0;
}

int ir_threshold_segment(const double *data, int length,
                          double threshold, int *binary_out) {
    if (!data || !binary_out || length <= 0) return -1;
    for (int i = 0; i < length; i++)
        binary_out[i] = (data[i] > threshold) ? 1 : 0;
    return 0;
}

/* ========================================================================
 * Otsu's Automatic Thresholding (1979)
 *
 * Maximizes between-class variance to find optimal binary threshold.
 * Widely used in IR image segmentation for automatic target detection.
 * ======================================================================== */

double ir_otsu_threshold(const double *data, int length, int n_bins) {
    if (!data || length <= 0 || n_bins < 2) return -1.0;

    double d_min = data[0], d_max = data[0];
    for (int i = 1; i < length; i++) {
        if (data[i] < d_min) d_min = data[i];
        if (data[i] > d_max) d_max = data[i];
    }
    if (fabs(d_max - d_min) < 1e-10) return d_min;

    int *hist = (int*)calloc(n_bins, sizeof(int));
    for (int i = 0; i < length; i++) {
        int bin = (int)((data[i] - d_min) / (d_max - d_min) * (n_bins - 1));
        if (bin < 0) bin = 0;
        if (bin >= n_bins) bin = n_bins - 1;
        hist[bin]++;
    }

    double total = length;
    double sum_all = 0.0;
    for (int i = 0; i < n_bins; i++) sum_all += i * hist[i];

    double best_thresh = 0.0, best_var = 0.0;
    double sum_b = 0.0, w_b = 0.0;

    for (int t = 0; t < n_bins; t++) {
        w_b += hist[t];
        if (w_b == 0.0) continue;
        double w_f = total - w_b;
        if (w_f == 0.0) break;
        sum_b += t * hist[t];
        double m_b = sum_b / w_b;
        double m_f = (sum_all - sum_b) / w_f;
        double var_between = w_b * w_f * (m_b - m_f) * (m_b - m_f);
        if (var_between > best_var) {
            best_var = var_between;
            best_thresh = d_min + (d_max - d_min) * t / (n_bins - 1.0);
        }
    }
    free(hist);
    return best_thresh;
}

/* ========================================================================
 * Matched Filter (Correlation Detector)
 * ======================================================================== */

double ir_matched_filter_1d(const double *signal, int sig_len,
                             const double *template, int temp_len,
                             double *correlation_out) {
    if (!signal || !template || sig_len <= 0 || temp_len <= 0) return -1.0;
    if (temp_len > sig_len) return -1.0;

    double max_corr = -1e308;
    for (int i = 0; i <= sig_len - temp_len; i++) {
        double corr = 0.0;
        for (int j = 0; j < temp_len; j++)
            corr += signal[i + j] * template[j];
        if (correlation_out) correlation_out[i] = corr;
        if (corr > max_corr) max_corr = corr;
    }
    return max_corr;
}

/* ========================================================================
 * Signal-to-Clutter Ratio (SCR) Estimation
 *
 * SCR = (mean_target - mean_background) / stddev_background
 *
 * SCR is a key performance metric for IR detection systems.
 * SCR > 5 is typically required for reliable detection.
 * ======================================================================== */

double ir_scr_estimate(const double *data, int length,
                        int target_start, int target_end) {
    if (!data || length <= 0) return -1.0;
    if (target_start < 0 || target_end > length || target_start >= target_end)
        return -1.0;

    double mu_t = 0.0;
    for (int i = target_start; i < target_end; i++) mu_t += data[i];
    mu_t /= (target_end - target_start);

    double mu_bg = 0.0, var_bg = 0.0;
    int n_bg = 0;
    for (int i = 0; i < target_start; i++) { mu_bg += data[i]; n_bg++; }
    for (int i = target_end; i < length; i++) { mu_bg += data[i]; n_bg++; }
    if (n_bg == 0) return -1.0;
    mu_bg /= n_bg;

    for (int i = 0; i < target_start; i++)
        var_bg += (data[i] - mu_bg) * (data[i] - mu_bg);
    for (int i = target_end; i < length; i++)
        var_bg += (data[i] - mu_bg) * (data[i] - mu_bg);
    double sigma_bg = sqrt(var_bg / n_bg);

    if (sigma_bg < 1e-10) return (mu_t > mu_bg) ? 1e308 : 0.0;
    return (mu_t - mu_bg) / sigma_bg;
}

/* ========================================================================
 * Detection Probability Estimation
 *
 * P_d = Q(Q^{-1}(P_fa) - sqrt(SCR))
 * where Q is the Q-function (complementary CDF of standard normal).
 * Simple approximation using erf for demonstration.
 * ======================================================================== */

double ir_detection_probability(double scr, double p_fa) {
    if (scr < 0.0 || p_fa <= 0.0 || p_fa >= 1.0) return -1.0;

    /* Q^{-1}(x) approx: sqrt(2) * erfinv(1 - 2x) */
    /* Simplified using erf/erfc */
    (void)p_fa; /* unused in simplified approx */
    /* Actually use the approximation via the inverse error function */
    /* Q^{-1}(P_fa) ~ sqrt(-2*ln(P_fa)) for P_fa << 1 */
    double arg = sqrt(scr) - sqrt(-2.0 * log(p_fa));
    /* P_d = Q(-arg) = 1 - Phi(-arg) = Phi(arg) where arg is the SNR margin */
    return 0.5 * erfc(-arg / sqrt(2.0));
}
