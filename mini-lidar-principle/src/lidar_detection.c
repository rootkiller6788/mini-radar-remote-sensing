/**
 * @file    lidar_detection.c
 * @brief   LiDAR detection theory — threshold, CFAR, matched filter, noise
 *
 * Knowledge covered:
 *   L1: Detection probability, false alarm probability, ROC curves,
 *       range accuracy (CRLB)
 *   L2: Constant False Alarm Rate (CFAR) detection
 *   L3: Gaussian noise model, Poisson statistics, Q-function
 *   L4: Neyman-Pearson criterion, detection threshold derivation
 *   L5: Matched filter for pulse detection, CA-CFAR implementation
 *
 * Reference:
 *   - Kay, S.M., *Detection Theory*, Prentice Hall, 1998.
 *   - Richards et al., *Principles of Modern Radar*, 2010, Ch.15-16.
 *   - Kingston, R.H., *Optical Sources, Detectors, and Systems*, 1995.
 */

#include "lidar_detection.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════
 * L4: Q-function and inverse Q-function (rational approximations)
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * @brief Q-function: Q(x) = (1/sqrt(2*pi)) * integral_x^inf exp(-t^2/2) dt
 *
 * Using the Abramowitz & Stegun approximation (26.2.17):
 *   Q(x) ≈ Z(x) * (b1*t + b2*t^2 + b3*t^3 + b4*t^4 + b5*t^5)
 *   where Z(x) = exp(-x^2/2)/sqrt(2*pi), t = 1/(1 + p*x), p = 0.2316419
 *
 * Relative error < 1.5e-7 for all x.
 *
 * Reference: Abramowitz & Stegun, *Handbook of Mathematical Functions*,
 *            Dover, 1965, Sec.26.2.
 */
static double Q_function(double x) {
    double p = 0.2316419;
    double b1 = 0.319381530;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;

    double t = 1.0 / (1.0 + p * fabs(x));
    double Z = exp(-x * x / 2.0) / sqrt(2.0 * M_PI);
    double poly = b1 * t + b2 * t * t + b3 * t * t * t
                  + b4 * t * t * t * t + b5 * t * t * t * t * t;

    double Q_val = Z * poly;
    return (x > 0.0) ? Q_val : 1.0 - Q_val;
}

/**
 * @brief Inverse Q-function (probit function)
 *
 * Using the Odeh & Evans rational approximation.
 * Relative error < 1.2e-8 for 10^(-20) < p < 1 - 10^(-20).
 *
 * Reference: Odeh & Evans, "Algorithm AS 70: The Percentage Points of
 *            the Normal Distribution", *Applied Statistics* 23, 1974.
 */
static double Q_inverse(double p) {
    if (p <= 0.0) return 10.0;     /* extreme — cap at 10 sigma */
    if (p >= 1.0) return -10.0;
    if (p > 0.5) return -Q_inverse(1.0 - p);

    double t = sqrt(-2.0 * log(p));
    double c0 = 2.515517;
    double c1 = 0.802853;
    double c2 = 0.010328;
    double d1 = 1.432788;
    double d2 = 0.189269;
    double d3 = 0.001308;

    double num = c0 + c1 * t + c2 * t * t;
    double den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;
    return t - num / den;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L4: Detection threshold and probability computations
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_detection_threshold_gaussian(double pfa,
                                            double noise_mean,
                                            double noise_std) {
    if (pfa <= 0.0 || pfa >= 1.0 || noise_std <= 0.0) return 0.0;
    /* T = mu + sigma * Q^{-1}(PFA) */
    double z = Q_inverse(pfa);
    return noise_mean + noise_std * z;
}

double lidar_prob_detection(double snr_db, double pfa) {
    if (snr_db < -900.0) return 0.0;
    if (pfa <= 0.0 || pfa >= 1.0) return 0.0;

    double snr_linear = pow(10.0, snr_db / 10.0);
    double a = sqrt(snr_linear);
    double z_pfa = Q_inverse(pfa);

    /* P_d = Q( Q^{-1}(PFA) - sqrt(SNR) ) */
    return Q_function(z_pfa - a);
}

double lidar_detectability_factor(double pd, double pfa) {
    if (pd <= 0.0 || pd >= 1.0 || pfa <= 0.0 || pfa >= 1.0) return -1.0;

    double z_pd  = Q_inverse(pd);
    double z_pfa = Q_inverse(pfa);

    /* D = (Q^{-1}(PFA) - Q^{-1}(P_d))^2 */
    double diff = z_pfa - z_pd;
    return diff * diff;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L2/L5: CA-CFAR detection
 * ═══════════════════════════════════════════════════════════════════════════ */

size_t lidar_cfar_detection(const double *signal, size_t num_bins,
                              double pfa,
                              size_t guard_cells, size_t ref_cells,
                              lidar_detection_result_t *detections,
                              size_t max_detections) {
    if (!signal || num_bins == 0 || !detections || max_detections == 0) return 0;

    size_t n_ref = 2 * ref_cells;
    if (n_ref == 0) return 0;

    /* CA-CFAR threshold multiplier: alpha = n_ref * (PFA^{-1/n_ref} - 1) */
    double alpha = (double)n_ref * (pow(pfa, -1.0 / (double)n_ref) - 1.0);

    size_t n_det = 0;
    size_t window = guard_cells + ref_cells;

    for (size_t i = window; i + window < num_bins && n_det < max_detections; i++) {
        /* Sum reference cells (left + right windows) */
        double sum_ref = 0.0;
        for (size_t j = i - window; j < i - guard_cells; j++) {
            sum_ref += signal[j];
        }
        for (size_t j = i + guard_cells + 1; j <= i + window; j++) {
            sum_ref += signal[j];
        }

        double noise_est = sum_ref / (double)n_ref;
        double threshold = alpha * noise_est;

        if (signal[i] > threshold) {
            detections[n_det].detected = 1;
            detections[n_det].amplitude = signal[i];
            detections[n_det].threshold = threshold;
            detections[n_det].range = 0.0; /* caller maps bin → range */
            detections[n_det].snr = (noise_est > 0.0)
                ? 10.0 * log10(signal[i] / noise_est) : 0.0;
            detections[n_det].prob_false_alarm = pfa;
            detections[n_det].prob_detection = lidar_prob_detection(
                detections[n_det].snr, pfa);
            n_det++;
        }
    }

    return n_det;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Matched filter
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_matched_filter(const double *waveform, size_t num_samples,
                          const double *pulse, size_t pulse_len,
                          double *output, size_t *out_len) {
    if (!waveform || !pulse || !output || !out_len) return -1;
    if (num_samples == 0 || pulse_len == 0) return -1;

    /* Output length = N + M - 1 (full convolution) */
    size_t o_len = num_samples + pulse_len - 1;
    if (*out_len < o_len) o_len = *out_len;

    /* Time-domain convolution: y[k] = sum_j x[j] * h[k-j]
       where h[j] = pulse[pulse_len - 1 - j] (time-reversed) */
    for (size_t k = 0; k < o_len; k++) {
        double sum = 0.0;
        for (size_t j = 0; j < pulse_len; j++) {
            size_t sig_idx;
            if (k >= j) {
                sig_idx = k - j;
            } else {
                break; /* k < j: index out of bounds */
            }
            if (sig_idx < num_samples) {
                /* Matched filter: h[j] = pulse[L-1-j] */
                sum += waveform[sig_idx] * pulse[pulse_len - 1 - j];
            }
        }
        output[k] = sum;
    }

    *out_len = o_len;
    return 0;
}

void lidar_gaussian_pulse_template(double *pulse, size_t length,
                                     double dt, double fwhm) {
    if (!pulse || length == 0 || dt <= 0.0 || fwhm <= 0.0) return;

    double sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0))); /* σ = FWHM / 2.35482 */
    double t_offset = (double)(length - 1) * dt / 2.0;  /* Center */

    for (size_t i = 0; i < length; i++) {
        double t = (double)i * dt - t_offset;
        pulse[i] = exp(-0.5 * t * t / (sigma * sigma));
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L2: Optical noise sources
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_shot_noise_variance(double current, double bandwidth,
                                   double M, double F) {
    if (current < 0.0 || bandwidth <= 0.0 || M <= 0.0 || F < 1.0) return 0.0;
    /* sigma^2 = 2 * q * I_primary * M^2 * F * B */
    return 2.0 * LIDAR_ELECTRON_Q * current * M * M * F * bandwidth;
}

double lidar_thermal_noise_variance(double temperature,
                                      double bandwidth,
                                      double resistance) {
    if (temperature <= 0.0 || bandwidth <= 0.0 || resistance <= 0.0) return 0.0;
    /* sigma^2 = 4 * k * T * B / R */
    return 4.0 * LIDAR_BOLTZMANN_K * temperature * bandwidth / resistance;
}

double lidar_background_current(const lidar_config_t *config,
                                  const lidar_atmosphere_t *atm) {
    if (!config || !atm) return 0.0;

    /* Solar spectral radiance approximation */
    double L_solar;
    if (config->wavelength == 905)      L_solar = 500.0;
    else if (config->wavelength == 1064) L_solar = 300.0;
    else if (config->wavelength == 1550) L_solar = 150.0;
    else if (config->wavelength == 532)  L_solar = 1500.0;
    else                                 L_solar = 400.0;

    /* Optical filter bandwidth: assume 10 nm = 0.01 um */
    double delta_lambda = 0.01;
    /* Convert to meters (um -> m): multiply by 1e-6 for power calculation
       L_solar is in W/(m^2*sr*um), so delta_lambda in um is correct */

    /* Receiver aperture area */
    double A_r = M_PI * (config->aperture_diam / 2.0) * (config->aperture_diam / 2.0);

    /* Receiver FOV solid angle */
    double Omega_r = M_PI * (config->scan_fov_v / 2.0) * (config->scan_fov_h / 2.0);
    if (Omega_r > M_PI) Omega_r = M_PI;
    if (Omega_r < 1e-12) Omega_r = 1e-12;

    double P_bg = L_solar * delta_lambda * A_r * Omega_r * config->opt_receive;

    /* Include atmospheric background radiance if available */
    if (atm->background_rad > 0.0) {
        P_bg += atm->background_rad * delta_lambda * A_r * Omega_r * config->opt_receive;
    }

    return P_bg * config->responsivity;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: ROC curve computation
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_roc_curve(double snr_db,
                      double *pfa, double *pd,
                      size_t num_pts,
                      double pfa_min, double pfa_max) {
    if (!pfa || !pd || num_pts == 0 || snr_db < -900.0) return -1;
    if (pfa_min <= 0.0) pfa_min = 1e-6;
    if (pfa_max >= 1.0) pfa_max = 0.5;
    if (pfa_min >= pfa_max) return -1;

    double snr_linear = pow(10.0, snr_db / 10.0);
    double a = sqrt(snr_linear);

    for (size_t i = 0; i < num_pts; i++) {
        /* Log-spaced PFA values */
        double log_min = log10(pfa_min);
        double log_max = log10(pfa_max);
        double log_pfa = log_min + (log_max - log_min) * (double)i / (double)(num_pts - 1);
        pfa[i] = pow(10.0, log_pfa);

        double z_pfa = Q_inverse(pfa[i]);
        pd[i] = Q_function(z_pfa - a);
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Range accuracy — Cramer-Rao Lower Bound
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_range_crlb(double snr_linear, double b_rms) {
    if (snr_linear <= 0.0 || b_rms <= 0.0) return 0.0;
    /* sigma_R^2 >= c^2 / (8 * SNR * B_rms^2) */
    double sigma_R_sq = (LIDAR_C * LIDAR_C) / (8.0 * snr_linear * b_rms * b_rms);
    return sqrt(sigma_R_sq);
}