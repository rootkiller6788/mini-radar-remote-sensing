/**
 * @file    ir_calibration.c
 * @brief   IR FPA Calibration - NUC, Bad Pixel Detection/Replacement
 *
 * Non-uniformity correction (NUC) is critical for IR focal plane arrays.
 * Each pixel has unique gain and offset; NUC compensates pixel-to-pixel
 * variation for radiometrically accurate imagery.
 *
 * References:
 *   Perry & Dereniak (1993) Optical Engineering, 32(8)
 *   Scribner et al. (1990) SPIE Vol. 1308
 */

#include "ir_calibration.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ========================================================================
 * Non-Uniformity Correction
 * ======================================================================== */

int ir_nuc_one_point_correction(ir_nuc_state_t *nuc,
                                 const double *raw_frame,
                                 double *corrected_out) {
    if (!nuc || !raw_frame || !corrected_out) return -1;
    if (!nuc->offset || nuc->rows <= 0 || nuc->cols <= 0) return -1;

    int N = nuc->rows * nuc->cols;
    for (int i = 0; i < N; i++) {
        corrected_out[i] = raw_frame[i] - nuc->offset[i];
    }
    return 0;
}

int ir_cal_two_point_compute(const double *raw_low, const double *raw_high,
                              double T_low_K, double T_high_K,
                              int rows, int cols,
                              ir_cal_pixel_coeff_t *coeffs_out) {
    if (!raw_low || !raw_high || !coeffs_out) return -1;
    if (rows <= 0 || cols <= 0) return -1;
    if (T_low_K <= 0.0 || T_high_K <= T_low_K) return -1;

    int N = rows * cols;
    double dT = T_high_K - T_low_K;

    for (int i = 0; i < N; i++) {
        double diff = raw_high[i] - raw_low[i];
        if (fabs(diff) < 1e-10) {
            /* Degenerate pixel: mark as bad */
            coeffs_out[i].gain = 0.0;
            coeffs_out[i].offset = raw_low[i];
            coeffs_out[i].is_bad = 1;
        } else {
            coeffs_out[i].gain = dT / diff;
            coeffs_out[i].offset = raw_low[i];
            coeffs_out[i].is_bad = 0;
        }
    }
    return 0;
}

int ir_cal_two_point_apply(const ir_cal_pixel_coeff_t *coeffs,
                            const double *raw_frame,
                            int rows, int cols,
                            double *corrected_out) {
    if (!coeffs || !raw_frame || !corrected_out) return -1;
    if (rows <= 0 || cols <= 0) return -1;

    int N = rows * cols;
    for (int i = 0; i < N; i++) {
        if (coeffs[i].is_bad) {
            corrected_out[i] = NAN;
        } else {
            corrected_out[i] = coeffs[i].gain * (raw_frame[i] - coeffs[i].offset);
        }
    }
    return 0;
}

int ir_cal_multi_point_apply(const ir_multi_point_cal_t *cal,
                              const double *raw_value,
                              int rows, int cols,
                              double *corrected_radiance_out) {
    if (!cal || !raw_value || !corrected_radiance_out) return -1;
    if (cal->n_points < 2 || rows <= 0 || cols <= 0) return -1;

    int N = rows * cols;
    for (int i = 0; i < N; i++) {
        double rv = raw_value[i];

        /* Find bracketing calibration points */
        int idx = -1;
        for (int k = 1; k < cal->n_points; k++) {
            if (cal->raw_values[k][i] >= rv) {
                idx = k - 1;
                break;
            }
        }
        if (idx < 0) idx = cal->n_points - 2;

        double r_lo = cal->raw_values[idx][i];
        double r_hi = cal->raw_values[idx + 1][i];
        double T_lo = cal->T_points_K[idx];
        double T_hi = cal->T_points_K[idx + 1];

        if (fabs(r_hi - r_lo) < 1e-10) {
            corrected_radiance_out[i] = NAN;
        } else {
            double frac = (rv - r_lo) / (r_hi - r_lo);
            corrected_radiance_out[i] = T_lo + frac * (T_hi - T_lo);
        }
    }
    return 0;
}

/* ========================================================================
 * Scene-Based NUC (Constant-Statistics)
 *
 * Iteratively adjusts gain/offset to equalize mean and variance across
 * all pixels, using the assumption that scene statistics are stationary.
 * ======================================================================== */

int ir_nuc_scene_based_update(ir_nuc_state_t *nuc,
                               const double *raw_frame,
                               int rows, int cols,
                               double learning_rate) {
    if (!nuc || !raw_frame) return -1;
    if (rows <= 0 || cols <= 0) return -1;
    if (learning_rate <= 0.0 || learning_rate > 1.0) return -1;

    int N = rows * cols;

    /* Compute frame mean */
    double mean = 0.0;
    for (int i = 0; i < N; i++) mean += raw_frame[i];
    mean /= N;

    /* Update offset map to drive each pixel toward the global mean */
    if (nuc->offset) {
        for (int i = 0; i < N; i++) {
            double error = raw_frame[i] - mean;
            nuc->offset[i] += learning_rate * error;
        }
    }

    return 0;
}

/* ========================================================================
 * Bad Pixel Detection and Replacement
 * ======================================================================== */

int ir_bad_pixel_detect(const double *raw_frame, int rows, int cols,
                         double n_sigma, ir_bad_pixel_map_t *map_out) {
    if (!raw_frame || !map_out) return -1;
    if (rows <= 0 || cols <= 0 || n_sigma <= 0.0) return -1;

    int N = rows * cols;

    /* Compute mean and stddev */
    double mean = 0.0, var = 0.0;
    for (int i = 0; i < N; i++) mean += raw_frame[i];
    mean /= N;
    for (int i = 0; i < N; i++) {
        double d = raw_frame[i] - mean;
        var += d * d;
    }
    double stddev = sqrt(var / N);

    /* Flag pixels outside n_sigma */
    map_out->rows = rows;
    map_out->cols = cols;
    map_out->n_bad = 0;
    double threshold = n_sigma * stddev;

    for (int i = 0; i < N; i++) {
        if (fabs(raw_frame[i] - mean) > threshold) {
            map_out->map[i] = 1;
            map_out->n_bad++;
        } else {
            map_out->map[i] = 0;
        }
    }
    return 0;
}

int ir_bad_pixel_replace(double *frame, int rows, int cols,
                          const ir_bad_pixel_map_t *map) {
    if (!frame || !map) return -1;
    if (rows <= 0 || cols <= 0) return -1;

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            int idx = r * cols + c;
            if (!map->map[idx]) continue;

            /* Median of nearest good neighbors (3x3 kernel) */
            double neighbors[8];
            int n_neigh = 0;
            for (int dr = -1; dr <= 1; dr++) {
                for (int dc = -1; dc <= 1; dc++) {
                    if (dr == 0 && dc == 0) continue;
                    int nr = r + dr, nc = c + dc;
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                    int nidx = nr * cols + nc;
                    if (map->map[nidx]) continue;
                    neighbors[n_neigh++] = frame[nidx];
                }
            }
            if (n_neigh > 0) {
                /* Simple bubble sort for median */
                for (int j = 0; j < n_neigh - 1; j++)
                    for (int k = j + 1; k < n_neigh; k++)
                        if (neighbors[j] > neighbors[k]) {
                            double tmp = neighbors[j];
                            neighbors[j] = neighbors[k];
                            neighbors[k] = tmp;
                        }
                frame[idx] = neighbors[n_neigh / 2];
            }
        }
    }
    return 0;
}

int ir_dead_pixel_detect(const ir_cal_pixel_coeff_t *coeffs,
                          int rows, int cols,
                          double gain_threshold,
                          ir_bad_pixel_map_t *map_out) {
    if (!coeffs || !map_out) return -1;
    if (rows <= 0 || cols <= 0 || gain_threshold <= 0.0) return -1;

    int N = rows * cols;
    map_out->rows = rows;
    map_out->cols = cols;
    map_out->n_bad = 0;

    for (int i = 0; i < N; i++) {
        if (coeffs[i].is_bad || fabs(coeffs[i].gain) < gain_threshold) {
            map_out->map[i] = 1;
            map_out->n_bad++;
        } else {
            map_out->map[i] = 0;
        }
    }
    return 0;
}

int ir_drifting_pixel_detect(const double *temporal_stddev,
                              int rows, int cols,
                              double noise_threshold,
                              ir_bad_pixel_map_t *map_out) {
    if (!temporal_stddev || !map_out) return -1;
    if (rows <= 0 || cols <= 0 || noise_threshold <= 0.0) return -1;

    int N = rows * cols;
    map_out->rows = rows;
    map_out->cols = cols;
    map_out->n_bad = 0;

    for (int i = 0; i < N; i++) {
        if (temporal_stddev[i] > noise_threshold) {
            map_out->map[i] = 1;
            map_out->n_bad++;
        } else {
            map_out->map[i] = 0;
        }
    }
    return 0;
}

/* ========================================================================
 * Temperature Calibration
 * ======================================================================== */

double ir_cal_raw_to_temperature(double raw_value,
                                  double gain_coeff,
                                  double offset_coeff,
                                  double lambda_um) {
    if (gain_coeff <= 0.0) return -1.0;

    /* Corrected signal = gain * (raw - offset) */
    double corrected = gain_coeff * (raw_value - offset_coeff);
    if (corrected <= 0.0) return 0.0;

    /* Invert Planck: T = c2 / (lambda * ln(1 + c1/(lambda^5 * L))) */
    double c1 = IR_FIRST_RADIATION_C1;
    double c2 = IR_SECOND_RADIATION_C2;
    double lm = lambda_um * 1e-6;
    double l5 = lm * lm * lm * lm * lm;

    return c2 / (lm * log(1.0 + c1 / (l5 * corrected)));
}

int ir_cal_shutter_correction(ir_nuc_state_t *nuc,
                               const double *shutter_frame,
                               int rows, int cols) {
    if (!nuc || !shutter_frame) return -1;
    if (rows <= 0 || cols <= 0) return -1;
    if (!nuc->offset) return -1;

    int N = rows * cols;
    for (int i = 0; i < N; i++) {
        nuc->offset[i] = shutter_frame[i];
    }
    return 0;
}

double ir_nuc_residual_nuc(const double *corrected, int rows, int cols) {
    if (!corrected || rows <= 0 || cols <= 0) return -1.0;

    int N = rows * cols;
    double mean = 0.0;
    for (int i = 0; i < N; i++) mean += corrected[i];
    mean /= N;
    if (fabs(mean) < 1e-10) return 0.0;

    double rms = 0.0;
    for (int i = 0; i < N; i++) {
        double d = corrected[i] - mean;
        rms += d * d;
    }
    return 100.0 * sqrt(rms / N) / fabs(mean);  /* percent */
}
