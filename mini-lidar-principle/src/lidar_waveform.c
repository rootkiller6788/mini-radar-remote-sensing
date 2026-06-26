/**
 * @file    lidar_waveform.c
 * @brief   LiDAR full-waveform processing implementation
 *
 * Knowledge covered:
 *   L1: Waveform digitization, pulse width measurement
 *   L3: Gaussian function, multi-Gaussian superposition
 *   L5: Peak detection (CFD, leading-edge, derivative), noise estimation,
 *       Gaussian decomposition via Levenberg-Marquardt,
 *       pulse characterization (FWHM, energy, skewness)
 *
 * Reference:
 *   - Wagner et al., *ISPRS JPRS* 60(2), pp.100-112, 2006.
 *   - Marquardt, D.W., *J. SIAM* 11(2), pp.431-441, 1963.
 */

#include "lidar_waveform.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ─── Internal helpers ──────────────────────────────────────────────────── */

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "lidar_waveform: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SQRT_2PI  2.5066282746310005024
#define SQRT_2LN2 1.1774100225154746910  /* sqrt(2*ln(2)) */

/* ─── Box-Muller Gaussian random number ──────────────────────────────────── */

static double rand_gauss(void) {
    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;
    if (u1 < 1e-10) u1 = 1e-10;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Waveform creation and management
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_waveform_init(lidar_waveform_t *wf, size_t num_samples,
                         double dt, double t_offset, double pulse_width) {
    if (!wf || num_samples == 0 || dt <= 0.0) return -1;
    wf->samples = (double*)safe_malloc(num_samples * sizeof(double));
    wf->num_samples = num_samples;
    wf->dt = dt;
    wf->t_offset = t_offset;
    wf->pulse_width = pulse_width;
    wf->noise_floor = 0.0;
    wf->noise_std = 0.0;
    wf->max_amplitude = 0.0;
    wf->max_index = 0;
    memset(wf->samples, 0, num_samples * sizeof(double));
    return 0;
}

void lidar_waveform_free(lidar_waveform_t *wf) {
    if (!wf) return;
    free(wf->samples);
    memset(wf, 0, sizeof(*wf));
}

int lidar_waveform_synthesize(lidar_waveform_t *wf,
                               const lidar_gaussian_component_t *components,
                               size_t num_comp,
                               double noise_std) {
    if (!wf || !components || num_comp == 0 || num_comp > LIDAR_WAVEFORM_MAX_RETURNS)
        return -1;

    for (size_t i = 0; i < wf->num_samples; i++) {
        double t = wf->t_offset + (double)i * wf->dt;
        double val = 0.0;
        for (size_t j = 0; j < num_comp; j++) {
            val += lidar_gaussian_eval(&components[j], t);
        }
        if (noise_std > 0.0) {
            val += noise_std * rand_gauss();
        }
        wf->samples[i] = val;
        if (val > wf->max_amplitude) {
            wf->max_amplitude = val;
            wf->max_index = i;
        }
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Peak detection algorithms
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_cfd_timing(const lidar_waveform_t *wf, double fraction) {
    if (!wf || !wf->samples || wf->num_samples == 0) return -1.0;
    if (fraction <= 0.0 || fraction > 1.0) return -1.0;

    double threshold = wf->noise_floor + fraction * (wf->max_amplitude - wf->noise_floor);
    size_t wi = wf->max_index;

    /* Search backward from peak for threshold crossing */
    while (wi > 0 && wf->samples[wi] > threshold) wi--;
    if (wi >= wf->num_samples - 1) return -1.0;

    /* Linear interpolation between samples wi and wi+1 */
    double y0 = wf->samples[wi];
    double y1 = wf->samples[wi + 1];
    if (fabs(y1 - y0) < 1e-15) {
        return wf->t_offset + (double)wi * wf->dt;
    }
    double frac = (threshold - y0) / (y1 - y0);
    return wf->t_offset + ((double)wi + frac) * wf->dt;
}

double lidar_leading_edge_timing(const lidar_waveform_t *wf, double k_sigma) {
    if (!wf || !wf->samples || wf->num_samples == 0) return -1.0;

    double threshold = wf->noise_floor + k_sigma * wf->noise_std;
    for (size_t i = 0; i < wf->num_samples - 1; i++) {
        if (wf->samples[i] < threshold && wf->samples[i + 1] >= threshold) {
            double y0 = wf->samples[i];
            double y1 = wf->samples[i + 1];
            double frac;
            if (fabs(y1 - y0) < 1e-15) {
                frac = 0.0;
            } else {
                frac = (threshold - y0) / (y1 - y0);
            }
            return wf->t_offset + ((double)i + frac) * wf->dt;
        }
    }
    return -1.0;
}

int lidar_detect_peaks_derivative(const lidar_waveform_t *wf,
                                   double *peak_times,
                                   size_t max_peaks,
                                   double min_sep,
                                   size_t smooth_win) {
    if (!wf || !peak_times || max_peaks == 0 || wf->num_samples < 3) return 0;

    size_t n = wf->num_samples;
    double *diff = (double*)malloc(n * sizeof(double));
    if (!diff) return 0;

    /* First derivative via central differences (smoothed) */
    for (size_t i = 0; i < n; i++) {
        size_t i0 = (i > smooth_win) ? i - smooth_win : 0;
        size_t i1 = (i + smooth_win < n) ? i + smooth_win : n - 1;
        diff[i] = (wf->samples[i1] - wf->samples[i0]) / (double)(i1 - i0);
    }

    size_t num_found = 0;
    for (size_t i = 1; i < n - 1 && num_found < max_peaks; i++) {
        /* Positive to negative zero-crossing of derivative = local maximum */
        if (diff[i - 1] > 0.0 && diff[i] <= 0.0) {
            /* Check that signal is above noise floor */
            if (wf->samples[i] > wf->noise_floor + 3.0 * wf->noise_std) {
                double peak_t = wf->t_offset + (double)i * wf->dt;

                /* Minimum separation check */
                int too_close = 0;
                for (size_t j = 0; j < num_found; j++) {
                    if (fabs(peak_t - peak_times[j]) < min_sep) {
                        too_close = 1;
                        /* Keep the larger peak */
                        if (wf->samples[i] > wf->samples[(size_t)((peak_times[j] - wf->t_offset) / wf->dt)]) {
                            peak_times[j] = peak_t;
                        }
                        break;
                    }
                }
                if (!too_close) {
                    peak_times[num_found++] = peak_t;
                }
            }
        }
    }

    free(diff);
    return (int)num_found;
}

int lidar_waveform_noise_estimate(lidar_waveform_t *wf, double tail_ratio) {
    if (!wf || tail_ratio <= 0.0 || tail_ratio >= 1.0) return -1;

    size_t tail_start = (size_t)((double)wf->num_samples * (1.0 - tail_ratio));
    if (tail_start >= wf->num_samples) tail_start = wf->num_samples - 1;
    size_t tail_len = wf->num_samples - tail_start;
    if (tail_len < 4) return -1;

    /* Mean */
    double sum = 0.0;
    for (size_t i = tail_start; i < wf->num_samples; i++) {
        sum += wf->samples[i];
    }
    wf->noise_floor = sum / (double)tail_len;

    /* Standard deviation */
    double sum_sq = 0.0;
    for (size_t i = tail_start; i < wf->num_samples; i++) {
        double d = wf->samples[i] - wf->noise_floor;
        sum_sq += d * d;
    }
    wf->noise_std = sqrt(sum_sq / (double)tail_len);

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Gaussian decomposition
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_gaussian_decompose(const lidar_waveform_t *wf,
                              lidar_waveform_decomp_t *result,
                              int max_iter,
                              double tol) {
    if (!wf || !result || wf->num_samples == 0) return -1;

    memset(result, 0, sizeof(*result));

    /* Step 1: Peak detection */
    double peak_times[LIDAR_WAVEFORM_MAX_RETURNS];
    int n_peaks = lidar_detect_peaks_derivative(wf, peak_times,
        LIDAR_WAVEFORM_MAX_RETURNS, wf->pulse_width, 2);
    if (n_peaks == 0) return 0;

    /* Step 2: Initial parameter estimation */
    size_t n_comp = (size_t)n_peaks;
    if (n_comp > LIDAR_WAVEFORM_MAX_RETURNS) n_comp = LIDAR_WAVEFORM_MAX_RETURNS;

    /* Parameter vector: {A_1,mu_1,sigma_1,...,A_n,mu_n,sigma_n,offset}
       3*n_comp + 1 parameters */
    size_t n_params = 3 * n_comp + 1;
    double *params = (double*)calloc(n_params, sizeof(double));
    if (!params) return -1;

    for (size_t j = 0; j < n_comp; j++) {
        size_t idx = (size_t)((peak_times[j] - wf->t_offset) / wf->dt);
        if (idx >= wf->num_samples) idx = wf->num_samples - 1;
        params[3*j]     = wf->samples[idx];           /* amplitude */
        params[3*j + 1] = peak_times[j];               /* center */
        params[3*j + 2] = wf->pulse_width / 2.35482;   /* sigma from FWHM */
    }
    params[n_params - 1] = wf->noise_floor;            /* offset */

    /* Step 3: Levenberg-Marquardt fitting */
    int n_iter = lidar_levenberg_marquardt_gaussian(wf, params, n_comp, max_iter, 0.001, tol);

    /* Step 4: Fill result */
    result->num_components = n_comp;
    result->converged = (n_iter < max_iter) ? 1 : 0;
    result->num_iterations = abs(n_iter);

    for (size_t j = 0; j < n_comp && j < LIDAR_WAVEFORM_MAX_RETURNS; j++) {
        result->components[j].amplitude = params[3*j];
        result->components[j].center    = params[3*j + 1];
        result->components[j].sigma     = params[3*j + 2];
        result->components[j].fwhm      = params[3*j + 2] * 2.35482;
        result->components[j].range     = params[3*j + 1] * LIDAR_C / 2000.0;
        result->components[j].energy    = params[3*j] * params[3*j + 2] * sqrt(2.0 * M_PI);
    }

    /* Residual RMSE */
    double rms = 0.0;
    double signal_var = 0.0;
    double signal_mean = 0.0;
    for (size_t i = 0; i < wf->num_samples; i++) {
        double t = wf->t_offset + (double)i * wf->dt;
        double model = lidar_multigaussian_eval(result->components, n_comp,
                                                 params[n_params - 1], t);
        double resid = wf->samples[i] - model;
        rms += resid * resid;
        signal_mean += wf->samples[i];
    }
    result->residual_rms = sqrt(rms / (double)wf->num_samples);
    signal_mean /= (double)wf->num_samples;
    for (size_t i = 0; i < wf->num_samples; i++) {
        double d = wf->samples[i] - signal_mean;
        signal_var += d * d;
    }
    if (signal_var > 1e-20) {
        result->r_squared = 1.0 - rms / signal_var;
    } else {
        result->r_squared = 0.0;
    }

    free(params);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Levenberg-Marquardt for Gaussian model
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_levenberg_marquardt_gaussian(const lidar_waveform_t *wf,
                                        double *params,
                                        size_t num_comp,
                                        int max_iter,
                                        double lambda,
                                        double tol) {
    if (!wf || !params || num_comp == 0) return -1;

    size_t n_params = 3 * num_comp + 1;
    size_t n_data   = wf->num_samples;
    double lam = lambda;

    /* Workspace */
    double *residual = (double*)malloc(n_data * sizeof(double));
    double *J        = (double*)malloc(n_data * n_params * sizeof(double));
    double *JTJ      = (double*)malloc(n_params * n_params * sizeof(double));
    double *JTr      = (double*)malloc(n_params * sizeof(double));
    double *delta    = (double*)malloc(n_params * sizeof(double));
    double *params_new = (double*)malloc(n_params * sizeof(double));
    if (!residual || !J || !JTJ || !JTr || !delta || !params_new) {
        free(residual); free(J); free(JTJ); free(JTr); free(delta); free(params_new);
        return -1;
    }

    double prev_cost = 1e100;
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* Compute residuals and Jacobian */
        double cost = 0.0;
        for (size_t i = 0; i < n_data; i++) {
            double t = wf->t_offset + (double)i * wf->dt;
            double model = params[n_params - 1];
            for (size_t j = 0; j < num_comp; j++) {
                double A   = params[3*j];
                double mu  = params[3*j + 1];
                double sig = params[3*j + 2];
                if (sig < 1e-30) sig = 1e-30;
                double dt  = (t - mu) / sig;
                double g   = exp(-0.5 * dt * dt);
                model     += A * g;
                /* Jacobian entries */
                J[i * n_params + 3*j]     = g;                           /* d/dA */
                J[i * n_params + 3*j + 1] = A * g * dt / sig;            /* d/dmu */
                J[i * n_params + 3*j + 2] = A * g * dt * dt / sig;      /* d/dsigma */
            }
            J[i * n_params + n_params - 1] = 1.0;                        /* d/doffset */
            residual[i] = wf->samples[i] - model;
            cost += residual[i] * residual[i];
        }

        if (fabs(cost - prev_cost) < tol * (prev_cost + 1.0)) {
            break;
        }
        prev_cost = cost;

        /* Form normal equations: (J^T J + lambda * diag(J^T J)) * delta = J^T r */
        memset(JTJ, 0, n_params * n_params * sizeof(double));
        memset(JTr, 0, n_params * sizeof(double));

        for (size_t i = 0; i < n_data; i++) {
            for (size_t a = 0; a < n_params; a++) {
                double Jia = J[i * n_params + a];
                JTr[a] += Jia * residual[i];
                for (size_t b = a; b < n_params; b++) {
                    JTJ[a * n_params + b] += Jia * J[i * n_params + b];
                }
            }
        }
        /* Symmetric fill */
        for (size_t a = 0; a < n_params; a++) {
            for (size_t b = 0; b < a; b++) {
                JTJ[a * n_params + b] = JTJ[b * n_params + a];
            }
        }

        /* Add damping */
        double *A_aug = JTJ; /* reuse JTJ as augmented matrix */
        for (size_t a = 0; a < n_params; a++) {
            A_aug[a * n_params + a] += lam * A_aug[a * n_params + a];
            if (A_aug[a * n_params + a] < 1e-50) A_aug[a * n_params + a] = 1.0;
        }

        /* Solve via Gaussian elimination (Cholesky would be better but this is portable) */
        /* Copy into workspace for solve */
        double *A_solve = (double*)malloc(n_params * n_params * sizeof(double));
        memcpy(A_solve, A_aug, n_params * n_params * sizeof(double));
        memcpy(delta, JTr, n_params * sizeof(double));

        /* Gaussian elimination with partial pivoting */
        for (size_t k = 0; k < n_params; k++) {
            /* Find pivot */
            size_t pivot = k;
            double max_val = fabs(A_solve[k * n_params + k]);
            for (size_t r = k + 1; r < n_params; r++) {
                if (fabs(A_solve[r * n_params + k]) > max_val) {
                    max_val = fabs(A_solve[r * n_params + k]);
                    pivot = r;
                }
            }
            /* Swap rows */
            if (pivot != k) {
                for (size_t c = 0; c < n_params; c++) {
                    double tmp = A_solve[k * n_params + c];
                    A_solve[k * n_params + c] = A_solve[pivot * n_params + c];
                    A_solve[pivot * n_params + c] = tmp;
                }
                double tmp = delta[k];
                delta[k] = delta[pivot];
                delta[pivot] = tmp;
            }
            /* Eliminate */
            double pivot_val = A_solve[k * n_params + k];
            if (fabs(pivot_val) < 1e-50) continue;
            for (size_t r = k + 1; r < n_params; r++) {
                double factor = A_solve[r * n_params + k] / pivot_val;
                for (size_t c = k; c < n_params; c++) {
                    A_solve[r * n_params + c] -= factor * A_solve[k * n_params + c];
                }
                delta[r] -= factor * delta[k];
            }
        }
        /* Back substitution */
        for (int k_int = (int)n_params - 1; k_int >= 0; k_int--) {
            size_t k = (size_t)k_int;
            double sum = delta[k];
            for (size_t c = k + 1; c < n_params; c++) {
                sum -= A_solve[k * n_params + c] * delta[c];
            }
            if (fabs(A_solve[k * n_params + k]) > 1e-50) {
                delta[k] = sum / A_solve[k * n_params + k];
            } else {
                delta[k] = 0.0;
            }
        }
        free(A_solve);

        /* Trial step */
        memcpy(params_new, params, n_params * sizeof(double));
        for (size_t a = 0; a < n_params; a++) {
            params_new[a] += delta[a];
        }
        /* Constrain sigmas > 0 */
        for (size_t j = 0; j < num_comp; j++) {
            if (params_new[3*j + 2] < wf->dt * 0.5) params_new[3*j + 2] = wf->dt * 0.5;
        }

        /* Compute new cost */
        double cost_new = 0.0;
        for (size_t i = 0; i < n_data; i++) {
            double t = wf->t_offset + (double)i * wf->dt;
            double model = params_new[n_params - 1];
            for (size_t j = 0; j < num_comp; j++) {
                double A   = params_new[3*j];
                double mu  = params_new[3*j + 1];
                double sig = params_new[3*j + 2];
                double dt2  = (t - mu) / sig;
                model += A * exp(-0.5 * dt2 * dt2);
            }
            double r = wf->samples[i] - model;
            cost_new += r * r;
        }

        if (cost_new < cost) {
            /* Accept step, reduce lambda */
            memcpy(params, params_new, n_params * sizeof(double));
            lam *= 0.1;
        } else {
            /* Reject step, increase lambda */
            lam *= 10.0;
        }
        if (lam > 1e10) lam = 1e10;
        if (lam < 1e-10) lam = 1e-10;
    }

    free(residual); free(J); free(JTJ); free(JTr); free(delta); free(params_new);
    return iter;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Pulse characterization
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_pulse_fwhm(const lidar_waveform_t *wf, size_t center) {
    if (!wf || !wf->samples || center >= wf->num_samples) return -1.0;

    double peak = wf->samples[center];
    double half_peak = (peak - wf->noise_floor) * 0.5 + wf->noise_floor;

    /* Search left */
    size_t left = center;
    while (left > 0 && wf->samples[left] > half_peak) left--;
    double t_left;
    if (left < wf->num_samples - 1 && wf->samples[left + 1] > wf->samples[left]) {
        double frac = (half_peak - wf->samples[left]) / (wf->samples[left + 1] - wf->samples[left]);
        t_left = wf->t_offset + ((double)left + frac) * wf->dt;
    } else {
        t_left = wf->t_offset + (double)left * wf->dt;
    }

    /* Search right */
    size_t right = center;
    while (right < wf->num_samples - 1 && wf->samples[right] > half_peak) right++;
    double t_right;
    if (right > 0 && wf->samples[right - 1] > wf->samples[right]) {
        double frac = (wf->samples[right - 1] - half_peak) / (wf->samples[right - 1] - wf->samples[right]);
        t_right = wf->t_offset + ((double)right - frac) * wf->dt;
    } else {
        t_right = wf->t_offset + (double)right * wf->dt;
    }

    return t_right - t_left;
}

double lidar_pulse_energy(const lidar_waveform_t *wf,
                           size_t center, size_t half_win) {
    if (!wf || !wf->samples) return 0.0;

    size_t start = (center > half_win) ? center - half_win : 0;
    size_t end   = (center + half_win < wf->num_samples) ? center + half_win : wf->num_samples - 1;
    if (end <= start) return 0.0;

    double sum = 0.0;
    for (size_t i = start; i <= end; i++) {
        double val = wf->samples[i] - wf->noise_floor;
        if (val > 0.0) sum += val;
    }
    return sum * wf->dt;
}

double lidar_pulse_skewness(const lidar_waveform_t *wf,
                             size_t center, size_t half_win) {
    if (!wf || !wf->samples || half_win == 0) return 0.0;

    size_t start = (center > half_win) ? center - half_win : 0;
    size_t end   = (center + half_win < wf->num_samples) ? center + half_win : wf->num_samples - 1;
    if (end <= start + 2) return 0.0;

    /* Compute mean */
    double sum = 0.0, sum_w = 0.0;
    for (size_t i = start; i <= end; i++) {
        double val = wf->samples[i] - wf->noise_floor;
        if (val > 0.0) {
            sum += (double)i * val;
            sum_w += val;
        }
    }
    if (sum_w < 1e-15) return 0.0;
    double mean = sum / sum_w;

    /* Compute std and skewness */
    double m2 = 0.0, m3 = 0.0;
    sum_w = 0.0;
    for (size_t i = start; i <= end; i++) {
        double val = wf->samples[i] - wf->noise_floor;
        if (val > 0.0) {
            double d = (double)i - mean;
            m2 += val * d * d;
            m3 += val * d * d * d;
            sum_w += val;
        }
    }
    if (sum_w < 1e-15 || m2 < 1e-30) return 0.0;
    double variance = m2 / sum_w;
    double std_dev = sqrt(variance);
    if (std_dev < 1e-30) return 0.0;
    return (m3 / sum_w) / (std_dev * std_dev * std_dev);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Gaussian model evaluation
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_gaussian_eval(const lidar_gaussian_component_t *comp, double t) {
    if (!comp || comp->sigma <= 0.0) return 0.0;
    double dt = (t - comp->center) / comp->sigma;
    return comp->amplitude * exp(-0.5 * dt * dt);
}

double lidar_multigaussian_eval(const lidar_gaussian_component_t *components,
                                 size_t num_comp, double offset, double t) {
    double val = offset;
    for (size_t j = 0; j < num_comp; j++) {
        val += lidar_gaussian_eval(&components[j], t);
    }
    return val;
}