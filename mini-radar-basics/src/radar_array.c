/**
 * @file    radar_array.c
 * @brief   Phased array basics for radar: beamforming, steering vectors, MIMO
 *
 * Knowledge: L1 ULA geometry, L5 beamforming/steering, L8 MIMO virtual array
 */
#include "radar_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "radar_array: alloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* L1: Uniform Linear Array geometry */
typedef struct {
    size_t n_elements;
    double element_spacing_m;
    double wavelength_m;
    double *element_positions;
} ula_array_t;

/* L5: Initialize ULA */
int ula_init(ula_array_t *arr, size_t n_el, double spacing, double lambda)
{
    if (!arr || n_el < 1 || spacing <= 0.0 || lambda <= 0.0) return -1;
    arr->n_elements = n_el;
    arr->element_spacing_m = spacing;
    arr->wavelength_m = lambda;
    arr->element_positions = safe_malloc(n_el * sizeof(double));
    for (size_t i = 0; i < n_el; i++)
        arr->element_positions[i] = (double)i * spacing;
    return 0;
}

/* L5: Array factor - complex sum of element contributions */
double complex ula_array_factor(const ula_array_t *arr,
                                 double theta_rad)
{
    if (!arr || arr->n_elements == 0) return CMPLX(0.0, 0.0);
    double d = arr->element_spacing_m;
    double lambda = arr->wavelength_m;
    double k = 2.0 * M_PI / lambda;
    double complex af = CMPLX(0.0, 0.0);
    for (size_t i = 0; i < arr->n_elements; i++) {
        double phase = k * (double)i * d * sin(theta_rad);
        af += CMPLX(cos(phase), sin(phase));
    }
    return af;
}

/* L5: Steering vector for beamforming */
void ula_steering_vector(const ula_array_t *arr,
                         double theta_rad,
                         double complex *sv, size_t sv_len)
{
    if (!arr || !sv || sv_len < arr->n_elements) return;
    double d = arr->element_spacing_m;
    double lambda = arr->wavelength_m;
    double k = 2.0 * M_PI / lambda;
    for (size_t i = 0; i < arr->n_elements; i++) {
        double phase = -k * (double)i * d * sin(theta_rad);
        sv[i] = CMPLX(cos(phase), sin(phase));
    }
}

/* L5: Normalized beam pattern (power) */
double ula_beam_pattern(const ula_array_t *arr, double theta_rad)
{
    if (!arr || arr->n_elements == 0) return 0.0;
    double d = arr->element_spacing_m;
    double lambda = arr->wavelength_m;
    double k = 2.0 * M_PI * d * sin(theta_rad) / lambda;
    if (fabs(k) < 1e-12) return 1.0;
    double N = (double)arr->n_elements;
    double val = sin(N * k / 2.0) / (N * sin(k / 2.0));
    return val * val;
}

/* L5: Half-power beamwidth for ULA (broadside) */
double ula_beamwidth_broadside(const ula_array_t *arr)
{
    if (!arr || arr->n_elements < 2) return M_PI;
    double d = arr->element_spacing_m;
    double L = (double)(arr->n_elements - 1) * d;
    double lambda = arr->wavelength_m;
    return 0.886 * lambda / L;
}

/* L8: MIMO virtual array */
size_t mimo_virtual_elements(size_t n_tx, size_t n_rx)
{
    return n_tx * n_rx;
}

void mimo_virtual_positions(const double *tx_pos, size_t n_tx,
                            const double *rx_pos, size_t n_rx,
                            double *virt_pos)
{
    if (!tx_pos || !rx_pos || !virt_pos) return;
    for (size_t t = 0; t < n_tx; t++)
        for (size_t r = 0; r < n_rx; r++)
            virt_pos[t * n_rx + r] = tx_pos[t] + rx_pos[r];
}

/* L5: Digital Beamforming (DBF) */
void ula_dbf(const double complex *elements, size_t n_el,
             const double complex *weights,
             double complex *output)
{
    if (!elements || !weights || !output) return;
    *output = CMPLX(0.0, 0.0);
    for (size_t i = 0; i < n_el; i++)
        *output += elements[i] * conj(weights[i]);
}

/* L5: Phase shift for electronic beam steering */
double ula_phase_shift(size_t element_idx, double theta_rad,
                       double spacing_m, double lambda_m)
{
    if (lambda_m <= 0.0) return 0.0;
    double k = 2.0 * M_PI / lambda_m;
    return -k * (double)element_idx * spacing_m * sin(theta_rad);
}

/* L7: Grating lobe angle condition */
double ula_grating_lobe_angle(double lambda_m, double spacing_m)
{
    if (spacing_m <= 0.0 || lambda_m >= spacing_m) return -1.0;
    return asin(lambda_m / spacing_m);
}

/* L1: Array gain from coherent combining */
double ula_array_gain_db(size_t n_elements)
{
    if (n_elements == 0) return 0.0;
    return 10.0 * log10((double)n_elements);
}

/* L5: Free ULA resources */
void ula_free(ula_array_t *arr)
{
    if (!arr) return;
    free(arr->element_positions);
    arr->element_positions = NULL;
    arr->n_elements = 0;
}

/* L7: Beam pattern in dB */
double ula_beam_pattern_db(const ula_array_t *arr, double theta_rad)
{
    double p = ula_beam_pattern(arr, theta_rad);
    if (p <= 0.0) return -100.0;
    return 10.0 * log10(p);
}

/* L7: Compute array factor for rectangular planar array */
typedef struct {
    size_t n_rows, n_cols;
    double dx_m, dy_m;
    double wavelength_m;
} planar_array_t;

int planar_array_init(planar_array_t *arr, size_t rows, size_t cols,
                      double dx, double dy, double lambda)
{
    if (!arr || rows < 1 || cols < 1 || dx <= 0.0 || dy <= 0.0 || lambda <= 0.0)
        return -1;
    arr->n_rows = rows;
    arr->n_cols = cols;
    arr->dx_m = dx;
    arr->dy_m = dy;
    arr->wavelength_m = lambda;
    return 0;
}

double complex planar_array_factor(const planar_array_t *arr,
                                    double theta_rad, double phi_rad)
{
    if (!arr) return CMPLX(0.0, 0.0);
    double k = 2.0 * M_PI / arr->wavelength_m;
    double complex af = CMPLX(0.0, 0.0);
    double u = sin(theta_rad) * cos(phi_rad);
    double v = sin(theta_rad) * sin(phi_rad);
    for (size_t r = 0; r < arr->n_rows; r++) {
        for (size_t c = 0; c < arr->n_cols; c++) {
            double phase = k * ((double)r * arr->dy_m * v +
                                (double)c * arr->dx_m * u);
            af += CMPLX(cos(phase), sin(phase));
        }
    }
    return af;
}

/* L8: Adaptive beamforming - MVDR (Capon) weight computation */
void mvdr_weights(const double complex *cov_matrix, size_t n_el,
                  double theta_rad, double spacing_m, double lambda_m,
                  double complex *weights)
{
    if (!cov_matrix || !weights || n_el == 0 || lambda_m <= 0.0) return;
    /* Simple: use steering vector directly (no interference nulling) */
    double k = 2.0 * M_PI / lambda_m;
    for (size_t i = 0; i < n_el; i++) {
        double phase = -k * (double)i * spacing_m * sin(theta_rad);
        weights[i] = CMPLX(cos(phase), sin(phase));
    }
}

/* L8: MUSIC pseudo-spectrum for DOA estimation */
double music_spectrum_sample(const double complex *steering_vec,
                              const double complex *noise_subspace,
                              size_t n_el, size_t n_signals)
{
    if (!steering_vec || !noise_subspace || n_el == 0 || n_signals >= n_el)
        return 0.0;
    size_t n_noise = n_el - n_signals;
    double num = 0.0;
    for (size_t i = 0; i < n_el; i++)
        num += cabs(steering_vec[i]) * cabs(steering_vec[i]);
    double den = 0.0;
    for (size_t j = 0; j < n_noise; j++) {
        double complex proj = CMPLX(0.0, 0.0);
        for (size_t i = 0; i < n_el; i++)
            proj += conj(noise_subspace[i * n_noise + j]) * steering_vec[i];
        den += cabs(proj) * cabs(proj);
    }
    if (den <= 0.0) return 1e6;
    return num / den;
}

/* L7: Compute array manifold (steering vectors over angle scan) */
void ula_array_manifold(const ula_array_t *arr, double theta_start_rad,
                         double theta_end_rad, size_t n_angles,
                         double complex *manifold, size_t manifold_stride)
{
    if (!arr || !manifold || n_angles < 2) return;
    double dtheta = (theta_end_rad - theta_start_rad) / (double)(n_angles - 1);
    size_t n_el = arr->n_elements;
    for (size_t a = 0; a < n_angles; a++) {
        double theta = theta_start_rad + (double)a * dtheta;
        double complex *sv = manifold + a * manifold_stride;
        ula_steering_vector(arr, theta, sv, n_el);
    }
}

/* L7: Element pattern multiplication (pattern multiplication principle) */
double array_total_pattern(double element_pattern, double array_factor_mag)
{
    return element_pattern * array_factor_mag;
}

/* L8: Computed directivity via numerical integration */
double ula_directivity_numerical(const ula_array_t *arr)
{
    if (!arr || arr->n_elements == 0) return 0.0;
    size_t n_theta = 360;
    double sum = 0.0;
    for (size_t i = 0; i < n_theta; i++) {
        double theta = -M_PI / 2.0 + M_PI * (double)i / (double)(n_theta - 1);
        double pat = ula_beam_pattern(arr, theta);
        sum += pat * cos(theta);
    }
    double avg = sum / (double)n_theta;
    if (avg <= 0.0) return 0.0;
    return 2.0 / avg;
}
