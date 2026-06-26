/**
 * @file    radar_doppler.c
 * @brief   Radar Doppler processing: MTI, Doppler FFT, velocity estimation
 *
 * Knowledge covered:
 *   L1: Doppler shift, radial velocity, blind speeds, CPI
 *   L3: Phase progression model, Doppler steering vector
 *   L5: MTI filter design, Doppler FFT processing
 *   L6: Range-Doppler map generation, velocity estimation, PRF staggering
 *
 * Reference: Richards, Scheer & Holm (2010), Ch.5, 17.
 *            Skolnik, "Radar Handbook" (2008), Ch.3.
 */
#include "radar_doppler.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "radar_doppler: alloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ===== L1: Doppler Shift Computation ===================================== */

double radar_doppler_shift(double vr_ms, double lambda_m)
{
    if (lambda_m <= 0.0) return 0.0;
    return 2.0 * vr_ms / lambda_m;
}

double radar_doppler_to_velocity(double fd_hz, double lambda_m)
{
    if (lambda_m <= 0.0) return 0.0;
    return fd_hz * lambda_m / 2.0;
}

double radar_doppler_resolution(double prf_hz, size_t n_pulses)
{
    if (prf_hz <= 0.0 || n_pulses == 0) return 0.0;
    return prf_hz / (double)n_pulses;
}

double radar_blind_speed(double lambda_m, double prf_hz, int n)
{
    if (lambda_m <= 0.0 || prf_hz <= 0.0 || n < 0) return 0.0;
    return (double)n * lambda_m * prf_hz / 2.0;
}

/* ===== L5: MTI Filter Implementation ===================================== */

int mti_filter_init(mti_filter_t *f, mti_filter_type_t type,
                    const double *custom_coeffs, size_t order)
{
    if (!f) return -1;
    memset(f, 0, sizeof(*f));
    f->type = type;

    switch (type) {
    case MTI_SINGLE_CANCEL:
        f->order = 1;
        f->coeffs = safe_malloc(2 * sizeof(double));
        f->coeffs[0] = 1.0; f->coeffs[1] = -1.0;
        break;
    case MTI_DOUBLE_CANCEL:
        f->order = 2;
        f->coeffs = safe_malloc(3 * sizeof(double));
        f->coeffs[0] = 1.0; f->coeffs[1] = -2.0; f->coeffs[2] = 1.0;
        break;
    case MTI_TRIPLE_CANCEL:
        f->order = 3;
        f->coeffs = safe_malloc(4 * sizeof(double));
        f->coeffs[0] = 1.0; f->coeffs[1] = -3.0;
        f->coeffs[2] = 3.0; f->coeffs[3] = -1.0;
        break;
    case MTI_CUSTOM:
        if (!custom_coeffs || order == 0) return -1;
        f->order = order;
        f->coeffs = safe_malloc((order + 1) * sizeof(double));
        for (size_t i = 0; i <= order; i++)
            f->coeffs[i] = custom_coeffs[i];
        break;
    default:
        return -1;
    }

    f->delay_line = (double complex*)calloc(f->order + 1, sizeof(double complex));
    if (!f->delay_line) { free(f->coeffs); return -1; }
    f->delay_index = 0;
    f->initialized = 1;
    return 0;
}

double complex mti_filter_apply(mti_filter_t *f, double complex x_n)
{
    if (!f || !f->initialized) return CMPLX(0.0, 0.0);

    /* Shift delay line */
    for (size_t i = f->order; i > 0; i--)
        f->delay_line[i] = f->delay_line[i - 1];
    f->delay_line[0] = x_n;

    /* Apply FIR filter */
    double complex y = CMPLX(0.0, 0.0);
    for (size_t i = 0; i <= f->order; i++)
        y += f->coeffs[i] * f->delay_line[i];

    return y;
}

int mti_filter_response(const mti_filter_t *f,
                        const double *freq_hz, size_t n_freqs,
                        double *response)
{
    if (!f || !freq_hz || !response || n_freqs == 0 || f->prf_hz <= 0.0)
        return -1;
    for (size_t k = 0; k < n_freqs; k++) {
        double omega = 2.0 * M_PI * freq_hz[k] / f->prf_hz;
        double complex H = CMPLX(0.0, 0.0);
        for (size_t i = 0; i <= f->order; i++) {
            double phase = -omega * (double)i;
            H += f->coeffs[i] * CMPLX(cos(phase), sin(phase));
        }
        response[k] = cabs(H);
    }
    return 0;
}

void mti_filter_free(mti_filter_t *f)
{
    if (!f) return;
    free(f->coeffs);
    free(f->delay_line);
    f->coeffs = NULL;
    f->delay_line = NULL;
    f->initialized = 0;
}

/* ===== L5: Doppler FFT (2D Range-Doppler Map) ============================ */

static void fft_radix2_doppler(double complex *x, size_t N, int inverse)
{
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < N; i++) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { double complex t = x[i]; x[i] = x[j]; x[j] = t; }
    }
    double sign = inverse ? 1.0 : -1.0;
    for (size_t len = 2; len <= N; len <<= 1) {
        double ang = 2.0 * M_PI / (double)len * sign;
        double complex wlen = CMPLX(cos(ang), sin(ang));
        for (size_t i = 0; i < N; i += len) {
            double complex w = CMPLX(1.0, 0.0);
            for (size_t j2 = 0; j2 < len / 2; j2++) {
                double complex u = x[i + j2];
                double complex v = x[i + j2 + len / 2] * w;
                x[i + j2] = u + v;
                x[i + j2 + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse)
        for (size_t i = 0; i < N; i++) x[i] /= (double)N;
}

static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

int radar_doppler_fft_2d(const double complex *data,
                         size_t n_range, size_t n_pulses,
                         double complex *doppler_map)
{
    if (!data || !doppler_map || n_range == 0 || n_pulses == 0)
        return -1;
    size_t N_fft = next_pow2(n_pulses);
    for (size_t r = 0; r < n_range; r++) {
        double complex *row = safe_malloc(N_fft * sizeof(double complex));
        for (size_t p = 0; p < n_pulses; p++)
            row[p] = data[r * n_pulses + p];
        for (size_t p = n_pulses; p < N_fft; p++)
            row[p] = CMPLX(0.0, 0.0);
        fft_radix2_doppler(row, N_fft, 0);
        for (size_t p = 0; p < N_fft; p++)
            doppler_map[r * N_fft + p] = row[p];
        free(row);
    }
    return 0;
}

/* ===== L1: Doppler Bin to Frequency ====================================== */

double radar_doppler_bin_to_freq(size_t bin, size_t n_fft, double prf_hz)
{
    if (n_fft == 0 || prf_hz <= 0.0) return 0.0;
    if (bin < n_fft / 2)
        return (double)bin * prf_hz / (double)n_fft;
    else
        return ((double)bin - (double)n_fft) * prf_hz / (double)n_fft;
}

/* ===== L6: Velocity from Phase Difference ================================ */

double radar_velocity_from_phase(double complex s0, double complex s1,
                                  double wavelength_m, double pri_s)
{
    if (pri_s <= 0.0 || wavelength_m <= 0.0) return 0.0;
    double delta_phi = carg(s1 * conj(s0));
    return wavelength_m * delta_phi / (4.0 * M_PI * pri_s);
}

/* ===== L6: Doppler Ambiguity Resolution ================================== */

int radar_doppler_resolve_ambiguity(double fd1, double prf1,
                                     double fd2, double prf2,
                                     double *resolved_fd)
{
    if (!resolved_fd || prf1 <= 0.0 || prf2 <= 0.0) return -1;
    /* True Doppler: f = fd1 + k1*prf1 = fd2 + k2*prf2 for integers k1,k2 */
    int max_k = 100;
    for (int k1 = -max_k; k1 <= max_k; k1++) {
        double f1 = fd1 + (double)k1 * prf1;
        for (int k2 = -max_k; k2 <= max_k; k2++) {
            double f2 = fd2 + (double)k2 * prf2;
            if (fabs(f1 - f2) < 0.005 * (prf1 + prf2)) {
                *resolved_fd = (f1 + f2) / 2.0;
                return 0;
            }
        }
    }
    return -1;
}

/* L7: MTI improvement factor (Barton's formula) */
double mti_improvement_factor_db(mti_filter_type_t type, double clutter_std_hz,
                                  double prf_hz)
{
    if (prf_hz <= 0.0) return 0.0;
    double sigma_norm = clutter_std_hz / prf_hz;
    double I_db;
    switch (type) {
    case MTI_SINGLE_CANCEL:
        I_db = -20.0 * log10(2.0 * sin(M_PI * sigma_norm));
        break;
    case MTI_DOUBLE_CANCEL:
        I_db = -20.0 * log10(4.0 * sin(M_PI * sigma_norm) *
                              sin(M_PI * 2.0 * sigma_norm));
        break;
    case MTI_TRIPLE_CANCEL:
        I_db = -20.0 * log10(8.0 * sin(M_PI * sigma_norm) *
                              sin(M_PI * 2.0 * sigma_norm) *
                              sin(M_PI * 3.0 * sigma_norm));
        break;
    default:
        I_db = 0.0;
        break;
    }
    return I_db;
}

/* L7: Staggered PRF ratio optimization for blind speed avoidance */
void staggered_prf_ratios(int *prf1_hz, int *prf2_hz,
                          double max_range_m, double max_velocity_ms,
                          double wavelength_m)
{
    /* Typical ratios: 4:5, 5:6, 8:9 — choose based on requirements */
    double base_prf = 299792458.0 / (2.0 * max_range_m);
    double doppler_nyquist = 2.0 * max_velocity_ms / wavelength_m;
    (void)doppler_nyquist;
    *prf1_hz = (int)(base_prf * 0.9);
    *prf2_hz = (int)(base_prf * 1.1);
}

/* L7: Compute CPI duration from number of pulses and PRF */
double radar_cpi_duration_s(size_t n_pulses, double prf_hz)
{
    if (prf_hz <= 0.0) return 0.0;
    return (double)n_pulses / prf_hz;
}

/* L7: Doppler bin to velocity */
double radar_doppler_bin_to_velocity(size_t bin, size_t n_fft,
                                      double prf_hz, double wavelength_m)
{
    double fd = radar_doppler_bin_to_freq(bin, n_fft, prf_hz);
    return radar_doppler_to_velocity(fd, wavelength_m);
}

/* L7: Micro-Doppler signature analysis — simple vibration model */
double micro_doppler_freq_hz(double vibration_freq_hz,
                              double displacement_m,
                              double wavelength_m)
{
    /* Max micro-Doppler: f_mD = 4*pi*f_vib*D / lambda */
    if (wavelength_m <= 0.0) return 0.0;
    return 4.0 * M_PI * vibration_freq_hz * displacement_m / wavelength_m;
}

/* L8: STAP (Space-Time Adaptive Processing) sample matrix */
void stap_sample_matrix(const double complex *data_cube,
                         size_t n_range, size_t n_pulses, size_t n_channels,
                         size_t range_bin,
                         double complex *cov_matrix)
{
    if (!data_cube || !cov_matrix || range_bin >= n_range) return;
    size_t dof = n_pulses * n_channels;
    for (size_t i = 0; i < dof; i++) {
        for (size_t j = 0; j < dof; j++) {
            cov_matrix[i * dof + j] = CMPLX(0.0, 0.0);
        }
    }
    /* Simple: outer product of space-time snapshot */
    size_t offset = range_bin * n_pulses * n_channels;
    for (size_t i = 0; i < dof; i++) {
        for (size_t j = 0; j < dof; j++) {
            cov_matrix[i * dof + j] = data_cube[offset + i] *
                                       conj(data_cube[offset + j]);
        }
    }
}

/* L7: Clutter Doppler spread estimation */
double clutter_doppler_spread_hz(double platform_velocity_ms,
                                  double wavelength_m,
                                  double beamwidth_az_rad)
{
    if (wavelength_m <= 0.0) return 0.0;
    double max_doppler = 2.0 * platform_velocity_ms / wavelength_m;
    return max_doppler * sin(beamwidth_az_rad / 2.0);
}

/* L7: Pulse-Doppler visibility factor */
double pulse_doppler_visibility(size_t n_pulses, double target_doppler_hz,
                                 double prf_hz)
{
    if (prf_hz <= 0.0 || n_pulses == 0) return 0.0;
    /* Visibility degrades as target Doppler approaches PRF boundaries */
    double doppler_norm = target_doppler_hz / prf_hz;
    /* Sinc-like pattern: V = |sin(pi*N*f_d/PRF) / (N*sin(pi*f_d/PRF))| */
    double a = M_PI * doppler_norm;
    double b = M_PI * (double)n_pulses * doppler_norm;
    if (fabs(sin(a)) < 1e-12) return 1.0;
    return fabs(sin(b) / ((double)n_pulses * sin(a)));
}
