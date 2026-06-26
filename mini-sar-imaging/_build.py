
import os
BASE = os.path.dirname(os.path.abspath(__file__))

def write_file(relpath, content):
    full = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  {relpath}: {len(content.splitlines())} lines")

# Files will be added here via subsequent appends
print("_build.py base ready")


# ==== Content strings ====

GEO_H = """
#ifndef SAR_GEOMETRY_H
#define SAR_GEOMETRY_H
#include "sar_core.h"

typedef enum { SAR_COORD_SLANT_RANGE=0, SAR_COORD_GROUND_RANGE=1, SAR_COORD_GEOGRAPHIC=2, SAR_COORD_MAP=3 } sar_coord_type_t;
typedef struct { double x_m, y_m, z_m; } sar_ecef_t;
typedef struct { double lat_deg, lon_deg, alt_m; } sar_geodetic_t;

double sar_range_hyperbolic(double R0, double v, double eta, double eta_0);
double sar_range_parabolic(double R0, double v, double eta, double eta_0);
double sar_range_quartic(double R0, double v, double eta, double eta_0);
double sar_doppler_centroid(double lambda, double v, double theta_sq);
double sar_doppler_rate(double lambda, double v, double R0, double theta_sq);
double sar_synthetic_aperture_time(double lambda, double R0, double L_a, double v, double theta_sq);

typedef struct {
    double R0, v, lambda, f_Dc, f_R;
    double rcm_linear_coeff, rcm_quadratic_coeff;
    double max_rcm_m, range_cell_size_m;
    size_t max_rcm_cells;
} sar_rcm_t;

void sar_rcm_init(sar_rcm_t *rcm, double R0, double v, double lambda, double theta_sq, double range_cell_size);
double sar_rcm_at_eta(const sar_rcm_t *rcm, double eta);
double sar_slant_to_ground(double R_slant, double H);
double sar_ground_to_slant(double G, double H);
double sar_incidence_angle(double H, double R_slant);
void sar_ecef_to_geodetic(const sar_ecef_t *ecef, sar_geodetic_t *geo);
void sar_geodetic_to_ecef(const sar_geodetic_t *geo, sar_ecef_t *ecef);

typedef struct {
    double antenna_length_m, wavelength_m;
    double beamwidth_rad, beamwidth_az_km;
} sar_antenna_t;

void sar_antenna_init(sar_antenna_t *ant, double L_a, double lambda, double R0);
double sar_antenna_gain(const sar_antenna_t *ant, double theta);
double sar_effective_velocity(double v_plat, double H, double Re);
int sar_resolution_check(const sar_params_t *sp, double tol);
int sar_prf_check(const sar_params_t *sp);

#endif
"""

print("Content strings defined")


# Write sar_geometry.h
write_file("include/sar_geometry.h", GEO_H)

# Write sar_algorithm.h
ALG_H = """
#ifndef SAR_ALGORITHM_H
#define SAR_ALGORITHM_H
#include "sar_core.h"
#include "sar_geometry.h"

typedef struct {
    int use_secondary_range_compression;
    int rcmc_interpolation_method;
    double range_oversample_factor;
    size_t rcmc_kernel_length;
} sar_rda_params_t;

int sar_rda_process(const sar_raw_data_t *raw, const sar_rda_params_t *rda_p, sar_image_t *image);
void sar_rda_range_compression(const sar_raw_data_t *raw, double **compressed_I, double **compressed_Q);
void sar_rda_rcmc(double **rda_domain_I, double **rda_domain_Q, size_t naz, size_t nrng, const sar_params_t *sp, const sar_rda_params_t *rda_p);
void sar_rda_azimuth_compression(double **rda_domain_I, double **rda_domain_Q, size_t naz, size_t nrng, const sar_params_t *sp);

typedef struct { double scaling_factor; int use_extended_csa; } sar_csa_params_t;
int sar_csa_process(const sar_raw_data_t *raw, const sar_csa_params_t *csa_p, sar_image_t *image);
void sar_csa_chirp_scaling(double **data_I, double **data_Q, size_t naz, size_t nrng, const sar_params_t *sp, const sar_csa_params_t *csa_p);

int sar_omegak_process(const sar_raw_data_t *raw, double R_ref, sar_image_t *image);
void sar_omegak_ref_multiply(double **data_I, double **data_Q, size_t naz, size_t nrng, const sar_params_t *sp, double R_ref);
void sar_omegak_stolt(double **data_I, double **data_Q, size_t naz, size_t nrng, const sar_params_t *sp);

typedef struct {
    double grid_x_min, grid_x_max, grid_y_min, grid_y_max;
    size_t grid_nx, grid_ny;
    double pixel_size_x, pixel_size_y;
} sar_bp_params_t;
int sar_bp_process(const sar_raw_data_t *raw, const sar_bp_params_t *bp_p, const double *traj_x, const double *traj_y, const double *traj_z, sar_image_t *image);

int sar_specan_process(const sar_raw_data_t *raw, size_t burst_len, sar_image_t *image);

int sar_pga_autofocus(sar_image_t *image, int niter);
int sar_mapdrift_autofocus(sar_image_t *image);

typedef struct {
    double range_resolution_3dB_m, azimuth_resolution_3dB_m;
    double range_pslr_db, azimuth_pslr_db;
    double range_islr_db, azimuth_islr_db;
    double peak_range_idx, peak_azimuth_idx;
    double peak_magnitude, peak_phase_rad;
} sar_impulse_response_t;
void sar_analyze_impulse_response(const sar_image_t *image, size_t r0, size_t c0, size_t win, sar_impulse_response_t *ir);
#endif
"""
write_file("include/sar_algorithm.h", ALG_H)

print("Headers 2-3 written")

# Write sar_interferometry.h + sar_advanced.h
IFR_H = """
#ifndef SAR_INTERFEROMETRY_H
#define SAR_INTERFEROMETRY_H
#include "sar_core.h"
#include "sar_geometry.h"

/* --- L1: Interferometric Definitions --- */
typedef struct {
    sar_image_t *master;
    sar_image_t *slave;
    double baseline_m;
    double baseline_parallel_m;
    double baseline_perp_m;
    double baseline_angle_rad;
    double height_ambiguity_m;
    double critical_baseline_m;
} sar_insar_pair_t;

/* --- L1: Coherence --- */
typedef struct {
    size_t nrows, ncols;
    double **coherence;
    double mean_coherence;
    double **interferometric_phase;
    double **phase_stddev;
} sar_coherence_map_t;

sar_insar_pair_t *sar_insar_pair_alloc(void);
void sar_insar_pair_free(sar_insar_pair_t *pair);
void sar_insar_set_baseline(sar_insar_pair_t *pair, double B, double alpha, double H, double R0, double lambda);

sar_coherence_map_t *sar_coherence_alloc(size_t nrows, size_t ncols);
void sar_coherence_free(sar_coherence_map_t *cmap);

/* L5: Coherence estimation via spatial averaging */
void sar_coherence_estimate(const sar_image_t *master, const sar_image_t *slave, sar_coherence_map_t *cmap, size_t win_rows, size_t win_cols);

/* L5: Interferometric phase computation */
void sar_interferogram_compute(const sar_image_t *master, const sar_image_t *slave, double **ifgram, size_t nrows, size_t ncols);

/* L5: Flat-earth phase removal */
void sar_flat_earth_removal(double **ifgram, size_t nrows, size_t ncols, double range_spacing, double lambda, double H, double look_angle);

/* L5: Phase unwrapping - Goldstein branch cut */
int sar_phase_unwrap_goldstein(double **wrapped_phase, double **unwrapped_phase, size_t nrows, size_t ncols);

/* L5: Phase unwrapping - quality guided */
int sar_phase_unwrap_quality(double **wrapped_phase, double **unwrapped_phase, size_t nrows, size_t ncols);

/* L5: Phase to height conversion */
void sar_phase_to_height(const double **unwrapped_phase, size_t nrows, size_t ncols, double B_perp, double R0, double lambda, double H, double theta, double **dem);

/* L7: DInSAR - differential interferometry for displacement */
int sar_dinsar_displacement(const sar_image_t *master, const sar_image_t *slave, double **displacement_m, size_t nrows, size_t ncols, double lambda, const double **topo_phase);

/* L7: Phase filtering - Goldstein filter */
void sar_goldstein_filter(double **ifgram, size_t nrows, size_t ncols, double alpha, int patch_size);

#endif
"""
write_file("include/sar_interferometry.h", IFR_H)

ADV_H = """
#ifndef SAR_ADVANCED_H
#define SAR_ADVANCED_H
#include "sar_core.h"

/* --- L8: Compressive Sensing SAR --- */
typedef struct {
    double undersampling_ratio;
    double regularization_lambda;
    int max_iterations;
    double convergence_tol;
    int algorithm_type;
} sar_cs_params_t;

int sar_cs_reconstruct(const double **measurements, size_t n_meas, size_t n_pixels, const sar_cs_params_t *csp, double *reconstructed);

/* --- L8: MIMO-SAR --- */
typedef struct {
    size_t n_tx;
    size_t n_rx;
    double *tx_positions_x;
    double *rx_positions_x;
    double *tx_waveforms;
    size_t waveform_length;
} sar_mimo_config_t;

int sar_mimo_image(const sar_mimo_config_t *cfg, double **image, size_t nrows, size_t ncols);

/* --- L8: Bistatic SAR --- */
typedef struct {
    double tx_x, tx_y, tx_z;
    double rx_x, rx_y, rx_z;
    double tx_vel_x, tx_vel_y, tx_vel_z;
    double rx_vel_x, rx_vel_y, rx_vel_z;
} sar_bistatic_config_t;

double sar_bistatic_range(const sar_bistatic_config_t *bc, double eta, double x_tgt, double y_tgt, double z_tgt);
int sar_bistatic_processor(const sar_raw_data_t *raw, const sar_bistatic_config_t *bc, sar_image_t *image);

/* --- L8: Polarimetric SAR --- */
typedef struct {
    double **HH_I, **HH_Q;
    double **HV_I, **HV_Q;
    double **VH_I, **VH_Q;
    double **VV_I, **VV_Q;
    size_t nrows, ncols;
} sar_polarimetric_t;

typedef struct {
    double T11, T22, T33;
    double T12_real, T12_imag;
    double T13_real, T13_imag;
    double T23_real, T23_imag;
} sar_coherency_matrix_t;

/* Freeman-Durden 3-component decomposition */
void sar_freeman_durden(const sar_coherency_matrix_t *T, size_t n, double *Ps, double *Pd, double *Pv);

/* H-Alpha decomposition */
void sar_h_alpha_decomp(const sar_coherency_matrix_t *T, size_t n, double *entropy, double *alpha, double *anisotropy);

#endif
"""
write_file("include/sar_advanced.h", ADV_H)

print("Headers 4-5 written")


# === Write C source files ===

print("Writing C source files...")

# ============================================================
# src/sar_core.c - Core SAR functionality implementation
# ============================================================
CORE_C = """
/**
 * @file    sar_core.c
 * @brief   SAR Core Implementations -- L1 Definitions + L3 Math + L4 Laws
 *
 * @details Implements: SAR parameter initialization, chirp generation and
 *          matched filtering, raw echo data simulation (point target),
 *          SAR image management (multi-look, calibration, magnitude/phase),
 *          2D FFT operations, wavenumber domain mapping, and Stolt interpolation.
 *
 * Knowledge Mapping:
 *   L1 - sar_params_init, sar_chirp_alloc, sar_chirp_autocorrelation,
 *        sar_pulse_compression, sar_raw_data_point_target,
 *        sar_image_magnitude/phase, sar_multilook, sar_calibrate_sigma0
 *   L3 - sar_fft2d, sar_fft_range, sar_fft_azimuth, sar_wavenumber_init,
 *        sar_stolt_interpolation
 *   L4 - Range resolution rho_r = c/(2*B_r) verified in params_init
 *        Nyquist PRF constraint in azimuth bandwidth computation
 *
 * Reference:
 *   - Cumming & Wong, "Digital Processing of SAR Data" (2005)
 *   - Richards, Scheer & Holm, "Principles of Modern Radar" (2010)
 *   - Cooley & Tukey, "An Algorithm for the Machine Calculation of
 *     Complex Fourier Series" (1965), Math. Comp.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sar_core.h"

/* ==================================================================
 * L1: SAR Parameter Initialization
 * ================================================================== */

/**
 * Initialize all SAR system parameters from fundamental values.
 *
 * Derived quantities:
 *   lambda          = c / f0                 (wavelength)
 *   chirp_rate      = B_r / tau_p            (Hz/s)
 *   range_resolution = c / (2 * B_r)          (m)
 *   azimuth_resolution = L_a / 2              (m, focused stripmap)
 *   swath_width     = R_far*cos(theta_i_far) - R_near*cos(theta_i_near)
 *
 * The azimuth resolution rho_a = L_a/2 is the fundamental SAR result:
 * finer resolution is achieved with a smaller antenna, counter-intuitively,
 * because a smaller antenna creates a wider beam, hence a longer synthetic
 * aperture and greater Doppler bandwidth.
 */
void sar_params_init(sar_params_t *p,
                     double freq, double bw, double pw, double prf,
                     double ant, double vel, double alt,
                     double squint, double look,
                     double near, double far)
{
    if (!p) return;
    memset(p, 0, sizeof(*p));

    p->carrier_freq_hz     = freq;
    p->wavelength_m        = SAR_C / freq;
    p->bandwidth_hz        = bw;
    p->pulse_width_s       = pw;
    p->prf_hz              = prf;
    p->chirp_rate          = (pw > 0.0) ? bw / pw : 0.0;
    p->sample_rate_hz      = bw * 1.2;
    p->range_resolution_m  = (bw > 0.0) ? SAR_C / (2.0 * bw) : 0.0;
    p->azimuth_resolution_m = ant / 2.0;
    p->antenna_length_m    = ant;
    p->platform_velocity_ms = vel;
    p->squint_angle_rad    = squint;
    p->look_angle_rad      = look;
    p->incidence_angle_rad = asin(alt / near);
    p->platform_altitude_m = alt;
    p->near_range_m        = near;
    p->far_range_m         = far;

    /* Swath width in ground range */
    double theta_i_near = asin(alt / near);
    double theta_i_far  = asin(alt / far);
    double G_near = near * cos(theta_i_near);
    double G_far  = far  * cos(theta_i_far);
    p->swath_width_m = G_far - G_near;
    if (p->swath_width_m < 0.0) p->swath_width_m = 0.0;
}

/* ==================================================================
 * L1: Chirp (Linear FM) Generation
 * ================================================================== */

/**
 * Allocate and generate a complex baseband chirp:
 *
 *   s[n] = exp(j * pi * K * (n*dt)^2)  for n = 0..N-1
 *   I[n] = cos(pi * K * (n*dt)^2)
 *   Q[n] = sin(pi * K * (n*dt)^2)
 *
 * where K = B/tau is the chirp rate and dt = 1/fs.
 *
 * Time-bandwidth product TBP = B * tau = K * tau^2 characterizes
 * the pulse compression gain: compressed pulse width ~ 1/B.
 *
 * @param N    Number of samples
 * @param tau  Pulse width [s]
 * @param B    Bandwidth [Hz]
 * @param fs   Sampling rate [Hz]
 * @return     Allocated chirp structure, or NULL on failure
 */
sar_chirp_t *sar_chirp_alloc(size_t N, double tau, double B, double fs)
{
    if (N == 0 || tau <= 0.0 || B <= 0.0 || fs <= 0.0) return NULL;

    sar_chirp_t *c = (sar_chirp_t *)malloc(sizeof(sar_chirp_t));
    if (!c) return NULL;

    c->I = (double *)malloc(N * sizeof(double));
    c->Q = (double *)malloc(N * sizeof(double));
    if (!c->I || !c->Q) {
        free(c->I); free(c->Q); free(c);
        return NULL;
    }

    c->num_samples  = N;
    c->pulse_width_s = tau;
    c->bandwidth_hz  = B;
    c->chirp_rate    = B / tau;
    c->sample_rate_hz = fs;
    c->dt            = 1.0 / fs;
    c->t_min         = -tau / 2.0;
    c->t_max         = +tau / 2.0;
    c->time_bandwidth_product = B * tau;

    /* Generate I/Q samples */
    double K  = c->chirp_rate;
    double dt = c->dt;
    for (size_t n = 0; n < N; n++) {
        double t = (double)((int64_t)n - (int64_t)(N/2)) * dt;
        double phase = M_PI * K * t * t;
        c->I[n] = cos(phase);
        c->Q[n] = sin(phase);
    }
    return c;
}

void sar_chirp_free(sar_chirp_t *c)
{
    if (!c) return;
    free(c->I);
    free(c->Q);
    free(c);
}

/**
 * Compute chirp autocorrelation (matched filter response):
 *
 *   R[k] = sum_{n} s[n] * conj(s[n-k])
 *
 * For a chirp, R[k] approximates a sinc function:
 *   |R(tau)| ~ |sinc(B * tau)|  for |tau| <= tau_p
 *
 * Peak-to-Sidelobe Ratio (PSLR) for rectangular window:
 *   PSLR ~ -13.26 dB (first sidelobe of sinc)
 *
 * The -3dB width of the mainlobe is approximately 0.886/B,
 * giving the compressed range resolution rho_r = c/(2B).
 *
 * @param c      Chirp signal
 * @param R_out  Output autocorrelation magnitude [2*N-1]
 */
void sar_chirp_autocorrelation(const sar_chirp_t *c, double *R_out)
{
    if (!c || !R_out) return;
    size_t N = c->num_samples;
    size_t M = 2 * N - 1;

    for (size_t k = 0; k < M; k++) {
        double sum_I = 0.0, sum_Q = 0.0;
        int64_t lag = (int64_t)k - (int64_t)(N - 1);

        for (size_t n = 0; n < N; n++) {
            int64_t nk = (int64_t)n - lag;
            if (nk >= 0 && (size_t)nk < N) {
                /* conj multiply: (I1+jQ1)*(I2-jQ2) = (I1*I2+Q1*Q2) + j(I2*Q1-I1*Q2) */
                sum_I += c->I[n] * c->I[nk] + c->Q[n] * c->Q[nk];
                sum_Q += c->I[nk] * c->Q[n] - c->I[n] * c->Q[nk];
            }
        }
        R_out[k] = sqrt(sum_I * sum_I + sum_Q * sum_Q) / (double)N;
    }
}
"""

write_file("src/sar_core.c", CORE_C)
print("sar_core.c part 1 written")
