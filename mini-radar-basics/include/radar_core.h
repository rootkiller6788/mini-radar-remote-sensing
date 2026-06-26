/**
 * @file    radar_core.h
 * @brief   Core Radar Definitions �� L1 Definitions + L3 Mathematical Structures
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Radar range equation (monostatic, bistatic)
 *     - Radar cross section (RCS), Swerling fluctuation models
 *     - Pulse width (tau), PRF, PRI, duty cycle
 *     - Range resolution, unambiguous range, unambiguous velocity
 *     - Antenna gain G, effective aperture A_e, beamwidth
 *     - Peak power P_t, average power P_avg
 *     - SNR, noise figure F, noise temperature T_sys
 *     - Radar wavelength lambda, operating frequency
 *   L2 - Core Concepts:
 *     - Time-of-flight ranging principle
 *     - Pulse radar vs CW radar vs FMCW
 *     - Power-aperture product
 *     - System noise model and loss budget
 *   L3 - Mathematical Structures:
 *     - Radar range equation as Friis-like power budget
 *     - Complex baseband representation
 *     - Point target scattering model
 *     - Free-space path loss
 *   L4 - Fundamental Laws:
 *     - Radar range equation: P_r = P_t*G^2*lambda^2*sigma / ((4*pi)^3*R^4*L)
 *     - Maximum detection range (fourth-root dependence)
 *     - Unambiguous range/velocity constraints
 *
 * Reference: Richards, Scheer & Holm, "Principles of Modern Radar" (2010), Ch.1-2.
 *            Skolnik, "Radar Handbook" (2008), Ch.1.
 */
#ifndef RADAR_CORE_H
#define RADAR_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef CMPLX
#define CMPLX(r, i) ((double complex)((double)(r) + I * (double)(i)))
#endif

/* Physical constants */
#define RADAR_C               299792458.0
#define RADAR_K_BOLTZMANN     1.380649e-23
#define RADAR_T0              290.0

/* ������ L1: Radar Operating Mode ������������������������������������������������������������������������������ */

typedef enum {
    RADAR_MODE_PULSE,
    RADAR_MODE_CW,
    RADAR_MODE_FMCW,
    RADAR_MODE_PULSE_DOPPLER,
    RADAR_MODE_STEPPED_FREQ,
    RADAR_MODE_NOISE
} radar_mode_t;

typedef enum {
    RADAR_POL_HH,
    RADAR_POL_VV,
    RADAR_POL_HV,
    RADAR_POL_VH,
    RADAR_POL_CIRCULAR_LL,
    RADAR_POL_CIRCULAR_LR
} radar_polarization_t;

/* ������ L1: Radar System Parameters ������������������������������������������������������������������������ */

typedef struct {
    double center_freq_hz;
    double wavelength_m;
    double peak_power_w;
    double antenna_gain_db;
    double antenna_gain_linear;
    double pulse_width_s;
    double prf_hz;
    double pri_s;
    double bandwidth_hz;
    double noise_figure_db;
    double noise_figure_linear;
    double system_loss_db;
    double system_loss_linear;
    double antenna_temp_k;
    radar_mode_t mode;
    radar_polarization_t pol;
} radar_params_t;

int radar_params_init(radar_params_t *params,
                      double freq_hz, double peak_power_w, double ant_gain_db,
                      double tau_s, double prf_hz, double bw_hz,
                      double nf_db, double loss_db, double ant_temp_k,
                      radar_mode_t mode);

/* ������ L1: Radar Cross Section (RCS) �������������������������������������������������������������������� */

typedef enum {
    RCS_CONSTANT,
    RCS_SWERLING_I,
    RCS_SWERLING_II,
    RCS_SWERLING_III,
    RCS_SWERLING_IV
} rcs_model_t;

typedef struct {
    double mean_rcs_dbsm;
    double mean_rcs_m2;
    rcs_model_t fluctuation;
    double decorrelation_time_s;
    size_t n_samples;
} radar_rcs_t;

int radar_rcs_init(radar_rcs_t *rcs, double rcs_dbsm, rcs_model_t model,
                   double decorr_time_s);
double radar_rcs_sample(radar_rcs_t *rcs);

/* ������ L1: Range Resolution and Unambiguous Range ���������������������������������������� */

double radar_range_resolution(double pulse_width_or_bw, int use_bandwidth);
double radar_unambiguous_range(double prf_hz);
double radar_unambiguous_velocity(double wavelength_m, double prf_hz);

/* ������ L3/L4: Radar Range Equation ������������������������������������������������������������������������ */

double radar_received_power(const radar_params_t *params,
                            double rcs_m2, double range_m);
double radar_snr(const radar_params_t *params, double rcs_m2, double range_m);
double radar_snr_db(const radar_params_t *params, double rcs_m2, double range_m);
double radar_max_range(const radar_params_t *params, double rcs_m2,
                       double snr_min);
double radar_bistatic_power(double pt_w, double gt_linear, double gr_linear,
                            double lambda_m, double rcs_m2,
                            double rt_m, double rr_m, double loss_linear);

/* ������ L2: System Noise ���������������������������������������������������������������������������������������������� */

double radar_system_noise_temp(double ant_temp_k, double noise_fig_lin);
double radar_noise_power(double tsys_k, double bw_hz);

/* ������ L3: Propagation ������������������������������������������������������������������������������������������������ */

double radar_path_loss_2way(double range_m, double lambda_m, double gain_linear);

/* ������ L1: Antenna Model �������������������������������������������������������������������������������������������� */

double radar_antenna_gain_from_aperture(double a_e_m2, double lambda_m);
double radar_beamwidth_circular(double lambda_m, double diameter_m);
double radar_directivity_kraus(double beamwidth_az_rad, double beamwidth_el_rad);

/* ������ L2: Duty Cycle and Average Power �������������������������������������������������������������� */

double radar_duty_cycle(double tau_s, double prf_hz);
double radar_average_power(double peak_power_w, double duty_cycle);

/* ������ L5: Pulse Integration Gain �������������������������������������������������������������������������� */

double radar_coherent_integration_gain(size_t n_pulses, double snr_single);
double radar_ncoherent_integration_gain(size_t n_pulses);

/* ������ Utility ���������������������������������������������������������������������������������������������������������������� */

static inline double lin2db(double x) {
    return (x > 0.0) ? 10.0 * log10(x) : -INFINITY;
}
static inline double db2lin(double db) {
    return pow(10.0, db / 10.0);
}

/* ─── L7: Application-Level Helpers (declarations matching new src functions) ── */

double radar_power_aperture_product(const radar_params_t *params);
double radar_search_frame_time(double solid_angle_sr,
                               const radar_params_t *params,
                               double rcs_m2, double snr_min);
double radar_integration_gain_db(size_t n_pulses, int coherent);
double radar_min_detectable_rcs(const radar_params_t *params,
                                 double range_m, double snr_min);
double radar_atmospheric_loss_db(double freq_ghz, double range_km);
double radar_horizon_range_km(double h_radar_m, double h_target_m);

/* ─── L8: Advanced Radar Concepts ────────────────────────────────────── */

/**
 * Radar waveform agility modes.
 *
 * Modern radars adapt waveform parameters pulse-to-pulse
 * to optimize detection and tracking performance.
 */
typedef enum {
    AGILITY_NONE,
    AGILITY_FREQUENCY,       /**< Frequency hopping (ECCM) */
    AGILITY_PRF,             /**< PRF staggering (blind speed mitigation) */
    AGILITY_WAVEFORM,        /**< Waveform diversity (LFM + phase-coded) */
    AGILITY_POLARIZATION,    /**< Polarization agility */
    AGILITY_FULL             /**< Full cognitive radar adaptation */
} radar_agility_mode_t;

/**
 * Cognitive radar state.
 *
 * Perception-action cycle: sense environment, learn from
 * observations, adapt waveform and processing accordingly.
 */
typedef struct {
    radar_agility_mode_t agility;
    double current_snr_db;
    double interference_power_db;
    int target_present;
    size_t waveform_index;
    double adaptation_rate;
} cognitive_radar_state_t;

/**
 * Initialize cognitive radar state.
 */
int cognitive_radar_init(cognitive_radar_state_t *state);

/**
 * Adapt cognitive radar parameters based on observed SNR.
 */
int cognitive_radar_adapt(cognitive_radar_state_t *state,
                          double measured_snr_db);

/* ─── L7: Radar Operating Bands ──────────────────────────────────────── */

/**
 * IEEE radar band designations and typical applications.
 *
 * Band    Frequency (GHz)    Applications
 * ----    ---------------    ------------
 * HF      0.003-0.03         OTH (over-the-horizon) radar
 * VHF     0.03-0.3           Early warning, foliage penetration
 * UHF     0.3-1.0            Air surveillance, ballistic missile warning
 * L       1.0-2.0            Long-range air surveillance, ATC
 * S       2.0-4.0            Medium-range surveillance, weather radar
 * C       4.0-8.0            Fire control, SAR, weather
 * X       8.0-12.0           Fire control, maritime, airborne radar
 * Ku      12.0-18.0          High-resolution SAR, satellite altimetry
 * K       18.0-27.0          Automotive radar (24 GHz), police radar
 * Ka      27.0-40.0          Automotive radar (35/38 GHz), missile seeker
 * V       40.0-75.0          Automotive radar (60/77 GHz), imaging
 * W       75.0-110.0         High-resolution imaging, security screening
 * mm      110.0-300.0        Experimental, THz imaging
 */
typedef enum {
    RADAR_BAND_HF,
    RADAR_BAND_VHF,
    RADAR_BAND_UHF,
    RADAR_BAND_L,
    RADAR_BAND_S,
    RADAR_BAND_C,
    RADAR_BAND_X,
    RADAR_BAND_KU,
    RADAR_BAND_K,
    RADAR_BAND_KA,
    RADAR_BAND_V,
    RADAR_BAND_W,
    RADAR_BAND_MM
} radar_band_t;

/**
 * Get center frequency for a given radar band.
 */
double radar_band_center_freq(radar_band_t band);

/**
 * Get wavelength for a given radar band.
 */
double radar_band_wavelength(radar_band_t band);

/**
 * Get typical antenna gain range for a given radar band (dB).
 */
double radar_band_typical_gain_db(radar_band_t band);

#endif /* RADAR_CORE_H */
