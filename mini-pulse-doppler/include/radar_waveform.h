#ifndef RADAR_WAVEFORM_H
#define RADAR_WAVEFORM_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

/** L1: Waveform modulation type enum */
typedef enum {
    WAVEFORM_RECT_PULSE = 0,
    WAVEFORM_LFM_UP,
    WAVEFORM_LFM_DOWN,
    WAVEFORM_BARKER,
    WAVEFORM_NLFM,
    WAVEFORM_COSTAS,
    WAVEFORM_STEPPED_FREQ
} waveform_type_t;

/** L1: Core radar waveform parameter definitions */
typedef struct {
    double pulse_width;
    double bandwidth;
    double center_frequency;
    double prf;
    double pri;
    double duty_cycle;
    double sampling_rate;
    uint32_t num_samples;
    waveform_type_t type;
} radar_waveform_params_t;

/** L1: Barker code parameters */
typedef struct {
    uint32_t code_length;
    double chip_rate;
    double chip_duration;
    int8_t *code_sequence;
} barker_code_t;

/** L2: LFM chirp parameters */
typedef struct {
    double chirp_rate;
    double time_bw_product;
    int    is_upchirp;
} lfm_params_t;

/** L8: Non-linear FM parameters */
typedef struct {
    double *freq_profile;
    double *amp_profile;
    uint32_t num_profile_pts;
    double taylor_nbar;
    double peak_sidelobe_db;
} nlfm_params_t;

/* L5: Waveform Generation Algorithms */

int rect_pulse_generate(const radar_waveform_params_t *params,
                        double complex *buffer);

int lfm_chirp_generate(const radar_waveform_params_t *params,
                       const lfm_params_t *chirp,
                       double complex *buffer);

int barker_code_generate(const radar_waveform_params_t *params,
                         const barker_code_t *barker,
                         double complex *buffer);

int pulse_train_generate(const radar_waveform_params_t *params,
                         uint32_t num_pulses,
                         const double complex *single_pulse,
                         double complex *buffer);

int waveform_autocorrelation(const double complex *signal,
                             size_t N,
                             double *acorr);

int waveform_cross_ambiguity_1d(const double complex *tx,
                                const double complex *rx,
                                size_t length,
                                double fd,
                                double *result);

int waveform_instantaneous_freq(const double complex *signal,
                                size_t N,
                                double *freq,
                                double fs);

int waveform_apply_window(double complex *signal,
                          size_t N,
                          int window_type,
                          double taylor_nbar,
                          double taylor_sll_db);

int waveform_rms_bandwidth(const double complex *signal,
                           size_t N,
                           double fs,
                           double *bw_rms);

/** L5: Peak sidelobe level from autocorrelation. PSL = 20*log10(max_sidelobe/peak) dB */
int compute_peak_sidelobe_level(const double *autocorr, size_t N, double *psl_db);

/** L5: Integrated sidelobe level ratio. ISLR = 10*log10(E_sidelobes/E_mainlobe) dB */
int compute_integrated_sidelobe_level(const double *autocorr, size_t N, double *islr_db);

/** L1: Rayleigh range resolution: Delta_R = c/(2*B) */
double compute_range_resolution(double bandwidth_hz);

/** L1: Doppler resolution: Delta_fd = PRF/M */
double compute_doppler_resolution(double prf_hz, uint32_t num_pulses);

/** L1: Velocity resolution: Delta_v = lambda*PRF/(2*M) */
double compute_velocity_resolution(double wavelength_m, double prf_hz, uint32_t num_pulses);

/** L2: Pulse compression gain: G = 10*log10(tau*B) dB */
double compute_pulse_compression_gain(double pulse_width_s, double bandwidth_hz);

/** L1: Duty cycle: d = tau * PRF */
double compute_duty_cycle(double pulse_width_s, double prf_hz);

/** L1: Peak power from average power: P_peak = P_avg / duty */
double compute_peak_power(double average_power_w, double duty_cycle);

/** L1: Average power from peak: P_avg = P_peak * duty */
double compute_average_power(double peak_power_w, double duty_cycle);

/** L1: Energy per pulse: E = P_peak * tau */
double compute_energy_per_pulse(double peak_power_w, double pulse_width_s);

/** L3: DFT-based magnitude spectrum of a complex waveform */
int waveform_spectrum_magnitude(const double complex *signal, size_t N,
                                double *magnitude_spectrum);

/** L2: Costas frequency-hopping code waveform generation */
int costas_code_generate(const radar_waveform_params_t *params,
                         uint32_t code_length,
                         double complex *buffer);

/** L7: Generate stepped-frequency CW waveform for high-range-resolution profiling */
int stepped_frequency_generate(const radar_waveform_params_t *params,
                               uint32_t num_steps,
                               double frequency_step_hz,
                               double complex *buffer);

/** L7: Compute signal-to-clutter ratio for ground clutter */
double compute_signal_to_clutter_ratio(double target_rcs_dbsm,
                                       double clutter_rcs_per_area_dbsm,
                                       double range_resolution_m,
                                       double azimuth_resolution_m);

/** L1: Compute unambiguous range from PRF: Rmax = c/(2*PRF) */
double compute_unambiguous_range(double prf_hz);

/** L2: Compute blind speed: v_blind = k*lambda*PRF/2 */
double compute_blind_speed(double prf_hz, double wavelength_m, uint32_t k);

#endif /* RADAR_WAVEFORM_H */
