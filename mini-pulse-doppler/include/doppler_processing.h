#ifndef DOPPLER_PROCESSING_H
#define DOPPLER_PROCESSING_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include "pulse_doppler.h"

/* ---------------------------------------------------------------------------
 * L1: Doppler Processing Definitions
 * --------------------------------------------------------------------------- */

/** Doppler processing mode */
typedef enum {
    DOPPLER_MODE_STANDARD = 0,
    DOPPLER_MODE_STAGGERED_PRF,
    DOPPLER_MODE_BURST,
    DOPPLER_MODE_CONTINUOUS_WAVE
} doppler_mode_t;

/** Window function type for Doppler sidelobe control */
typedef enum {
    WINDOW_RECTANGULAR = 0,
    WINDOW_HAMMING     = 1,
    WINDOW_HANNING     = 2,
    WINDOW_BLACKMAN    = 3,
    WINDOW_BLACKMAN_HARRIS = 4,
    WINDOW_FLATTOP     = 5,
    WINDOW_CHEBYSHEV   = 6,
    WINDOW_KAISER      = 7,
    WINDOW_TAYLOR      = 8
} window_type_t;

/** Window function descriptor */
typedef struct {
    window_type_t type;
    size_t length;
    double *coefficients;
    double coherent_gain;
    double scalloping_loss;
    double worst_case_processing_loss;
    double peak_sidelobe_db;
    double _6db_bandwidth_bins;
    double _3db_bandwidth_bins;
} window_func_t;

/** Doppler FFT processing parameters */
typedef struct {
    uint32_t fft_size;
    uint32_t num_pulses;
    double prf_hz;
    double doppler_resolution_hz;
    double max_velocity_mps;
    double velocity_resolution_mps;
    double wavelength_m;
    int zero_pad_enabled;
    window_type_t window;
    int remove_center_clutter;
} doppler_fft_params_t;

/** Doppler velocity spectrum */
typedef struct {
    double *velocity_bins;
    double *magnitude_spectrum;
    double *phase_spectrum;
    uint32_t num_bins;
    double velocity_resolution;
    double max_velocity;
} doppler_spectrum_t;

/** MTI filter structure (Moving Target Indication) */
typedef enum {
    MTI_SINGLE_DELAY = 0,
    MTI_DOUBLE_DELAY,
    MTI_TRIPLE_DELAY,
    MTI_FEEDBACK,
    MTI_OPTIMAL
} mti_filter_type_t;

typedef struct {
    mti_filter_type_t type;
    double complex *coefficients;
    size_t filter_order;
    double clutter_attenuation_db;
    double improvement_factor_db;
    double notch_width_hz;
} mti_filter_t;

/** Staggered PRF parameters for Doppler ambiguity resolution */
typedef struct {
    double *prf_values;
    uint32_t num_prfs;
    double *velocity_estimates;
    double unambiguous_velocity;
} staggered_prf_t;

/* ---------------------------------------------------------------------------
 * L5: Doppler Processing Algorithms
 * --------------------------------------------------------------------------- */

int window_function_create(window_type_t type, size_t length,
                           double param, window_func_t *win);

void window_function_free(window_func_t *win);

int window_function_apply(const window_func_t *win,
                          double complex *data,
                          size_t length);

/** L5: FFT-based Doppler processing along slow-time dimension */
int doppler_fft_process(const double complex *pulse_matrix,
                        uint32_t num_range_bins,
                        uint32_t num_pulses,
                        const doppler_fft_params_t *params,
                        range_doppler_map_t *rdmap);

int doppler_spectrum_compute(const double complex *slow_time_data,
                             uint32_t num_pulses,
                             double prf_hz,
                             double wavelength_m,
                             window_type_t window,
                             doppler_spectrum_t *spectrum);

void doppler_spectrum_free(doppler_spectrum_t *spectrum);

/** L1: Fundamental Doppler shift formula: fd = 2*v/λ */
int velocity_from_doppler(double doppler_hz,
                          double wavelength_m,
                          double *velocity_mps);

/** L1: Inverse: v → fd */
int doppler_from_velocity(double velocity_mps,
                          double wavelength_m,
                          double *doppler_hz);

/** L2: Maximum unambiguous velocity: vmax = λ*PRF/4 */
int max_unambiguous_velocity(double prf_hz,
                             double wavelength_m,
                             double *vmax_mps);

/** L2: Blind speed computation: v_blind = k*λ*PRF/2, k=1,2,3,... */
int blind_speed_compute(double prf_hz,
                        double wavelength_m,
                        uint32_t k,
                        double *blind_speed_mps);

/** L7: Staggered PRF design for Doppler ambiguity resolution */
int stagger_prf_design(double wavelength_m,
                       double max_velocity_mps,
                       uint32_t num_prfs,
                       staggered_prf_t *stagger);

void stagger_prf_free(staggered_prf_t *stagger);

/** L7: Chinese Remainder Theorem for resolving Doppler ambiguities */
int resolve_doppler_ambiguity_crt(const double *ambiguous_frequencies,
                                  const double *prf_values,
                                  uint32_t num_prfs,
                                  double *true_frequency_hz);

/** L2: MTI filter design for clutter suppression */
int mti_filter_design(mti_filter_type_t type,
                      double prf_hz,
                      double clutter_width_hz,
                      mti_filter_t *filter);

int mti_filter_apply(const mti_filter_t *filter,
                     double complex *slow_time_data,
                     size_t length,
                     double complex *output);

void mti_filter_free(mti_filter_t *filter);

/** L8: Adaptive Doppler processing (STAP simplification) */
int adaptive_doppler_filter(const double complex *data_matrix,
                            uint32_t num_range_bins,
                            uint32_t num_pulses,
                            uint32_t num_guard_cells,
                            double complex *filtered);

#endif /* DOPPLER_PROCESSING_H */
