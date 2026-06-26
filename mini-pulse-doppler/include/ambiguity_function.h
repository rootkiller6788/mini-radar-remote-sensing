#ifndef AMBIGUITY_FUNCTION_H
#define AMBIGUITY_FUNCTION_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

/* ---------------------------------------------------------------------------
 * L1/L3: Ambiguity Function Definitions
 *
 * The ambiguity function χ(τ, fd) characterizes a radar waveform's
 * ability to resolve targets in both range and Doppler simultaneously.
 *
 * Definition (Woodward, 1953):
 *   χ(τ, fd) = ∫_{-∞}^{∞} s(t) s*(t - τ) exp(-j2π fd t) dt
 *
 *   |χ(τ, fd)|² is the ambiguity surface.
 *
 * Key Properties:
 *   1. |χ(0,0)|² = (2E)², the maximum value
 *   2. ∫∫ |χ(τ,fd)|² dτ dfd = (2E)² (volume invariance)
 *   3. |χ(-τ, -fd)| = |χ(τ, fd)| (symmetry)
 *
 * Reference: Woodward, P.M. "Probability and Information Theory with
 *            Applications to Radar" (1953)
 * --------------------------------------------------------------------------- */

/** Ambiguity function result structure */
typedef struct {
    double complex *surface;
    uint32_t num_delays;
    uint32_t num_dopplers;
    double delay_resolution_s;
    double doppler_resolution_hz;
    double max_delay_s;
    double max_doppler_hz;
} ambiguity_surface_t;

/** Ambiguity function metrics */
typedef struct {
    double peak_value;
    double peak_sidelobe_ratio_db;
    double integrated_sidelobe_ratio_db;
    double _3db_delay_width_s;
    double _3db_doppler_width_hz;
    double delay_doppler_volume;
    double q_function_value;
} ambiguity_metrics_t;

/** Range-Doppler coupling analysis for LFM */
typedef struct {
    double coupling_factor;
    double range_error_per_doppler_hz;
    double doppler_error_per_delay_s;
    int is_range_doppler_coupled;
} rd_coupling_t;

/* ---------------------------------------------------------------------------
 * L5: Ambiguity Function Computation Algorithms
 * --------------------------------------------------------------------------- */

/**
 * L5: Compute the full 2-D ambiguity function.
 *
 * Uses direct computation with FFT for the Doppler dimension:
 *   For each delay τ, the inner integral over t is a Fourier transform
 *   evaluated at fd.
 *
 * Complexity: O(N_delay × N_doppler × N_signal) for direct method,
 *             O(N_delay × N_signal × log N_signal) with FFT optimization
 */
int ambiguity_function_compute(const double complex *waveform,
                               size_t signal_length,
                               double sampling_rate,
                               uint32_t num_delays,
                               uint32_t num_dopplers,
                               double max_doppler_hz,
                               ambiguity_surface_t *result);

/** L5: Range cut (zero-Doppler) of the ambiguity function */
int ambiguity_cut_range(const ambiguity_surface_t *surface,
                        double *range_cut);

/** L5: Doppler cut (zero-delay) of the ambiguity function */
int ambiguity_cut_doppler(const ambiguity_surface_t *surface,
                          double *doppler_cut);

/** L5: Compute ambiguity function metrics */
int ambiguity_metrics_compute(const ambiguity_surface_t *surface,
                              double signal_energy,
                              ambiguity_metrics_t *metrics);

/** L5: Compare ambiguity performance of two waveforms */
int ambiguity_compare(const ambiguity_surface_t *af1,
                      const ambiguity_surface_t *af2,
                      double *similarity_metric);

/** L3: Range-Doppler coupling strength for LFM waveforms */
int rd_coupling_analyze(double chirp_rate,
                        double carrier_frequency_hz,
                        rd_coupling_t *coupling);

/** L3: Q-function: inverse of the ambiguity function volume property */
int q_function_compute(const ambiguity_surface_t *surface,
                       double *q_value);

/** Memory management */
int ambiguity_surface_alloc(uint32_t num_delays,
                            uint32_t num_dopplers,
                            ambiguity_surface_t *surface);

void ambiguity_surface_free(ambiguity_surface_t *surface);

#endif /* AMBIGUITY_FUNCTION_H */
