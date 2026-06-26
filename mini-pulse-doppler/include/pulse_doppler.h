#ifndef PULSE_DOPPLER_H
#define PULSE_DOPPLER_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include "radar_waveform.h"

/* ---------------------------------------------------------------------------
 * L1: Core Definitions
 * --------------------------------------------------------------------------- */

/** Detection decision: H0 = noise only, H1 = target present */
typedef enum {
    DETECTION_H0 = 0,
    DETECTION_H1 = 1
} detection_decision_t;

/** Target model: Swerling cases for RCS fluctuation statistics */
typedef enum {
    SWERLING_0 = 0,
    SWERLING_1 = 1,
    SWERLING_2 = 2,
    SWERLING_3 = 3,
    SWERLING_4 = 4
} swerling_case_t;

/**
 * L1: Radar detection parameters
 *
 * Probability of detection (Pd) and probability of false alarm (Pfa)
 * are the fundamental metrics of radar detection performance.
 */
typedef struct {
    double pfa;
    double pd;
    double snr_db;
    double threshold;
    uint32_t num_samples;
    uint32_t num_pulses;
} radar_detection_params_t;

/**
 * L1: Range-Doppler cell
 *
 * The fundamental unit of radar processing: a specific range bin
 * and Doppler frequency bin pair.
 */
typedef struct {
    uint32_t range_bin;
    uint32_t doppler_bin;
    double  range_meters;
    double  doppler_hz;
    double  velocity_mps;
    double  magnitude;
    double  phase_rad;
} range_doppler_cell_t;

/**
 * L1: Range-Doppler map
 *
 * The primary output of pulse Doppler processing: a 2-D matrix
 * where rows = range bins, columns = Doppler bins.
 */
typedef struct {
    double complex *data;
    uint32_t num_range_bins;
    uint32_t num_doppler_bins;
    double range_resolution_m;
    double doppler_resolution_hz;
    double velocity_resolution_mps;
    double max_unambiguous_range_m;
    double max_unambiguous_velocity_mps;
} range_doppler_map_t;

/**
 * L1: Coherent Processing Interval (CPI) parameters
 */
typedef struct {
    uint32_t num_pulses;
    double pri_seconds;
    double prf_hz;
    double center_frequency_hz;
    double wavelength_m;
    double coherent_gain_db;
    double integration_time_s;
} cpi_params_t;

/**
 * L1: Target detection report
 */
typedef struct {
    uint32_t target_id;
    double range_m;
    double velocity_mps;
    double doppler_hz;
    double snr_db;
    double azimuth_deg;
    double elevation_deg;
    double rcs_dbsm;
    detection_decision_t detected;
} target_detection_t;

/**
 * L2: Matched filter structure
 */
typedef struct {
    double complex *coefficients;
    size_t filter_length;
    double processing_gain_db;
    int is_normalized;
} matched_filter_t;

/* ---------------------------------------------------------------------------
 * L5: Core Pulse Doppler Processing Algorithms
 * --------------------------------------------------------------------------- */

int matched_filter_init(matched_filter_t *mf,
                        const double complex *reference,
                        size_t length);

void matched_filter_free(matched_filter_t *mf);

int matched_filter_apply(const matched_filter_t *mf,
                         const double complex *input,
                         size_t input_len,
                         double complex *output);

int matched_filter_apply_freq_domain(const double complex *reference,
                                     const double complex *input,
                                     size_t length,
                                     double complex *output);

int pulse_compress_lfm(const double complex *chirp_reference,
                       const double complex *received_signal,
                       size_t signal_length,
                       size_t ref_length,
                       double complex *compressed);

int range_gate_extract(const double complex *pulse_data,
                       uint32_t num_range_bins,
                       uint32_t num_pulses,
                       uint32_t gate_start,
                       uint32_t gate_width,
                       double complex *gate_data);

int coherent_integration(double complex *pulse_matrix,
                         uint32_t num_range_bins,
                         uint32_t num_pulses,
                         double complex *integrated);

int noncoherent_integration(const double *magnitude_matrix,
                            uint32_t num_range_bins,
                            uint32_t num_pulses,
                            double *integrated);

int compute_coherent_gain(uint32_t num_pulses,
                          double *gain_db);

int compute_noncoherent_gain(uint32_t num_pulses,
                             double *gain_db);

int compute_integration_loss(uint32_t num_pulses,
                             double pfa,
                             double pd,
                             double *loss_db);

int detection_threshold_compute(double pfa,
                                uint32_t num_samples,
                                int noise_power_known,
                                double noise_power,
                                double *threshold);

/** L4: Albersheim's equation for Pd/Pfa/SNR relationship */
int albersheim_equation(double pfa, double pd,
                        uint32_t num_pulses,
                        double *required_snr_db);

/** L4: Shnidman's equation (improved Albersheim for Swerling targets) */
int shnidman_equation(double pfa, double pd,
                      uint32_t num_pulses,
                      swerling_case_t swerling,
                      double *required_snr_db);

int target_detection_report_init(target_detection_t *report);

int cpi_params_from_waveform(const radar_waveform_params_t *wf,
                             uint32_t num_pulses,
                             cpi_params_t *cpi);

int range_doppler_map_allocate(range_doppler_map_t *rdmap,
                               uint32_t num_range,
                               uint32_t num_doppler);

void range_doppler_map_free(range_doppler_map_t *rdmap);

int range_doppler_map_find_peaks(const range_doppler_map_t *rdmap,
                                 double threshold,
                                 uint32_t max_peaks,
                                 range_doppler_cell_t *peaks,
                                 uint32_t *num_found);

#endif /* PULSE_DOPPLER_H */
