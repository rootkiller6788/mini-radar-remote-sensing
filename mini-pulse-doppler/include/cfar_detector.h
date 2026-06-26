#ifndef CFAR_DETECTOR_H
#define CFAR_DETECTOR_H

#include <stddef.h>
#include <stdint.h>

/* ---------------------------------------------------------------------------
 * L1: CFAR Detector Definitions
 *
 * Constant False Alarm Rate (CFAR) detection maintains a constant
 * probability of false alarm by adaptively estimating the local
 * noise/clutter power from surrounding range/Doppler cells.
 *
 * Reference: Rohling, "Radar CFAR Thresholding in Clutter and Multiple
 * Target Situations," IEEE Trans. AES, 1983.
 * --------------------------------------------------------------------------- */

/** CFAR algorithm type */
typedef enum {
    CFAR_CA = 0,
    CFAR_GO = 1,
    CFAR_SO = 2,
    CFAR_OS = 3,
    CFAR_OSCA = 4,
    CFAR_CAGO = 5,
    CFAR_CASO = 6,
    CFAR_TM = 7,
    CFAR_VI = 8
} cfar_type_t;

/** CFAR detector configuration */
typedef struct {
    cfar_type_t type;
    uint32_t num_reference_cells;
    uint32_t num_guard_cells;
    double pfa;
    double threshold_factor;
    double os_rank;
    double vi_threshold;
    double mr_threshold;
} cfar_config_t;

/** CFAR detection result for a single cell */
typedef struct {
    uint32_t cell_index;
    double cell_value;
    double estimated_noise_power;
    double threshold_value;
    int detected;
} cfar_detection_t;

/** 2-D CFAR detection map (range-Doppler plane) */
typedef struct {
    uint8_t *detection_map;
    double *threshold_map;
    double *noise_power_map;
    uint32_t num_range_bins;
    uint32_t num_doppler_bins;
    uint32_t num_detections;
} cfar_detection_map_t;

/* ---------------------------------------------------------------------------
 * L5: CFAR Detection Algorithms
 * --------------------------------------------------------------------------- */

/**
 * Configure CFAR detector with standard parameters.
 *
 * For CA-CFAR, the threshold factor α is computed as:
 *   α = N × (Pfa^(-1/N) - 1)
 * where N is the number of reference cells.
 *
 * For OS-CFAR, the k-th ordered statistic is used.
 */
int cfar_config_init(cfar_type_t type,
                     uint32_t num_reference,
                     uint32_t num_guard,
                     double pfa,
                     cfar_config_t *config);

/**
 * L5: Cell-Averaging CFAR — the fundamental CFAR algorithm.
 *
 * Estimates noise power as the arithmetic mean of reference cells.
 * Optimal in homogeneous Rayleigh clutter but degrades near clutter edges
 * and in multiple-target situations.
 *
 * Threshold: T = α × (1/N) × Σ_{i=1}^{N} x_i
 *
 * Complexity: O(N_cells × N_ref)
 */
int ca_cfar_detect(const double *data,
                   size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections,
                   size_t *num_detections);

/**
 * L5: Greatest-Of CFAR — mitigates clutter edge false alarms.
 *
 * Splits reference window into leading/lagging halves, takes max of the two.
 * Reduces false alarms at clutter edges but may mask closely-spaced targets.
 *
 * T = α × max(μ_leading, μ_lagging)
 */
int go_cfar_detect(const double *data,
                   size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections,
                   size_t *num_detections);

/**
 * L5: Smallest-Of CFAR — better detection of closely-spaced targets.
 *
 * Uses min of leading/lagging window averages. Better for resolving
 * multiple nearby targets but has higher false alarm rate at edges.
 *
 * T = α × min(μ_leading, μ_lagging)
 */
int so_cfar_detect(const double *data,
                   size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections,
                   size_t *num_detections);

/**
 * L5: Ordered-Statistic CFAR — robust in non-homogeneous clutter.
 *
 * Sorts reference cells and uses the k-th value as noise estimate.
 * Inherently immune to up to (N - k) interfering targets.
 *
 * T = α × x_(k)  where x_(k) is the k-th order statistic
 *
 * Complexity: O(N_cells × N_ref × log(N_ref)) due to sorting
 */
int os_cfar_detect(const double *data,
                   size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections,
                   size_t *num_detections);

/**
 * L8: Variability Index CFAR — adaptive algorithm selection.
 *
 * Uses VI statistic to detect non-homogeneity and automatically
 * switches between CA/GO/SO strategies per cell.
 *
 * VI = 1 + (σ² / μ²)
 */
int vi_cfar_detect(const double *data,
                   size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections,
                   size_t *num_detections);

/** 2-D CFAR scan over a range-Doppler map */
int cfar_2d_scan(const double *rd_map,
                 uint32_t num_range,
                 uint32_t num_doppler,
                 const cfar_config_t *config,
                 cfar_detection_map_t *result);

/** Compute the CFAR threshold factor α from Pfa and N */
int cfar_threshold_factor_compute(cfar_type_t type,
                                  double pfa,
                                  uint32_t num_reference,
                                  double os_rank,
                                  double *alpha);

/** Evaluate detection performance: Pd for given SNR and CFAR params */
int cfar_pd_compute(cfar_type_t type,
                    double snr_db,
                    double pfa,
                    uint32_t num_reference,
                    double *pd);

/** CFAR detection map management */
int cfar_detection_map_alloc(uint32_t num_range,
                              uint32_t num_doppler,
                              cfar_detection_map_t *map);

void cfar_detection_map_free(cfar_detection_map_t *map);

#endif /* CFAR_DETECTOR_H */
