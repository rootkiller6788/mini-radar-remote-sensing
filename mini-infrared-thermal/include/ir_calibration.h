/**
 * @file    ir_calibration.h
 * @brief   Infrared Detector Calibration - NUC, Bad Pixel, Temperature Cal
 *
 * Knowledge: L1 (NUC, bad pixel), L5 (two-point calibration, shutter correction)
 *
 * Calibration is essential for converting detector output to radiometric
 * quantities. IR focal plane arrays (FPAs) have pixel-to-pixel non-uniformity
 * that must be corrected for accurate temperature measurement and imaging.
 *
 * References:
 *   Perry, D.L. & Dereniak, E.L. (1993) "Linear theory of nonuniformity
 *     correction in IR staring sensors", Optical Engineering, 32(8)
 *   Scribner, D.A. et al. (1990) "IR FPA nonuniformity correction",
 *     SPIE Vol. 1308
 */

#ifndef IR_CALIBRATION_H
#define IR_CALIBRATION_H

#include <stddef.h>
#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "ir_core.h"

/* Maximum FPA dimensions */
#define IR_CAL_MAX_ROWS    2048
#define IR_CAL_MAX_COLS    2048

/* ========================================================================
 * L1: Calibration Data Structures
 * ======================================================================== */

/** Gain and offset correction coefficients per pixel */
typedef struct {
    double gain;        /**< multiplicative correction (ideal: 1.0) */
    double offset;      /**< additive correction (ideal: 0.0) */
    int    is_bad;       /**< 1 if pixel is non-responsive */
} ir_cal_pixel_coeff_t;

/** Two-point calibration data */
typedef struct {
    int     rows;
    int     cols;
    double  T_low_K;           /**< low reference temperature [K] */
    double  T_high_K;          /**< high reference temperature [K] */
    double  *raw_low;          /**< raw values at T_low [rows*cols] */
    double  *raw_high;         /**< raw values at T_high [rows*cols] */
    ir_cal_pixel_coeff_t *coeffs; /**< computed coefficients [rows*cols] */
    double *corrected;         /**< corrected output buffer [rows*cols] */
} ir_two_point_cal_t;

/** Multi-point calibration with piecewise linear correction */
typedef struct {
    int      rows;
    int      cols;
    int      n_points;         /**< number of calibration temperatures */
    double   *T_points_K;      /**< calibration temperatures [n_points] */
    double  **raw_values;      /**< raw data [n_points][rows*cols] */
} ir_multi_point_cal_t;

/** Bad pixel map */
typedef struct {
    int     rows;
    int     cols;
    uint8_t *map;              /**< 1 = bad, 0 = good [rows*cols] */
    int     n_bad;             /**< total number of bad pixels */
} ir_bad_pixel_map_t;

/** Non-uniformity correction (NUC) state */
typedef struct {
    int     rows;
    int     cols;
    double *gain;              /**< gain correction map [rows*cols] */
    double *offset;            /**< offset correction map [rows*cols] */
    double *current_frame;     /**< most recent corrected frame */
} ir_nuc_state_t;

/* ========================================================================
 * L5: Non-Uniformity Correction Algorithms
 * ======================================================================== */

/** One-point offset correction.
 *  V_corrected[i] = V_raw[i] - offset[i]
 *  Offset is computed from a uniform reference (shutter). */
int ir_nuc_one_point_correction(ir_nuc_state_t *nuc,
                                 const double *raw_frame,
                                 double *corrected_out);

/** Two-point gain and offset correction.
 *  V_corrected = gain * (V_raw - offset)
 *  gain = (V_high_ref - V_low_ref) / (V_high_raw - V_low_raw)
 *  offset = V_low_raw */
int ir_cal_two_point_compute(const double *raw_low, const double *raw_high,
                              double T_low_K, double T_high_K,
                              int rows, int cols,
                              ir_cal_pixel_coeff_t *coeffs_out);

/** Apply two-point correction to a frame */
int ir_cal_two_point_apply(const ir_cal_pixel_coeff_t *coeffs,
                            const double *raw_frame,
                            int rows, int cols,
                            double *corrected_out);

/** Multi-point piecewise linear calibration.
 *  Uses linear interpolation between nearest calibration points. */
int ir_cal_multi_point_apply(const ir_multi_point_cal_t *cal,
                              const double *raw_value,
                              int rows, int cols,
                              double *corrected_radiance_out);

/** Scene-based NUC using constant-statistics assumption.
 *  Adjusts gain/offset to equalize mean and variance across pixels. */
int ir_nuc_scene_based_update(ir_nuc_state_t *nuc,
                               const double *raw_frame,
                               int rows, int cols,
                               double learning_rate);

/* ========================================================================
 * L5: Bad Pixel Detection and Replacement
 * ======================================================================== */

/** Detect bad pixels from uniformity data.
 *  Pixels deviating > n_sigma * stddev from the mean are flagged. */
int ir_bad_pixel_detect(const double *raw_frame, int rows, int cols,
                         double n_sigma, ir_bad_pixel_map_t *map_out);

/** Replace bad pixels by median of nearest good neighbors.
 *  Uses a 3x3 kernel; skips other bad pixels in neighborhood. */
int ir_bad_pixel_replace(double *frame, int rows, int cols,
                          const ir_bad_pixel_map_t *map);

/** Dead pixel detection: pixels with responsivity below threshold */
int ir_dead_pixel_detect(const ir_cal_pixel_coeff_t *coeffs,
                          int rows, int cols,
                          double gain_threshold,
                          ir_bad_pixel_map_t *map_out);

/** Drifting pixel detection: pixels with excessive temporal noise */
int ir_drifting_pixel_detect(const double *temporal_stddev,
                              int rows, int cols,
                              double noise_threshold,
                              ir_bad_pixel_map_t *map_out);

/* ========================================================================
 * L5: Temperature Calibration
 * ======================================================================== */

/** Convert raw digital numbers to scene temperature [K].
 *  Uses Planck curve fitting with calibration coefficients. */
double ir_cal_raw_to_temperature(double raw_value,
                                  double gain_coeff,
                                  double offset_coeff,
                                  double lambda_um);

/** Shutter correction: update offset map from shutter frame.
 *  Removes temporal drift of the FPA offset. */
int ir_cal_shutter_correction(ir_nuc_state_t *nuc,
                               const double *shutter_frame,
                               int rows, int cols);

/** Compute residual non-uniformity after correction.
 *  Returns RMS of (corrected - mean) / mean as percentage. */
double ir_nuc_residual_nuc(const double *corrected, int rows, int cols);

#endif /* IR_CALIBRATION_H */
