/**
 * @file    ir_detection.h
 * @brief   Infrared Target Detection Algorithms
 *
 * L5 - Algorithms: CA-CFAR, Otsu thresholding, matched filter
 * L7 - Applications: IRST (Infrared Search and Track), missile warning
 */

#ifndef IR_DETECTION_H
#define IR_DETECTION_H

#include <stddef.h>

/* CA-CFAR (Cell-Averaging Constant False Alarm Rate) */
int ir_cfar_detect_1d(const double *signal, int length,
                       int test_idx, int guard_cells,
                       int reference_cells, double p_fa,
                       double *threshold_out);

/* Background subtraction */
int ir_background_subtract(const double *frame, const double *background,
                            int length, double *diff_out);
int ir_background_update(const double *frame, double *background,
                          int length, double alpha);
int ir_threshold_segment(const double *data, int length,
                          double threshold, int *binary_out);

/* Otsu automatic thresholding */
double ir_otsu_threshold(const double *data, int length, int n_bins);

/* Matched filter */
double ir_matched_filter_1d(const double *signal, int sig_len,
                             const double *template, int temp_len,
                             double *correlation_out);

/* Signal-to-Clutter Ratio */
double ir_scr_estimate(const double *data, int length,
                        int target_start, int target_end);

/* Detection probability */
double ir_detection_probability(double scr, double p_fa);

#endif /* IR_DETECTION_H */
