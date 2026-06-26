/**
 * @file    ir_image.h
 * @brief   Thermal Image Representation and Processing
 *
 * Knowledge: L1 (thermal image, temperature matrix), L2 (AGC, DDE),
 *            L5 (histogram, spatial filtering)
 *
 * Thermal images encode scene temperature as pixel intensity.
 * Processing includes contrast enhancement, dynamic range compression,
 * and detail enhancement for optimal visualization.
 *
 * References:
 *   Gonzalez, R.C. & Woods, R.E. (2018) "Digital Image Processing", 4th Ed.
 *   Poynton, C. (2012) "Digital Video and HD: Algorithms and Interfaces"
 */

#ifndef IR_IMAGE_H
#define IR_IMAGE_H

#include <stddef.h>
#include <stdint.h>

/** Thermal image frame */
typedef struct {
    int     rows;
    int     cols;
    double *data;          /**< temperature or radiance values [rows*cols] */
    double  min_val;
    double  max_val;
    double  mean_val;
    double  std_val;
} ir_thermal_image_t;

/** Display image (8-bit grayscale) */
typedef struct {
    int      rows;
    int      cols;
    uint8_t *pixels;       /**< 0-255 grayscale [rows*cols] */
    double   temp_min;     /**< temperature mapped to 0 */
    double   temp_max;     /**< temperature mapped to 255 */
} ir_display_image_t;

/** Image statistics */
typedef struct {
    double min;
    double max;
    double mean;
    double median;
    double stddev;
    int    n_pixels;
} ir_image_stats_t;

/** Histogram */
typedef struct {
    int    n_bins;
    int   *bins;           /**< count per bin [n_bins] */
    double bin_min;
    double bin_max;
    double bin_width;
} ir_histogram_t;

/* Image creation and management */
ir_thermal_image_t* ir_image_create(int rows, int cols);
void ir_image_destroy(ir_thermal_image_t *img);
int ir_image_copy(const ir_thermal_image_t *src, ir_thermal_image_t *dst);
int ir_image_set_constant(ir_thermal_image_t *img, double value);

/* Statistics and analysis */
ir_image_stats_t ir_image_compute_stats(const ir_thermal_image_t *img);
ir_histogram_t* ir_image_histogram(const ir_thermal_image_t *img, int n_bins);
void ir_histogram_destroy(ir_histogram_t *hist);
double ir_image_median(const ir_thermal_image_t *img);

/* Contrast enhancement */
/** Automatic Gain Control: linear mapping from [min,max] to [0,255] */
int ir_image_agc_linear(const ir_thermal_image_t *img,
                         ir_display_image_t *disp,
                         double clip_percent_low,
                         double clip_percent_high);

/** Histogram equalization for enhanced contrast */
int ir_image_histogram_equalization(const ir_thermal_image_t *img,
                                     ir_display_image_t *disp);

/** Plateau histogram equalization (CLAHE-like) */
int ir_image_plateau_equalization(const ir_thermal_image_t *img,
                                   ir_display_image_t *disp,
                                   double plateau_limit);

/* Detail enhancement */
/** Digital Detail Enhancement (DDE) using bilateral-filter-like approach */
int ir_image_dde(const ir_thermal_image_t *img,
                  ir_display_image_t *disp,
                  double sigma_range,
                  double sigma_spatial);

/** High-pass sharpen filter */
int ir_image_sharpen(const ir_thermal_image_t *img,
                      ir_thermal_image_t *out,
                      double strength);

/* Spatial filtering */
int ir_image_median_filter(const ir_thermal_image_t *img,
                            ir_thermal_image_t *out,
                            int kernel_size);

int ir_image_gaussian_filter(const ir_thermal_image_t *img,
                              ir_thermal_image_t *out,
                              double sigma);

/** Temperature difference between two images (for change detection) */
int ir_image_difference(const ir_thermal_image_t *img1,
                         const ir_thermal_image_t *img2,
                         ir_thermal_image_t *diff);

/** Extract temperature matrix from raw data */
int ir_image_extract_temperature(const double *raw_data,
                                  int rows, int cols,
                                  double gain, double offset,
                                  ir_thermal_image_t *temp_img);

/** Region of Interest (ROI) statistics */
ir_image_stats_t ir_image_roi_stats(const ir_thermal_image_t *img,
                                     int row_start, int col_start,
                                     int row_end, int col_end);

/** Maximum temperature location (hot spot detection) */
int ir_image_find_hotspot(const ir_thermal_image_t *img,
                           int *row_out, int *col_out, double *temp_out);

#endif /* IR_IMAGE_H */
