/**
 * @file    ir_image.c
 * @brief   Thermal Image Processing - Implementation
 *
 * Implements thermal image representation, statistics, contrast enhancement
 * (AGC, histogram equalization), detail enhancement (DDE), and spatial filtering.
 *
 * References:
 *   Gonzalez & Woods (2018) "Digital Image Processing", 4th Ed.
 *   FLIR Systems (2015) "DDE - Digital Detail Enhancement Technical Note"
 */

#include "ir_image.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ========================================================================
 * Image Creation and Management
 * ======================================================================== */

ir_thermal_image_t* ir_image_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    ir_thermal_image_t *img = (ir_thermal_image_t*)malloc(sizeof(*img));
    if (!img) return NULL;
    img->rows = rows;
    img->cols = cols;
    img->data = (double*)calloc(rows * cols, sizeof(double));
    if (!img->data) { free(img); return NULL; }
    img->min_val = 0.0;
    img->max_val = 0.0;
    img->mean_val = 0.0;
    img->std_val = 0.0;
    return img;
}

void ir_image_destroy(ir_thermal_image_t *img) {
    if (img) { free(img->data); free(img); }
}

int ir_image_copy(const ir_thermal_image_t *src, ir_thermal_image_t *dst) {
    if (!src || !dst) return -1;
    if (src->rows != dst->rows || src->cols != dst->cols) return -1;
    int N = src->rows * src->cols;
    memcpy(dst->data, src->data, N * sizeof(double));
    dst->min_val = src->min_val;
    dst->max_val = src->max_val;
    dst->mean_val = src->mean_val;
    dst->std_val = src->std_val;
    return 0;
}

int ir_image_set_constant(ir_thermal_image_t *img, double value) {
    if (!img) return -1;
    int N = img->rows * img->cols;
    for (int i = 0; i < N; i++) img->data[i] = value;
    img->min_val = img->max_val = value;
    img->mean_val = value;
    img->std_val = 0.0;
    return 0;
}

/* ========================================================================
 * Image Statistics
 * ======================================================================== */

ir_image_stats_t ir_image_compute_stats(const ir_thermal_image_t *img) {
    ir_image_stats_t stats = {0};
    if (!img || img->rows <= 0 || img->cols <= 0) return stats;

    int N = img->rows * img->cols;
    double min_v = img->data[0], max_v = img->data[0];
    double sum = 0.0, sum2 = 0.0;

    for (int i = 0; i < N; i++) {
        double v = img->data[i];
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
        sum += v;
        sum2 += v * v;
    }

    stats.min = min_v;
    stats.max = max_v;
    stats.mean = sum / N;
    stats.stddev = sqrt(sum2 / N - stats.mean * stats.mean);
    stats.n_pixels = N;
    return stats;
}

double ir_image_median(const ir_thermal_image_t *img) {
    if (!img || img->rows <= 0 || img->cols <= 0) return 0.0;

    int N = img->rows * img->cols;
    double *sorted = (double*)malloc(N * sizeof(double));
    if (!sorted) return 0.0;
    memcpy(sorted, img->data, N * sizeof(double));

    /* Simple bubble sort for median (use qsort in production) */
    for (int i = 0; i < N - 1; i++)
        for (int j = i + 1; j < N; j++)
            if (sorted[i] > sorted[j]) {
                double tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
    double med = sorted[N / 2];
    free(sorted);
    return med;
}

/* ========================================================================
 * Automatic Gain Control (AGC) - Linear Mapping
 * ======================================================================== */

int ir_image_agc_linear(const ir_thermal_image_t *img,
                         ir_display_image_t *disp,
                         double clip_percent_low,
                         double clip_percent_high) {
    if (!img || !disp || img->rows != disp->rows || img->cols != disp->cols)
        return -1;

    ir_image_stats_t s = ir_image_compute_stats(img);
    double t_min = s.min + clip_percent_low * (s.max - s.min) / 100.0;
    double t_max = s.max - clip_percent_high * (s.max - s.min) / 100.0;
    if (fabs(t_max - t_min) < 1e-10) return -1;

    disp->temp_min = t_min;
    disp->temp_max = t_max;

    int N = img->rows * img->cols;
    for (int i = 0; i < N; i++) {
        double norm = (img->data[i] - t_min) / (t_max - t_min);
        if (norm < 0.0) norm = 0.0;
        if (norm > 1.0) norm = 1.0;
        disp->pixels[i] = (uint8_t)(norm * 255.0);
    }
    return 0;
}

/* ========================================================================
 * Histogram Equalization
 * ======================================================================== */

int ir_image_histogram_equalization(const ir_thermal_image_t *img,
                                     ir_display_image_t *disp) {
    if (!img || !disp || img->rows != disp->rows || img->cols != disp->cols)
        return -1;

    int N = img->rows * img->cols;
    int hist[256] = {0};

    ir_image_stats_t s = ir_image_compute_stats(img);
    double t_min = s.min, t_max = s.max;
    if (fabs(t_max - t_min) < 1e-10) return -1;

    /* Build histogram (256 bins) */
    for (int i = 0; i < N; i++) {
        int bin = (int)((img->data[i] - t_min) / (t_max - t_min) * 255.0);
        if (bin < 0) bin = 0;
        if (bin > 255) bin = 255;
        hist[bin]++;
    }

    /* Cumulative distribution function */
    int cdf[256];
    cdf[0] = hist[0];
    for (int i = 1; i < 256; i++) cdf[i] = cdf[i - 1] + hist[i];

    /* Map to equalized values */
    double scale = 255.0 / N;
    for (int i = 0; i < N; i++) {
        int bin = (int)((img->data[i] - t_min) / (t_max - t_min) * 255.0);
        if (bin < 0) bin = 0;
        if (bin > 255) bin = 255;
        disp->pixels[i] = (uint8_t)(cdf[bin] * scale);
    }
    disp->temp_min = t_min;
    disp->temp_max = t_max;
    return 0;
}

/* ========================================================================
 * Spatial Filtering
 * ======================================================================== */

int ir_image_median_filter(const ir_thermal_image_t *img,
                            ir_thermal_image_t *out,
                            int kernel_size) {
    if (!img || !out || img->rows != out->rows || img->cols != out->cols)
        return -1;
    if (kernel_size < 3 || kernel_size % 2 == 0) return -1;

    int rows = img->rows, cols = img->cols;
    int half = kernel_size / 2;

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            /* Gather kernel values */
            double *vals = (double*)malloc(kernel_size * kernel_size * sizeof(double));
            int nv = 0;
            for (int dr = -half; dr <= half; dr++) {
                for (int dc = -half; dc <= half; dc++) {
                    int nr = r + dr, nc = c + dc;
                    if (nr >= 0 && nr < rows && nc >= 0 && nc < cols)
                        vals[nv++] = img->data[nr * cols + nc];
                }
            }
            /* Sort for median */
            for (int i = 0; i < nv - 1; i++)
                for (int j = i + 1; j < nv; j++)
                    if (vals[i] > vals[j]) {
                        double t = vals[i]; vals[i] = vals[j]; vals[j] = t;
                    }
            out->data[r * cols + c] = vals[nv / 2];
            free(vals);
        }
    }
    return 0;
}

int ir_image_gaussian_filter(const ir_thermal_image_t *img,
                              ir_thermal_image_t *out,
                              double sigma) {
    if (!img || !out || img->rows != out->rows || img->cols != out->cols)
        return -1;
    if (sigma <= 0.0) return -1;
    int rows = img->rows, cols = img->cols;
    int ks = (int)(6.0 * sigma + 1.0);
    if (ks % 2 == 0) ks++;
    int half = ks / 2;
    double *kernel = (double*)malloc(ks * ks * sizeof(double));
    double sum = 0.0;
    for (int i = 0; i < ks; i++) {
        for (int j = 0; j < ks; j++) {
            double dx = i - half, dy = j - half;
            double v = exp(-(dx*dx + dy*dy) / (2.0 * sigma * sigma));
            kernel[i * ks + j] = v;
            sum += v;
        }
    }
    for (int i = 0; i < ks * ks; i++) kernel[i] /= sum;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double v = 0.0;
            for (int dr = -half; dr <= half; dr++) {
                for (int dc = -half; dc <= half; dc++) {
                    int nr = r + dr, nc = c + dc;
                    if (nr >= 0 && nr < rows && nc >= 0 && nc < cols)
                        v += img->data[nr*cols+nc] * kernel[(dr+half)*ks + (dc+half)];
                }
            }
            out->data[r * cols + c] = v;
        }
    }
    free(kernel);
    return 0;
}

/* ========================================================================
 * Hot Spot Detection
 * ======================================================================== */

int ir_image_find_hotspot(const ir_thermal_image_t *img,
                           int *row_out, int *col_out, double *temp_out) {
    if (!img || !row_out || !col_out || !temp_out) return -1;
    int N = img->rows * img->cols;
    if (N <= 0) return -1;
    int max_idx = 0;
    for (int i = 1; i < N; i++)
        if (img->data[i] > img->data[max_idx]) max_idx = i;
    *row_out = max_idx / img->cols;
    *col_out = max_idx % img->cols;
    *temp_out = img->data[max_idx];
    return 0;
}

int ir_image_extract_temperature(const double *raw_data,
                                  int rows, int cols,
                                  double gain, double offset,
                                  ir_thermal_image_t *temp_img) {
    if (!raw_data || !temp_img) return -1;
    if (rows <= 0 || cols <= 0) return -1;
    if (temp_img->rows != rows || temp_img->cols != cols) return -1;
    int N = rows * cols;
    for (int i = 0; i < N; i++)
        temp_img->data[i] = gain * raw_data[i] + offset;
    return 0;
}

ir_image_stats_t ir_image_roi_stats(const ir_thermal_image_t *img,
                                     int row_start, int col_start,
                                     int row_end, int col_end) {
    ir_image_stats_t stats = {0};
    if (!img || row_start < 0 || col_start < 0) return stats;
    if (row_end > img->rows || col_end > img->cols) return stats;
    if (row_start >= row_end || col_start >= col_end) return stats;

    double min_v = img->data[row_start * img->cols + col_start];
    double max_v = min_v, sum = 0.0, sum2 = 0.0;
    int count = 0;

    for (int r = row_start; r < row_end; r++) {
        for (int c = col_start; c < col_end; c++) {
            double v = img->data[r * img->cols + c];
            if (v < min_v) min_v = v;
            if (v > max_v) max_v = v;
            sum += v;
            sum2 += v * v;
            count++;
        }
    }
    stats.min = min_v;
    stats.max = max_v;
    stats.mean = sum / count;
    stats.stddev = sqrt(sum2 / count - stats.mean * stats.mean);
    stats.n_pixels = count;
    return stats;
}

int ir_image_difference(const ir_thermal_image_t *img1,
                         const ir_thermal_image_t *img2,
                         ir_thermal_image_t *diff) {
    if (!img1 || !img2 || !diff) return -1;
    if (img1->rows != img2->rows || img1->cols != img2->cols) return -1;
    if (diff->rows != img1->rows || diff->cols != img1->cols) return -1;
    int N = img1->rows * img1->cols;
    for (int i = 0; i < N; i++)
        diff->data[i] = img1->data[i] - img2->data[i];
    return 0;
}
