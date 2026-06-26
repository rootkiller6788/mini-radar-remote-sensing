/**
 * @file    ir_enhancement.c
 * @brief   Thermal Image Enhancement - DDE, AGC, Sharpening
 *
 * Implements Digital Detail Enhancement (DDE) and other visualization
 * algorithms for thermal imagery. DDE separates the image into a low-frequency
 * (contrast) component and a high-frequency (detail) component, enhancing
 * each independently for optimal display on 8-bit monitors.
 *
 * References:
 *   FLIR Systems (2015) "DDE - Digital Detail Enhancement"
 *   Gonzalez & Woods (2018) "Digital Image Processing", 4th Ed.
 *   Peli, T. (1990) "Multiscale image contrast enhancement", Opt. Eng.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ========================================================================
 * Digital Detail Enhancement (DDE)
 *
 * DDE separates the thermal image into:
 *   1. Base layer (low-pass filtered): contains large-scale contrast
 *   2. Detail layer (original - base): contains fine details
 *
 * The two layers are processed independently and recombined:
 *   output = gain_base * base + gain_detail * detail
 *
 * This allows simultaneous display of both hot/cold objects and fine
 * texture details in the same 8-bit image.
 * ======================================================================== */

int ir_dde_process(const double *input, int rows, int cols,
                    double sigma_filter, double gain_detail,
                    double *base_out, double *detail_out,
                    double *enhanced_out) {
    if (!input || rows <= 0 || cols <= 0 || sigma_filter <= 0.0)
        return -1;

    /* int N = rows * cols; */ (void)rows; (void)cols;
    int ks = (int)(6.0 * sigma_filter + 1.0);
    if (ks % 2 == 0) ks++;
    int half = ks / 2;

    /* Build Gaussian kernel */
    double *kernel = (double*)malloc(ks * ks * sizeof(double));
    double ksum = 0.0;
    for (int i = 0; i < ks; i++) {
        for (int j = 0; j < ks; j++) {
            double dx = i - half, dy = j - half;
            double v = exp(-(dx*dx + dy*dy) / (2.0 * sigma_filter * sigma_filter));
            kernel[i * ks + j] = v;
            ksum += v;
        }
    }
    for (int i = 0; i < ks * ks; i++) kernel[i] /= ksum;

    /* Low-pass filter to get base layer */
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double v = 0.0;
            for (int dr = -half; dr <= half; dr++) {
                for (int dc = -half; dc <= half; dc++) {
                    int nr = r + dr, nc = c + dc;
                    if (nr >= 0 && nr < rows && nc >= 0 && nc < cols)
                        v += input[nr * cols + nc]
                             * kernel[(dr + half) * ks + (dc + half)];
                }
            }
            int idx = r * cols + c;
            if (base_out) base_out[idx] = v;
            if (detail_out) detail_out[idx] = input[idx] - v;
            if (enhanced_out) {
                enhanced_out[idx] = v + gain_detail * (input[idx] - v);
            }
        }
    }

    free(kernel);
    return 0;
}

/* ========================================================================
 * Contrast Enhancement
 * ======================================================================== */

/** Linear contrast stretch: maps [in_min, in_max] to [out_min, out_max] */
int ir_contrast_stretch(const double *input, int length,
                         double in_min, double in_max,
                         double out_min, double out_max,
                         double *output) {
    if (!input || !output || length <= 0) return -1;
    if (fabs(in_max - in_min) < 1e-10) return -1;

    double scale = (out_max - out_min) / (in_max - in_min);
    for (int i = 0; i < length; i++) {
        double v = input[i];
        if (v < in_min) v = in_min;
        if (v > in_max) v = in_max;
        output[i] = out_min + scale * (v - in_min);
    }
    return 0;
}

/** Gamma correction for display optimization.
 *  output = input^gamma  (normalized to [0,1]) */
int ir_gamma_correction(const double *input, int length,
                         double gamma, double *output) {
    if (!input || !output || length <= 0 || gamma <= 0.0) return -1;
    for (int i = 0; i < length; i++) {
        double v = input[i];
        if (v < 0.0) v = 0.0;
        if (v > 1.0) v = 1.0;
        output[i] = pow(v, 1.0 / gamma);
    }
    return 0;
}

/* ========================================================================
 * Temporal Noise Reduction (simple frame averaging)
 * ======================================================================== */

int ir_temporal_average(const double *frame, double *accumulator,
                         int length, int frame_count, double *output) {
    if (!frame || !accumulator || !output || length <= 0 || frame_count < 1)
        return -1;

    double alpha = 1.0 / frame_count;
    for (int i = 0; i < length; i++) {
        accumulator[i] += frame[i];
        output[i] = accumulator[i] * alpha;
    }
    return 0;
}

/* ========================================================================
 * Super-Resolution by Pixel Interpolation (Bilinear)
 * ======================================================================== */

int ir_bilinear_upscale(const double *input, int in_rows, int in_cols,
                         int out_rows, int out_cols, double *output) {
    if (!input || !output || in_rows <= 0 || in_cols <= 0) return -1;
    if (out_rows <= 0 || out_cols <= 0) return -1;

    double scale_r = (double)in_rows / out_rows;
    double scale_c = (double)in_cols / out_cols;

    for (int r = 0; r < out_rows; r++) {
        for (int c = 0; c < out_cols; c++) {
            double src_r = r * scale_r;
            double src_c = c * scale_c;
            int r0 = (int)src_r;
            int c0 = (int)src_c;
            int r1 = (r0 + 1 < in_rows) ? r0 + 1 : r0;
            int c1 = (c0 + 1 < in_cols) ? c0 + 1 : c0;
            double dr = src_r - r0;
            double dc = src_c - c0;

            double v00 = input[r0 * in_cols + c0];
            double v01 = input[r0 * in_cols + c1];
            double v10 = input[r1 * in_cols + c0];
            double v11 = input[r1 * in_cols + c1];

            output[r * out_cols + c] = (1.0 - dr) * (1.0 - dc) * v00
                                       + (1.0 - dr) * dc * v01
                                       + dr * (1.0 - dc) * v10
                                       + dr * dc * v11;
        }
    }
    return 0;
}
