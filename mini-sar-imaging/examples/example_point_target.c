/**
 * @file    example_point_target.c
 * @brief   SAR Point Target Simulation and Focusing Demo
 *
 * Demonstrates: chirp generation, raw data simulation for a point target,
 * range compression (matched filtering), and impulse response analysis.
 * This is the canonical L6 problem: simulate and focus a single point.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sar_core.h"
#include "sar_geometry.h"
#include "sar_algorithm.h"

void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward);

int main(void)
{
    printf("=== SAR Point Target Simulation & Focusing ===\n\n");

    /* Radar parameters: X-band spaceborne SAR */
    double f0 = 9.65e9;       /* Carrier frequency [Hz] */
    double B  = 150e6;        /* Range bandwidth [Hz] */
    double tau = 10e-6;       /* Pulse width [s] */
    double PRF = 3000.0;      /* PRF [Hz] */
    double L_a = 4.8;         /* Antenna length [m] */
    double v   = 7500.0;      /* Platform velocity [m/s] */
    double H   = 514e3;       /* Altitude [m] */
    double R0  = 850e3;       /* Slant range [m] */

    /* Derived parameters */
    double lambda = SAR_C / f0;
    double rho_r = SAR_C / (2.0 * B);
    double rho_a = L_a / 2.0;
    double T_sa = lambda * R0 / (L_a * v);

    printf("System Parameters:\n");
    printf("  Frequency:      %.2f GHz\n", f0/1e9);
    printf("  Wavelength:     %.3f m\n", lambda);
    printf("  Bandwidth:      %.0f MHz\n", B/1e6);
    printf("  Range res:      %.2f m\n", rho_r);
    printf("  Azimuth res:    %.2f m\n", rho_a);
    printf("  PRF:            %.0f Hz\n", PRF);
    printf("  Synth aperture: %.3f s\n", T_sa);

    /* Generate chirp */
    size_t N_chirp = 1024;
    double fs = 1.2 * B;
    sar_chirp_t *chirp = sar_chirp_alloc(N_chirp, tau, B, fs);
    if (!chirp) { printf("ERROR: chirp allocation failed\n"); return 1; }

    printf("\nChirp:\n");
    printf("  Samples:        %zu\n", chirp->num_samples);
    printf("  TBP:            %.0f\n", chirp->time_bandwidth_product);
    printf("  Chirp rate:     %.2e Hz/s\n", chirp->chirp_rate);

    /* Create raw data with single point target */
    size_t naz = 256, nrng = 512;
    sar_raw_data_t *raw = sar_raw_data_alloc(naz, nrng);
    if (!raw) { printf("ERROR: raw data alloc failed\n"); sar_chirp_free(chirp); return 1; }

    /* Initialize SAR parameters */
    sar_params_t sp;
    sar_params_init(&sp, f0, B, tau, PRF, L_a, v, H, 0.0, 0.35, R0-5000, R0+5000);
    raw->params = sp;
    raw->range_sampling_interval = SAR_C / (2.0 * fs);
    raw->azimuth_sampling_interval = 1.0 / PRF;

    /* Place point target at center */
    sar_raw_data_point_target(raw, nrng/2, naz/2, 1.0);

    printf("\nRaw data:\n");
    printf("  Azimuth lines:  %zu\n", raw->naz);
    printf("  Range bins:     %zu\n", raw->nrng);

    /* Range compression */
    printf("\nApplying range compression...\n");
    double *h_I = (double *)malloc(N_chirp * sizeof(double));
    double *h_Q = (double *)malloc(N_chirp * sizeof(double));
    sar_matched_filter_coeff(chirp, h_I, h_Q);

    double **rc_I = (double **)malloc(naz * sizeof(double *));
    double **rc_Q = (double **)malloc(naz * sizeof(double *));
    for (size_t i = 0; i < naz; i++) {
        rc_I[i] = (double *)malloc(nrng * sizeof(double));
        rc_Q[i] = (double *)malloc(nrng * sizeof(double));
    }

    /* Simple range compression along each azimuth line */
    for (size_t a = 0; a < naz; a++) {
        sar_pulse_compression_fft(raw->data_I[a], raw->data_Q[a], nrng,
                                   h_I, h_Q, N_chirp, rc_I[a], rc_Q[a]);
    }

    printf("Range compression complete.\n");

    /* Create image from range-compressed data */
    sar_image_t *image = sar_image_alloc(naz, nrng);
    for (size_t a = 0; a < naz; a++) {
        memcpy(image->data_I[a], rc_I[a], nrng * sizeof(double));
        memcpy(image->data_Q[a], rc_Q[a], nrng * sizeof(double));
    }
    image->range_pixel_spacing_m = rho_r;
    image->azimuth_pixel_spacing_m = rho_a;

    /* Analyze impulse response */
    sar_impulse_response_t ir;
    sar_analyze_impulse_response(image, naz/2, nrng/2, 30, &ir);

    printf("\nImpulse Response Analysis:\n");
    printf("  Peak location:    (%.0f, %.0f)\n", ir.peak_azimuth_idx, ir.peak_range_idx);
    printf("  Peak magnitude:   %.4f\n", ir.peak_magnitude);
    printf("  Range -3dB width: %.3f m (expected: %.3f m)\n",
           ir.range_resolution_3dB_m, rho_r);
    printf("  Azim -3dB width:  %.3f m (expected: %.3f m)\n",
           ir.azimuth_resolution_3dB_m, rho_a);
    printf("  Range PSLR:       %.1f dB\n", ir.range_pslr_db);
    printf("  Azimuth PSLR:     %.1f dB\n", ir.azimuth_pslr_db);

    /* Cleanup */
    for (size_t i = 0; i < naz; i++) { free(rc_I[i]); free(rc_Q[i]); }
    free(rc_I); free(rc_Q);
    free(h_I); free(h_Q);
    sar_image_free(image);
    sar_raw_data_free(raw);
    sar_chirp_free(chirp);

    printf("\n=== Demo Complete ===\n");
    return 0;
}