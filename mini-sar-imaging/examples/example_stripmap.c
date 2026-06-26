/**
 * @file    example_stripmap.c
 * @brief   Stripmap SAR Processing Demo -- RDA simplified
 *
 * Demonstrates the stripmap SAR mode processing chain:
 * raw data simulation -> range compression -> simple azimuth FFT focusing.
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
    printf("=== Stripmap SAR Processing Demo ===\n\n");

    /* C-band SAR parameters (like Radarsat-2) */
    double f0 = 5.405e9, B = 100e6, tau = 42e-6, PRF = 2000.0;
    double L_a = 7.5, v = 7500.0, H = 798e3;
    double lambda = SAR_C / f0;
    double rho_r = SAR_C / (2.0 * B);
    double rho_a = L_a / 2.0;

    printf("C-Band Stripmap SAR:\n");
    printf("  Range resolution:  %.2f m\n", rho_r);
    printf("  Azimuth resolution: %.2f m (focused)\n", rho_a);
    printf("  PRF:               %.0f Hz\n", PRF);
    printf("  Swath width:       ~%.0f km\n", 50.0);

    /* Create small raw dataset */
    size_t naz = 256, nrng = 512;
    sar_raw_data_t *raw = sar_raw_data_alloc(naz, nrng);
    if (!raw) return 1;

    sar_params_t sp;
    sar_params_init(&sp, f0, B, tau, PRF, L_a, v, H, 0.0, 0.35, 1020e3, 1070e3);
    raw->params = sp;
    raw->range_sampling_interval = SAR_C / (2.0 * 1.2 * B);
    raw->azimuth_sampling_interval = 1.0 / PRF;

    /* Simulate 3 point targets at different ranges */
    printf("\nSimulating 3 point targets...\n");
    sar_raw_data_point_target(raw, nrng/4, naz/3, 1.0);
    sar_raw_data_point_target(raw, nrng/2, naz/2, 0.8);
    sar_raw_data_point_target(raw, 3*nrng/4, 2*naz/3, 0.6);

    /* Add noise */
    printf("Adding thermal noise...\n");
    sar_raw_data_add_noise(raw, 0.05);

    /* Range compression */
    printf("Range compression via FFT...\n");
    size_t Nc = 512;
    sar_chirp_t *chirp = sar_chirp_alloc(Nc, tau, B, 1.2*B);
    double *h_I = malloc(Nc*8), *h_Q = malloc(Nc*8);
    sar_matched_filter_coeff(chirp, h_I, h_Q);

    sar_image_t *image = sar_image_alloc(naz, nrng);
    for (size_t a = 0; a < naz; a++) {
        double *tmpI = malloc(nrng*8), *tmpQ = malloc(nrng*8);
        memcpy(tmpI, raw->data_I[a], nrng*8);
        memcpy(tmpQ, raw->data_Q[a], nrng*8);
        sar_fft_range_1d(tmpI, tmpQ, nrng, 1);
        /* Multiply by conj(H) in freq domain */
        double *HI = malloc(nrng*8), *HQ = malloc(nrng*8);
        memset(HI,0,nrng*8); memset(HQ,0,nrng*8);
        memcpy(HI,h_I,Nc*8); memcpy(HQ,h_Q,Nc*8);
        sar_fft_range_1d(HI, HQ, nrng, 1);
        for (size_t r=0; r<nrng; r++) {
            double pI=tmpI[r]*HI[r]+tmpQ[r]*HQ[r];
            double pQ=tmpQ[r]*HI[r]-tmpI[r]*HQ[r];
            tmpI[r]=pI; tmpQ[r]=pQ;
        }
        sar_fft_range_1d(tmpI, tmpQ, nrng, 0);
        memcpy(image->data_I[a], tmpI, nrng*8);
        memcpy(image->data_Q[a], tmpQ, nrng*8);
        free(tmpI); free(tmpQ); free(HI); free(HQ);
    }
    image->range_pixel_spacing_m = rho_r;
    image->azimuth_pixel_spacing_m = rho_a;

    /* RDA simplified: azimuth FFT for each range */
    printf("Azimuth focusing (FFT-based)...\n");
    for (size_t c = 0; c < nrng; c++) {
        double *colI = malloc(naz*8), *colQ = malloc(naz*8);
        for (size_t r=0; r<naz; r++) { colI[r]=image->data_I[r][c]; colQ[r]=image->data_Q[r][c]; }
        sar_fft_range_1d(colI, colQ, naz, 1);
        for (size_t r=0; r<naz; r++) { image->data_I[r][c]=colI[r]; image->data_Q[r][c]=colQ[r]; }
        free(colI); free(colQ);
    }

    /* Find peaks */
    double max_mag = 0; size_t max_r=0, max_c=0;
    for (size_t r=0; r<naz; r++)
        for (size_t c=0; c<nrng; c++) {
            double mag = image->data_I[r][c]*image->data_I[r][c]
                       + image->data_Q[r][c]*image->data_Q[r][c];
            if (mag > max_mag) { max_mag=mag; max_r=r; max_c=c; }
        }

    printf("\nStrongest target detected at (az=%zu, rng=%zu), mag=%.2f\n",
           max_r, max_c, sqrt(max_mag));

    /* Cleanup */
    free(h_I); free(h_Q);
    sar_chirp_free(chirp);
    sar_image_free(image);
    sar_raw_data_free(raw);

    printf("\n=== Stripmap Demo Complete ===\n");
    return 0;
}