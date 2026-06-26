/**
 * @file    example_insar.c
 * @brief   SAR Interferometry Demo -- InSAR Phase and Coherence
 *
 * Demonstrates: interferogram generation, coherence estimation,
 * flat-earth phase removal, and phase-to-height conversion.
 * Simulates two SAR acquisitions with a baseline for InSAR.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sar_core.h"
#include "sar_geometry.h"
#include "sar_interferometry.h"

int main(void)
{
    printf("=== SAR Interferometry (InSAR) Demo ===\n\n");

    /* Create master and slave images */
    size_t nr = 128, nc = 128;
    sar_image_t *master = sar_image_alloc(nr, nc);
    sar_image_t *slave  = sar_image_alloc(nr, nc);

    /* Simulate master: constant magnitude + phase ramp (flat earth) */
    for (size_t r = 0; r < nr; r++) {
        for (size_t c = 0; c < nc; c++) {
            double phase = 0.02 * (double)c;
            master->data_I[r][c] = cos(phase);
            master->data_Q[r][c] = sin(phase);
        }
    }

    /* Simulate slave: master + topographic phase + noise
     * Topo phase: delta_phi = 4*pi*B_perp*h/(lambda*R0*sin(theta))
     * For a simple Gaussian hill at center */
    double B_perp = 200.0, lambda = 0.0566, R0 = 850e3, theta = 0.65;
    double H = 514e3;
double scale = 4.0 * M_PI * B_perp / (lambda * R0 * sin(theta));

    for (size_t r = 0; r < nr; r++) {
        for (size_t c = 0; c < nc; c++) {
            double flat_phase = 0.02 * (double)c;
            /* Gaussian hill: h = 500 * exp(-((r-64)^2+(c-64)^2)/400) meters */
            double dr = (double)((int64_t)r - 64), dc = (double)((int64_t)c - 64);
            double h = 500.0 * exp(-(dr*dr + dc*dc) / 400.0);
            double topo_phase = scale * h;
            double total_phase = flat_phase + topo_phase;
            slave->data_I[r][c] = cos(total_phase);
            slave->data_Q[r][c] = sin(total_phase);
        }
    }

    /* Compute interferogram */
    printf("Computing interferogram...\n");
    double **ifgram = (double **)malloc(nr * sizeof(double *));
    double *ifgram_flat = (double *)calloc(nr * nc, sizeof(double));
    for (size_t i = 0; i < nr; i++) ifgram[i] = ifgram_flat + i * nc;

    sar_interferogram_compute(master, slave, ifgram, nr, nc);

    printf("Interferogram statistics:\n");
    double minp = 1e9, maxp = -1e9, sump = 0;
    for (size_t i = 0; i < nr*nc; i++) {
        if (ifgram_flat[i] < minp) minp = ifgram_flat[i];
        if (ifgram_flat[i] > maxp) maxp = ifgram_flat[i];
        sump += ifgram_flat[i];
    }
    printf("  Phase range: [%.3f, %.3f] rad\n", minp, maxp);
    printf("  Mean phase:  %.3f rad\n", sump / (double)(nr*nc));

    /* Coherence estimation */
    printf("\nEstimating coherence...\n");
    sar_coherence_map_t *cm = sar_coherence_alloc(nr, nc);
    sar_coherence_estimate(master, slave, cm, 5, 5);
    printf("  Mean coherence: %.4f\n", cm->mean_coherence);

    /* Flat-earth removal */
    printf("\nRemoving flat-earth phase...\n");
    sar_flat_earth_removal(ifgram, nr, nc, 1.5, lambda, 514e3, theta);
    minp = 1e9; maxp = -1e9;
    for (size_t i = 0; i < nr*nc; i++) {
        if (ifgram_flat[i] < minp) minp = ifgram_flat[i];
        if (ifgram_flat[i] > maxp) maxp = ifgram_flat[i];
    }
    printf("  After flat-earth removal: [%.3f, %.3f] rad\n", minp, maxp);

    /* Phase unwrapping */
    printf("\nPhase unwrapping (Goldstein)...\n");
    double **unwrapped = (double **)malloc(nr * sizeof(double *));
    double *unwrapped_flat = (double *)calloc(nr * nc, sizeof(double));
    for (size_t i = 0; i < nr; i++) unwrapped[i] = unwrapped_flat + i * nc;

    sar_phase_unwrap_goldstein(ifgram, unwrapped, nr, nc);

    /* Phase to height */
    printf("Converting phase to height...\n");
    double **dem = (double **)malloc(nr * sizeof(double *));
    double *dem_flat = (double *)calloc(nr * nc, sizeof(double));
    for (size_t i = 0; i < nr; i++) dem[i] = dem_flat + i * nc;

    sar_phase_to_height((const double **)unwrapped, nr, nc,
                        B_perp, R0, lambda, H, theta, dem);

    /* Find maximum height */
    double max_h = 0; size_t max_hr = 0, max_hc = 0;
    for (size_t r = 0; r < nr; r++)
        for (size_t c = 0; c < nc; c++)
            if (dem[r][c] > max_h) { max_h = dem[r][c]; max_hr = r; max_hc = c; }

    printf("\nResults:\n");
    printf("  Maximum height: %.1f m at (%zu, %zu)\n", max_h, max_hr, max_hc);
    printf("  Expected peak:  ~500 m at (64, 64)\n");

    /* Cleanup */
    free(ifgram_flat); free(ifgram);
    free(unwrapped_flat); free(unwrapped);
    free(dem_flat); free(dem);
    sar_coherence_free(cm);
    sar_image_free(master);
    sar_image_free(slave);

    printf("\n=== InSAR Demo Complete ===\n");
    return 0;
}