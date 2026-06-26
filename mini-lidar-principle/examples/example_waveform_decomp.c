/**
 * @file    example_waveform_decomp.c
 * @brief   End-to-end example: Full-waveform LiDAR decomposition
 *
 * Simulates a digitized LiDAR waveform with multiple overlapping returns
 * from forest canopy and ground, then decomposes using Gaussian fitting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lidar_core.h"
#include "lidar_waveform.h"

int main(void) {
    printf("=== LiDAR Full-Waveform Decomposition Example ===\n\n");

    /* Simulate waveform: 5 ns FWHM pulse, sampled at 0.25 ns intervals.
       Returns at: 30 ns (canopy top), 34 ns (mid-canopy), 42 ns (ground) */
    lidar_waveform_t wf;
    lidar_waveform_init(&wf, 300, 0.25, -5.0, 5.0);

    /* Define ground-truth components */
    lidar_gaussian_component_t comps[3];
    memset(comps, 0, sizeof(comps));
    /* Canopy top: early return, moderate amplitude */
    comps[0].amplitude = 0.7;
    comps[0].center    = 30.0;  /* ns */
    comps[0].sigma     = 2.5;
    /* Mid-canopy: partial penetration */
    comps[1].amplitude = 0.4;
    comps[1].center    = 34.0;
    comps[1].sigma     = 2.5;
    /* Ground: strongest return, later */
    comps[2].amplitude = 1.0;
    comps[2].center    = 42.0;
    comps[2].sigma     = 1.5;

    /* Synthesize waveform with noise */
    lidar_waveform_synthesize(&wf, comps, 3, 0.02);
    lidar_waveform_noise_estimate(&wf, 0.2);

    printf("Synthesized waveform: %zu samples, dt=%.2f ns\n",
           wf.num_samples, wf.dt);
    printf("Ground truth: 3 returns at %.1f, %.1f, %.1f ns\n",
           comps[0].center, comps[1].center, comps[2].center);
    printf("Max amplitude: %.3f at sample %zu\n",
           wf.max_amplitude, wf.max_index);

    /* Peak detection */
    double peak_times[10];
    int n_peaks = lidar_detect_peaks_derivative(&wf, peak_times, 10, 5.0, 2);
    printf("\nDerivative peak detection: %d peaks found\n", n_peaks);
    for (int i = 0; i < n_peaks; i++) {
        printf("  Peak %d at t = %.2f ns\n", i + 1, peak_times[i]);
    }

    /* CFD timing for strongest peak */
    double cfd_t = lidar_cfd_timing(&wf, 0.5);
    printf("CFD timing (50%%): %.2f ns\n", cfd_t);

    /* Leading-edge timing */
    double le_t = lidar_leading_edge_timing(&wf, 3.0);
    printf("Leading-edge timing (3-sigma): %.2f ns\n", le_t);

    /* FWHM of main peak */
    double fwhm = lidar_pulse_fwhm(&wf, wf.max_index);
    printf("Main peak FWHM: %.2f ns\n", fwhm);

    /* Energy of main peak */
    double energy = lidar_pulse_energy(&wf, wf.max_index, 30);
    printf("Main peak energy: %.4f\n", energy);

    /* Gaussian decomposition */
    printf("\nPerforming Gaussian decomposition...\n");
    lidar_waveform_decomp_t result;
    if (lidar_gaussian_decompose(&wf, &result, 100, 1e-6) == 0) {
        printf("Decomposition result:\n");
        printf("  Components found: %zu\n", result.num_components);
        printf("  Converged: %s\n", result.converged ? "yes" : "no");
        printf("  Iterations: %d\n", result.num_iterations);
        printf("  Residual RMS: %.4f\n", result.residual_rms);
        printf("  R-squared: %.4f\n\n", result.r_squared);

        for (size_t j = 0; j < result.num_components; j++) {
            double range = result.components[j].center * LIDAR_C / 2000.0;
            printf("  Component %zu:\n", j + 1);
            printf("    Amplitude: %.3f\n", result.components[j].amplitude);
            printf("    Center: %.2f ns\n", result.components[j].center);
            printf("    FWHM: %.2f ns\n", result.components[j].fwhm);
            printf("    Range: %.2f m\n", range);
            printf("    Energy: %.4f\n", result.components[j].energy);
        }

        if (result.r_squared > 0.9) {
            printf("\nExcellent fit (R^2 > 0.9) — decomposition successful!\n");
        } else if (result.r_squared > 0.7) {
            printf("\nAcceptable fit (R^2 > 0.7).\n");
        } else {
            printf("\nPoor fit — try more iterations or adjust initial estimates.\n");
        }
    } else {
        printf("Decomposition failed.\n");
    }

    printf("\nWaveform decomposition example completed.\n");
    lidar_waveform_free(&wf);
    return 0;
}