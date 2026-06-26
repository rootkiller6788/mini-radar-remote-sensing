/* Linear Array Beam Steering Demo
 * Demonstrates beam steering with a 16-element uniform linear array at 10 GHz.
 * Computes array factor patterns for broadside and 30° scan.
 */

#include "phased_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(void)
{
    printf("=== Linear Phased Array Beam Steering Demo ===\n");
    printf("16-element ULA, d=lambda/2, f=10 GHz (X-band)\n\n");

    /* Configure array */
    pa_array_config_t config;
    pa_init_linear_array(&config, 16, 0.015, 10.0e9);

    uint32_t N = config.num_elements;
    pa_element_t *elements = pa_allocate_elements(&config);
    if (!elements) { printf("Error: allocation failed\n"); return 1; }

    /* Set all elements to half-wave dipole type */
    for (uint32_t i = 0; i < N; i++) {
        elements[i].element_type = PA_ELEMENT_HALF_WAVE_DIPOLE;
    }

    double lambda = pa_wavelength(config.frequency_hz);
    printf("Wavelength: %.3f mm\n", lambda * 1000.0);
    printf("Electrical spacing (d/lambda): %.3f\n",
           config.element_spacing_x / lambda);

    /* Demo 1: Broadside beam (theta = 90°) */
    printf("\n--- Broadside Beam (theta = 90 deg) ---\n");
    double complex *sv_broad = pa_steering_vector(&config, elements,
                                                   M_PI/2.0, 0.0);
    if (!sv_broad) { printf("Error: steering vector\n"); return 1; }

    printf("Steering vector (first 4 elements):\n");
    for (uint32_t i = 0; i < 4 && i < N; i++) {
        printf("  w[%u] = %.3f + j%.3f  (|w|=%.3f, arg=%.1f deg)\n",
               i, creal(sv_broad[i]), cimag(sv_broad[i]),
               cabs(sv_broad[i]), carg(sv_broad[i]) * 180.0/M_PI);
    }

    /* Compute array factor at broadside */
    pa_af_result_t r = pa_array_factor(&config, elements, sv_broad,
                                        M_PI/2.0, 0.0);
    printf("AF at broadside: |AF|=%.1f (%.1f dB)\n",
           r.af_magnitude, r.af_magnitude_db);

    /* Directivity */
    double D = pa_directivity(&config, elements, sv_broad);
    printf("Peak directivity: %.1f (%.1f dBi)\n", D, 10.0 * log10(D));

    /* Beamwidth */
    double hp_deg;
    pa_beamwidth_3db(&config, elements, sv_broad, M_PI/2.0, 0.0, &hp_deg);
    printf("3-dB beamwidth: %.2f deg\n", hp_deg);

    /* Demo 2: Steered beam to 30° */
    printf("\n--- Steered Beam (theta = 60 deg, 30 deg from broadside) ---\n");
    double theta_steer = M_PI / 2.0 - 30.0 * M_PI / 180.0;  /* 60° = 30° from broadside */

    double complex *sv_steer = pa_steering_vector(&config, elements,
                                                   theta_steer, 0.0);
    if (!sv_steer) { printf("Error: steered steering vector\n"); return 1; }

    printf("Phase shifts (first 8 elements):\n");
    double phase[16];
    pa_phase_shifts(&config, elements, theta_steer, 0.0, phase);
    for (uint32_t i = 0; i < 8 && i < N; i++) {
        printf("  elem[%2u]: phi = %+7.1f deg\n",
               i, phase[i] * 180.0 / M_PI);
    }

    /* Compute AF in steered direction */
    pa_af_result_t r_steer = pa_array_factor(&config, elements, sv_steer,
                                              theta_steer, 0.0);
    printf("AF at steer angle: |AF|=%.1f (%.1f dB)\n",
           r_steer.af_magnitude, r_steer.af_magnitude_db);

    /* Beamwidth broadens with scan */
    double hp_steer_deg;
    pa_beamwidth_3db(&config, elements, sv_steer, theta_steer, 0.0,
                     &hp_steer_deg);
    printf("3-dB beamwidth (scanned): %.2f deg\n", hp_steer_deg);

    /* Demo 3: Amplitude taper for sidelobe suppression */
    printf("\n--- Sidelobe Suppression with Taylor Taper (-30 dB) ---\n");
    double *taper = pa_taylor_weights(N, -30.0, 5);
    if (!taper) { printf("Error: taper computation\n"); return 1; }

    double complex *weights_tapered = (double complex *)malloc(
        N * sizeof(double complex));
    for (uint32_t i = 0; i < N; i++) {
        weights_tapered[i] = sv_broad[i] * taper[i];
    }

    printf("Amplitude taper (first 4 elements):\n");
    for (uint32_t i = 0; i < 4; i++) {
        printf("  taper[%u] = %.3f\n", i, taper[i]);
    }

    /* Compute pattern metrics with taper */
    pa_pattern_t *pat = pa_allocate_pattern(181, 1, 0.0, M_PI, 0.0, 0.0);
    if (!pat) { printf("Error: pattern allocation\n"); return 1; }

    pa_compute_pattern(&config, elements, weights_tapered, pat);
    printf("With Taylor (-30 dB): SLL = %.1f dB, HPBW = %.2f deg\n",
           pat->sidelobe_level_db, pat->half_power_beamwidth_deg);

    /* Gear ratio: SLL / HPBW trade-off */
    printf("\nTrade-off: Lower sidelobes → wider beamwidth\n");

    /* Demo 4: Grating lobe check */
    printf("\n--- Grating Lobe Analysis ---\n");
    pa_grating_lobe_t lobes[10];
    int n_lobes = pa_find_grating_lobes(&config, elements, weights_tapered,
                                         theta_steer, 0.0, 10, lobes);
    if (n_lobes > 0) {
        printf("Warning: %d grating lobe(s) detected!\n", n_lobes);
        for (int i = 0; i < n_lobes; i++) {
            printf("  Lobe %d: theta=%.1f deg, level=%.1f dB, order=m=%d\n",
                   i, lobes[i].theta_grating_rad * 180.0/M_PI,
                   lobes[i].relative_level_db, lobes[i].order_m);
        }
    } else {
        printf("No grating lobes in visible space (d < lambda).\n");
    }

    printf("\n=== Demo Complete ===\n");

    /* Cleanup */
    free(sv_broad);
    free(sv_steer);
    free(taper);
    free(weights_tapered);
    pa_free_pattern(pat);
    pa_free_elements(elements);

    return 0;
}
