/* Planar Array 2D Beam Steering Demo
 * 8x8 rectangular grid array at 10 GHz, steering in azimuth and elevation.
 * Demonstrates 2D pattern computation and AESA radar performance.
 */

#include "phased_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(void)
{
    printf("=== Planar Phased Array (8x8 AESA) Demo ===\n\n");

    /* 8x8 planar array, λ/2 spacing at 10 GHz */
    pa_array_config_t config;
    pa_init_planar_array(&config, 8, 8, 0.015, 0.015, 10.0e9);

    uint32_t R = config.rows;
    uint32_t C = config.cols;
    uint32_t N = config.num_elements;

    printf("Array: %ux%u = %u elements\n", R, C, N);
    printf("Spacing: dx=%.2f mm, dy=%.2f mm\n",
           config.element_spacing_x * 1000, config.element_spacing_y * 1000);
    printf("Frequency: %.1f GHz (X-band)\n", config.frequency_hz / 1e9);

    pa_element_t *elements = pa_allocate_elements(&config);
    if (!elements) { printf("Error: allocation\n"); return 1; }

    /* Assign realistic T/R module parameters */
    for (uint32_t i = 0; i < N; i++) {
        elements[i].tr_module_id       = (int)i;
        elements[i].tr_module_power_watt = 10.0;
        elements[i].tr_module_gain_db   = 20.0;
        elements[i].tr_module_noise_figure_db = 3.0;
    }

    /* Uniform weights (no taper) — compute pattern */
    double complex *weights = (double complex *)malloc(N * sizeof(double complex));
    for (uint32_t i = 0; i < N; i++) weights[i] = 1.0 + 0.0 * I;

    /* Scan in azimuth cut (phi=0, vary theta) */
    uint32_t N_theta = 181;
    pa_pattern_t *pat = pa_allocate_pattern(N_theta, 1, 0.0, M_PI, 0.0, 0.0);
    if (!pat) { printf("Error: pattern alloc\n"); return 1; }

    pa_compute_pattern(&config, elements, weights, pat);

    printf("\n--- Pattern Metrics (uniform, broadside) ---\n");
    printf("Half-power beamwidth: %.2f deg\n", pat->half_power_beamwidth_deg);
    printf("First-null beamwidth: %.2f deg\n", pat->first_null_beamwidth_deg);
    printf("Peak sidelobe level: %.1f dB\n", pat->sidelobe_level_db);
    printf("Peak directivity: %.1f dBi\n", pat->max_directivity_dbi);

    /* Print pattern cut at selected angles */
    printf("\nPattern cut (phi=0, theta sweep):\n");
    printf("  theta(deg)  |  AF(dB)\n");
    printf("  ------------+---------\n");
    for (int i = 0; i < 10; i++) {
        int idx = 90 + i * 5;  /* Near broadside (theta ~= 90 deg) */
        if (idx < (int)N_theta) {
            double theta_deg = (double)idx * 180.0 / (double)(N_theta - 1);
            printf("  %6.1f      |  %+6.1f\n", theta_deg, pat->af_db[idx]);
        }
    }

    /* Radar range equation for this AESA */
    printf("\n--- Radar Performance (AESA) ---\n");
    double range_1m2 = pa_radar_range_equation(&config, elements, 1.0, 13.0, 5.0);
    double range_01m2 = pa_radar_range_equation(&config, elements, 0.1, 13.0, 5.0);

    printf("Detection range (sigma=1.0 m^2, SNR=13 dB): %.1f km\n",
           range_1m2 / 1000.0);
    printf("Detection range (sigma=0.1 m^2, SNR=13 dB): %.1f km\n",
           range_01m2 / 1000.0);

    /* Effective aperture */
    double Ae = pa_effective_aperture(&config, elements, 0.7);
    printf("Effective aperture (eta=0.7): %.4f m^2\n", Ae);
    printf("Physical aperture: %.4f m x %.4f m = %.4f m^2\n",
           (double)C * config.element_spacing_x,
           (double)R * config.element_spacing_y,
           ((double)C * config.element_spacing_x) * ((double)R * config.element_spacing_y));

    /* G/T */
    double gt = pa_array_gt_fom(&config, elements, 290.0);
    printf("Array G/T: %.1f dB/K\n", gt);

    /* Steered beam at Az=30°, El=10° */
    printf("\n--- Steered Beam (Az=30 deg, El=10 deg) ---\n");
    double theta_s = M_PI/2.0 - 10.0 * M_PI/180.0;  /* El=10 -> theta=80° */
    double phi_s   = 30.0 * M_PI / 180.0;

    double complex *sv_steer = pa_steering_vector(&config, elements,
                                                   theta_s, phi_s);
    if (!sv_steer) { printf("Error: steering vector\n"); return 1; }

    /* Compute pattern with steered weights */
    pa_pattern_t *pat_steer = pa_allocate_pattern(181, 1,
                                                    theta_s - 0.3, theta_s + 0.3,
                                                    phi_s, phi_s);
    if (!pat_steer) { printf("Error: pattern alloc 2\n"); return 1; }

    pa_compute_pattern(&config, elements, sv_steer, pat_steer);
    printf("Steered beam HPBW: %.2f deg\n",
           pat_steer->half_power_beamwidth_deg);

    /* With Dolph-Chebyshev taper */
    double *cheb = pa_dolph_chebyshev_weights(C, -40.0);
    if (cheb) {
        printf("\n--- With Dolph-Chebyshev Taper (SLL=-40 dB, cols) ---\n");
        printf("Chebyshev taper (columns):\n");
        for (uint32_t c = 0; c < C; c++) {
            printf("  [%u] = %.3f\n", c, cheb[c]);
        }
        free(cheb);
    }

    /* T/R module health check */
    printf("\n--- T/R Module Health ---\n");
    int healthy = 0, total = 0;
    for (uint32_t i = 0; i < N; i++) {
        if (elements[i].tr_module_id >= 0) {
            total++;
            pa_tr_module_t mod = pa_tr_module_init(elements[i].tr_module_id);
            if (mod.is_healthy) healthy++;
        }
    }
    printf("T/R modules: %d/%d healthy\n", healthy, total);

    printf("\n=== Planar Array Demo Complete ===\n");

    /* Cleanup */
    free(weights);
    free(sv_steer);
    pa_free_pattern(pat);
    pa_free_pattern(pat_steer);
    pa_free_elements(elements);

    return 0;
}
