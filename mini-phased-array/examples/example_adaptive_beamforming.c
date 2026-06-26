/* Adaptive Beamforming Demo — MVDR/Capon Interference Rejection
 *
 * 8-element ULA at 10 GHz. Desired signal from broadside, interferer
 * from 30° off broadside. Demonstrates:
 *   - Covariance matrix estimation from snapshots
 *   - MVDR weight computation (adaptive nulling)
 *   - Capon spatial spectrum (high-resolution DOA)
 *   - SINR improvement over conventional beamforming
 *
 * Scenario: Comms link with co-channel interference
 *   (5G mmWave beamforming / AESA radar jammer rejection)
 */

#include "phased_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(void)
{
    printf("=== Adaptive Beamforming (MVDR/Capon) Demo ===\n\n");
    printf("Scenario: 8-element ULA, f=10 GHz, d=lambda/2\n");
    printf("  Desired signal: broadside (theta=90 deg)\n");
    printf("  Interferer: 30 deg off broadside (theta=60 deg)\n");
    printf("  SNR = 20 dB, INR = 30 dB (strong jammer)\n\n");

    /* Array setup */
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    if (!elements) { printf("Error: allocation\n"); return 1; }

    uint32_t M = config.num_elements;
    uint32_t K = 200;  /* Number of snapshots */

    /* Generate synthetic snapshot data */
    double sigma_s = pow(10.0, 20.0/20.0);  /* Signal amplitude (SNR=20dB) */
    double sigma_i = pow(10.0, 30.0/20.0);  /* Interferer amplitude (INR=30dB) */
    double sigma_n = 1.0;                    /* Noise std */

    /* Steering vectors */
    double complex *a_signal   = pa_steering_vector(&config, elements,
                                                     M_PI/2.0, 0.0);
    double complex *a_interfer = pa_steering_vector(&config, elements,
                                                     M_PI/2.0 - 0.5236, 0.0);
    if (!a_signal || !a_interfer) { printf("Error: steering\n"); return 1; }

    /* Allocate snapshot array */
    double complex **snapshots = (double complex **)malloc(
        K * sizeof(double complex *));
    for (uint32_t k = 0; k < K; k++) {
        snapshots[k] = (double complex *)calloc(M, sizeof(double complex));
    }

    /* Generate snapshots: x[k] = s[k]*a_s + i[k]*a_i + noise[k] */
    srand(42);  /* Fixed seed for reproducibility */
    for (uint32_t k = 0; k < K; k++) {
        /* Random signal and interferer symbols (QPSK-like) */
        double s_phase = 2.0 * M_PI * (double)rand() / RAND_MAX;
        double i_phase = 2.0 * M_PI * (double)rand() / RAND_MAX;
        double complex s_sym = sigma_s * (cos(s_phase) + sin(s_phase) * I);
        double complex i_sym = sigma_i * (cos(i_phase) + sin(i_phase) * I);

        for (uint32_t m = 0; m < M; m++) {
            double nr = ((double)rand() / RAND_MAX - 0.5) * 2.0 * sigma_n;
            double ni = ((double)rand() / RAND_MAX - 0.5) * 2.0 * sigma_n;
            snapshots[k][m] = s_sym * a_signal[m]
                            + i_sym * a_interfer[m]
                            + (nr + ni * I);
        }
    }

    /* Conventional beamformer (delay-and-sum) */
    printf("--- Conventional Beamformer (uniform weights) ---\n");
    double complex *w_conv = (double complex *)malloc(M * sizeof(double complex));
    for (uint32_t m = 0; m < M; m++) w_conv[m] = 1.0 / sqrt((double)M);

    /* Compute conventional output power */
    double conv_output_power = 0.0;
    for (uint32_t k = 0; k < K; k++) {
        double complex y = 0.0 + 0.0 * I;
        for (uint32_t m = 0; m < M; m++) {
            y += conj(w_conv[m]) * snapshots[k][m];
        }
        conv_output_power += cabs(y) * cabs(y);
    }
    conv_output_power /= (double)K;
    printf("Conventional output power: %.2f\n", conv_output_power);

    /* MVDR adaptive beamformer */
    printf("\n--- MVDR Adaptive Beamformer ---\n");
    pa_mvdr_beamformer_t *mvdr = pa_mvdr_beamformer_init(M, K);
    if (!mvdr) { printf("Error: MVDR init\n"); return 1; }

    /* Estimate covariance */
    pa_smi_estimate_covariance(&mvdr->base, (const double complex **)snapshots);

    /* Compute MVDR weights for broadside look direction */
    pa_mvdr_compute_weights(mvdr, &config, elements, M_PI/2.0, 0.0);

    printf("MVDR weights (first 4 elements):\n");
    for (uint32_t m = 0; m < 4; m++) {
        printf("  w[%u] = %.3f + j%.3f  (|w|=%.3f, arg=%.1f deg)\n",
               m, creal(mvdr->base.weight_vector[m]),
               cimag(mvdr->base.weight_vector[m]),
               cabs(mvdr->base.weight_vector[m]),
               carg(mvdr->base.weight_vector[m]) * 180.0 / M_PI);
    }

    /* Compute MVDR output power */
    double mvdr_output_power = 0.0;
    for (uint32_t k = 0; k < K; k++) {
        double complex y = 0.0 + 0.0 * I;
        for (uint32_t m = 0; m < M; m++) {
            y += conj(mvdr->base.weight_vector[m]) * snapshots[k][m];
        }
        mvdr_output_power += cabs(y) * cabs(y);
    }
    mvdr_output_power /= (double)K;
    printf("MVDR output power: %.2f\n", mvdr_output_power);

    /* Output SINR */
    double sinr = pa_mvdr_output_sinr(mvdr);
    printf("MVDR output SINR: %.1f dB\n", sinr);

    /* Capon spatial spectrum */
    printf("\n--- Capon Spatial Spectrum ---\n");
    uint32_t N_angles = 91;  /* 0° to 90° from broadside */
    double *theta_vals = (double *)malloc(N_angles * sizeof(double));
    double *phi_vals   = (double *)malloc(N_angles * sizeof(double));
    double *spectrum   = (double *)malloc(N_angles * sizeof(double));

    for (uint32_t i = 0; i < N_angles; i++) {
        theta_vals[i] = M_PI/2.0 - (double)i * (M_PI/180.0);
        phi_vals[i]   = 0.0;
    }

    pa_mvdr_capon_spectrum(mvdr, &config, elements,
                            spectrum, theta_vals, phi_vals, N_angles);

    printf("Angle(deg)  |  P_Capon\n");
    printf("------------+-----------\n");
    /* Show selected angles: broadside, interferer, and nearby */
    int angles_to_show[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    for (int a = 0; a < 10; a++) {
        int ia = angles_to_show[a];
        if (ia < (int)N_angles) {
            printf("  %3d       |  %.2e\n", ia, spectrum[ia]);
        }
    }

    /* Find peak in spectrum (should be near 0°, the signal DOA) */
    double max_spectrum = spectrum[0];
    int max_idx = 0;
    for (uint32_t i = 1; i < N_angles; i++) {
        if (spectrum[i] > max_spectrum) {
            max_spectrum = spectrum[i];
            max_idx = (int)i;
        }
    }
    printf("Capon spectrum peak at %.0f deg: %.2e\n",
           (double)max_idx, max_spectrum);

    /* Performance summary */
    printf("\n--- Performance Summary ---\n");
    printf("Conventional beamformer: P_out = %.2f (signal+interferer+noise)\n",
           conv_output_power);
    printf("MVDR (adaptive): P_out = %.2f (nulls interferer)\n",
           mvdr_output_power);
    double improvement = 10.0 * log10(conv_output_power / mvdr_output_power);
    printf("Interference rejection: %.1f dB (ratio conv/mvdr power)\n",
           improvement);
    printf("\nMVDR places a spatial null at the interferer direction\n");
    printf("while maintaining unity gain in the look direction.\n");

    printf("\n=== Adaptive Beamforming Demo Complete ===\n");

    /* Cleanup */
    for (uint32_t k = 0; k < K; k++) free(snapshots[k]);
    free(snapshots);
    free(w_conv);
    free(a_signal);
    free(a_interfer);
    free(theta_vals);
    free(phi_vals);
    free(spectrum);
    pa_mvdr_beamformer_free(mvdr);
    pa_free_elements(elements);

    return 0;
}
