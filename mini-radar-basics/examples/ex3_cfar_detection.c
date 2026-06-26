/**
 * @file    ex3_cfar_detection.c
 * @brief   Example: CFAR target detection in noise with ROC analysis
 *
 * Demonstrates: CA-CFAR, OS-CFAR, GO-CFAR detection on noisy range profile,
 *               ROC curve generation, and Pd/Pfa trade-off analysis.
 *
 * Knowledge: L5 CFAR, L6 Target detection, L4 Detection theory
 * Reference: Richards, Scheer & Holm (2010), Ch.16
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "radar_detection.h"

int main(void)
{
    printf("=== Example 3: CFAR Detection and ROC Analysis ===\n\n");

    /* Generate synthetic range profile with noise + target */
    size_t nbins = 256;
    double *rp = calloc(nbins, sizeof(double));

    /* Noise floor: unit-mean exponential (square-law detector output) */
    srand(12345);
    for (size_t i = 0; i < nbins; i++) {
        double u = (double)rand() / RAND_MAX;
        if (u < 1e-12) u = 1e-12;
        rp[i] = -log(u);  /* Exponential with mean = 1.0 */
    }

    /* Inject targets */
    rp[50] = 25.0;   /* Strong target at bin 50 */
    rp[100] = 12.0;  /* Medium target at bin 100 */
    rp[150] = 8.0;   /* Weak target at bin 150 */

    printf("Generated range profile: %zu bins\n", nbins);
    printf("  Noise floor mean = 1.0 (exponential)\n");
    printf("  Targets at bins 50 (SNR~14dB), 100 (SNR~11dB), 150 (SNR~9dB)\n\n");

    /* CA-CFAR detection */
    cfar_config_t cfg_ca;
    cfar_config_init(&cfg_ca, CFAR_CA, 4, 16, 1e-4);
    int *det_ca = calloc(nbins, sizeof(int));
    cfar_detect(rp, nbins, &cfg_ca, det_ca);

    printf("CA-CFAR Results (Pfa=1e-4, guard=4, ref=16):\n");
    for (size_t i = 0; i < nbins; i++) {
        if (det_ca[i]) printf("  Detection at bin %zu: power=%.1f\n", i, rp[i]);
    }

    /* OS-CFAR detection */
    cfar_config_t cfg_os;
    cfar_config_init(&cfg_os, CFAR_OS, 4, 16, 1e-4);
    int *det_os = calloc(nbins, sizeof(int));
    cfar_os_detect(rp, nbins, &cfg_os, det_os);

    printf("\nOS-CFAR Results:\n");
    for (size_t i = 0; i < nbins; i++) {
        if (det_os[i]) printf("  Detection at bin %zu: power=%.1f\n", i, rp[i]);
    }

    /* GO-CFAR (best for clutter edges) */
    cfar_config_t cfg_go;
    cfar_config_init(&cfg_go, CFAR_GO, 4, 16, 1e-4);
    int *det_go = calloc(nbins, sizeof(int));
    cfar_detect(rp, nbins, &cfg_go, det_go);

    printf("\nGO-CFAR Results:\n");
    for (size_t i = 0; i < nbins; i++) {
        if (det_go[i]) printf("  Detection at bin %zu: power=%.1f\n", i, rp[i]);
    }

    /* ROC curve analysis */
    printf("\nROC Curve (10 dB SNR, N=1, Swerling 0):\n");
    printf("  P_fa          P_d\n");
    printf("  ----          ---\n");

    double pfa_arr[10], pd_arr[10];
    radar_roc_curve(10.0, pfa_arr, pd_arr, 10, 1, 0);
    for (int i = 0; i < 10; i++) {
        printf("  %.1e     %.4f\n", pfa_arr[i], pd_arr[i]);
    }

    /* Albersheim: required SNR for P_d=0.9, P_fa=1e-6 */
    double snr_req = radar_albersheim_snr(0.9, 1e-6, 1);
    printf("\nAlbersheim: Required SNR for Pd=0.9, Pfa=1e-6, N=1: %.1f dB\n",
           snr_req);

    /* Fluctuation loss: compare Swerling 0 vs Swerling I */
    double pd_sw0 = radar_pd_marcum(10.0, 1e-4, 1);
    double pd_sw1 = radar_pd_swerling1(10.0, 1e-4);
    printf("\nFluctuation loss at SNR=10dB, Pfa=1e-4:\n");
    printf("  Swerling 0 (non-fluctuating): Pd = %.4f\n", pd_sw0);
    printf("  Swerling I (slow Rayleigh):   Pd = %.4f\n", pd_sw1);
    printf("  Fluctuation loss: Pd ratio = %.2f\n", pd_sw0 / pd_sw1);

    free(rp); free(det_ca); free(det_os); free(det_go);

    printf("\nCFAR example complete.\n");
    return 0;
}
