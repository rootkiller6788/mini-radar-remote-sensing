/**
 * @file    demo_sar.c
 * @brief   SAR Visualization Demo
 * Demonstrates point target simulation and ASCII-art range profile.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sar_core.h"
#include "sar_geometry.h"

int main(void) {
    printf("SAR Point Target Visualization\n");
    printf("==============================\n\n");

    /* Generate chirp and show its autocorrelation */
    sar_chirp_t *c = sar_chirp_alloc(128, 20e-6, 50e6, 100e6);
    size_t M = 2*c->num_samples-1;
    double *R = malloc(M*8);
    sar_chirp_autocorrelation(c, R);

    printf("Chirp Autocorrelation (matched filter output):\n");
    for (size_t i=0; i<M; i+=4) {
        int bar = (int)(R[i]*50.0);
        if (bar>50) bar=50;
        printf("%3zu: ", i);
        for (int j=0; j<bar; j++) printf("#");
        printf("\n");
    }

    /* SAR parameter summary */
    printf("\nTypical Spaceborne SAR Parameters:\n");
    printf("  TerraSAR-X:  X-band, 9.65 GHz, 1m res\n");
    printf("  Radarsat-2:  C-band, 5.405 GHz, 3m res\n");
    printf("  ALOS-2:      L-band, 1.27 GHz, 10m res\n");
    printf("  Sentinel-1:  C-band, 5.405 GHz, 5m res\n");

    free(R); sar_chirp_free(c);
    printf("\nDemo complete.\n");
    return 0;
}