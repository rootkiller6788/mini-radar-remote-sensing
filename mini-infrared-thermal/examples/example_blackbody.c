/**
 * @file    example_blackbody.c
 * @brief   Example: Blackbody radiation spectrum computation
 *
 * Demonstrates Planck's law, Wien's displacement, and Stefan-Boltzmann
 * by computing spectral radiance curves for bodies at different temperatures
 * (human body, boiling water, fire, industrial furnace, Sun).
 *
 * L6 - Canonical Problem: Blackbody spectrum computation and verification
 * L7 - Application: Non-contact temperature measurement
 */

#include <stdio.h>
#include <math.h>
#include "ir_core.h"

int main(void) {
    printf("=== Blackbody Radiation Spectrum ===");

    /* Compute spectra for various temperatures */
    double temps[] = {310.0, 373.0, 800.0, 1500.0, 5778.0};
    const char *labels[] = {
        "Human body (310K)",
        "Boiling water (373K)",
        "Fire (800K)",
        "Industrial furnace (1500K)",
        "Sun surface (5778K)"
    };

    for (int t = 0; t < 5; t++) {
        double T = temps[t];
        double lambda_peak = ir_wien_peak_wavelength(T);
        double M = ir_stefan_boltzmann_exitance(T);
        double B_peak = ir_planck_peak_radiance(T);

        printf("=== %s ===", labels[t]);
        printf("  Peak wavelength: %.2f um", lambda_peak * 1e6);
        printf("  Radiant exitance: %.2f W/m2", M);
        printf("  Peak spectral radiance: %.2e W/(sr*m3)", B_peak);
        printf("  Total power from 1 m2 blackbody: %.2f kW", M / 1000.0);

        /* Print spectral radiance at key wavelengths */
        printf("  Spectral radiance at key wavelengths:");
        double lambdas[] = {0.5e-6, 1e-6, 3e-6, 5e-6, 8e-6, 10e-6, 14e-6};
        for (int i = 0; i < 7; i++) {
            double B = ir_planck_spectral_radiance_wavelength(lambdas[i], T);
            printf("    %.1f um: %.3e W/(sr*m3)", lambdas[i] * 1e6, B);
        }
        printf("");
    }

    /* Verify Wien's law */
    printf("=== Wien's Law Verification ===");
    for (double T = 200.0; T <= 1000.0; T += 200.0) {
        double lambda_max = ir_wien_peak_wavelength(T);
        printf("  T=%.0fK: lambda_max*T = %.2f um*K (expected ~2898)",               T, lambda_max * T * 1e6);
    }

    /* Verify Stefan-Boltzmann by numerical integration */
    printf("=== Stefan-Boltzmann Numerical Verification ===");
    double err = ir_stefan_boltzmann_verify(500.0, 1000);
    printf("  Numerical integration error at 500K: %.4f%%", err * 100.0);

    return 0;
}
