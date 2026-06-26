/**
 * @file    example_fever_screening.c
 * @brief   Example: Fever Screening using IR Thermography
 *
 * Demonstrates the use of brightness temperature, NETD, and thermal
 * contrast for medical fever screening applications.
 *
 * L6 - Canonical Problem: Non-contact human body temperature measurement
 * L7 - Application: Public health screening (e.g., airport thermal scanners)
 *
 * Reference: ISO/TR 13154:2017 "Medical electrical equipment - Deployment,
 *   implementation and operational guidelines for screening of humans
 *   for fever using a screening thermograph"
 */

#include <stdio.h>
#include <math.h>
#include "ir_core.h"
#include "ir_detector.h"
#include "ir_radiometry.h"

int main(void) {
    printf("=== Fever Screening with IR Thermography ===");

    /* Scenario: Airport thermal screening */
    double T_normal = 310.0;      /* Normal body temp (37C, 98.6F) */
    double T_fever = 311.5;       /* Fever threshold (38.5C, 101.3F) */
    double T_bg = 295.0;          /* Background (22C, 72F) */

    printf("Scenario: Airport thermal screening system");
    printf("  Normal body temp: %.1f K (%.1f C)", T_normal, T_normal - 273.15);
    printf("  Fever threshold:  %.1f K (%.1f C)", T_fever, T_fever - 273.15);
    printf("  Background temp:  %.1f K (%.1f C)", T_bg, T_bg - 273.15);

    /* Thermal contrast */
    double contrast_normal = ir_thermal_contrast(T_normal, T_bg);
    double contrast_fever = ir_thermal_contrast(T_fever, T_bg);
    printf("Thermal contrast:");
    printf("  Normal vs bg: %.4f", contrast_normal);
    printf("  Fever vs bg:  %.4f", contrast_fever);

    /* Radiance in LWIR band */
    double L_normal = ir_in_band_radiance(8e-6, 14e-6, T_normal, 500);
    double L_fever  = ir_in_band_radiance(8e-6, 14e-6, T_fever, 500);
    printf("In-band radiance (8-14 um):");
    printf("  Normal: %.2f W/(sr*m2)", L_normal);
    printf("  Fever:  %.2f W/(sr*m2)", L_fever);
    printf("  Difference: %.2f W/(sr*m2) (%.1f%%)",           L_fever - L_normal, 100.0 * (L_fever - L_normal) / L_normal);

    /* NETD requirement */
    ir_detector_params_t det;
    ir_detector_params_init(&det, IR_DET_THERMAL_BOLOMETER);
    double dLdT = ir_radiance_derivative_wrt_temperature(T_normal, 8.0, 14.0, 500);
    double I_noise = ir_noise_total_current(&det, 30.0, 30.0);
    double R_i = ir_detector_responsivity_current(det.quantum_efficiency, 10.0);
    double nep = I_noise / R_i;
    double netd = ir_netd(1.0, 0.9, det.active_area_cm2, nep, dLdT);

    printf("System sensitivity:");
    printf("  dL/dT at 310K: %.3f W/(sr*m2*K)", dLdT);
    printf("  Detector NEP:  %.2e W", nep);
    printf("  NETD:          %.1f mK", netd * 1000.0);

    if (netd * 1000.0 < (T_fever - T_normal) * 1000.0 / 2.0) {
        printf("CONCLUSION: System is SUITABLE for fever screening.");
        printf("  NETD (%.1f mK) < required (%.1f mK)",               netd * 1000.0, (T_fever - T_normal) * 1000.0 / 2.0);
    } else {
        printf("CONCLUSION: System NOT suitable. NETD too large.");
    }

    /* Atmospheric effects */
    printf("Atmospheric attenuation effects:");
    for (double tau = 0.7; tau <= 1.0; tau += 0.1) {
        double dT_app = ir_apparent_temperature_difference(
            T_fever - T_normal, tau);
        printf("  tau=%.1f: apparent dT = %.3f K (%.1f C)",               tau, dT_app, dT_app);
    }

    return 0;
}
