/**
 * @file    example_building_inspection.c
 * @brief   Example: Building Thermal Inspection
 *
 * Demonstrates thermal image analysis for building energy audits:
 * detecting heat loss, insulation defects, and thermal bridging.
 *
 * L6 - Canonical Problem: Thermal pattern analysis for defect detection
 * L7 - Application: Building energy efficiency (ISO 6781, EN 13187)
 *
 * Reference: ISO 6781:1983 "Thermal insulation - Qualitative detection of
 *   thermal irregularities in building envelopes - Infrared method"
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ir_core.h"
#include "ir_image.h"
#include "ir_detection.h"

int main(void) {
    printf("=== Building Thermal Inspection ===");

    /* Simulated thermal image of a wall section (8x8 grid) */
    /* Temperatures in Celsius -> convert to Kelvin */
    double wall_data_8x8[] = {
        18.5, 18.3, 18.1, 17.9, 17.8, 17.6, 17.5, 17.4,
        18.4, 18.2, 22.5, 23.1, 22.8, 17.5, 17.4, 17.3,
        18.3, 18.1, 23.2, 24.5, 23.9, 17.6, 17.5, 17.4,
        18.2, 18.0, 22.8, 23.8, 23.2, 17.7, 17.6, 17.5,
        18.1, 17.9, 17.8, 17.7, 17.6, 17.8, 17.7, 17.6,
        18.0, 17.8, 17.7, 17.6, 17.5, 17.9, 18.0, 17.9,
        17.9, 17.7, 17.6, 17.5, 17.4, 17.8, 17.9, 18.0,
        17.8, 17.6, 17.5, 17.4, 17.3, 17.7, 17.8, 17.9
    };

    int rows = 8, cols = 8;
    int N = rows * cols;

    /* Convert C to K */
    double *wall_K = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
        wall_K[i] = wall_data_8x8[i] + 273.15;

    /* Compute statistics */
    double sum = 0.0, sum2 = 0.0;
    double min_v = wall_K[0], max_v = wall_K[0];
    for (int i = 0; i < N; i++) {
        double v = wall_K[i];
        sum += v;
        sum2 += v * v;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    double mean = sum / N;
    double std = sqrt(sum2 / N - mean * mean);

    printf("Wall thermal image statistics:");
    printf("  Min: %.1f C  Max: %.1f C", min_v - 273.15, max_v - 273.15);
    printf("  Mean: %.1f C  StdDev: %.2f C", mean - 273.15, std);
    printf("  Temperature span: %.1f C", max_v - min_v);

    /* Detect thermal anomalies using Otsu's method */
    double thresh = ir_otsu_threshold(wall_K, N, 64);
    printf("Otsu threshold: %.1f C", thresh - 273.15);

    /* Segment the image */
    int *binary = (int*)malloc(N * sizeof(int));
    ir_threshold_segment(wall_K, N, thresh, binary);

    printf("Thermal anomaly map (X = anomaly):");
    for (int r = 0; r < rows; r++) {
        printf("  ");
        for (int c = 0; c < cols; c++) {
            int idx = r * cols + c;
            printf("%s ", binary[idx] ? "X" : ".");
        }
        printf("");
    }

    /* Check for thermal bridging (linear patterns) */
    int anomaly_count = 0;
    for (int i = 0; i < N; i++)
        if (binary[i]) anomaly_count++;

    double anomaly_pct = 100.0 * anomaly_count / N;
    printf("Anomaly coverage: %d/%d pixels (%.1f%%)",           anomaly_count, N, anomaly_pct);

    /* Radiometric analysis */
    double L_wall = ir_in_band_radiance(8e-6, 14e-6, mean, 200);
    double L_anomaly = ir_in_band_radiance(8e-6, 14e-6, max_v, 200);
    double contrast = ir_thermal_contrast(max_v, mean);

    printf("Radiometric analysis (8-14 um LWIR):");
    printf("  Wall radiance: %.2f W/(sr*m2)", L_wall);
    printf("  Anomaly radiance: %.2f W/(sr*m2)", L_anomaly);
    printf("  Thermal contrast: %.2f", contrast);

    double delta_T = max_v - mean;
    printf("  Temperature anomaly: %.1f C", delta_T);
    if (delta_T > 3.0) {
        printf("  SEVERE: Significant heat loss detected! (>3C anomaly)");
    } else if (delta_T > 1.0) {
        printf("  MODERATE: Insulation defect suspected. (>1C anomaly)");
    } else {
        printf("  NORMAL: Within acceptable range.");
    }

    /* Energy loss estimate */
    double sigma_T4_diff = IR_STEFAN_BOLTZMANN
                           * (max_v * max_v * max_v * max_v
                              - mean * mean * mean * mean);
    printf("Radiative heat loss from anomaly: %.1f W/m2",           sigma_T4_diff);

    free(wall_K);
    free(binary);
    return 0;
}
