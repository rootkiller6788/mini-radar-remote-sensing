/**
 * @file    example_object_detection.c
 * @brief   End-to-end example: Object detection for autonomous driving
 *
 * Simulates an automotive LiDAR scan with vehicles, pedestrians, and clutter,
 * then performs Euclidean clustering to detect objects.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lidar_core.h"
#include "lidar_geometry.h"
#include "lidar_detection.h"
#include "lidar_applications.h"

int main(void) {
    printf("=== Automotive LiDAR Object Detection Example ===\n\n");

    /* Configure an automotive LiDAR */
    lidar_config_t config;
    lidar_config_init_automotive(&config);
    printf("LiDAR config: %d nm, max range %.0f m, PRF %.0f kHz\n",
           config.wavelength, config.range_max, config.prf / 1000.0);
    printf("Range resolution: %.3f m\n", lidar_range_resolution(&config));

    /* Create a simulated scan with objects */
    lidar_scan_t scan;
    lidar_scan_init(&scan, 5000);

    /* Ground surface */
    for (int i = 0; i < 1000; i++) {
        double x = (double)(rand() % 1000) / 10.0 - 20.0;
        double y = (double)(rand() % 2000) / 10.0;
        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = -1.5;
        pt.intensity = 0.1;
        pt.class_label = 2;
        lidar_scan_add_point(&scan, &pt);
    }

    /* Vehicle 1: car at (10, 5), 4.5m x 1.8m x 1.5m */
    for (int i = 0; i < 200; i++) {
        double x = 10.0 + (double)(rand() % 45) / 10.0 - 2.25;
        double y = 5.0 + (double)(rand() % 18) / 10.0 - 0.9;
        double z = (double)(rand() % 15) / 10.0 - 1.5;
        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = z;
        pt.intensity = 0.4;
        lidar_scan_add_point(&scan, &pt);
    }

    /* Vehicle 2: truck at (5, 25), 8m x 2.5m x 3m */
    for (int i = 0; i < 300; i++) {
        double x = 5.0 + (double)(rand() % 80) / 10.0 - 4.0;
        double y = 25.0 + (double)(rand() % 25) / 10.0 - 1.25;
        double z = (double)(rand() % 30) / 10.0 - 1.5;
        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = z;
        pt.intensity = 0.5;
        lidar_scan_add_point(&scan, &pt);
    }

    /* Pedestrian at (15, 15) */
    for (int i = 0; i < 30; i++) {
        double x = 15.0 + (double)(rand() % 4) / 10.0 - 0.2;
        double y = 15.0 + (double)(rand() % 4) / 10.0 - 0.2;
        double z = (double)(rand() % 18) / 10.0 - 1.5;
        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = z;
        pt.intensity = 0.3;
        lidar_scan_add_point(&scan, &pt);
    }

    printf("Generated %zu points\n", scan.num_points);

    /* Detect objects using Euclidean clustering */
    printf("\nDetecting objects (Euclidean clustering, 2.0m threshold)...\n");
    lidar_object_detection_t detections[20];
    size_t n_det = lidar_detect_objects_euclidean(&scan, 2.0, 5, 500,
                                                     detections, 20);

    printf("Objects detected: %zu\n\n", n_det);
    for (size_t i = 0; i < n_det; i++) {
        const char *class_name;
        double w = detections[i].width;
        double l = detections[i].length;
        double h = detections[i].height;

        /* Simple heuristic classification based on dimensions */
        if (h > 2.5 && l > 6.0) {
            class_name = "Truck";
        } else if (w > 1.5 && l > 3.0) {
            class_name = "Car";
        } else if (h > 1.0 && h < 2.5) {
            class_name = "Pedestrian";
        } else {
            class_name = "Unknown";
        }

        printf("  Object %zu: %s\n", i + 1, class_name);
        printf("    Position: (%.1f, %.1f)\n",
               detections[i].x_center, detections[i].y_center);
        printf("    Dimensions: %.1f x %.1f x %.1f m\n", w, l, h);
        printf("    Confidence: %.2f\n\n", detections[i].confidence);
    }

    /* Compute SNR for a 10% reflective target at 50m */
    lidar_target_t target;
    memset(&target, 0, sizeof(target));
    target.reflectivity = 0.1;
    target.area = 1.0;

    lidar_atmosphere_t atm;
    memset(&atm, 0, sizeof(atm));
    atm.visibility = 23000.0;
    atm.extinction_coeff = lidar_atmosphere_extinction(atm.visibility, config.wavelength);

    double snr = lidar_snr(&config, &target, &atm, 50.0);
    printf("SNR for 10%% reflective 1m^2 target at 50m: %.1f dB\n", snr);

    double pd = lidar_prob_detection(snr, 1e-4);
    printf("Probability of detection (PFA=1e-4): %.4f\n", pd);

    printf("\nObject detection example completed.\n");
    lidar_scan_free(&scan);
    return 0;
}