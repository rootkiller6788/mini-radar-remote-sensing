/**
 * @file    example_dem_generation.c
 * @brief   End-to-end example: DEM generation from simulated airborne LiDAR
 *
 * Simulates an airborne LiDAR scan over synthetic terrain with buildings,
 * then generates a Digital Terrain Model (DTM) and Digital Surface Model (DSM),
 * and computes the Canopy Height Model (CHM) difference.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lidar_core.h"
#include "lidar_geometry.h"
#include "lidar_applications.h"

int main(void) {
    printf("=== LiDAR DEM Generation Example ===\n\n");

    /* Create a simulated airborne scan over 100m x 100m terrain.
       Terrain: gentle hill with buildings at (20,20) and (60,70). */
    lidar_scan_t scan;
    lidar_scan_init(&scan, 10000);
    printf("Generating synthetic terrain with buildings...\n");

    for (int i = 0; i < 2000; i++) {
        double x = (double)(rand() % 1000) / 10.0;
        double y = (double)(rand() % 1000) / 10.0;
        double z_terrain = 10.0 + 2.0 * sin(x / 30.0) * cos(y / 25.0);

        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = z_terrain;
        pt.return_num = 1;
        pt.num_returns = 1;
        pt.class_label = 2;
        lidar_scan_add_point(&scan, &pt);
    }

    /* Building 1 at (20, 20): 15m x 10m, height ~12m */
    for (int i = 0; i < 500; i++) {
        double x = 20.0 + (double)(rand() % 150) / 10.0;
        double y = 20.0 + (double)(rand() % 100) / 10.0;
        double z = 10.0 + 2.0 * sin(x / 30.0) * cos(y / 25.0) + 12.0;
        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = z;
        pt.return_num = 1;
        pt.num_returns = 1;
        pt.class_label = 6;
        lidar_scan_add_point(&scan, &pt);
    }

    /* Building 2 at (60, 70): 10m x 15m, height ~8m */
    for (int i = 0; i < 400; i++) {
        double x = 60.0 + (double)(rand() % 100) / 10.0;
        double y = 70.0 + (double)(rand() % 150) / 10.0;
        double z = 10.0 + 2.0 * sin(x / 30.0) * cos(y / 25.0) + 8.0;
        lidar_point_t pt;
        memset(&pt, 0, sizeof(pt));
        pt.x = x; pt.y = y; pt.z = z;
        pt.return_num = 1;
        pt.num_returns = 1;
        pt.class_label = 6;
        lidar_scan_add_point(&scan, &pt);
    }

    printf("Generated %zu points\n", scan.num_points);

    /* Configure DEM */
    lidar_dem_config_t dem_cfg;
    lidar_dem_config_from_scan(&scan, 2.0, &dem_cfg);
    printf("DEM grid: %zu x %zu cells at %.1f m resolution\n",
           dem_cfg.cols, dem_cfg.rows, dem_cfg.cell_size);

    /* Generate DTM (ground points only) */
    lidar_dem_t dtm;
    if (lidar_dem_generate(&scan, &dem_cfg, &dtm) == 0) {
        printf("DTM generated: %zu x %zu cells\n", dtm.cols, dtm.rows);

        /* Compute elevation statistics */
        double min_el = 1e100, max_el = -1e100;
        size_t valid = 0;
        for (size_t i = 0; i < dtm.cols * dtm.rows; i++) {
            if (!isnan(dtm.data[i])) {
                if (dtm.data[i] < min_el) min_el = dtm.data[i];
                if (dtm.data[i] > max_el) max_el = dtm.data[i];
                valid++;
            }
        }
        printf("DTM elevation range: %.2f — %.2f m (%zu valid cells)\n",
               min_el, max_el, valid);
        lidar_dem_free(&dtm);
    }

    /* Generate DSM (all points) */
    lidar_dem_t dsm;
    if (lidar_dsm_generate(&scan, &dem_cfg, &dsm) == 0) {
        printf("DSM generated: %zu x %zu cells\n", dsm.cols, dsm.rows);
        lidar_dem_free(&dsm);
    }

    /* Extract buildings */
    lidar_building_footprint_t footprints[20];
    size_t n_bldg = lidar_extract_buildings(&scan, 2.0, 4.0, 3.0,
                                              footprints, 20);
    printf("\nBuildings detected: %zu\n", n_bldg);
    for (size_t i = 0; i < n_bldg; i++) {
        printf("  Building %zu: center=(%.1f, %.1f), "
               "area=%.1f m^2, height=%.1f m\n",
               i + 1,
               (footprints[i].x_min + footprints[i].x_max) / 2.0,
               (footprints[i].y_min + footprints[i].y_max) / 2.0,
               footprints[i].area,
               footprints[i].height);
    }

    printf("\nDEM generation example completed.\n");
    lidar_scan_free(&scan);
    return 0;
}