/**
 * example_automotive_radar.c — Automotive Radar Target Tracking
 *
 * Demonstrates:
 *   - DBSCAN clustering of radar point cloud
 *   - Object classification based on radar features
 *   - Ego-motion compensation for moving vehicle
 *   - Occupancy grid mapping (OGM)
 *   - Group detection of convoy vehicles
 *
 * Scenario: Ego vehicle at 20 m/s approaches 2 vehicles ahead.
 *           Front radar at 77 GHz, 1 GHz bandwidth, 10 Hz update.
 *           Vehicle 1: 40m ahead, 15 m/s (slower)
 *           Vehicle 2: 60m ahead, 20 m/s (same speed)
 *
 * Reference: TI AWR1642 SDK, Werber et al. (2015)
 */

#include "mmwave_tracker.h"
#include "motion_models.h"
#include "track_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

static double gauss_rand(double mean, double std)
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z = sqrt(-2.0 * log(u1 + 1e-10)) * cos(2.0 * M_PI * u2);
    return mean + std * z;
}

int main(void)
{
    srand((unsigned int)time(NULL));

    printf("=== Automotive Radar Target Tracking ===\n\n");
    printf("Configuration:\n");
    printf("  Radar: 77 GHz, BW=1 GHz, PRF=5 kHz, 10 Hz update\n");
    printf("  Ego speed: 20 m/s (72 km/h)\n");
    printf("  Vehicle 1: 40 m ahead, 15 m/s (slower)\n");
    printf("  Vehicle 2: 60 m ahead, 20 m/s (same speed)\n\n");

    /* --- Part 1: DBSCAN clustering --- */
    printf("--- Part 1: DBSCAN Point Cloud Clustering ---\n");

    /* Simulate radar point cloud with two vehicles + noise */
    mmwave_point_t points[200];
    int n_points = 0;

    /* Vehicle 1 at (40, 0): ~30 points */
    for (int i = 0; i < 30; i++) {
        points[n_points].x = (float)gauss_rand(40.0, 0.3);
        points[n_points].y = (float)gauss_rand(0.0, 0.5);
        points[n_points].velocity = (float)gauss_rand(-5.0, 0.1); /* relative = 20-15 = 5 approaching */
        points[n_points].snr = (float)gauss_rand(15.0, 3.0);
        n_points++;
    }

    /* Vehicle 2 at (60, 1): ~40 points (larger truck) */
    for (int i = 0; i < 40; i++) {
        points[n_points].x = (float)gauss_rand(60.0, 0.5);
        points[n_points].y = (float)gauss_rand(1.0, 1.0);
        points[n_points].velocity = (float)gauss_rand(0.0, 0.2); /* relative = 0 */
        points[n_points].snr = (float)gauss_rand(20.0, 4.0);
        n_points++;
    }

    /* Noise/clutter points */
    for (int i = 0; i < 20; i++) {
        points[n_points].x = (float)(10.0 + (rand() % 800) / 10.0);
        points[n_points].y = (float)((rand() % 200 - 100) / 10.0);
        points[n_points].velocity = (float)((rand() % 200 - 100) / 10.0);
        points[n_points].snr = (float)(2.0 + (rand() % 50) / 10.0);
        n_points++;
    }

    printf("  Total points: %d\n", n_points);

    /* DBSCAN clustering */
    double eps = 1.5; /* 1.5m neighborhood for vehicles */
    int minPts = 5;
    dbscan_result_t *cluster = dbscan_cluster(points, n_points, eps, minPts);

    if (cluster) {
        printf("  Clusters found: %d\n", cluster->num_clusters);

        /* Count points per cluster */
        for (int c = 1; c <= cluster->num_clusters; c++) {
            int count = 0;
            for (int i = 0; i < n_points; i++) {
                if (cluster->labels[i] == c) count++;
            }
            printf("    Cluster %d: %d points\n", c, count);
        }

        int noise_count = 0;
        for (int i = 0; i < n_points; i++) {
            if (cluster->labels[i] == 0) noise_count++;
        }
        printf("  Noise points: %d\n", noise_count);

        /* Compute centroids for each cluster */
        for (int c = 1; c <= cluster->num_clusters; c++) {
            double cx = 0.0, cy = 0.0;
            int count = 0;
            for (int i = 0; i < n_points; i++) {
                if (cluster->labels[i] == c) {
                    cx += points[i].x;
                    cy += points[i].y;
                    count++;
                }
            }
            cx /= count;
            cy /= count;
            printf("    Centroid %d: (%.1f, %.1f) m\n", c, cx, cy);

            /* Classify */
            object_class_t cls = classify_radar_object(points, n_points,
                                                       cluster->labels, cx, cy);
            const char *cls_names[] = {"Unknown", "Pedestrian", "Bicycle",
                                       "Car", "Truck", "Static"};
            printf("      Classification: %s\n", cls_names[cls]);
        }
    }
    printf("\n");

    /* --- Part 2: Ego-motion compensation --- */
    printf("--- Part 2: Ego-Motion Compensation ---\n");

    mmwave_point_t moving_test[3] = {
        {10.0f, 0.0f, 0.0f, -20.0f, 15.0f, 0},  /* stationary object ahead */
        {10.0f, 2.0f, 0.0f, -18.0f, 12.0f, 0},  /* should be stationary */
        {10.0f, 0.0f, 0.0f, -15.0f, 10.0f, 0}   /* actually moving away */
    };

    printf("  Before compensation:\n");
    for (int i = 0; i < 3; i++) {
        printf("    Point %d: v=%.1f m/s\n", i, moving_test[i].velocity);
    }

    ego_motion_compensate(moving_test, 3, 20.0, 0.0, 0.0, 0.0, 0.0);
    /* Ego moving at 20 m/s forward, stationary objects should show 0 radial velocity */

    printf("  After compensation (ego v=20 m/s):\n");
    for (int i = 0; i < 3; i++) {
        printf("    Point %d: v=%.1f m/s\n", i, moving_test[i].velocity);
    }

    int n_moving = detect_moving_points(moving_test, 3, 1.0);
    printf("  Moving points detected: %d\n\n", n_moving);

    /* --- Part 3: Occupancy Grid Mapping --- */
    printf("--- Part 3: Occupancy Grid Mapping ---\n");

    occupancy_grid_t *grid = ogm_create(50, 30, 2.0, 0.0, -15.0);
    if (grid) {
        printf("  Grid: %d x %d cells, %.1f m resolution\n",
               grid->width, grid->height, grid->resolution);
        printf("  Coverage: %.0f m x %.0f m\n",
               grid->width * grid->resolution,
               grid->height * grid->resolution);

        /* Update with radar scan */
        ogm_update_radar_scan(grid, points, n_points, 0.0, 0.0, 0.0);

        /* Query occupancy at vehicle locations */
        printf("  Occupancy at (40, 0): %.3f\n", ogm_query(grid, 40.0, 0.0));
        printf("  Occupancy at (60, 1): %.3f\n", ogm_query(grid, 60.0, 1.0));
        printf("  Occupancy at (20, 0): %.3f (should be low/free)\n",
               ogm_query(grid, 20.0, 0.0));

        /* Extract obstacles */
        double obs_x[100], obs_y[100];
        int n_obs = ogm_extract_obstacles(grid, 0.6, obs_x, obs_y, 100);
        printf("  Obstacles (p > 0.6): %d\n", n_obs);
        for (int i = 0; i < n_obs && i < 5; i++) {
            printf("    Obstacle at (%.1f, %.1f)\n", obs_x[i], obs_y[i]);
        }

        ogm_free(grid);
    }
    printf("\n");

    /* --- Part 4: Group Tracking --- */
    printf("--- Part 4: Group Detection ---\n");

    double object_positions[6 * 2] = {
        40.0, 0.0,   /* Vehicle 1 (group 1) */
        42.0, 1.5,   /* Vehicle 1 trailer (close, same speed) */
        60.0, 1.0,   /* Vehicle 2 (group 2) */
        61.0, -1.0,  /* Vehicle 2 escort */
        100.0, 10.0, /* Solo vehicle far away */
        101.0, 12.0  /* Another solo */
    };

    double object_velocities[6 * 2] = {
        15.0, 0.0,   /* group 1 */
        15.0, 0.5,
        20.0, 0.0,   /* group 2 */
        20.0, -0.5,
        25.0, 5.0,   /* solo */
        25.0, 6.0    /* solo (should form group 3) */
    };

    int groups[6];
    int n_groups = 0;

    detect_groups(object_positions, object_velocities, 6, 5.0, 2.0,
                   groups, &n_groups);

    printf("  Objects: 6, Groups found: %d\n", n_groups);
    for (int i = 0; i < 6; i++) {
        printf("    Object %d at (%.0f, %.0f): Group %d\n",
               i, object_positions[2*i], object_positions[2*i+1], groups[i]);
    }

    /* Compute group bounding boxes */
    for (int g = 0; g < n_groups; g++) {
        double cx, cy, xmin, xmax, ymin, ymax;
        compute_group_bounds(object_positions, 6, groups, g,
                              &cx, &cy, &xmin, &xmax, &ymin, &ymax);
        printf("    Group %d centroid: (%.1f, %.1f), bbox: [%.1f-%.1f, %.1f-%.1f]\n",
               g, cx, cy, xmin, xmax, ymin, ymax);
    }

    /* --- Summary --- */
    printf("\n=== Automotive Radar Example Summary ===\n");
    printf("DBSCAN: %d points in %d clusters\n", n_points,
           cluster ? cluster->num_clusters : 0);
    printf("Ego-motion: Compensated %d points\n", 3);
    printf("OGM: 50x30 grid at 2m resolution\n");
    printf("Groups: %d groups from %d objects\n", n_groups, 6);

    if (cluster) dbscan_free_result(cluster);
    printf("\n=== Example Complete ===\n");
    return 0;
}
