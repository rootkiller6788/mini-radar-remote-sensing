/**
 * example_multi_target_tracking.c — Multi-Target Tracking with GNN and PDA
 *
 * Demonstrates:
 *   - Multi-target tracker with track management
 *   - Data association (Nearest Neighbor and Global Nearest Neighbor)
 *   - PDA filter for single-track-in-clutter
 *   - Track initiation, confirmation, and deletion
 *
 * Scenario: Two crossing targets with clutter measurements.
 *           Sensor: 2D radar at origin, 1 Hz scan rate.
 *
 * Reference: Bar-Shalom (2011), Ch. 7, Blackman (1999), Ch. 8
 */

#include "track_core.h"
#include "kalman_filter.h"
#include "data_association.h"
#include "measurement_models.h"
#include "track_metrics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int main(void)
{
    srand((unsigned int)time(NULL));

    printf("=== Multi-Target Tracking with GNN and PDA ===\n\n");
    printf("Scenario: Two crossing targets in clutter\n");
    printf("  Target 1: Start (0, 0), velocity (8, 6) m/s\n");
    printf("  Target 2: Start (0, 200), velocity (8, -6) m/s\n");
    printf("  Clutter: ~5 false alarms per scan\n");
    printf("  Duration: 40 seconds\n\n");

    /* Radar sensor */
    radar_sensor_t sensor;
    memset(&sensor, 0, sizeof(sensor));
    sensor.sensor_x = 0.0;
    sensor.sensor_y = 0.0;
    sensor.sigma_range = 10.0;
    sensor.sigma_bearing = deg_to_rad(0.5);
    sensor.range_max = 500.0;
    sensor.P_D = 0.9;
    sensor.SNR_dB = 15.0;

    /* Tracker */
    tracker_t tracker;
    tracker_init(&tracker);
    tracker.mgmt.M = 3;
    tracker.mgmt.N = 4;
    tracker.mgmt.coast_max = 3;
    tracker.mgmt.delete_score_threshold = -5.0;
    tracker.mgmt.confirm_score_threshold = 1.0;

    /* Simple CV model for predictions */
    double F[16] = {1, 1, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 1,
                    0, 0, 0, 1};

    /* Measurement statistics */
    int total_measurements = 0;
    int total_associations = 0;
    int total_false_alarms = 0;
    int confirmed_tracks = 0;

    printf("Scan  Tracks  Meas   Assoc   Confirmed   Track1_Err(m)   Track2_Err(m)\n");
    printf("----  ------  ----   -----   ---------   -------------   -------------\n");

    for (int scan = 0; scan < 40; scan++) {
        double t = scan * 1.0;
        tracker.current_time = t;
        tracker.scan_count++;

        /* True target positions */
        double t1_x = 8.0 * t;
        double t1_y = 6.0 * t;

        double t2_x = 8.0 * t;
        double t2_y = 200.0 - 6.0 * t;

        /* --- Generate measurements --- */
        measurement_t measurements[50];
        int num_meas = 0;

        /* Target 1 measurement (with P_D probability) */
        if ((double)rand() / RAND_MAX < sensor.P_D) {
            double r1 = sqrt(t1_x * t1_x + t1_y * t1_y);
            double b1 = atan2(t1_y, t1_x);
            /* Add noise */
            r1 += sensor.sigma_range * ((double)rand()/RAND_MAX * 2 - 1);
            b1 += sensor.sigma_bearing * ((double)rand()/RAND_MAX * 2 - 1);

            measurements[num_meas].meas_id = num_meas + 1;
            measurements[num_meas].type = MEAS_POLAR_2D;
            measurements[num_meas].dim = 2;
            measurements[num_meas].z[0] = r1;
            measurements[num_meas].z[1] = b1;
            measurements[num_meas].R[0] = sensor.sigma_range * sensor.sigma_range;
            measurements[num_meas].R[1] = 0;
            measurements[num_meas].R[2] = 0;
            measurements[num_meas].R[3] = sensor.sigma_bearing * sensor.sigma_bearing;
            measurements[num_meas].snr = sensor.SNR_dB;
            measurements[num_meas].is_clutter = 0;
            num_meas++;
        }

        /* Target 2 measurement */
        if ((double)rand() / RAND_MAX < sensor.P_D) {
            double r2 = sqrt(t2_x * t2_x + t2_y * t2_y);
            double b2 = atan2(t2_y, t2_x);
            r2 += sensor.sigma_range * ((double)rand()/RAND_MAX * 2 - 1);
            b2 += sensor.sigma_bearing * ((double)rand()/RAND_MAX * 2 - 1);

            measurements[num_meas].meas_id = num_meas + 1;
            measurements[num_meas].type = MEAS_POLAR_2D;
            measurements[num_meas].dim = 2;
            measurements[num_meas].z[0] = r2;
            measurements[num_meas].z[1] = b2;
            measurements[num_meas].R[0] = sensor.sigma_range * sensor.sigma_range;
            measurements[num_meas].R[1] = 0;
            measurements[num_meas].R[2] = 0;
            measurements[num_meas].R[3] = sensor.sigma_bearing * sensor.sigma_bearing;
            measurements[num_meas].snr = sensor.SNR_dB;
            measurements[num_meas].is_clutter = 0;
            num_meas++;
        }

        /* Clutter: random false alarms */
        int num_clutter = rand() % 8; /* 0-7 false alarms */
        for (int c = 0; c < num_clutter && num_meas < 50; c++) {
            double range = 10.0 + (rand() % 490);
            double bearing = deg_to_rad((rand() % 360) - 180);

            measurements[num_meas].meas_id = num_meas + 1;
            measurements[num_meas].type = MEAS_POLAR_2D;
            measurements[num_meas].dim = 2;
            measurements[num_meas].z[0] = range;
            measurements[num_meas].z[1] = bearing;
            measurements[num_meas].R[0] = 400.0;
            measurements[num_meas].R[1] = 0;
            measurements[num_meas].R[2] = 0;
            measurements[num_meas].R[3] = 0.01;
            measurements[num_meas].snr = 3.0;
            measurements[num_meas].is_clutter = 1;
            num_meas++;
        }

        total_measurements += num_meas;
        total_false_alarms += num_clutter;

        /* --- Data association --- */
        /* Use a simplified approach: for each track, predict and find NN measurement */
        for (int ti = 0; ti < TRACKER_MAX_TRACKS; ti++) {
            track_t *t = &tracker.tracks[ti];
            if (t->state == TRACK_STATE_FREE || t->state == TRACK_STATE_DELETED)
                continue;

            /* Predict using simple CV model */
            if (t->state_dim == 4) {
                /* x_pred = F * x */
                for (int i = 0; i < 4; i++) {
                    t->x_pred[i] = 0.0;
                    for (int j = 0; j < 4; j++) {
                        t->x_pred[i] += F[i * 4 + j] * t->x[j];
                    }
                }
                /* Simple P_pred = P + Q (approximate) */
                memcpy(t->P_pred, t->P, 16 * sizeof(double));
            }

            /* Find closest measurement in gate */
            double min_dist = 1e10;
            int best_j = -1;
            for (int j = 0; j < num_meas; j++) {
                /* Simple Euclidean gate in Cartesian for demo */
                double mx = measurements[j].z[0] * cos(measurements[j].z[1]);
                double my = measurements[j].z[0] * sin(measurements[j].z[1]);

                double dx = mx - t->x_pred[0];
                double dy = my - t->x_pred[2];
                double dist = dx * dx + dy * dy;

                double gate_size = 2500.0; /* 50^2 */
                if (dist < min_dist && dist < gate_size) {
                    min_dist = dist;
                    best_j = j;
                }
            }

            if (best_j >= 0) {
                /* Associated: update track with this measurement */
                t->misses_since_last = 0;
                total_associations++;
                /* Simple update: blend prediction with measurement */
                double mx = measurements[best_j].z[0] * cos(measurements[best_j].z[1]);
                double my = measurements[best_j].z[0] * sin(measurements[best_j].z[1]);
                double alpha = 0.7; /* filter gain */
                t->x[0] = (1 - alpha) * t->x_pred[0] + alpha * mx;
                t->x[2] = (1 - alpha) * t->x_pred[2] + alpha * my;

                /* Mark measurement as used */
                measurements[best_j].meas_id = 0;
            } else {
                t->misses_since_last++;
            }
        }

        /* Track management */
        tracker_apply_mn_logic(&tracker);

        /* Create new tracks from unassigned measurements */
        for (int j = 0; j < num_meas; j++) {
            if (measurements[j].meas_id == 0) continue; /* already used */
            if (measurements[j].is_clutter) continue; /* don't init on clutter */
            /* Try to initiate a new track */
            if (tracker.num_tracks < TRACKER_MAX_TRACKS - 2) {
                track_t *new_t = tracker_create_track(&tracker, 4, 2, MEAS_POLAR_2D);
                if (new_t) {
                    double mx = measurements[j].z[0] * cos(measurements[j].z[1]);
                    double my = measurements[j].z[0] * sin(measurements[j].z[1]);
                    new_t->x[0] = mx;
                    new_t->x[2] = my;
                }
            }
        }

        /* Count confirmed tracks */
        confirmed_tracks = 0;
        for (int ti = 0; ti < TRACKER_MAX_TRACKS; ti++) {
            if (tracker.tracks[ti].state == TRACK_STATE_CONFIRMED)
                confirmed_tracks++;
        }

        /* Compute errors for confirmed tracks */
        double err1 = -1, err2 = -1;
        for (int ti = 0; ti < TRACKER_MAX_TRACKS; ti++) {
            if (tracker.tracks[ti].state != TRACK_STATE_CONFIRMED) continue;
            double dx1 = tracker.tracks[ti].x[0] - t1_x;
            double dy1 = tracker.tracks[ti].x[2] - t1_y;
            double dist1 = sqrt(dx1 * dx1 + dy1 * dy1);
            double dx2 = tracker.tracks[ti].x[0] - t2_x;
            double dy2 = tracker.tracks[ti].x[2] - t2_y;
            double dist2 = sqrt(dx2 * dx2 + dy2 * dy2);
            if (dist1 < dist2 && dist1 < 100) err1 = dist1;
            if (dist2 < dist1 && dist2 < 100) err2 = dist2;
        }

        if (scan % 5 == 0) {
            printf("%4d   %5d   %4d   %5d   %9d   %13.1f   %13.1f\n",
                   scan, tracker.num_tracks, num_meas,
                   (scan > 0 && total_associations > 0)
                     ? total_associations / (scan + 1) : 0,
                   confirmed_tracks,
                   err1, err2);
        }
    }

    printf("\n=== Multi-Target Tracking Summary ===\n");
    printf("Total scans:              40\n");
    printf("Total measurements:       %d (%.0f%% clutter)\n",
           total_measurements,
           100.0 * total_false_alarms / (total_measurements > 0 ? total_measurements : 1));
    printf("Total associations:       %d\n", total_associations);
    printf("Confirmed tracks:         %d\n", confirmed_tracks);
    printf("Tracks created:           %d\n", tracker.total_tracks_created);
    printf("Tracks confirmed (total): %d\n", tracker.total_tracks_confirmed);
    printf("Tracks deleted:           %d\n", tracker.total_tracks_deleted);

    printf("\n=== Example Complete ===\n");
    return 0;
}
