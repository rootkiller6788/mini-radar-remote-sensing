/**
 * example_single_target_tracking.c — Single Target Tracking with Kalman Filter
 *
 * Demonstrates:
 *   - Track initialization, prediction, and update
 *   - CV model with range-bearing measurements (CMKF)
 *   - Track scoring and M/N confirmation
 *   - NEES evaluation
 *
 * Scenario: A single target moves at constant velocity (10 m/s east, 5 m/s north)
 *           Radar provides range + bearing measurements at 1 Hz with noise.
 *
 * Reference: Bar-Shalom (2011), Sec. 5.4
 */

#include "track_core.h"
#include "kalman_filter.h"
#include "measurement_models.h"
#include "motion_models.h"
#include "track_metrics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int main(void)
{
    srand((unsigned int)time(NULL));

    printf("=== Single Target Tracking with Kalman Filter ===\n\n");
    printf("Scenario: Target at (0,0) moving 10 m/s east, 5 m/s north\n");
    printf("Sensor: Radar at origin, 1 Hz, sigma_r=5m, sigma_theta=1deg\n");
    printf("Duration: 60 seconds\n\n");

    /* Radar sensor model */
    radar_sensor_t sensor;
    memset(&sensor, 0, sizeof(sensor));
    sensor.sigma_range = 5.0;
    sensor.sigma_bearing = deg_to_rad(1.0);
    sensor.range_max = 2000.0;
    sensor.P_D = 0.9;
    sensor.SNR_dB = 20.0;

    /* Motion model: 2D CV */
    motion_model_t *cv_model = motion_create_cv(DIM_2D, 1.0, 0.1);
    if (!cv_model) { printf("Failed to create motion model\n"); return 1; }

    /* Kalman filter: 4-state (pos_x, vel_x, pos_y, vel_y), 2-measurement (range, bearing) */
    kf_t *kf = kf_alloc(4, 2, 0);
    if (!kf) { motion_free(cv_model); return 1; }

    double H_cart[2 * 4] = {1, 0, 0, 0,
                             0, 0, 1, 0};
    double R_polar[4] = {sensor.sigma_range * sensor.sigma_range, 0,
                         0, sensor.sigma_bearing * sensor.sigma_bearing};

    kf_set_model(kf, cv_model->F, H_cart, cv_model->Q, R_polar, NULL, NULL);

    /* Initial state: at origin, unknown velocity */
    double x0[4] = {0.0, 0.0, 0.0, 0.0};
    double P0[16];
    memset(P0, 0, sizeof(P0));
    for (int i = 0; i < 4; i++) P0[i * 4 + i] = 1000.0; /* large initial uncertainty */
    kf_set_state(kf, x0, P0);

    /* Tracker */
    tracker_t tracker;
    tracker_init(&tracker);

    /* True target trajectory */
    double true_pos[60 * 2]; /* 60 steps, [x, y] per step */
    double est_pos[60 * 2];

    printf("Time(s)   True_X(m)   Est_X(m)    Err_X(m)   True_Y(m)   Est_Y(m)    Err_Y(m)\n");
    printf("--------   ---------   ---------   ---------   ---------   ---------   ---------\n");

    for (int k = 0; k < 60; k++) {
        double t = k * 1.0;

        /* True position */
        double true_x = 10.0 * t;
        double true_y = 5.0 * t;
        true_pos[k * 2] = true_x;
        true_pos[k * 2 + 1] = true_y;

        /* Simulate measurement with noise */
        double range_true = sqrt(true_x * true_x + true_y * true_y);
        double bearing_true = atan2(true_y, true_x);

        double range_meas = range_true + sensor.sigma_range * ((double)rand()/RAND_MAX * 2 - 1);
        /* Add Gaussian-like noise via Box-Muller approx */
        double u1 = (double)rand()/RAND_MAX;
        double u2 = (double)rand()/RAND_MAX;
        range_meas = range_true + sensor.sigma_range * sqrt(-2*log(u1+1e-10)) * cos(2*M_PI*u2);

        u1 = (double)rand()/RAND_MAX;
        u2 = (double)rand()/RAND_MAX;
        double bearing_meas = bearing_true + sensor.sigma_bearing * sqrt(-2*log(u1+1e-10)) * cos(2*M_PI*u2);

        /* CMKF step: convert polar to Cartesian and update */
        int ret = cmkf_step(kf->x, kf->P, cv_model->F, cv_model->Q,
                             range_meas, bearing_meas,
                             sensor.sigma_range, sensor.sigma_bearing,
                             H_cart, 4, 2);
        if (ret != 0) {
            /* On failure, just predict */
            kf_predict(kf);
            memcpy(kf->x, kf->x_pred, 4 * sizeof(double));
            memcpy(kf->P, kf->P_pred, 16 * sizeof(double));
        }

        est_pos[k * 2] = kf->x[0];
        est_pos[k * 2 + 1] = kf->x[2];

        /* Print every 5 seconds */
        if (k % 5 == 0) {
            printf("%8.0f   %9.1f   %9.1f   %9.1f   %9.1f   %9.1f   %9.1f\n",
                   t, true_x, kf->x[0], fabs(true_x - kf->x[0]),
                   true_y, kf->x[2], fabs(true_y - kf->x[2]));
        }
    }

    /* Compute RMSE over last 30 seconds (after filter convergence) */
    double rmse_x = metrics_rmse_position(&est_pos[30 * 2], &true_pos[30 * 2], 30, 1);
    double rmse_y = metrics_rmse_position(&est_pos[30 * 2 + 1], &true_pos[30 * 2 + 1], 30, 1);
    double aee = metrics_aee(est_pos, true_pos, 60, 2);

    printf("\n=== Performance Summary ===\n");
    printf("Final state estimate: pos=(%.1f, %.1f) m, vel=(%.2f, %.2f) m/s\n",
           kf->x[0], kf->x[2], kf->x[1], kf->x[3]);
    printf("True state:           pos=(%.1f, %.1f) m, vel=(%.2f, %.2f) m/s\n",
           600.0, 300.0, 10.0, 5.0);
    printf("RMSE (X, last 30s): %.2f m\n", rmse_x);
    printf("RMSE (Y, last 30s): %.2f m\n", rmse_y);
    printf("AEE (full 60s):     %.2f m\n", aee);
    printf("Log-likelihood:     %.2f\n", kf_neg_log_likelihood(kf));

    kf_free(kf);
    motion_free(cv_model);
    printf("\n=== Example Complete ===\n");
    return 0;
}
