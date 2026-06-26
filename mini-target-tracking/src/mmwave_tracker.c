/**
 * mmwave_tracker.c — Automotive mmWave Radar Target Tracking Implementation
 *
 * Implements: DBSCAN clustering, object classification, shape computation,
 *             ego-motion compensation, occupancy grid mapping, extended
 *             object tracking, group detection.
 *
 * References:
 *   - Ester et al. "DBSCAN" (1996)
 *   - Koch, J.W. "Extended Object Tracking using Random Matrices" (2008)
 *   - Thrun, Burgard, Fox "Probabilistic Robotics" (2005)
 *   - Werber et al. "Automotive Radar Gridmap Representations" (2015)
 */

#include "mmwave_tracker.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * DBSCAN clustering
 * ============================================================================
 */

dbscan_result_t *dbscan_cluster(const mmwave_point_t *points, int n,
                                 double eps, int minPts)
{
    if (!points || n <= 0 || eps <= 0.0 || minPts < 1) return NULL;

    dbscan_result_t *result = (dbscan_result_t *)calloc(1, sizeof(dbscan_result_t));
    if (!result) return NULL;

    result->num_points = n;
    result->labels = (int *)malloc(n * sizeof(int));
    if (!result->labels) { free(result); return NULL; }

    /* Initialize all points as unvisited (label = -1) */
    for (int i = 0; i < n; i++) result->labels[i] = -1;

    int cluster_id = 0;

    for (int p = 0; p < n; p++) {
        if (result->labels[p] != -1) continue; /* already processed */

        /* Find neighbors within eps */
        int *neighbors = (int *)malloc(n * sizeof(int));
        int nn = 0;

        for (int q = 0; q < n; q++) {
            double dx = points[p].x - points[q].x;
            double dy = points[p].y - points[q].y;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist <= eps) {
                neighbors[nn++] = q;
            }
        }

        if (nn < minPts) {
            /* Noise point (may be reassigned later) */
            result->labels[p] = 0; /* noise = 0 */
            free(neighbors);
        } else {
            /* New cluster */
            cluster_id++;
            result->labels[p] = cluster_id;

            /* Expand cluster */
            int head = 0;
            while (head < nn) {
                int q = neighbors[head];
                head++;

                if (result->labels[q] == 0) {
                    /* Reassign noise to border point */
                    result->labels[q] = cluster_id;
                }
                if (result->labels[q] != -1) continue;

                result->labels[q] = cluster_id;

                /* Find neighbors of q */
                int *neighbors_q = (int *)malloc(n * sizeof(int));
                int nq = 0;
                for (int r = 0; r < n; r++) {
                    double dx = points[q].x - points[r].x;
                    double dy = points[q].y - points[r].y;
                    double dist = sqrt(dx * dx + dy * dy);
                    if (dist <= eps) {
                        neighbors_q[nq++] = r;
                    }
                }

                if (nq >= minPts) {
                    /* Add new neighbors to the list */
                    int *new_neighbors = (int *)realloc(neighbors,
                                                         (nn + nq) * sizeof(int));
                    if (new_neighbors) {
                        neighbors = new_neighbors;
                        for (int k = 0; k < nq; k++) {
                            int already = 0;
                            for (int m = 0; m < nn; m++) {
                                if (neighbors[m] == neighbors_q[k]) {
                                    already = 1;
                                    break;
                                }
                            }
                            if (!already) {
                                neighbors[nn++] = neighbors_q[k];
                            }
                        }
                    }
                }
                free(neighbors_q);
            }
            free(neighbors);
        }
    }

    result->num_clusters = cluster_id;
    return result;
}

void dbscan_free_result(dbscan_result_t *result)
{
    if (!result) return;
    free(result->labels);
    free(result);
}

dbscan_result_t *dbscan_cluster_adaptive(const mmwave_point_t *points, int n,
                                           double eps_min, double eps_max,
                                           double r_max, int minPts)
{
    if (!points || n <= 0) return NULL;

    /* For adaptive DBSCAN, we use a simple approach:
     * partition points into range bins, run DBSCAN per bin with range-dependent eps,
     * then merge clusters across bins.
     */
    dbscan_result_t *result = (dbscan_result_t *)calloc(1, sizeof(dbscan_result_t));
    if (!result) return NULL;

    result->num_points = n;
    result->labels = (int *)malloc(n * sizeof(int));
    if (!result->labels) { free(result); return NULL; }

    for (int i = 0; i < n; i++) result->labels[i] = 0; /* noise */

    int cluster_id = 0;

    /* Process per range bin */
    const int num_bins = 5;
    for (int bin = 0; bin < num_bins; bin++) {
        double r_low = (bin / (double)num_bins) * r_max;
        double r_high = ((bin + 1) / (double)num_bins) * r_max;
        double r_mid = (r_low + r_high) / 2.0;

        /* Adaptive eps for this range */
        double eps = eps_min + (r_mid / r_max) * (eps_max - eps_min);
        if (eps < eps_min) eps = eps_min;
        if (eps > eps_max) eps = eps_max;

        /* Collect points in this range bin */
        int *bin_points = (int *)malloc(n * sizeof(int));
        int nb = 0;
        for (int i = 0; i < n; i++) {
            double r = sqrt(points[i].x * points[i].x + points[i].y * points[i].y);
            if (r >= r_low && r < r_high && result->labels[i] == 0) {
                bin_points[nb++] = i;
            }
        }

        if (nb < minPts) { free(bin_points); continue; }

        /* Build temporary point array for this bin */
        mmwave_point_t *temp = (mmwave_point_t *)malloc(nb * sizeof(mmwave_point_t));
        for (int k = 0; k < nb; k++) temp[k] = points[bin_points[k]];

        dbscan_result_t *bin_result = dbscan_cluster(temp, nb, eps, minPts);
        if (bin_result && bin_result->num_clusters > 0) {
            for (int k = 0; k < nb; k++) {
                if (bin_result->labels[k] > 0) {
                    int orig_idx = bin_points[k];
                    result->labels[orig_idx] = cluster_id + bin_result->labels[k];
                }
            }
            cluster_id += bin_result->num_clusters;
            dbscan_free_result(bin_result);
        } else {
            dbscan_free_result(bin_result);
        }

        free(bin_points);
        free(temp);
    }

    result->num_clusters = cluster_id;
    return result;
}

/* ============================================================================
 * Object classification
 * ============================================================================
 */

object_class_t classify_radar_object(const mmwave_point_t *points, int n_points,
                                      const int *labels,
                                      double centroid_x, double centroid_y)
{
    if (!points || n_points < 1 || !labels) return OBJ_UNKNOWN;

    (void)centroid_x;
    (void)centroid_y;

    /* Compute features */
    double sum_rcs_db = 0.0, sum_speed = 0.0;
    double max_speed = 0.0, min_speed = 1e10;
    int count = 0;

    for (int i = 0; i < n_points; i++) {
        if (labels[i] <= 0) continue; /* not in class */
        sum_rcs_db += points[i].snr;
        double speed = fabs(points[i].velocity);
        sum_speed += speed;
        if (speed > max_speed) max_speed = speed;
        if (speed < min_speed) min_speed = speed;
        count++;
    }

    if (count == 0) return OBJ_UNKNOWN;

    double mean_rcs = sum_rcs_db / count;
    double mean_speed = sum_speed / count;
    double speed_spread = (count > 1) ? max_speed - min_speed : 0.0;

    /* Rule-based classifier */
    /* Low point count, low RCS, modest speed → pedestrian */
    if (count < 10 && mean_rcs < -5.0 && mean_speed < 2.0) {
        return OBJ_PEDESTRIAN;
    }
    /* Very low speed, spread out → static */
    if (mean_speed < 0.5 && count >= 3) {
        return OBJ_STATIC;
    }
    /* High point count, high RCS, low speed → truck */
    if (count > 30 && mean_rcs > 0.0 && mean_speed < 25.0) {
        return OBJ_TRUCK;
    }
    /* Medium count, medium speed → car */
    if (count >= 10 && mean_speed > 1.0 && mean_speed < 50.0) {
        return OBJ_CAR;
    }
    /* Modest count, low speed, medium RCS → bicycle */
    if (count >= 3 && count < 15 && mean_speed < 10.0 && speed_spread > 1.0) {
        return OBJ_BICYCLE;
    }

    return OBJ_UNKNOWN;
}

void compute_object_shape(const mmwave_point_t *points, int n_points,
                           double centroid_x, double centroid_y,
                           double *length, double *width, double *orientation)
{
    if (!points || !length || !width || !orientation || n_points < 2) {
        if (length) *length = 1.0;
        if (width) *width = 1.0;
        if (orientation) *orientation = 0.0;
        return;
    }

    /* Compute spatial covariance matrix */
    double C_xx = 0.0, C_yy = 0.0, C_xy = 0.0;

    for (int i = 0; i < n_points; i++) {
        double dx = points[i].x - centroid_x;
        double dy = points[i].y - centroid_y;
        C_xx += dx * dx;
        C_yy += dy * dy;
        C_xy += dx * dy;
    }
    C_xx /= n_points;
    C_yy /= n_points;
    C_xy /= n_points;

    /* Eigenvalues of 2×2 matrix [C_xx, C_xy; C_xy, C_yy] */
    double trace = C_xx + C_yy;
    double det = C_xx * C_yy - C_xy * C_xy;
    double disc = sqrt(trace * trace - 4.0 * det);

    double lambda_max = 0.5 * (trace + disc);
    double lambda_min = 0.5 * (trace - disc);
    if (lambda_min < 0.0) lambda_min = 0.0;
    if (lambda_max < lambda_min) { double tmp = lambda_max; lambda_max = lambda_min; lambda_min = tmp; }

    /* Orientation: angle of eigenvector for largest eigenvalue */
    *orientation = 0.5 * atan2(2.0 * C_xy, C_xx - C_yy);

    /* Semi-axes (2-sigma for ~95% containment) */
    *length = 2.0 * sqrt(lambda_max);
    *width  = 2.0 * sqrt(lambda_min);

    if (*length < *width) {
        double tmp = *length; *length = *width; *width = tmp;
        *orientation += M_PI / 2.0;
        while (*orientation > M_PI) *orientation -= M_PI;
    }
}

/* ============================================================================
 * Ego-motion compensation
 * ============================================================================
 */

void ego_motion_compensate(mmwave_point_t *points, int n,
                            double v_ego_x, double v_ego_y,
                            double yaw_rate_ego,
                            double sensor_x, double sensor_y)
{
    if (!points || n <= 0) return;

    for (int i = 0; i < n; i++) {
        /* Position of point relative to vehicle center */
        double px = points[i].x + sensor_x;
        double py = points[i].y + sensor_y;

        /* Ego velocity at this point: v_ego + ω × r */
        double v_ego_at_point_x = v_ego_x - yaw_rate_ego * py;
        double v_ego_at_point_y = v_ego_y + yaw_rate_ego * px;

        /* Project ego velocity onto radial direction to the point */
        double dist = sqrt(px * px + py * py);
        if (dist < 1e-6) continue;
        double ux = px / dist;
        double uy = py / dist;

        double v_ego_radial = v_ego_at_point_x * ux + v_ego_at_point_y * uy;

        /* Compensate: remove ego-motion contribution */
        points[i].velocity -= (float)v_ego_radial;
    }
}

int detect_moving_points(const mmwave_point_t *points, int n, double threshold)
{
    if (!points || n <= 0) return 0;

    int moving = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(points[i].velocity) > threshold) {
            moving++;
        }
    }
    return moving;
}

/* ============================================================================
 * Occupancy grid mapping
 * ============================================================================
 */

occupancy_grid_t *ogm_create(int width, int height, double resolution,
                              double origin_x, double origin_y)
{
    if (width <= 0 || height <= 0 || resolution <= 0.0) return NULL;

    occupancy_grid_t *grid = (occupancy_grid_t *)calloc(1, sizeof(occupancy_grid_t));
    if (!grid) return NULL;

    grid->width = width;
    grid->height = height;
    grid->resolution = resolution;
    grid->origin_x = origin_x;
    grid->origin_y = origin_y;

    grid->cells = (occupancy_cell_t *)calloc(width * height, sizeof(occupancy_cell_t));
    if (!grid->cells) { free(grid); return NULL; }

    grid->prob_prior = 0.5;
    grid->prob_occupied_given_hit = 0.7;
    grid->prob_occupied_given_miss = 0.4;

    /* Initialize all cells to prior */
    double l0 = log(grid->prob_prior / (1.0 - grid->prob_prior));
    for (int i = 0; i < width * height; i++) {
        grid->cells[i].log_odds = l0;
        grid->cells[i].prob_occupied = grid->prob_prior;
    }

    return grid;
}

void ogm_free(occupancy_grid_t *grid)
{
    if (!grid) return;
    free(grid->cells);
    free(grid);
}

void ogm_update_radar_scan(occupancy_grid_t *grid,
                            const mmwave_point_t *points, int n_points,
                            double sensor_x, double sensor_y, double sensor_yaw)
{
    (void)sensor_yaw;
    if (!grid || !points || n_points <= 0) return;

    double l_hit = log(grid->prob_occupied_given_hit
                       / (1.0 - grid->prob_occupied_given_hit));
    double l_miss = log(grid->prob_occupied_given_miss
                        / (1.0 - grid->prob_occupied_given_miss));
    double l_prior = log(grid->prob_prior / (1.0 - grid->prob_prior));

    for (int i = 0; i < n_points; i++) {
        /* Point in global frame (already in sensor frame, add sensor offset) */
        double gx = points[i].x + sensor_x;
        double gy = points[i].y + sensor_y;

        /* Convert to grid cell */
        int cx = (int)((gx - grid->origin_x) / grid->resolution);
        int cy = (int)((gy - grid->origin_y) / grid->resolution);

        if (cx < 0 || cx >= grid->width || cy < 0 || cy >= grid->height)
            continue;

        /* Ray casting: cells from sensor to detection (free) */
        double sx = sensor_x;
        double sy = sensor_y;

        /* Bresenham-like ray casting */
        double dx = gx - sx;
        double dy = gy - sy;
        double dist = sqrt(dx * dx + dy * dy);
        if (dist < 1e-6) continue;

        int steps = (int)(dist / grid->resolution) + 1;
        for (int s = 0; s < steps; s++) {
            double t = (double)s / steps;
            double rx = sx + t * dx;
            double ry = sy + t * dy;
            int rcx = (int)((rx - grid->origin_x) / grid->resolution);
            int rcy = (int)((ry - grid->origin_y) / grid->resolution);
            if (rcx >= 0 && rcx < grid->width && rcy >= 0 && rcy < grid->height) {
                int idx = rcy * grid->width + rcx;
                grid->cells[idx].log_odds += l_miss - l_prior;
                grid->cells[idx].log_odds = fmax(-20.0, fmin(20.0, grid->cells[idx].log_odds));
                grid->cells[idx].num_misses++;
            }
        }

        /* Detection cell (occupied) */
        int idx = cy * grid->width + cx;
        grid->cells[idx].log_odds += l_hit - l_prior;
        grid->cells[idx].log_odds = fmax(-20.0, fmin(20.0, grid->cells[idx].log_odds));
        grid->cells[idx].num_hits++;
        grid->cells[idx].mean_velocity = 0.9 * grid->cells[idx].mean_velocity
                                         + 0.1 * points[i].velocity;
    }

    /* Update probabilities from log-odds */
    for (int i = 0; i < grid->width * grid->height; i++) {
        double l = grid->cells[i].log_odds;
        grid->cells[i].prob_occupied = 1.0 / (1.0 + exp(-l));
    }
}

double ogm_query(const occupancy_grid_t *grid, double x, double y)
{
    if (!grid) return -1.0;

    int cx = (int)((x - grid->origin_x) / grid->resolution);
    int cy = (int)((y - grid->origin_y) / grid->resolution);

    if (cx < 0 || cx >= grid->width || cy < 0 || cy >= grid->height)
        return -1.0;

    return grid->cells[cy * grid->width + cx].prob_occupied;
}

int ogm_extract_obstacles(const occupancy_grid_t *grid, double threshold,
                           double *obstacles_x, double *obstacles_y,
                           int max_obs)
{
    if (!grid || !obstacles_x || !obstacles_y) return 0;

    int count = 0;
    for (int cy = 0; cy < grid->height && count < max_obs; cy++) {
        for (int cx = 0; cx < grid->width && count < max_obs; cx++) {
            int idx = cy * grid->width + cx;
            if (grid->cells[idx].prob_occupied > threshold) {
                obstacles_x[count] = grid->origin_x
                                     + (cx + 0.5) * grid->resolution;
                obstacles_y[count] = grid->origin_y
                                     + (cy + 0.5) * grid->resolution;
                count++;
            }
        }
    }
    return count;
}

/* ============================================================================
 * Extended object tracking (random matrix)
 * ============================================================================
 */

void extended_object_random_matrix(const double *measurements, int n_meas,
                                    int meas_dim, double *x_kin, double *P_kin,
                                    double *V_extent, double *nu,
                                    const double *H, int state_dim)
{
    if (!measurements || !x_kin || !P_kin || !V_extent || !nu || !H)
        return;
    if (n_meas < 2 || meas_dim <= 0) return;

    /* Compute sample mean of measurements */
    double z_bar[TRACK_MAX_MEAS_DIM];
    memset(z_bar, 0, meas_dim * sizeof(double));
    for (int j = 0; j < n_meas; j++) {
        for (int d = 0; d < meas_dim; d++) {
            z_bar[d] += measurements[j * meas_dim + d];
        }
    }
    for (int d = 0; d < meas_dim; d++) z_bar[d] /= n_meas;

    /* Kinematic state update: simple centroid tracking */
    /* x_kin = pseudo-inverse(H) * z_bar */
    double HTH[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(HTH, 0, sizeof(HTH));
    double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_transpose(H, HT, meas_dim, state_dim);
    mat_mat_mul(HT, H, HTH, state_dim, meas_dim, state_dim);

    double H_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_inv_cholesky(HTH, H_inv, state_dim);
    double HT_z[TRACK_MAX_STATE_DIM];
    mat_vec_mul(HT, z_bar, HT_z, state_dim, meas_dim);

    /* Least-squares: x = (H'*H)^{-1} * H' * z_bar */
    double HHT_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_inv_cholesky(HTH, HHT_inv, state_dim);
    mat_vec_mul(HHT_inv, HT_z, x_kin, state_dim, state_dim);

    /* Compute scatter matrix of measurements around centroid */
    memset(V_extent, 0, meas_dim * meas_dim * sizeof(double));
    for (int j = 0; j < n_meas; j++) {
        double dz[TRACK_MAX_MEAS_DIM];
        vec_sub(&measurements[j * meas_dim], z_bar, dz, meas_dim);
        for (int a = 0; a < meas_dim; a++) {
            for (int b = 0; b < meas_dim; b++) {
                V_extent[a * meas_dim + b] += dz[a] * dz[b];
            }
        }
    }

    /* Scale factor: V = V / (n_meas - 1) -> unbiased scatter estimate */
    if (n_meas > 1) {
        for (int i = 0; i < meas_dim * meas_dim; i++) {
            V_extent[i] /= (n_meas - 1);
        }
    }

    /* Degrees of freedom update */
    *nu += n_meas;

    /* Simple kinematic covariance based on measurement spread */
    for (int i = 0; i < state_dim * state_dim; i++) P_kin[i] = 0.0;
    for (int i = 0; i < meas_dim && i < state_dim; i++) {
        for (int j = 0; j < meas_dim && j < state_dim; j++) {
            P_kin[i * state_dim + j] = V_extent[i * meas_dim + j] / n_meas;
        }
    }
}

void extent_to_ellipse(const double *X, double c,
                        double *semi_major, double *semi_minor,
                        double *orientation)
{
    if (!X || !semi_major || !semi_minor || !orientation) return;

    double trace = X[0] + X[3]; /* X is 2x2: [x11, x12; x21, x22] */
    double det = X[0] * X[3] - X[1] * X[2];
    double disc = sqrt(trace * trace - 4.0 * det);
    double lambda1 = 0.5 * (trace + disc);
    double lambda2 = 0.5 * (trace - disc);
    if (lambda2 < 0.0) lambda2 = 0.0;
    if (lambda1 < lambda2) { double t = lambda1; lambda1 = lambda2; lambda2 = t; }

    *semi_major = c * sqrt(lambda1);
    *semi_minor = c * sqrt(lambda2);
    *orientation = 0.5 * atan2(2.0 * X[1], X[0] - X[3]);
}

/* ============================================================================
 * Group tracking
 * ============================================================================
 */

void detect_groups(const double *positions, const double *velocities,
                    int N, double d_max, double v_max,
                    int *groups, int *n_groups)
{
    if (!positions || !velocities || !groups || !n_groups || N <= 0) return;

    for (int i = 0; i < N; i++) groups[i] = -1; /* unassigned */

    /* Build adjacency matrix: A[i][j] = 1 if within distance & velocity threshold */
    int *adj = (int *)calloc(N * N, sizeof(int));
    if (!adj) return;

    for (int i = 0; i < N; i++) {
        adj[i * N + i] = 1; /* self */
        for (int j = i + 1; j < N; j++) {
            double dx = positions[2 * i] - positions[2 * j];
            double dy = positions[2 * i + 1] - positions[2 * j + 1];
            double dist = sqrt(dx * dx + dy * dy);

            double dvx = velocities[2 * i] - velocities[2 * j];
            double dvy = velocities[2 * i + 1] - velocities[2 * j + 1];
            double dvel = sqrt(dvx * dvx + dvy * dvy);

            if (dist <= d_max && dvel <= v_max) {
                adj[i * N + j] = 1;
                adj[j * N + i] = 1;
            }
        }
    }

    /* Find connected components (DFS) */
    int group_id = 0;
    int *visited = (int *)calloc(N, sizeof(int));
    if (!visited) { free(adj); return; }

    for (int i = 0; i < N; i++) {
        if (visited[i] || groups[i] >= 0) continue;

        /* Check if i has any connections */
        int has_connections = 0;
        for (int j = 0; j < N; j++) {
            if (i != j && adj[i * N + j]) { has_connections = 1; break; }
        }

        if (!has_connections) {
            groups[i] = -1; /* solo target */
            visited[i] = 1;
            continue;
        }

        /* DFS from i */
        int *stack = (int *)malloc(N * sizeof(int));
        int stack_top = 0;
        stack[stack_top++] = i;
        groups[i] = group_id;
        visited[i] = 1;

        while (stack_top > 0) {
            int v = stack[--stack_top];
            for (int w = 0; w < N; w++) {
                if (adj[v * N + w] && !visited[w]) {
                    visited[w] = 1;
                    groups[w] = group_id;
                    stack[stack_top++] = w;
                }
            }
        }
        free(stack);
        group_id++;
    }

    *n_groups = group_id;

    free(adj);
    free(visited);
}

void compute_group_bounds(const double *positions, int N,
                           const int *group_ids, int target_group,
                           double *centroid_x, double *centroid_y,
                           double *bbox_xmin, double *bbox_xmax,
                           double *bbox_ymin, double *bbox_ymax)
{
    if (!positions || !group_ids || N <= 0) return;

    double sum_x = 0.0, sum_y = 0.0;
    int count = 0;
    double xmin = 1e10, xmax = -1e10, ymin = 1e10, ymax = -1e10;

    for (int i = 0; i < N; i++) {
        if (group_ids[i] != target_group) continue;
        double x = positions[2 * i];
        double y = positions[2 * i + 1];
        sum_x += x; sum_y += y; count++;
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }

    if (count > 0) {
        if (centroid_x) *centroid_x = sum_x / count;
        if (centroid_y) *centroid_y = sum_y / count;
        if (bbox_xmin) *bbox_xmin = xmin;
        if (bbox_xmax) *bbox_xmax = xmax;
        if (bbox_ymin) *bbox_ymin = ymin;
        if (bbox_ymax) *bbox_ymax = ymax;
    }
}
