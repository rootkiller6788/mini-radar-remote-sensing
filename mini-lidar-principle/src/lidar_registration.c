/**
 * @file    lidar_registration.c
 * @brief   Point cloud registration — SVD transform estimation, ICP
 *
 * Knowledge covered:
 *   L3: SVD decomposition, orthogonal Procrustes problem
 *   L5: Rigid transform estimation from correspondences,
 *       Iterative Closest Point (ICP), nearest-neighbor correspondence,
 *       registration quality metrics
 *   L6: Multi-view registration
 *
 * Reference:
 *   - Besl & McKay, *IEEE TPAMI* 14(2), pp.239-256, 1992.
 *   - Arun, Huang & Blostein, *IEEE TPAMI* 9(5), pp.698-700, 1987.
 *   - Rusinkiewicz & Levoy, *3DIM*, 2001.
 */

#include "lidar_registration.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ─── Internal ──────────────────────────────────────────────────────────── */

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "lidar_registration: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: SVD-based rigid transform estimation (Arun et al. 1987)
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * @brief SVD of a 3x3 matrix (simple 3x3 Jacobi SVD)
 *
 * Computes A = U * S * V^T for a 3x3 matrix.
 * S stored as diagonal elements, U and V as rotation matrices.
 */
static int svd_3x3(const double A[9], double U[9], double S[3], double V[9]) {
    /* Compute A^T * A (3x3 symmetric) */
    double ATA[9];
    for (int r = 0; r < 3; r++) {
        for (int c = r; c < 3; c++) {
            ATA[r * 3 + c] = 0.0;
            for (int k = 0; k < 3; k++) {
                ATA[r * 3 + c] += A[k * 3 + r] * A[k * 3 + c];
            }
            ATA[c * 3 + r] = ATA[r * 3 + c];
        }
    }

    /* Jacobi eigenvalue decomposition of ATA to get V and S^2 */
    double ATA_copy[3][3] = {
        {ATA[0], ATA[1], ATA[2]},
        {ATA[3], ATA[4], ATA[5]},
        {ATA[6], ATA[7], ATA[8]}
    };
    double Vcopy[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

    for (int iter = 0; iter < 20; iter++) {
        int p = 0, q = 1;
        double max_off = fabs(ATA_copy[0][1]);
        if (fabs(ATA_copy[0][2]) > max_off) { max_off = fabs(ATA_copy[0][2]); p = 0; q = 2; }
        if (fabs(ATA_copy[1][2]) > max_off) { max_off = fabs(ATA_copy[1][2]); p = 1; q = 2; }
        if (max_off < 1e-15) break;

        double theta = 0.5 * atan2(2.0 * ATA_copy[p][q],
                                    ATA_copy[q][q] - ATA_copy[p][p]);
        double c = cos(theta), s = sin(theta);

        double App = c*c * ATA_copy[p][p] - 2*c*s * ATA_copy[p][q] + s*s * ATA_copy[q][q];
        double Aqq = s*s * ATA_copy[p][p] + 2*c*s * ATA_copy[p][q] + c*c * ATA_copy[q][q];
        ATA_copy[p][q] = ATA_copy[q][p] = 0.0;

        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                double Arp = c * ATA_copy[r][p] - s * ATA_copy[r][q];
                double Arq = s * ATA_copy[r][p] + c * ATA_copy[r][q];
                ATA_copy[r][p] = ATA_copy[p][r] = Arp;
                ATA_copy[r][q] = ATA_copy[q][r] = Arq;
            }
        }
        ATA_copy[p][p] = App;
        ATA_copy[q][q] = Aqq;

        for (int r = 0; r < 3; r++) {
            double Vrp = c * Vcopy[r][p] - s * Vcopy[r][q];
            double Vrq = s * Vcopy[r][p] + c * Vcopy[r][q];
            Vcopy[r][p] = Vrp;
            Vcopy[r][q] = Vrq;
        }
    }

    /* Eigenvalues (sorted descending) */
    double eig[3] = {ATA_copy[0][0], ATA_copy[1][1], ATA_copy[2][2]};
    int order[3] = {0, 1, 2};
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            if (eig[order[j]] > eig[order[i]]) {
                int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
            }
        }
    }

    /* Fill V with sorted eigenvectors */
    for (int r = 0; r < 3; r++) {
        V[r * 3 + 0] = Vcopy[r][order[0]];
        V[r * 3 + 1] = Vcopy[r][order[1]];
        V[r * 3 + 2] = Vcopy[r][order[2]];
    }

    /* Compute singular values: S = sqrt(|eigenvalues|) */
    for (int i = 0; i < 3; i++) {
        S[i] = sqrt(fabs(eig[order[i]]));
        if (S[i] < 1e-15) S[i] = 0.0;
    }

    /* Compute U = A * V * S^{-1} */
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += A[r * 3 + k] * V[k * 3 + c];
            }
            U[r * 3 + c] = (S[c] > 1e-15) ? sum / S[c] : 0.0;
        }
    }

    return 0;
}

int lidar_estimate_rigid_transform(const lidar_vec3_t *src,
                                     const lidar_vec3_t *dst,
                                     size_t num_pts,
                                     lidar_mat3_t *R,
                                     lidar_vec3_t *t) {
    if (!src || !dst || !R || !t || num_pts < 3) return -1;

    /* 1. Centroids */
    lidar_vec3_t src_c = {0, 0, 0}, dst_c = {0, 0, 0};
    for (size_t i = 0; i < num_pts; i++) {
        src_c.x += src[i].x; src_c.y += src[i].y; src_c.z += src[i].z;
        dst_c.x += dst[i].x; dst_c.y += dst[i].y; dst_c.z += dst[i].z;
    }
    double inv_n = 1.0 / (double)num_pts;
    src_c.x *= inv_n; src_c.y *= inv_n; src_c.z *= inv_n;
    dst_c.x *= inv_n; dst_c.y *= inv_n; dst_c.z *= inv_n;

    /* 2. Cross-covariance matrix H = sum(src'_i * dst'_i^T) */
    double H[9] = {0};
    for (size_t i = 0; i < num_pts; i++) {
        double sx = src[i].x - src_c.x;
        double sy = src[i].y - src_c.y;
        double sz = src[i].z - src_c.z;
        double dx = dst[i].x - dst_c.x;
        double dy = dst[i].y - dst_c.y;
        double dz = dst[i].z - dst_c.z;

        H[0] += sx * dx; H[1] += sx * dy; H[2] += sx * dz;
        H[3] += sy * dx; H[4] += sy * dy; H[5] += sy * dz;
        H[4] += sz * dx; H[7] += sz * dy; H[8] += sz * dz;
    }
    /* Fix the H[4] double-add bug above: correct computation */
    /* Actually let me recompute properly */
    memset(H, 0, 9 * sizeof(double));
    for (size_t i = 0; i < num_pts; i++) {
        double sx = src[i].x - src_c.x;
        double sy = src[i].y - src_c.y;
        double sz = src[i].z - src_c.z;
        double dx = dst[i].x - dst_c.x;
        double dy = dst[i].y - dst_c.y;
        double dz = dst[i].z - dst_c.z;
        H[0] += sx * dx; H[1] += sx * dy; H[2] += sx * dz;
        H[3] += sy * dx; H[4] += sy * dy; H[5] += sy * dz;
        H[6] += sz * dx; H[7] += sz * dy; H[8] += sz * dz;
    }

    /* 3. SVD: H = U * S * V^T */
    double U[9], S[3], VT[9];
    svd_3x3(H, U, S, VT);

    /* 4. Rotation: R = V * U^T, corrected for reflection */
    double Rm[9];
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            /* R = V * U^T */
            Rm[r * 3 + c] = VT[0 * 3 + r] * U[0 * 3 + c]
                          + VT[1 * 3 + r] * U[1 * 3 + c]
                          + VT[2 * 3 + r] * U[2 * 3 + c];
        }
    }

    /* Ensure det(R) = +1 (correct for reflection) */
    double det_R = Rm[0]*(Rm[4]*Rm[8] - Rm[5]*Rm[7])
                 - Rm[1]*(Rm[3]*Rm[8] - Rm[5]*Rm[6])
                 + Rm[2]*(Rm[3]*Rm[7] - Rm[4]*Rm[6]);

    if (det_R < 0.0) {
        /* Flip sign of last column of V (last row of V^T) */
        for (int r = 0; r < 3; r++) {
            Rm[r * 3 + 0] -= 2.0 * VT[2 * 3 + r] * U[2 * 3 + 0];
            Rm[r * 3 + 1] -= 2.0 * VT[2 * 3 + r] * U[2 * 3 + 1];
            Rm[r * 3 + 2] -= 2.0 * VT[2 * 3 + r] * U[2 * 3 + 2];
        }
        /* Recalculate properly:
           Set V[:, 2] = -V[:, 2], then R = V * U^T */
        double V_corrected[9];
        memcpy(V_corrected, VT, 9 * sizeof(double));
        V_corrected[2 * 3 + 0] = -VT[2 * 3 + 0];
        V_corrected[2 * 3 + 1] = -VT[2 * 3 + 1];
        V_corrected[2 * 3 + 2] = -VT[2 * 3 + 2];
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                Rm[r * 3 + c] = V_corrected[0 * 3 + r] * U[0 * 3 + c]
                              + V_corrected[1 * 3 + r] * U[1 * 3 + c]
                              + V_corrected[2 * 3 + r] * U[2 * 3 + c];
            }
        }
    }

    memcpy(R->m, Rm, 9 * sizeof(double));

    /* 5. Translation: t = dst_c - R * src_c */
    t->x = dst_c.x - (Rm[0]*src_c.x + Rm[1]*src_c.y + Rm[2]*src_c.z);
    t->y = dst_c.y - (Rm[3]*src_c.x + Rm[4]*src_c.y + Rm[5]*src_c.z);
    t->z = dst_c.z - (Rm[6]*src_c.x + Rm[7]*src_c.y + Rm[8]*src_c.z);

    return 0;
}

void lidar_apply_transform(lidar_vec3_t *points, size_t num_pts,
                            lidar_mat3_t R, lidar_vec3_t t) {
    if (!points || num_pts == 0) return;
    for (size_t i = 0; i < num_pts; i++) {
        double x = points[i].x, y = points[i].y, z = points[i].z;
        points[i].x = R.m[0]*x + R.m[1]*y + R.m[2]*z + t.x;
        points[i].y = R.m[3]*x + R.m[4]*y + R.m[5]*z + t.y;
        points[i].z = R.m[6]*x + R.m[7]*y + R.m[8]*z + t.z;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: ICP algorithm
 * ═══════════════════════════════════════════════════════════════════════════ */

void lidar_icp_config_default(lidar_icp_config_t *config) {
    if (!config) return;
    config->max_iterations           = 100;
    config->translation_tol         = 1e-6;
    config->rotation_tol            = 1e-6;
    config->max_correspondence_dist = 1.0;
    config->outlier_rejection_std   = 2.5;
    config->use_point_to_plane      = 0;
    config->verbose                 = 0;
}

size_t lidar_brute_force_correspondences(const lidar_vec3_t *src, size_t num_src,
                                           const lidar_vec3_t *tgt, size_t num_tgt,
                                           double max_dist,
                                           int *correspondences,
                                           double *distances) {
    if (!src || !tgt || !correspondences || !distances) return 0;

    size_t n_valid = 0;
    for (size_t i = 0; i < num_src; i++) {
        double best_dist = max_dist * max_dist;
        int best_idx = -1;
        for (size_t j = 0; j < num_tgt; j++) {
            double dx = src[i].x - tgt[j].x;
            double dy = src[i].y - tgt[j].y;
            double dz = src[i].z - tgt[j].z;
            double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best_dist) {
                best_dist = d2;
                best_idx = (int)j;
            }
        }
        correspondences[i] = best_idx;
        if (best_idx >= 0) {
            distances[i] = sqrt(best_dist);
            n_valid++;
        } else {
            distances[i] = max_dist * 2.0;
        }
    }
    return n_valid;
}

int lidar_icp_point_to_point(lidar_scan_t *source,
                               const lidar_scan_t *target,
                               const lidar_icp_config_t *config,
                               lidar_icp_result_t *result) {
    lidar_mat4_t identity = lidar_mat4_identity();
    return lidar_icp_with_initial(source, target, identity, config, result);
}

int lidar_icp_with_initial(lidar_scan_t *source,
                             const lidar_scan_t *target,
                             lidar_mat4_t init_transform,
                             const lidar_icp_config_t *config,
                             lidar_icp_result_t *result) {
    if (!source || !target || !config || !result) return -1;
    if (source->num_points == 0 || target->num_points == 0) return -1;

    memset(result, 0, sizeof(*result));
    result->rmse_history = (double*)safe_malloc(config->max_iterations * sizeof(double));

    /* Extract initial R and t from 4x4 matrix */
    lidar_mat3_t R_accum;
    lidar_vec3_t t_accum;
    for (int i = 0; i < 9; i++) R_accum.m[i] = init_transform.m[(i/3)*4 + (i%3)];
    t_accum.x = init_transform.m[3];
    t_accum.y = init_transform.m[7];
    t_accum.z = init_transform.m[11];

    /* Working copy of source points */
    lidar_vec3_t *src_pts = (lidar_vec3_t*)safe_malloc(source->num_points * sizeof(lidar_vec3_t));
    lidar_vec3_t *tgt_pts = (lidar_vec3_t*)safe_malloc(target->num_points * sizeof(lidar_vec3_t));
    int *correspondences = (int*)malloc(source->num_points * sizeof(int));
    double *distances = (double*)malloc(source->num_points * sizeof(double));

    if (!src_pts || !tgt_pts || !correspondences || !distances) {
        free(src_pts); free(tgt_pts); free(correspondences); free(distances);
        free(result->rmse_history);
        return -1;
    }

    for (size_t i = 0; i < source->num_points; i++) {
        src_pts[i].x = source->points[i].x;
        src_pts[i].y = source->points[i].y;
        src_pts[i].z = source->points[i].z;
    }
    for (size_t i = 0; i < target->num_points; i++) {
        tgt_pts[i].x = target->points[i].x;
        tgt_pts[i].y = target->points[i].y;
        tgt_pts[i].z = target->points[i].z;
    }

    double prev_rmse = 1e100;
    int iter;
    for (iter = 0; iter < config->max_iterations; iter++) {
        /* 1. Find correspondences */
        size_t n_corr = lidar_brute_force_correspondences(
            src_pts, source->num_points,
            tgt_pts, target->num_points,
            config->max_correspondence_dist,
            correspondences, distances);

        if (n_corr < 3) break;

        /* 2. Outlier rejection */
        /* Compute mean and std of valid distances */
        double sum_d = 0.0, sum_d2 = 0.0;
        for (size_t i = 0; i < source->num_points; i++) {
            if (correspondences[i] >= 0) {
                sum_d += distances[i];
                sum_d2 += distances[i] * distances[i];
            }
        }
        if (n_corr == 0) break;
        double mean_d = sum_d / (double)n_corr;
        double std_d = sqrt(sum_d2/(double)n_corr - mean_d*mean_d);
        if (std_d < 1e-15) std_d = 1e-15;

        double max_d = mean_d + config->outlier_rejection_std * std_d;

        /* Collect inlier correspondences */
        lidar_vec3_t *src_inl = (lidar_vec3_t*)malloc(n_corr * sizeof(lidar_vec3_t));
        lidar_vec3_t *tgt_inl = (lidar_vec3_t*)malloc(n_corr * sizeof(lidar_vec3_t));
        size_t n_inl = 0;
        for (size_t i = 0; i < source->num_points; i++) {
            if (correspondences[i] >= 0
                && distances[i] <= max_d
                && distances[i] <= config->max_correspondence_dist) {
                src_inl[n_inl] = src_pts[i];
                tgt_inl[n_inl] = tgt_pts[correspondences[i]];
                n_inl++;
            }
        }

        if (n_inl < 3) { free(src_inl); free(tgt_inl); break; }

        /* 3. Estimate transform */
        lidar_mat3_t R_inc;
        lidar_vec3_t t_inc;
        int est_ok = lidar_estimate_rigid_transform(src_inl, tgt_inl, n_inl,
                                                       &R_inc, &t_inc);
        free(src_inl); free(tgt_inl);

        if (est_ok != 0) break;

        /* 4. Apply incremental transform to source points */
        lidar_apply_transform(src_pts, source->num_points, R_inc, t_inc);

        /* 5. Update accumulated transform: T_i = T_inc * T_{i-1} */
        lidar_mat3_t R_new = lidar_mat3_mul(R_inc, R_accum);
        R_accum = R_new;
        t_accum.x += t_inc.x;
        t_accum.y += t_inc.y;
        t_accum.z += t_inc.z;

        /* 6. Compute RMSE */
        double rmse = 0.0;
        n_corr = lidar_brute_force_correspondences(
            src_pts, source->num_points,
            tgt_pts, target->num_points,
            config->max_correspondence_dist,
            correspondences, distances);
        for (size_t i = 0; i < source->num_points; i++) {
            if (correspondences[i] >= 0) {
                rmse += distances[i] * distances[i];
            }
        }
        if (n_corr > 0) rmse = sqrt(rmse / (double)n_corr);
        else rmse = 1e100;

        result->rmse_history[iter] = rmse;
        result->num_corresp = n_corr;

        if (config->verbose) {
            printf("ICP iter %d: RMSE=%.6f m, corr=%zu\n", iter, rmse, n_corr);
        }

        /* 7. Convergence check */
        double t_change = lidar_vec3_norm(t_inc);
        /* Rotation change: angle of R_inc */
        double trace_R = R_inc.m[0] + R_inc.m[4] + R_inc.m[8];
        double cos_angle = (trace_R - 1.0) / 2.0;
        if (cos_angle > 1.0) cos_angle = 1.0;
        if (cos_angle < -1.0) cos_angle = -1.0;
        double r_change = fabs(acos(cos_angle));

        if (t_change < config->translation_tol && r_change < config->rotation_tol) {
            result->converged = 1;
            iter++;
            break;
        }

        if (fabs(rmse - prev_rmse) < 1e-10) {
            result->converged = 1;
            iter++;
            break;
        }
        prev_rmse = rmse;
    }

    /* Fill result */
    result->transform = lidar_mat4_transform(R_accum, t_accum);
    result->rmse = prev_rmse;
    result->num_iterations = iter;

    /* Update the source point cloud in-place */
    for (size_t i = 0; i < source->num_points; i++) {
        lidar_vec3_t p = { source->points[i].x, source->points[i].y, source->points[i].z };
        p = lidar_mat4_transform_point(result->transform, p);
        source->points[i].x = p.x;
        source->points[i].y = p.y;
        source->points[i].z = p.z;
    }

    free(src_pts); free(tgt_pts); free(correspondences); free(distances);
    return 0;
}

void lidar_icp_result_free(lidar_icp_result_t *result) {
    if (!result) return;
    free(result->rmse_history);
    memset(result, 0, sizeof(*result));
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Registration fitness evaluation
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_registration_rmse(const lidar_scan_t *source,
                                 const lidar_scan_t *target,
                                 double max_dist) {
    if (!source || !target || source->num_points == 0 || target->num_points == 0)
        return -1.0;

    /* Extract point arrays */
    lidar_vec3_t *src_pts = (lidar_vec3_t*)malloc(source->num_points * sizeof(lidar_vec3_t));
    lidar_vec3_t *tgt_pts = (lidar_vec3_t*)malloc(target->num_points * sizeof(lidar_vec3_t));
    int *corr = (int*)malloc(source->num_points * sizeof(int));
    double *dist = (double*)malloc(source->num_points * sizeof(double));

    if (!src_pts || !tgt_pts || !corr || !dist) {
        free(src_pts); free(tgt_pts); free(corr); free(dist);
        return -1.0;
    }

    for (size_t i = 0; i < source->num_points; i++) {
        src_pts[i].x = source->points[i].x;
        src_pts[i].y = source->points[i].y;
        src_pts[i].z = source->points[i].z;
    }
    for (size_t i = 0; i < target->num_points; i++) {
        tgt_pts[i].x = target->points[i].x;
        tgt_pts[i].y = target->points[i].y;
        tgt_pts[i].z = target->points[i].z;
    }

    lidar_brute_force_correspondences(src_pts, source->num_points,
                                                   tgt_pts, target->num_points,
                                                   max_dist, corr, dist);

    double rmse = 0.0;
    size_t n_valid = 0;
    for (size_t i = 0; i < source->num_points; i++) {
        if (corr[i] >= 0) {
            rmse += dist[i] * dist[i];
            n_valid++;
        }
    }

    free(src_pts); free(tgt_pts); free(corr); free(dist);

    if (n_valid == 0) return -1.0;
    return sqrt(rmse / (double)n_valid);
}

double lidar_registration_overlap(const lidar_scan_t *source,
                                    const lidar_scan_t *target,
                                    double max_dist) {
    if (!source || !target || source->num_points == 0) return 0.0;

    lidar_vec3_t *src_pts = (lidar_vec3_t*)malloc(source->num_points * sizeof(lidar_vec3_t));
    lidar_vec3_t *tgt_pts = (lidar_vec3_t*)malloc(target->num_points * sizeof(lidar_vec3_t));
    int *corr = (int*)malloc(source->num_points * sizeof(int));
    double *dist = (double*)malloc(source->num_points * sizeof(double));

    if (!src_pts || !tgt_pts || !corr || !dist) {
        free(src_pts); free(tgt_pts); free(corr); free(dist);
        return 0.0;
    }

    for (size_t i = 0; i < source->num_points; i++) {
        src_pts[i].x = source->points[i].x;
        src_pts[i].y = source->points[i].y;
        src_pts[i].z = source->points[i].z;
    }
    for (size_t i = 0; i < target->num_points; i++) {
        tgt_pts[i].x = target->points[i].x;
        tgt_pts[i].y = target->points[i].y;
        tgt_pts[i].z = target->points[i].z;
    }

    lidar_brute_force_correspondences(src_pts, source->num_points,
                                        tgt_pts, target->num_points,
                                        max_dist, corr, dist);

    size_t n_overlap = 0;
    for (size_t i = 0; i < source->num_points; i++) {
        if (corr[i] >= 0) n_overlap++;
    }

    free(src_pts); free(tgt_pts); free(corr); free(dist);

    return (double)n_overlap / (double)source->num_points;
}