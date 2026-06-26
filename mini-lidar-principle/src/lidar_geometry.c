/**
 * @file    lidar_geometry.c
 * @brief   Point cloud geometry — vectors, matrices, transforms, normals, filtering
 *
 * Knowledge covered:
 *   L1: 3D vector operations, coordinate transforms
 *   L3: Matrix operations (3x3, 4x4), rotation representations,
 *       PCA, singular value decomposition, covariance
 *   L5: Normal estimation (PCA-based), RANSAC plane fitting,
 *       voxel grid downsampling, brute-force KNN
 *   L6: Ground filtering (height threshold, progressive morphological)
 *
 * Reference:
 *   - Hoppe et al., *SIGGRAPH*, 1992.
 *   - Fischler & Bolles, *Comm. ACM* 24(6), pp.381-395, 1981.
 *   - Zhang et al., *IEEE TGRS* 41(4), pp.872-882, 2003.
 */

#include "lidar_geometry.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "lidar_geometry: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L3: Vector operations
 * ═══════════════════════════════════════════════════════════════════════════ */

lidar_vec3_t lidar_vec3_add(lidar_vec3_t a, lidar_vec3_t b) {
    lidar_vec3_t result = { a.x + b.x, a.y + b.y, a.z + b.z };
    return result;
}

lidar_vec3_t lidar_vec3_sub(lidar_vec3_t a, lidar_vec3_t b) {
    lidar_vec3_t result = { a.x - b.x, a.y - b.y, a.z - b.z };
    return result;
}

lidar_vec3_t lidar_vec3_scale(lidar_vec3_t v, double s) {
    lidar_vec3_t result = { v.x * s, v.y * s, v.z * s };
    return result;
}

double lidar_vec3_dot(lidar_vec3_t a, lidar_vec3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

lidar_vec3_t lidar_vec3_cross(lidar_vec3_t a, lidar_vec3_t b) {
    lidar_vec3_t result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

double lidar_vec3_norm(lidar_vec3_t v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

lidar_vec3_t lidar_vec3_normalize(lidar_vec3_t v) {
    double n = lidar_vec3_norm(v);
    if (n < 1e-30) {
        lidar_vec3_t zero = {0.0, 0.0, 0.0};
        return zero;
    }
    return lidar_vec3_scale(v, 1.0 / n);
}

double lidar_vec3_distance(lidar_vec3_t a, lidar_vec3_t b) {
    return lidar_vec3_norm(lidar_vec3_sub(a, b));
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L3: Matrix operations
 * ═══════════════════════════════════════════════════════════════════════════ */

lidar_mat3_t lidar_mat3_identity(void) {
    lidar_mat3_t M;
    memset(M.m, 0, sizeof(M.m));
    M.m[0] = 1.0; M.m[4] = 1.0; M.m[8] = 1.0;
    return M;
}

lidar_mat3_t lidar_mat3_mul(lidar_mat3_t a, lidar_mat3_t b) {
    lidar_mat3_t C;
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            C.m[row * 3 + col] = a.m[row * 3 + 0] * b.m[0 * 3 + col]
                               + a.m[row * 3 + 1] * b.m[1 * 3 + col]
                               + a.m[row * 3 + 2] * b.m[2 * 3 + col];
        }
    }
    return C;
}

lidar_vec3_t lidar_mat3_vec_mul(lidar_mat3_t a, lidar_vec3_t x) {
    lidar_vec3_t y;
    y.x = a.m[0] * x.x + a.m[1] * x.y + a.m[2] * x.z;
    y.y = a.m[3] * x.x + a.m[4] * x.y + a.m[5] * x.z;
    y.z = a.m[6] * x.x + a.m[7] * x.y + a.m[8] * x.z;
    return y;
}

lidar_mat3_t lidar_mat3_transpose(lidar_mat3_t a) {
    lidar_mat3_t T;
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            T.m[r * 3 + c] = a.m[c * 3 + r];
        }
    }
    return T;
}

double lidar_mat3_det(lidar_mat3_t a) {
    return a.m[0] * (a.m[4] * a.m[8] - a.m[5] * a.m[7])
         - a.m[1] * (a.m[3] * a.m[8] - a.m[5] * a.m[6])
         + a.m[2] * (a.m[3] * a.m[7] - a.m[4] * a.m[6]);
}

lidar_mat3_t lidar_mat3_inverse(lidar_mat3_t a) {
    double det = lidar_mat3_det(a);
    lidar_mat3_t inv;
    if (fabs(det) < 1e-30) {
        return lidar_mat3_identity();
    }
    double inv_det = 1.0 / det;
    inv.m[0] = (a.m[4] * a.m[8] - a.m[5] * a.m[7]) * inv_det;
    inv.m[1] = (a.m[2] * a.m[7] - a.m[1] * a.m[8]) * inv_det;
    inv.m[2] = (a.m[1] * a.m[5] - a.m[2] * a.m[4]) * inv_det;
    inv.m[3] = (a.m[5] * a.m[6] - a.m[3] * a.m[8]) * inv_det;
    inv.m[4] = (a.m[0] * a.m[8] - a.m[2] * a.m[6]) * inv_det;
    inv.m[5] = (a.m[2] * a.m[3] - a.m[0] * a.m[5]) * inv_det;
    inv.m[6] = (a.m[3] * a.m[7] - a.m[4] * a.m[6]) * inv_det;
    inv.m[7] = (a.m[1] * a.m[6] - a.m[0] * a.m[7]) * inv_det;
    inv.m[8] = (a.m[0] * a.m[4] - a.m[1] * a.m[3]) * inv_det;
    return inv;
}

lidar_mat4_t lidar_mat4_identity(void) {
    lidar_mat4_t M;
    memset(M.m, 0, sizeof(M.m));
    M.m[0] = 1.0; M.m[5] = 1.0; M.m[10] = 1.0; M.m[15] = 1.0;
    return M;
}

lidar_mat4_t lidar_mat4_transform(lidar_mat3_t R, lidar_vec3_t t) {
    lidar_mat4_t T;
    memset(T.m, 0, sizeof(T.m));
    T.m[0] = R.m[0]; T.m[1] = R.m[1]; T.m[2] = R.m[2]; T.m[3] = t.x;
    T.m[4] = R.m[3]; T.m[5] = R.m[4]; T.m[6] = R.m[5]; T.m[7] = t.y;
    T.m[8] = R.m[6]; T.m[9] = R.m[7]; T.m[10]= R.m[8]; T.m[11]= t.z;
    T.m[12]= 0.0;    T.m[13]= 0.0;    T.m[14]= 0.0;    T.m[15]= 1.0;
    return T;
}

lidar_vec3_t lidar_mat4_transform_point(lidar_mat4_t T, lidar_vec3_t p) {
    double x = T.m[0]*p.x + T.m[1]*p.y + T.m[2]*p.z + T.m[3];
    double y = T.m[4]*p.x + T.m[5]*p.y + T.m[6]*p.z + T.m[7];
    double z = T.m[8]*p.x + T.m[9]*p.y + T.m[10]*p.z + T.m[11];
    double w = T.m[12]*p.x + T.m[13]*p.y + T.m[14]*p.z + T.m[15];
    if (fabs(w) < 1e-15) w = 1.0;
    lidar_vec3_t result = { x / w, y / w, z / w };
    return result;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L3: Rotation representations
 * ═══════════════════════════════════════════════════════════════════════════ */

lidar_mat3_t lidar_rotation_axis_angle(lidar_vec3_t axis, double angle) {
    double n = lidar_vec3_norm(axis);
    if (n < 1e-15) return lidar_mat3_identity();
    lidar_vec3_t k = lidar_vec3_scale(axis, 1.0 / n);

    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;

    lidar_mat3_t R;
    R.m[0] = c + k.x * k.x * t;
    R.m[1] = k.x * k.y * t - k.z * s;
    R.m[2] = k.x * k.z * t + k.y * s;
    R.m[3] = k.y * k.x * t + k.z * s;
    R.m[4] = c + k.y * k.y * t;
    R.m[5] = k.y * k.z * t - k.x * s;
    R.m[6] = k.z * k.x * t - k.y * s;
    R.m[7] = k.z * k.y * t + k.x * s;
    R.m[8] = c + k.z * k.z * t;
    return R;
}

lidar_mat3_t lidar_rotation_euler(double roll, double pitch, double yaw) {
    double cr = cos(roll),  sr = sin(roll);
    double cp = cos(pitch), sp = sin(pitch);
    double cy = cos(yaw),   sy = sin(yaw);

    /* R = Rz(yaw) * Ry(pitch) * Rx(roll) */
    lidar_mat3_t R;
    R.m[0] = cy * cp;
    R.m[1] = cy * sp * sr - sy * cr;
    R.m[2] = cy * sp * cr + sy * sr;
    R.m[3] = sy * cp;
    R.m[4] = sy * sp * sr + cy * cr;
    R.m[5] = sy * sp * cr - cy * sr;
    R.m[6] = -sp;
    R.m[7] = cp * sr;
    R.m[8] = cp * cr;
    return R;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Coordinate transforms
 * ═══════════════════════════════════════════════════════════════════════════ */

lidar_point_t lidar_spherical_to_cartesian(double azimuth, double elevation,
                                            double range) {
    lidar_point_t pt;
    memset(&pt, 0, sizeof(pt));
    double cos_el = cos(elevation);
    pt.x = range * cos_el * sin(azimuth);
    pt.y = range * cos_el * cos(azimuth);
    pt.z = range * sin(elevation);
    return pt;
}

int lidar_transform_scan(lidar_scan_t *scan,
                          lidar_mat3_t R_boresight, lidar_vec3_t t_boresight,
                          lidar_mat3_t R_platform, lidar_vec3_t t_platform) {
    if (!scan || scan->num_points == 0) return -1;

    for (size_t i = 0; i < scan->num_points; i++) {
        lidar_vec3_t p_scanner = { scan->points[i].x,
                                    scan->points[i].y,
                                    scan->points[i].z };
        /* Scanner frame → platform frame */
        lidar_vec3_t p_platform = lidar_mat3_vec_mul(R_boresight, p_scanner);
        p_platform = lidar_vec3_add(p_platform, t_boresight);
        /* Platform frame → global frame */
        lidar_vec3_t p_global = lidar_mat3_vec_mul(R_platform, p_platform);
        p_global = lidar_vec3_add(p_global, t_platform);

        scan->points[i].x = p_global.x;
        scan->points[i].y = p_global.y;
        scan->points[i].z = p_global.z;
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Normal estimation via PCA
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_knn_brute_force(const lidar_scan_t *scan, size_t point_idx,
                           size_t k, size_t *neighbors) {
    if (!scan || point_idx >= scan->num_points || !neighbors || k == 0) return 0;

    size_t n = scan->num_points;
    lidar_vec3_t query = { scan->points[point_idx].x,
                           scan->points[point_idx].y,
                           scan->points[point_idx].z };

    /* Simple O(N*k) sort of distances — adequate for small clouds */
    double *dists = (double*)malloc(n * sizeof(double));
    size_t *indices = (size_t*)malloc(n * sizeof(size_t));
    if (!dists || !indices) { free(dists); free(indices); return 0; }

    for (size_t i = 0; i < n; i++) {
        lidar_vec3_t pi = { scan->points[i].x, scan->points[i].y, scan->points[i].z };
        dists[i] = lidar_vec3_distance(query, pi);
        indices[i] = i;
    }

    /* Bubble-sort distances (small n, simplicity over speed) */
    for (size_t i = 0; i < n && i < k + 10; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (dists[j] < dists[i]) {
                double tmp_d = dists[i];
                dists[i] = dists[j];
                dists[j] = tmp_d;
                size_t tmp_i = indices[i];
                indices[i] = indices[j];
                indices[j] = tmp_i;
            }
        }
    }

    size_t found = 0;
    for (size_t i = 0; i < n && found < k; i++) {
        if (indices[i] != point_idx && dists[i] < 1e-6) {
            /* exclude self with zero distance, but not others */
            if (dists[i] > 0.0) {
                neighbors[found++] = indices[i];
            }
        }
    }

    /* Fill remaining with actual nearest (including self if needed) */
    size_t i = 0;
    while (found < k && i < n) {
        if (indices[i] != point_idx) {
            neighbors[found++] = indices[i];
        }
        i++;
    }

    free(dists);
    free(indices);
    return (int)found;
}

lidar_normal_t lidar_estimate_normal_knn(const lidar_scan_t *scan,
                                           size_t point_idx, size_t k) {
    lidar_normal_t result;
    memset(&result, 0, sizeof(result));

    if (!scan || point_idx >= scan->num_points || k < 3) return result;

    size_t *neighbors = (size_t*)malloc(k * sizeof(size_t));
    if (!neighbors) return result;

    int n_neighbors = lidar_knn_brute_force(scan, point_idx, k, neighbors);
    if (n_neighbors < 3) { free(neighbors); return result; }
    size_t nn = (size_t)n_neighbors;

    /* Compute centroid */
    lidar_vec3_t centroid = {0, 0, 0};
    for (size_t j = 0; j < nn; j++) {
        centroid.x += scan->points[neighbors[j]].x;
        centroid.y += scan->points[neighbors[j]].y;
        centroid.z += scan->points[neighbors[j]].z;
    }
    centroid.x /= (double)nn;
    centroid.y /= (double)nn;
    centroid.z /= (double)nn;
    result.centroid = centroid;

    /* Compute 3x3 covariance matrix */
    double C[9] = {0};
    for (size_t j = 0; j < nn; j++) {
        double dx = scan->points[neighbors[j]].x - centroid.x;
        double dy = scan->points[neighbors[j]].y - centroid.y;
        double dz = scan->points[neighbors[j]].z - centroid.z;
        C[0] += dx * dx; C[1] += dx * dy; C[2] += dx * dz;
        C[3] += dy * dx; C[4] += dy * dy; C[5] += dy * dz;
        C[6] += dz * dx; C[7] += dz * dy; C[8] += dz * dz;
    }
    for (int i = 0; i < 9; i++) C[i] /= (double)nn;
    for (int i = 0; i < 9; i++) result.covariance.m[i] = C[i];

    /* Jacobi eigenvalue decomposition for 3x3 symmetric matrix */
    double A[3][3] = {{C[0], C[1], C[2]},
                       {C[3], C[4], C[5]},
                       {C[6], C[7], C[8]}};
    double V[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    double lambda[3] = {0};

    /* Jacobi rotation iterations */
    for (int iter = 0; iter < 20; iter++) {
        /* Find largest off-diagonal */
        int p = 0, q = 1;
        double max_off = fabs(A[0][1]);
        if (fabs(A[0][2]) > max_off) { max_off = fabs(A[0][2]); p = 0; q = 2; }
        if (fabs(A[1][2]) > max_off) { max_off = fabs(A[1][2]); p = 1; q = 2; }
        if (max_off < 1e-15) break;

        double theta = 0.5 * atan2(2.0 * A[p][q], A[q][q] - A[p][p]);
        double c = cos(theta), s = sin(theta);

        /* Apply rotation to A */
        double App = c*c*A[p][p] - 2*c*s*A[p][q] + s*s*A[q][q];
        double Aqq = s*s*A[p][p] + 2*c*s*A[p][q] + c*c*A[q][q];
        A[p][q] = A[q][p] = 0.0;

        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                double Arp = c*A[r][p] - s*A[r][q];
                double Arq = s*A[r][p] + c*A[r][q];
                A[r][p] = A[p][r] = Arp;
                A[r][q] = A[q][r] = Arq;
            }
        }
        A[p][p] = App;
        A[q][q] = Aqq;

        /* Accumulate eigenvectors */
        for (int r = 0; r < 3; r++) {
            double Vrp = c*V[r][p] - s*V[r][q];
            double Vrq = s*V[r][p] + c*V[r][q];
            V[r][p] = Vrp;
            V[r][q] = Vrq;
        }
    }

    /* Extract eigenvalues */
    lambda[0] = A[0][0]; lambda[1] = A[1][1]; lambda[2] = A[2][2];

    /* Sort eigenvalues and eigenvectors in descending order */
    int order[3] = {0, 1, 2};
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            if (lambda[order[j]] > lambda[order[i]]) {
                int tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }
        }
    }
    result.eigenvalues[0] = lambda[order[0]];
    result.eigenvalues[1] = lambda[order[1]];
    result.eigenvalues[2] = lambda[order[2]];

    /* Normal = eigenvector corresponding to smallest eigenvalue */
    int min_idx = order[2];
    result.normal.x = V[0][min_idx];
    result.normal.y = V[1][min_idx];
    result.normal.z = V[2][min_idx];

    /* Ensure normal points toward sensor (positive Z direction by convention) */
    if (result.normal.z < 0.0) {
        result.normal.x = -result.normal.x;
        result.normal.y = -result.normal.y;
        result.normal.z = -result.normal.z;
    }

    /* Curvature = lambda_min / (lambda_1 + lambda_2 + lambda_3) */
    double sum_lambda = result.eigenvalues[0] + result.eigenvalues[1] + result.eigenvalues[2];
    if (sum_lambda > 1e-15) {
        result.curvature = result.eigenvalues[2] / sum_lambda;
    }
    result.valid = 1;

    free(neighbors);
    return result;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: RANSAC plane fitting
 * ═══════════════════════════════════════════════════════════════════════════ */

lidar_ransac_plane_t lidar_ransac_plane(const lidar_scan_t *scan,
                                          double distance_threshold,
                                          int max_iter,
                                          double min_inlier_ratio) {
    lidar_ransac_plane_t result;
    memset(&result, 0, sizeof(result));

    if (!scan || scan->num_points < 3 || max_iter < 1) return result;

    size_t N = scan->num_points;
    size_t best_n_inliers = 0;
    lidar_vec3_t best_normal = {0, 0, 1};
    double best_d = 0.0;

    /* RANSAC iterations */
    for (int iter = 0; iter < max_iter; iter++) {
        /* Randomly select 3 points */
        size_t i1 = (size_t)rand() % N;
        size_t i2 = (size_t)rand() % N;
        size_t i3 = (size_t)rand() % N;
        if (i1 == i2 || i2 == i3 || i1 == i3) continue;

        lidar_vec3_t p1 = { scan->points[i1].x, scan->points[i1].y, scan->points[i1].z };
        lidar_vec3_t p2 = { scan->points[i2].x, scan->points[i2].y, scan->points[i2].z };
        lidar_vec3_t p3 = { scan->points[i3].x, scan->points[i3].y, scan->points[i3].z };

        /* Plane normal = (p2 - p1) x (p3 - p1) */
        lidar_vec3_t v1 = lidar_vec3_sub(p2, p1);
        lidar_vec3_t v2 = lidar_vec3_sub(p3, p1);
        lidar_vec3_t n = lidar_vec3_cross(v1, v2);
        double n_norm = lidar_vec3_norm(n);
        if (n_norm < 1e-15) continue; /* Collinear points */
        n = lidar_vec3_scale(n, 1.0 / n_norm);

        /* Plane equation: n·x + d = 0 => d = -n·p1 */
        double d = -lidar_vec3_dot(n, p1);

        /* Count inliers */
        size_t n_inliers = 0;
        for (size_t j = 0; j < N; j++) {
            lidar_vec3_t pj = { scan->points[j].x, scan->points[j].y, scan->points[j].z };
            double dist = fabs(lidar_vec3_dot(n, pj) + d);
            if (dist < distance_threshold) n_inliers++;
        }

        if (n_inliers > best_n_inliers) {
            best_n_inliers = n_inliers;
            best_normal = n;
            best_d = d;
        }

        double ratio = (double)best_n_inliers / (double)N;
        if (ratio > 0.95) break; /* Early exit for clear planes */
    }

    double ratio = (double)best_n_inliers / (double)N;
    if (ratio < min_inlier_ratio || best_n_inliers < 3) return result;

    result.normal = best_normal;
    result.d = best_d;
    result.num_inliers = best_n_inliers;
    result.inlier_ratio = ratio;
    result.converged = 1;
    result.inliers = (size_t*)safe_malloc(best_n_inliers * sizeof(size_t));

    /* Fill inlier indices */
    size_t idx = 0;
    for (size_t j = 0; j < N; j++) {
        lidar_vec3_t pj = { scan->points[j].x, scan->points[j].y, scan->points[j].z };
        double dist = fabs(lidar_vec3_dot(best_normal, pj) + best_d);
        if (dist < distance_threshold) {
            result.inliers[idx++] = j;
        }
    }

    return result;
}

void lidar_ransac_plane_free(lidar_ransac_plane_t *result) {
    if (!result) return;
    free(result->inliers);
    memset(result, 0, sizeof(*result));
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Voxel grid downsampling
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_voxel_downsample(const lidar_scan_t *input,
                            double voxel_size,
                            lidar_scan_t *output) {
    if (!input || !output || voxel_size <= 0.0 || input->num_points == 0) return -1;

    /* Compute grid bounds */
    double min_xyz[3], max_xyz[3];
    if (lidar_scan_bounding_box(input, min_xyz, max_xyz) != 0) return -1;

    size_t nx = (size_t)ceil((max_xyz[0] - min_xyz[0]) / voxel_size) + 1;
    size_t ny = (size_t)ceil((max_xyz[1] - min_xyz[1]) / voxel_size) + 1;
    size_t nz = (size_t)ceil((max_xyz[2] - min_xyz[2]) / voxel_size) + 1;
    if (nx == 0) nx = 1;
    if (ny == 0) ny = 1;
    if (nz == 0) nz = 1;
    if (nx * ny * nz > 10000000) return -1; /* Too large */

    /* Per-voxel accumulators */
    typedef struct { double sx, sy, sz; size_t count; } voxel_acc_t;
    size_t n_voxels = nx * ny * nz;
    voxel_acc_t *voxels = (voxel_acc_t*)calloc(n_voxels, sizeof(voxel_acc_t));
    if (!voxels) return -1;

    for (size_t i = 0; i < input->num_points; i++) {
        size_t ix = (size_t)((input->points[i].x - min_xyz[0]) / voxel_size);
        size_t iy = (size_t)((input->points[i].y - min_xyz[1]) / voxel_size);
        size_t iz = (size_t)((input->points[i].z - min_xyz[2]) / voxel_size);
        if (ix >= nx) ix = nx - 1;
        if (iy >= ny) iy = ny - 1;
        if (iz >= nz) iz = nz - 1;
        size_t vi = ix + iy * nx + iz * nx * ny;
        voxels[vi].sx += input->points[i].x;
        voxels[vi].sy += input->points[i].y;
        voxels[vi].sz += input->points[i].z;
        voxels[vi].count++;
    }

    /* Output centroid of each populated voxel */
    lidar_scan_init(output, n_voxels);
    for (size_t vi = 0; vi < n_voxels; vi++) {
        if (voxels[vi].count > 0) {
            lidar_point_t pt;
            memset(&pt, 0, sizeof(pt));
            double inv = 1.0 / (double)voxels[vi].count;
            pt.x = voxels[vi].sx * inv;
            pt.y = voxels[vi].sy * inv;
            pt.z = voxels[vi].sz * inv;
            lidar_scan_add_point(output, &pt);
        }
    }

    free(voxels);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L6: Ground filtering
 * ═══════════════════════════════════════════════════════════════════════════ */

size_t lidar_ground_filter_height(const lidar_scan_t *scan, double height_thr) {
    if (!scan || scan->num_points == 0) return 0;

    /* Find minimum Z */
    double z_min = scan->points[0].z;
    for (size_t i = 1; i < scan->num_points; i++) {
        if (scan->points[i].z < z_min) z_min = scan->points[i].z;
    }

    double ground_z = z_min + height_thr;
    size_t count = 0;
    for (size_t i = 0; i < scan->num_points; i++) {
        if (scan->points[i].z <= ground_z) {
            scan->points[i].class_label = 2; /* Ground */
            count++;
        } else {
            scan->points[i].class_label = 1; /* Non-ground */
        }
    }
    return count;
}

/* Progressive Morphological Filter — simplified implementation */
size_t lidar_ground_filter_pmf(lidar_scan_t *scan,
                                 double max_window,
                                 double slope_threshold,
                                 double init_distance) {
    if (!scan || scan->num_points == 0) return 0;

    /* This is a simplified PMF: iterative erosion-dilation with increasing window.
       Full PMF requires gridding the point cloud first.  Here we use a
       cell-based approximation. */

    double min_xyz[3], max_xyz[3];
    if (lidar_scan_bounding_box(scan, min_xyz, max_xyz) != 0) return 0;

    double grid_size = 1.0; /* 1m grid */
    size_t nx = (size_t)((max_xyz[0] - min_xyz[0]) / grid_size) + 2;
    size_t ny = (size_t)((max_xyz[1] - min_xyz[1]) / grid_size) + 2;
    if (nx > 10000 || ny > 10000) return 0;

    /* Minimum Z per grid cell */
    double *grid_zmin = (double*)malloc(nx * ny * sizeof(double));
    if (!grid_zmin) return 0;
    for (size_t i = 0; i < nx * ny; i++) grid_zmin[i] = 1e100;

    for (size_t i = 0; i < scan->num_points; i++) {
        size_t ix = (size_t)((scan->points[i].x - min_xyz[0]) / grid_size);
        size_t iy = (size_t)((scan->points[i].y - min_xyz[1]) / grid_size);
        if (ix < nx && iy < ny) {
            size_t idx = ix + iy * nx;
            if (scan->points[i].z < grid_zmin[idx]) {
                grid_zmin[idx] = scan->points[i].z;
            }
        }
    }

    /* Initialize ground surface = cell minima */
    double *surface = (double*)malloc(nx * ny * sizeof(double));
    memcpy(surface, grid_zmin, nx * ny * sizeof(double));

    /* Iterative PMF */
    double window = grid_size;
    double dh_max = init_distance;
    while (window <= max_window) {
        size_t half_win = (size_t)(window / grid_size / 2.0);
        if (half_win < 1) half_win = 1;

        /* Morphological opening: erosion then dilation */
        double *eroded = (double*)malloc(nx * ny * sizeof(double));

        /* Erosion: minimum in window */
        for (size_t iy = 0; iy < ny; iy++) {
            for (size_t ix = 0; ix < nx; ix++) {
                double min_val = 1e100;
                size_t x0 = (ix > half_win) ? ix - half_win : 0;
                size_t x1 = (ix + half_win < nx) ? ix + half_win : nx - 1;
                size_t y0 = (iy > half_win) ? iy - half_win : 0;
                size_t y1 = (iy + half_win < ny) ? iy + half_win : ny - 1;
                for (size_t wy = y0; wy <= y1; wy++) {
                    for (size_t wx = x0; wx <= x1; wx++) {
                        double v = surface[wx + wy * nx];
                        if (v < min_val) min_val = v;
                    }
                }
                eroded[ix + iy * nx] = (min_val < 1e99) ? min_val : surface[ix + iy * nx];
            }
        }

        /* Dilation: maximum in window */
        double *dilated = (double*)malloc(nx * ny * sizeof(double));
        for (size_t iy = 0; iy < ny; iy++) {
            for (size_t ix = 0; ix < nx; ix++) {
                double max_val = -1e100;
                size_t x0 = (ix > half_win) ? ix - half_win : 0;
                size_t x1 = (ix + half_win < nx) ? ix + half_win : nx - 1;
                size_t y0 = (iy > half_win) ? iy - half_win : 0;
                size_t y1 = (iy + half_win < ny) ? iy + half_win : ny - 1;
                for (size_t wy = y0; wy <= y1; wy++) {
                    for (size_t wx = x0; wx <= x1; wx++) {
                        double v = eroded[wx + wy * nx];
                        if (v > max_val) max_val = v;
                    }
                }
                dilated[ix + iy * nx] = (max_val > -1e99) ? max_val : surface[ix + iy * nx];
            }
        }

        /* Update surface: keep points close to the opened surface */
        for (size_t iy = 0; iy < ny; iy++) {
            for (size_t ix = 0; ix < nx; ix++) {
                size_t idx = ix + iy * nx;
                double diff = surface[idx] - dilated[idx];
                if (diff <= dh_max) {
                    surface[idx] = dilated[idx];
                }
            }
        }

        free(eroded);
        free(dilated);

        /* Increase window size and dh_max for next iteration */
        window *= 2.0;
        dh_max += slope_threshold * window * 0.5;
    }

    /* Classify points */
    size_t count = 0;
    for (size_t i = 0; i < scan->num_points; i++) {
        size_t ix = (size_t)((scan->points[i].x - min_xyz[0]) / grid_size);
        size_t iy = (size_t)((scan->points[i].y - min_xyz[1]) / grid_size);
        int is_ground = 0;
        if (ix < nx && iy < ny) {
            double surf_z = surface[ix + iy * nx];
            if (scan->points[i].z <= surf_z + init_distance) {
                is_ground = 1;
            }
        } else {
            /* Points outside grid: use height-based fallback */
            if (scan->points[i].z <= min_xyz[2] + init_distance) is_ground = 1;
        }
        scan->points[i].class_label = is_ground ? 2 : 1;
        if (is_ground) count++;
    }

    free(grid_zmin);
    free(surface);
    return count;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L3: Point cloud statistics
 * ═══════════════════════════════════════════════════════════════════════════ */

lidar_vec3_t lidar_scan_centroid(const lidar_scan_t *scan) {
    lidar_vec3_t c = {0, 0, 0};
    if (!scan || scan->num_points == 0) return c;
    for (size_t i = 0; i < scan->num_points; i++) {
        c.x += scan->points[i].x;
        c.y += scan->points[i].y;
        c.z += scan->points[i].z;
    }
    double inv = 1.0 / (double)scan->num_points;
    c.x *= inv; c.y *= inv; c.z *= inv;
    return c;
}

double lidar_scan_density(const lidar_scan_t *scan, double z_min, double z_max) {
    if (!scan || scan->num_points == 0) return 0.0;
    double bmin[3], bmax[3];
    if (lidar_scan_bounding_box(scan, bmin, bmax) != 0) return 0.0;
    double area = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]);
    if (area < 1e-10) return 0.0;

    size_t count = 0;
    for (size_t i = 0; i < scan->num_points; i++) {
        double z = scan->points[i].z;
        if (z >= z_min && z <= z_max) count++;
    }
    return (double)count / area;
}

void lidar_scan_elevation_stats(const lidar_scan_t *scan,
                                 double *z_min, double *z_max,
                                 double *z_mean, double *z_std) {
    if (!scan || scan->num_points == 0) {
        if (z_min) *z_min = 0;
        if (z_max) *z_max = 0;
        if (z_mean) *z_mean = 0;
        if (z_std) *z_std = 0;
        return;
    }
    double zmin = scan->points[0].z, zmax = scan->points[0].z;
    double sum = 0.0;
    for (size_t i = 0; i < scan->num_points; i++) {
        double z = scan->points[i].z;
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
        sum += z;
    }
    double mean = sum / (double)scan->num_points;
    double sum_sq = 0.0;
    for (size_t i = 0; i < scan->num_points; i++) {
        double d = scan->points[i].z - mean;
        sum_sq += d * d;
    }
    double std_dev = sqrt(sum_sq / (double)scan->num_points);
    if (z_min) *z_min = zmin;
    if (z_max) *z_max = zmax;
    if (z_mean) *z_mean = mean;
    if (z_std) *z_std = std_dev;
}