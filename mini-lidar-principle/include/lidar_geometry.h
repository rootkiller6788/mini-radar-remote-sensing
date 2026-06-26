/**
 * @file    lidar_geometry.h
 * @brief   Point cloud geometry — coordinate transforms, normals, features
 *
 * Knowledge covered:
 *   L1: Coordinate systems (sensor, local, global), point cloud representation
 *   L3: Rigid-body transforms, rotation matrices, quaternions,
 *       PCA for normal estimation, covariance structure
 *   L5: Normal estimation, curvature computation, RANSAC plane fitting,
 *       voxel grid downsampling, KD-tree spatial indexing
 *   L6: Ground filtering, feature extraction
 *
 * Reference:
 *   - Hoppe et al., "Surface Reconstruction from Unorganized Points",
 *     *SIGGRAPH*, 1992 (PCA normal estimation).
 *   - Rusu, R.B., *Semantic 3D Object Maps for Everyday Manipulation*,
 *     PhD Thesis, TU Munich, 2009 (PCL design).
 *   - Fischler & Bolles, "Random Sample Consensus: A Paradigm for Model
 *     Fitting", *Comm. ACM* 24(6), pp.381-395, 1981.
 */

#ifndef LIDAR_GEOMETRY_H
#define LIDAR_GEOMETRY_H

#include "lidar_core.h"
#include <stddef.h>

/* ─── L1: 3D vector operations ──────────────────────────────────────────── */

typedef struct {
    double x, y, z;
} lidar_vec3_t;

typedef struct {
    double m[9];  /**< Row-major 3x3 matrix: m[row*3+col] */
} lidar_mat3_t;

typedef struct {
    double m[16]; /**< Row-major 4x4 homogeneous transform matrix */
} lidar_mat4_t;

/* ─── L3: Vector math ──────────────────────────────────────────────────── */

/** Vector addition: c = a + b */
lidar_vec3_t lidar_vec3_add(lidar_vec3_t a, lidar_vec3_t b);

/** Vector subtraction: c = a - b */
lidar_vec3_t lidar_vec3_sub(lidar_vec3_t a, lidar_vec3_t b);

/** Scalar multiplication: v * s */
lidar_vec3_t lidar_vec3_scale(lidar_vec3_t v, double s);

/** Dot product: a · b */
double lidar_vec3_dot(lidar_vec3_t a, lidar_vec3_t b);

/** Cross product: a × b */
lidar_vec3_t lidar_vec3_cross(lidar_vec3_t a, lidar_vec3_t b);

/** Vector norm (Euclidean length) */
double lidar_vec3_norm(lidar_vec3_t v);

/** Normalize vector to unit length — returns zero vector if input is zero */
lidar_vec3_t lidar_vec3_normalize(lidar_vec3_t v);

/** Distance between two points */
double lidar_vec3_distance(lidar_vec3_t a, lidar_vec3_t b);

/* ─── L3: Matrix operations ─────────────────────────────────────────────── */

/** 3x3 identity matrix */
lidar_mat3_t lidar_mat3_identity(void);

/** 3x3 matrix multiplication: C = A * B */
lidar_mat3_t lidar_mat3_mul(lidar_mat3_t a, lidar_mat3_t b);

/** Matrix-vector multiply: y = A * x */
lidar_vec3_t lidar_mat3_vec_mul(lidar_mat3_t a, lidar_vec3_t x);

/** Matrix transpose */
lidar_mat3_t lidar_mat3_transpose(lidar_mat3_t a);

/** Matrix determinant */
double lidar_mat3_det(lidar_mat3_t a);

/** Matrix inverse — returns identity on singular matrix */
lidar_mat3_t lidar_mat3_inverse(lidar_mat3_t a);

/** 4x4 identity matrix */
lidar_mat4_t lidar_mat4_identity(void);

/** 4x4 homogeneous transform: translate then rotate (R|t form) */
lidar_mat4_t lidar_mat4_transform(lidar_mat3_t R, lidar_vec3_t t);

/** Transform a 3D point by 4x4 homogeneous matrix */
lidar_vec3_t lidar_mat4_transform_point(lidar_mat4_t T, lidar_vec3_t p);

/* ─── L3: Rotation representations ─────────────────────────────────────── */

/**
 * @brief Rotation matrix from axis-angle representation (Rodrigues formula)
 *
 * R = I + sin(θ)·[k]× + (1-cos(θ))·[k]×²
 *
 * where k is the unit rotation axis and θ is the angle.
 *
 * @param axis   Rotation axis (must be non-zero, will be normalized internally)
 * @param angle  Rotation angle [rad]
 * @return       3x3 rotation matrix
 */
lidar_mat3_t lidar_rotation_axis_angle(lidar_vec3_t axis, double angle);

/**
 * @brief Rotation matrix from Euler angles (ZYX convention: yaw-pitch-roll)
 *
 * R = Rz(yaw) * Ry(pitch) * Rx(roll)
 *
 * @param roll   Rotation about X axis [rad]
 * @param pitch  Rotation about Y axis [rad]
 * @param yaw    Rotation about Z axis [rad]
 * @return       3x3 rotation matrix
 */
lidar_mat3_t lidar_rotation_euler(double roll, double pitch, double yaw);

/* ─── L5: Coordinate transforms ─────────────────────────────────────────── */

/**
 * @brief Convert spherical (scan) coordinates to Cartesian
 *
 * Given scanner azimuth φ, elevation θ, and slant range R:
 *   x = R · cos(θ) · sin(φ)
 *   y = R · cos(θ) · cos(φ)
 *   z = R · sin(θ)
 *
 * This is the fundamental transform for any scanning LiDAR.
 * Convention: φ=0 along +Y axis, increasing clockwise from above (North-East-Down).
 *
 * @param azimuth    Scanner azimuth angle [rad]
 * @param elevation  Scanner elevation angle [rad] (0 = horizontal)
 * @param range      Measured slant range [m]
 * @return           Cartesian (x, y, z) in scanner coordinate frame [m]
 */
lidar_point_t lidar_spherical_to_cartesian(double azimuth, double elevation,
                                            double range);

/**
 * @brief Convert point cloud from scanner frame to global frame
 *
 * Applies boresight alignment and platform pose to transform the
 * entire scan into a geo-referenced coordinate system.
 *
 * p_global = R_platform * (R_boresight * p_scanner + t_boresight) + t_platform
 *
 * @param scan          Input scan (modified in-place)
 * @param R_boresight   Boresight rotation matrix (scanner → platform)
 * @param t_boresight   Boresight translation [m]
 * @param R_platform    Platform attitude rotation matrix
 * @param t_platform    Platform position in global frame [m]
 * @return              0 on success
 */
int lidar_transform_scan(lidar_scan_t *scan,
                          lidar_mat3_t R_boresight, lidar_vec3_t t_boresight,
                          lidar_mat3_t R_platform, lidar_vec3_t t_platform);

/* ─── L5: Normal estimation via PCA ─────────────────────────────────────── */

/**
 * @brief Covariance matrix of point neighborhoods for normal estimation
 */
typedef struct {
    lidar_vec3_t centroid;     /**< Mean of neighborhood points */
    lidar_mat3_t covariance;   /**< 3x3 covariance matrix */
    lidar_vec3_t normal;       /**< Estimated surface normal (unit vector) */
    double       curvature;    /**< Surface curvature: λ_min / (λ_1 + λ_2 + λ_3) */
    double       eigenvalues[3]; /**< Sorted eigenvalues: λ_1 ≥ λ_2 ≥ λ_3 */
    int          valid;        /**< 1 if sufficient points for reliable estimate */
} lidar_normal_t;

/**
 * @brief Estimate surface normal using PCA on k-nearest neighbors
 *
 * The normal is the eigenvector corresponding to the smallest eigenvalue
 * of the 3x3 covariance matrix of the k-nearest neighbors.
 *
 *   C = (1/k) · Σ (p_i - p_centroid) · (p_i - p_centroid)^T
 *
 *   C · v_j = λ_j · v_j    (j = 1,2,3; λ_1 ≥ λ_2 ≥ λ_3)
 *   n = ± v_3              (oriented toward sensor for consistency)
 *
 * Curvature estimate:
 *   σ = λ_3 / (λ_1 + λ_2 + λ_3)
 *
 * Reference: Hoppe et al. (1992); Pauly et al., "Efficient Simplification
 *            of Point-Sampled Surfaces", *IEEE Vis*, 2002.
 *
 * @param scan        Input point cloud
 * @param point_idx   Index of query point
 * @param k           Number of nearest neighbors
 * @return            Normal estimation result
 */
lidar_normal_t lidar_estimate_normal_knn(const lidar_scan_t *scan,
                                           size_t point_idx, size_t k);

/**
 * @brief Brute-force k-nearest neighbors search (for small point clouds)
 *
 * Finds the k nearest points to the query point using linear scan.
 * O(N) per query.  For production use, a KD-tree is preferred.
 *
 * @param scan        Input point cloud
 * @param point_idx   Query point index
 * @param k           Number of neighbors
 * @param neighbors   Output: indices of k nearest neighbors
 * @return            Number of neighbors found (may be < k for tiny clouds)
 */
int lidar_knn_brute_force(const lidar_scan_t *scan, size_t point_idx,
                           size_t k, size_t *neighbors);

/* ─── L5: RANSAC plane fitting ──────────────────────────────────────────── */

/**
 * @brief RANSAC plane fitting result
 */
typedef struct {
    lidar_vec3_t normal;       /**< Plane normal (unit vector) */
    double       d;            /**< Plane offset: n·x + d = 0 */
    size_t      *inliers;      /**< Array of inlier point indices */
    size_t       num_inliers;  /**< Number of inliers */
    double       inlier_ratio; /**< inliers / total_points */
    int          converged;    /**< 1 if a valid plane was found */
} lidar_ransac_plane_t;

/**
 * @brief Fit a plane to point cloud using RANSAC
 *
 * Algorithm:
 *   1. Randomly select 3 non-collinear points
 *   2. Compute plane parameters (normal + offset)
 *   3. Count inliers (points within distance_threshold of plane)
 *   4. Repeat for max_iter, keep best model
 *   5. Optional: re-fit plane to all inliers using total least squares
 *
 * Complexity: O(max_iter · N)
 *
 * The plane equation is: n · x + d = 0, where n is the unit normal
 * and d = -n · p_0 for any point p_0 on the plane.
 *
 * Reference: Fischler & Bolles (1981).
 *
 * @param scan                Input point cloud
 * @param distance_threshold  Maximum point-to-plane distance for inlier [m]
 * @param max_iter            Maximum RANSAC iterations
 * @param min_inlier_ratio    Minimum ratio of inliers to accept (e.g., 0.3)
 * @return                    Plane fitting result (caller must free inliers)
 */
lidar_ransac_plane_t lidar_ransac_plane(const lidar_scan_t *scan,
                                          double distance_threshold,
                                          int max_iter,
                                          double min_inlier_ratio);

/**
 * @brief Free memory allocated by RANSAC result
 */
void lidar_ransac_plane_free(lidar_ransac_plane_t *result);

/* ─── L5: Voxel grid downsampling ──────────────────────────────────────── */

/**
 * @brief Voxel grid downsampling (3D box filter)
 *
 * Reduces point cloud density by replacing points within each voxel
 * with their centroid.  Preserves surface structure while reducing
 * data volume for downstream processing (registration, classification).
 *
 * @param input      Input scan
 * @param voxel_size Voxel edge length [m]
 * @param output     Output downsampled scan (must be initialized)
 * @return           0 on success
 */
int lidar_voxel_downsample(const lidar_scan_t *input,
                            double voxel_size,
                            lidar_scan_t *output);

/* ─── L6: Ground filtering ─────────────────────────────────────────────── */

/**
 * @brief Simple ground filtering by height threshold
 *
 * Classifies points below a height threshold (relative to min Z)
 * as ground (class_label = 2).  Points above are non-ground (class_label = 1).
 *
 * This is the simplest form of ground filtering; for production use,
 * progressive morphological filters or cloth simulation are preferred.
 *
 * @param scan       Point cloud (modified in-place: class_label updated)
 * @param height_thr Height above local minimum to classify as ground [m]
 * @return           Number of points classified as ground
 */
size_t lidar_ground_filter_height(const lidar_scan_t *scan, double height_thr);

/**
 * @brief Progressive Morphological Filter (PMF) for ground classification
 *
 * Iteratively applies morphological opening (erosion then dilation) with
 * increasing window sizes, progressively removing non-ground objects
 * (buildings, vegetation) while preserving terrain.
 *
 * Reference: Zhang et al., "A Progressive Morphological Filter for
 *            Removing Nonground Measurements from Airborne LIDAR Data",
 *            *IEEE TGRS* 41(4), pp.872-882, 2003.
 *
 * @param scan            Point cloud (modified in-place)
 * @param max_window      Maximum morphological window size [m]
 * @param slope_threshold Terrain slope threshold for window size increase
 * @param init_distance   Initial elevation difference threshold [m]
 * @return                Number of ground points found
 */
size_t lidar_ground_filter_pmf(lidar_scan_t *scan,
                                 double max_window,
                                 double slope_threshold,
                                 double init_distance);

/* ─── L3: Point cloud statistics ───────────────────────────────────────── */

/** Compute centroid of point cloud */
lidar_vec3_t lidar_scan_centroid(const lidar_scan_t *scan);

/** Compute point density [points/m²] in a horizontal slice */
double lidar_scan_density(const lidar_scan_t *scan, double z_min, double z_max);

/** Compute elevation statistics: min, max, mean, std */
void lidar_scan_elevation_stats(const lidar_scan_t *scan,
                                 double *z_min, double *z_max,
                                 double *z_mean, double *z_std);

#endif /* LIDAR_GEOMETRY_H */