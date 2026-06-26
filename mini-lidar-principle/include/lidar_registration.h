/**
 * @file    lidar_registration.h
 * @brief   Point cloud registration — ICP and variants
 *
 * Knowledge covered:
 *   L1: Point cloud registration, rigid transform, correspondence
 *   L3: Singular Value Decomposition (SVD), orthogonal Procrustes problem
 *   L5: Iterative Closest Point (ICP), point-to-point and point-to-plane
 *       variants, transformation estimation from correspondences
 *   L6: Multi-view registration, scan alignment
 *
 * Reference:
 *   - Besl & McKay, "A Method for Registration of 3-D Shapes",
 *     *IEEE TPAMI* 14(2), pp.239-256, 1992 (original ICP).
 *   - Chen & Medioni, "Object Modelling by Registration of Multiple
 *     Range Images", *Image and Vision Computing* 10(3), pp.145-155, 1992.
 *   - Rusinkiewicz & Levoy, "Efficient Variants of the ICP Algorithm",
 *     *3DIM*, 2001.
 *   - Arun, Huang & Blostein, "Least-Squares Fitting of Two 3-D Point Sets",
 *     *IEEE TPAMI* 9(5), pp.698-700, 1987 (SVD method).
 */

#ifndef LIDAR_REGISTRATION_H
#define LIDAR_REGISTRATION_H

#include "lidar_core.h"
#include "lidar_geometry.h"
#include <stddef.h>

/* ─── L1: ICP configuration and result ──────────────────────────────────── */

/**
 * @brief ICP algorithm configuration
 */
typedef struct {
    int    max_iterations;      /**< Maximum ICP iterations (typical: 50-100) */
    double translation_tol;    /**< Convergence tolerance: translation change [m] */
    double rotation_tol;       /**< Convergence tolerance: rotation change [rad] */
    double max_correspondence_dist; /**< Maximum distance for point correspondence [m] */
    double outlier_rejection_std;   /**< Reject correspondences > N * std dev (0 = no rejection) */
    int    use_point_to_plane; /**< 0 = point-to-point, 1 = point-to-plane */
    int    verbose;            /**< Print iteration progress */
} lidar_icp_config_t;

/**
 * @brief ICP registration result
 */
typedef struct {
    lidar_mat4_t transform;       /**< Final 4x4 rigid transform (source → target) */
    double       rmse;            /**< Root-mean-square error of correspondences [m] */
    int          num_iterations;  /**< Number of iterations performed */
    int          converged;       /**< 1 if convergence criterion met */
    size_t       num_corresp;     /**< Number of correspondences in final iteration */
    double       *rmse_history;   /**< RMSE per iteration (allocated internally) */
} lidar_icp_result_t;

/* ─── L5: SVD-based rigid transform estimation ──────────────────────────── */

/**
 * @brief Estimate rigid transform from two sets of corresponding points
 *
 * Given source points {p_i} and target points {q_i} with known correspondence,
 * find the rotation R and translation t that minimize:
 *
 *   min Σ || R·p_i + t - q_i ||²
 *
 * Solution via the orthogonal Procrustes problem:
 *   1. Compute centroids: p̄ = (1/N)·Σ p_i,  q̄ = (1/N)·Σ q_i
 *   2. Center points: p'_i = p_i - p̄,  q'_i = q_i - q̄
 *   3. Cross-covariance: H = Σ p'_i · q'_i^T
 *   4. SVD: H = U · Σ · V^T
 *   5. Rotation: R = V · U^T  (corrected for reflection: det(R) must be +1)
 *   6. Translation: t = q̄ - R · p̄
 *
 * Complexity: O(N) (SVD of 3x3 is O(27) constant)
 *
 * Reference: Arun, Huang & Blostein (1987).
 *
 * @param src       Array of source points
 * @param dst       Array of target (destination) points
 * @param num_pts   Number of corresponding point pairs
 * @param R         Output: 3x3 rotation matrix
 * @param t         Output: translation vector
 * @return          0 on success, -1 on degenerate input
 */
int lidar_estimate_rigid_transform(const lidar_vec3_t *src,
                                     const lidar_vec3_t *dst,
                                     size_t num_pts,
                                     lidar_mat3_t *R,
                                     lidar_vec3_t *t);

/**
 * @brief Apply rigid transform to a set of points (in-place)
 *
 * p'_i = R · p_i + t
 *
 * @param points   Array of points (modified in-place)
 * @param num_pts  Number of points
 * @param R        Rotation matrix
 * @param t        Translation vector
 */
void lidar_apply_transform(lidar_vec3_t *points, size_t num_pts,
                            lidar_mat3_t R, lidar_vec3_t t);

/* ─── L5: Iterative Closest Point (ICP) ─────────────────────────────────── */

/**
 * @brief Point-to-point ICP registration
 *
 * Classic ICP algorithm:
 *   while (not converged):
 *     1. For each source point, find nearest neighbor in target
 *     2. Reject outlier correspondences
 *     3. Estimate optimal R, t from correspondences
 *     4. Apply transform to source
 *     5. Check convergence
 *
 * Complexity per iteration: O(N_src · N_tgt) for brute-force NN,
 * or O(N_src · log N_tgt) with KD-tree.
 *
 * Reference: Besl & McKay (1992).
 *
 * @param source   Source point cloud (moved to align with target)
 * @param target   Target (reference) point cloud
 * @param config   ICP configuration parameters
 * @param result   Output registration result
 * @return         0 on success, -1 on error
 */
int lidar_icp_point_to_point(lidar_scan_t *source,
                               const lidar_scan_t *target,
                               const lidar_icp_config_t *config,
                               lidar_icp_result_t *result);

/**
 * @brief ICP with specified initial transform
 *
 * Same as lidar_icp_point_to_point but starts from a given initial
 * transform (coarse alignment from GPS/INS or feature matching).
 *
 * @param source        Source point cloud
 * @param target        Target point cloud
 * @param init_transform Initial rigid transform (applied before ICP loop)
 * @param config        ICP configuration
 * @param result        Output result
 * @return              0 on success
 */
int lidar_icp_with_initial(lidar_scan_t *source,
                             const lidar_scan_t *target,
                             lidar_mat4_t init_transform,
                             const lidar_icp_config_t *config,
                             lidar_icp_result_t *result);

/**
 * @brief Initialize ICP configuration with sensible defaults
 *
 * Defaults:
 *   max_iterations = 100
 *   translation_tol = 1e-6 m
 *   rotation_tol = 1e-6 rad
 *   max_correspondence_dist = 1.0 m
 *   outlier_rejection_std = 2.5
 *   use_point_to_plane = 0
 *   verbose = 0
 *
 * @param config  ICP config to initialize
 */
void lidar_icp_config_default(lidar_icp_config_t *config);

/**
 * @brief Free memory associated with ICP result
 */
void lidar_icp_result_free(lidar_icp_result_t *result);

/* ─── L5: Nearest-neighbor correspondence ───────────────────────────────── */

/**
 * @brief Compute closest-point correspondences (brute-force)
 *
 * For each point in source, find the nearest neighbor in target.
 * O(N_src · N_tgt).
 *
 * @param src           Source points
 * @param num_src       Number of source points
 * @param tgt           Target points
 * @param num_tgt       Number of target points
 * @param max_dist      Maximum correspondence distance [m]
 * @param correspondences Output: indices of corresponding target points
 *                        (length = num_src, -1 if no correspondence within max_dist)
 * @param distances     Output: distances to corresponding points
 * @return              Number of valid correspondences found
 */
size_t lidar_brute_force_correspondences(const lidar_vec3_t *src, size_t num_src,
                                           const lidar_vec3_t *tgt, size_t num_tgt,
                                           double max_dist,
                                           int *correspondences,
                                           double *distances);

/* ─── L5: ICP fitness evaluation ────────────────────────────────────────── */

/**
 * @brief Compute RMSE between two aligned point clouds
 *
 * RMSE = sqrt( (1/N) · Σ ||p_i - q_closest(i)||² )
 *
 * @param source   Aligned source point cloud
 * @param target   Target point cloud
 * @param max_dist Maximum distance for valid correspondence
 * @return         RMSE [m], -1 if no correspondences
 */
double lidar_registration_rmse(const lidar_scan_t *source,
                                 const lidar_scan_t *target,
                                 double max_dist);

/**
 * @brief Compute overlap ratio: fraction of source points with
 *        a corresponding target point within max_dist
 *
 * @param source   Source point cloud
 * @param target   Target point cloud
 * @param max_dist Maximum correspondence distance
 * @return         Overlap ratio [0-1]
 */
double lidar_registration_overlap(const lidar_scan_t *source,
                                    const lidar_scan_t *target,
                                    double max_dist);

#endif /* LIDAR_REGISTRATION_H */