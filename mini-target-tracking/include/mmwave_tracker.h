/**
 * mmwave_tracker.h — Automotive mmWave Radar Target Tracking Applications
 *
 * Covers: L7 (Applications) — TI AWR1642/1843 processing chain, DBSCAN clustering
 *         L7 — Automotive radar object classification, ego-motion compensation
 *         L8 (Advanced Topics) — Occupancy grid mapping, extended object tracking
 *         L8 — Random Finite Set introduction, group tracking concepts
 *
 * References:
 *   - Texas Instruments "AWR1642/AWR1843 Technical Reference Manual" (2020)
 *   - Ester, M., Kriegel, H.P., et al. "A density-based algorithm for
 *     discovering clusters" (DBSCAN, 1996)
 *   - Granström, K., Baum, M., Reuter, S. "Extended Object Tracking:
 *     Introduction, Overview, Applications" JAIF (2017)
 *   - Mahler, R.P.S. "Statistical Multisource-Multitarget Information Fusion" (2007)
 *   - Werber, K. et al. "Automotive Radar Gridmap Representations" (2015)
 *
 * Curriculum:
 *   - MIT 6.450, Stanford EE359, Michigan EECS 411 (Automotive)
 *   - Georgia Tech ECE 6601, TU Munich High-Frequency Eng
 */

#ifndef MMWAVE_TRACKER_H
#define MMWAVE_TRACKER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"

/* ============================================================================
 * L7 — TI mmWave radar data structures
 * ============================================================================
 */

/** TI mmWave SDK detection object (emulating DPIF_PointCloudCartesian) */
typedef struct {
    float x;            /**< X coordinate (m) */
    float y;            /**< Y coordinate (m) */
    float z;            /**< Z coordinate (m), 0 for 2D */
    float velocity;     /**< Radial velocity (m/s) */
    float snr;          /**< SNR (dB) */
    uint16_t noise;     /**< Noise level */
} mmwave_point_t;

/** Cluster label assignments */
typedef struct {
    int     num_points;  /**< Number of input points */
    int    *labels;      /**< Cluster label for each point (0 = noise) */
    int     num_clusters; /**< Number of clusters found (excluding noise) */
} dbscan_result_t;

/** Object classification type (automotive) */
typedef enum {
    OBJ_UNKNOWN   = 0,  /**< Unclassified */
    OBJ_PEDESTRIAN = 1, /**< Pedestrian / VRU */
    OBJ_BICYCLE   = 2,  /**< Bicycle / motorcycle */
    OBJ_CAR       = 3,  /**< Passenger vehicle */
    OBJ_TRUCK     = 4,  /**< Truck / bus */
    OBJ_STATIC    = 5   /**< Static infrastructure */
} object_class_t;

/** Automotive object descriptor */
typedef struct {
    int             id;              /**< Object ID */
    object_class_t  class_type;      /**< Classification */
    double          x;               /**< Position x (m) */
    double          y;               /**< Position y (m) */
    double          vx;              /**< Velocity x (m/s) */
    double          vy;              /**< Velocity y (m/s) */
    double          length;          /**< Estimated object length (m) */
    double          width;           /**< Estimated object width (m) */
    double          orientation;     /**< Estimated orientation (rad) */
    double          rcs_mean;        /**< Mean RCS (dBsm) */
    double          existence_prob;  /**< Track existence probability */
    int             num_points;      /**< Number of radar points in this object */
    double          age;             /**< Object track age (s) */
} automotive_object_t;

/* ============================================================================
 * L7 — DBSCAN point cloud clustering
 * ============================================================================
 */

/**
 * DBSCAN clustering of radar point cloud.
 *
 * Density-Based Spatial Clustering of Applications with Noise.
 *
 * Definition: A point p is a core point if at least minPts points
 * (including p) are within distance eps of p.
 *
 * Algorithm:
 *   1. For each unvisited point p:
 *      a. Find all points within eps (neighborhood N(p))
 *      b. If |N(p)| ≥ minPts: p is a core point, start new cluster
 *         i.  For each point q in cluster:
 *             Find N(q). If |N(q)| ≥ minPts, add N(q) to cluster.
 *      c. Else: label p as noise (for now)
 *
 * Key parameters:
 *   - eps: neighborhood radius (m), based on object size and sensor spread
 *   - minPts: minimum points for core point, typically 3-5 for mmWave radar
 *
 * Reference: Ester et al. (1996)
 * Complexity: O(N·log(N)) with spatial index, O(N²) naive
 *
 * @param points  Array of radar points
 * @param n       Number of points
 * @param eps     Neighborhood radius (e.g., 0.5 m for pedestrian, 1.5 m for car)
 * @param minPts  Minimum points per cluster (e.g., 3)
 * @return Clustering result (free with dbscan_free_result)
 */
dbscan_result_t *dbscan_cluster(const mmwave_point_t *points, int n,
                                 double eps, int minPts);

/**
 * Free DBSCAN result.
 */
void dbscan_free_result(dbscan_result_t *result);

/**
 * Adaptive DBSCAN: eps is chosen per-point based on range-dependent
 * radar resolution. Points at longer range have larger eps.
 *
 * eps(r) = eps_min + (r / r_max) · (eps_max − eps_min)
 *
 * This reflects the increasing measurement uncertainty with range.
 *
 * Complexity: O(N²) naive, O(N·log(N)) with spatial index
 */
dbscan_result_t *dbscan_cluster_adaptive(const mmwave_point_t *points, int n,
                                           double eps_min, double eps_max,
                                           double r_max, int minPts);

/* ============================================================================
 * L7 — Object classification based on radar features
 * ============================================================================
 */

/**
 * Classify a clustered object based on radar features.
 *
 * Features used:
 *   - Number of points (spatial extent)
 *   - Velocity profile (dispersion, mean)
 *   - RCS statistics (mean, variance)
 *   - Shape (eigenvalues of point cloud)
 *   - Micro-Doppler signature (if available)
 *
 * Simple rule-based classifier:
 *   - Low point count + low RCS + low speed → likely pedestrian
 *   - Medium point count + medium RCS + high speed → likely car
 *   - High point count + high RCS + low speed → likely truck
 *   - Zero velocity + spread points → static object
 *
 * Reference: Heuel, S., Rohling, H. "Pedestrian recognition in automotive
 *            radar sensors" (2012)
 * Complexity: O(n_points) per object
 *
 * @param points     Points belonging to a single cluster
 * @param n_points   Number of points in cluster
 * @param labels     Cluster point indices
 * @param centroid_x Cluster centroid x (m)
 * @param centroid_y Cluster centroid y (m)
 * @return Classified object type
 */
object_class_t classify_radar_object(const mmwave_point_t *points, int n_points,
                                      const int *labels,
                                      double centroid_x, double centroid_y);

/**
 * Compute point cloud shape features (eigenvalues, orientation).
 *
 * From the spatial covariance of the point cloud:
 *   C = (1/N) Σ (p_i − μ)(p_i − μ)ᵀ
 *
 * Eigenvalues λ₁ ≥ λ₂ give semi-major/minor axis lengths:
 *   length ≈ 2·√(λ₁), width ≈ 2·√(λ₂)
 *
 * Orientation angle: θ = 0.5·atan2(2·C_xy, C_xx − C_yy)
 *
 * Complexity: O(N)
 *
 * @param points      Cluster points
 * @param n_points    Number of points
 * @param centroid_x  Cluster centroid x
 * @param centroid_y  Cluster centroid y
 * @param length      Output: estimated length (m)
 * @param width       Output: estimated width (m)
 * @param orientation Output: orientation angle (rad)
 */
void compute_object_shape(const mmwave_point_t *points, int n_points,
                           double centroid_x, double centroid_y,
                           double *length, double *width, double *orientation);

/* ============================================================================
 * L7 — Ego-motion compensation
 * ============================================================================
 */

/**
 * Compensate radar measurements for ego-vehicle motion.
 *
 * When the ego vehicle moves, stationary objects appear to have relative
 * velocity. Ego-motion compensation removes this artifact.
 *
 * Given ego velocity (v_ego_x, v_ego_y) and yaw rate ω_ego:
 *   v_compensated = v_measured − v_ego_projected
 *
 * where v_ego_projected = (v_ego − ω×r)·û_r
 * (projection of ego velocity onto the radial direction to the target)
 *
 * This separates truly moving targets from stationary clutter.
 *
 * Reference: Kellner, D. et al. "Instantaneous ego-motion estimation
 *            using Doppler radar" (2013)
 * Complexity: O(N) for N points
 *
 * @param points        Radar point array, velocities updated in place
 * @param n             Number of points
 * @param v_ego_x       Ego longitudinal velocity (m/s)
 * @param v_ego_y       Ego lateral velocity (m/s)
 * @param yaw_rate_ego  Ego yaw rate (rad/s)
 * @param sensor_x      Sensor x-offset from vehicle center (m)
 * @param sensor_y      Sensor y-offset from vehicle center (m)
 */
void ego_motion_compensate(mmwave_point_t *points, int n,
                            double v_ego_x, double v_ego_y,
                            double yaw_rate_ego,
                            double sensor_x, double sensor_y);

/**
 * Detect moving objects after ego-motion compensation.
 *
 * Points with |v_compensated| > threshold are classified as moving.
 *
 * Complexity: O(N)
 *
 * @param points    Radar points (already ego-motion compensated)
 * @param n         Number of points
 * @param threshold Velocity threshold (m/s), e.g., 0.5
 * @return Number of moving points detected
 */
int detect_moving_points(const mmwave_point_t *points, int n, double threshold);

/* ============================================================================
 * L8 — Occupancy grid mapping (OGM)
 * ============================================================================
 */

/** Occupancy grid cell */
typedef struct {
    double prob_occupied;   /**< Probability of occupancy p(m_i | z_{1:t}) */
    double log_odds;        /**< Log-odds representation l_i = log(p/(1−p)) */
    int    num_hits;        /**< Number of times this cell was observed occupied */
    int    num_misses;      /**< Number of times this cell was observed free */
    double mean_velocity;   /**< Mean Doppler velocity in this cell */
    int    last_update_scan; /**< Last scan index when this cell was updated */
} occupancy_cell_t;

/** Occupancy grid map */
typedef struct {
    int             width;   /**< Grid width (cells) */
    int             height;  /**< Grid height (cells) */
    double          resolution; /**< Cell size (m) */
    double          origin_x;/**< Origin x in global coordinates */
    double          origin_y;/**< Origin y in global coordinates */
    occupancy_cell_t *cells; /**< Cell array [width × height] (row-major) */
    double          prob_prior; /**< Prior occupancy probability */
    double          prob_occupied_given_hit;  /**< P(m=1|z=hit) */
    double          prob_occupied_given_miss; /**< P(m=1|z=miss) */
} occupancy_grid_t;

/**
 * Initialize an occupancy grid.
 *
 * Complexity: O(W·H)
 *
 * @param width      Grid width in cells
 * @param height     Grid height in cells
 * @param resolution Cell size in meters
 * @param origin_x   Global x of grid origin (lower-left corner)
 * @param origin_y   Global y of grid origin
 */
occupancy_grid_t *ogm_create(int width, int height, double resolution,
                              double origin_x, double origin_y);

/**
 * Free occupancy grid.
 */
void ogm_free(occupancy_grid_t *grid);

/**
 * Update occupancy grid with a radar scan using the binary Bayes filter.
 *
 * Log-odds update:
 *   l_{i,k} = l_{i,k-1} + log( P(z_k|m_i=1) / P(z_k|m_i=0) ) − l₀
 *
 * where l₀ = log(p₀/(1−p₀)) is the log-odds prior.
 *
 * Ray-casting: cells along the line from sensor to detection are marked
 * as free; cells at detection range are marked as occupied.
 *
 * Reference: Thrun, S., Burgard, W., Fox, D. "Probabilistic Robotics" (2005)
 * Complexity: O(N·R) where N = points, R = range/resolution
 *
 * @param grid            Occupancy grid to update
 * @param points          Radar points in sensor frame
 * @param n_points        Number of points
 * @param sensor_x        Sensor x in global frame
 * @param sensor_y        Sensor y in global frame
 * @param sensor_yaw      Sensor yaw (rad)
 */
void ogm_update_radar_scan(occupancy_grid_t *grid,
                            const mmwave_point_t *points, int n_points,
                            double sensor_x, double sensor_y, double sensor_yaw);

/**
 * Query occupancy probability at a location.
 * Complexity: O(1)
 *
 * @return Occupancy probability [0, 1], or -1 if out of bounds
 */
double ogm_query(const occupancy_grid_t *grid, double x, double y);

/**
 * Extract the set of occupied cells (above threshold) as point-like obstacles.
 *
 * Complexity: O(W·H)
 *
 * @param threshold Probability threshold (e.g., 0.7)
 * @param obstacles Output: array of (x, y) occupied positions
 * @param max_obs   Max obstacles to extract
 * @return Number of obstacles extracted
 */
int ogm_extract_obstacles(const occupancy_grid_t *grid, double threshold,
                           double *obstacles_x, double *obstacles_y,
                           int max_obs);

/* ============================================================================
 * L8 — Extended object tracking concepts
 * ============================================================================
 */

/**
 * Compute the random matrix estimate for an extended object.
 *
 * The extended object is modeled as an ellipse with kinematic state
 * (position, velocity) and extent state (shape matrix X).
 *
 * The random matrix approach models the extent as an inverse-Wishart
 * distributed random matrix:
 *
 *   p(X_k | z_{1:k}) ≈ IW(X_k; v_k, V_k)
 *
 * where v_k is the degrees of freedom and V_k is the scale matrix.
 *
 * Measurement likelihood (spatial distribution model):
 *   z_k ~ N(H·x_k, s·X_k + R)
 *
 * where H is the measurement matrix, x_k is the kinematic state,
 * s is a scaling factor, X_k is the extent matrix, R is sensor noise.
 *
 * Reference: Koch, J.W. "Bayesian approach to extended object and cluster
 *            tracking using random matrices" IEEE T-AES (2008)
 *            Granström et al. (2017)
 * Complexity: O(d³) for d-dimensional extent
 *
 * @param measurements Measurements assigned to this object [n_meas × meas_dim]
 * @param n_meas       Number of measurements
 * @param meas_dim     Measurement dimension (2 for 2D)
 * @param x_kin        Kinematic state [state_dim], updated in place
 * @param P_kin        Kinematic covariance [state_dim²], updated in place
 * @param V_extent     Extent scale matrix [meas_dim²], updated in place
 * @param nu           Degrees of freedom for inverse-Wishart, updated in place
 * @param H            Measurement matrix [meas_dim × state_dim]
 * @param state_dim    State dimension
 */
void extended_object_random_matrix(const double *measurements, int n_meas,
                                    int meas_dim, double *x_kin, double *P_kin,
                                    double *V_extent, double *nu,
                                    const double *H, int state_dim);

/**
 * Compute ellipse parameters (semi-axes, orientation) from extent matrix.
 *
 * For extent matrix X (2×2 positive definite):
 *   λ₁, λ₂ = eigenvalues of X
 *   semi_major = c·√(λ₁), semi_minor = c·√(λ₂)
 *   orientation = 0.5·atan2(2·X₁₂, X₁₁−X₂₂)
 *
 * where c is a scaling constant (e.g., c = 2 for 95% containment).
 *
 * Complexity: O(1) for 2×2 matrix
 */
void extent_to_ellipse(const double *X, double c,
                        double *semi_major, double *semi_minor,
                        double *orientation);

/* ============================================================================
 * L8 — Group tracking
 * ============================================================================
 */

/**
 * Detect if a set of targets form a group (move together).
 *
 * Group definition: targets are considered a group if:
 *   1. Pairwise distance < d_max
 *   2. Velocity difference < v_max
 *   3. Targets stay together for ≥ T_group scans
 *
 * Uses a connectivity graph approach: build adjacency matrix, find
 * connected components, each component is a candidate group.
 *
 * Complexity: O(N²) for N targets
 *
 * @param positions  Target positions [N × 2] (x, y in columns)
 * @param velocities Target velocities [N × 2] (vx, vy in columns)
 * @param N          Number of targets
 * @param d_max      Maximum inter-target distance for grouping (m)
 * @param v_max      Maximum velocity difference for grouping (m/s)
 * @param groups     Output: group_id for each target [N], -1 = solo
 * @param n_groups   Output: number of distinct groups found
 */
void detect_groups(const double *positions, const double *velocities,
                    int N, double d_max, double v_max,
                    int *groups, int *n_groups);

/**
 * Compute group centroid and bounding box.
 *
 * Centroid is the mean of member positions.
 * Bounding box is the axis-aligned rectangle containing all members.
 *
 * Complexity: O(N_group)
 */
void compute_group_bounds(const double *positions, int N,
                           const int *group_ids, int target_group,
                           double *centroid_x, double *centroid_y,
                           double *bbox_xmin, double *bbox_xmax,
                           double *bbox_ymin, double *bbox_ymax);

#ifdef __cplusplus
}
#endif

#endif /* MMWAVE_TRACKER_H */
