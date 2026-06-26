/**
 * track_core.h — Core Track Management Data Structures and API
 *
 * Covers: L1 (Definitions) — Track state machine, scoring, gating
 *         L2 (Core Concepts) — Track lifecycle, M/N logic, maintenance
 *         L4 (Fundamental Laws) — Reid's track score, Bar-Shalom's MHT
 *
 * References:
 *   - Bar-Shalom, Y., Willett, P.K., Tian, X. "Tracking and Data Fusion" (2011)
 *   - Blackman, S., Popoli, R. "Design and Analysis of Modern Tracking Systems" (1999)
 *   - Reid, D.B. "An Algorithm for Tracking Multiple Targets" IEEE TAC (1979)
 *   - Richards, Scheer, Holm "Principles of Modern Radar" (2010)
 *
 * Curriculum:
 *   - MIT 6.450 (Digital Communications), Stanford EE359 (Wireless), Michigan EECS 411
 *   - Georgia Tech ECE 6350 (EM), TU Munich High-Frequency Eng, ETH 227-0455
 */

#ifndef TRACK_CORE_H
#define TRACK_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L1 — Core Definitions: Track states, measurement types, gate shapes
 * ============================================================================
 */

/** Track lifecycle states (IEC 62388 radar tracking standard alignment) */
typedef enum {
    TRACK_STATE_FREE       = 0,  /**< Unused track slot */
    TRACK_STATE_TENTATIVE  = 1,  /**< Initialized, awaiting confirmation */
    TRACK_STATE_CONFIRMED  = 2,  /**< Confirmed, output to user */
    TRACK_STATE_COAST      = 3,  /**< Missing measurements, coasting */
    TRACK_STATE_DELETED    = 4   /**< Marked for deletion */
} track_state_t;

/** Measurement types supported by tracker */
typedef enum {
    MEAS_NONE         = 0,  /**< No measurement */
    MEAS_CARTESIAN_2D = 1,  /**< [x, y] in meters */
    MEAS_CARTESIAN_3D = 2,  /**< [x, y, z] in meters */
    MEAS_POLAR_2D     = 3,  /**< [range(m), bearing(rad)] */
    MEAS_POLAR_3D     = 4,  /**< [range, azimuth, elevation] */
    MEAS_DOPPLER_ONLY = 5,  /**< [range-rate (m/s)] */
    MEAS_RADAR_CUBE   = 6   /**< Range-Doppler-Angle cube point */
} measurement_type_t;

/** Gate shape for measurement-to-track association gating */
typedef enum {
    GATE_RECTANGULAR   = 0,  /**< Independent rectangular gates per dimension */
    GATE_ELLIPSOIDAL   = 1,  /**< Mahalanobis-distance-based ellipsoid */
    GATE_ADAPTIVE      = 2   /**< Adaptive gate based on innovation covariance */
} gate_shape_t;

/** Gating parameters structure */
typedef struct {
    gate_shape_t shape;          /**< Gate shape selection */
    double       gate_size;      /**< Gate probability threshold (e.g., 0.99) */
    double       gate_threshold; /**< Chi-squared threshold for gate_size, dof */
    int          gate_dof;       /**< Degrees of freedom for chi-squared */
    double       max_volume;     /**< Maximum gate volume constraint */
} gating_params_t;

/** M/N track confirmation logic parameters */
typedef struct {
    int M;  /**< Required number of detections (e.g., 3) */
    int N;  /**< Within this many scans (e.g., 5) */
    int   coast_max;  /**< Max coast cycles before deletion */
    double delete_score_threshold; /**< LLR score below which track is deleted */
    double confirm_score_threshold; /**< LLR score above which tentative → confirmed */
} track_mgmt_params_t;

/* ============================================================================
 * L1 — Track state vector and covariance
 * ============================================================================
 */

/** Maximum state dimension supported */
#define TRACK_MAX_STATE_DIM  9
/** Maximum measurement dimension supported */
#define TRACK_MAX_MEAS_DIM   6
/** Default track ID for uninitialized tracks */
#define TRACK_INVALID_ID     0

/**
 * Track structure — central entity representing one target hypothesis.
 *
 * State vector x is ordered: [pos_x, vel_x, acc_x, pos_y, vel_y, acc_y, pos_z, vel_z, acc_z]
 * for 3D CA model. Sub-models use leading subset.
 */
typedef struct {
    /* Identity */
    uint32_t      track_id;         /**< Unique track identifier */
    track_state_t state;             /**< Current lifecycle state */
    double        birth_time;        /**< Track creation time (seconds) */

    /* State estimate */
    int     state_dim;               /**< Dimensionality of state vector (e.g., 6 for 2D CA) */
    double  x[TRACK_MAX_STATE_DIM]; /**< State vector estimate */
    double  P[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM]; /**< Covariance matrix (row-major) */

    /* Measurement */
    int     meas_dim;                /**< Dimensionality of measurement */
    measurement_type_t meas_type;   /**< Type of measurement expected */

    /* Motion model identity */
    int motion_model_id; /**< Index into motion model table */

    /* Track quality */
    double score;           /**< Log-likelihood ratio (LLR) track score (Reid 1979) */
    int    updates_total;   /**< Total number of measurement updates */
    int    misses_since_last; /**< Consecutive scans without measurement */
    int    confirm_counter; /**< M/N logic: detections in recent N scans */
    int    history_window[16]; /**< Recent N-scan detection bitmask for M/N */

    /* Gating */
    gating_params_t gate;    /**< Gating parameters for this track */

    /* Prediction */
    double  x_pred[TRACK_MAX_STATE_DIM];  /**< Predicted state */
    double  P_pred[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM]; /**< Predicted covariance */
    double  innovation[TRACK_MAX_MEAS_DIM]; /**< Latest innovation vector */
    double  S[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM]; /**< Innovation covariance */
    double  mahalanobis_dist; /**< Latest Mahalanobis distance */

    /* Lifecycle */
    double time_since_birth; /**< Seconds since track creation */
    double last_update_time; /**< Timestamp of last valid update */

    /* Metadata */
    int    label;             /**< Object class label (-1 = unknown) */
    double existence_prob;    /**< Track existence probability (IPDA-style) */
    void  *user_data;         /**< User-defined pointer hook */
} track_t;

/* ============================================================================
 * L1 — Track cluster and multi-track container
 * ============================================================================
 */

/** Maximum tracks in a single tracker instance */
#define TRACKER_MAX_TRACKS 512
/** Maximum measurements per scan */
#define TRACKER_MAX_MEAS   2048

/** Measurement point (radar detection) */
typedef struct {
    uint32_t         meas_id;       /**< Measurement identifier */
    measurement_type_t type;        /**< Measurement type */
    double           time;          /**< Timestamp of measurement */
    int              dim;           /**< Actual measurement dimension */
    double           z[TRACK_MAX_MEAS_DIM]; /**< Measurement vector */
    double           R[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM]; /**< Measurement noise covariance */
    double           snr;           /**< SNR in dB */
    double           rcs;           /**< Radar cross section in dBsm */
    double           range_rate;    /**< Doppler-derived range rate if available */
    int              is_clutter;    /**< Flag: this detection may be false alarm */
    double           amplitude;     /**< Detection amplitude */
    void            *raw_data;      /**< Pointer to raw radar data */
} measurement_t;

/** Association pairing (measurement → track) */
typedef struct {
    uint32_t meas_id;   /**< Measurement ID */
    uint32_t track_id;  /**< Assigned track ID, 0 = unassigned/clutter */
    double   distance;  /**< Association distance metric (Mahalanobis) */
    double   likelihood; /**< Association likelihood */
    int      in_gate;   /**< Boolean: within validation gate */
} association_t;

/** Multi-object tracker container */
typedef struct {
    track_t      tracks[TRACKER_MAX_TRACKS]; /**< Track array */
    int          num_tracks;                 /**< Current active track count */
    uint32_t     next_track_id;              /**< Monotonic ID allocator */
    double       current_time;              /**< Current scan time */
    int          scan_count;                /**< Total scans processed */
    track_mgmt_params_t mgmt;               /**< Track management parameters */
    gating_params_t     default_gate;       /**< Default gating configuration */
    int          total_tracks_created;      /**< Cumulative tracks created */
    int          total_tracks_deleted;      /**< Cumulative tracks deleted */
    int          total_tracks_confirmed;    /**< Cumulative confirmed tracks */
} tracker_t;

/* ============================================================================
 * L2 — Core Concepts: Track lifecycle management functions
 * ============================================================================
 */

/**
 * Initialize the tracker container with default parameters.
 *
 * Complexity: O(1)
 */
void tracker_init(tracker_t *tracker);

/**
 * Initialize a single track slot to FREE state.
 *
 * Complexity: O(1)
 */
void track_init(track_t *track, uint32_t id, int state_dim, int meas_dim);

/**
 * Create a new track in the first free slot. Returns pointer to new track,
 * or NULL if tracker is full.
 *
 * Complexity: O(N_tracks)
 */
track_t *tracker_create_track(tracker_t *tracker, int state_dim, int meas_dim,
                               measurement_type_t meas_type);

/**
 * Delete a track by setting state to DELETED. Actual removal may be deferred
 * for output purposes. Marked tracks are reclaimed in tracker_cleanup().
 *
 * Complexity: O(1)
 */
void track_delete(track_t *track);

/**
 * Clean up all DELETED tracks, compacting the track array. Returns number
 * of tracks removed.
 *
 * Complexity: O(N_tracks)
 */
int tracker_cleanup(tracker_t *tracker);

/**
 * Apply M/N confirmation logic to all tentative tracks. Uses sliding window
 * of detections within the last N scans.
 *
 * Reference: Farina & Studer, "Radar Data Processing", Vol. I (1985)
 * Complexity: O(N_tracks * N_window)
 */
void tracker_apply_mn_logic(tracker_t *tracker);

/**
 * Prune tracks whose score falls below the deletion threshold.
 * Returns number of tracks pruned.
 *
 * Reference: Reid (1979) — track score pruning
 * Complexity: O(N_tracks)
 */
int tracker_prune_low_score(tracker_t *tracker);

/**
 * Coast tracks that have no measurement association in the current scan.
 * If coast_max is exceeded, the track is deleted.
 *
 * Complexity: O(N_tracks)
 */
void tracker_coast_orphan_tracks(tracker_t *tracker);

/* ============================================================================
 * L2 — Gating functions
 * ============================================================================
 */

/**
 * Compute Mahalanobis distance between measurement z and track prediction.
 *
 * d² = νᵀ S⁻¹ ν, where ν = z - Hx_pred is the innovation.
 *
 * Reference: Mahalanobis, P.C. "On the generalized distance in statistics" (1936)
 * Complexity: O(d³) for Cholesky-based inverse (d = measurement dimension)
 *
 * @param z     Measurement vector [meas_dim]
 * @param R     Measurement noise covariance [meas_dim * meas_dim]
 * @param x_pred Predicted state [state_dim]
 * @param P_pred Predicted covariance [state_dim * state_dim]
 * @param H     Measurement matrix [meas_dim * state_dim]
 * @param state_dim State dimension
 * @param meas_dim  Measurement dimension
 * @return Mahalanobis distance squared, or -1.0 on failure
 */
double mahalanobis_distance_sq(const double *z, const double *R,
                                const double *x_pred, const double *P_pred,
                                const double *H, int state_dim, int meas_dim);

/**
 * Check if a measurement falls within the ellipsoidal validation gate.
 *
 * Gate condition: d² ≤ γ, where γ = chi2inv(P_G, n_z) from chi-squared
 * distribution with n_z degrees of freedom at gate probability P_G.
 *
 * @param mahal_sq  Squared Mahalanobis distance
 * @param threshold Chi-squared threshold for given P_G and n_z
 * @return 1 if in gate, 0 otherwise
 */
int ellipsoidal_gate_check(double mahal_sq, double threshold);

/**
 * Compute chi-squared threshold for specified gate probability.
 * Uses the Wilson-Hilferty approximation for computational efficiency.
 *
 * γ ≈ n_z · (1 − 2/(9·n_z) + z_P · √(2/(9·n_z)))³
 * where z_P is the standard normal quantile for probability P_G.
 *
 * Reference: Wilson, E.B., Hilferty, M.M. (1931)
 * Complexity: O(1)
 *
 * @param P_G Gate probability (e.g., 0.99, 0.999)
 * @param n_z Degrees of freedom (measurement dimension)
 * @return Chi-squared threshold γ
 */
double chi2_threshold_approx(double P_G, int n_z);

/**
 * Compute the standard normal quantile (probit function) via
 * the Abramowitz-Stegun approximation (26.2.23).
 *
 * Complexity: O(1)
 *
 * @param p Probability (0 < p < 1)
 * @return z-score
 */
double normal_quantile(double p);

/* ============================================================================
 * L2 — Track scoring (Log-Likelihood Ratio)
 * ============================================================================
 */

/**
 * Compute track score increment from a measurement association.
 *
 * ΔLLR = log(P_D / (β_F · (2π)^{n_z/2} · √|S|)) − d²/2
 *
 * where P_D = probability of detection, β_F = false alarm density,
 * S = innovation covariance, d² = Mahalanobis distance.
 *
 * Reference: Reid, D.B. "An Algorithm for Tracking Multiple Targets" (1979)
 * Complexity: O(n_z³) for determinant
 *
 * @param P_D    Probability of detection
 * @param beta_F False alarm spatial density
 * @param S      Innovation covariance matrix [n_z × n_z]
 * @param d2     Squared Mahalanobis distance
 * @param n_z    Measurement dimension
 * @return Track score increment
 */
double track_score_increment(double P_D, double beta_F, const double *S,
                              double d2, int n_z);

/**
 * Compute track score penalty for missed detection.
 *
 * ΔLLR_miss = log(1 − P_D)
 *
 * Complexity: O(1)
 */
double track_score_miss_penalty(double P_D);

/**
 * Update track score — adds increment or penalty based on association status.
 *
 * Complexity: O(n_z³)
 */
void track_update_score(track_t *track, int associated, double P_D,
                         double beta_F, double d2);

/* ============================================================================
 * L2 — Utility: matrix operations for small fixed-size matrices
 * ============================================================================
 */

/**
 * Matrix-vector multiply: y = A * x.
 * Complexity: O(m * n)
 */
void mat_vec_mul(const double *A, const double *x, double *y, int m, int n);

/**
 * Matrix-matrix multiply: C = A * B (row-major storage).
 * Complexity: O(m * k * n)
 */
void mat_mat_mul(const double *A, const double *B, double *C,
                  int m, int k, int n);

/**
 * Matrix transpose: B = Aᵀ.
 * Complexity: O(m * n)
 */
void mat_transpose(const double *A, double *B, int m, int n);

/**
 * Matrix inverse via Cholesky decomposition for symmetric positive definite.
 * A must be SPD. Returns determinant or 0.0 on failure.
 * Complexity: O(n³)
 */
double mat_inv_cholesky(const double *A, double *A_inv, int n);

/**
 * Matrix determinant of symmetric positive definite via Cholesky.
 * Complexity: O(n³)
 */
double mat_det_cholesky(const double *A, int n);

/**
 * Solve A x = b for SPD A via Cholesky.
 * Complexity: O(n³)
 * Returns 0 on success, -1 on failure.
 */
int mat_solve_cholesky(const double *A, const double *b, double *x, int n);

/**
 * Vector dot product.
 * Complexity: O(n)
 */
double vec_dot(const double *a, const double *b, int n);

/**
 * Vector 2-norm.
 * Complexity: O(n)
 */
double vec_norm2(const double *v, int n);

/**
 * Vector subtraction: c = a - b.
 * Complexity: O(n)
 */
void vec_sub(const double *a, const double *b, double *c, int n);

/**
 * Vector addition: c = a + b.
 * Complexity: O(n)
 */
void vec_add(const double *a, const double *b, double *c, int n);

/**
 * Scalar-vector multiply: y = alpha * x.
 * Complexity: O(n)
 */
void vec_scale(const double *x, double alpha, double *y, int n);

/**
 * Allocate and zero a double array.
 * Complexity: O(n)
 */
void vec_zero(double *v, int n);

/**
 * Copy vector: dst = src.
 * Complexity: O(n)
 */
void vec_copy(const double *src, double *dst, int n);

/**
 * Convert bearing from degrees to radians.
 */
double deg_to_rad(double deg);

/**
 * Convert bearing from radians to degrees.
 */
double rad_to_deg(double rad);

/**
 * Wrap angle to [-pi, pi].
 */
double wrap_angle(double angle);

/**
 * Wrap angle to [0, 2*pi].
 */
double wrap_angle_2pi(double angle);

#ifdef __cplusplus
}
#endif

#endif /* TRACK_CORE_H */
