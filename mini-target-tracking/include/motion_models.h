/**
 * motion_models.h — Target Motion Models for State Prediction
 *
 * Covers: L1 (Definitions) — Motion model types, process noise models
 *         L2 (Core Concepts) — White-noise acceleration, Markov acceleration
 *         L3 (Mathematical Structures) — State transition matrices, discretization
 *         L4 (Fundamental Laws) — Wiener process, Singer model, coordinated turn
 *         L5 (Algorithms/Methods) — IMM model set design, discretization
 *
 * References:
 *   - Li, X.R., Jilkov, V.P. "Survey of Maneuvering Target Tracking" (2003)
 *   - Singer, R.A. "Estimating Optimal Tracking Filter Performance" (1970)
 *   - Bar-Shalom, Willett, Tian "Tracking and Data Fusion" (2011), Ch. 6
 *   - Blackman, Popoli "Design and Analysis of Modern Tracking Systems" (1999)
 *
 * Curriculum:
 *   - MIT 6.003, Stanford EE102A, ETH 227-0427, 清华 (信号与系统)
 *   - Michigan EECS 351, Georgia Tech ECE 4270
 */

#ifndef MOTION_MODELS_H
#define MOTION_MODELS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"

/* ============================================================================
 * L1 — Motion model type enumeration
 * ============================================================================
 */

/** Supported motion models */
typedef enum {
    MOTION_CV   = 0,  /**< Constant Velocity (second-order) */
    MOTION_CA   = 1,  /**< Constant Acceleration (third-order) */
    MOTION_CT   = 2,  /**< Coordinated Turn (constant turn rate) */
    MOTION_SINGER = 3, /**< Singer model (exponentially correlated acceleration) */
    MOTION_WPA  = 4,  /**< Wiener Process Acceleration (white jerk) */
    MOTION_STOP = 5,  /**< Stopped target model (zero velocity) */
    MOTION_BROWN = 6, /**< Brownian motion (random walk, first-order) */
    MOTION_NCV  = 7,  /**< Nearly Constant Velocity (small white acceleration) */
    MOTION_NCT  = 8,  /**< Nearly Coordinated Turn */
    MOTION_CURVILINEAR = 9 /**< General 2D curvilinear motion */
} motion_model_type_t;

/** Motion model dimensionality */
typedef enum {
    DIM_1D = 1,  /**< 1D: [x, vx] or [x, vx, ax] */
    DIM_2D = 2,  /**< 2D: [x, vx, ax, y, vy, ay] */
    DIM_3D = 3   /**< 3D: [x, vx, ax, y, vy, ay, z, vz, az] */
} motion_dim_t;

/* ============================================================================
 * L1 — Motion model parameters structure
 * ============================================================================
 */

/** Motion model descriptor with pre-computed transition matrix and noise */
typedef struct {
    motion_model_type_t type;        /**< Model type */
    motion_dim_t        dim;         /**< Dimensionality */
    int                 state_dim;   /**< Total state vector dimension */
    int                 num_params;  /**< Number of continuous-time parameters */

    /* Discretized model matrices */
    double *F;  /**< State transition matrix [state_dim × state_dim] */
    double *Q;  /**< Process noise covariance [state_dim × state_dim] */

    /* Continuous-time parameters */
    double   T;        /**< Sampling interval (seconds) */
    double   q;        /**< Process noise intensity (power spectral density) */
    double   alpha;    /**< Singer model: maneuver frequency (1/tau) */
    double   omega;    /**< Coordinated turn: turn rate (rad/s) */
    double   sigma_m;  /**< Singer model: maneuver standard deviation */

    /* Pre-computed helper quantities */
    double  *Gamma;    /**< Noise gain matrix [state_dim × noise_dim] */
    int      noise_dim; /**< Dimension of continuous-time process noise */

    /* IMM transition probability from this model */
    double   imm_prior; /**< Prior probability of this model */
} motion_model_t;

/* ============================================================================
 * L2+L3 — Model construction functions
 * ============================================================================
 */

/**
 * Allocate and construct a Constant Velocity (CV) model.
 *
 * Continuous-time: ẍ(t) = w(t),  w(t) ~ N(0, q)
 *
 * State vector (2D): [x, ẋ, y, ẏ]
 * Discrete-time F: piecewise-constant white acceleration model
 *
 * F = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1]
 * Q = q · [T³/3 T²/2 0   0  ;
 *          T²/2 T    0   0  ;
 *          0    0    T³/3 T²/2;
 *          0    0    T²/2 T   ]
 *
 * Reference: Bar-Shalom (2011), Eq. 6.3.2-2
 * Complexity: O(state_dim²) allocation
 *
 * @param dim Motion dimensionality (DIM_1D, DIM_2D, DIM_3D)
 * @param T   Sampling interval (seconds)
 * @param q   Process noise PSD for acceleration
 * @return Allocated model, or NULL on failure
 */
motion_model_t *motion_create_cv(motion_dim_t dim, double T, double q);

/**
 * Allocate and construct a Constant Acceleration (CA) model.
 *
 * Continuous-time: x⃛(t) = w(t), w(t) ~ N(0, q)
 *
 * State vector (2D): [x, ẋ, ẍ, y, ẏ, ÿ]
 * Discrete-time F: third-order integration of white jerk
 *
 * F = [1 T T²/2 0 0 0   ;
 *      0 1 T    0 0 0   ;
 *      0 0 1    0 0 0   ;
 *      0 0 0    1 T T²/2;
 *      0 0 0    0 1 T   ;
 *      0 0 0    0 0 1   ]
 *
 * Q = q · [T⁵/20 T⁴/8 T³/6 0     0     0    ;
 *          T⁴/8  T³/3 T²/2 0     0     0    ;
 *          T³/6  T²/2 T    0     0     0    ;
 *          0     0    0    T⁵/20 T⁴/8  T³/6;
 *          0     0    0    T⁴/8  T³/3  T²/2;
 *          0     0    0    T³/6  T²/2  T   ]
 *
 * Reference: Li & Jilkov (2003), Sec. III-B
 * Complexity: O(state_dim²)
 *
 * @param dim Motion dimensionality
 * @param T   Sampling interval (seconds)
 * @param q   Process noise PSD for jerk
 */
motion_model_t *motion_create_ca(motion_dim_t dim, double T, double q);

/**
 * Allocate and construct a Coordinated Turn (CT) model.
 *
 * Continuous-time (2D known turn rate ω):
 *   ẋ(t) = v_x(t)
 *   v̇_x(t) = −ω·v_y(t)
 *   ẏ(t) = v_y(t)
 *   v̇_y(t) = ω·v_x(t)
 *
 * State vector: [x, ẋ, y, ẏ, ω] (or [x, ẋ, y, ẏ] for known ω)
 *
 * For known ω:
 * F = [1 sin(ωT)/ω 0 −(1−cos(ωT))/ω;
 *      0 cos(ωT)    0 −sin(ωT)       ;
 *      0 (1−cos(ωT))/ω 1 sin(ωT)/ω   ;
 *      0 sin(ωT)    0 cos(ωT)        ]
 *
 * Process noise adds white acceleration:
 * Q = q · [T³/3 T²/2 0   0  ;
 *          T²/2 T    0   0  ;
 *          0    0    T³/3 T²/2;
 *          0    0    T²/2 T   ]
 *
 * Reference: Li & Jilkov (2003), Sec. IV; Bar-Shalom (2011), Sec. 11.7.2
 * Complexity: O(state_dim²)
 *
 * @param omega Turn rate (rad/s), positive = left turn
 * @param T     Sampling interval (seconds)
 * @param q     Process noise PSD
 */
motion_model_t *motion_create_ct(double omega, double T, double q);

/**
 * Allocate and construct a Singer maneuver model.
 *
 * The target acceleration a(t) is modeled as exponentially correlated
 * (first-order Markov) process:
 *
 *   R_a(τ) = E[a(t)·a(t+τ)] = σ_m² · exp(−α·|τ|)
 *
 * where σ_m² is the maneuver variance and α = 1/τ_m is the reciprocal
 * of the maneuver time constant.
 *
 * State vector (1D): [x, ẋ, ẍ]
 * Discrete-time transition:
 *   F = [1 T    (αT − 1 + e^{−αT})/α²;
 *        0 1    (1 − e^{−αT})/α;
 *        0 0     e^{−αT}                  ]
 *
 * Process noise Q is given by:
 *   Q = 2ασ_m²·q_singer(α, T)  (11-element matrix with α/T correlations)
 *
 * Key parameters:
 *   - α → 0: CA model (constant acceleration, long correlation)
 *   - α → ∞: CV model (white acceleration, zero correlation)
 *   - 0 < α < ∞: realistic maneuver between the extremes
 *
 * Reference: Singer, R.A. "Estimating Optimal Tracking Filter Performance
 *            for Piloted Maneuvering Targets" IEEE T-AES (1970)
 * Complexity: O(state_dim²)
 *
 * @param dim     Dimensionality
 * @param T       Sampling interval
 * @param alpha   Maneuver frequency (rad/s), = 1/τ_m
 * @param sigma_m Maneuver standard deviation (m/s²)
 */
motion_model_t *motion_create_singer(motion_dim_t dim, double T,
                                      double alpha, double sigma_m);

/**
 * Allocate and construct a Wiener Process Acceleration (WPA) model.
 *
 * Also called the "white-noise jerk" model:
 *   Continuous-time: a(t) = Wiener process (integrated white noise)
 *
 * This is a special case of the Singer model when α → 0 and σ_m² → ∞
 * such that 2ασ_m² → q (finite white noise spectrum).
 *
 * Reference: Bar-Shalom (2011), Sec. 6.3
 * Complexity: O(state_dim²)
 */
motion_model_t *motion_create_wpa(motion_dim_t dim, double T, double q);

/**
 * Allocate a Brownian (random walk) model.
 *
 * State: [x] (1D), [x, y] (2D), [x, y, z] (3D)
 * F = I (identity)
 * Q = q·T·I
 *
 * Used for stationary or very slowly moving targets.
 * Complexity: O(state_dim²)
 */
motion_model_t *motion_create_brownian(motion_dim_t dim, double T, double q);

/**
 * Free motion model and all associated matrices.
 * Complexity: O(1)
 */
void motion_free(motion_model_t *model);

/* ============================================================================
 * L3 — State propagation using motion model
 * ============================================================================
 */

/**
 * Propagate state forward by one time step: x_next = F·x.
 * Complexity: O(state_dim²)
 */
void motion_predict_state(const motion_model_t *model,
                           const double *x, double *x_next);

/**
 * Propagate covariance forward: P_next = F·P·Fᵀ + Q.
 * Complexity: O(state_dim³)
 */
void motion_predict_covariance(const motion_model_t *model,
                                const double *P, double *P_next);

/**
 * Combined state and covariance prediction.
 * Complexity: O(state_dim³)
 */
void motion_predict(const motion_model_t *model,
                     const double *x, const double *P,
                     double *x_next, double *P_next);

/* ============================================================================
 * L3 — Measurement matrix construction
 * ============================================================================
 */

/**
 * Build measurement matrix H for a given motion model and measurement type.
 *
 * For Cartesian position measurement of CV 2D model:
 *   H = [1 0 0 0; 0 0 1 0] (measures position only)
 *
 * For full-state measurement: H = I
 * For velocity-only: H selects velocity rows
 *
 * Complexity: O(meas_dim × state_dim)
 *
 * @param model      Motion model
 * @param meas_type  What is being measured
 * @param H          Output matrix [meas_dim × state_dim], pre-allocated
 * @param meas_dim   Measurement dimension
 * @return 0 on success, -1 on incompatibility
 */
int motion_build_H(const motion_model_t *model, measurement_type_t meas_type,
                    double *H, int meas_dim);

/* ============================================================================
 * L5 — IMM model set construction
 * ============================================================================
 */

/**
 * Build a standard IMM model set with mixing matrix.
 *
 * Standard 3-model set for air targets:
 *   Model 1: CV (non-maneuver)
 *   Model 2: CT with ω > 0 (left turn)
 *   Model 3: CT with ω < 0 (right turn)
 *
 * Markov transition probabilities:
 *   p_ii = 0.95 (stay probability)
 *   p_ij = 0.025 (transition to each other model)
 *
 * Complexity: O(r·n²) for r models
 *
 * @param models       Output: array of r model pointers
 * @param mixing       Output: Markov transition matrix [r×r], row-stochastic
 * @param max_models   Input: max models to create / output: actual count
 * @param T            Sampling interval
 * @param q_cv         CV process noise
 * @param q_ct         CT process noise
 * @param omega        CT turn rate magnitude
 * @return Number of models created
 */
int imm_build_standard_3model(motion_model_t **models, double *mixing,
                               int *max_models, double T,
                               double q_cv, double q_ct, double omega);

/**
 * Build a 2-model IMM for road vehicles (stop & go):
 *   Model 1: CV (moving)
 *   Model 2: Stopped (zero velocity, small position noise)
 *
 * Complexity: O(2·n²)
 */
int imm_build_vehicle_2model(motion_model_t **models, double *mixing,
                              double T, double q_moving, double q_stopped);

/* ============================================================================
 * L2 — Discretization of continuous-time models
 * ============================================================================
 */

/**
 * Discretize a continuous-time linear system:
 *
 *   ẋ(t) = A·x(t) + G·w(t),  w(t) ~ N(0, Q_c)
 *
 * to discrete-time:
 *
 *   x_{k+1} = F·x_k + w_k,  w_k ~ N(0, Q_d)
 *
 * where:
 *   F = exp(A·T) = Σ_{i=0}^{∞} (A·T)^i / i!
 *   Q_d = ∫_0^T exp(A·τ)·G·Q_c·Gᵀ·exp(Aᵀ·τ) dτ
 *
 * Truncate series at sufficient accuracy.
 *
 * Reference: Bar-Shalom (2011), Appendix D
 * Complexity: O(state_dim³ · series_terms)
 *
 * @param A          Continuous-time dynamics [state_dim × state_dim]
 * @param G          Noise gain [state_dim × noise_dim]
 * @param Q_c        Continuous process noise PSD [noise_dim × noise_dim]
 * @param state_dim  State dimension
 * @param noise_dim  Noise dimension
 * @param T          Sampling interval
 * @param F          Output discrete transition [state_dim × state_dim]
 * @param Q_d        Output discrete noise covariance [state_dim × state_dim]
 * @param terms      Number of series terms (≥8 recommended)
 */
void motion_discretize(const double *A, const double *G, const double *Q_c,
                        int state_dim, int noise_dim, double T,
                        double *F, double *Q_d, int terms);

/**
 * Compute matrix exponential via truncated series.
 * Complexity: O(state_dim³ · terms)
 */
void matrix_exponential(const double *A, double *expA, int n, int terms);

#ifdef __cplusplus
}
#endif

#endif /* MOTION_MODELS_H */
