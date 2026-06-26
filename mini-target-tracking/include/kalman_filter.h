/**
 * kalman_filter.h — Kalman Filter Variants for Target Tracking
 *
 * Covers: L3 (Mathematical Structures) — Gaussian state-space, Riccati equation
 *         L4 (Fundamental Laws) — Kalman optimality, Cramér-Rao lower bound
 *         L5 (Algorithms/Methods) — KF, EKF, UKF, SR-KF, Information Filter
 *         L6 (Canonical Problems) — Filter initialization, tuning, divergence
 *
 * References:
 *   - Kalman, R.E. "A New Approach to Linear Filtering and Prediction" (1960)
 *   - Bar-Shalom, Willett, Tian "Tracking and Data Fusion" (2011), Ch. 3-5
 *   - Julier, Uhlmann "Unscented Filtering and Nonlinear Estimation" (2004)
 *   - Bierman, G.J. "Factorization Methods for Discrete Sequential Estimation" (1977)
 *   - Maybeck, P.S. "Stochastic Models, Estimation and Control" Vol. 1 (1979)
 *
 * Curriculum:
 *   - MIT 6.003 (Signal Processing), Stanford EE102A, Berkeley EE123 (DSP)
 *   - Michigan EECS 351 (DSP), Georgia Tech ECE 4270
 *   - ETH 227-0427 (Signal Processing), 清华 (信号与系统)
 */

#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"

/* ============================================================================
 * L1 — Filter type enumeration
 * ============================================================================
 */

/** Kalman filter variants supported */
typedef enum {
    KF_LINEAR      = 0,  /**< Standard linear discrete-time KF */
    KF_EXTENDED    = 1,  /**< Extended Kalman Filter (Jacobian linearization) */
    KF_UNSCENTED   = 2,  /**< Unscented Kalman Filter (sigma-point) */
    KF_SQUARE_ROOT = 3,  /**< Square-Root KF (UDU or Cholesky factors) */
    KF_INFO        = 4,  /**< Information Filter (inverse covariance form) */
    KF_IMM         = 5,  /**< Interacting Multiple Model (uses multiple sub-filters) */
    KF_ADAPTIVE    = 6   /**< Adaptive KF with online Q/R estimation */
} kf_type_t;

/* ============================================================================
 * L1 — Kalman filter structure
 * ============================================================================
 */

/** Linear Kalman Filter instance */
typedef struct {
    kf_type_t type;        /**< Filter variant */
    int       state_dim;   /**< State vector dimension n */
    int       meas_dim;    /**< Measurement vector dimension m */
    int       control_dim; /**< Control input dimension c (0 = none) */

    /* State estimate */
    double *x;             /**< State estimate [n] */
    double *P;             /**< Error covariance [n × n], row-major */

    /* Predicted state for current step */
    double *x_pred;        /**< Predicted state [n] */
    double *P_pred;        /**< Predicted covariance [n × n] */

    /* Innovation */
    double *nu;            /**< Innovation/residual [m] */
    double *S;             /**< Innovation covariance [m × m] */

    /* Kalman gain */
    double *K;             /**< Kalman gain [n × m] */

    /* Model matrices */
    double *F;             /**< State transition [n × n] */
    double *H;             /**< Measurement matrix [m × n] */
    double *Q;             /**< Process noise covariance [n × n] */
    double *R;             /**< Measurement noise covariance [m × m] */

    /* Control */
    double *B;             /**< Control input matrix [n × c] (NULL if no control) */
    double *u;             /**< Control input [c] (NULL if no control) */

    /* Workspace for intermediate computations */
    double *work_n;        /**< Scratch buffer size n² */
    double *work_m;        /**< Scratch buffer size m² */

    /* Diagnostics */
    double  log_likelihood;  /**< Negative log-likelihood of measurements */
    double  nees;             /**< Normalized Estimation Error Squared */
    int     steps;            /**< Number of filter steps executed */
    int     diverged;         /**< Divergence flag (1 = diverged) */
    double  consistency_ratio; /**< NEES / n for consistency check */

    /* Adaptive filter fields */
    double *Q_est;          /**< Estimated process noise [n × n] (adaptive) */
    double *R_est;          /**< Estimated measurement noise [m × m] (adaptive) */
    double  forgetting_factor; /**< Exponential forgetting factor for adaptation */

    /* UKF-specific fields */
    double *sigma_points;   /**< Sigma points [(2n+1) × n] */
    double *sigma_weights;  /**< Sigma weights [2n+1] */
    double *sigma_pred;     /**< Predicted sigma points [(2n+1) × n] */
    double *sigma_meas;     /**< Measurement sigma points [(2n+1) × m] */
    int     num_sigma;       /**< Number of sigma points (2n+1) */
    double  alpha_ukf;       /**< UKF alpha parameter */
    double  beta_ukf;        /**< UKF beta parameter */
    double  kappa_ukf;       /**< UKF kappa parameter */
    double  lambda_ukf;      /**< UKF lambda = alpha²(n+κ) − n */

    /* Square-Root specific fields */
    double *S_chol;         /**< Cholesky factor of P [n × n] lower-triangular */
    int     use_ud;          /**< 1 = UDU factor, 0 = Cholesky */

    /* Function pointers for nonlinear models (EKF/UKF) */
    /** State transition function (nonlinear): x_{k+1} = f(x_k, u_k) */
    void (*f_nonlinear)(const double *x, const double *u, double *x_next,
                         int state_dim, int control_dim);
    /** Measurement function (nonlinear): z_k = h(x_k) */
    void (*h_nonlinear)(const double *x, double *z_pred,
                         int state_dim, int meas_dim);
    /** Linearization: F = ∂f/∂x at current estimate */
    void (*compute_F_jacobian)(const double *x, const double *u,
                                double *F, int state_dim, int control_dim);
    /** Linearization: H = ∂h/∂x at predicted state */
    void (*compute_H_jacobian)(const double *x_pred,
                                double *H, int state_dim, int meas_dim);

    /* Multi-model (IMM) sub-filter pointers */
    void **sub_filters;     /**< Array of sub-filter pointers */
    int    num_models;      /**< Number of IMM models */
    double *model_probs;    /**< Model probabilities [num_models] */
    double *mixing_matrix;  /**< Markov transition matrix [num_models²] */
} kf_t;

/* ============================================================================
 * L3 — Mathematical Structures: Gaussian state-space model definition
 * ============================================================================
 */

/**
 * Allocate and initialize a Kalman filter.
 *
 * The state-space model:
 *   x_{k+1} = F·x_k + B·u_k + w_k,   w_k ~ N(0, Q)
 *   z_k     = H·x_k + v_k,           v_k ~ N(0, R)
 *
 * Complexity: O(n² + m²) allocation
 *
 * @param state_dim   State dimension n
 * @param meas_dim    Measurement dimension m
 * @param control_dim Control dimension c (0 if none)
 * @return Initialized filter, or NULL on allocation failure
 */
kf_t *kf_alloc(int state_dim, int meas_dim, int control_dim);

/**
 * Deallocate filter and all internal buffers.
 * Complexity: O(1) free
 */
void kf_free(kf_t *kf);

/**
 * Set the linear model matrices by copying provided values.
 * F, H, Q, R are mandatory. B, u are optional (set NULL if unused).
 * Complexity: O(n² + m² + n·m)
 */
void kf_set_model(kf_t *kf,
                   const double *F, const double *H,
                   const double *Q, const double *R,
                   const double *B, const double *u);

/**
 * Initialize the state estimate.
 * Complexity: O(n²)
 */
void kf_set_state(kf_t *kf, const double *x0, const double *P0);

/* ============================================================================
 * L4 — Fundamental Laws: Kalman optimality
 * ============================================================================
 */

/**
 * Kalman filter prediction step.
 *
 * x_pred = F·x
 * P_pred = F·P·Fᵀ + Q
 *
 * For EKF: x_pred = f(x), P_pred = F·P·Fᵀ + Q (where F = ∂f/∂x)
 *
 * This implements the minimum-variance unbiased predictor under Gaussian noise.
 * Theorem (Kalman 1960): Among all linear estimators, the KF minimizes
 * E[||x - x̂||²] for linear Gaussian systems.
 *
 * Complexity: O(n³) for matrix multiplications
 *
 * @return 0 on success, -1 on numerical failure
 */
int kf_predict(kf_t *kf);

/**
 * Kalman filter update/correction step.
 *
 * ν = z − H·x_pred                                  (innovation)
 * S = H·P_pred·Hᵀ + R                                (innovation covariance)
 * K = P_pred·Hᵀ·S⁻¹                                  (Kalman gain)
 * x = x_pred + K·ν                                   (updated state)
 * P = (I − K·H)·P_pred                               (updated covariance)
 *
 * The gain K minimizes tr(P) — the minimum-variance optimal gain.
 * This is the discrete-time algebraic Riccati equation solution.
 *
 * Complexity: O(n³ + m³) for matrix inversions
 *
 * @param z Measurement vector [m]
 * @return 0 on success, -1 on failure
 */
int kf_update(kf_t *kf, const double *z);

/**
 * Combined predict + update in one call.
 * Complexity: O(n³ + m³)
 */
int kf_step(kf_t *kf, const double *z);

/* ============================================================================
 * L5 — Algorithms: EKF, UKF, SR-KF, Info Filter
 * ============================================================================
 */

/**
 * Extended Kalman Filter update — uses Jacobian linearization of h(x).
 *
 * H = ∂h/∂x evaluated at x_pred.
 * Then standard KF update with this H.
 *
 * The EKF approximates the nonlinear system by first-order Taylor expansion.
 * Consistency is not guaranteed — EKF may diverge for highly nonlinear systems.
 *
 * Reference: Jazwinski, A.H. "Stochastic Processes and Filtering Theory" (1970)
 * Complexity: O(n³ + m³)
 *
 * @param kf Filter with h_nonlinear and compute_H_jacobian function pointers set
 * @param z  Measurement vector [m]
 * @return 0 on success, -1 on failure
 */
int ekf_update(kf_t *kf, const double *z);

/**
 * Extended Kalman Filter predict — uses nonlinear state transition.
 *
 * x_pred = f(x, u)
 * P_pred = F·P·Fᵀ + Q, F = ∂f/∂x at current state
 *
 * Complexity: O(n³)
 */
int ekf_predict(kf_t *kf);

/**
 * Initialize UKF parameters.
 *
 * Recommended values:
 *   α = 1.0 (10⁻⁴ ≤ α ≤ 1), β = 2.0 (Gaussian optimal), κ = 3 − n
 *
 * Reference: Wan, E.A., van der Merwe, R. "The Unscented Kalman Filter" (2001)
 * Complexity: O(1)
 *
 * @param kf    Already allocated KF (type = KF_UNSCENTED)
 * @param alpha Spread of sigma points
 * @param beta  Prior knowledge of distribution (2 = Gaussian optimal)
 * @param kappa Secondary scaling parameter
 */
void ukf_init_params(kf_t *kf, double alpha, double beta, double kappa);

/**
 * Generate sigma points from mean and covariance.
 *
 * χ₀ = x̄
 * χᵢ = x̄ + √((n+λ)P)ᵢ    for i = 1,...,n
 * χᵢ₊ₙ = x̄ − √((n+λ)P)ᵢ  for i = 1,...,n
 *
 * Weights: W₀ᵐ = λ/(n+λ), W₀ᶜ = λ/(n+λ) + (1−α²+β)
 *          Wᵢᵐ = Wᵢᶜ = 1/(2(n+λ)) for i > 0
 *
 * Complexity: O(n³) for matrix square root (Cholesky)
 *
 * @param kf UKF with initialized state (x, P) and UKF params
 * @return 0 on success, -1 on failure
 */
int ukf_generate_sigma_points(kf_t *kf);

/**
 * Unscented Kalman Filter predict step.
 *
 * Transform sigma points through nonlinear f():
 *   χ_predᵢ = f(χᵢ, u_k)
 * Then compute predicted mean and covariance from transformed points.
 *
 * Complexity: O(n³ + (2n+1)·cost(f))
 */
int ukf_predict(kf_t *kf);

/**
 * Unscented Kalman Filter update step.
 *
 * Transform predicted sigma points through h():
 *   Zᵢ = h(χ_predᵢ)
 * Compute measurement prediction, innovation covariance, and cross-covariance
 * from empirical sigma-point statistics.
 *
 * Complexity: O(n³ + m³ + (2n+1)·cost(h))
 */
int ukf_update(kf_t *kf, const double *z);

/**
 * Square-Root Kalman filter via Potter's algorithm (scalar measurement update).
 *
 * Updates the Cholesky factor S (P = S·Sᵀ) directly, avoiding explicit P
 * and ensuring numerical stability for ill-conditioned problems.
 *
 * Reference: Potter, J.E., Stern, R.G. "Statistical filtering of space
 *            navigation measurements" (1963); Bierman (1977)
 * Complexity: O(n²) per scalar measurement
 *
 * @param kf  Filter with S_chol set as Cholesky factor of P
 * @param H_row Single row of H for scalar measurement
 * @param z_scalar Scalar measurement
 * @param R_scalar Measurement noise variance
 * @return 0 on success, -1 on failure
 */
int srkf_scalar_update(kf_t *kf, const double *H_row,
                        double z_scalar, double R_scalar);

/**
 * Information Filter — operates on information matrix Y = P⁻¹ and
 * information state ŷ = P⁻¹x.
 *
 * Prediction: Y_{k|k-1} = (F·Y⁻¹·Fᵀ + Q)⁻¹
 * Update:     Y_k = Y_{k|k-1} + Hᵀ·R⁻¹·H
 *             ŷ_k = ŷ_{k|k-1} + Hᵀ·R⁻¹·z
 *
 * Key advantage: infinite initial uncertainty (P = ∞ ⇒ Y = 0)
 * does not require special initialization.
 *
 * Reference: Maybeck, P.S. (1979), Vol. 1, Ch. 5
 * Complexity: O(n³ + m³)
 *
 * @param Y_info Information matrix Y [n × n], updated in place
 * @param y_info Information state ŷ [n], updated in place
 * @return 0 on success
 */
int info_filter_update(double *Y_info, double *y_info,
                        const double *H, const double *R_inv,
                        const double *z, int n, int m);

/* ============================================================================
 * L5 — Adaptive Kalman filtering (innovation-based)
 * ============================================================================
 */

/**
 * Adaptively estimate measurement noise covariance R using innovation sequence.
 *
 * R̂_k = (1−α)R̂_{k-1} + α(ν_k ν_kᵀ − H·P_pred·Hᵀ)
 *
 * Where α is the forgetting factor. The subtraction term removes the
 * contribution of state prediction uncertainty from the innovation covariance.
 *
 * Reference: Mehra, R.K. "On the identification of variances and adaptive
 *            Kalman filtering" IEEE TAC (1970)
 * Complexity: O(m² + n·m²)
 *
 * @param kf Kalman filter with forgetting_factor set
 */
void kf_adapt_R(kf_t *kf);

/**
 * Adaptively estimate process noise covariance Q using innovation.
 *
 * Q̂ = K·ν·νᵀ·Kᵀ  (approximation based on smoothed state changes)
 *
 * A more robust approach uses the Myers-Tapley scheme.
 * Complexity: O(n²·m)
 */
void kf_adapt_Q(kf_t *kf);

/* ============================================================================
 * L5 — Extended algorithms: IMM base, filter bank
 * ============================================================================
 */

/**
 * Interacting Multiple Model (IMM) mixing step.
 *
 * Computes mixed initial conditions for each sub-filter:
 *   x̂₀ⱼ = Σᵢ x̂ᵢ · μ_{i|j}
 *   P₀ⱼ = Σᵢ [Pᵢ + (x̂ᵢ − x̂₀ⱼ)(x̂ᵢ − x̂₀ⱼ)ᵀ] · μ_{i|j}
 *
 * where μ_{i|j} = p_{ij}·μᵢ / Σₖ p_{kj}·μₖ are the mixing probabilities.
 *
 * Reference: Blom, H.A.P., Bar-Shalom, Y. "IMM Estimator" IEEE TAC (1988)
 * Complexity: O(r²·n²) where r = number of models
 *
 * @param kf IMM filter with sub_filters and mixing_matrix set
 * @return 0 on success
 */
int imm_mix(kf_t *kf);

/**
 * Update IMM model probabilities using likelihood of each sub-filter.
 *
 * μⱼ = (1/c) · Λⱼ · Σᵢ p_{ij}·μᵢ
 *
 * where Λⱼ = N(νⱼ; 0, Sⱼ) is the Gaussian likelihood.
 *
 * Complexity: O(r² + r·m³)
 */
int imm_update_probs(kf_t *kf);

/**
 * Compute overall IMM estimate as weighted sum of sub-filter estimates.
 *
 * x̂ = Σⱼ μⱼ · x̂ⱼ
 * P  = Σⱼ μⱼ · [Pⱼ + (x̂ⱼ − x̂)(x̂ⱼ − x̂)ᵀ]
 *
 * Complexity: O(r·n²)
 */
int imm_combine(kf_t *kf);

/* ============================================================================
 * L4 — Filter consistency and diagnostics
 * ============================================================================
 */

/**
 * Compute Normalized Innovation Squared (NIS).
 *
 * ε_ν = νᵀ·S⁻¹·ν
 *
 * Under consistent filter, ε_ν ~ χ²(m).
 *
 * Complexity: O(m³)
 */
double kf_nis(const kf_t *kf);

/**
 * Compute Normalized Estimation Error Squared (NEES).
 *
 * ε_x = (x_true − x̂)ᵀ·P⁻¹·(x_true − x̂)
 *
 * Under consistent filter, ε_x ~ χ²(n).
 * Requires ground truth; for simulation validation only.
 *
 * Complexity: O(n³)
 */
double kf_nees(const kf_t *kf, const double *x_true);

/**
 * Chi-squared test for filter consistency.
 * Returns 1 if consistent at given significance level α, 0 otherwise.
 *
 * H₀: filter is consistent
 * Reject H₀ if ε exceeds χ²_{1−α}(dof).
 *
 * Complexity: O(1)
 */
int kf_consistency_test(double epsilon, int dof, double alpha);

/**
 * Detect filter divergence by monitoring innovation magnitude.
 * Returns 1 if divergence detected.
 *
 * Criterion: NIS > threshold for k consecutive steps.
 *
 * Complexity: O(1)
 */
int kf_divergence_detect(kf_t *kf, double threshold, int consecutive);

/**
 * Get the negative log-likelihood of measurements for model selection.
 *
 * −log p(z₁,...,zₖ) = ½ Σ [m·log(2π) + log|Sₖ| + νₖᵀ·Sₖ⁻¹·νₖ]
 *
 * Complexity: O(1) — uses running accumulator
 */
double kf_neg_log_likelihood(const kf_t *kf);

#ifdef __cplusplus
}
#endif

#endif /* KALMAN_FILTER_H */
