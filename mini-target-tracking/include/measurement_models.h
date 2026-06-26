/**
 * measurement_models.h — Radar Measurement Models and Coordinate Conversions
 *
 * Covers: L1 (Definitions) — Range/bearing/Doppler, polar-to-Cartesian
 *         L2 (Core Concepts) — Converted measurements, debiasing
 *         L3 (Mathematical Structures) — Jacobians, coordinate transforms
 *         L4 (Fundamental Laws) — CRLB for radar, resolution bounds
 *         L5 (Algorithms/Methods) — CMKF, debiased conversion, unscented transform
 *
 * References:
 *   - Lerro, D., Bar-Shalom, Y. "Tracking with debiased consistent converted
 *     measurements versus EKF" IEEE T-AES (1993)
 *   - Longbin, M. et al. "Unbiased converted measurements for tracking" (1998)
 *   - Suchomski, P. "Explicit expressions for debiased statistics of 3D
 *     converted measurements" IEEE T-AES (1999)
 *   - Richards, Scheer, Holm "Principles of Modern Radar" (2010), Ch. 18
 *
 * Curriculum:
 *   - MIT 6.630 (EM Waves), Stanford EE359, Michigan EECS 411
 *   - Georgia Tech ECE 6350, TU Munich High-Frequency Eng
 */

#ifndef MEASUREMENT_MODELS_H
#define MEASUREMENT_MODELS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"

/* ============================================================================
 * L1 — Radar measurement parameters
 * ============================================================================
 */

/** Radar sensor parameters for measurement model */
typedef struct {
    /* Position */
    double sensor_x;        /**< Sensor x-position in global frame (m) */
    double sensor_y;        /**< Sensor y-position in global frame (m) */
    double sensor_z;        /**< Sensor z-position in global frame (m) */
    double yaw;             /**< Sensor yaw angle (rad) */
    double pitch;           /**< Sensor pitch angle (rad) */
    double roll;            /**< Sensor roll angle (rad) */

    /* Measurement noise standard deviations */
    double sigma_range;     /**< Range measurement noise std (m) */
    double sigma_bearing;   /**< Bearing measurement noise std (rad) */
    double sigma_elevation; /**< Elevation measurement noise std (rad) */
    double sigma_doppler;   /**< Doppler/range-rate noise std (m/s) */

    /* Detection */
    double P_D;             /**< Probability of detection */
    double P_FA;            /**< False alarm probability */
    double SNR_dB;          /**< Signal-to-noise ratio (dB) */
    double range_min;       /**< Minimum detection range (m) */
    double range_max;       /**< Maximum unambiguous range (m) */
    double v_max;           /**< Maximum unambiguous velocity (m/s) */

    /* Resolution */
    double range_resolution;   /**< Range resolution (m) */
    double bearing_resolution; /**< Bearing resolution (rad) */
    double doppler_resolution; /**< Doppler resolution (m/s) */

    /* Beam */
    double beamwidth_az;    /**< Azimuth beamwidth (rad) */
    double beamwidth_el;    /**< Elevation beamwidth (rad) */
} radar_sensor_t;

/* ============================================================================
 * L1+L2 — Coordinate conversion functions
 * ============================================================================
 */

/**
 * Convert polar coordinates (range, bearing) to Cartesian (x, y).
 *
 * x = r·cos(θ),  y = r·sin(θ)
 *
 * Complexity: O(1)
 *
 * @param r     Range (meters)
 * @param theta Bearing angle (radians, 0 = +x axis)
 * @param x     Output x-coordinate
 * @param y     Output y-coordinate
 */
void polar_to_cartesian_2d(double r, double theta, double *x, double *y);

/**
 * Convert Cartesian (x, y) to polar (range, bearing).
 *
 * r = √(x² + y²), θ = atan2(y, x)
 *
 * Complexity: O(1)
 */
void cartesian_to_polar_2d(double x, double y, double *r, double *theta);

/**
 * Convert spherical (range, azimuth, elevation) to Cartesian 3D.
 *
 * x = r·cos(el)·cos(az), y = r·cos(el)·sin(az), z = r·sin(el)
 *
 * Complexity: O(1)
 */
void spherical_to_cartesian_3d(double r, double az, double el,
                                double *x, double *y, double *z);

/**
 * Convert Cartesian 3D to spherical.
 *
 * r = √(x²+y²+z²), az = atan2(y, x), el = atan2(z, √(x²+y²))
 *
 * Complexity: O(1)
 */
void cartesian_to_spherical_3d(double x, double y, double z,
                                double *r, double *az, double *el);

/* ============================================================================
 * L2+L3 — Measurement Jacobians for EKF
 * ============================================================================
 */

/**
 * Compute the Jacobian of the 2D polar-to-Cartesian measurement function.
 *
 * h(x) = [√((x−sx)² + (y−sy)²); atan2(y−sy, x−sx)]
 *
 * H = ∂h/∂x = [∂r/∂x ∂r/∂y ∂r/∂vx ∂r/∂vy;
 *               ∂θ/∂x ∂θ/∂y ∂θ/∂vx ∂θ/∂vy]
 *
 * where:
 *   ∂r/∂x = (x−sx)/d,  ∂r/∂y = (y−sy)/d,  ∂r/∂vx = 0,  ∂r/∂vy = 0
 *   ∂θ/∂x = −(y−sy)/d², ∂θ/∂y = (x−sx)/d², ∂θ/∂vx = 0, ∂θ/∂vy = 0
 *   d = √((x−sx)² + (y−sy)²)
 *
 * Complexity: O(state_dim · meas_dim)
 *
 * @param state     State vector [state_dim] (positions at appropriate indices)
 * @param sensor    Sensor parameters
 * @param H         Output Jacobian [2 × state_dim]
 * @param state_dim State dimension
 * @param pos_indices Array of position indices in state vector [2] (e.g., [0, 3] for CV)
 */
void measurement_jacobian_polar2d(const double *state, const radar_sensor_t *sensor,
                                   double *H, int state_dim, const int *pos_indices);

/**
 * Compute the Jacobian of the 3D spherical measurement function.
 *
 * h(x) = [r; az; el] with r, az, el as defined in spherical_to_cartesian_3d.
 *
 * Complexity: O(state_dim · 3)
 */
void measurement_jacobian_spherical(const double *state, const radar_sensor_t *sensor,
                                     double *H, int state_dim, const int *pos_indices);

/**
 * Compute the Jacobian of the Doppler measurement (range-rate projection).
 *
 * h_doppler(x) = ((x−sx)·vx + (y−sy)·vy + (z−sz)·vz) / d
 *
 * where d = √((x−sx)² + (y−sy)² + (z−sz)²)
 *
 * This is the radial component of the velocity vector.
 *
 * Complexity: O(state_dim)
 *
 * @param state     State vector with position and velocity
 * @param sensor    Sensor parameters
 * @param H_dop     Output Jacobian [1 × state_dim]
 * @param state_dim State dimension
 * @param pos_indices Position indices [3]
 * @param vel_indices Velocity indices [3]
 */
void measurement_jacobian_doppler(const double *state, const radar_sensor_t *sensor,
                                   double *H_dop, int state_dim,
                                   const int *pos_indices, const int *vel_indices);

/* ============================================================================
 * L2+L5 — Measurement prediction functions
 * ============================================================================
 */

/**
 * Predict measurement from state (2D polar).
 *
 * z_pred = [√((x−sx)² + (y−sy)²); atan2(y−sy, x−sx)]
 *
 * Complexity: O(1)
 */
void measurement_predict_polar2d(const double *state, const radar_sensor_t *sensor,
                                  double *z_pred, const int *pos_indices);

/**
 * Predict measurement from state (3D spherical).
 * Complexity: O(1)
 */
void measurement_predict_spherical(const double *state, const radar_sensor_t *sensor,
                                    double *z_pred,
                                    const int *pos_indices);

/**
 * Predict range-rate (Doppler) measurement.
 * Complexity: O(1)
 */
double measurement_predict_doppler(const double *state, const radar_sensor_t *sensor,
                                    const int *pos_indices, const int *vel_indices);

/* ============================================================================
 * L5 — Converted Measurement Kalman Filter (CMKF)
 * ============================================================================
 */

/**
 * Convert a polar measurement to Cartesian coordinates with debiasing.
 *
 * Standard conversion:
 *   x_c = r_m·cos(θ_m),  y_c = r_m·sin(θ_m)
 *
 * This is biased because:
 *   E[x_c] = r·cos(θ)·exp(−σ²_θ/2) ≠ x_true
 *   E[y_c] = r·sin(θ)·exp(−σ²_θ/2) ≠ y_true
 *
 * Debiased conversion (Lerro & Bar-Shalom 1993):
 *   x_d = r_m·cos(θ_m) / exp(−σ²_θ/2)
 *   y_d = r_m·sin(θ_m) / exp(−σ²_θ/2)
 *
 * Or equivalently (Longbin 1998):
 *   x_d = r_m·cos(θ_m) − μ_x,  y_d = r_m·sin(θ_m) − μ_y
 *
 * where μ_x, μ_y are the bias correction terms.
 *
 * Complexity: O(1)
 *
 * @param r_meas    Measured range (m)
 * @param theta_meas Measured bearing (rad)
 * @param sigma_r   Range noise std
 * @param sigma_theta Bearing noise std
 * @param x_conv    Output debiased x
 * @param y_conv    Output debiased y
 */
void polar_to_cartesian_debiased(double r_meas, double theta_meas,
                                  double sigma_r, double sigma_theta,
                                  double *x_conv, double *y_conv);

/**
 * Compute the converted measurement error covariance (Cartesian).
 *
 * For the debiased conversion (Longbin 1998), elements of R_c:
 *
 * R_xx = r_m²·σ²_θ·sin²(θ_m) + σ²_r·cos²(θ_m) − bias_xx_term
 * R_yy = r_m²·σ²_θ·cos²(θ_m) + σ²_r·sin²(θ_m) − bias_yy_term
 * R_xy = (σ²_r − r_m²·σ²_θ)·sin(θ_m)·cos(θ_m) − bias_xy_term
 *
 * Reference: Longbin et al. (1998)
 * Complexity: O(1)
 *
 * @param r_meas      Measured range
 * @param theta_meas  Measured bearing
 * @param sigma_r     Range noise std
 * @param sigma_theta Bearing noise std
 * @param R_conv      Output covariance [2×2] = [R_xx, R_xy; R_xy, R_yy] (row-major)
 */
void converted_measurement_covariance(double r_meas, double theta_meas,
                                       double sigma_r, double sigma_theta,
                                       double *R_conv);

/**
 * Converted Measurement Kalman Filter — one full step.
 *
 * Converts polar measurement to Cartesian, computes converted covariance,
 * then runs standard KF (linear update with Cartesian measurement).
 *
 * This avoids EKF linearization errors for range-bearing measurements
 * at the cost of approximate noise covariance.
 *
 * Reference: Lerro & Bar-Shalom (1993)
 * Complexity: O(n³) for KF step
 *
 * @param x       State estimate [state_dim], updated in place
 * @param P       Covariance [state_dim²], updated in place
 * @param F       Transition matrix [state_dim²]
 * @param Q       Process noise [state_dim²]
 * @param r_meas  Range measurement
 * @param theta_meas Bearing measurement
 * @param sigma_r Range std
 * @param sigma_theta Bearing std
 * @param H_cart  Cartesian measurement matrix (e.g., H=[1 0 0 0; 0 0 1 0])
 * @param state_dim State dimension
 * @param meas_dim  Measurement dimension (2 for 2D Cartesian)
 * @return 0 on success
 */
int cmkf_step(double *x, double *P,
               const double *F, const double *Q,
               double r_meas, double theta_meas,
               double sigma_r, double sigma_theta,
               const double *H_cart, int state_dim, int meas_dim);

/* ============================================================================
 * L5 — Unscented transform for measurement prediction
 * ============================================================================
 */

/**
 * Propagate state distribution through nonlinear measurement function
 * using the Unscented Transform.
 *
 * Given x ~ N(x̄, P), compute z = h(x) ~ N(ẑ, P_zz)
 * and cross-covariance P_xz.
 *
 * Algorithm:
 *   1. Generate 2n+1 sigma points from (x̄, P)
 *   2. Transform each: Z_i = h(X_i)
 *   3. Compute weighted mean and covariance of transformed points
 *
 * Reference: Julier & Uhlmann (2004)
 * Complexity: O(n³ + (2n+1)·cost(h))
 *
 * @param x_bar     Mean of input distribution [n]
 * @param P         Covariance of input [n×n]
 * @param n         State dimension
 * @param h_func    Nonlinear measurement function
 * @param h_context User context for h_func
 * @param m         Measurement dimension
 * @param z_bar     Output: predicted measurement mean [m]
 * @param P_zz      Output: measurement covariance [m×m]
 * @param P_xz      Output: cross-covariance [n×m]
 * @param alpha     UT alpha parameter
 * @param beta      UT beta parameter
 * @param kappa     UT kappa parameter
 */
void unscented_transform(const double *x_bar, const double *P, int n,
                          void (*h_func)(const double *x, double *z,
                                         int n, int m, void *context),
                          void *h_context, int m,
                          double *z_bar, double *P_zz, double *P_xz,
                          double alpha, double beta, double kappa);

/* ============================================================================
 * L2 — Radar resolution and measurement ambiguity
 * ============================================================================
 */

/**
 * Check range ambiguity: is the target beyond the unambiguous range?
 *
 * R_unamb = c / (2·PRF)  (pulse radar)
 *
 * Returns: 0 = unambiguous, n > 0 = nth ambiguity fold
 * Complexity: O(1)
 */
int range_ambiguity_fold(double range, double R_unamb);

/**
 * Check Doppler (velocity) ambiguity.
 *
 * v_unamb = λ·PRF / 2
 *
 * Returns: folded velocity into [-v_unamb/2, v_unamb/2]
 * Complexity: O(1)
 */
double doppler_unwrap(double v_measured, double v_unamb);

/**
 * Compute Cramér-Rao Lower Bound for range estimation.
 *
 * σ²_r ≥ c² / (8·π²·B²·SNR)
 *
 * where B is the signal bandwidth, SNR is the signal-to-noise ratio.
 * This is the fundamental limit on range estimation accuracy.
 *
 * Reference: Richards (2010), Sec. 16.3
 * Complexity: O(1)
 *
 * @param bandwidth Signal bandwidth (Hz)
 * @param SNR       Linear SNR (not dB)
 * @return CRLB for range variance (m²)
 */
double crlb_range(double bandwidth, double SNR);

/**
 * Compute Cramér-Rao Lower Bound for bearing (angle) estimation.
 *
 * σ²_θ ≥ λ² / (8·π²·D²·SNR·cos²(θ))
 *
 * where D is the antenna aperture, λ is the wavelength.
 *
 * Complexity: O(1)
 */
double crlb_bearing(double wavelength, double aperture, double SNR, double theta);

/**
 * Compute Cramér-Rao Lower Bound for Doppler estimation.
 *
 * σ²_v ≥ λ² / (8·π²·T_cpi²·SNR)
 *
 * where T_cpi is the coherent processing interval.
 *
 * Complexity: O(1)
 */
double crlb_doppler(double wavelength, double T_cpi, double SNR);

#ifdef __cplusplus
}
#endif

#endif /* MEASUREMENT_MODELS_H */
