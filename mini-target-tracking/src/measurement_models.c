/**
 * measurement_models.c — Radar Measurement Models Implementation
 *
 * Implements: Coordinate transformations, Jacobians, debiased conversion,
 *             CMKF, Unscented Transform, CRLB computations.
 *
 * References:
 *   - Lerro & Bar-Shalom (1993), Longbin et al. (1998)
 *   - Julier & Uhlmann (2004)
 *   - Richards, Scheer, Holm "Principles of Modern Radar" (2010)
 */

#include "measurement_models.h"
#include "kalman_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Coordinate conversions
 * ============================================================================
 */

void polar_to_cartesian_2d(double r, double theta, double *x, double *y)
{
    if (!x || !y) return;
    *x = r * cos(theta);
    *y = r * sin(theta);
}

void cartesian_to_polar_2d(double x, double y, double *r, double *theta)
{
    if (!r || !theta) return;
    *r = sqrt(x * x + y * y);
    *theta = atan2(y, x);
}

void spherical_to_cartesian_3d(double r, double az, double el,
                                double *x, double *y, double *z)
{
    if (!x || !y || !z) return;
    double cos_el = cos(el);
    *x = r * cos_el * cos(az);
    *y = r * cos_el * sin(az);
    *z = r * sin(el);
}

void cartesian_to_spherical_3d(double x, double y, double z,
                                double *r, double *az, double *el)
{
    if (!r || !az || !el) return;
    *r = sqrt(x * x + y * y + z * z);
    *az = atan2(y, x);
    double xy_dist = sqrt(x * x + y * y);
    *el = atan2(z, xy_dist);
}

/* ============================================================================
 * Measurement Jacobians (for EKF linearization)
 * ============================================================================
 */

void measurement_jacobian_polar2d(const double *state, const radar_sensor_t *sensor,
                                   double *H, int state_dim, const int *pos_indices)
{
    if (!state || !sensor || !H || !pos_indices || state_dim <= 0) return;

    double dx = state[pos_indices[0]] - sensor->sensor_x;
    double dy = state[pos_indices[1]] - sensor->sensor_y;
    double d = sqrt(dx * dx + dy * dy);

    if (d < 1e-10) d = 1e-10; /* avoid division by zero */

    /* Zero H */
    memset(H, 0, 2 * state_dim * sizeof(double));

    /* ∂r/∂x = dx/d, ∂r/∂y = dy/d */
    H[0 * state_dim + pos_indices[0]] = dx / d;
    H[0 * state_dim + pos_indices[1]] = dy / d;

    /* ∂θ/∂x = −dy/d², ∂θ/∂y = dx/d² */
    double d2 = d * d;
    H[1 * state_dim + pos_indices[0]] = -dy / d2;
    H[1 * state_dim + pos_indices[1]] = dx / d2;
}

void measurement_jacobian_spherical(const double *state, const radar_sensor_t *sensor,
                                     double *H, int state_dim, const int *pos_indices)
{
    if (!state || !sensor || !H || !pos_indices || state_dim <= 0) return;

    double dx = state[pos_indices[0]] - sensor->sensor_x;
    double dy = state[pos_indices[1]] - sensor->sensor_y;
    double dz = state[pos_indices[2]] - sensor->sensor_z;

    double d_xy = sqrt(dx * dx + dy * dy);
    double d = sqrt(d_xy * d_xy + dz * dz);

    if (d < 1e-10) d = 1e-10;
    if (d_xy < 1e-10) d_xy = 1e-10;

    memset(H, 0, 3 * state_dim * sizeof(double));

    /* ∂r/∂x = dx/d, ∂r/∂y = dy/d, ∂r/∂z = dz/d */
    H[0 * state_dim + pos_indices[0]] = dx / d;
    H[0 * state_dim + pos_indices[1]] = dy / d;
    H[0 * state_dim + pos_indices[2]] = dz / d;

    /* ∂az/∂x = −dy/d_xy², ∂az/∂y = dx/d_xy² */
    double d_xy2 = d_xy * d_xy;
    H[1 * state_dim + pos_indices[0]] = -dy / d_xy2;
    H[1 * state_dim + pos_indices[1]] = dx / d_xy2;

    /* ∂el/∂x = −dx·dz/(d²·d_xy), ∂el/∂y = −dy·dz/(d²·d_xy), ∂el/∂z = d_xy/d² */
    double d2 = d * d;
    H[2 * state_dim + pos_indices[0]] = -dx * dz / (d2 * d_xy);
    H[2 * state_dim + pos_indices[1]] = -dy * dz / (d2 * d_xy);
    H[2 * state_dim + pos_indices[2]] = d_xy / d2;
}

void measurement_jacobian_doppler(const double *state, const radar_sensor_t *sensor,
                                   double *H_dop, int state_dim,
                                   const int *pos_indices, const int *vel_indices)
{
    if (!state || !sensor || !H_dop || !pos_indices || !vel_indices) return;

    double dx = state[pos_indices[0]] - sensor->sensor_x;
    double dy = state[pos_indices[1]] - sensor->sensor_y;
    double dz = (state_dim >= 9 && pos_indices[2] < state_dim)
                ? state[pos_indices[2]] - sensor->sensor_z : 0.0;
    double d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < 1e-10) d = 1e-10;

    double vx = state[vel_indices[0]];
    double vy = state[vel_indices[1]];
    double vz = (state_dim >= 9 && vel_indices[2] < state_dim)
                ? state[vel_indices[2]] : 0.0;

    double vr = (dx * vx + dy * vy + dz * vz) / d;

    memset(H_dop, 0, state_dim * sizeof(double));

    double d2 = d * d;
    /* d3 = d2 * d — used in theoretical formula, optimized out */

    /* ∂vr/∂x = vx/d − dx·(dx·vx+dy·vy+dz·vz)/d³ */
    H_dop[pos_indices[0]] = vx / d - dx * vr / d2;
    H_dop[pos_indices[1]] = vy / d - dy * vr / d2;
    if (state_dim >= 9) {
        H_dop[pos_indices[2]] = vz / d - dz * vr / d2;
    }

    /* ∂vr/∂vx = dx/d */
    H_dop[vel_indices[0]] = dx / d;
    H_dop[vel_indices[1]] = dy / d;
    if (state_dim >= 9) {
        H_dop[vel_indices[2]] = dz / d;
    }
}

/* ============================================================================
 * Measurement prediction
 * ============================================================================
 */

void measurement_predict_polar2d(const double *state, const radar_sensor_t *sensor,
                                  double *z_pred, const int *pos_indices)
{
    if (!state || !sensor || !z_pred || !pos_indices) return;
    double dx = state[pos_indices[0]] - sensor->sensor_x;
    double dy = state[pos_indices[1]] - sensor->sensor_y;
    cartesian_to_polar_2d(dx, dy, &z_pred[0], &z_pred[1]);
}

void measurement_predict_spherical(const double *state, const radar_sensor_t *sensor,
                                    double *z_pred, const int *pos_indices)
{
    if (!state || !sensor || !z_pred || !pos_indices) return;
    double dx = state[pos_indices[0]] - sensor->sensor_x;
    double dy = state[pos_indices[1]] - sensor->sensor_y;
    double dz = state[pos_indices[2]] - sensor->sensor_z;
    cartesian_to_spherical_3d(dx, dy, dz, &z_pred[0], &z_pred[1], &z_pred[2]);
}

double measurement_predict_doppler(const double *state, const radar_sensor_t *sensor,
                                    const int *pos_indices, const int *vel_indices)
{
    if (!state || !sensor || !pos_indices || !vel_indices) return 0.0;
    double dx = state[pos_indices[0]] - sensor->sensor_x;
    double dy = state[pos_indices[1]] - sensor->sensor_y;
    double d = sqrt(dx * dx + dy * dy);
    if (d < 1e-10) return 0.0;
    return (dx * state[vel_indices[0]] + dy * state[vel_indices[1]]) / d;
}

/* ============================================================================
 * Converted Measurement Kalman Filter (CMKF)
 * ============================================================================
 */

void polar_to_cartesian_debiased(double r_meas, double theta_meas,
                                  double sigma_r, double sigma_theta,
                                  double *x_conv, double *y_conv)
{
    (void)sigma_r;
    if (!x_conv || !y_conv) return;

    /* Debiased conversion (Longbin 1998 method).
     *
     * Bias correction factors:
     *   μ_x = cos(θ_m)·(r_m·(1 − exp(−σ²_θ)) − ...) ≈ 0 for small σ_θ
     *   Actually: E[x_conv] = r·cos(θ)·exp(−σ²_θ/2)
     *   So debiased: x_d = r_m·cos(θ_m)·exp(σ²_θ/2)
     */

    double exp_corr = exp(sigma_theta * sigma_theta / 2.0);
    *x_conv = r_meas * cos(theta_meas) * exp_corr;
    *y_conv = r_meas * sin(theta_meas) * exp_corr;
}

void converted_measurement_covariance(double r_meas, double theta_meas,
                                       double sigma_r, double sigma_theta,
                                       double *R_conv)
{
    if (!R_conv) return;

    double cos_t = cos(theta_meas);
    double sin_t = sin(theta_meas);
    double cos2 = cos_t * cos_t;
    double sin2 = sin_t * sin_t;
    double sincos = sin_t * cos_t;
    double r2 = r_meas * r_meas;
    double sr2 = sigma_r * sigma_r;
    double st2 = sigma_theta * sigma_theta;

    /* Longbin 1998 unbiased conversion with consistent covariance:
     *
     * R_xx = r_m²·σ²_θ·sin²(θ_m) + σ²_r·cos²(θ_m)
     *        − (exp(−2σ²_θ)·r²·cos²(θ)·sinh(2σ²_θ) ...)
     *
     * Simplified first-order:
     */
    R_conv[0] = r2 * st2 * sin2 + sr2 * cos2;
    R_conv[1] = (sr2 - r2 * st2) * sincos;
    R_conv[2] = R_conv[1]; /* symmetric */
    R_conv[3] = r2 * st2 * cos2 + sr2 * sin2;

    /* Add small regularization */
    R_conv[0] += 1e-6;
    R_conv[3] += 1e-6;
}

int cmkf_step(double *x, double *P,
               const double *F, const double *Q,
               double r_meas, double theta_meas,
               double sigma_r, double sigma_theta,
               const double *H_cart, int state_dim, int meas_dim)
{
    if (!x || !P || !F || !Q || !H_cart) return -1;
    if (state_dim <= 0 || meas_dim != 2) return -1;

    int n = state_dim;
    int m = meas_dim;

    /* Predict: x_pred = F*x, P_pred = F*P*F' + Q */
    double x_pred[TRACK_MAX_STATE_DIM];
    mat_vec_mul(F, x, x_pred, n, n);

    double FP[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(F, P, FP, n, n, n);
    double FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(F, FT, n, n);
    double P_pred[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(FP, FT, P_pred, n, n, n);
    for (int i = 0; i < n * n; i++) P_pred[i] += Q[i];

    /* Convert measurement to Cartesian */
    double z_cart[2];
    double R_cart[4];
    polar_to_cartesian_debiased(r_meas, theta_meas, sigma_r, sigma_theta,
                                 &z_cart[0], &z_cart[1]);
    converted_measurement_covariance(r_meas, theta_meas, sigma_r, sigma_theta,
                                      R_cart);

    /* Linear KF update with Cartesian measurement */
    double nu[2];
    double Hx_pred[2];
    mat_vec_mul(H_cart, x_pred, Hx_pred, m, n);
    vec_sub(z_cart, Hx_pred, nu, m);

    /* S = H*P_pred*H' + R */
    double PHT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_transpose(H_cart, HT, m, n);
    mat_mat_mul(P_pred, HT, PHT, n, n, m);

    double S[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    mat_mat_mul(H_cart, PHT, S, m, n, m);
    for (int i = 0; i < m * m; i++) S[i] += R_cart[i];

    /* K = P_pred*H' * S^{-1} */
    double S_inv[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    if (mat_inv_cholesky(S, S_inv, m) <= 0.0) return -1;

    double K[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_mat_mul(PHT, S_inv, K, n, m, m);

    /* x = x_pred + K*nu */
    double Knu[TRACK_MAX_STATE_DIM];
    mat_vec_mul(K, nu, Knu, n, m);
    vec_add(x_pred, Knu, x, n);

    /* P = (I - K*H)*P_pred */
    double KH[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(K, H_cart, KH, n, m, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double val = -KH[i * n + j];
            if (i == j) val += 1.0;
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += val * P_pred[k * n + j]; /* I_KH[i][k] * P_pred[k][j] */
            }
            P[i * n + j] = sum;
        }
    }

    /* Recompute more carefully: P = P_pred - K*S*K' (Joseph form) */
    double KSK[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_mat_mul(K, S, KSK, n, m, m);
    double KT[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(K, KT, n, m);
    double KSKT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(KSK, KT, KSKT, n, m, n);

    for (int i = 0; i < n * n; i++) {
        P[i] = P_pred[i] - KSKT[i];
        if (P[i] < 0.0) P[i] = 1e-10; /* ensure PSD */
    }

    return 0;
}

/* ============================================================================
 * Unscented Transform for measurement prediction
 * ============================================================================
 */

void unscented_transform(const double *x_bar, const double *P, int n,
                          void (*h_func)(const double *x, double *z,
                                         int n, int m, void *context),
                          void *h_context, int m,
                          double *z_bar, double *P_zz, double *P_xz,
                          double alpha, double beta, double kappa)
{
    if (!x_bar || !P || !h_func || !z_bar || !P_zz || !P_xz) return;

    double lambda = alpha * alpha * (n + kappa) - n;
    int N_s = 2 * n + 1;

    /* Generate sigma points */
    double *X = (double *)malloc(N_s * n * sizeof(double));
    double *Wm = (double *)malloc(N_s * sizeof(double));
    double *Wc = (double *)malloc(N_s * sizeof(double));
    if (!X || !Wm || !Wc) { free(X); free(Wm); free(Wc); return; }

    /* Cholesky of (n+lambda)*P */
    double P_scaled[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    double scale = n + lambda;
    if (scale <= 0) scale = 0.01;
    for (int i = 0; i < n * n; i++) P_scaled[i] = P[i] * scale;

    double L[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(L, 0, sizeof(L));
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) sum += L[j * n + k] * L[j * n + k];
        double diag = P_scaled[j * n + j] - sum;
        if (diag <= 0.0) diag = 1e-10;
        L[j * n + j] = sqrt(diag);
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) sum += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (P_scaled[i * n + j] - sum) / L[j * n + j];
        }
    }

    /* Sigma point 0 */
    for (int j = 0; j < n; j++) X[0 * n + j] = x_bar[j];
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < n; j++) {
            X[i * n + j] = x_bar[j] + L[j * n + (i - 1)];
            X[(i + n) * n + j] = x_bar[j] - L[j * n + (i - 1)];
        }
    }

    /* Weights */
    Wm[0] = lambda / (n + lambda);
    Wc[0] = Wm[0] + (1.0 - alpha * alpha + beta);
    for (int i = 1; i < N_s; i++) {
        Wm[i] = 0.5 / (n + lambda);
        Wc[i] = Wm[i];
    }

    /* Transform and compute statistics */
    memset(z_bar, 0, m * sizeof(double));
    double *Z = (double *)malloc(N_s * m * sizeof(double));
    if (!Z) { free(X); free(Wm); free(Wc); return; }

    for (int i = 0; i < N_s; i++) {
        h_func(&X[i * n], &Z[i * m], n, m, h_context);
        for (int j = 0; j < m; j++) {
            z_bar[j] += Wm[i] * Z[i * m + j];
        }
    }

    memset(P_zz, 0, m * m * sizeof(double));
    memset(P_xz, 0, n * m * sizeof(double));

    for (int i = 0; i < N_s; i++) {
        double dz[TRACK_MAX_MEAS_DIM];
        vec_sub(&Z[i * m], z_bar, dz, m);

        double dx[TRACK_MAX_STATE_DIM];
        vec_sub(&X[i * n], x_bar, dx, n);

        for (int a = 0; a < m; a++) {
            for (int b = 0; b < m; b++) {
                P_zz[a * m + b] += Wc[i] * dz[a] * dz[b];
            }
        }
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < m; b++) {
                P_xz[a * m + b] += Wc[i] * dx[a] * dz[b];
            }
        }
    }

    free(X); free(Wm); free(Wc); free(Z);
}

/* ============================================================================
 * Radar resolution and ambiguity
 * ============================================================================
 */

int range_ambiguity_fold(double range, double R_unamb)
{
    if (R_unamb <= 0.0) return 0;
    return (int)(range / R_unamb);
}

double doppler_unwrap(double v_measured, double v_unamb)
{
    if (v_unamb <= 0.0) return v_measured;
    /* Fold velocity into [-v_unamb/2, v_unamb/2] */
    double v = v_measured;
    while (v > v_unamb / 2.0)  v -= v_unamb;
    while (v < -v_unamb / 2.0) v += v_unamb;
    return v;
}

double crlb_range(double bandwidth, double SNR)
{
    if (bandwidth <= 0.0 || SNR <= 0.0) return 1e10;
    double c = 299792458.0; /* speed of light */
    /* σ²_r ≥ c² / (8·π²·B²·SNR) */
    return (c * c) / (8.0 * M_PI * M_PI * bandwidth * bandwidth * SNR);
}

double crlb_bearing(double wavelength, double aperture, double SNR, double theta)
{
    if (wavelength <= 0.0 || aperture <= 0.0 || SNR <= 0.0) return 1e10;
    double cos_theta = cos(theta);
    if (fabs(cos_theta) < 1e-10) return 1e10;
    /* σ²_θ ≥ λ² / (8·π²·D²·SNR·cos²(θ)) */
    return (wavelength * wavelength)
           / (8.0 * M_PI * M_PI * aperture * aperture * SNR * cos_theta * cos_theta);
}

double crlb_doppler(double wavelength, double T_cpi, double SNR)
{
    if (wavelength <= 0.0 || T_cpi <= 0.0 || SNR <= 0.0) return 1e10;
    /* σ²_v ≥ λ² / (8·π²·T_cpi²·SNR) */
    return (wavelength * wavelength)
           / (8.0 * M_PI * M_PI * T_cpi * T_cpi * SNR);
}
