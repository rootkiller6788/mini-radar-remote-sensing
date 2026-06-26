/**
 * kalman_filter.c — Kalman Filter Variants Implementation
 *
 * Implements: KF predict/update, EKF, UKF, SR-KF, Information Filter,
 *             Adaptive KF, IMM mixing/combination.
 *
 * References:
 *   - Kalman, R.E. (1960)
 *   - Bar-Shalom, Willett, Tian (2011)
 *   - Julier & Uhlmann (2004)
 *   - Bierman, G.J. "Factorization Methods for Discrete Sequential Estimation" (1977)
 */

#include "kalman_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Filter allocation and setup
 * ============================================================================
 */

kf_t *kf_alloc(int state_dim, int meas_dim, int control_dim)
{
    if (state_dim <= 0 || meas_dim <= 0) return NULL;
    if (state_dim > TRACK_MAX_STATE_DIM || meas_dim > TRACK_MAX_MEAS_DIM)
        return NULL;

    kf_t *kf = (kf_t *)calloc(1, sizeof(kf_t));
    if (!kf) return NULL;

    kf->type = KF_LINEAR;
    kf->state_dim = state_dim;
    kf->meas_dim = meas_dim;
    kf->control_dim = control_dim;

    int n = state_dim;
    int m = meas_dim;
    int c = control_dim;

    kf->x     = (double *)calloc(n, sizeof(double));
    kf->P     = (double *)calloc(n * n, sizeof(double));
    kf->x_pred = (double *)calloc(n, sizeof(double));
    kf->P_pred = (double *)calloc(n * n, sizeof(double));
    kf->nu    = (double *)calloc(m, sizeof(double));
    kf->S     = (double *)calloc(m * m, sizeof(double));
    kf->K     = (double *)calloc(n * m, sizeof(double));
    kf->F     = (double *)calloc(n * n, sizeof(double));
    kf->H     = (double *)calloc(m * n, sizeof(double));
    kf->Q     = (double *)calloc(n * n, sizeof(double));
    kf->R     = (double *)calloc(m * m, sizeof(double));

    kf->work_n = (double *)calloc(n * n, sizeof(double));
    kf->work_m = (double *)calloc(m * m, sizeof(double));

    if (c > 0) {
        kf->B = (double *)calloc(n * c, sizeof(double));
        kf->u = (double *)calloc(c, sizeof(double));
    }

    /* Initialize identity transition */
    for (int i = 0; i < n; i++) {
        kf->F[i * n + i] = 1.0;
    }
    /* Initialize identity measurement (first m states measured) */
    for (int i = 0; i < m && i < n; i++) {
        kf->H[i * n + i] = 1.0;
    }

    kf->log_likelihood = 0.0;
    kf->steps = 0;
    kf->diverged = 0;
    kf->forgetting_factor = 0.98;

    kf->num_sigma = 2 * n + 1;
    kf->alpha_ukf = 0.001; /* default: tight sigma spread */
    kf->beta_ukf = 2.0;
    kf->kappa_ukf = 0.0;
    kf->lambda_ukf = kf->alpha_ukf * kf->alpha_ukf * (n + kf->kappa_ukf) - n;

    return kf;
}

void kf_free(kf_t *kf)
{
    if (!kf) return;
    free(kf->x);
    free(kf->P);
    free(kf->x_pred);
    free(kf->P_pred);
    free(kf->nu);
    free(kf->S);
    free(kf->K);
    free(kf->F);
    free(kf->H);
    free(kf->Q);
    free(kf->R);
    free(kf->work_n);
    free(kf->work_m);
    free(kf->B);
    free(kf->u);
    free(kf->Q_est);
    free(kf->R_est);
    free(kf->sigma_points);
    free(kf->sigma_weights);
    free(kf->sigma_pred);
    free(kf->sigma_meas);
    free(kf->S_chol);
    free(kf->model_probs);
    free(kf->mixing_matrix);
    free(kf);
}

void kf_set_model(kf_t *kf,
                   const double *F, const double *H,
                   const double *Q, const double *R,
                   const double *B, const double *u)
{
    if (!kf) return;
    int n = kf->state_dim;
    int m = kf->meas_dim;
    int c = kf->control_dim;

    if (F) memcpy(kf->F, F, n * n * sizeof(double));
    if (H) memcpy(kf->H, H, m * n * sizeof(double));
    if (Q) memcpy(kf->Q, Q, n * n * sizeof(double));
    if (R) memcpy(kf->R, R, m * m * sizeof(double));
    if (B && c > 0) memcpy(kf->B, B, n * c * sizeof(double));
    if (u && c > 0) memcpy(kf->u, u, c * sizeof(double));
}

void kf_set_state(kf_t *kf, const double *x0, const double *P0)
{
    if (!kf) return;
    int n = kf->state_dim;
    if (x0) memcpy(kf->x, x0, n * sizeof(double));
    if (P0) memcpy(kf->P, P0, n * n * sizeof(double));
}

/* ============================================================================
 * Linear KF predict
 * ============================================================================
 */

int kf_predict(kf_t *kf)
{
    if (!kf) return -1;
    int n = kf->state_dim;
    int c = kf->control_dim;

    double *F = kf->F;
    double *x = kf->x;
    double *P = kf->P;
    double *Q = kf->Q;
    double *x_pred = kf->x_pred;
    double *P_pred = kf->P_pred;

    /* x_pred = F*x (+ B*u) */
    mat_vec_mul(F, x, x_pred, n, n);
    if (c > 0 && kf->B && kf->u) {
        double Bu[TRACK_MAX_STATE_DIM];
        mat_vec_mul(kf->B, kf->u, Bu, n, c);
        for (int i = 0; i < n; i++) x_pred[i] += Bu[i];
    }

    /* P_pred = F*P*F' + Q */
    double FP[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(F, P, FP, n, n, n);

    double FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(F, FT, n, n);

    mat_mat_mul(FP, FT, P_pred, n, n, n);

    for (int i = 0; i < n * n; i++) {
        P_pred[i] += Q[i];
    }

    kf->steps++;
    return 0;
}

int kf_update(kf_t *kf, const double *z)
{
    if (!kf || !z) return -1;
    int n = kf->state_dim;
    int m = kf->meas_dim;

    double *H = kf->H;
    double *R = kf->R;
    double *x_pred = kf->x_pred;
    double *P_pred = kf->P_pred;
    double *nu = kf->nu;
    double *S  = kf->S;
    double *K  = kf->K;

    /* Innovation: nu = z - H*x_pred */
    double Hx[TRACK_MAX_MEAS_DIM];
    mat_vec_mul(H, x_pred, Hx, m, n);
    vec_sub(z, Hx, nu, m);

    /* S = H*P_pred*H' + R */
    double PHT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_transpose(H, HT, m, n);
    mat_mat_mul(P_pred, HT, PHT, n, n, m);
    mat_mat_mul(H, PHT, S, m, n, m);
    for (int i = 0; i < m * m; i++) S[i] += R[i];

    /* K = P_pred*H' * S^{-1} */
    double S_inv[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    if (mat_inv_cholesky(S, S_inv, m) <= 0.0) return -1;
    mat_mat_mul(PHT, S_inv, K, n, m, m);

    /* x = x_pred + K*nu */
    double Knu[TRACK_MAX_STATE_DIM];
    mat_vec_mul(K, nu, Knu, n, m);
    vec_add(x_pred, Knu, kf->x, n);

    /* P = (I - K*H)*P_pred  (Joseph form for symmetry) */
    /* First: P = P_pred */
    memcpy(kf->P, P_pred, n * n * sizeof(double));

    /* Then the standard form: P = (I-KH)*P_pred */
    double KH[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(K, H, KH, n, m, n);

    double I_KH[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            I_KH[i * n + j] = -KH[i * n + j];
            if (i == j) I_KH[i * n + j] += 1.0;
        }
    }

    double P_new[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(I_KH, P_pred, P_new, n, n, n);
    /* Ensure symmetry */
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double avg = 0.5 * (P_new[i * n + j] + P_new[j * n + i]);
            kf->P[i * n + j] = avg;
            kf->P[j * n + i] = avg;
        }
    }

    /* Log-likelihood: -(1/2) * [m*log(2π) + log|S| + nu'*S^{-1}*nu] */
    double det_S = mat_det_cholesky(S, m);
    double S_inv_nu[TRACK_MAX_MEAS_DIM];
    memset(S_inv_nu, 0, sizeof(S_inv_nu));
    /* Use solve: S * t = nu, then nu' * t */
    mat_solve_cholesky(S, nu, S_inv_nu, m);
    double nll = 0.5 * (m * log(2.0 * M_PI) + log(det_S > 0 ? det_S : 1e-300)
                        + vec_dot(nu, S_inv_nu, m));
    kf->log_likelihood += nll;

    return 0;
}

int kf_step(kf_t *kf, const double *z)
{
    if (kf_predict(kf) != 0) return -1;
    return kf_update(kf, z);
}

/* ============================================================================
 * Extended Kalman Filter
 * ============================================================================
 */

int ekf_predict(kf_t *kf)
{
    if (!kf) return -1;
    int n = kf->state_dim;
    int c = kf->control_dim;

    /* Nonlinear state transition */
    if (kf->f_nonlinear) {
        kf->f_nonlinear(kf->x, kf->u, kf->x_pred, n, c);
    } else {
        /* Fall back to linear F*x */
        mat_vec_mul(kf->F, kf->x, kf->x_pred, n, n);
        if (c > 0 && kf->B && kf->u) {
            double Bu[TRACK_MAX_STATE_DIM];
            mat_vec_mul(kf->B, kf->u, Bu, n, c);
            for (int i = 0; i < n; i++) kf->x_pred[i] += Bu[i];
        }
    }

    /* Jacobian F = ∂f/∂x */
    if (kf->compute_F_jacobian) {
        kf->compute_F_jacobian(kf->x, kf->u, kf->F, n, c);
    }

    /* P_pred = F*P*F' + Q */
    double FP[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(kf->F, kf->P, FP, n, n, n);

    double FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(kf->F, FT, n, n);

    mat_mat_mul(FP, FT, kf->P_pred, n, n, n);
    for (int i = 0; i < n * n; i++) kf->P_pred[i] += kf->Q[i];

    kf->steps++;
    return 0;
}

int ekf_update(kf_t *kf, const double *z)
{
    if (!kf || !z) return -1;

    /* Jacobian H = ∂h/∂x at predicted state */
    if (kf->compute_H_jacobian) {
        kf->compute_H_jacobian(kf->x_pred, kf->H, kf->state_dim, kf->meas_dim);
    }

    /* Standard KF update with linearized H */
    return kf_update(kf, z);
}

/* ============================================================================
 * Unscented Kalman Filter
 * ============================================================================
 */

void ukf_init_params(kf_t *kf, double alpha, double beta, double kappa)
{
    if (!kf) return;
    kf->alpha_ukf = alpha;
    kf->beta_ukf = beta;
    kf->kappa_ukf = kappa;
    int n = kf->state_dim;
    kf->lambda_ukf = alpha * alpha * (n + kappa) - n;
}

int ukf_generate_sigma_points(kf_t *kf)
{
    if (!kf) return -1;
    int n = kf->state_dim;
    int N_s = kf->num_sigma; /* 2n+1 */
    double lambda = kf->lambda_ukf;
    double alpha = kf->alpha_ukf;

    /* Allocate sigma point storage */
    double *X = (double *)realloc(kf->sigma_points, N_s * n * sizeof(double));
    double *Wm = (double *)realloc(kf->sigma_weights, N_s * sizeof(double));
    double *Wc = (double *)realloc(kf->sigma_pred, N_s * sizeof(double));
    if (!X || !Wm || !Wc) return -1;
    kf->sigma_points = X;
    kf->sigma_weights = Wm;
    /* Wc stored in sigma_pred temporarily */

    /* Compute sqrt of (n+lambda)*P via Cholesky */
    double sqrtP[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(sqrtP, 0, sizeof(sqrtP));

    /* Scale P by (n+lambda), then Cholesky factor */
    double scale = n + lambda;
    if (scale < 0) scale = 0.1; /* ensure positive */

    double P_scaled[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    for (int i = 0; i < n * n; i++) P_scaled[i] = kf->P[i] * scale;

    /* Cholesky: P_scaled = L*L' */
    double L[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(L, 0, sizeof(L));

    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) sum += L[j * n + k] * L[j * n + k];
        double diag = P_scaled[j * n + j] - sum;
        if (diag <= 0.0) { /* add small regularization */
            diag = 1e-6;
        }
        L[j * n + j] = sqrt(diag);
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) sum += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (P_scaled[i * n + j] - sum) / L[j * n + j];
        }
    }

    /* Sigma point 0: mean */
    for (int j = 0; j < n; j++) X[0 * n + j] = kf->x[j];

    /* Sigma points 1..n: mean + L column */
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < n; j++) {
            X[i * n + j] = kf->x[j] + L[j * n + (i - 1)];
        }
    }

    /* Sigma points n+1..2n: mean - L column */
    for (int i = n + 1; i <= 2 * n; i++) {
        for (int j = 0; j < n; j++) {
            X[i * n + j] = kf->x[j] - L[j * n + (i - n - 1)];
        }
    }

    /* Weights for mean */
    Wm[0] = lambda / (n + lambda);
    for (int i = 1; i < N_s; i++) Wm[i] = 0.5 / (n + lambda);

    /* Weights for covariance */
    Wc[0] = lambda / (n + lambda) + (1.0 - alpha * alpha + kf->beta_ukf);
    for (int i = 1; i < N_s; i++) Wc[i] = 0.5 / (n + lambda);

    return 0;
}

int ukf_predict(kf_t *kf)
{
    if (!kf) return -1;
    int n = kf->state_dim;
    int N_s = kf->num_sigma;
    int c = kf->control_dim;

    if (!kf->sigma_points) {
        if (ukf_generate_sigma_points(kf) != 0) return -1;
    }

    /* Transform sigma points through f() */
    if (!kf->sigma_pred) {
        kf->sigma_pred = (double *)calloc(N_s * n, sizeof(double));
        if (!kf->sigma_pred) return -1;
    }

    for (int i = 0; i < N_s; i++) {
        if (kf->f_nonlinear) {
            kf->f_nonlinear(&kf->sigma_points[i * n], kf->u,
                             &kf->sigma_pred[i * n], n, c);
        } else {
            /* Linear: x_next = F*x */
            mat_vec_mul(kf->F, &kf->sigma_points[i * n],
                        &kf->sigma_pred[i * n], n, n);
        }
    }

    /* Predicted mean: x_pred = Σ Wm_i * χ_pred_i */
    double *Wm = kf->sigma_weights;
    memset(kf->x_pred, 0, n * sizeof(double));
    for (int i = 0; i < N_s; i++) {
        for (int j = 0; j < n; j++) {
            kf->x_pred[j] += Wm[i] * kf->sigma_pred[i * n + j];
        }
    }

    /* Predicted covariance: P_pred = Σ Wc_i*(χ_pred_i - x_pred)*(χ_pred_i - x_pred)' + Q */
    double *Wc = kf->sigma_pred; /* Wc was stored here in ukf_generate_sigma_points */
    memset(kf->P_pred, 0, n * n * sizeof(double));
    for (int i = 0; i < N_s; i++) {
        double dx[TRACK_MAX_STATE_DIM];
        vec_sub(&kf->sigma_pred[i * n], kf->x_pred, dx, n);
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                kf->P_pred[a * n + b] += Wc[i] * dx[a] * dx[b];
            }
        }
    }
    for (int i = 0; i < n * n; i++) kf->P_pred[i] += kf->Q[i];

    kf->steps++;
    return 0;
}

int ukf_update(kf_t *kf, const double *z)
{
    if (!kf || !z) return -1;
    int n = kf->state_dim;
    int m = kf->meas_dim;
    int N_s = kf->num_sigma;

    if (!kf->sigma_meas) {
        kf->sigma_meas = (double *)calloc(N_s * m, sizeof(double));
        if (!kf->sigma_meas) return -1;
    }

    /* Transform predicted sigma points through h() */
    for (int i = 0; i < N_s; i++) {
        if (kf->h_nonlinear) {
            kf->h_nonlinear(&kf->sigma_pred[i * n],
                             &kf->sigma_meas[i * m], n, m);
        } else {
            mat_vec_mul(kf->H, &kf->sigma_pred[i * n],
                        &kf->sigma_meas[i * m], m, n);
        }
    }

    /* Predicted measurement mean */
    double *Wm = kf->sigma_weights;
    double z_pred[TRACK_MAX_MEAS_DIM];
    memset(z_pred, 0, m * sizeof(double));
    for (int i = 0; i < N_s; i++) {
        for (int j = 0; j < m; j++) {
            z_pred[j] += Wm[i] * kf->sigma_meas[i * m + j];
        }
    }

    /* Innovation covariance S */
    double *Wc = kf->sigma_pred; /* Wc stored here */
    memset(kf->S, 0, m * m * sizeof(double));
    for (int i = 0; i < N_s; i++) {
        double dz[TRACK_MAX_MEAS_DIM];
        vec_sub(&kf->sigma_meas[i * m], z_pred, dz, m);
        for (int a = 0; a < m; a++) {
            for (int b = 0; b < m; b++) {
                kf->S[a * m + b] += Wc[i] * dz[a] * dz[b];
            }
        }
    }
    for (int i = 0; i < m * m; i++) kf->S[i] += kf->R[i];

    /* Cross-covariance Pxz */
    double Pxz[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    memset(Pxz, 0, sizeof(Pxz));
    for (int i = 0; i < N_s; i++) {
        double dx[TRACK_MAX_STATE_DIM];
        double dz2[TRACK_MAX_MEAS_DIM];
        vec_sub(&kf->sigma_pred[i * n], kf->x_pred, dx, n);
        vec_sub(&kf->sigma_meas[i * m], z_pred, dz2, m);
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < m; b++) {
                Pxz[a * m + b] += Wc[i] * dx[a] * dz2[b];
            }
        }
    }

    /* Kalman gain K = Pxz * S^{-1} */
    double S_inv[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    if (mat_inv_cholesky(kf->S, S_inv, m) <= 0.0) return -1;
    mat_mat_mul(Pxz, S_inv, kf->K, n, m, m);

    /* Update state */
    double nu[TRACK_MAX_MEAS_DIM];
    vec_sub(z, z_pred, nu, m);
    double Knu[TRACK_MAX_STATE_DIM];
    mat_vec_mul(kf->K, nu, Knu, n, m);
    vec_add(kf->x_pred, Knu, kf->x, n);

    /* Update covariance: P = P_pred - K*S*K' */
    double KSKT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    double KSK[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_mat_mul(kf->K, kf->S, KSK, n, m, m);
    double KT[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(kf->K, KT, n, m);
    mat_mat_mul(KSK, KT, KSKT, n, m, n);
    for (int i = 0; i < n * n; i++) {
        kf->P[i] = kf->P_pred[i] - KSKT[i];
    }

    return 0;
}

/* ============================================================================
 * Square-Root Kalman (Potter scalar update)
 * ============================================================================
 */

int srkf_scalar_update(kf_t *kf, const double *H_row,
                        double z_scalar, double R_scalar)
{
    if (!kf || !H_row) return -1;
    int n = kf->state_dim;

    if (!kf->S_chol) {
        /* Initialize Cholesky factor from P if not set */
        kf->S_chol = (double *)calloc(n * n, sizeof(double));
        if (!kf->S_chol) return -1;
        /* Compute Cholesky of P */
        double P_copy[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        memcpy(P_copy, kf->P, n * n * sizeof(double));

        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += kf->S_chol[j * n + k] * kf->S_chol[j * n + k];
            double diag = P_copy[j * n + j] - sum;
            if (diag <= 0.0) diag = 1e-10;
            kf->S_chol[j * n + j] = sqrt(diag);
            for (int i = j + 1; i < n; i++) {
                sum = 0.0;
                for (int k = 0; k < j; k++)
                    sum += kf->S_chol[i * n + k] * kf->S_chol[j * n + k];
                kf->S_chol[i * n + j] = (P_copy[i * n + j] - sum)
                                        / kf->S_chol[j * n + j];
            }
        }
    }

    /* Potter's algorithm for scalar measurement update.
     *
     * f = S' * H'
     * α = f'*f + R
     * γ = 1 / (α + √(α*R))
     * K = γ * S * f
     * x = x_pred + K * (z - H*x_pred)
     * S = S - γ * (S*f) * f' / (1 + √(R/α))
     */
    double f[TRACK_MAX_STATE_DIM];
    memset(f, 0, n * sizeof(double));

    /* f = S' * H' */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            f[i] += kf->S_chol[j * n + i] * H_row[j];
        }
    }

    double alpha = vec_dot(f, f, n) + R_scalar;
    double sqrt_R_alpha = sqrt(R_scalar * alpha);
    double gamma = 1.0 / (alpha + sqrt_R_alpha);

    /* Innovation */
    double Hx_pred = vec_dot(H_row, kf->x_pred, n);
    double nu = z_scalar - Hx_pred;

    /* State update */
    double Sf[TRACK_MAX_STATE_DIM];
    mat_vec_mul(kf->S_chol, f, Sf, n, n); /* S*f — note S is lower triangular */
    for (int i = 0; i < n; i++) {
        kf->x[i] = kf->x_pred[i] + gamma * Sf[i] * nu;
    }

    /* Cholesky factor update: S_new = S - gamma' * (S*f) * f' */
    double gamma_prime = 1.0 / (1.0 + sqrt(R_scalar / alpha));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double correction = gamma_prime * Sf[i] * f[j];
            kf->S_chol[i * n + j] -= correction;
            if (i == j && kf->S_chol[i * n + j] < 0.0) {
                kf->S_chol[i * n + j] = 1e-10;
            }
        }
    }

    /* Propagate to P */
    double ST[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(kf->S_chol, ST, n, n);
    mat_mat_mul(kf->S_chol, ST, kf->P, n, n, n);

    return 0;
}

/* ============================================================================
 * Information Filter
 * ============================================================================
 */

int info_filter_update(double *Y_info, double *y_info,
                        const double *H, const double *R_inv,
                        const double *z, int n, int m)
{
    if (!Y_info || !y_info || !H || !R_inv || !z) return -1;

    /* Update: Y += H' * R^{-1} * H */
    /*          y += H' * R^{-1} * z */

    /* Compute H' * R_inv */
    double HT_Rinv[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < m; k++) {
                sum += H[k * n + i] * R_inv[k * m + j]; /* H'[i][k] * R_inv[k][j] */
            }
            HT_Rinv[i * m + j] = sum;
        }
    }

    /* Y += HT_Rinv * H */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < m; k++) {
                sum += HT_Rinv[i * m + k] * H[k * n + j];
            }
            Y_info[i * n + j] += sum;
        }
    }

    /* y += HT_Rinv * z */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < m; k++) {
            sum += HT_Rinv[i * m + k] * z[k];
        }
        y_info[i] += sum;
    }

    return 0;
}

/* ============================================================================
 * Adaptive Kalman filtering
 * ============================================================================
 */

void kf_adapt_R(kf_t *kf)
{
    if (!kf) return;
    int n = kf->state_dim;
    int m = kf->meas_dim;
    double alpha = 1.0 - kf->forgetting_factor;

    if (!kf->R_est) {
        kf->R_est = (double *)calloc(m * m, sizeof(double));
        if (!kf->R_est) return;
        memcpy(kf->R_est, kf->R, m * m * sizeof(double));
    }

    /* E[nu*nu'] = S = H*P_pred*H' + R
     * R_est = (1-alpha)*R_est + alpha*(nu*nu' - H*P_pred*H')
     */
    double HPH[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    double PHT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_transpose(kf->H, HT, m, n);
    mat_mat_mul(kf->P_pred, HT, PHT, n, n, m);
    mat_mat_mul(kf->H, PHT, HPH, m, n, m);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double nu_nu = kf->nu[i] * kf->nu[j];
            double raw_R = nu_nu - HPH[i * m + j];
            if (raw_R < 0.0) raw_R = 0.0; /* clip negative */
            kf->R_est[i * m + j] = (1.0 - alpha) * kf->R_est[i * m + j]
                                   + alpha * raw_R;
        }
    }
}

void kf_adapt_Q(kf_t *kf)
{
    if (!kf) return;
    int n = kf->state_dim;
    int m = kf->meas_dim;
    double alpha = 1.0 - kf->forgetting_factor;

    if (!kf->Q_est) {
        kf->Q_est = (double *)calloc(n * n, sizeof(double));
        if (!kf->Q_est) return;
        memcpy(kf->Q_est, kf->Q, n * n * sizeof(double));
    }

    /* Q_est = (1-alpha)*Q_est + alpha*(K*nu*nu'*K') */
    double Knu[TRACK_MAX_STATE_DIM];
    mat_vec_mul(kf->K, kf->nu, Knu, n, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dQ = Knu[i] * Knu[j];
            kf->Q_est[i * n + j] = (1.0 - alpha) * kf->Q_est[i * n + j]
                                   + alpha * dQ;
        }
    }
}

/* ============================================================================
 * IMM
 * ============================================================================
 */

int imm_mix(kf_t *kf)
{
    if (!kf || !kf->sub_filters || !kf->mixing_matrix || !kf->model_probs)
        return -1;
    int r = kf->num_models;
    if (r <= 0) return -1;

    kf_t **filters = (kf_t **)kf->sub_filters;
    int n = filters[0]->state_dim;

    /* Mixing probabilities: mu_{i|j} = p_{ij} * mu_i / cbar_j */
    double *mu = kf->model_probs;
    double *pij = kf->mixing_matrix;

    double cbar[16]; /* r ≤ 16 assumed */
    for (int j = 0; j < r; j++) {
        cbar[j] = 0.0;
        for (int i = 0; i < r; i++) {
            cbar[j] += pij[i * r + j] * mu[i];
        }
    }

    /* For each sub-filter j, compute mixed initial state and covariance */
    for (int j = 0; j < r; j++) {
        double *x_mixed = kf->work_n; /* temporary buffer */
        memset(x_mixed, 0, n * sizeof(double));

        for (int i = 0; i < r; i++) {
            double mu_ij = pij[i * r + j] * mu[i] / cbar[j];
            for (int d = 0; d < n; d++) {
                x_mixed[d] += mu_ij * filters[i]->x[d];
            }
        }

        double P_mixed[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        memset(P_mixed, 0, sizeof(P_mixed));

        for (int i = 0; i < r; i++) {
            double mu_ij = pij[i * r + j] * mu[i] / cbar[j];
            double dx[TRACK_MAX_STATE_DIM];
            vec_sub(filters[i]->x, x_mixed, dx, n);

            for (int a = 0; a < n; a++) {
                for (int b = 0; b < n; b++) {
                    P_mixed[a * n + b] += mu_ij * (
                        filters[i]->P[a * n + b] + dx[a] * dx[b]);
                }
            }
        }

        /* Copy mixed state to sub-filter j */
        memcpy(filters[j]->x, x_mixed, n * sizeof(double));
        memcpy(filters[j]->P, P_mixed, n * n * sizeof(double));
    }

    return 0;
}

int imm_update_probs(kf_t *kf)
{
    if (!kf || !kf->sub_filters || !kf->mixing_matrix || !kf->model_probs)
        return -1;
    int r = kf->num_models;
    int m = kf->meas_dim;

    kf_t **filters = (kf_t **)kf->sub_filters;
    double *mu = kf->model_probs;
    double *pij = kf->mixing_matrix;

    /* Compute likelihood for each model:
     * Lambda_j = N(nu_j; 0, S_j) = exp(-0.5*nu_j'*S_j^{-1}*nu_j) / sqrt((2π)^m * |S_j|)
     */
    double Lambda[16];
    for (int j = 0; j < r; j++) {
        double det_S = mat_det_cholesky(filters[j]->S, m);
        if (det_S <= 0.0) det_S = 1e-300;
        double S_inv_nu[TRACK_MAX_MEAS_DIM];
        double nu_copy[TRACK_MAX_MEAS_DIM];
        vec_copy(filters[j]->nu, nu_copy, m);
        mat_solve_cholesky(filters[j]->S, nu_copy, S_inv_nu, m);
        double nis = vec_dot(filters[j]->nu, S_inv_nu, m);
        Lambda[j] = exp(-0.5 * nis) / sqrt(pow(2.0 * M_PI, m) * det_S);
        if (Lambda[j] < 1e-300) Lambda[j] = 1e-300;
    }

    /* cbar_j first */
    double cbar[16];
    for (int j = 0; j < r; j++) {
        cbar[j] = 0.0;
        for (int i = 0; i < r; i++) {
            cbar[j] += pij[i * r + j] * mu[i];
        }
    }

    /* Updated probabilities */
    double mu_new[16];
    double c_total = 0.0;
    for (int j = 0; j < r; j++) {
        mu_new[j] = Lambda[j] * cbar[j];
        c_total += mu_new[j];
    }
    if (c_total > 0.0) {
        for (int j = 0; j < r; j++) {
            mu[j] = mu_new[j] / c_total;
        }
    }

    return 0;
}

int imm_combine(kf_t *kf)
{
    if (!kf || !kf->sub_filters || !kf->model_probs) return -1;
    int r = kf->num_models;
    int n = kf->state_dim;

    kf_t **filters = (kf_t **)kf->sub_filters;
    double *mu = kf->model_probs;

    /* Combined state: x = Σ μ_j * x_j */
    memset(kf->x, 0, n * sizeof(double));
    for (int j = 0; j < r; j++) {
        for (int d = 0; d < n; d++) {
            kf->x[d] += mu[j] * filters[j]->x[d];
        }
    }

    /* Combined covariance: P = Σ μ_j * [P_j + (x_j - x)(x_j - x)'] */
    memset(kf->P, 0, n * n * sizeof(double));
    for (int j = 0; j < r; j++) {
        double dx[TRACK_MAX_STATE_DIM];
        vec_sub(filters[j]->x, kf->x, dx, n);
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                kf->P[a * n + b] += mu[j] * (
                    filters[j]->P[a * n + b] + dx[a] * dx[b]);
            }
        }
    }

    return 0;
}

/* ============================================================================
 * Filter diagnostics
 * ============================================================================
 */

double kf_nis(const kf_t *kf)
{
    if (!kf) return 0.0;
    int m = kf->meas_dim;

    double S_inv[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    double nu_copy[TRACK_MAX_MEAS_DIM];
    vec_copy(kf->nu, nu_copy, m);

    if (mat_inv_cholesky(kf->S, S_inv, m) <= 0.0) return 0.0;

    double nis = 0.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            nis += nu_copy[i] * S_inv[i * m + j] * nu_copy[j];
        }
    }
    return nis;
}

double kf_nees(const kf_t *kf, const double *x_true)
{
    if (!kf || !x_true) return 0.0;
    int n = kf->state_dim;

    double dx[TRACK_MAX_STATE_DIM];
    vec_sub(kf->x, x_true, dx, n);

    double P_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    if (mat_inv_cholesky(kf->P, P_inv, n) <= 0.0) return 0.0;

    double nees = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            nees += dx[i] * P_inv[i * n + j] * dx[j];
        }
    }
    return nees;
}

int kf_consistency_test(double epsilon, int dof, double alpha)
{
    /* Compute chi-squared threshold and compare */
    double threshold = chi2_threshold_approx(1.0 - alpha, dof);
    return (epsilon <= threshold) ? 1 : 0;
}

int kf_divergence_detect(kf_t *kf, double threshold, int consecutive)
{
    if (!kf) return 0;
    double nis = kf_nis(kf);
    if (nis > threshold) {
        static int count = 0;
        count++;
        if (count >= consecutive) {
            count = 0;
            kf->diverged = 1;
            return 1;
        }
        return 0;
    }
    /* Reset counter on in-bound NIS */
    return 0;
}

double kf_neg_log_likelihood(const kf_t *kf)
{
    if (!kf) return 0.0;
    return kf->log_likelihood;
}
