/**
 * motion_models.c — Target Motion Models Implementation
 *
 * Implements: CV, CA, CT, Singer, WPA, Brownian motion models
 *             with discretized state transition and process noise matrices.
 *
 * References:
 *   - Li, X.R., Jilkov, V.P. "Survey of Maneuvering Target Tracking" (2003)
 *   - Singer, R.A. "Estimating Optimal Tracking Filter Performance" (1970)
 *   - Bar-Shalom, Willett, Tian (2011), Ch. 6
 */

#include "motion_models.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Constant Velocity (CV) model
 * ============================================================================
 */

motion_model_t *motion_create_cv(motion_dim_t dim, double T, double q)
{
    if (T <= 0.0 || q < 0.0) return NULL;

    /* state_dim = dim * 2 (position + velocity per axis) */
    int state_dim = (int)dim * 2;
    int n = state_dim;

    motion_model_t *model = (motion_model_t *)calloc(1, sizeof(motion_model_t));
    if (!model) return NULL;
    model->type = MOTION_CV;
    model->dim = dim;
    model->state_dim = n;
    model->T = T;
    model->q = q;

    model->F = (double *)calloc(n * n, sizeof(double));
    model->Q = (double *)calloc(n * n, sizeof(double));
    if (!model->F || !model->Q) { motion_free(model); return NULL; }

    /* Build F: per-axis block [1 T; 0 1] */
    for (int d = 0; d < (int)dim; d++) {
        int i = d * 2; /* position index */
        model->F[i * n + i]     = 1.0;
        model->F[i * n + i + 1] = T;
        model->F[(i + 1) * n + i + 1] = 1.0;
    }

    /* Build Q = q * block([T³/3, T²/2; T²/2, T]) per axis.
     * q is the continuous-time white acceleration PSD.
     */
    double T3_3 = T * T * T / 3.0;
    double T2_2 = T * T / 2.0;

    for (int d = 0; d < (int)dim; d++) {
        int i = d * 2;
        model->Q[i * n + i]         = q * T3_3;
        model->Q[i * n + i + 1]     = q * T2_2;
        model->Q[(i + 1) * n + i]   = q * T2_2;
        model->Q[(i + 1) * n + i + 1] = q * T;
    }

    return model;
}

/* ============================================================================
 * Constant Acceleration (CA) model
 * ============================================================================
 */

motion_model_t *motion_create_ca(motion_dim_t dim, double T, double q)
{
    if (T <= 0.0 || q < 0.0) return NULL;

    int n = (int)dim * 3; /* [pos, vel, acc] per axis */

    motion_model_t *model = (motion_model_t *)calloc(1, sizeof(motion_model_t));
    if (!model) return NULL;
    model->type = MOTION_CA;
    model->dim = dim;
    model->state_dim = n;
    model->T = T;
    model->q = q;

    model->F = (double *)calloc(n * n, sizeof(double));
    model->Q = (double *)calloc(n * n, sizeof(double));
    if (!model->F || !model->Q) { motion_free(model); return NULL; }

    /* F = [1 T T²/2; 0 1 T; 0 0 1] per axis */
    for (int d = 0; d < (int)dim; d++) {
        int i = d * 3;
        model->F[i * n + i]     = 1.0;
        model->F[i * n + i + 1] = T;
        model->F[i * n + i + 2] = T * T / 2.0;
        model->F[(i + 1) * n + i + 1] = 1.0;
        model->F[(i + 1) * n + i + 2] = T;
        model->F[(i + 2) * n + i + 2] = 1.0;
    }

    /* Q = q * [T⁵/20, T⁴/8, T³/6;
     *          T⁴/8,  T³/3, T²/2;
     *          T³/6,  T²/2, T   ] per axis
     */
    double q00 = q * T * T * T * T * T / 20.0;
    double q01 = q * T * T * T * T / 8.0;
    double q02 = q * T * T * T / 6.0;
    double q11 = q * T * T * T / 3.0;
    double q12 = q * T * T / 2.0;
    double q22 = q * T;

    for (int d = 0; d < (int)dim; d++) {
        int i = d * 3;
        model->Q[i * n + i]         = q00;
        model->Q[i * n + i + 1]     = q01;
        model->Q[i * n + i + 2]     = q02;
        model->Q[(i + 1) * n + i]   = q01;
        model->Q[(i + 1) * n + i + 1] = q11;
        model->Q[(i + 1) * n + i + 2] = q12;
        model->Q[(i + 2) * n + i]   = q02;
        model->Q[(i + 2) * n + i + 1] = q12;
        model->Q[(i + 2) * n + i + 2] = q22;
    }

    return model;
}

/* ============================================================================
 * Coordinated Turn (CT) model
 * ============================================================================
 */

motion_model_t *motion_create_ct(double omega, double T, double q)
{
    if (T <= 0.0 || q < 0.0) return NULL;

    int n = 4; /* [x, vx, y, vy] — 2D with constant turn rate */
    /* Note: CT with known omega; 5-state version adds omega estimation */

    motion_model_t *model = (motion_model_t *)calloc(1, sizeof(motion_model_t));
    if (!model) return NULL;
    model->type = MOTION_CT;
    model->dim = DIM_2D;
    model->state_dim = n;
    model->T = T;
    model->q = q;
    model->omega = omega;

    model->F = (double *)calloc(n * n, sizeof(double));
    model->Q = (double *)calloc(n * n, sizeof(double));
    if (!model->F || !model->Q) { motion_free(model); return NULL; }

    if (fabs(omega) < 1e-10) {
        /* Degenerate to CV when omega ≈ 0 */
        model->F[0] = 1.0; model->F[1] = T;
        model->F[2*4+2] = 1.0; model->F[2*4+3] = T;
        model->F[1*4+1] = 1.0; model->F[3*4+3] = 1.0;
    } else {
        double sin_wT = sin(omega * T);
        double cos_wT = cos(omega * T);

        /* Row 0: [1 sin(ωT)/ω 0 −(1−cos(ωT))/ω] */
        model->F[0 * n + 0] = 1.0;
        model->F[0 * n + 1] = sin_wT / omega;
        model->F[0 * n + 2] = 0.0;
        model->F[0 * n + 3] = -(1.0 - cos_wT) / omega;

        /* Row 1: [0 cos(ωT) 0 −sin(ωT)] */
        model->F[1 * n + 1] = cos_wT;
        model->F[1 * n + 3] = -sin_wT;

        /* Row 2: [0 (1−cos(ωT))/ω 1 sin(ωT)/ω] */
        model->F[2 * n + 1] = (1.0 - cos_wT) / omega;
        model->F[2 * n + 2] = 1.0;
        model->F[2 * n + 3] = sin_wT / omega;

        /* Row 3: [0 sin(ωT) 0 cos(ωT)] */
        model->F[3 * n + 1] = sin_wT;
        model->F[3 * n + 3] = cos_wT;
    }

    /* White acceleration process noise (same as CV Q blocks) */
    double T3_3 = T * T * T / 3.0;
    double T2_2 = T * T / 2.0;

    model->Q[0 * n + 0] = q * T3_3;
    model->Q[0 * n + 1] = q * T2_2;
    model->Q[1 * n + 0] = q * T2_2;
    model->Q[1 * n + 1] = q * T;
    model->Q[2 * n + 2] = q * T3_3;
    model->Q[2 * n + 3] = q * T2_2;
    model->Q[3 * n + 2] = q * T2_2;
    model->Q[3 * n + 3] = q * T;

    return model;
}

/* ============================================================================
 * Singer maneuver model
 * ============================================================================
 */

motion_model_t *motion_create_singer(motion_dim_t dim, double T,
                                      double alpha, double sigma_m)
{
    if (T <= 0.0) return NULL;

    /* Singer model is 1D: [x, v, a] per axis (but applied per-axis) */
    int state_per_axis = 3;
    int n = (int)dim * state_per_axis;

    motion_model_t *model = (motion_model_t *)calloc(1, sizeof(motion_model_t));
    if (!model) return NULL;
    model->type = MOTION_SINGER;
    model->dim = dim;
    model->state_dim = n;
    model->T = T;
    model->alpha = alpha;
    model->sigma_m = sigma_m;

    model->F = (double *)calloc(n * n, sizeof(double));
    model->Q = (double *)calloc(n * n, sizeof(double));
    if (!model->F || !model->Q) { motion_free(model); return NULL; }

    double exp_neg_aT = exp(-alpha * T);
    double aT = alpha * T;

    /* F per axis = [1 T (aT − 1 + e^{−αT})/α²;
     *                0 1 (1 − e^{−αT})/α;
     *                0 0 e^{−αT}              ]
     */
    for (int d = 0; d < (int)dim; d++) {
        int i = d * 3;
        model->F[i * n + i]     = 1.0;
        model->F[i * n + i + 1] = T;
        model->F[i * n + i + 2] = (aT - 1.0 + exp_neg_aT) / (alpha * alpha);
        model->F[(i + 1) * n + i + 1] = 1.0;
        model->F[(i + 1) * n + i + 2] = (1.0 - exp_neg_aT) / alpha;
        model->F[(i + 2) * n + i + 2] = exp_neg_aT;
    }

    /* Q = 2ασ_m² * q_singer (approximate as σ_m² * Q_simple) */
    double sig2 = sigma_m * sigma_m;
    double q_scale = 2.0 * alpha * sig2;

    /* Simplify: use a reasonable Q approximation */
    double T5_20 = T * T * T * T * T / 20.0;
    double T4_8  = T * T * T * T / 8.0;
    double T3_6  = T * T * T / 6.0;
    double T3_3  = T3_6 * 2.0;
    double T2_2  = T * T / 2.0;

    for (int d = 0; d < (int)dim; d++) {
        int i = d * 3;
        /* Scaled by q_scale in the acceleration dimension primarily */
        model->Q[i * n + i]         = q_scale * T5_20;
        model->Q[i * n + i + 1]     = q_scale * T4_8;
        model->Q[i * n + i + 2]     = q_scale * T3_6;
        model->Q[(i + 1) * n + i]   = q_scale * T4_8;
        model->Q[(i + 1) * n + i + 1] = q_scale * T3_3;
        model->Q[(i + 1) * n + i + 2] = q_scale * T2_2;
        model->Q[(i + 2) * n + i]   = q_scale * T3_6;
        model->Q[(i + 2) * n + i + 1] = q_scale * T2_2;
        model->Q[(i + 2) * n + i + 2] = q_scale * T;
    }

    return model;
}

/* ============================================================================
 * Wiener Process Acceleration (WPA / white jerk)
 * ============================================================================
 */

motion_model_t *motion_create_wpa(motion_dim_t dim, double T, double q)
{
    /* WPA is similar to CA but with white jerk driving the acceleration.
     * We implement as a 3rd-order (CA) model.
     */
    return motion_create_ca(dim, T, q);
}

/* ============================================================================
 * Brownian motion / random walk
 * ============================================================================
 */

motion_model_t *motion_create_brownian(motion_dim_t dim, double T, double q)
{
    if (T <= 0.0 || q < 0.0) return NULL;

    int n = (int)dim; /* position only per axis */

    motion_model_t *model = (motion_model_t *)calloc(1, sizeof(motion_model_t));
    if (!model) return NULL;
    model->type = MOTION_BROWN;
    model->dim = dim;
    model->state_dim = n;
    model->T = T;
    model->q = q;

    model->F = (double *)calloc(n * n, sizeof(double));
    model->Q = (double *)calloc(n * n, sizeof(double));
    if (!model->F || !model->Q) { motion_free(model); return NULL; }

    /* F = I (identity) */
    for (int i = 0; i < n; i++) model->F[i * n + i] = 1.0;

    /* Q = q*T*I */
    for (int i = 0; i < n; i++) model->Q[i * n + i] = q * T;

    return model;
}

void motion_free(motion_model_t *model)
{
    if (!model) return;
    free(model->F);
    free(model->Q);
    free(model->Gamma);
    free(model);
}

/* ============================================================================
 * State and covariance prediction
 * ============================================================================
 */

void motion_predict_state(const motion_model_t *model,
                           const double *x, double *x_next)
{
    if (!model || !x || !x_next) return;
    int n = model->state_dim;
    mat_vec_mul(model->F, x, x_next, n, n);
}

void motion_predict_covariance(const motion_model_t *model,
                                const double *P, double *P_next)
{
    if (!model || !P || !P_next) return;
    int n = model->state_dim;

    double FP[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(model->F, P, FP, n, n, n);

    double FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(model->F, FT, n, n);

    mat_mat_mul(FP, FT, P_next, n, n, n);

    for (int i = 0; i < n * n; i++) {
        P_next[i] += model->Q[i];
    }
}

void motion_predict(const motion_model_t *model,
                     const double *x, const double *P,
                     double *x_next, double *P_next)
{
    motion_predict_state(model, x, x_next);
    motion_predict_covariance(model, P, P_next);
}

/* ============================================================================
 * Measurement matrix construction
 * ============================================================================
 */

int motion_build_H(const motion_model_t *model, measurement_type_t meas_type,
                    double *H, int meas_dim)
{
    if (!model || !H) return -1;
    int n = model->state_dim;

    memset(H, 0, meas_dim * n * sizeof(double));

    switch (meas_type) {
    case MEAS_CARTESIAN_2D:
        /* Measure position x and y */
        if (model->dim >= DIM_2D && meas_dim >= 2) {
            H[0 * n + 0] = 1.0; /* x position */
            H[1 * n + 2] = 1.0; /* y position */
        }
        break;
    case MEAS_CARTESIAN_3D:
        if (model->dim >= DIM_3D && meas_dim >= 3) {
            H[0 * n + 0] = 1.0;
            H[1 * n + 3] = 1.0;
            H[2 * n + 6] = 1.0;
        }
        break;
    case MEAS_POLAR_2D:
        /* Position measurement, but note: this H is used after conversion */
        if (model->dim >= DIM_2D && meas_dim >= 2) {
            H[0 * n + 0] = 1.0;
            H[1 * n + 2] = 1.0;
        }
        break;
    case MEAS_DOPPLER_ONLY:
        /* Measure velocity in x (or radial, depending on setup) */
        if (model->state_dim >= 2 && meas_dim >= 1) {
            H[0 * n + 1] = 1.0; /* x velocity */
        }
        break;
    default:
        return -1;
    }
    return 0;
}

/* ============================================================================
 * IMM model set construction
 * ============================================================================
 */

int imm_build_standard_3model(motion_model_t **models, double *mixing,
                               int *max_models, double T,
                               double q_cv, double q_ct, double omega)
{
    if (!models || !mixing || !max_models || *max_models < 3) return -1;

    models[0] = motion_create_cv(DIM_2D, T, q_cv);
    models[1] = motion_create_ct(omega, T, q_ct);
    models[2] = motion_create_ct(-omega, T, q_ct);

    if (!models[0] || !models[1] || !models[2]) {
        motion_free(models[0]); motion_free(models[1]); motion_free(models[2]);
        return -1;
    }

    /* Markov transition matrix (3×3):
     *   [0.95  0.025 0.025;
     *    0.025 0.95  0.025;
     *    0.025 0.025 0.95 ]
     */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mixing[i * 3 + j] = (i == j) ? 0.95 : 0.025;
        }
    }

    *max_models = 3;
    return 3;
}

int imm_build_vehicle_2model(motion_model_t **models, double *mixing,
                              double T, double q_moving, double q_stopped)
{
    if (!models || !mixing) return -1;

    models[0] = motion_create_cv(DIM_2D, T, q_moving);
    models[1] = motion_create_brownian(DIM_2D, T, q_stopped);

    if (!models[0] || !models[1]) {
        motion_free(models[0]); motion_free(models[1]);
        return -1;
    }

    mixing[0] = 0.9; mixing[1] = 0.1; /* stay moving / stop */
    mixing[2] = 0.3; mixing[3] = 0.7; /* start moving / stay stopped */

    return 2;
}

/* ============================================================================
 * Continuous-time discretization
 * ============================================================================
 */

void motion_discretize(const double *A, const double *G, const double *Q_c,
                        int state_dim, int noise_dim, double T,
                        double *F, double *Q_d, int terms)
{
    if (!A || !F || !Q_d) return;
    int n = state_dim;

    if (terms < 1) terms = 8;

    /* F = exp(A*T) via truncated Taylor series */
    matrix_exponential(A, F, n, terms);

    /* Scale: exp(A*T) with T scaling */
    double A_scaled[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    for (int i = 0; i < n * n; i++) A_scaled[i] = A[i] * T;
    matrix_exponential(A_scaled, F, n, terms);

    /* Q_d = ∫₀ᵀ exp(A·τ)·G·Q_c·Gᵀ·exp(Aᵀ·τ) dτ
     * Approximate via Van Loan's method or discretized integration.
     *
     * Simplified: Q_d ≈ G·Q_c·Gᵀ·T + (A·G·Q_c·Gᵀ + G·Q_c·Gᵀ·Aᵀ)·T²/2
     * This is a first-order approximation.
     */
    if (G && Q_c && noise_dim > 0) {
        double GQG[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        memset(GQG, 0, sizeof(GQG));

        /* GQG = G * Q_c * G' */
        double GQ[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
        mat_mat_mul(G, Q_c, GQ, n, noise_dim, noise_dim);
        double GT[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
        mat_transpose(G, GT, n, noise_dim);
        mat_mat_mul(GQ, GT, GQG, n, noise_dim, n);

        /* First term */
        for (int i = 0; i < n * n; i++) Q_d[i] = GQG[i] * T;

        /* Higher-order terms (second order) */
        double AGQG[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        mat_mat_mul(A_scaled, GQG, AGQG, n, n, n);
        double GQGAT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        double AT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        mat_transpose(A_scaled, AT, n, n);
        mat_mat_mul(GQG, AT, GQGAT, n, n, n);

        for (int i = 0; i < n * n; i++) {
            Q_d[i] += (AGQG[i] + GQGAT[i]) * T / 2.0;
        }
    }
}

void matrix_exponential(const double *A, double *expA, int n, int terms)
{
    /* exp(A) = Σ_{k=0}^{terms-1} A^k / k!
     *
     * Use scaling-and-squaring for better numerical accuracy when ||A|| is large.
     * For small ||A||, the Taylor series converges quickly.
     */

    /* Initialize as identity */
    memset(expA, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) expA[i * n + i] = 1.0;

    double A_pow[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memcpy(A_pow, A, n * n * sizeof(double));

    double factorial = 1.0;

    for (int k = 1; k < terms; k++) {
        factorial *= k;
        for (int i = 0; i < n * n; i++) {
            expA[i] += A_pow[i] / factorial;
        }

        /* A_pow = A_pow * A for next iteration */
        double A_next[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        mat_mat_mul(A_pow, A, A_next, n, n, n);
        memcpy(A_pow, A_next, n * n * sizeof(double));
    }
}
