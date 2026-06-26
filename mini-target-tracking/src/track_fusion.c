/**
 * track_fusion.c — Multi-Sensor Track-to-Track Fusion Implementation
 *
 * Implements: Track association, CI, WCF, Info Matrix Fusion,
 *             Tracklet fusion, cross-covariance, decorrelation.
 *
 * References:
 *   - Bar-Shalom (1981), Julier & Uhlmann (1997)
 *   - Chang, Saha, Bar-Shalom (1997)
 *   - Chong, Mori, Chang (1990)
 */

#include "track_fusion.h"
#include "kalman_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Track-to-track association
 * ============================================================================
 */

t2t_association_t *t2t_associate(const track_t *local_tracks, int n_local,
                                  const track_t *remote_tracks, int n_remote,
                                  double threshold)
{
    if (!local_tracks || !remote_tracks || n_local <= 0 || n_remote <= 0)
        return NULL;

    t2t_association_t *assoc = (t2t_association_t *)calloc(1, sizeof(t2t_association_t));
    if (!assoc) return NULL;

    assoc->n_local = n_local;
    assoc->n_remote = n_remote;
    assoc->assignment = (int *)malloc(n_local * sizeof(int));
    assoc->association_cost = (double *)malloc(n_local * sizeof(double));
    if (!assoc->assignment || !assoc->association_cost) {
        t2t_association_free(assoc); return NULL;
    }

    for (int i = 0; i < n_local; i++) assoc->assignment[i] = -1;

    /* Greedy NN association: for each local track, find closest remote */
    int *remote_taken = (int *)calloc(n_remote, sizeof(int));
    if (!remote_taken) { t2t_association_free(assoc); return NULL; }

    for (int i = 0; i < n_local; i++) {
        const track_t *lt = &local_tracks[i];
        if (lt->state == TRACK_STATE_FREE || lt->state == TRACK_STATE_DELETED)
            continue;

        double min_dist = 1e15;
        int best_j = -1;

        for (int j = 0; j < n_remote; j++) {
            if (remote_taken[j]) continue;
            const track_t *rt = &remote_tracks[j];
            if (rt->state == TRACK_STATE_FREE || rt->state == TRACK_STATE_DELETED)
                continue;
            if (lt->state_dim != rt->state_dim) continue;

            int n = lt->state_dim;
            /* d² = (x1 - x2)' * (P1 + P2)^{-1} * (x1 - x2) */
            double dx[TRACK_MAX_STATE_DIM];
            vec_sub(lt->x, rt->x, dx, n);

            double P_sum[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
            for (int k = 0; k < n * n; k++) P_sum[k] = lt->P[k] + rt->P[k];

            /* Solve P_sum * t = dx */
            double t[TRACK_MAX_STATE_DIM];
            vec_copy(dx, t, n);
            if (mat_solve_cholesky(P_sum, dx, t, n) == 0) {
                double d2 = vec_dot(dx, t, n); /* dx' * P_sum^{-1} * dx */
                if (d2 < min_dist) {
                    min_dist = d2;
                    best_j = j;
                }
            }
        }

        if (best_j >= 0 && min_dist <= threshold) {
            assoc->assignment[i] = best_j;
            assoc->association_cost[i] = min_dist;
            remote_taken[best_j] = 1;
            assoc->n_pairs++;
        }
    }

    free(remote_taken);
    return assoc;
}

void t2t_association_free(t2t_association_t *assoc)
{
    if (!assoc) return;
    free(assoc->assignment);
    free(assoc->association_cost);
    free(assoc);
}

/* ============================================================================
 * Covariance Intersection (CI)
 * ============================================================================
 */

int fusion_covariance_intersection(const double *x1, const double *P1,
                                    const double *x2, const double *P2,
                                    int n, double *x_fused, double *P_fused,
                                    double *omega)
{
    if (!x1 || !P1 || !x2 || !P2 || !x_fused || !P_fused) return -1;
    if (n <= 0 || n > TRACK_MAX_STATE_DIM) return -1;

    /* Compute P1^{-1} and P2^{-1} */
    double P1_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    double P2_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    if (mat_inv_cholesky(P1, P1_inv, n) <= 0.0) return -1;
    if (mat_inv_cholesky(P2, P2_inv, n) <= 0.0) return -1;

    /* Optimize omega ∈ [0, 1] by golden section search minimizing det(P_CI) */
    double w_best = 0.5;
    double det_best = 1e30;

    for (int k = 0; k <= 20; k++) {
        double w = k / 20.0; /* coarse grid search */

        /* P_CI^{-1} = w * P1^{-1} + (1-w) * P2^{-1} */
        double P_CI_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        for (int i = 0; i < n * n; i++) {
            P_CI_inv[i] = w * P1_inv[i] + (1.0 - w) * P2_inv[i];
        }

        /* Invert to get P_CI */
        double P_CI[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        double det = mat_inv_cholesky(P_CI_inv, P_CI, n);

        if (det > 0.0 && det < det_best) {
            det_best = det;
            w_best = w;

            /* Store this P_CI */
            memcpy(P_fused, P_CI, n * n * sizeof(double));
        }
    }

    /* More precision near optimum */
    for (int k = 0; k <= 100; k++) {
        double w = (k / 100.0) * (w_best + 0.05) + (w_best - 0.05);
        if (w < 0.0) w = 0.0;
        if (w > 1.0) w = 1.0;

        double P_CI_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        for (int i = 0; i < n * n; i++) {
            P_CI_inv[i] = w * P1_inv[i] + (1.0 - w) * P2_inv[i];
        }

        double P_CI[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        double det = mat_inv_cholesky(P_CI_inv, P_CI, n);
        if (det > 0.0 && det < det_best) {
            det_best = det;
            w_best = w;
            memcpy(P_fused, P_CI, n * n * sizeof(double));
        }
    }

    /* Fused state: x_CI = P_CI * (w * P1^{-1} * x1 + (1-w) * P2^{-1} * x2) */
    double P1_inv_x1[TRACK_MAX_STATE_DIM];
    double P2_inv_x2[TRACK_MAX_STATE_DIM];
    mat_vec_mul(P1_inv, x1, P1_inv_x1, n, n);
    mat_vec_mul(P2_inv, x2, P2_inv_x2, n, n);

    double combined[TRACK_MAX_STATE_DIM];
    for (int i = 0; i < n; i++) {
        combined[i] = w_best * P1_inv_x1[i] + (1.0 - w_best) * P2_inv_x2[i];
    }
    mat_vec_mul(P_fused, combined, x_fused, n, n);

    if (omega) *omega = w_best;
    return 0;
}

/* ============================================================================
 * Weighted Covariance Fusion (WCF)
 * ============================================================================
 */

int fusion_weighted_covariance(const double *x1, const double *P1,
                                const double *x2, const double *P2,
                                const double *P12,
                                int n, double *x_fused, double *P_fused)
{
    if (!x1 || !P1 || !x2 || !P2 || !P12 || !x_fused || !P_fused) return -1;

    /* P_sum = P1 + P2 - P12 - P12' */
    double P_sum[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(P_sum, 0, sizeof(P_sum));
    for (int i = 0; i < n * n; i++) {
        double P12_T;
        /* P12 is n×n, P12'[i] = P12[j*n+i] where i = row, j = col */
        int row = i / n, col = i % n;
        P12_T = P12[col * n + row];
        P_sum[i] = P1[i] + P2[i] - P12[i] - P12_T;
    }

    /* P_diff = P1 - P12 */
    double P_diff[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    for (int i = 0; i < n * n; i++) P_diff[i] = P1[i] - P12[i];

    /* Invert P_sum */
    double P_sum_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    if (mat_inv_cholesky(P_sum, P_sum_inv, n) <= 0.0) return -1;

    /* x_fused = x1 + P_diff * P_sum^{-1} * (x2 - x1) */
    double dx[TRACK_MAX_STATE_DIM];
    vec_sub(x2, x1, dx, n);

    double P_sum_inv_dx[TRACK_MAX_STATE_DIM];
    mat_vec_mul(P_sum_inv, dx, P_sum_inv_dx, n, n);

    double correction[TRACK_MAX_STATE_DIM];
    mat_vec_mul(P_diff, P_sum_inv_dx, correction, n, n);

    vec_add(x1, correction, x_fused, n);

    /* P_fused = P1 - P_diff * P_sum^{-1} * P_diff' */
    double P_diff_T[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(P_diff, P_diff_T, n, n);

    double temp[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(P_diff, P_sum_inv, temp, n, n, n);
    double temp2[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(temp, P_diff_T, temp2, n, n, n);

    for (int i = 0; i < n * n; i++) P_fused[i] = P1[i] - temp2[i];

    return 0;
}

/* ============================================================================
 * Information Matrix Fusion
 * ============================================================================
 */

int fusion_information_matrix(const double *Y1, const double *y1,
                               const double *Y2, const double *y2,
                               const double *Y_common, const double *y_common,
                               int n,
                               double *Y_fused, double *y_fused,
                               double *x_fused, double *P_fused)
{
    if (!Y1 || !y1 || !Y2 || !y2 || !Y_fused || !y_fused) return -1;

    /* Y_fused = Y1 + Y2 - Y_common */
    for (int i = 0; i < n * n; i++) {
        Y_fused[i] = Y1[i] + Y2[i];
        if (Y_common) Y_fused[i] -= Y_common[i];
    }

    /* y_fused = y1 + y2 - y_common */
    for (int i = 0; i < n; i++) {
        y_fused[i] = y1[i] + y2[i];
        if (y_common) y_fused[i] -= y_common[i];
    }

    /* Convert to Cartesian if requested: x = Y^{-1} * y, P = Y^{-1} */
    if (x_fused || P_fused) {
        double Y_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        if (mat_inv_cholesky(Y_fused, Y_inv, n) <= 0.0) return -1;

        if (P_fused) memcpy(P_fused, Y_inv, n * n * sizeof(double));
        if (x_fused) mat_vec_mul(Y_inv, y_fused, x_fused, n, n);
    }

    return 0;
}

/* ============================================================================
 * Tracklet fusion
 * ============================================================================
 */

int fusion_tracklet(double *x1, double *P1,
                     const double *z2, const double *R2,
                     const double *F, const double *Q, const double *H,
                     int n, int m, int k2, double dt)
{
    (void)dt; /* time step reserved for future variable-rate support */
    if (!x1 || !P1 || !z2 || !F || !Q || !H) return -1;

    /* Process each measurement from tracklet 2 sequentially */
    double x[TRACK_MAX_STATE_DIM];
    double P[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memcpy(x, x1, n * sizeof(double));
    memcpy(P, P1, n * n * sizeof(double));

    for (int k = 0; k < k2; k++) {
        /* Predict */
        double x_pred[TRACK_MAX_STATE_DIM];
        mat_vec_mul(F, x, x_pred, n, n);

        double FP[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        mat_mat_mul(F, P, FP, n, n, n);
        double FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        mat_transpose(F, FT, n, n);
        double P_pred[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        mat_mat_mul(FP, FT, P_pred, n, n, n);
        for (int i = 0; i < n * n; i++) P_pred[i] += Q[i];

        /* Update with measurement z2[k*m ... k*m+m-1] */
        const double *zk = &z2[k * m];
        const double *Rk = &R2[k * m * m]; /* or same R for all */

        double nu[TRACK_MAX_MEAS_DIM];
        double Hx_pred[TRACK_MAX_MEAS_DIM];
        mat_vec_mul(H, x_pred, Hx_pred, m, n);
        vec_sub(zk, Hx_pred, nu, m);

        double PHT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
        double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
        mat_transpose(H, HT, m, n);
        mat_mat_mul(P_pred, HT, PHT, n, n, m);

        double S[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
        mat_mat_mul(H, PHT, S, m, n, m);
        for (int i = 0; i < m * m; i++) S[i] += Rk[i];

        double S_inv[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
        if (mat_inv_cholesky(S, S_inv, m) <= 0.0) continue;

        double K[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
        mat_mat_mul(PHT, S_inv, K, n, m, m);

        double Knu[TRACK_MAX_STATE_DIM];
        mat_vec_mul(K, nu, Knu, n, m);
        vec_add(x_pred, Knu, x, n);

        double KT[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
        mat_transpose(K, KT, n, m);
        double KSKT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
        double KS[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
        mat_mat_mul(K, S, KS, n, m, m);
        mat_mat_mul(KS, KT, KSKT, n, m, n);

        for (int i = 0; i < n * n; i++) {
            P[i] = P_pred[i] - KSKT[i];
        }
    }

    memcpy(x1, x, n * sizeof(double));
    memcpy(P1, P, n * n * sizeof(double));
    return 0;
}

/* ============================================================================
 * Cross-covariance and decorrelation
 * ============================================================================
 */

void fusion_cross_covariance_update(double *P12,
                                     const double *K1, const double *K2,
                                     const double *F, const double *Q,
                                     const double *H, int n, int m)
{
    if (!P12 || !K1 || !K2 || !F || !Q || !H) return;

    /* I - K1*H and I - K2*H */
    double I_K1H[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    double I_K2H[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(I_K1H, 0, sizeof(I_K1H));
    memset(I_K2H, 0, sizeof(I_K2H));

    for (int i = 0; i < n; i++) {
        I_K1H[i * n + i] = 1.0;
        I_K2H[i * n + i] = 1.0;
    }

    double K1H[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(K1, H, K1H, n, m, n);
    double K2H[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(K2, H, K2H, n, m, n);

    for (int i = 0; i < n * n; i++) {
        I_K1H[i] -= K1H[i];
        I_K2H[i] -= K2H[i];
    }

    /* P12_new = I_K1H * F * P12_current * F' * I_K2H' + I_K1H * Q * I_K2H' */
    double F_P12[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(F, P12, F_P12, n, n, n);
    double FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(F, FT, n, n);
    double F_P12_FT[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(F_P12, FT, F_P12_FT, n, n, n);

    double I1_FPF[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(I_K1H, F_P12_FT, I1_FPF, n, n, n);

    double I_K2H_T[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(I_K2H, I_K2H_T, n, n);
    double term1[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(I1_FPF, I_K2H_T, term1, n, n, n);

    /* I_K1H * Q * I_K2H' */
    double I1_Q[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(I_K1H, Q, I1_Q, n, n, n);
    double term2[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(I1_Q, I_K2H_T, term2, n, n, n);

    for (int i = 0; i < n * n; i++) P12[i] = term1[i] + term2[i];
}

void fusion_decorrelate(const double *x1, const double *P1,
                         const double *x2, const double *P2,
                         const double *P12, int n,
                         double *z_dec, double *R_dec)
{
    if (!x1 || !P1 || !x2 || !P2 || !P12 || !z_dec || !R_dec) return;

    /* z_dec = x2 - P12' * P1^{-1} * x1 */
    double P1_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    if (mat_inv_cholesky(P1, P1_inv, n) <= 0.0) return;

    double P1_inv_x1[TRACK_MAX_STATE_DIM];
    mat_vec_mul(P1_inv, x1, P1_inv_x1, n, n);

    double P12_T[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_transpose(P12, P12_T, n, n);

    double P12T_P1inv_x1[TRACK_MAX_STATE_DIM];
    mat_vec_mul(P12_T, P1_inv_x1, P12T_P1inv_x1, n, n);

    vec_sub(x2, P12T_P1inv_x1, z_dec, n);

    /* R_dec = P2 - P12' * P1^{-1} * P12 */
    double P1_inv_P12[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(P1_inv, P12, P1_inv_P12, n, n, n);

    double temp[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(P12_T, P1_inv_P12, temp, n, n, n);

    for (int i = 0; i < n * n; i++) R_dec[i] = P2[i] - temp[i];
}

/* ============================================================================
 * Rumor-robust novelty detection
 * ============================================================================
 */

int fusion_is_novel_info(const double *Y_local, const double *Y_remote,
                          int n, double threshold)
{
    if (!Y_local || !Y_remote) return 0;

    /* d(Y1, Y2) = tr(Y1^{-1} * Y2) - log|Y1^{-1} * Y2| - n
     * = tr(Y1^{-1} * Y2) + log|Y1| - log|Y2| - n
     */

    double Y1_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    if (mat_inv_cholesky(Y_local, Y1_inv, n) <= 0.0) return 1; /* assume novel */

    double Y1_inv_Y2[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(Y1_inv, Y_remote, Y1_inv_Y2, n, n, n);

    /* Trace */
    double trace = 0.0;
    for (int i = 0; i < n; i++) trace += Y1_inv_Y2[i * n + i];

    /* Log-determinants */
    double log_det_Y1 = log(mat_det_cholesky(Y_local, n));
    double log_det_Y2 = log(mat_det_cholesky(Y_remote, n));

    double d = trace + log_det_Y1 - log_det_Y2 - n;
    if (d < 0.0) d = 0.0;

    return (d > threshold) ? 1 : 0;
}
