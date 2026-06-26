/**
 * track_core.c — Track Management, Gating, Scoring Implementation
 *
 * References:
 *   - Bar-Shalom, Willett, Tian "Tracking and Data Fusion" (2011)
 *   - Reid, D.B. "An Algorithm for Tracking Multiple Targets" (1979)
 *   - Blackman, Popoli "Design and Analysis of Modern Tracking Systems" (1999)
 */

#include "track_core.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Forward reference to motion model ID (defined in motion_models.h) */
#ifndef MOTION_CV
#define MOTION_CV 0
#endif

/* ============================================================================
 * Tracker lifecycle management
 * ============================================================================
 */

void tracker_init(tracker_t *tracker)
{
    if (!tracker) return;
    memset(tracker, 0, sizeof(tracker_t));

    /* Default M/N logic: 3 detections in 5 scans to confirm */
    tracker->mgmt.M = 3;
    tracker->mgmt.N = 5;
    tracker->mgmt.coast_max = 5;
    tracker->mgmt.delete_score_threshold = -10.0;
    tracker->mgmt.confirm_score_threshold = 2.0;

    /* Default gating: 99% ellipsoidal gate */
    tracker->default_gate.shape = GATE_ELLIPSOIDAL;
    tracker->default_gate.gate_size = 0.99;
    tracker->default_gate.gate_dof = 2;

    tracker->next_track_id = 1;
    tracker->num_tracks = 0;
    tracker->current_time = 0.0;
    tracker->scan_count = 0;

    /* Initialize all track slots as FREE */
    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        tracker->tracks[i].state = TRACK_STATE_FREE;
        tracker->tracks[i].track_id = TRACK_INVALID_ID;
    }
}

void track_init(track_t *track, uint32_t id, int state_dim, int meas_dim)
{
    if (!track) return;
    memset(track, 0, sizeof(track_t));

    track->track_id = id;
    track->state = TRACK_STATE_FREE;
    track->state_dim = state_dim;
    track->meas_dim = meas_dim;
    track->meas_type = MEAS_CARTESIAN_2D;
    track->motion_model_id = MOTION_CV;
    track->score = 0.0;

    /* Default gate */
    track->gate.shape = GATE_ELLIPSOIDAL;
    track->gate.gate_size = 0.99;
    track->gate.gate_dof = meas_dim;
    track->gate.gate_threshold = chi2_threshold_approx(0.99, meas_dim);

    vec_zero(track->x, TRACK_MAX_STATE_DIM);
    vec_zero(track->P, TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM);
    /* Initialize P as large identity (high initial uncertainty) */
    for (int i = 0; i < state_dim; i++) {
        track->P[i * TRACK_MAX_STATE_DIM + i] = 1e6;
    }
}

track_t *tracker_create_track(tracker_t *tracker, int state_dim, int meas_dim,
                               measurement_type_t meas_type)
{
    if (!tracker) return NULL;
    if (state_dim > TRACK_MAX_STATE_DIM || meas_dim > TRACK_MAX_MEAS_DIM)
        return NULL;

    /* Find first free slot */
    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        if (tracker->tracks[i].state == TRACK_STATE_FREE ||
            tracker->tracks[i].state == TRACK_STATE_DELETED) {
            track_t *t = &tracker->tracks[i];
            track_init(t, tracker->next_track_id++, state_dim, meas_dim);
            t->state = TRACK_STATE_TENTATIVE;
            t->meas_type = meas_type;
            t->birth_time = tracker->current_time;
            t->last_update_time = tracker->current_time;
            t->gate = tracker->default_gate;
            t->gate.gate_dof = meas_dim;
            t->gate.gate_threshold = chi2_threshold_approx(
                tracker->default_gate.gate_size, meas_dim);

            tracker->num_tracks++;
            tracker->total_tracks_created++;
            return t;
        }
    }
    return NULL; /* Tracker full */
}

void track_delete(track_t *track)
{
    if (!track) return;
    track->state = TRACK_STATE_DELETED;
}

int tracker_cleanup(tracker_t *tracker)
{
    if (!tracker) return 0;
    int removed = 0;

    /* Compact: move non-DELETED tracks to front */
    int write_idx = 0;
    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        if (tracker->tracks[i].state != TRACK_STATE_DELETED &&
            tracker->tracks[i].state != TRACK_STATE_FREE) {
            if (write_idx != i) {
                tracker->tracks[write_idx] = tracker->tracks[i];
                /* Clear the old slot */
                tracker->tracks[i].state = TRACK_STATE_FREE;
                tracker->tracks[i].track_id = TRACK_INVALID_ID;
            }
            write_idx++;
        } else if (tracker->tracks[i].state == TRACK_STATE_DELETED) {
            tracker->tracks[i].state = TRACK_STATE_FREE;
            tracker->tracks[i].track_id = TRACK_INVALID_ID;
            removed++;
            tracker->total_tracks_deleted++;
        }
    }
    tracker->num_tracks = write_idx;
    return removed;
}

void tracker_apply_mn_logic(tracker_t *tracker)
{
    if (!tracker) return;
    int N = tracker->mgmt.N;
    int M = tracker->mgmt.M;

    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        track_t *t = &tracker->tracks[i];
        if (t->state != TRACK_STATE_TENTATIVE && t->state != TRACK_STATE_CONFIRMED)
            continue;

        /* Update sliding window: shift left, add current bit */
        for (int j = N - 1; j > 0; j--) {
            t->history_window[j] = t->history_window[j - 1];
        }
        t->history_window[0] = (t->misses_since_last == 0) ? 1 : 0;

        /* Count detections in window */
        int detections = 0;
        int scans = 0;
        for (int j = 0; j < N; j++) {
            if (t->history_window[j] >= 0) {
                detections += t->history_window[j];
                scans++;
            }
        }

        /* M/N logic */
        if (t->state == TRACK_STATE_TENTATIVE && detections >= M) {
            t->state = TRACK_STATE_CONFIRMED;
            tracker->total_tracks_confirmed++;
        }

        /* Store for external use */
        t->confirm_counter = detections;
    }
}

int tracker_prune_low_score(tracker_t *tracker)
{
    if (!tracker) return 0;
    int pruned = 0;
    double threshold = tracker->mgmt.delete_score_threshold;

    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        if (tracker->tracks[i].track_id == TRACK_INVALID_ID) continue;
        if (tracker->tracks[i].state == TRACK_STATE_FREE ||
            tracker->tracks[i].state == TRACK_STATE_DELETED) continue;

        if (tracker->tracks[i].score < threshold) {
            track_delete(&tracker->tracks[i]);
            pruned++;
        }
    }
    if (pruned > 0) {
        tracker_cleanup(tracker);
    }
    return pruned;
}

void tracker_coast_orphan_tracks(tracker_t *tracker)
{
    if (!tracker) return;
    int coast_max = tracker->mgmt.coast_max;

    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        track_t *t = &tracker->tracks[i];
        if (t->track_id == TRACK_INVALID_ID) continue;
        if (t->state == TRACK_STATE_FREE || t->state == TRACK_STATE_DELETED)
            continue;

        if (t->misses_since_last > 0) {
            /* Coast: track was NOT associated in this scan */
            if (t->state == TRACK_STATE_CONFIRMED) {
                t->state = TRACK_STATE_COAST;
            }
            /* Check if coasting too long */
            if (t->misses_since_last > coast_max) {
                t->state = TRACK_STATE_DELETED;
            }
        }
    }
    /* Increment miss counters for all active tracks */
    /* Caller should reset misses_since_last for associated tracks */
}

/* ============================================================================
 * Mahalanobis distance and gating
 * ============================================================================
 */

double mahalanobis_distance_sq(const double *z, const double *R,
                                const double *x_pred, const double *P_pred,
                                const double *H, int state_dim, int meas_dim)
{
    if (!z || !R || !x_pred || !P_pred || !H) return -1.0;
    if (state_dim <= 0 || meas_dim <= 0) return -1.0;
    if (state_dim > TRACK_MAX_STATE_DIM || meas_dim > TRACK_MAX_MEAS_DIM)
        return -1.0;

    /* Innovation: nu = z - H*x_pred */
    double nu[TRACK_MAX_MEAS_DIM];
    double Hx[TRACK_MAX_MEAS_DIM];
    mat_vec_mul(H, x_pred, Hx, meas_dim, state_dim);
    vec_sub(z, Hx, nu, meas_dim);

    /* Innovation covariance: S = H*P_pred*H' + R */
    double PHT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    double S[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];

    /* PHT = P_pred * H' */
    double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    for (int i = 0; i < state_dim; i++) {
        for (int j = 0; j < meas_dim; j++) {
            HT[j * state_dim + i] = H[i * state_dim + j]; /* H' */
        }
    }
    mat_mat_mul(P_pred, HT, PHT, state_dim, state_dim, meas_dim);
    /* S = H * PHT */
    mat_mat_mul(H, PHT, S, meas_dim, state_dim, meas_dim);

    /* S += R */
    for (int i = 0; i < meas_dim * meas_dim; i++) {
        S[i] += R[i];
    }

    /* d^2 = nu' * S^{-1} * nu via Cholesky solve */
    double nu_copy[TRACK_MAX_MEAS_DIM];
    vec_copy(nu, nu_copy, meas_dim);

    double S_inv_nu[TRACK_MAX_MEAS_DIM];
    if (mat_solve_cholesky(S, nu_copy, S_inv_nu, meas_dim) != 0) {
        /* Fallback: use pseudo-inverse via direct solve attempt */
        return -1.0;
    }

    return vec_dot(nu, S_inv_nu, meas_dim);
}

int ellipsoidal_gate_check(double mahal_sq, double threshold)
{
    return (mahal_sq >= 0 && mahal_sq <= threshold) ? 1 : 0;
}

double chi2_threshold_approx(double P_G, int n_z)
{
    /* Wilson-Hilferty approximation for chi-squared quantile.
     *
     * γ ≈ n_z · (1 − 2/(9·n_z) + z_P · √(2/(9·n_z)))³
     *
     * where z_P = Φ⁻¹(P_G) is the standard normal quantile.
     */
    if (P_G <= 0.0 || P_G >= 1.0 || n_z <= 0) return 0.0;

    double z_P = normal_quantile(P_G);
    double a = 2.0 / (9.0 * n_z);
    double inner = 1.0 - a + z_P * sqrt(a);
    double gamma = n_z * inner * inner * inner;

    return (gamma > 0) ? gamma : 0.0;
}

double normal_quantile(double p)
{
    /* Abramowitz-Stegun approximation 26.2.23 for the standard normal quantile.
     *
     * For 0 < p < 0.5: use the rational approximation for the lower tail.
     * For p ≥ 0.5: use symmetry: z_p = −z_{1−p}.
     */
    if (p <= 0.0) return -10.0; /* practical limit */
    if (p >= 1.0) return 10.0;

    double p_use = (p < 0.5) ? p : 1.0 - p;
    double t = sqrt(-2.0 * log(p_use));

    /* Coefficients from AS 26.2.23 */
    const double c0 = 2.515517;
    const double c1 = 0.802853;
    const double c2 = 0.010328;
    const double d1 = 1.432788;
    const double d2 = 0.189269;
    const double d3 = 0.001308;

    double num = c0 + c1 * t + c2 * t * t;
    double den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;
    double z = t - num / den;

    return (p < 0.5) ? -z : z;
}

/* ============================================================================
 * Track scoring (LLR)
 * ============================================================================
 */

double track_score_increment(double P_D, double beta_F, const double *S,
                              double d2, int n_z)
{
    /* ΔLLR = log(P_D / (β_F · (2π)^(n_z/2) · √|S|)) − d²/2
     *
     * Simplification: Use the determinant of S for volume penalty.
     */
    if (n_z <= 0 || !S) return 0.0;

    double det_S = mat_det_cholesky(S, n_z);
    if (det_S <= 0.0) det_S = 1e-10; /* avoid log(0) */

    double vol_term = beta_F * pow(2.0 * M_PI, n_z / 2.0) * sqrt(det_S);

    double increment = log(P_D / vol_term) - 0.5 * d2;
    return increment;
}

double track_score_miss_penalty(double P_D)
{
    if (P_D >= 1.0) return -1e10; /* impossible to miss */
    if (P_D <= 0.0) return 0.0;    /* never detects */
    return log(1.0 - P_D);
}

void track_update_score(track_t *track, int associated, double P_D,
                         double beta_F, double d2)
{
    if (!track) return;
    if (associated) {
        track->score += track_score_increment(P_D, beta_F,
                                               track->S, d2, track->meas_dim);
        track->updates_total++;
    } else {
        track->score += track_score_miss_penalty(P_D);
    }
}

/* ============================================================================
 * Matrix and vector utilities
 * ============================================================================
 */

void mat_vec_mul(const double *A, const double *x, double *y, int m, int n)
{
    for (int i = 0; i < m; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i * n + j] * x[j];
        }
    }
}

void mat_mat_mul(const double *A, const double *B, double *C,
                  int m, int k, int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
            for (int p = 0; p < k; p++) {
                C[i * n + j] += A[i * k + p] * B[p * n + j];
            }
        }
    }
}

void mat_transpose(const double *A, double *B, int m, int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[j * m + i] = A[i * n + j];
        }
    }
}

double mat_inv_cholesky(const double *A, double *A_inv, int n)
{
    /* Cholesky decomposition: A = L * L', L lower-triangular.
     * Then A^{-1} = (L^{-1})' * L^{-1}.
     * Returns determinant of A, or 0.0 on failure.
     */
    if (!A || !A_inv || n <= 0) return 0.0;

    double L[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(L, 0, sizeof(double) * n * n);

    double det = 1.0;

    /* Cholesky: L computed in-place */
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += L[j * n + k] * L[j * n + k];
        }
        double diag = A[j * n + j] - sum;
        if (diag <= 0.0) return 0.0; /* not SPD */
        L[j * n + j] = sqrt(diag);
        det *= L[j * n + j] * L[j * n + j];

        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = (A[i * n + j] - sum) / L[j * n + j];
        }
    }

    /* Forward substitution to find L^{-1} column by column */
    double Linv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(Linv, 0, sizeof(double) * n * n);

    for (int col = 0; col < n; col++) {
        double e[TRACK_MAX_STATE_DIM];
        memset(e, 0, sizeof(double) * n);
        e[col] = 1.0;

        /* Solve L * x = e (forward substitution) */
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[i * n + k] * Linv[k * n + col];
            }
            Linv[i * n + col] = (e[i] - sum) / L[i * n + i];
        }
    }

    /* A^{-1} = (L^{-1})' * L^{-1} */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += Linv[k * n + i] * Linv[k * n + j];
            }
            A_inv[i * n + j] = sum;
        }
    }

    return det;
}

double mat_det_cholesky(const double *A, int n)
{
    if (!A || n <= 0) return 0.0;

    double L[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(L, 0, sizeof(double) * n * n);
    double det = 1.0;

    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += L[j * n + k] * L[j * n + k];
        }
        double diag = A[j * n + j] - sum;
        if (diag <= 0.0) return 0.0;
        L[j * n + j] = sqrt(diag);
        det *= diag;

        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = (A[i * n + j] - sum) / L[j * n + j];
        }
    }
    return det;
}

int mat_solve_cholesky(const double *A, const double *b, double *x, int n)
{
    /* Solve A*x = b for symmetric positive definite A via Cholesky.
     * 1. Compute L such that A = L*L'.
     * 2. Solve L*y = b (forward).
     * 3. Solve L'*x = y (backward).
     */
    if (!A || !b || !x || n <= 0) return -1;

    double L[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(L, 0, sizeof(double) * n * n);

    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += L[j * n + k] * L[j * n + k];
        }
        double diag = A[j * n + j] - sum;
        if (diag <= 0.0) return -1;
        L[j * n + j] = sqrt(diag);

        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = (A[i * n + j] - sum) / L[j * n + j];
        }
    }

    /* Forward: solve L*y = b */
    double y[TRACK_MAX_STATE_DIM];
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += L[i * n + k] * y[k];
        }
        y[i] = (b[i] - sum) / L[i * n + i];
    }

    /* Backward: solve L'*x = y */
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int k = i + 1; k < n; k++) {
            sum += L[k * n + i] * x[k];
        }
        x[i] = (y[i] - sum) / L[i * n + i];
    }

    return 0;
}

double vec_dot(const double *a, const double *b, int n)
{
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += a[i] * b[i];
    }
    return result;
}

double vec_norm2(const double *v, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

void vec_sub(const double *a, const double *b, double *c, int n)
{
    for (int i = 0; i < n; i++) {
        c[i] = a[i] - b[i];
    }
}

void vec_add(const double *a, const double *b, double *c, int n)
{
    for (int i = 0; i < n; i++) {
        c[i] = a[i] + b[i];
    }
}

void vec_scale(const double *x, double alpha, double *y, int n)
{
    for (int i = 0; i < n; i++) {
        y[i] = alpha * x[i];
    }
}

void vec_zero(double *v, int n)
{
    for (int i = 0; i < n; i++) {
        v[i] = 0.0;
    }
}

void vec_copy(const double *src, double *dst, int n)
{
    for (int i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}

double deg_to_rad(double deg) { return deg * M_PI / 180.0; }
double rad_to_deg(double rad) { return rad * 180.0 / M_PI; }

double wrap_angle(double angle)
{
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

double wrap_angle_2pi(double angle)
{
    while (angle >= 2.0 * M_PI) angle -= 2.0 * M_PI;
    while (angle < 0.0)         angle += 2.0 * M_PI;
    return angle;
}
