/**
 * data_association.c — Data Association Algorithms Implementation
 *
 * Implements: GNN via Hungarian algorithm, NN, PDA, JPDA, MHT, Auction.
 *
 * References:
 *   - Kuhn, H.W. "The Hungarian Method" (1955)
 *   - Munkres, J. "Algorithms for Assignment and Transportation Problems" (1957)
 *   - Bar-Shalom, Tse "Tracking in Clutter" Automatica (1975)
 *   - Fortmann, Bar-Shalom, Scheffe "JPDA" IEEE JOE (1983)
 *   - Reid, D.B. "MHT" IEEE TAC (1979)
 *   - Bertsekas, D.P. "Auction Algorithm" (1988)
 */

#include "data_association.h"
#include "kalman_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Cost matrix construction
 * ============================================================================
 */

void association_build_cost_matrix(const track_t *tracks, int num_tracks,
                                    const measurement_t *measurements, int num_meas,
                                    cost_matrix_t *matrix)
{
    if (!tracks || !measurements || !matrix) return;

    matrix->num_tracks = num_tracks;
    matrix->num_measurements = num_meas;

    /* Allocate cost and gate mask */
    int size = num_tracks * num_meas;
    matrix->cost = (double *)realloc(matrix->cost, size * sizeof(double));
    matrix->gate_mask = (int *)realloc(matrix->gate_mask, size * sizeof(int));
    if (!matrix->cost || !matrix->gate_mask) return;

    for (int i = 0; i < num_tracks; i++) {
        const track_t *t = &tracks[i];
        if (t->state == TRACK_STATE_FREE || t->state == TRACK_STATE_DELETED) {
            /* Inactive tracks: infinite cost */
            for (int j = 0; j < num_meas; j++) {
                matrix->cost[i * num_meas + j] = 1e15;
                matrix->gate_mask[i * num_meas + j] = 0;
            }
            continue;
        }

        for (int j = 0; j < num_meas; j++) {
            const measurement_t *m = &measurements[j];

            /* Compute Mahalanobis distance using track's prediction */
            /* Build measurement matrix for this track */
            double H_cart[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
            memset(H_cart, 0, sizeof(H_cart));
            int st_dim = t->state_dim;
            int ms_dim = t->meas_dim;

            /* Simple position-measurement H: measure first ms_dim states */
            for (int d = 0; d < ms_dim && d < st_dim; d++) {
                H_cart[d * st_dim + d] = 1.0;
            }

            double d2 = mahalanobis_distance_sq(m->z, m->R,
                                                  t->x_pred, t->P_pred,
                                                  H_cart, st_dim, ms_dim);
            if (d2 < 0) d2 = 1e15;

            matrix->cost[i * num_meas + j] = d2;

            /* Gate check */
            double thr = matrix->gate_threshold > 0
                         ? matrix->gate_threshold
                         : chi2_threshold_approx(0.99, ms_dim);
            matrix->gate_mask[i * num_meas + j] = (d2 <= thr) ? 1 : 0;

            /* Infinite cost for out-of-gate measurements */
            if (!matrix->gate_mask[i * num_meas + j]) {
                matrix->cost[i * num_meas + j] = 1e15;
            }
        }
    }
}

void association_build_gate_mask(cost_matrix_t *matrix, double gate_threshold)
{
    if (!matrix) return;
    matrix->gate_threshold = gate_threshold;
    int nt = matrix->num_tracks;
    int nm = matrix->num_measurements;
    for (int i = 0; i < nt; i++) {
        for (int j = 0; j < nm; j++) {
            matrix->gate_mask[i * nm + j] =
                (matrix->cost[i * nm + j] <= gate_threshold) ? 1 : 0;
        }
    }
}

/* ============================================================================
 * Nearest Neighbor (NN)
 * ============================================================================
 */

assignment_result_t *nn_associate(const cost_matrix_t *matrix,
                                   const double *track_scores)
{
    if (!matrix) return NULL;

    /* Sort track indices by score (descending: best tracks first) */
    int nt = matrix->num_tracks;
    int nm = matrix->num_measurements;

    int *track_order = (int *)malloc(nt * sizeof(int));
    double *track_keys = (double *)malloc(nt * sizeof(double));
    if (!track_order || !track_keys) {
        free(track_order); free(track_keys); return NULL;
    }

    for (int i = 0; i < nt; i++) {
        track_order[i] = i;
        track_keys[i] = (track_scores) ? track_scores[i] : 0.0;
    }

    /* Simple insertion sort by descending score */
    for (int i = 1; i < nt; i++) {
        int key_i = track_order[i];
        double key_s = track_keys[i];
        int j = i - 1;
        while (j >= 0 && track_keys[j] < key_s) {
            track_order[j + 1] = track_order[j];
            track_keys[j + 1] = track_keys[j];
            j--;
        }
        track_order[j + 1] = key_i;
        track_keys[j + 1] = key_s;
    }
    free(track_keys);

    assignment_result_t *result = assignment_alloc(nt, nm);
    if (!result) { free(track_order); return NULL; }

    int *meas_taken = (int *)calloc(nm, sizeof(int));
    if (!meas_taken) { free(track_order); assignment_free(result); return NULL; }

    for (int idx = 0; idx < nt; idx++) {
        int i = track_order[idx];
        double min_cost = 1e15;
        int best_j = -1;

        for (int j = 0; j < nm; j++) {
            if (!meas_taken[j] && matrix->cost[i * nm + j] < min_cost) {
                min_cost = matrix->cost[i * nm + j];
                best_j = j;
            }
        }

        if (best_j >= 0 && min_cost < 1e14) {
            result->assignment[i] = best_j;
            result->total_cost += min_cost;
            meas_taken[best_j] = 1;
        }
    }

    /* Collect unassigned measurements */
    result->num_unassigned = 0;
    for (int j = 0; j < nm; j++) {
        if (!meas_taken[j]) {
            result->unassigned_meas[result->num_unassigned++] = j;
        }
    }

    free(track_order);
    free(meas_taken);
    return result;
}

/* ============================================================================
 * Hungarian algorithm for optimal assignment (GNN)
 * ============================================================================
 */

/* Internal Hungarian state */
#define HUNG_MAX 256
static double hung_cost[HUNG_MAX * HUNG_MAX];
static double hung_u[HUNG_MAX], hung_v[HUNG_MAX];
static int hung_p[HUNG_MAX], hung_way[HUNG_MAX];
static int hung_n, hung_m;

assignment_result_t *gnn_hungarian_associate(const cost_matrix_t *matrix)
{
    if (!matrix) return NULL;
    int n = matrix->num_tracks;
    int m = matrix->num_measurements;
    if (n > HUNG_MAX || m > HUNG_MAX) {
        /* Fall back to NN for large problems */
        return nn_associate(matrix, NULL);
    }

    hung_n = n;
    hung_m = m;
    int N = (n > m) ? n : m; /* square matrix size */

    /* Copy cost matrix, padding to square with zeros */
    memset(hung_cost, 0, sizeof(hung_cost));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            hung_cost[i * N + j] = matrix->cost[i * m + j];
        }
    }

    /* Hungarian algorithm (O(N³)) */
    memset(hung_u, 0, N * sizeof(double));
    memset(hung_v, 0, N * sizeof(double));
    memset(hung_p, 0, N * sizeof(int));
    memset(hung_way, 0, N * sizeof(int));

    for (int i = 1; i <= N; i++) {
        hung_p[0] = i;
        int j0 = 0;
        double minv[HUNG_MAX];
        for (int k = 0; k <= N; k++) minv[k] = 1e18;
        int used[HUNG_MAX] = {0};

        do {
            used[j0] = 1;
            int i0 = hung_p[j0];
            double delta = 1e18;
            int j1 = 0;

            for (int j = 1; j <= N; j++) {
                if (!used[j]) {
                    double cur = hung_cost[(i0 - 1) * N + (j - 1)]
                                 - hung_u[i0] - hung_v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        hung_way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }

            for (int j = 0; j <= N; j++) {
                if (used[j]) {
                    hung_u[hung_p[j]] += delta;
                    hung_v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (hung_p[j0] != 0);

        do {
            int j1 = hung_way[j0];
            hung_p[j0] = hung_p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    /* Extract assignment */
    assignment_result_t *result = assignment_alloc(n, m);
    if (!result) return NULL;

    int *assigned = (int *)calloc(N + 1, sizeof(int));
    for (int j = 1; j <= N; j++) {
        if (hung_p[j] > 0) {
            assigned[hung_p[j]] = j;
        }
    }

    for (int i = 1; i <= n; i++) {
        int j = assigned[i] - 1;
        if (j >= 0 && j < m && matrix->cost[(i - 1) * m + j] < 1e14) {
            result->assignment[i - 1] = j;
            result->total_cost += matrix->cost[(i - 1) * m + j];
        }
    }

    /* Unassigned measurements */
    int *meas_assigned = (int *)calloc(m, sizeof(int));
    for (int i = 0; i < n; i++) {
        if (result->assignment[i] >= 0) {
            meas_assigned[result->assignment[i]] = 1;
        }
    }
    result->num_unassigned = 0;
    for (int j = 0; j < m; j++) {
        if (!meas_assigned[j]) {
            result->unassigned_meas[result->num_unassigned++] = j;
        }
    }

    free(assigned);
    free(meas_assigned);
    return result;
}

assignment_result_t *gnn_gated_associate(const cost_matrix_t *matrix)
{
    return gnn_hungarian_associate(matrix);
}

/* ============================================================================
 * Probabilistic Data Association (PDA)
 * ============================================================================
 */

int pda_update(track_t *track, const measurement_t *measurements,
               int num_validated, const int *valid_indices,
               double P_D, double beta_F,
               double *x_updated, double *P_updated)
{
    if (!track || !measurements || !valid_indices || !x_updated || !P_updated)
        return -1;
    if (num_validated < 1) return -1;

    int n = track->state_dim;
    int m = track->meas_dim;

    /* Compute association weight for each validated measurement */
    double *e = (double *)calloc(num_validated + 1, sizeof(double));
    if (!e) return -1;

    double det_S = mat_det_cholesky(track->S, m);
    if (det_S <= 0.0) det_S = 1e-10;
    double vol_S = sqrt(det_S);

    double sum_e = 0.0;
    for (int k = 0; k < num_validated; k++) {
        int idx = valid_indices[k];
        const measurement_t *meas = &measurements[idx];

        /* Compute Mahalanobis distance for this measurement */
        double H_c[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
        memset(H_c, 0, sizeof(H_c));
        for (int d = 0; d < m && d < n; d++) H_c[d * n + d] = 1.0;

        double nu_k[TRACK_MAX_MEAS_DIM];
        double Hx[TRACK_MAX_MEAS_DIM];
        mat_vec_mul(H_c, track->x_pred, Hx, m, n);
        vec_sub(meas->z, Hx, nu_k, m);

        double S_inv_nu[TRACK_MAX_MEAS_DIM];
        vec_copy(nu_k, S_inv_nu, m);
        if (mat_solve_cholesky(track->S, nu_k, S_inv_nu, m) != 0) continue;

        double d2 = vec_dot(nu_k, S_inv_nu, m);

        /* Association weight numerator */
        e[k] = exp(-0.5 * d2);
        sum_e += e[k];
    }

    /* Weight for "no detection" hypothesis */
    double b = beta_F * pow(2.0 * M_PI, m / 2.0) * vol_S * (1.0 - P_D) / P_D;
    e[num_validated] = b;
    sum_e += b;

    /* Normalize to get β_k */
    double *beta = (double *)malloc((num_validated + 1) * sizeof(double));
    if (!beta) { free(e); return -1; }

    if (sum_e > 0.0) {
        for (int k = 0; k <= num_validated; k++) {
            beta[k] = e[k] / sum_e;
        }
    } else {
        beta[num_validated] = 1.0;
        for (int k = 0; k < num_validated; k++) beta[k] = 0.0;
    }

    /* Combined innovation */
    double nu_combined[TRACK_MAX_MEAS_DIM];
    memset(nu_combined, 0, m * sizeof(double));
    for (int k = 0; k < num_validated; k++) {
        int idx = valid_indices[k];
        double H_c[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
        memset(H_c, 0, sizeof(H_c));
        for (int d = 0; d < m && d < n; d++) H_c[d * n + d] = 1.0;

        double Hx[TRACK_MAX_MEAS_DIM];
        mat_vec_mul(H_c, track->x_pred, Hx, m, n);
        double nu_k[TRACK_MAX_MEAS_DIM];
        vec_sub(measurements[idx].z, Hx, nu_k, m);

        for (int d = 0; d < m; d++) {
            nu_combined[d] += beta[k] * nu_k[d];
        }
    }

    /* PDA Kalman gain: K = P_pred*H'*S^{-1} */
    double H_simple[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
    memset(H_simple, 0, sizeof(H_simple));
    for (int d = 0; d < m && d < n; d++) H_simple[d * n + d] = 1.0;

    double HT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_transpose(H_simple, HT, m, n);
    double PHT[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_mat_mul(track->P_pred, HT, PHT, n, n, m);
    double S_inv[TRACK_MAX_MEAS_DIM * TRACK_MAX_MEAS_DIM];
    mat_inv_cholesky(track->S, S_inv, m);
    double K_pda[TRACK_MAX_STATE_DIM * TRACK_MAX_MEAS_DIM];
    mat_mat_mul(PHT, S_inv, K_pda, n, m, m);

    /* State update */
    double Knu[TRACK_MAX_STATE_DIM];
    mat_vec_mul(K_pda, nu_combined, Knu, n, m);
    vec_add(track->x_pred, Knu, x_updated, n);

    /* Covariance update (PDA form) */
    double I_KH[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(I_KH, 0, sizeof(I_KH));
    double KH[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(K_pda, H_simple, KH, n, m, n);
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            I_KH[a * n + b] = -KH[a * n + b];
            if (a == b) I_KH[a * n + b] += 1.0;
        }
    }

    double P_corrected[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    mat_mat_mul(I_KH, track->P_pred, P_corrected, n, n, n);

    /* Spread of innovations term P_tilde */
    double P_tilde[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    memset(P_tilde, 0, sizeof(P_tilde));
    for (int k = 0; k < num_validated; k++) {
        int idx = valid_indices[k];
        double nu_k[TRACK_MAX_MEAS_DIM];
        double H_c[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
        memset(H_c, 0, sizeof(H_c));
        for (int d = 0; d < m && d < n; d++) H_c[d * n + d] = 1.0;
        double Hx[TRACK_MAX_MEAS_DIM];
        mat_vec_mul(H_c, track->x_pred, Hx, m, n);
        vec_sub(measurements[idx].z, Hx, nu_k, m);

        /* u_k = K * (nu_k - nu_combined) */
        double diff[TRACK_MAX_MEAS_DIM];
        vec_sub(nu_k, nu_combined, diff, m);
        double u_k[TRACK_MAX_STATE_DIM];
        mat_vec_mul(K_pda, diff, u_k, n, m);

        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                P_tilde[a * n + b] += beta[k] * u_k[a] * u_k[b];
            }
        }
    }
    /* Subtract self term */
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            P_tilde[a * n + b] -= Knu[a] * Knu[b];
        }
    }

    /* P = beta0 * P_pred + (1-beta0) * P_corrected + P_tilde */
    double beta0 = beta[num_validated];
    for (int a = 0; a < n * n; a++) {
        P_updated[a] = beta0 * track->P_pred[a]
                       + (1.0 - beta0) * P_corrected[a]
                       + P_tilde[a];
    }

    free(e);
    free(beta);
    return 0;
}

/* ============================================================================
 * JPDA
 * ============================================================================
 */

int jpda_compute_beta(track_t *tracks, int num_tracks,
                      const measurement_t *measurements, int num_meas,
                      double P_D, double beta_F,
                      double **beta_matrix, double *beta0)
{
    if (!tracks || !measurements || !beta_matrix || !beta0) return -1;
    if (num_tracks <= 0 || num_meas <= 0) return -1;

    /* Simplified JPDA: for each track t, compute β_tj by approximate method.
     *
     * βₜⱼ = eₜⱼ / (bₜ + Σ_{k∈gₜ} eₜₖ)
     *
     * This is a simplified version that avoids enumeration of all feasible
     * joint events. Full JPDA enumeration is O(2^{NT·NM}) and impractical
     * without gating and efficient search (Murty's algorithm).
     */

    for (int t = 0; t < num_tracks; t++) {
        if (!beta_matrix[t]) {
            beta_matrix[t] = (double *)calloc(num_meas + 1, sizeof(double));
            if (!beta_matrix[t]) return -1;
        }

        double e_sum = 0.0;
        int m_dim = tracks[t].meas_dim;

        double det_S = mat_det_cholesky(tracks[t].S, m_dim);
        if (det_S <= 0.0) det_S = 1e-10;
        double b_t = beta_F * pow(2.0 * M_PI, m_dim / 2.0)
                     * sqrt(det_S) * (1.0 - P_D) / P_D;

        for (int j = 0; j < num_meas; j++) {
            double H_c[TRACK_MAX_MEAS_DIM * TRACK_MAX_STATE_DIM];
            memset(H_c, 0, sizeof(H_c));
            for (int d = 0; d < m_dim && d < tracks[t].state_dim; d++)
                H_c[d * tracks[t].state_dim + d] = 1.0;

            double d2 = mahalanobis_distance_sq(measurements[j].z,
                                                  measurements[j].R,
                                                  tracks[t].x_pred,
                                                  tracks[t].P_pred,
                                                  H_c, tracks[t].state_dim,
                                                  m_dim);
            if (d2 < 0) d2 = 1e15;

            double gate_thr = chi2_threshold_approx(0.99, m_dim);
            if (d2 <= gate_thr) {
                beta_matrix[t][j + 1] = exp(-0.5 * d2);
                e_sum += beta_matrix[t][j + 1];
            } else {
                beta_matrix[t][j + 1] = 0.0;
            }
        }

        beta_matrix[t][0] = b_t;
        e_sum += b_t;

        if (e_sum > 0.0) {
            for (int j = 0; j <= num_meas; j++) {
                beta_matrix[t][j] /= e_sum;
            }
        } else {
            beta_matrix[t][0] = 1.0;
        }
        beta0[t] = beta_matrix[t][0];
    }

    return 0;
}

int jpda_update_tracks(track_t *tracks, int num_tracks,
                       const measurement_t *measurements, int num_meas,
                       const double *const *beta_matrix, const double *beta0)
{
    if (!tracks || !beta_matrix || !beta0) return -1;

    for (int t = 0; t < num_tracks; t++) {
        /* Use PDA update for each track with JPDA β probabilities */
        int in_gate[2048];
        int n_in_gate = 0;
        for (int j = 0; j < num_meas; j++) {
            if (beta_matrix[t] && beta_matrix[t][j + 1] > 1e-10) {
                in_gate[n_in_gate++] = j;
            }
        }
        if (n_in_gate > 0) {
            pda_update(&tracks[t], measurements, n_in_gate, in_gate,
                       0.9, 1e-6, tracks[t].x, tracks[t].P);
        }
    }
    return 0;
}

/* ============================================================================
 * MHT
 * ============================================================================
 */

void mht_init(mht_tree_t *tree, int max_hypotheses, int max_depth)
{
    if (!tree) return;
    memset(tree, 0, sizeof(mht_tree_t));
    tree->max_hypotheses = max_hypotheses;
    tree->max_depth = max_depth;
    tree->hypotheses = (mht_hypothesis_t *)calloc(max_hypotheses, sizeof(mht_hypothesis_t));
    tree->num_hypotheses = 0;
    tree->best_index = -1;
}

int mht_generate_hypotheses(mht_tree_t *tree,
                             const cost_matrix_t *matrix,
                             const track_t *tracks, int num_tracks,
                             double P_D, double beta_F)
{
    (void)tracks;
    (void)P_D;
    (void)beta_F;
    if (!tree || !matrix) return 0;

    /* Generate new hypotheses from the best N existing ones */
    int n_new = 0;
    int nt = num_tracks;

    /* Simplified: create one hypothesis with GNN assignment */
    assignment_result_t *assign = gnn_hungarian_associate(matrix);
    if (!assign) return 0;

    if (tree->num_hypotheses < tree->max_hypotheses) {
        mht_hypothesis_t *h = &tree->hypotheses[tree->num_hypotheses];
        h->scan_id = tree->num_hypotheses;
        h->assignment = (int *)malloc(nt * sizeof(int));
        if (h->assignment) {
            memcpy(h->assignment, assign->assignment, nt * sizeof(int));
            h->num_tracks = nt;
            h->hypothesis_score = 0.0; /* computed externally */
            h->parent_hypothesis = tree->best_index;
            tree->best_index = tree->num_hypotheses;
            tree->num_hypotheses++;
            n_new++;
        }
    }

    assignment_free(assign);
    return n_new;
}

int mht_prune(mht_tree_t *tree, int keep_top_k)
{
    if (!tree) return 0;
    if (tree->num_hypotheses <= keep_top_k) return 0;
    int pruned = tree->num_hypotheses - keep_top_k;
    tree->num_hypotheses = keep_top_k;
    if (tree->best_index >= keep_top_k) tree->best_index = 0;
    return pruned;
}

const mht_hypothesis_t *mht_get_best(const mht_tree_t *tree)
{
    if (!tree || tree->best_index < 0) return NULL;
    return &tree->hypotheses[tree->best_index];
}

void mht_free(mht_tree_t *tree)
{
    if (!tree) return;
    for (int i = 0; i < tree->num_hypotheses; i++) {
        free(tree->hypotheses[i].assignment);
    }
    free(tree->hypotheses);
    memset(tree, 0, sizeof(mht_tree_t));
}

/* ============================================================================
 * Auction algorithm
 * ============================================================================
 */

int auction_assign(const double *cost, int n_tracks, int n_meas,
                   int *assignment, double epsilon)
{
    if (!cost || !assignment || n_tracks <= 0 || n_meas <= 0) return -1;
    if (epsilon <= 0.0) epsilon = 0.01;

    /* Auction: each track "bids" on a measurement.
     * Profit = -cost (maximizing profit = minimizing cost).
     *
     * For proper auction with track count ≠ meas count, add dummy
     * measurements with zero profit.
     */
    int N = (n_tracks > n_meas) ? n_tracks : n_meas;
    double *profit = (double *)calloc(N * N, sizeof(double));
    double *price = (double *)calloc(N, sizeof(double));
    int *assigned_meas = (int *)malloc(N * sizeof(int));

    if (!profit || !price || !assigned_meas) {
        free(profit); free(price); free(assigned_meas); return -1;
    }

    /* Build profit matrix (convert cost to profit) */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i < n_tracks && j < n_meas) {
                double c = cost[i * n_meas + j];
                if (c >= 1e14) c = 1e6; /* cap infinite cost */
                profit[i * N + j] = -c; /* profit = -cost */
            } else {
                profit[i * N + j] = -1e6; /* dummy: very low profit */
            }
        }
    }

    /* Initialize: no assignments, zero prices */
    memset(price, 0, N * sizeof(double));
    for (int i = 0; i < N; i++) assigned_meas[i] = -1;

    /* Auction iterations */
    for (int iter = 0; iter < N * 10; iter++) {
        int any_change = 0;

        for (int i = 0; i < n_tracks; i++) {
            /* Find best and second-best measurement for track i */
            double best_val = -1e18, second_best = -1e18;
            int best_j = -1;

            for (int j = 0; j < N; j++) {
                double val = profit[i * N + j] - price[j];
                if (val > best_val) {
                    second_best = best_val;
                    best_val = val;
                    best_j = j;
                } else if (val > second_best) {
                    second_best = val;
                }
            }

            if (best_j >= 0) {
                /* Bid: increase price by (best_val - second_best + epsilon) */
                double bid = best_val - second_best + epsilon;
                price[best_j] += bid;
                any_change = 1;

                /* Reassign: if best_j was assigned to someone else, free them */
                for (int k = 0; k < N; k++) {
                    if (assigned_meas[k] == best_j) {
                        assigned_meas[k] = -1;
                    }
                }
                assigned_meas[i] = best_j;
            }
        }

        if (!any_change) break;
    }

    /* Extract assignment */
    for (int i = 0; i < n_tracks; i++) {
        int j = assigned_meas[i];
        if (j >= 0 && j < n_meas) {
            assignment[i] = j;
        } else {
            assignment[i] = -1;
        }
    }

    free(profit); free(price); free(assigned_meas);
    return 0;
}

/* ============================================================================
 * Assignment utilities
 * ============================================================================
 */

assignment_result_t *assignment_alloc(int num_tracks, int num_measurements)
{
    assignment_result_t *r = (assignment_result_t *)calloc(1, sizeof(assignment_result_t));
    if (!r) return NULL;
    r->num_tracks = num_tracks;
    r->num_measurements = num_measurements;
    r->assignment = (int *)malloc(num_tracks * sizeof(int));
    r->unassigned_meas = (int *)malloc(num_measurements * sizeof(int));
    if (!r->assignment || !r->unassigned_meas) {
        assignment_free(r); return NULL;
    }
    for (int i = 0; i < num_tracks; i++) r->assignment[i] = -1;
    r->num_unassigned = 0;
    r->total_cost = 0.0;
    return r;
}

void assignment_free(assignment_result_t *result)
{
    if (!result) return;
    free(result->assignment);
    free(result->unassigned_meas);
    free(result);
}

int assignment_count_pairs(const assignment_result_t *result)
{
    if (!result) return 0;
    int count = 0;
    for (int i = 0; i < result->num_tracks; i++) {
        if (result->assignment[i] >= 0) count++;
    }
    return count;
}

void assignment_print(const assignment_result_t *result)
{
    if (!result) return;
    printf("Assignment: %d tracks, %d measurements, %d pairs\n",
           result->num_tracks, result->num_measurements,
           assignment_count_pairs(result));
    printf("Total cost: %.4f\n", result->total_cost);
    for (int i = 0; i < result->num_tracks; i++) {
        if (result->assignment[i] >= 0) {
            printf("  Track %d -> Meas %d\n", i, result->assignment[i]);
        }
    }
    if (result->num_unassigned > 0) {
        printf("  Unassigned measurements:");
        for (int j = 0; j < result->num_unassigned; j++) {
            printf(" %d", result->unassigned_meas[j]);
        }
        printf("\n");
    }
}
