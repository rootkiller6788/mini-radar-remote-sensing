/**
 * track_metrics.c — Tracking Performance Metrics Implementation
 *
 * Implements: RMSE, AEE, NEES/ANEES, track purity/continuity/fragmentation,
 *             OSPA distance, CLEAR MOT metrics.
 *
 * References:
 *   - Bar-Shalom, Willett, Tian (2011), Ch. 11
 *   - Schuhmacher, Vo, Vo (2008) — OSPA
 *   - Bernardin & Stiefelhagen (2008) — CLEAR MOT
 */

#include "track_metrics.h"
#include "kalman_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Point estimation errors
 * ============================================================================
 */

double metrics_rmse_position(const double *estimates, const double *truths,
                              int K, int dim)
{
    if (!estimates || !truths || K <= 0 || dim <= 0) return 0.0;

    double sum_sq = 0.0;
    for (int k = 0; k < K; k++) {
        double dist_sq = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = estimates[k * dim + d] - truths[k * dim + d];
            dist_sq += diff * diff;
        }
        sum_sq += dist_sq;
    }
    return sqrt(sum_sq / K);
}

double metrics_rmse_velocity(const double *estimates, const double *truths,
                              int K, int dim)
{
    /* Same computation as position RMSE, applied to velocity components */
    return metrics_rmse_position(estimates, truths, K, dim);
}

double metrics_aee(const double *estimates, const double *truths,
                    int K, int dim)
{
    if (!estimates || !truths || K <= 0 || dim <= 0) return 0.0;

    double sum = 0.0;
    for (int k = 0; k < K; k++) {
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = estimates[k * dim + d] - truths[k * dim + d];
            dist += diff * diff;
        }
        sum += sqrt(dist);
    }
    return sum / K;
}

double metrics_nees(const double *estimate, const double *truth,
                     const double *P, int n)
{
    if (!estimate || !truth || !P || n <= 0) return 0.0;

    double dx[TRACK_MAX_STATE_DIM];
    vec_sub(estimate, truth, dx, n);

    double P_inv[TRACK_MAX_STATE_DIM * TRACK_MAX_STATE_DIM];
    if (mat_inv_cholesky(P, P_inv, n) <= 0.0) return -1.0;

    double nees = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            nees += dx[i] * P_inv[i * n + j] * dx[j];
        }
    }
    return nees;
}

double metrics_anees(const double *estimates, const double *truths,
                      const double *P_seq, int N, int K, int n)
{
    if (!estimates || !truths || !P_seq || N <= 0 || K <= 0 || n <= 0)
        return 0.0;

    double sum_nees = 0.0;
    int count = 0;

    for (int run = 0; run < N; run++) {
        for (int k = 0; k < K; k++) {
            const double *x = &estimates[(run * K + k) * n];
            /* Truth is same across runs */
            const double *x_t = &truths[k * n];

            /* Use per-step covariance or shared covariance */
            const double *P;
            if (P_seq) {
                /* Check if P_seq is per-step (K × n²) or single (n²) */
                /* Assume per-step if expanded */
                P = &P_seq[k * n * n];
            } else {
                P = NULL;
                return -1.0;
            }

            double nees = metrics_nees(x, x_t, P, n);
            if (nees >= 0.0) {
                sum_nees += nees;
                count++;
            }
        }
    }

    return (count > 0) ? sum_nees / count : 0.0;
}

void metrics_anees_interval(int n, int N, int K, double alpha,
                              double *lower, double *upper)
{
    double z = normal_quantile(1.0 - alpha / 2.0);
    double std = sqrt(2.0 * n / (N * K));
    *lower = n - z * std;
    *upper = n + z * std;
}

int metrics_is_consistent(double anees, int n, int N, int K, double alpha)
{
    double lower, upper;
    metrics_anees_interval(n, N, K, alpha, &lower, &upper);
    return (anees >= lower && anees <= upper) ? 1 : 0;
}

/* ============================================================================
 * Multi-target tracking metrics
 * ============================================================================
 */

void metrics_track_purity(const int *const *truth_assignments,
                           const int *track_lengths, int T,
                           double *purities)
{
    if (!truth_assignments || !track_lengths || !purities) return;

    for (int t = 0; t < T; t++) {
        int len = track_lengths[t];
        if (len <= 0) { purities[t] = 0.0; continue; }

        /* Count assignments to each truth ID */
        int *truth_counts = (int *)calloc(len + 1, sizeof(int));
        int max_count = 0;

        for (int k = 0; k < len; k++) {
            int tid = truth_assignments[t][k];
            if (tid >= 0 && tid < len) {
                truth_counts[tid]++;
                if (truth_counts[tid] > max_count) max_count = truth_counts[tid];
            }
        }

        purities[t] = (double)max_count / len;
        free(truth_counts);
    }
}

double metrics_average_purity(const double *purities, int T)
{
    if (!purities || T <= 0) return 0.0;
    double sum = 0.0;
    for (int t = 0; t < T; t++) sum += purities[t];
    return sum / T;
}

void metrics_track_continuity(const int *const *truth_assignments,
                               const int *track_lengths, int T,
                               int num_truths, double *continuity)
{
    if (!truth_assignments || !track_lengths || !continuity) return;

    /* For each truth, count how many tracks covered it */
    for (int tr = 0; tr < num_truths; tr++) {
        int track_count = 0;
        for (int t = 0; t < T; t++) {
            int covered = 0;
            for (int k = 0; k < track_lengths[t] && !covered; k++) {
                if (truth_assignments[t][k] == tr) covered = 1;
            }
            if (covered) track_count++;
        }
        continuity[tr] = (track_count > 0) ? 1.0 / track_count : 0.0;
    }
}

void metrics_track_fragmentation(const int *const *truth_assignments,
                                  const int *track_lengths, int T,
                                  int num_truths, int *frag_count)
{
    if (!truth_assignments || !track_lengths || !frag_count) return;

    for (int tr = 0; tr < num_truths; tr++) {
        int switches = 0;
        int prev_track = -1;

        for (int t = 0; t < T; t++) {
            for (int k = 0; k < track_lengths[t]; k++) {
                if (truth_assignments[t][k] == tr) {
                    if (prev_track >= 0 && prev_track != t) {
                        switches++;
                    }
                    prev_track = t;
                }
            }
        }
        frag_count[tr] = switches;
    }
}

double metrics_false_track_rate(const int *const *truth_assignments,
                                 const int *track_lengths, int T,
                                 double total_time)
{
    if (!truth_assignments || !track_lengths || T <= 0 || total_time <= 0.0)
        return 0.0;

    int false_tracks = 0;
    for (int t = 0; t < T; t++) {
        int has_truth = 0;
        for (int k = 0; k < track_lengths[t] && !has_truth; k++) {
            if (truth_assignments[t][k] >= 0) has_truth = 1;
        }
        if (!has_truth) false_tracks++;
    }
    return (double)false_tracks / total_time;
}

double metrics_initiation_delay(const double *truth_appear_time,
                                 const double *track_start_time,
                                 const int *track_to_truth, int T)
{
    if (!truth_appear_time || !track_start_time || !track_to_truth) return 0.0;

    double total_delay = 0.0;
    int count = 0;

    for (int t = 0; t < T; t++) {
        int tr = track_to_truth[t];
        if (tr >= 0 && track_start_time[t] > truth_appear_time[tr]) {
            total_delay += track_start_time[t] - truth_appear_time[tr];
            count++;
        }
    }
    return (count > 0) ? total_delay / count : 0.0;
}

double metrics_termination_delay(const double *truth_disappear_time,
                                  const double *track_end_time,
                                  const int *track_to_truth, int T)
{
    if (!truth_disappear_time || !track_end_time || !track_to_truth) return 0.0;

    double total_delay = 0.0;
    int count = 0;

    for (int t = 0; t < T; t++) {
        int tr = track_to_truth[t];
        if (tr >= 0 && track_end_time[t] > truth_disappear_time[tr]) {
            total_delay += track_end_time[t] - truth_disappear_time[tr];
            count++;
        }
    }
    return (count > 0) ? total_delay / count : 0.0;
}

/* ============================================================================
 * OSPA distance
 * ============================================================================
 */

double metrics_ospa(const double *X_est, int n_est,
                     const double *X_true, int n_true,
                     int dim, double c, int p)
{
    if (!X_est || !X_true) return c; /* cardinality penalty only */
    if (n_est == 0 && n_true == 0) return 0.0;

    int m = (n_est > n_true) ? n_est : n_true;
    if (m == 0) return 0.0;

    /* Build distance matrix with cutoff c */
    double *D = (double *)calloc(m * m, sizeof(double));
    if (!D) return c;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i < n_est && j < n_true) {
                double dist = 0.0;
                for (int d = 0; d < dim; d++) {
                    double diff = X_est[i * dim + d] - X_true[j * dim + d];
                    dist += diff * diff;
                }
                dist = sqrt(dist);
                D[i * m + j] = (dist < c) ? dist : c;
            } else {
                D[i * m + j] = c; /* dummy assignment cost */
            }
        }
    }

    /* Hungarian algorithm for optimal assignment */
    int *assignment = (int *)malloc(m * sizeof(int));
    double *u = (double *)calloc(m, sizeof(double));
    double *v = (double *)calloc(m, sizeof(double));
    int *p_row = (int *)calloc(m, sizeof(int));
    int *way = (int *)calloc(m, sizeof(int));

    for (int i = 1; i <= m; i++) {
        p_row[0] = i;
        int j0 = 0;
        double *minv = (double *)calloc(m + 1, sizeof(double));
        int *used = (int *)calloc(m + 1, sizeof(int));
        for (int jj = 0; jj <= m; jj++) minv[jj] = 1e18;

        do {
            used[j0] = 1;
            int i0 = p_row[j0];
            double delta = 1e18;
            int j1 = 0;
            for (int j = 1; j <= m; j++) {
                if (!used[j]) {
                    double cur = D[(i0 - 1) * m + (j - 1)] - u[i0] - v[j];
                    if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                    if (minv[j] < delta) { delta = minv[j]; j1 = j; }
                }
            }
            for (int j = 0; j <= m; j++) {
                if (used[j]) { u[p_row[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p_row[j0] != 0);

        do {
            int j1 = way[j0];
            p_row[j0] = p_row[j1];
            j0 = j1;
        } while (j0 != 0);
        free(minv);
        free(used);
    }

    for (int j = 1; j <= m; j++) {
        assignment[p_row[j] - 1] = (p_row[j] > 0 && p_row[j] <= m) ? j - 1 : -1;
    }

    /* Compute OSPA */
    double sum_p = 0.0;
    for (int i = 0; i < n_est; i++) {
        int j = (i < m) ? assignment[i] : -1;
        double cost;
        if (j >= 0 && j < n_true) {
            double dist = 0.0;
            for (int d = 0; d < dim; d++) {
                double diff = X_est[i * dim + d] - X_true[j * dim + d];
                dist += diff * diff;
            }
            dist = sqrt(dist);
            cost = (dist < c) ? dist : c;
        } else {
            cost = c;
        }
        sum_p += pow(cost, p);
    }

    /* Cardinality penalty */
    sum_p += pow(c, p) * abs(n_true - n_est);

    double ospa = pow(sum_p / m, 1.0 / p);

    free(D); free(assignment); free(u); free(v); free(p_row);
    return ospa;
}

double metrics_ospa_avg(const double *X_est_seq, const int *n_est_seq,
                         const double *X_true_seq, const int *n_true_seq,
                         int K, int dim, double c, int p)
{
    if (!X_est_seq || !n_est_seq || !X_true_seq || !n_true_seq || K <= 0)
        return 0.0;

    double sum = 0.0;
    int est_offset = 0, true_offset = 0;

    for (int k = 0; k < K; k++) {
        int ne = n_est_seq[k];
        int nt = n_true_seq[k];

        sum += metrics_ospa(X_est_seq + est_offset, ne,
                             X_true_seq + true_offset, nt, dim, c, p);

        est_offset += ne * dim;
        true_offset += nt * dim;
    }

    return sum / K;
}

/* ============================================================================
 * CLEAR MOT metrics
 * ============================================================================
 */

void metrics_clear_mot(const int *const *assignments,
                        const int *track_lengths,
                        const double *const *distances,
                        int K, int T,
                        double *mota, double *motp, int *id_switches)
{
    (void)K; /* frame count available for per-frame analysis extension */
    if (!assignments || !track_lengths || !distances) return;

    int total_FN = 0, total_FP = 0, total_IDSW = 0;
    int total_GT = 0;
    double total_dist = 0.0;
    int total_matches = 0;

    /* Simplified CLEAR MOT computation.
     * Full computation requires frame-by-frame matching with ground truth.
     */
    for (int t = 0; t < T; t++) {
        int len = track_lengths[t];
        total_GT += len;

        for (int k = 0; k < len; k++) {
            if (assignments[t][k] < 0) {
                total_FN++; /* no truth for this track point */
            } else {
                total_matches++;
                if (distances && distances[t]) {
                    total_dist += distances[t][k];
                }
            }
        }
    }

    /* MOTA = 1 - (FN + FP + IDSW) / GT */
    if (total_GT > 0) {
        *mota = 1.0 - (double)(total_FN + total_FP + total_IDSW) / total_GT;
        if (*mota < 0.0) *mota = 0.0;
    } else {
        *mota = 0.0;
    }

    /* MOTP = total_dist / total_matches */
    *motp = (total_matches > 0) ? total_dist / total_matches : 0.0;

    if (id_switches) *id_switches = total_IDSW;
}
