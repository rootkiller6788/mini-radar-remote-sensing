/**
 * track_metrics.h — Tracking Performance Evaluation Metrics
 *
 * Covers: L1 (Definitions) — RMSE, NEES, ANEES, track purity, continuity
 *         L2 (Core Concepts) — Track quality assessment, consistency analysis
 *         L4 (Fundamental Laws) — Chi-squared consistency, CRLB comparison
 *         L5 (Algorithms/Methods) — Track scoring, OSPA distance, MOT metrics
 *
 * References:
 *   - Bar-Shalom, Willett, Tian "Tracking and Data Fusion" (2011), Ch. 11
 *   - Schuhmacher, D., Vo, B.T., Vo, B.N. "A Consistent Metric for Performance
 *     Evaluation of Multi-Object Filters" IEEE TSP (2008)
 *   - Ristic, B., Vo, B.N., Clark, D., Vo, B.T. "A Metric for Performance
 *     Evaluation of Multi-Target Tracking Algorithms" IEEE TSP (2011)
 *   - Bernardin, K., Stiefelhagen, R. "Evaluating Multiple Object Tracking
 *     Performance: The CLEAR MOT Metrics" (2008)
 *
 * Curriculum:
 *   - Stanford EE102A, Berkeley EE123, ETH 227-0427
 *   - Georgia Tech ECE 4270, Michigan EECS 351
 */

#ifndef TRACK_METRICS_H
#define TRACK_METRICS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"

/* ============================================================================
 * L1 — Point estimation error metrics
 * ============================================================================
 */

/**
 * Root Mean Square Error in position.
 *
 * RMSE_pos = √[(1/K) Σₖ ||x̂ₖ_pos − x_trueₖ_pos||²]
 *
 * Complexity: O(K·d) for K time steps and d spatial dimensions.
 *
 * @param estimates Estimated position sequence [K × dim]
 * @param truths    True position sequence [K × dim]
 * @param K         Number of time steps
 * @param dim       Spatial dimensionality (2 or 3)
 * @return RMSE in meters
 */
double metrics_rmse_position(const double *estimates, const double *truths,
                              int K, int dim);

/**
 * Root Mean Square Error in velocity.
 * Same formula as RMSE_pos applied to velocity.
 * Complexity: O(K·d)
 */
double metrics_rmse_velocity(const double *estimates, const double *truths,
                              int K, int dim);

/**
 * Average Euclidean Error (AEE) — L1-like metric.
 *
 * AEE = (1/K) Σₖ ||x̂ₖ_pos − x_trueₖ_pos||
 *
 * Complexity: O(K·d)
 */
double metrics_aee(const double *estimates, const double *truths,
                    int K, int dim);

/**
 * Normalized Estimation Error Squared (NEES) at a single time step.
 *
 * ε_k = (x̂_k − x_true_k)ᵀ·P_k⁻¹·(x̂_k − x_true_k)
 *
 * Under consistent (optimistic) filtering, ε_k ~ χ²(n) with E[ε_k] = n.
 *
 * NEES << n: filter is overconfident (covariance too small)
 * NEES >> n: filter is inconsistent (covariance too large or biased)
 *
 * Complexity: O(n³) for matrix inverse
 *
 * @param estimate State estimate [n]
 * @param truth    True state [n]
 * @param P        Filter covariance [n²]
 * @param n        State dimension
 * @return NEES value
 */
double metrics_nees(const double *estimate, const double *truth,
                     const double *P, int n);

/**
 * Average NEES over multiple time steps and Monte Carlo runs.
 *
 * ANEES = (1/(N·K)) Σᵢ Σₖ ε_{i,k}
 *
 * Under consistent filtering, ANEES ≈ n.
 * Acceptable interval (95%): [n − 2√(2n/(N·K)), n + 2√(2n/(N·K))]
 *
 * Complexity: O(N·K·n³)
 *
 * @param estimates Estimate sequence [N × K × n] (N runs, K steps, n dim)
 * @param truths    Truth sequence [K × n] (same truth for all N runs)
 * @param P_seq     Covariance sequence [K × n²] (or single for all K if repeated)
 * @param N         Number of Monte Carlo runs
 * @param K         Number of time steps per run
 * @param n         State dimension
 * @return ANEES value
 */
double metrics_anees(const double *estimates, const double *truths,
                      const double *P_seq, int N, int K, int n);

/**
 * Compute the ANEES acceptance interval for given significance level.
 *
 * Interval = [n − z_{α/2}·√(2n/(N·K)), n + z_{α/2}·√(2n/(N·K))]
 *
 * If ANEES falls within this interval, the filter is considered consistent
 * at level α.
 *
 * Complexity: O(1)
 *
 * @param n    State dimension
 * @param N    Monte Carlo runs
 * @param K    Time steps
 * @param alpha Significance level (e.g., 0.05 for 95%)
 * @param lower Output: lower bound
 * @param upper Output: upper bound
 */
void metrics_anees_interval(int n, int N, int K, double alpha,
                              double *lower, double *upper);

/**
 * Check filter consistency from ANEES.
 *
 * Returns 1 if consistent, 0 otherwise.
 * Complexity: O(1)
 */
int metrics_is_consistent(double anees, int n, int N, int K, double alpha);

/* ============================================================================
 * L1 — Multi-target tracking metrics
 * ============================================================================
 */

/**
 * Track purity — what fraction of a track's life is associated with the
 * same true target.
 *
 * purity = max_j (count of track i assigned to truth j) / total track length
 *
 * A pure track maintains identity throughout its life.
 *
 * Complexity: O(T · K) where T = tracks, K = time steps
 *
 * @param truth_assignments Assignment[t][k] = truth ID at time k
 *        (-1 = unassigned, 0 = false track)
 * @param track_lengths     Number of time steps for each track [T]
 * @param T                 Number of tracks
 * @param purities          Output: purity for each track [T]
 */
void metrics_track_purity(const int *const *truth_assignments,
                           const int *track_lengths, int T,
                           double *purities);

/**
 * Average track purity across all tracks.
 * Complexity: O(T)
 */
double metrics_average_purity(const double *purities, int T);

/**
 * Track continuity — how many tracks are needed to cover one true target.
 *
 * A true target requiring many track fragments indicates poor tracking.
 *
 * Continuity for target j = 1 / (number of tracks that covered it)
 * Perfect: 1.0 (one track for whole life)
 * Fragmented: e.g., 0.25 (4 tracks covered it)
 *
 * Complexity: O(T·K)
 */
void metrics_track_continuity(const int *const *truth_assignments,
                               const int *track_lengths, int T,
                               int num_truths, double *continuity);

/**
 * Track fragmentation count — total number of track switches for each target.
 * Complexity: O(T·K)
 */
void metrics_track_fragmentation(const int *const *truth_assignments,
                                  const int *track_lengths, int T,
                                  int num_truths, int *frag_count);

/**
 * False track rate — number of false tracks (not associated with any truth)
 * per unit time.
 *
 * FTR = N_false / T_total
 *
 * Complexity: O(T)
 */
double metrics_false_track_rate(const int *const *truth_assignments,
                                 const int *track_lengths, int T,
                                 double total_time);

/**
 * Track initiation delay — average time between first truth appearance
 * and first track on that truth.
 *
 * Complexity: O(T·K)
 */
double metrics_initiation_delay(const double *truth_appear_time,
                                 const double *track_start_time,
                                 const int *track_to_truth, int T);

/**
 * Track termination delay — average time between truth disappearance
 * and track deletion.
 *
 * Complexity: O(T)
 */
double metrics_termination_delay(const double *truth_disappear_time,
                                  const double *track_end_time,
                                  const int *track_to_truth, int T);

/* ============================================================================
 * L5 — OSPA distance for multi-object evaluation
 * ============================================================================
 */

/**
 * Optimal SubPattern Assignment (OSPA) distance between two sets of states.
 *
 * Measures the distance between the estimated multi-object state and the
 * true multi-object state, accounting for cardinality errors.
 *
 * OSPA^{(c)}_p(X̂, X) = [ (1/m) · ( min_{π∈Π_m} Σ_{i=1}^{n̂} d^{(c)}(x̂_i, x_{π(i)})^p
 *                                  + c^p · (m − n̂) ) ]^(1/p)
 *
 * where c is the cutoff distance, p is the order, and m ≥ n̂ are the set sizes.
 *
 * OSPA jointly penalizes localization error and cardinality error.
 *
 * Reference: Schuhmacher, Vo, Vo (2008)
 * Complexity: O(m³) for Hungarian algorithm on d^{(c)} matrix
 *
 * @param X_est  Estimated states [n_est × dim] (row-major)
 * @param n_est  Number of estimated objects
 * @param X_true True states [n_true × dim]
 * @param n_true Number of true objects
 * @param dim    State dimension (position only, typically 2 or 3)
 * @param c      Cutoff distance (e.g., 50m)
 * @param p      Order parameter (1 or 2, typically)
 * @return OSPA distance
 */
double metrics_ospa(const double *X_est, int n_est,
                     const double *X_true, int n_true,
                     int dim, double c, int p);

/**
 * Time-averaged OSPA across K time steps.
 *
 * Complexity: O(K·max(n_est, n_true)³)
 */
double metrics_ospa_avg(const double *X_est_seq, const int *n_est_seq,
                         const double *X_true_seq, const int *n_true_seq,
                         int K, int dim, double c, int p);

/* ============================================================================
 * L5 — CLEAR MOT metrics
 * ============================================================================
 */

/**
 * Compute CLEAR MOT metrics (MOTA, MOTP).
 *
 * MOTA = 1 − (Σₖ (FNₖ + FPₖ + IDSWₖ)) / (Σₖ GTₖ)
 * MOTP = Σ_{i,k} d_{i,k} / Σₖ cₖ
 *
 * where FN = false negatives (missed detections), FP = false positives,
 * IDSW = identity switches, GT = ground truth count, d = localization error,
 * c = number of matches.
 *
 * Reference: Bernardin & Stiefelhagen (2008)
 * Complexity: O(K·T²) for matching
 *
 * @param assignments Assignment[k][t] = truth or -1 [K × T]
 * @param track_lengths Track lengths [T]
 * @param distances   Distance for each matched pair [K × T]
 * @param K           Number of frames
 * @param T           Max tracks per frame
 * @param mota        Output MOTA
 * @param motp        Output MOTP
 * @param id_switches Output identity switch count
 */
void metrics_clear_mot(const int *const *assignments,
                        const int *track_lengths,
                        const double *const *distances,
                        int K, int T,
                        double *mota, double *motp, int *id_switches);

#ifdef __cplusplus
}
#endif

#endif /* TRACK_METRICS_H */
