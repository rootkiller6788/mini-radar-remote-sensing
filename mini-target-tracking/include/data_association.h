/**
 * data_association.h — Measurement-to-Track Data Association Algorithms
 *
 * Covers: L2 (Core Concepts) — Gating, assignment problem
 *         L5 (Algorithms/Methods) — NN, GNN (Hungarian), PDA, JPDA, MHT
 *         L6 (Canonical Problems) — Multi-target tracking data association
 *         L8 (Advanced) — Auction algorithm, multi-scan MHT
 *
 * References:
 *   - Bar-Shalom, Fortmann "Tracking and Data Association" (1988)
 *   - Munkres, J. "Algorithms for Assignment and Transportation Problems" (1957)
 *   - Kuhn, H.W. "The Hungarian Method" (1955)
 *   - Reid, D.B. "An Algorithm for Tracking Multiple Targets" (1979)
 *   - Blackman, Popoli "Design and Analysis of Modern Tracking Systems" (1999)
 *
 * Curriculum:
 *   - MIT 6.450, Stanford EE359, Berkeley EE123, ETH 227-0436
 */

#ifndef DATA_ASSOCIATION_H
#define DATA_ASSOCIATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"
#include <stddef.h>

/* ============================================================================
 * L1 — Cost matrix structure
 * ============================================================================
 */

/** Cost matrix for measurement-to-track association */
typedef struct {
    int     num_tracks;       /**< Number of tracks */
    int     num_measurements; /**< Number of measurements */
    double *cost;             /**< Cost matrix [num_tracks × num_measurements], row-major */
    int    *gate_mask;        /**< Gate mask [num_tracks × num_measurements], 1 = in gate */
    double  gate_threshold;   /**< Mahalanobis distance gate threshold */
    double  beta_F;           /**< False alarm spatial density */
    double  P_D;              /**< Probability of detection */
} cost_matrix_t;

/**
 * Assignment result structure — which measurement goes to which track.
 */
typedef struct {
    int     num_tracks;
    int     num_measurements;
    int    *assignment;       /**< assignment[i] = j means track i → measurement j
                                   -1 means track i is unassigned (missed detection) */
    int    *unassigned_meas;  /**< List of unassigned measurement indices (new tracks/clutter) */
    int     num_unassigned;
    double  total_cost;       /**< Total assignment cost */
} assignment_result_t;

/* ============================================================================
 * L2 — Gate computation for cost matrix
 * ============================================================================
 */

/**
 * Compute the Mahalanobis-based cost matrix between all tracks and measurements.
 *
 * cost(i, j) = d_M²(track_i, meas_j) if measurement j is in track_i's gate,
 *              ∞ otherwise.
 *
 * Complexity: O(N_tracks × N_meas × m³) where m = measurement dimension
 *
 * @param tracks        Array of track_t pointers
 * @param num_tracks    Number of tracks
 * @param measurements  Array of measurement_t pointers
 * @param num_meas      Number of measurements
 * @param matrix        Output cost matrix (will be allocated or resized)
 */
void association_build_cost_matrix(const track_t *tracks, int num_tracks,
                                    const measurement_t *measurements, int num_meas,
                                    cost_matrix_t *matrix);

/**
 * Fill the gate mask: mask[i][j] = 1 if meas_j falls within track_i's gate.
 * Complexity: O(N_tracks × N_meas)
 */
void association_build_gate_mask(cost_matrix_t *matrix, double gate_threshold);

/* ============================================================================
 * L5 — Nearest Neighbor (NN) association
 * ============================================================================
 */

/**
 * Nearest Neighbor (NN) data association.
 *
 * For each track, assign the nearest (lowest cost) measurement that is in gate.
 * Each measurement can be assigned to at most one track (greedy):
 * tracks are processed in order, bypassing already-assigned measurements.
 *
 * Algorithm (greedy NN):
 *   1. Sort tracks by track score (better tracks get first pick)
 *   2. For each track: find min-cost measurement in gate, not yet assigned
 *   3. If cost < gate_threshold, assign; otherwise leave unassigned
 *
 * Complexity: O(N_tracks × N_meas × log(N_tracks)) with sorting
 *
 * Reference: Blackman (1999), Ch. 7
 */
assignment_result_t *nn_associate(const cost_matrix_t *matrix,
                                   const double *track_scores);

/**
 * Global Nearest Neighbor (GNN) via Hungarian/Munkres algorithm.
 *
 * Solves the optimal assignment problem minimizing total cost:
 *   min Σᵢ Σⱼ C_{ij} · x_{ij}
 *   s.t. Σᵢ x_{ij} ≤ 1 (each measurement assigned at most once)
 *        Σⱼ x_{ij} ≤ 1 (each track assigned at most once)
 *        x_{ij} ∈ {0, 1}
 *
 * The Hungarian algorithm solves this in polynomial time.
 *
 * Reference: Kuhn (1955), Munkres (1957)
 * Complexity: O(max(N_t, N_m)³)
 *
 * @param matrix Cost matrix (modified in place as working buffer)
 */
assignment_result_t *gnn_hungarian_associate(const cost_matrix_t *matrix);

/**
 * GNN with gating: only consider assignments within the gate.
 * Entries outside the gate are set to an infinite cost.
 * The Hungarian algorithm then naturally avoids them.
 *
 * Complexity: same as Hungarian
 */
assignment_result_t *gnn_gated_associate(const cost_matrix_t *matrix);

/* ============================================================================
 * L5 — Probabilistic Data Association (PDA)
 * ============================================================================
 */

/**
 * Probabilistic Data Association for a single track in clutter.
 *
 * Computes association probabilities β_i for each validated measurement:
 *
 * β₀ = b / (b + Σⱼ eⱼ)     (probability that no measurement is from target)
 * βⱼ = eⱼ / (b + Σⱼ eⱼ)    (probability meas_j is from target)
 *
 * where eⱼ = exp(−dⱼ²/2) / (P_G · √|S| · (2π)^{n_z/2})
 *       b  = (1 − P_D·P_G) · β_F · (2π)^{n_z/2} · √|S| / P_D
 *
 * The PDA state update uses the combined innovation:
 *   ν = Σⱼ βⱼ·νⱼ
 *   P = β₀·P_pred + (1−β₀)·P_corrected + P̃
 *
 * where P_corrected is the standard KF updated covariance and
 * P̃ is the "spread of the innovations" term.
 *
 * Reference: Bar-Shalom, Y., Tse, E. "Tracking in a cluttered environment
 *            with probabilistic data association" Automatica (1975)
 * Complexity: O(m_k³ + m_k²·n) where m_k = validated measurements
 *
 * @param track          Track to update
 * @param measurements   Array of validated measurements
 * @param num_validated  Number of measurements in gate
 * @param valid_indices  Indices into the full measurement array
 * @param P_D            Probability of detection
 * @param beta_F         False alarm density
 * @param x_updated      Output: updated state estimate [state_dim]
 * @param P_updated      Output: updated covariance [state_dim²]
 * @return 0 on success, -1 on failure
 */
int pda_update(track_t *track, const measurement_t *measurements,
               int num_validated, const int *valid_indices,
               double P_D, double beta_F,
               double *x_updated, double *P_updated);

/* ============================================================================
 * L5 — Joint Probabilistic Data Association (JPDA)
 * ============================================================================
 */

/**
 * Joint Probabilistic Data Association filter.
 *
 * Unlike PDA which processes tracks independently, JPDA computes joint
 * association probabilities considering all tracks simultaneously.
 *
 * Feasible joint event θ: a hypothesis assigning each measurement to at most
 * one track, and each track to at most one measurement.
 *
 * β_{jt} = P(θ_{jt} | Z^k) = Σ_{θ: θ_{jt}∈θ} P(θ | Z^k)
 *
 * The probability of each feasible joint event:
 *   P(θ | Z^k) = (1/c) · ∏ⱼ (P_D)^{δ_j} · (1−P_D)^{1−δ_j}
 *                        · ∏_{t:τ(t,j)>0} p_{t,τ(t,j)}
 *
 * Reference: Fortmann, Bar-Shalom, Scheffe "Sonar tracking of multiple
 *            targets using JPDA" IEEE JOE (1983)
 * Complexity: worst-case O(N!·M!) but approximated with efficient gating
 *
 * @param tracks        Array of track pointers
 * @param num_tracks    Number of tracks
 * @param measurements  Array of measurements
 * @param num_meas      Number of measurements
 * @param P_D           Detection probability
 * @param beta_F        False alarm density
 * @param beta_matrix   Output: β[t][m] association probabilities [num_tracks × num_meas]
 * @param beta0         Output: β[t][0] = prob track t had no detection [num_tracks]
 * @return 0 on success, -1 on memory failure
 */
int jpda_compute_beta(track_t *tracks, int num_tracks,
                      const measurement_t *measurements, int num_meas,
                      double P_D, double beta_F,
                      double **beta_matrix, double *beta0);

/**
 * Apply JPDA association probabilities to update all tracks.
 * Uses the combined innovation weighted by β probabilities.
 *
 * Complexity: O(N_tracks × N_meas × n³)
 */
int jpda_update_tracks(track_t *tracks, int num_tracks,
                       const measurement_t *measurements, int num_meas,
                       const double *const *beta_matrix, const double *beta0);

/* ============================================================================
 * L5 — Multiple Hypothesis Tracking (MHT)
 * ============================================================================
 */

/** Single hypothesis: a set of track-measurement assignments for one scan */
typedef struct {
    int     scan_id;            /**< Scan number */
    int    *assignment;         /**< Track-to-measurement mapping */
    int     num_tracks;
    double  hypothesis_score;   /**< Cumulative LLR score */
    int     parent_hypothesis;  /**< Index of parent hypothesis (-1 = root) */
} mht_hypothesis_t;

/** MHT hypothesis tree */
typedef struct {
    mht_hypothesis_t *hypotheses;     /**< Array of hypotheses */
    int               num_hypotheses; /**< Current count */
    int               max_hypotheses; /**< Maximum capacity */
    int               best_index;     /**< Index of best (highest score) hypothesis */
    int               max_depth;      /**< Maximum scan depth for pruning */
} mht_tree_t;

/**
 * Initialize MHT hypothesis tree.
 * Complexity: O(1)
 */
void mht_init(mht_tree_t *tree, int max_hypotheses, int max_depth);

/**
 * Generate new hypotheses from the current best hypotheses.
 *
 * For each surviving hypothesis and each feasible association,
 * create a child hypothesis with updated score.
 *
 * Complexity: O(N_hypotheses × N_associations)
 *
 * @return Number of new hypotheses generated
 */
int mht_generate_hypotheses(mht_tree_t *tree,
                             const cost_matrix_t *matrix,
                             const track_t *tracks, int num_tracks,
                             double P_D, double beta_F);

/**
 * Prune MHT hypotheses: keep top-K by score, delete the rest.
 * Returns number of hypotheses pruned.
 *
 * Complexity: O(N_hypotheses × log(N_hypotheses))
 */
int mht_prune(mht_tree_t *tree, int keep_top_k);

/**
 * Get the best hypothesis assignment.
 * Complexity: O(N_tracks)
 */
const mht_hypothesis_t *mht_get_best(const mht_tree_t *tree);

/**
 * Free MHT tree memory.
 * Complexity: O(1)
 */
void mht_free(mht_tree_t *tree);

/* ============================================================================
 * L8 — Auction algorithm for assignment
 * ============================================================================
 */

/**
 * Auction algorithm for the assignment problem.
 *
 * An alternative to the Hungarian algorithm with better average-case
 * performance for large sparse problems.
 *
 * Each "person" (track) bids on its preferred "item" (measurement).
 * Prices (thresholds) are updated dynamically.
 *
 * Reference: Bertsekas, D.P. "The auction algorithm: A distributed
 *            relaxation method for the assignment problem" (1988)
 * Complexity: potentially O(N³) worst-case, often O(N²) average
 *
 * @param cost   Cost matrix [n_tracks × n_meas], row-major
 *                (note: auction maximizes profit = −cost, so smaller cost = better)
 * @param n_tracks Number of tracks
 * @param n_meas   Number of measurements
 * @param assignment Output: assignment[i] = j or -1 if unassigned
 * @param epsilon   Auction precision parameter (smaller = more accurate, slower)
 * @return 0 on success
 */
int auction_assign(const double *cost, int n_tracks, int n_meas,
                   int *assignment, double epsilon);

/* ============================================================================
 * L2 — Assignment result utilities
 * ============================================================================
 */

/**
 * Allocate and initialize an assignment result.
 * Complexity: O(N_tracks + N_meas)
 */
assignment_result_t *assignment_alloc(int num_tracks, int num_measurements);

/**
 * Free assignment result.
 * Complexity: O(1)
 */
void assignment_free(assignment_result_t *result);

/**
 * Count the number of associated pairs in the result.
 * Complexity: O(N_tracks)
 */
int assignment_count_pairs(const assignment_result_t *result);

/**
 * Print assignment result for debugging.
 */
void assignment_print(const assignment_result_t *result);

#ifdef __cplusplus
}
#endif

#endif /* DATA_ASSOCIATION_H */
