/**
 * track_fusion.h — Multi-Sensor Track-to-Track Fusion
 *
 * Covers: L2 (Core Concepts) — Track-to-track association, fusion architectures
 *         L3 (Mathematical Structures) — Covariance intersection, information fusion
 *         L5 (Algorithms/Methods) — CI, WCF, information matrix fusion, tracklet
 *         L6 (Canonical Problems) — Distributed vs centralized fusion
 *         L8 (Advanced Topics) — Decentralized fusion, rumor propagation
 *
 * References:
 *   - Bar-Shalom, Y. "On the track-to-track correlation problem" IEEE TAC (1981)
 *   - Julier, S.J., Uhlmann, J.K. "A non-divergent estimation algorithm in the
 *     presence of unknown correlations" Proc. ACC (1997)
 *   - Chang, K.C., Saha, R.K., Bar-Shalom, Y. "On optimal track-to-track fusion"
 *     IEEE T-AES (1997)
 *   - Chong, C.Y., Mori, S., Chang, K.C. "Distributed multitarget multisensor
 *     tracking" in Bar-Shalom (ed.) Multitarget-Multisensor Tracking (1990)
 *
 * Curriculum:
 *   - MIT 6.450, Stanford EE359, ETH 227-0436
 *   - Georgia Tech ECE 6601, Michigan EECS 455
 */

#ifndef TRACK_FUSION_H
#define TRACK_FUSION_H

#ifdef __cplusplus
extern "C" {
#endif

#include "track_core.h"

/* ============================================================================
 * L1 — Fusion architecture enumeration
 * ============================================================================
 */

/** Fusion architecture types */
typedef enum {
    FUSION_CENTRALIZED     = 0,  /**< Raw measurements sent to fusion center */
    FUSION_DISTRIBUTED_T2T = 1,  /**< Track-to-track fusion (tracks fused) */
    FUSION_HYBRID          = 2,  /**< Mixed: some tracks, some measurements */
    FUSION_DECENTRALIZED   = 3   /**< No central node, peer-to-peer fusion */
} fusion_architecture_t;

/** Track-to-track association result */
typedef struct {
    int     n_local;        /**< Number of local tracks */
    int     n_remote;       /**< Number of remote tracks */
    int    *assignment;     /**< assignment[i] = j: local track i → remote track j */
    double *association_cost; /**< Per-pair association cost */
    int     n_pairs;        /**< Number of associated pairs */
} t2t_association_t;

/* ============================================================================
 * L2 — Track-to-track association
 * ============================================================================
 */

/**
 * Associate local and remote tracks based on state similarity.
 *
 * Distance metric: d² = (x̂₁ − x̂₂)ᵀ·(P₁ + P₂)⁻¹·(x̂₁ − x̂₂)
 *
 * This is the Mahalanobis distance between two track estimates, accounting
 * for the fact that tracks from different sensors have independent errors
 * (assuming no common process noise correlation).
 *
 * When common process noise creates correlation, the distance metric should
 * use the cross-covariance P₁₂ instead of P₁ + P₂.
 *
 * Reference: Bar-Shalom (1981)
 * Complexity: O(N₁·N₂·n³) for n-dimensional state
 *
 * @param local_tracks  Local track array
 * @param n_local       Number of local tracks
 * @param remote_tracks Remote track array
 * @param n_remote      Number of remote tracks
 * @param threshold     Chi-squared gate threshold (e.g., 16.27 for p=0.999, dof=4)
 * @return Association result (free with t2t_association_free)
 */
t2t_association_t *t2t_associate(const track_t *local_tracks, int n_local,
                                  const track_t *remote_tracks, int n_remote,
                                  double threshold);

/**
 * Free T2T association result.
 * Complexity: O(1)
 */
void t2t_association_free(t2t_association_t *assoc);

/* ============================================================================
 * L3+L5 — Fusion algorithms
 * ============================================================================
 */

/**
 * Covariance Intersection (CI) fusion.
 *
 * Fuses two estimates (x̂₁, P₁) and (x̂₂, P₂) when the cross-correlation
 * between them is unknown. Produces a consistent (conservative) fused estimate.
 *
 * CI update:
 *   P⁻¹_CI = ω·P₁⁻¹ + (1−ω)·P₂⁻¹
 *   x̂_CI = P_CI · (ω·P₁⁻¹·x̂₁ + (1−ω)·P₂⁻¹·x̂₂)
 *
 * where ω ∈ [0, 1] minimizes det(P_CI) or tr(P_CI).
 *
 * Key property: Always consistent regardless of unknown correlation.
 * P_CI ≥ P_true (in the positive semi-definite sense) for any correlation.
 *
 * Reference: Julier & Uhlmann (1997)
 * Complexity: O(n³) for matrix inversions, O(k·n³) for ω optimization
 *
 * @param x1, P1 First estimate (state [n], covariance [n²])
 * @param x2, P2 Second estimate
 * @param n      State dimension
 * @param x_fused Output fused state [n]
 * @param P_fused Output fused covariance [n²]
 * @param omega  Output: optimal weight (pass NULL if not needed)
 * @return 0 on success
 */
int fusion_covariance_intersection(const double *x1, const double *P1,
                                    const double *x2, const double *P2,
                                    int n, double *x_fused, double *P_fused,
                                    double *omega);

/**
 * Weighted Covariance Fusion (WCF) — optimal when cross-covariance is known.
 *
 * Assumes: P₁₂ = E[(x̂₁−x)(x̂₂−x)ᵀ] is known.
 *
 * Fused estimate:
 *   x̂_WCF = x̂₁ + (P₁ − P₁₂)·(P₁ + P₂ − P₁₂ − P₁₂ᵀ)⁻¹·(x̂₂ − x̂₁)
 *   P_WCF = P₁ − (P₁ − P₁₂)·(P₁ + P₂ − P₁₂ − P₁₂ᵀ)⁻¹·(P₁ − P₁₂ᵀ)
 *
 * This is the optimal (minimum variance) linear fusion when P₁₂ is known.
 *
 * Reference: Bar-Shalom (1981), Chang et al. (1997)
 * Complexity: O(n³)
 *
 * @param x1, P1 First estimate
 * @param x2, P2 Second estimate
 * @param P12    Cross-covariance E[(x̂₁−x)(x̂₂−x)ᵀ] [n²]
 * @param n      State dimension
 * @param x_fused Output fused state
 * @param P_fused Output fused covariance
 * @return 0 on success, -1 if fusion is degenerate
 */
int fusion_weighted_covariance(const double *x1, const double *P1,
                                const double *x2, const double *P2,
                                const double *P12,
                                int n, double *x_fused, double *P_fused);

/**
 * Information Matrix Fusion — optimal for distributed estimation without
 * double-counting of common information.
 *
 * Information form:
 *   Y_fused = Y₁ + Y₂ − Y_common
 *   ŷ_fused = ŷ₁ + ŷ₂ − ŷ_common
 *
 * where Y = P⁻¹ is the information matrix and ŷ = P⁻¹·x̂ is the information state.
 *
 * When common information is zero (independent sensors with independent
 * process noise), this reduces to simple addition:
 *   Y_fused = Y₁ + Y₂
 *   ŷ_fused = ŷ₁ + ŷ₂
 *
 * Reference: Chong, Mori, Chang (1990)
 * Complexity: O(n³)
 *
 * @param Y1, y1 First estimate in information form
 * @param Y2, y2 Second estimate in information form
 * @param Y_common, y_common Common information to subtract (NULL = zero)
 * @param n      State dimension
 * @param Y_fused Output information matrix [n²]
 * @param y_fused Output information state [n]
 * @param x_fused Output fused state in Cartesian form [n] (NULL if not needed)
 * @param P_fused Output fused covariance [n²] (NULL if not needed)
 * @return 0 on success
 */
int fusion_information_matrix(const double *Y1, const double *y1,
                               const double *Y2, const double *y2,
                               const double *Y_common, const double *y_common,
                               int n,
                               double *Y_fused, double *y_fused,
                               double *x_fused, double *P_fused);

/* ============================================================================
 * L5 — Tracklet fusion
 * ============================================================================
 */

/**
 * Fuse two tracklets (short track segments) into a single track.
 *
 * A tracklet is a segment of track with a start and end time.
 * Tracklet fusion differs from track fusion because the two tracklets
 * may cover different time intervals.
 *
 * This function retrofits tracklet 2's information into tracklet 1's
 * estimate by re-running the Kalman filter forward/backward.
 *
 * Reference: Drummond, O.E. "Tracklet fusion for netted sensors" (1995)
 * Complexity: O(n³·k) for k measurements in tracklet 2
 *
 * @param x1, P1 Tracklet 1 state/covariance at time t1, updated in place
 * @param z2      Tracklet 2 measurement sequence
 * @param R2      Tracklet 2 measurement covariances (may repeat same R)
 * @param F       State transition matrix
 * @param Q       Process noise covariance
 * @param H       Measurement matrix
 * @param n       State dimension
 * @param m       Measurement dimension
 * @param k2      Number of measurements in tracklet 2
 * @param dt      Time step per measurement in tracklet 2
 * @return 0 on success
 */
int fusion_tracklet(double *x1, double *P1,
                     const double *z2, const double *R2,
                     const double *F, const double *Q, const double *H,
                     int n, int m, int k2, double dt);

/* ============================================================================
 * L6 — Distributed fusion consistency maintenance
 * ============================================================================
 */

/**
 * Compute the cross-covariance P₁₂ between two local estimates that share
 * common process noise history.
 *
 * P₁₂ evolves as:
 *   P₁₂(k+1) = (I − K₁·H)·F·P₁₂(k)·Fᵀ·(I − K₂·H)ᵀ + (I − K₁·H)·Q·(I − K₂·H)ᵀ
 *
 * This is needed for optimal WCF fusion between communicating local trackers.
 *
 * Reference: Bar-Shalom (1981)
 * Complexity: O(n³)
 *
 * @param P12 Current cross-covariance [n²], updated in place
 * @param K1  Kalman gain of tracker 1 [n×m]
 * @param K2  Kalman gain of tracker 2 [n×m]
 * @param F   Common transition matrix [n²]
 * @param Q   Common process noise [n²]
 * @param H   Common measurement matrix [m×n]
 * @param n   State dimension
 * @param m   Measurement dimension
 */
void fusion_cross_covariance_update(double *P12,
                                     const double *K1, const double *K2,
                                     const double *F, const double *Q,
                                     const double *H, int n, int m);

/**
 * Compute the decorrelated pseudo-measurement for track fusion.
 *
 * Given two local tracks (x̂₁,P₁), (x̂₂,P₂) and their cross-covariance P₁₂,
 * form the decorrelated pair:
 *
 *   z_dec = x̂₂ − P₁₂ᵀ·P₁⁻¹·x̂₁
 *   R_dec = P₂ − P₁₂ᵀ·P₁⁻¹·P₁₂
 *
 * Then (z_dec, R_dec) can be treated as an independent measurement
 * with no cross-correlation, fed into tracker 1's filter.
 *
 * Complexity: O(n³)
 *
 * @param x1, P1 First track estimate
 * @param x2, P2 Second track estimate
 * @param P12    Cross-covariance [n²]
 * @param n      State dimension
 * @param z_dec  Output decorrelated measurement [n]
 * @param R_dec  Output decorrelated covariance [n²]
 */
void fusion_decorrelate(const double *x1, const double *P1,
                         const double *x2, const double *P2,
                         const double *P12, int n,
                         double *z_dec, double *R_dec);

/* ============================================================================
 * L8 — Rumor-robust decentralized fusion
 * ============================================================================
 */

/**
 * Check if information from a remote node is "novel" (not already incorporated).
 * Prevents "rumor propagation" — the same information cycling through the
 * network and getting double-counted.
 *
 * Uses information-theoretic distance:
 *   d(Y₁, Y₂) = tr(Y₁⁻¹·Y₂) − log|Y₁⁻¹·Y₂| − n
 *
 * If d < threshold, the remote information is considered to be already known.
 *
 * Reference: Chong et al. (1990), Sec. 5
 * Complexity: O(n³)
 *
 * @param Y_local  Local information matrix [n²]
 * @param Y_remote Remote information matrix [n²]
 * @param n        State dimension
 * @param threshold Novelty threshold
 * @return 1 if novel (should fuse), 0 if redundant
 */
int fusion_is_novel_info(const double *Y_local, const double *Y_remote,
                          int n, double threshold);

#ifdef __cplusplus
}
#endif

#endif /* TRACK_FUSION_H */
