#ifndef PHASED_ARRAY_ADAPTIVE_H
#define PHASED_ARRAY_ADAPTIVE_H

#include "phased_array.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L8: Adaptive Beamforming — Data-Dependent Weight Optimization
 *
 * Adaptive beamformers adjust weights based on received signal statistics
 * to suppress interference while preserving the desired signal.
 *
 * Algorithms:
 *   - SMI (Sample Matrix Inversion): Covariance estimation from snapshots
 *   - MVDR/Capon: Minimum Variance Distortionless Response
 *   - LCMV: Linearly Constrained Minimum Variance (generalized MVDR)
 *
 * Reference: Van Trees (2002) Optimum Array Processing, Part IV.
 * ============================================================================ */

/**
 * Initialize an SMI-based adaptive beamformer.
 *
 * Allocates: weight vector, steering vector, covariance matrix,
 *            inverse covariance matrix.
 *
 * @param num_elements M = number of array elements.
 * @param num_snapshots K = number of temporal snapshots.
 * @return Initialized beamformer (caller frees via pa_smi_beamformer_free).
 */
pa_smi_beamformer_t *pa_smi_beamformer_init(uint32_t num_elements,
                                             uint32_t num_snapshots);

/**
 * Free an SMI beamformer.
 */
void pa_smi_beamformer_free(pa_smi_beamformer_t *bf);

/**
 * Estimate sample covariance matrix from snapshot data.
 *
 * R̂ = (1/K) Σ_{k=0}^{K-1} x[k] x[k]^H   [MLE for i.i.d. Gaussian]
 *
 * Applies diagonal loading for regularization:
 *   R̂_loaded = R̂ + γ · (tr(R̂)/M) · I
 *
 * Convergence requirement: K ≥ 2M for invertibility.
 *
 * @param bf Beamformer state (updated in place).
 * @param snapshots 2D array: snapshots[k][m], k=0..K-1, m=0..M-1.
 */
void pa_smi_estimate_covariance(pa_smi_beamformer_t *bf,
                                 const double complex **snapshots);

/**
 * Initialize an MVDR (Capon) beamformer.
 *
 * @param num_elements Number of array elements.
 * @param num_snapshots Number of temporal snapshots.
 * @return Initialized MVDR beamformer.
 */
pa_mvdr_beamformer_t *pa_mvdr_beamformer_init(uint32_t num_elements,
                                               uint32_t num_snapshots);

/**
 * Free an MVDR beamformer.
 */
void pa_mvdr_beamformer_free(pa_mvdr_beamformer_t *bf);

/**
 * Compute MVDR weight vector.
 *
 * Optimization problem:
 *   minimize  w^H R w      (minimize output power)
 *   subject to  w^H a = 1  (unity gain in look direction)
 *
 * Closed-form solution:
 *   w_MVDR = R^{-1} a / (a^H R^{-1} a)
 *
 * The denominator a^H R^{-1} a is a real scalar (R is Hermitian).
 *
 * MVDR places spatial nulls in interferer directions while maintaining
 * a fixed response toward the desired signal.
 *
 * @param bf MVDR beamformer (updated in place).
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Look direction θ (rad).
 * @param phi_steer Look direction φ (rad).
 *
 * Reference: Capon (1969) Proc. IEEE; Van Trees (2002) §6.2.
 */
void pa_mvdr_compute_weights(pa_mvdr_beamformer_t *bf,
                              const pa_array_config_t *config,
                              const pa_element_t *elements,
                              double theta_steer, double phi_steer);

/**
 * Compute Capon spatial power spectrum.
 *
 * P_Capon(θ,φ) = 1 / (a(θ,φ)^H R^{-1} a(θ,φ))
 *
 * The Capon spectrum provides high-resolution DOA estimation:
 * peaks correspond to source directions with peak height ∝ source power.
 * Resolution exceeds the Rayleigh (beamwidth) limit.
 *
 * @param bf MVDR beamformer (must have R^{-1} computed).
 * @param config Array config.
 * @param elements Array elements.
 * @param spectrum Output: P_Capon values, size = num_angles.
 * @param theta_vals Theta evaluation angles.
 * @param phi_vals Phi evaluation angles.
 * @param num_angles Number of angular bins.
 */
void pa_mvdr_capon_spectrum(pa_mvdr_beamformer_t *bf,
                             const pa_array_config_t *config,
                             const pa_element_t *elements,
                             double *spectrum,
                             const double *theta_vals,
                             const double *phi_vals,
                             uint32_t num_angles);

/**
 * Compute MVDR output SINR.
 *
 * SINR = |w^H a|² / (w^H R w − |w^H a|² · σ_s²)
 *
 * At optimality: SINR_opt = σ_s² · a^H R_{i+n}^{-1} a
 *
 * The number of interferers that can be nulled ≤ M−1
 * (M degrees of freedom minus the look-direction constraint).
 *
 * @param bf MVDR beamformer with computed weights.
 * @return SINR in dB.
 */
double pa_mvdr_output_sinr(const pa_mvdr_beamformer_t *bf);

#ifdef __cplusplus
}
#endif

#endif /* PHASED_ARRAY_ADAPTIVE_H */
