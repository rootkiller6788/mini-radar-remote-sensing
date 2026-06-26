#ifndef PHASED_ARRAY_SYNTHESIS_H
#define PHASED_ARRAY_SYNTHESIS_H

#include "phased_array.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L5: Array Pattern Synthesis — Sidelobe Control via Amplitude Tapering
 *
 * Synthesis Methods:
 *   1. Dolph-Chebyshev: Constant sidelobe level, minimal beamwidth
 *   2. Taylor one-parameter: Near-in constant SLL, far-out monotonic decay
 *   3. Binomial: Zero sidelobes, maximum beamwidth
 *
 * Reference: Balanis (2016) §6.8
 * ============================================================================ */

/**
 * Compute Dolph-Chebyshev amplitude coefficients.
 *
 * Produces an excitation with all sidelobes at exactly sll_db.
 * Uses Chebyshev polynomial T_{N-1}(x) mapping.
 *
 * Theory:
 *   x0 = cosh(acosh(R0)/(N-1))  where R0 = 10^{|sll_db|/20}
 *   AF(ψ) = T_{N-1}(x0 · cos(ψ/2))
 *
 * For N odd (N=2M+1):
 *   a_{M±k} = (1/N)[R0 + 2 Σ T_{N-1}(x0·cos(pπ/N))·cos(2k pπ/N)]
 *
 * For N even (N=2M):
 *   a_{M-1-k} = a_{M+k} = (2/N) Σ T_{N-1}(x0·cos((p-½)π/N))
 *                          · cos((2k+1)(p-½)π/N)
 *
 * @param num_elements N ≥ 2.
 * @param sll_db Desired sidelobe level in dB (negative, e.g. -30.0).
 * @return Amplitude weights (caller frees). NULL on error.
 *
 * Reference: Dolph (1946) Proc. IRE; Balanis §6.8.3.
 */
double *pa_dolph_chebyshev_weights(uint32_t num_elements, double sll_db);

/**
 * Compute Taylor one-parameter line-source distribution.
 *
 * w(z) = I₀(πB √(1−(2z/L)²)) / I₀(πB)
 * where B = acosh(R₀)/π, R₀ = 10^{|sll_db|/20}
 *
 * Benefits over Dolph-Chebyshev:
 *   - Far-out sidelobes decay monotonically (~1/u)
 *   - More practical edge taper (no extreme peaking)
 *   - Suitable for large arrays (N > 20)
 *
 * @param num_elements N ≥ 3.
 * @param sll_db Desired peak sidelobe level (dB, negative).
 * @param nbar Number of equal-level sidelobes before decay (≥ 2, typical 3–8).
 * @return Amplitude weights (caller frees).
 *
 * Reference: Taylor (1955) IRE Trans. AP-3; Balanis §6.8.4.
 */
double *pa_taylor_weights(uint32_t num_elements, double sll_db, uint32_t nbar);

/**
 * Compute binomial (maximally flat) array coefficients.
 *
 * a_n = C(N−1, n) / 2^{N-1}
 *
 * AF(ψ) = (2 cos(ψ/2))^{N-1}
 *
 * Properties:
 *   - Exactly zero sidelobes (all zeros degenerate at ψ=π)
 *   - Widest possible beamwidth for given N
 *   - Equivalent to maximally flat filter in signal processing
 *
 * Used primarily as a theoretical reference for the beamwidth-SLL trade-off.
 *
 * @param num_elements N ≥ 1.
 * @return Amplitude weights (caller frees).
 *
 * Reference: Balanis §6.8.2.
 */
double *pa_binomial_weights(uint32_t num_elements);

#ifdef __cplusplus
}
#endif

#endif /* PHASED_ARRAY_SYNTHESIS_H */
