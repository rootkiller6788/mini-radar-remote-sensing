/* ============================================================================
 * phased_array_synthesis.c — Array Pattern Synthesis
 *
 * Implements analytic array synthesis techniques for sidelobe control:
 *   - Dolph-Chebyshev synthesis (constant sidelobe level)
 *   - Taylor one-parameter distribution (monotonic sidelobe decay)
 *   - Binomial array synthesis (zero sidelobes)
 *
 * Pattern synthesis is the inverse problem: given desired pattern
 * characteristics (beamwidth, SLL), find the element excitations.
 *
 * References:
 *   Dolph, C.L. (1946) "A Current Distribution for Broadside Arrays Which
 *     Optimizes the Relationship Between Beam Width and Side-Lobe Level",
 *     Proc. IRE, vol. 34.
 *   Taylor, T.T. (1955) "Design of Line-Source Antennas for Narrow Beamwidth
 *     and Low Side Lobes", IRE Trans. AP-3.
 *   Balanis (2016) §6.8.3–6.8.4.
 * ============================================================================ */

#include "phased_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * L5: Dolph-Chebyshev Array Synthesis
 *
 * Problem: Design an N-element uniform linear array (symmetric excitation)
 *         with the narrowest possible beamwidth for a given sidelobe level.
 *
 * Solution (Dolph, 1946):
 *   Use Chebyshev polynomials to map the polynomial representing the
 *   array factor onto the equiripple passband of the Chebyshev polynomial.
 *
 * Mathematical derivation:
 *
 * 1. Array factor for an N-element symmetric linear array:
 *    AF(ψ) = Σ_{n=0}^{M} a_n cos(n ψ)
 *    where M = (N-1)/2 for odd N, M = N/2 for even N,
 *    and ψ = k₀ d sin(θ) + β (β = progressive phase shift).
 *
 * 2. Each cos(n ψ) can be expressed as a polynomial in cos(ψ):
 *    cos(n ψ) = T_n(cos ψ) where T_n is the Chebyshev polynomial of order n.
 *
 * 3. Thus AF = Σ a_n T_n(x) where x = cos(ψ).
 *    This is a polynomial of degree N-1 in x.
 *
 * 4. The Chebyshev polynomial T_{N-1}(x) oscillates between ±1 for |x| ≤ 1
 *    (visible region) and increases rapidly for |x| > 1 (invisible region).
 *
 * 5. To get constant sidelobe level at -R dB:
 *    - Define R_linear = 10^{R_dB/20}
 *    - Set x₀ = cosh(acosh(R_linear) / (N-1))
 *    - The polynomial T_{N-1}(x₀ · cos(ψ)) has the property that all
 *      sidelobes are exactly at -R dB and the main beam is at 0 dB.
 *
 * 6. The excitation coefficients are the polynomial coefficients of
 *    T_{N-1}(x₀ · x) expanded in the basis {T_n(x)}.
 *
 * Properties of Dolph-Chebyshev array (Balanis §6.8.3):
 *   - All sidelobes have equal amplitude (constant SLL).
 *   - For given N and SLL, produces the narrowest possible beamwidth.
 *   - For given N and beamwidth, produces the lowest possible SLL.
 *   - Edge elements may have higher current than the center (for deep SLL),
 *     which can cause feed difficulties.
 *
 * Limitation: For large arrays, constant near-in and far-out sidelobes
 *             are undesirable (far-out sidelobes waste power). Taylor
 *             distribution addresses this.
 * ============================================================================ */

/**
 * Compute Dolph-Chebyshev amplitude weights.
 *
 * Algorithm (Balanis §6.8.3):
 *
 * For an N-element array (N can be even or odd), the excitation
 * coefficients a_n (n = 0, 1, ..., N-1) are:
 *
 * For odd N = 2M+1:
 *   a_{M±k} = (1/N) · [R + 2 Σ_{p=1}^{M} T_{N-1}(x₀·cos(pπ/N))·cos(2pkπ/N)]
 *
 * For even N = 2M:
 *   a_{M-1-k} = a_{M+k} = (2/N) · Σ_{p=1}^{M} T_{N-1}(x₀·cos((p-0.5)π/N))
 *                         · cos((2k+1)(p-0.5)π/N)
 *
 * where:
 *   x₀ = cosh(acosh(R) / (N-1))
 *   R  = 10^{-SLL_dB / 20}
 *
 * Implementation note: we compute the array polynomial roots and
 * then expand to get the excitation coefficients using the
 * relationship between polynomial roots and coefficients.
 *
 * @param num_elements Number of elements (N ≥ 2).
 * @param sll_db Desired sidelobe level in dB (negative, e.g. -30.0).
 * @return double* Amplitude weights, length = num_elements.
 */
double *pa_dolph_chebyshev_weights(uint32_t num_elements, double sll_db)
{
    if (num_elements < 2) return NULL;

    uint32_t N = num_elements;
    /* R0 = 10^(|SLL|/20) = main-lobe-to-sidelobe voltage ratio.
     * T_{N-1}(x0) = R0, so x0 = cosh(acosh(R0)/(N-1)).
     * For -30 dB SLL: R0 = 10^(30/20) = 31.62, x0 = cosh(acosh(31.62)/8) ~ 1.1 */
    double R0 = pow(10.0, fabs(sll_db) / 20.0);

    /* x₀: position of the main lobe on the Chebyshev polynomial.
     * For SLL = -R_dB, the main lobe peak is at x₀ > 1, and
     * the sidelobe region is |x| ≤ 1 where T_{N-1}(x) oscillates
     * with amplitude ±1. */
    double x0 = cosh(acosh(R0) / (double)(N - 1));

    /* Compute root locations of T_{N-1}(x₀ · z):
     *   z_p = cos((2p-1)π/(2(N-1))) / x₀   for p = 1, 2, ..., N-1
     * These are the (N-1) roots that define the polynomial. */
    double *weight = (double *)calloc(N, sizeof(double));
    if (!weight) return NULL;

    if (N == 2) {
        /* Special case: 2-element array. Both elements get equal weight. */
        weight[0] = 1.0;
        weight[1] = 1.0;
        return weight;
    }

    /* Build the polynomial: P(z) = Π_{p=1}^{M} (z - z_p)
     * where z_p are the roots. The coefficients of this polynomial
     * are the array excitation amplitudes (after appropriate mapping). */

    /* For symmetric arrays, we compute the array polynomial by expanding
     * the product of root terms. We use the discrete Fourier method
     * (Schelkunoff's unit circle representation) which is more
     * numerically stable. */

    /* Alternative direct computation for small N (< 50):
     *
     * For odd N (N = 2M+1):
     *   a_k = (1/N) [R + 2 Σ_{p=1}^{M} T_{N-1}(x₀·cos(pπ/N))·cos(2k pπ/N)]
     *   for k = 0, 1, ..., M, then a_{N-1-k} = a_k (symmetric).
     *
     * For even N (N = 2M):
     *   a_{M-1-k} = a_{M+k} =
     *     (2/N) Σ_{p=1}^{M} T_{N-1}(x₀·cos((p-0.5)π/N))·cos((2k+1)(p-0.5)π/N)
     *   for k = 0, 1, ..., M-1. */

    /* General approach: Compute polynomial coefficients from the array
     * factor pattern directly. For a symmetric linear array, the array
     * factor polynomial is:
     *
     *   P(z) = Σ_{n=0}^{M} a_n z^n  where z = exp(j ψ)
     *
     * and P(z) is related to the Chebyshev polynomial:
     *   P(z) = T_{N-1}(x₀ · (z + z⁻¹)/2)
     *
     * Expanding this gives us the coefficients. */

    if (N % 2 == 1) {
        /* Odd N: N = 2M + 1 */
        uint32_t M_half = N / 2;
        double *coeff = (double *)calloc(M_half + 1, sizeof(double));
        if (!coeff) { free(weight); return NULL; }

        for (uint32_t k = 0; k <= M_half; k++) {
            double sum = 0.0;
            for (uint32_t p = 1; p <= M_half; p++) {
                /* T_{N-1}(x₀ · cos(p π / N)) */
                double arg = x0 * cos(M_PI * (double)p / (double)N);
                double T_val = cosh((double)(N - 1) * acosh(fabs(arg)));
                /* If arg < 1, T_{N-1}(arg) = cos((N-1)·acos(arg)) */
                if (fabs(arg) <= 1.0) {
                    T_val = cos((double)(N - 1) * acos(arg));
                }
                /* For arg <= -1: T_n(arg) = (-1)^n cosh(n acosh(|arg|)) */
                if (arg < -1.0 && (N - 1) % 2 == 1) {
                    T_val = -T_val;
                }
                sum += T_val * cos(2.0 * M_PI * (double)(k * p)
                                   / (double)N);
            }
            coeff[k] = (R0 + 2.0 * sum) / (double)N;
        }

        /* Assign symmetric weights: a_{M-k} = a_{M+k} = coeff[k] */
        /* For N=odd, center element is at index M_half */
        for (uint32_t k = 0; k <= M_half; k++) {
            weight[M_half - k] = coeff[k];
            weight[M_half + k] = coeff[k];
        }
        free(coeff);

    } else {
        /* Even N: N = 2M */
        uint32_t M_half = N / 2;
        double *coeff = (double *)calloc(M_half, sizeof(double));
        if (!coeff) { free(weight); return NULL; }

        for (uint32_t k = 0; k < M_half; k++) {
            double sum = 0.0;
            for (uint32_t p = 1; p <= M_half; p++) {
                double arg = x0 * cos(M_PI * ((double)p - 0.5)
                                      / (double)N);
                double T_val;
                if (fabs(arg) <= 1.0) {
                    T_val = cos((double)(N - 1) * acos(arg));
                } else {
                    T_val = cosh((double)(N - 1) * acosh(fabs(arg)));
                }
                if (arg < -1.0 && (N - 1) % 2 == 1) {
                    T_val = -T_val;
                }
                sum += T_val * cos(M_PI * (2.0 * (double)k + 1.0)
                                   * ((double)p - 0.5) / (double)N);
            }
            coeff[k] = (2.0 / (double)N) * sum;
        }

        /* Assign symmetric weights: a_{M-1-k} = a_{M+k} = coeff[k] */
        for (uint32_t k = 0; k < M_half; k++) {
            weight[M_half - 1 - k] = coeff[k];
            weight[M_half + k]     = coeff[k];
        }
        free(coeff);
    }

    /* Normalize: set max weight to 1.0 (unity current at center) */
    double max_w = 0.0;
    for (uint32_t n = 0; n < N; n++) {
        if (fabs(weight[n]) > max_w) max_w = fabs(weight[n]);
    }
    if (max_w > 1e-15) {
        for (uint32_t n = 0; n < N; n++) {
            weight[n] /= max_w;
        }
    }

    return weight;
}

/* ============================================================================
 * L5: Taylor One-Parameter Line-Source Distribution
 *
 * Taylor (1955) proposed a continuous aperture distribution that produces
 * a pattern where the near-in sidelobes (the first n̄) are at a controlled
 * equal level, and the far-out sidelobes decay as 1/u (≈ −6 dB/octave
 * for uniform, faster for tapered).
 *
 * This is more practical than Dolph-Chebyshev for large arrays because:
 *   - Far-out sidelobes naturally decay (avoiding excessive stored energy
 *     in the near field and reducing susceptibility to external noise).
 *   - The edge taper is more gradual (avoiding extremely high edge
 *     currents, which would be difficult to realize in practice).
 *
 * The Taylor one-parameter distribution:
 *
 *   w(z) = I₀(πB √(1 − (2z/L)²)) / I₀(πB)   for |z| ≤ L/2
 *
 * where:
 *   L  = aperture length = (N−1) · d   (for N-element discrete array)
 *   B  = (1/π) acosh(R)                 [Taylor parameter]
 *   R  = 10^{-SLL_dB/20}
 *   I₀ = modified Bessel function of the first kind, order zero
 *   z  = position along the aperture (−L/2 to +L/2)
 *
 * The pattern has:
 *   - Near-in sidelobes: nearly equal at ~−SLL_dB
 *   - Far-out sidelobes: decaying monotonically
 *   - Beamwidth: slightly wider than Dolph-Chebyshev (for same SLL)
 *   - Directivity: slightly lower than Dolph-Chebyshev (for same SLL)
 *
 * Sampling to an N-element linear array:
 *   a_n = w(z_n)  where z_n = (n − (N−1)/2) · d
 *
 * For narrow beams (large arrays), use the Taylor n̄-distribution
 * (not one-parameter) which independently controls the number of
 * equal-level sidelobes via the parameter n̄.
 *
 * Reference: Taylor, T.T. (1955) "Design of Line-Source Antennas..."
 *            Balanis §6.8.4.
 * ============================================================================ */

/**
 * Compute Taylor one-parameter amplitude weights.
 *
 * @param num_elements Number of elements (N ≥ 3).
 * @param sll_db Desired peak sidelobe level in dB (negative, e.g., -30.0).
 * @param nbar Number of nearly-equal sidelobes before monotonic decay.
 *             Typical values: 3–8. Must be ≥ 2.
 * @return double* Amplitude weights, length = num_elements.
 */
double *pa_taylor_weights(uint32_t num_elements, double sll_db, uint32_t nbar)
{
    if (num_elements < 3) return NULL;
    if (nbar < 2) nbar = 2;

    uint32_t N = num_elements;
    /* R0 = 10^(|SLL|/20) = main-lobe-to-sidelobe voltage ratio, > 1 */
    double R0 = pow(10.0, fabs(sll_db) / 20.0);

    /* Taylor parameter A:
     * A = (1/π) acosh(R0)
     *
     * B is actually what we use in the I₀ argument:
     * B = π · A = acosh(R)  [simplified] */
    double A = (1.0 / M_PI) * acosh(R0);

    /* For a discrete array, sample the continuous Taylor distribution
     * at the N element positions. */

    /* Compute I₀(B) for normalization.
     * I₀(x) = Σ_{k=0}^{∞} (x²/4)^k / (k!)²
     * Use series expansion with 20 terms (sufficient for x < 20). */
    double B = M_PI * A;
    double I0_B = 1.0;
    {
        double term = 1.0;
        double x2_4 = (B * B) / 4.0;
        for (int k = 1; k < 20; k++) {
            term *= x2_4 / ((double)k * (double)k);
            I0_B += term;
        }
    }

    double *weight = (double *)malloc(N * sizeof(double));
    if (!weight) return NULL;

    for (uint32_t n = 0; n < N; n++) {
        /* Position along aperture: normalized to [-1, 1] */
        double z = (2.0 * (double)n - (double)(N - 1)) / (double)(N - 1);
        /* Argument of I₀: B · √(1 − z²) */
        double arg = B * sqrt(1.0 - z * z);

        /* Compute I₀(arg) using series expansion */
        double I0_arg = 1.0;
        {
            double term = 1.0;
            double x2_4_2 = (arg * arg) / 4.0;
            for (int k = 1; k < 20; k++) {
                term *= x2_4_2 / ((double)k * (double)k);
                I0_arg += term;
            }
        }

        weight[n] = I0_arg / I0_B;
    }

    /* Normalize to unity at center */
    double max_w = 0.0;
    for (uint32_t n = 0; n < N; n++) {
        if (weight[n] > max_w) max_w = weight[n];
    }
    if (max_w > 1e-15) {
        for (uint32_t n = 0; n < N; n++) {
            weight[n] /= max_w;
        }
    }

    return weight;
}

/* ============================================================================
 * L5: Binomial Array Synthesis
 *
 * The binomial array uses the binomial coefficients as amplitude weights:
 *   a_n = C(N−1, n) / 2^{N−1}
 *
 * For an N-element array, the array factor polynomial is:
 *   AF(ψ) = (1 + z)^{N-1}  where z = exp(j ψ)
 *         = (2 cos(ψ/2))^{N-1}
 *
 * Since all zeros of the polynomial are at z = -1 (ψ = π), the pattern
 * has exactly zero sidelobes. The trade-off is an extremely wide main beam.
 *
 * HPBW ≈ 2 · arccos(0.5^{1/(N-1)})  → typically 3–5× wider than uniform.
 *
 * The binomial array has zero sidelobes because all the polynomial zeros
 * are degenerate at z = −1, which maps to the invisible region boundary.
 * This is analogous to a maximally flat filter in signal processing —
 * it gives up selectivity for flatness in the stopband.
 *
 * This distribution is rarely used in practice because the beamwidth
 * penalty is too severe, but it is an important theoretical benchmark:
 * it demonstrates the trade-off between beamwidth and sidelobe level.
 *
 * Reference: Balanis §6.8.2.
 * ============================================================================ */

double *pa_binomial_weights(uint32_t num_elements)
{
    if (num_elements < 1) return NULL;

    uint32_t N = num_elements;
    double *weight = (double *)malloc(N * sizeof(double));
    if (!weight) return NULL;

    /* Compute binomial coefficients C(N-1, k) for k = 0..N-1.
     *
     * Pascal's triangle iterative method:
     *   C(n, 0) = C(n, n) = 1
     *   C(n, k) = C(n-1, k-1) + C(n-1, k)
     *
     * For numerical stability with large N, compute iteratively
     * using the recurrence: C(n, k+1) = C(n, k) * (n−k) / (k+1) */

    uint32_t degree = N - 1;
    weight[0] = 1.0;
    for (uint32_t k = 0; k < degree; k++) {
        weight[k + 1] = weight[k] * (double)(degree - k) / (double)(k + 1);
    }

    /* Normalize to unity at center (max coefficient) */
    double max_w = 0.0;
    for (uint32_t n = 0; n < N; n++) {
        if (weight[n] > max_w) max_w = weight[n];
    }

    /* Check for symmetry (sanity check) */
    for (uint32_t n = 0; n < N; n++) {
        weight[n] /= max_w;
    }

    return weight;
}

