#ifndef PHASED_ARRAY_BEAMFORMING_H
#define PHASED_ARRAY_BEAMFORMING_H

#include "phased_array.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L5: Beam Steering — Phase and Time-Delay Beamforming
 *
 * Core beamforming functions for directing the array beam.
 * Covers both narrowband (phase) and wideband (TTD) steering.
 *
 * Reference: Van Trees (2002) §2.4; Mailloux (2005) §5.1-5.2
 * ============================================================================ */

/**
 * Compute steering vector a(θ_s, φ_s) for the look direction.
 *
 * a_n = exp(-j k₀ (r_n · û_steer))
 *
 * Applying w = a* as excitation weights produces a beam in direction
 * (θ_s, φ_s) — this is delay-and-sum (DAS) beamforming.
 *
 * @param config Array geometry.
 * @param elements Array elements (positions used).
 * @param theta_steer Target θ (rad).
 * @param phi_steer Target φ (rad).
 * @return Complex steering vector, length = num_elements. Caller frees.
 */
double complex *pa_steering_vector(const pa_array_config_t *config,
                                    const pa_element_t *elements,
                                    double theta_steer, double phi_steer);

/**
 * Compute per-element phase shifts φ_n [−π, +π] for beam steering.
 *
 * φ_n = k₀ · (r_n · û_steer)  (wrapped to [−π, +π])
 *
 * Phase slope (linear array): Δφ = k₀ · d · sin(θ_s) per element.
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Target θ (rad).
 * @param phi_steer Target φ (rad).
 * @param phase_shifts Output: per-element phases in radians.
 *                      Caller allocates size = num_elements.
 */
void pa_phase_shifts(const pa_array_config_t *config,
                     const pa_element_t *elements,
                     double theta_steer, double phi_steer,
                     double *phase_shifts);

/**
 * Quantize phase shifts to finite-bit resolution.
 *
 * Models the effect of digital phase shifters in T/R modules.
 *
 * Quantization effects:
 *   1. RMS beam pointing error ≈ BW / (√12 · 2^N)
 *   2. Quantization lobes: ΔSLL ≈ −5 − 6·N_bits (dB)
 *   3. Gain loss ≈ 1 − sinc²(π/2^N)
 *
 * @param phase_shifts I/O: unquantized → quantized (rad).
 * @param num_elements Number of elements.
 * @param bits Phase shifter resolution (e.g., PA_PHASE_BITS_4 for 22.5° steps).
 *
 * Reference: Mailloux (2005) §3.4.
 */
void pa_quantize_phase_shifts(double *phase_shifts, uint32_t num_elements,
                              pa_phase_bits_t bits);

/**
 * Compute per-element true-time-delay values for wideband beamforming.
 *
 * τ_n = (r_n · û_steer) / c₀   [seconds]
 *
 * TTD eliminates beam squint — the delay is independent of frequency,
 * so the beam direction remains correct across the full bandwidth.
 *
 * Delays are referenced to the minimum delay element for causal implementation
 * (all delays ≥ 0).
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Steer θ (rad).
 * @param phi_steer Steer φ (rad).
 * @param delays_ps Output: per-element delays in picoseconds.
 *                  Caller allocates size = num_elements.
 */
void pa_ttd_compute_delays(const pa_array_config_t *config,
                            const pa_element_t *elements,
                            double theta_steer, double phi_steer,
                            double *delays_ps);

/**
 * Analyze beam squint for phase-steered arrays.
 *
 * For a phase-steered array at f₀, the beam direction at frequency f is:
 *   sin(θ(f)) = (f₀/f) · sin(θ_s)
 *
 * Beam squint = θ(f) − θ_s at band edges.
 *
 * @param config Array config (center frequency).
 * @param theta_steer_deg Nominal steer angle (degrees from broadside).
 * @param squint Output: squint analysis results.
 */
void pa_beam_squint_analysis(const pa_array_config_t *config,
                              double theta_steer_deg,
                              pa_beam_squint_t *squint);

/**
 * Generate amplitude taper coefficients for sidelobe control.
 *
 * Supported windows:
 *   Uniform, Hamming, Hann, Blackman, Gaussian, Kaiser,
 *   Raised Cosine, Dolph-Chebyshev, Taylor, Binomial
 *
 * @param num_elements Number of elements.
 * @param window_type Taper type.
 * @param window_param Window-specific parameter (SLL for Chebyshev/Taylor,
 *                      beta for Kaiser, sigma for Gaussian).
 * @return Taper coefficients (caller frees). NULL on error.
 */
double *pa_amplitude_taper(uint32_t num_elements, pa_window_t window_type,
                           double window_param);

#ifdef __cplusplus
}
#endif

#endif /* PHASED_ARRAY_BEAMFORMING_H */
