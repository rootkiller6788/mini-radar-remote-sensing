/* ============================================================================
 * phased_array_steering.c — Beam Steering Implementation
 *
 * Implements:
 *   - Steering vector computation for arbitrary array geometries
 *   - Per-element phase shift calculation
 *   - Phase shifter quantization effects
 *   - True time delay (TTD) computation for wideband arrays
 *   - Beam squint analysis for phase-steered arrays
 *
 * Core concept (Van Trees §2.4):
 *   To steer a beam to direction (theta_s, phi_s), the complex weight
 *   for element n is:
 *     w_n = A_n * exp(j * alpha_n)
 *   where alpha_n = -k0 * (r_n · u_steer) cancels the spatial phase
 *   delay, producing constructive interference in the steer direction.
 *
 * References:
 *   Van Trees, H.L. (2002) "Optimum Array Processing", Part IV.
 *   Mailloux, R.J. (2005) "Phased Array Antenna Handbook" §3.1-3.4.
 *   Balanis, C.A. (2016) "Antenna Theory" §6.4.
 * ============================================================================ */

#include "phased_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#define C0 299792458.0

/* ============================================================================
 * L5: Steering Vector Computation
 *
 * The steering vector a(theta, phi) for an N-element array is:
 *   a_n(theta, phi) = exp(-j * k0 * (r_n · u(theta, phi)))
 *
 * where u = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
 * is the unit vector in the observation direction.
 *
 * Beamforming weights that point the beam to (theta_s, phi_s):
 *   w = a*(theta_s, phi_s)   [complex conjugate]
 *
 * This is the simplest beamformer: delay-and-sum (DAS).
 * The conjugate cancels the spatial propagation phase, aligning
 * all element contributions in the desired direction.
 *
 * For a uniform linear array along x-axis with spacing d:
 *   a_n = exp(-j * k0 * x_n * sin(theta) * cos(phi))
 *       = exp(-j * k0 * n*d * sin(theta))        [if phi=0, elements on x-axis]
 *
 * Phase difference between adjacent elements:
 *   Delta_phi = k0 * d * sin(theta) = 2*pi * d/lambda * sin(theta)
 *
 * ============================================================================ */

/**
 * Compute steering vector for the full array.
 *
 * Memory: caller must free() the returned array.
 * Length: config->num_elements complex doubles.
 */
double complex *pa_steering_vector(const pa_array_config_t *config,
                                    const pa_element_t *elements,
                                    double theta_steer, double phi_steer)
{
    if (!config || !elements || config->num_elements == 0)
        return NULL;

    uint32_t N = config->num_elements;
    double complex *sv = (double complex *)malloc(N * sizeof(double complex));
    if (!sv) return NULL;

    /* Unit direction vector for the steer direction */
    double u, v, w;
    pa_spherical_to_uvw(theta_steer, phi_steer, &u, &v, &w);

    double k0 = pa_wavenumber(config->frequency_hz);
    if (k0 <= 0.0) {
        /* DC case: return all-ones steering vector */
        for (uint32_t n = 0; n < N; n++)
            sv[n] = 1.0 + 0.0 * I;
        return sv;
    }

    for (uint32_t n = 0; n < N; n++) {
        /* Projection of element position onto the steer direction */
        double dot = elements[n].position.x * u
                   + elements[n].position.y * v
                   + elements[n].position.z * w;

        /* Phase = k0 * (r_n · u_steer)
         * Steering vector element a_n = exp(-j * phase)
         * The negative sign means we add phase to compensate for
         * the propagation delay, which cancels the spatial phase
         * when signals arrive from the desired direction. */
        double phase = -k0 * dot;
        sv[n] = cos(phase) + sin(phase) * I;
    }

    return sv;
}

/* ============================================================================
 * L5: Per-Element Phase Shift Computation
 *
 * For each element n, compute the phase shift phi_n needed to steer
 * the beam toward direction (theta_s, phi_s).
 *
 * Phase shift (radians):
 *   phi_n = k0 * (r_n · u_steer)
 *
 * For an N-element uniform linear array along x-axis:
 *   phi_n = k0 · n · d · sin(theta_s) · cos(phi_s)
 *
 * Phase slope (linear array):
 *   delta_phi/deg = 2*pi*d/lambda * cos(theta_s)  [rad per degree of scan]
 *
 * The phase shift is wrapped to [-pi, +pi] range, which is how
 * digital phase shifters operate (modulo 2*pi).
 * ============================================================================ */

/**
 * Compute per-element phase shifts for beam steering.
 *
 * @param config Array config
 * @param elements Array elements
 * @param theta_steer Target polar angle (rad)
 * @param phi_steer Target azimuth (rad)
 * @param phase_shifts Output array: size = num_elements (caller allocates)
 */
void pa_phase_shifts(const pa_array_config_t *config,
                     const pa_element_t *elements,
                     double theta_steer, double phi_steer,
                     double *phase_shifts)
{
    if (!config || !elements || !phase_shifts || config->num_elements == 0)
        return;

    uint32_t N = config->num_elements;
    double k0 = pa_wavenumber(config->frequency_hz);
    if (k0 <= 0.0) {
        for (uint32_t n = 0; n < N; n++)
            phase_shifts[n] = 0.0;
        return;
    }

    double u, v, w;
    pa_spherical_to_uvw(theta_steer, phi_steer, &u, &v, &w);

    for (uint32_t n = 0; n < N; n++) {
        double dot = elements[n].position.x * u
                   + elements[n].position.y * v
                   + elements[n].position.z * w;
        /* Phase shift = k0 · (r_n · u_steer), wrapped to [-pi, pi] */
        double phi_n = k0 * dot;

        /* Wrap to [-pi, +pi] for digital phase shifter compatibility */
        phi_n = fmod(phi_n + M_PI, 2.0 * M_PI);
        if (phi_n < 0.0) phi_n += 2.0 * M_PI;
        phi_n -= M_PI;
        phase_shifts[n] = phi_n;
    }
}

/* ============================================================================
 * L5: Phase Shifter Quantization
 *
 * Real T/R modules use digital phase shifters with finite resolution
 * (typically 3-6 bits). Quantization causes:
 *
 * 1. Beam pointing error: the beam peak shifts slightly because the
 *    ideal linear phase progression is quantized to discrete steps.
 *    RMS pointing error ≈ beamwidth / (sqrt(12) * 2^N).
 *
 * 2. Quantization lobes: periodic phase errors create periodic
 *    grating-lobe-like artifacts at angles where the phase error
 *    pattern constructively adds. RMS sidelobe level increase:
 *      Delta_SLL ≈ -5 - 6*N_bits  (dB, from Mailloux §3.4)
 *
 * 3. Gain loss: peak gain reduced by ~(1 - sinc²(π/2^N)).
 *    For 4-bit (N=4): gain loss ≈ 0.1 dB
 *    For 2-bit (N=2): gain loss ≈ 0.9 dB
 *
 * Phase step size:
 *   delta_phi = 2*pi / 2^{bits}
 *   For 4-bit: delta_phi = 22.5°
 *   For 6-bit: delta_phi = 5.625°
 *
 * Reference: Mailloux (2005) §3.4 "Quantization Effects".
 * ============================================================================ */

/**
 * Quantize phase shifts to finite-bit resolution.
 *
 * Rounds each phase to nearest available phase state.
 *
 * @param phase_shifts Input/output array (radians, size = N).
 * @param num_elements Number of elements.
 * @param bits Phase shifter resolution (1-8, or PA_PHASE_BITS_CONTINUOUS).
 */
void pa_quantize_phase_shifts(double *phase_shifts, uint32_t num_elements,
                              pa_phase_bits_t bits)
{
    if (!phase_shifts || num_elements == 0) return;

    /* Continuous phase (ideal analog phase shifter) — no quantization */
    if (bits == PA_PHASE_BITS_CONTINUOUS || bits == 0)
        return;

    uint32_t num_states = (uint32_t)(1u << (uint32_t)bits);
    double step = 2.0 * M_PI / (double)num_states;

    for (uint32_t n = 0; n < num_elements; n++) {
        /* Shift to [0, 2*pi) for quantization */
        double phi = phase_shifts[n];
        while (phi < 0.0) phi += 2.0 * M_PI;
        while (phi >= 2.0 * M_PI) phi -= 2.0 * M_PI;

        /* Round to nearest phase state:
         * state = round(phi / step)
         * quantized = state * step */
        double state = round(phi / step);
        if (state >= (double)num_states) state = 0.0;
        phi = state * step;

        /* Shift back to [-pi, pi] for consistency */
        if (phi > M_PI) phi -= 2.0 * M_PI;
        phase_shifts[n] = phi;
    }
}

/* ============================================================================
 * L6: True Time Delay (TTD) Beamforming
 *
 * Phase steering has a fundamental limitation: it works correctly only
 * at a single frequency. The phase shift phi_n = k0 * r_n · u_s is
 * frequency-dependent because k0 = 2*pi*f/c.
 *
 * TTD beamforming replaces phase shifters with variable time delays:
 *   tau_n = (r_n · u_steer) / c
 *
 * This produces a frequency-independent beam direction because:
 *   Phase at frequency f:  phi_n(f) = 2*pi*f * tau_n
 *                            = 2*pi*f * (r_n · u_s)/c
 *                            = k0 * r_n · u_s   [same as phase steering]
 *
 * The key difference: tau_n is constant across frequency, so the beam
 * direction remains correct over the entire instantaneous bandwidth.
 *
 * TTD is essential for:
 *   - Wideband radars (>10% fractional bandwidth)
 *   - Large arrays where aperture fill time > pulse width
 *   - High-resolution imaging with large instantaneous bandwidth
 *
 * TTD implementation: switched delay lines, photonic delay lines,
 * or (in digital beamforming) fractional-sample time delays.
 * ============================================================================ */

/**
 * Compute per-element true-time-delay values.
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Steer direction theta (rad).
 * @param phi_steer Steer direction phi (rad).
 * @param delays_ps Output: per-element delays in picoseconds.
 *                   Caller allocates size = num_elements.
 */
void pa_ttd_compute_delays(const pa_array_config_t *config,
                            const pa_element_t *elements,
                            double theta_steer, double phi_steer,
                            double *delays_ps)
{
    if (!config || !elements || !delays_ps || config->num_elements == 0)
        return;

    uint32_t N = config->num_elements;
    double u, v, w;
    pa_spherical_to_uvw(theta_steer, phi_steer, &u, &v, &w);

    /* tau_n = (r_n · u_steer) / c0
     * Convert to picoseconds: t_ps = t_seconds * 1e12 */
    for (uint32_t n = 0; n < N; n++) {
        double dot = elements[n].position.x * u
                   + elements[n].position.y * v
                   + elements[n].position.z * w;
        delays_ps[n] = (dot / C0) * 1.0e12;

        /* Reference delay to the element with the largest positive delay:
         * All delays should be non-negative for causal implementation.
         * Later, we subtract the minimum delay to make all >= 0. */
    }

    /* Find minimum delay and subtract to make all delays non-negative.
     * Real TTD systems can only delay, not advance in time. */
    double min_delay = delays_ps[0];
    for (uint32_t n = 1; n < N; n++) {
        if (delays_ps[n] < min_delay) min_delay = delays_ps[n];
    }
    for (uint32_t n = 0; n < N; n++) {
        delays_ps[n] -= min_delay;
    }
}

/* ============================================================================
 * L6: Beam Squint Analysis
 *
 * When a phased array uses phase shifters (not TTD), the beam direction
 * changes with frequency — this is "beam squint."
 *
 * For a linear array steered to angle theta_s at center frequency f0:
 *   sin(theta(f)) = (f0/f) * sin(theta_s)
 *
 * The beam squints toward broadside as frequency increases above f0,
 * and toward endfire as frequency decreases below f0.
 *
 * Beam squint angle:
 *   theta_squint = theta(f) - theta_s
 *
 * This limits the instantaneous bandwidth of a phased array:
 *   B_max / f0 = sin(theta_s) * Delta_theta_3dB   [Mailloux §5.2]
 *
 * For a large array with narrow beamwidth, even small frequency
 * offsets cause significant squint relative to the beamwidth.
 *
 * TTD solves this: tau_n is independent of frequency, so the beam
 * direction is the same at all frequencies within the bandwidth.
 *
 * Reference: Mailloux (2005) §5.2 "Bandwidth and Squint Effects".
 * ============================================================================ */

/**
 * Analyze beam squint for a phase-steered array.
 *
 * @param config Array config (frequency_hz is the center frequency).
 * @param theta_steer_deg Nominal steer angle (deg from broadside).
 * @param squint Output: beam squint analysis results.
 */
void pa_beam_squint_analysis(const pa_array_config_t *config,
                              double theta_steer_deg,
                              pa_beam_squint_t *squint)
{
    if (!config || !squint) return;

    memset(squint, 0, sizeof(*squint));
    squint->center_frequency_hz = config->frequency_hz;
    squint->steer_angle_nominal_deg = theta_steer_deg;

    /* Assume 10% bandwidth by default */
    squint->bandwidth_hz = config->frequency_hz * 0.10;

    double f0 = config->frequency_hz;
    if (f0 <= 0.0) return;

    double theta_s = theta_steer_deg * M_PI / 180.0;
    double sin_theta_s = sin(theta_s);

    /* Beam angle at the band edge (f_max = f0 + B/2):
     * sin(theta(f_max)) = (f0 / f_max) * sin(theta_s) */
    double f_upper = f0 + squint->bandwidth_hz / 2.0;
    double f_lower = f0 - squint->bandwidth_hz / 2.0;
    if (f_lower <= 0.0) f_lower = f0 * 0.01;

    /* Compute theta at upper and lower band edges */
    double sin_theta_upper = (f0 / f_upper) * sin_theta_s;
    double sin_theta_lower = (f0 / f_lower) * sin_theta_s;

    /* Clamp to [-1, 1] */
    if (sin_theta_upper > 1.0)  sin_theta_upper = 1.0;
    if (sin_theta_upper < -1.0) sin_theta_upper = -1.0;
    if (sin_theta_lower > 1.0)  sin_theta_lower = 1.0;
    if (sin_theta_lower < -1.0) sin_theta_lower = -1.0;

    double theta_upper = asin(sin_theta_upper);
    double theta_lower = asin(sin_theta_lower);

    /* Squint is the worst-case deviation */
    double squint_upper = (theta_upper - theta_s) * 180.0 / M_PI;
    double squint_lower = (theta_s - theta_lower) * 180.0 / M_PI;
    squint->squint_at_band_edge_deg =
        (fabs(squint_upper) > fabs(squint_lower)) ? squint_upper : squint_lower;

    /* Compute the maximum allowable bandwidth: when squint equals
     * half the 3dB beamwidth (commonly used criterion).
     * HPBW ≈ 0.886 * lambda / (N*d) => in degrees.
     * Then B_max/f0 = sin(theta_s) * HPBW (Mailloux Eq. 5.16) */
    double hp_deg;
    pa_element_t *elem = pa_allocate_elements(config);
    double complex *w = (double complex *)calloc(config->num_elements,
                                                  sizeof(double complex));
    if (w) {
        for (uint32_t i = 0; i < config->num_elements; i++)
            w[i] = 1.0 + 0.0 * I;
    }
    pa_beamwidth_3db(config, elem, w, theta_s, 0.0, &hp_deg);
    free(w);
    pa_free_elements(elem);

    if (hp_deg > 0.0) {
        double hp_rad = hp_deg * M_PI / 180.0;
        double f_bw = f0 * hp_rad / sin_theta_s;
        if (sin_theta_s < 0.01) f_bw = f0; /* No squint at broadside */
        squint->max_allowable_bandwidth_hz = f_bw;
        squint->fractional_bandwidth_limit = f_bw / f0;
    } else {
        squint->max_allowable_bandwidth_hz = squint->bandwidth_hz;
        squint->fractional_bandwidth_limit = 0.10;
    }
}

/* ============================================================================
 * L5: Amplitude Taper Generation
 *
 * Amplitude tapering reduces sidelobes at the expense of wider beamwidth
 * and reduced directivity.
 *
 * Window Type      | Peak SLL (dB) | HPBW broadening | Directivity loss
 * -----------------+---------------+-----------------+-----------------
 * Uniform           | -13.3         | 1.00x           | 0.0 dB
 * Hamming           | -42.6         | 1.46x           | 1.34 dB
 * Hann              | -31.5         | 1.62x           | 1.76 dB
 * Blackman          | -58.1         | 1.86x           | 2.36 dB
 * Dolph-Chebyshev   | user-defined  | varies          | varies
 * Taylor            | user-defined  | varies          | varies
 * Binomial          | -∞ (none)     | very wide       | large loss
 *
 * Reference: Harris (1978) "On the Use of Windows for Harmonic Analysis"
 *            Balanis §6.8 "Amplitude Taper Concepts"
 * ============================================================================ */

/**
 * Generate amplitude taper coefficients.
 *
 * @param num_elements Number of elements in the array.
 * @param window_type Type of amplitude taper.
 * @param window_param Window-specific parameter:
 *        Chebyshev: desired SLL in negative dB (e.g., -40.0)
 *        Kaiser: beta parameter (e.g., 3.0)
 *        Gaussian: sigma (e.g., 0.3)
 * @return double* Taper coefficients (caller frees).
 */
double *pa_amplitude_taper(uint32_t num_elements, pa_window_t window_type,
                           double window_param)
{
    if (num_elements == 0) return NULL;

    double *coeff = (double *)malloc(num_elements * sizeof(double));
    if (!coeff) return NULL;

    double N = (double)num_elements;

    switch (window_type) {
    case PA_WINDOW_UNIFORM:
        for (uint32_t n = 0; n < num_elements; n++)
            coeff[n] = 1.0;
        break;

    case PA_WINDOW_HAMMING:
        /* w[n] = 0.54 - 0.46 * cos(2*pi*n/(N-1)) */
        for (uint32_t n = 0; n < num_elements; n++)
            coeff[n] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)n / (N - 1.0));
        break;

    case PA_WINDOW_HANN:
        /* w[n] = 0.5 - 0.5 * cos(2*pi*n/(N-1)) */
        for (uint32_t n = 0; n < num_elements; n++)
            coeff[n] = 0.5 - 0.5 * cos(2.0 * M_PI * (double)n / (N - 1.0));
        break;

    case PA_WINDOW_BLACKMAN:
        /* w[n] = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1)) */
        for (uint32_t n = 0; n < num_elements; n++) {
            double a = 2.0 * M_PI * (double)n / (N - 1.0);
            coeff[n] = 0.42 - 0.5 * cos(a) + 0.08 * cos(2.0 * a);
        }
        break;

    case PA_WINDOW_GAUSSIAN: {
        /* w[n] = exp(-0.5 * ((n - (N-1)/2) / (sigma*(N-1)/2))^2) */
        double sigma = (window_param > 0.0) ? window_param : 0.4;
        double denom = sigma * (N - 1.0) / 2.0;
        for (uint32_t n = 0; n < num_elements; n++) {
            double x = ((double)n - (N - 1.0) / 2.0) / denom;
            coeff[n] = exp(-0.5 * x * x);
        }
        break;
    }

    case PA_WINDOW_KAISER: {
        /* Kaiser-Bessel window: w[n] = I0(beta*sqrt(1-(2n/(N-1)-1)^2))/I0(beta)
         * Approximate with sinh formulation for simplicity. */
        double beta = (window_param > 0.0) ? window_param : 3.0;
        for (uint32_t n = 0; n < num_elements; n++) {
            double x = (2.0 * (double)n / (N - 1.0)) - 1.0;
            double arg = beta * sqrt(1.0 - x * x);
            /* sinh(arg) approximation for I0 (zeroth-order modified Bessel) */
            coeff[n] = sinh(arg) / sinh(beta);
        }
        break;
    }

    case PA_WINDOW_RAISED_COSINE:
        /* w[n] = cos(pi * (n - (N-1)/2) / (N-1)) */
        for (uint32_t n = 0; n < num_elements; n++) {
            double x = M_PI * ((double)n - (N - 1.0) / 2.0) / (N - 1.0);
            coeff[n] = cos(x);
            if (coeff[n] < 0.0) coeff[n] = 0.0;
        }
        break;

    case PA_WINDOW_DOLPH_CHEBYSHEV: {
        double *w = pa_dolph_chebyshev_weights(num_elements, window_param);
        if (w) {
            memcpy(coeff, w, num_elements * sizeof(double));
            free(w);
        } else {
            for (uint32_t n = 0; n < num_elements; n++) coeff[n] = 1.0;
        }
        break;
    }

    case PA_WINDOW_TAYLOR: {
        double *w = pa_taylor_weights(num_elements, window_param, 5);
        if (w) {
            memcpy(coeff, w, num_elements * sizeof(double));
            free(w);
        } else {
            for (uint32_t n = 0; n < num_elements; n++) coeff[n] = 1.0;
        }
        break;
    }

    case PA_WINDOW_BINOMIAL: {
        double *w = pa_binomial_weights(num_elements);
        if (w) {
            memcpy(coeff, w, num_elements * sizeof(double));
            free(w);
        } else {
            for (uint32_t n = 0; n < num_elements; n++) coeff[n] = 1.0;
        }
        break;
    }

    default:
        /* Default to uniform */
        for (uint32_t n = 0; n < num_elements; n++)
            coeff[n] = 1.0;
        break;
    }

    /* Normalize so that the maximum coefficient is 1.0 (unity at center) */
    double max_val = 0.0;
    for (uint32_t n = 0; n < num_elements; n++)
        if (coeff[n] > max_val) max_val = coeff[n];
    if (max_val > 1e-15) {
        for (uint32_t n = 0; n < num_elements; n++)
            coeff[n] /= max_val;
    }

    return coeff;
}