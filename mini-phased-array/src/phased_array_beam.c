/* ============================================================================
 * phased_array_beam.c — Beam Pattern Computation and Analysis
 *
 * Implements:
 *   - Full 2D beam pattern computation over angular range
 *   - Grating lobe analysis and visibility detection
 *   - Monopulse sum/difference pattern generation
 *   - Sidelobe level analysis
 *   - Pattern metrics extraction (peak, nulls, beamwidth)
 *
 * The beam pattern P(theta,phi) = |AF(theta,phi)|^2 is the fundamental
 * metric for array performance. It determines:
 *   - Angular resolution (beamwidth)
 *   - Interference rejection (sidelobe level)
 *   - Field of view (scan limits, grating lobe onset)
 *   - Target discrimination (sidelobe structure)
 *
 * References:
 *   Balanis (2016) §6.5-6.9
 *   Sherman & Barton (2011) "Monopulse Principles and Techniques"
 *   Skolnik (2008) "Radar Handbook" §8.7
 * ============================================================================ */

#include "phased_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#define DEG2RAD(d) ((d) * M_PI / 180.0)

/* ============================================================================
 * L6: Full 2D Beam Pattern Computation
 *
 * Computes the array factor magnitude pattern P(theta,phi) over a
 * specified angular grid. This is the most complete characterization
 * of array radiation behavior.
 *
 * Algorithm:
 *   For each (theta, phi) in the scan grid:
 *     1. Compute unit direction vector u = (sinθ cosφ, sinθ sinφ, cosθ)
 *     2. Sum contributions: AF = Σ w_n exp(j k0 r_n · u)
 *     3. Store |AF| in dB (normalized)
 *
 * After the scan:
 *     4. Find peak (main lobe)
 *     5. Find highest sidelobe
 *     6. Find 3dB points → compute beamwidth
 *     7. Find first null → compute FNBW
 *
 * Complexity: O(N · K_theta · K_phi) where N = number of elements
 *
 * Memory: pattern->af_db must be pre-allocated (K_theta × K_phi doubles).
 * ============================================================================ */

void pa_compute_pattern(const pa_array_config_t *config,
                        const pa_element_t *elements,
                        const double complex *weights,
                        pa_pattern_t *pattern)
{
    if (!config || !elements || !weights || !pattern) return;
    if (pattern->num_theta == 0 || pattern->num_phi == 0) return;

    uint32_t Nt = pattern->num_theta;
    uint32_t Np = pattern->num_phi;
    double th_start = pattern->theta_start;
    double th_end   = pattern->theta_end;
    double ph_start = pattern->phi_start;
    double ph_end   = pattern->phi_end;

    double th_step = (Nt > 1) ? (th_end - th_start) / (double)(Nt - 1) : 0.0;
    double ph_step = (Np > 1) ? (ph_end - ph_start) / (double)(Np - 1) : 0.0;

    double max_af_mag = -1.0;
    uint32_t peak_ti = 0, peak_pi = 0;
    double sum_sq_af = 0.0;  /* For directivity integration */

    /* Phase 1: Scan all angles, compute AF, track peak */
    for (uint32_t ti = 0; ti < Nt; ti++) {
        double theta = th_start + (double)ti * th_step;
        for (uint32_t pi = 0; pi < Np; pi++) {
            double phi = ph_start + (double)pi * ph_step;

            pa_af_result_t r = pa_array_factor(config, elements, weights,
                                                theta, phi);
            double mag = r.af_magnitude;

            size_t idx = (size_t)ti * Np + pi;
            pattern->af_db[idx] = mag;

            /* Track peak */
            if (mag > max_af_mag) {
                max_af_mag = mag;
                peak_ti = ti;
                peak_pi = pi;
            }
            sum_sq_af += mag * mag;
        }
    }

    /* Phase 2: Normalize to peak = 0 dB, compute metrics */
    double max_sll_mag = 0.0; /* Highest sidelobe magnitude */
    double peak_theta = th_start + (double)peak_ti * th_step;

    /* Simple main lobe detection: search outward from peak until first null,
     * then everything beyond is a sidelobe. This uses theta direction
     * (at the peak phi cut). */

    /* Find first null to the right of the peak (in theta direction) */
    uint32_t pi_peak = peak_pi;
    int null_found_right = 0;
    double fnbw_right = peak_theta;
    for (uint32_t ti = peak_ti + 1; ti < Nt; ti++) {
        size_t idx = (size_t)ti * Np + pi_peak;
        /* A null: magnitude drops and then rises */
        if (ti + 1 < Nt) {
            size_t idx_next = (size_t)(ti + 1) * Np + pi_peak;
            if (pattern->af_db[idx] < pattern->af_db[idx_next]
                && pattern->af_db[idx] < max_af_mag * 0.01) {
                fnbw_right = th_start + (double)ti * th_step;
                null_found_right = 1;
                break;
            }
        }
    }

    double fnbw_left = peak_theta;
    int null_found_left = 0;
    if (peak_ti > 0) {
        for (int32_t ti = (int32_t)peak_ti - 1; ti >= 0; ti--) {
            size_t idx = (size_t)ti * Np + pi_peak;
            if (ti > 0) {
                size_t idx_prev = (size_t)(ti - 1) * Np + pi_peak;
                if (pattern->af_db[idx] < pattern->af_db[idx_prev]
                    && pattern->af_db[idx] < max_af_mag * 0.01) {
                    fnbw_left = th_start + (double)ti * th_step;
                    null_found_left = 1;
                    break;
                }
            }
        }
    }

    /* Compute FNBW */
    if (null_found_right && null_found_left) {
        double fnbw_rad = fnbw_right - fnbw_left;
        pattern->first_null_beamwidth_deg = fabs(fnbw_rad) * 180.0 / M_PI;
    }

    /* Find 3dB beamwidth: search for half-power points (-3dB = 0.707*peak) */
    double hp_level = max_af_mag / sqrt(2.0);  /* -3dB in linear */
    double hp_right = peak_theta, hp_left = peak_theta;

    /* Search right from peak */
    for (uint32_t ti = peak_ti + 1; ti < Nt; ti++) {
        size_t idx = (size_t)ti * Np + pi_peak;
        size_t idx_prev = (size_t)(ti - 1) * Np + pi_peak;
        if (pattern->af_db[idx] <= hp_level && pattern->af_db[idx_prev] >= hp_level) {
            /* Linear interpolation for more accurate -3dB point */
            double mag_prev = pattern->af_db[idx_prev];
            double mag_curr = pattern->af_db[idx];
            double frac = (hp_level - mag_prev) / (mag_curr - mag_prev);
            hp_right = th_start + ((double)(ti - 1) + frac) * th_step;
            break;
        }
    }

    /* Search left from peak */
    if (peak_ti > 0) {
        for (int32_t ti = (int32_t)peak_ti - 1; ti >= 0; ti--) {
            size_t idx = (size_t)ti * Np + pi_peak;
            size_t idx_next = (size_t)(ti + 1) * Np + pi_peak;
            if (pattern->af_db[idx] <= hp_level && pattern->af_db[idx_next] >= hp_level) {
                double mag_prev = pattern->af_db[idx];
                double mag_curr = pattern->af_db[idx_next];
                double frac = (hp_level - mag_prev) / (mag_curr - mag_prev);
                hp_left = th_start + ((double)ti + frac) * th_step;
                break;
            }
        }
    }

    double hpbw_rad = hp_right - hp_left;
    pattern->half_power_beamwidth_deg = fabs(hpbw_rad) * 180.0 / M_PI;

    /* Find worst-case sidelobe level */
    /* Everything in the pattern not in the main lobe region (between first
     * nulls or ±2*HPBW from peak) is considered sidelobe. */
    double lobe_half_width = pattern->half_power_beamwidth_deg * M_PI / 180.0 * 2.0;
    if (null_found_right && null_found_left) {
        lobe_half_width = fabs(fnbw_right - peak_theta);
    }
    double lobe_lo = peak_theta - lobe_half_width;
    double lobe_hi = peak_theta + lobe_half_width;

    for (uint32_t ti = 0; ti < Nt; ti++) {
        double theta = th_start + (double)ti * th_step;
        /* Skip main lobe region */
        if (theta > lobe_lo && theta < lobe_hi) continue;
        for (uint32_t pi = 0; pi < Np; pi++) {
            size_t idx = (size_t)ti * Np + pi;
            if (pattern->af_db[idx] > max_sll_mag) {
                max_sll_mag = pattern->af_db[idx];
            }
        }
    }

    /* SLL in dB relative to peak */
    if (max_sll_mag > 1e-15 && max_af_mag > 1e-15) {
        pattern->sidelobe_level_db = 20.0 * log10(max_sll_mag / max_af_mag);
    } else {
        pattern->sidelobe_level_db = -300.0;
    }

    /* Normalize all values to peak = 0 dB */
    for (size_t i = 0; i < (size_t)Nt * Np; i++) {
        if (pattern->af_db[i] < 1e-15) {
            pattern->af_db[i] = -300.0;
        } else {
            pattern->af_db[i] = 20.0 * log10(pattern->af_db[i] / max_af_mag);
        }
    }

    /* Peak directivity estimate */
    pattern->max_directivity_dbi = 10.0 * log10(
        pa_directivity(config, elements, weights));
}

/* ============================================================================
 * L6: Grating Lobe Analysis
 *
 * Grating lobes are replicas of the main beam that occur when the
 * element spacing exceeds a critical threshold. They are a spatial
 * aliasing phenomenon analogous to aliasing in time-domain sampling.
 *
 * Condition for grating lobe visibility (linear array, Balanis §6.9):
 *
 *   sin(theta_g) - sin(theta_s) = ± m * lambda/d   (m = ±1, ±2, ...)
 *
 * For a grating lobe to exist in visible space: |sin(theta_g)| ≤ 1.
 *
 * Maximum spacing to avoid grating lobes:
 *   d_max = lambda / (1 + |sin(theta_s_max)|)
 *
 * For broadside (theta_s = 0):   d < lambda      to avoid all grating lobes
 * For scan to ±60° (theta_s = 60°): d < 0.536*lambda
 * For ±30° scan: d < 0.667*lambda
 *
 * Standard practical choice: d = 0.5*lambda (λ/2), which allows
 * scan to ±90° without grating lobes, but with some scan loss.
 *
 * The element pattern suppresses grating lobes that fall near
 * or beyond the scan limits, which is why practical arrays can
 * sometimes use d > λ/2 if the element pattern is narrow enough.
 *
 * Reference: Balanis (2016) §6.9 "Grating Lobes".
 * ============================================================================ */

int pa_find_grating_lobes(const pa_array_config_t *config,
                          const pa_element_t *elements,
                          const double complex *weights,
                          double theta_steer, double phi_steer,
                          uint32_t max_lobes, pa_grating_lobe_t *lobes)
{
    if (!config || !elements || !weights || !lobes || max_lobes == 0)
        return 0;

    double lambda = pa_wavelength(config->frequency_hz);
    if (lambda <= 0.0) return 0;

    int count = 0;

    switch (config->geometry) {

    case PA_GEOMETRY_LINEAR: {
        /* For linear array along x-axis, grating lobe condition:
         * sin(theta_g) = sin(theta_s) + m*lambda/d   for m = ±1, ±2, ...
         *
         * Also at phi = phi_steer (in the scan plane). */
        double d = config->element_spacing_x;
        double sin_theta_s = sin(theta_steer);

        /* Search for grating lobes at order m = ±1, ±2, ... */
        for (int m = -5; m <= 5; m++) {
            if (m == 0) continue;  /* m=0 is the main lobe */
            if (count >= (int)max_lobes) break;

            double sin_theta_g = sin_theta_s + (double)m * lambda / d;

            /* Check if visible: |sin(theta_g)| <= 1 */
            if (fabs(sin_theta_g) <= 1.0) {
                double theta_g = asin(sin_theta_g);

                /* Skip if theta is outside [0, pi] */
                if (theta_g < 0.0 || theta_g > M_PI) continue;

                /* Compute the array factor at this angle to get
                 * the relative level of the grating lobe */
                pa_af_result_t r = pa_array_factor(config, elements, weights,
                                                    theta_g, phi_steer);
                pa_af_result_t r_peak = pa_array_factor(config, elements,
                                                         weights,
                                                         theta_steer, phi_steer);

                lobes[count].theta_grating_rad  = theta_g;
                lobes[count].phi_grating_rad    = phi_steer;
                lobes[count].order_m            = m;

                if (r_peak.af_magnitude > 1e-15) {
                    lobes[count].relative_level_db =
                        20.0 * log10(r.af_magnitude / r_peak.af_magnitude);
                } else {
                    lobes[count].relative_level_db = -300.0;
                }
                count++;
            }
        }
        break;
    }

    case PA_GEOMETRY_PLANAR_RECT: {
        /* Planar array has grating lobes in both dimensions.
         * Condition: sinθ cosφ = sinθ_s cosφ_s + m*lambda/dx
         *            sinθ sinφ = sinθ_s sinφ_s + n*lambda/dy */
        double dx = config->element_spacing_x;
        double dy = config->element_spacing_y;
        double u_s = sin(theta_steer) * cos(phi_steer);
        double v_s = sin(theta_steer) * sin(phi_steer);

        for (int m = -3; m <= 3; m++) {
            for (int n = -3; n <= 3; n++) {
                if (m == 0 && n == 0) continue;
                if (count >= (int)max_lobes) break;

                double u_g = u_s + (double)m * lambda / dx;
                double v_g = v_s + (double)n * lambda / dy;

                /* Check visibility: u² + v² ≤ 1 */
                double r2 = u_g*u_g + v_g*v_g;
                if (r2 <= 1.0) {
                    double theta_g = asin(sqrt(r2));
                    double phi_g = atan2(v_g, u_g);
                    if (phi_g < 0.0) phi_g += 2.0 * M_PI;

                    pa_af_result_t r = pa_array_factor(config, elements,
                                                        weights,
                                                        theta_g, phi_g);
                    pa_af_result_t r_peak = pa_array_factor(config, elements,
                                                             weights,
                                                             theta_steer,
                                                             phi_steer);

                    lobes[count].theta_grating_rad  = theta_g;
                    lobes[count].phi_grating_rad    = phi_g;
                    lobes[count].order_m            = m * 10 + n;

                    if (r_peak.af_magnitude > 1e-15) {
                        lobes[count].relative_level_db =
                            20.0 * log10(r.af_magnitude / r_peak.af_magnitude);
                    } else {
                        lobes[count].relative_level_db = -300.0;
                    }
                    count++;
                }
            }
        }
        break;
    }

    default:
        /* Grating lobe analysis for other geometries not implemented.
         * General principle: check spatial aliasing condition for
         * the specific lattice type. */
        break;
    }

    return count;
}

/* ============================================================================
 * L6: Monopulse Angle Estimation
 *
 * Monopulse is a radar technique for estimating target angle from a
 * single pulse (hence "mono-pulse"). It compares signals from two
 * overlapping beams to derive angular error.
 *
 * Principle (amplitude comparison monopulse):
 *   1. Generate two squinted beams: Beam 1 at θ₀ − Δθ/2,
 *                                   Beam 2 at θ₀ + Δθ/2
 *   2. Form Sum (Σ) = Beam1 + Beam2  — used for detection
 *   3. Form Difference (Δ) = Beam1 − Beam2 — used for angle error
 *   4. Error signal ε = Re{Δ/Σ} is zero at boresight and approximately
 *      linear in a small region around it.
 *
 * For a linear array (phase comparison):
 *   Split the array into left (indices 0..N/2-1) and right (N/2..N-1) halves.
 *   Σ(θ) = AF_left(θ) + AF_right(θ)
 *   Δ(θ) = AF_left(θ) − AF_right(θ)
 *
 * The monopulse slope k_m = d/dθ (Δ/Σ) at θ = θ₀ determines the
 * sensitivity of the angle estimate:
 *   θ_est = θ₀ + ε/k_m
 *
 * Monopulse accuracy: σ_θ ≈ θ_beamwidth / sqrt(2 · SNR)
 * (the Cramer-Rao bound for angle estimation).
 *
 * Applications: fire-control radar (F-35 AN/APG-81), missile seekers,
 * air traffic control, and satellite tracking.
 *
 * Reference: Sherman & Barton (2011) "Monopulse Principles and Techniques".
 *            Skolnik (2008) §8.7.
 * ============================================================================ */

void pa_monopulse_patterns(const pa_array_config_t *config,
                           const pa_element_t *elements,
                           const double complex *weights,
                           const double *theta_vals, uint32_t num_theta,
                           double *sum_pattern, double *diff_pattern,
                           double *error_signal)
{
    if (!config || !elements || !weights || !theta_vals
        || !sum_pattern || !diff_pattern || !error_signal)
        return;

    uint32_t N = config->num_elements;
    if (N < 2) return;

    /* Split the array into two halves.
     * For even N: halves are N/2 each.
     * For odd N: center element is shared (contributes to both halves). */
    uint32_t half = N / 2;

    /* Compute sum and difference patterns for each theta */
    for (uint32_t k = 0; k < num_theta; k++) {
        double theta = theta_vals[k];

        /* AF for left half */
        double complex af_left = 0.0 + 0.0 * I;
        {
            double u, v, w;
            pa_spherical_to_uvw(theta, 0.0, &u, &v, &w);
            double k0 = pa_wavenumber(config->frequency_hz);
            for (uint32_t n = 0; n < half; n++) {
                double dot = elements[n].position.x * u
                           + elements[n].position.y * v
                           + elements[n].position.z * w;
                double phase = k0 * dot;
                af_left += weights[n] * (cos(phase) + sin(phase) * I);
            }
            /* For odd N, share center element */
            if (N % 2 == 1) {
                uint32_t nc = N / 2;
                double dot = elements[nc].position.x * u
                           + elements[nc].position.y * v
                           + elements[nc].position.z * w;
                double phase = k0 * dot;
                af_left += 0.5 * weights[nc] * (cos(phase) + sin(phase) * I);
            }
        }

        /* AF for right half */
        double complex af_right = 0.0 + 0.0 * I;
        {
            double u, v, w;
            pa_spherical_to_uvw(theta, 0.0, &u, &v, &w);
            double k0 = pa_wavenumber(config->frequency_hz);
            uint32_t start = (N % 2 == 1) ? half + 1 : half;
            for (uint32_t n = start; n < N; n++) {
                double dot = elements[n].position.x * u
                           + elements[n].position.y * v
                           + elements[n].position.z * w;
                double phase = k0 * dot;
                af_right += weights[n] * (cos(phase) + sin(phase) * I);
            }
            /* For odd N, share center element */
            if (N % 2 == 1) {
                uint32_t nc = N / 2;
                double dot = elements[nc].position.x * u
                           + elements[nc].position.y * v
                           + elements[nc].position.z * w;
                double phase = k0 * dot;
                af_right += 0.5 * weights[nc] * (cos(phase) + sin(phase) * I);
            }
        }

        double complex af_sum  = af_left + af_right;
        double complex af_diff = af_left - af_right;

        sum_pattern[k]  = cabs(af_sum);
        diff_pattern[k] = cabs(af_diff);

        /* Error signal: ε = Re{Δ/Σ}
         * For small angles: ε ≈ k_m · (θ − θ₀), proportional to angle offset */
        double mag_sum = cabs(af_sum);
        if (mag_sum > 1e-15) {
            /* Real part of Δ/Σ = (Re{Δ}·Re{Σ} + Im{Δ}·Im{Σ}) / |Σ|² */
            error_signal[k] = (creal(af_diff) * creal(af_sum)
                               + cimag(af_diff) * cimag(af_sum))
                              / (mag_sum * mag_sum);
        } else {
            error_signal[k] = 0.0;
        }
    }
}

/**
 * Estimate DOA angle from monopulse error signal.
 *
 * θ_est = θ_boresight + ε / k_m
 *
 * where k_m = dε/dθ|_{θ=θ₀} is the monopulse slope, typically on the
 * order of 1-2 rad⁻¹ (i.e., ε ≈ 0.02–0.04 per degree).
 *
 * The slope can be calibrated by measuring ε at a known offset angle:
 *   k_m = ε(θ_offset) / (θ_offset − θ_boresight)
 *
 * @param error_signal Measured ε = Re{Δ/Σ}.
 * @param slope_factor Monopulse slope k_m (rad⁻¹).
 * @param boresight_rad Boresight angle (rad).
 * @return Estimated angle in radians.
 */
double pa_monopulse_estimate_angle(double error_signal,
                                    double slope_factor,
                                    double boresight_rad)
{
    if (fabs(slope_factor) < 1e-15) return boresight_rad;
    return boresight_rad + error_signal / slope_factor;
}

/* ============================================================================
 * L5: Mutual Coupling Model
 *
 * Mutual coupling is the electromagnetic interaction between array elements.
 * It distorts the element patterns, changes the input impedance, and degrades
 * adaptive beamformer performance if uncorrected.
 *
 * Coupling model using S-parameters:
 *   For an N-element array, the S-matrix S (N×N) relates incident waves a
 *   to reflected waves b:  b = S · a
 *
 * The coupling-corrected excitation is:
 *   w_corrected = C^{-1} · w_ideal
 * where C = I − S (the mutual coupling matrix).
 *
 * Simplified physical model for uniform linear array:
 *   S_{mn} = s₀ · exp(−α · |m−n|) · exp(−j k₀ d |m−n|)
 *
 * where:
 *   s₀    = nearest-neighbor coupling coefficient (~0.1–0.3)
 *   α     = exponential decay rate
 *   k₀ d  = electrical spacing between adjacent elements
 *
 * Reference: Gupta & Ksienski (1983) IEEE Trans. AP-31.
 * ============================================================================ */

void pa_mutual_coupling_correct(uint32_t num_elements,
                                 const double complex *s_matrix,
                                 double complex *weights)
{
    if (!s_matrix || !weights || num_elements < 2) return;

    /* The mutual coupling correction inverts the coupling matrix.
     * For weak coupling (|S_{mn}| << 1 for m≠n), approximate:
     *   w_corrected ≈ (I + S) · w  (first-order Neumann series)
     *
     * Full inversion would require O(N³) LAPACK — beyond this library.
     * For moderate coupling, the first-order correction is sufficient.
     * We apply: w_corrected[n] = w[n] + Σ_{m≠n} S_{nm} · w[m] */

    double complex *w_copy = (double complex *)malloc(num_elements
                                                      * sizeof(double complex));
    if (!w_copy) return;

    memcpy(w_copy, weights, num_elements * sizeof(double complex));

    for (uint32_t n = 0; n < num_elements; n++) {
        double complex correction = 0.0 + 0.0 * I;
        for (uint32_t m = 0; m < num_elements; m++) {
            if (m == n) continue;
            /* S_{nm} is the coupled signal from element m into element n */
            size_t idx = (size_t)n * num_elements + m;
            correction += s_matrix[idx] * w_copy[m];
        }
        /* Subtract the coupled contribution (cancellation)
         * Note: sign convention depends on the reference plane definition */
        weights[n] -= correction;
    }

    free(w_copy);
}

/**
 * Generate mutual coupling S-matrix for a uniform linear array.
 *
 * Uses an exponential decay model for coupling between identical elements.
 *
 * @param num_elements Number of elements.
 * @param spacing_m Element spacing (m).
 * @param frequency_hz Operating frequency (Hz).
 * @param coupling_coeff |S_{n,n+1}| nearest-neighbor coupling magnitude.
 * @param decay_rate α (per-element decay exponent).
 * @param s_matrix Output: N×N complex S-matrix (row-major, caller allocates).
 */
void pa_compute_coupling_matrix(uint32_t num_elements,
                                 double spacing_m, double frequency_hz,
                                 double coupling_coefficient,
                                 double decay_rate,
                                 double complex *s_matrix)
{
    if (!s_matrix || num_elements == 0) return;

    double k0 = 2.0 * M_PI * frequency_hz / 299792458.0;

    for (uint32_t m = 0; m < num_elements; m++) {
        for (uint32_t n = 0; n < num_elements; n++) {
            size_t idx = (size_t)m * num_elements + n;

            if (m == n) {
                /* Self-coupling: S_{nn} = 0 (matched element) */
                s_matrix[idx] = 0.0 + 0.0 * I;
            } else {
                int32_t dist = (int32_t)m - (int32_t)n;
                double abs_dist = (double)(dist > 0 ? dist : -dist);
                /* Exponential decay with electrical distance */
                double s_mag = coupling_coefficient
                               * exp(-decay_rate * abs_dist);
                /* Phase shift from propagation between elements */
                double s_phase = -k0 * spacing_m * abs_dist;
                s_matrix[idx] = s_mag * (cos(s_phase) + sin(s_phase) * I);
            }
        }
    }
}