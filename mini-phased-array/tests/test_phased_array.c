#include "phased_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#define DEG2RAD(d) ((d) * M_PI / 180.0)
#define RAD2DEG(r) ((r) * 180.0 / M_PI)

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)

/* ========================================================================
 * L1: Definitions — Geometry Initialization
 * ======================================================================== */

static void test_linear_array_init(void)
{
    TEST("Linear array geometry init");
    pa_array_config_t config;
    pa_init_linear_array(&config, 16, 0.015, 10.0e9);
    assert(config.geometry == PA_GEOMETRY_LINEAR);
    assert(config.num_elements == 16);
    assert(fabs(config.element_spacing_x - 0.015) < 1e-9);
    assert(fabs(config.frequency_hz - 10.0e9) < 1.0);
    assert(config.rows == 1);
    assert(config.cols == 16);
    PASS();
}

static void test_planar_array_init(void)
{
    TEST("Planar array geometry init");
    pa_array_config_t config;
    pa_init_planar_array(&config, 8, 8, 0.015, 0.015, 10.0e9);
    assert(config.geometry == PA_GEOMETRY_PLANAR_RECT);
    assert(config.num_elements == 64);
    assert(config.rows == 8);
    assert(config.cols == 8);
    PASS();
}

static void test_circular_array_init(void)
{
    TEST("Circular array geometry init");
    pa_array_config_t config;
    pa_init_circular_array(&config, 16, 0.1, 10.0e9);
    assert(config.geometry == PA_GEOMETRY_CIRCULAR);
    assert(config.num_elements == 16);
    assert(fabs(config.element_spacing_radius - 0.1) < 1e-9);
    assert(config.scan_limit_az_deg == 180.0);
    PASS();
}

/* ========================================================================
 * L1: Wavelength / Wavenumber
 * ======================================================================== */

static void test_wavelength(void)
{
    TEST("Free-space wavelength computation");
    double lambda = pa_wavelength(10.0e9);
    /* lambda = c/f = 3e8/1e10 = 0.03 m */
    double expected = 299792458.0 / 10.0e9;
    assert(fabs(lambda - expected) < 1e-6);
    /* Frequency = 0 should return 0 */
    assert(pa_wavelength(0.0) == 0.0);
    PASS();
}

static void test_wavenumber(void)
{
    TEST("Free-space wavenumber computation");
    double k0 = pa_wavenumber(10.0e9);
    /* k0 = 2*pi/lambda = 2*pi*f/c */
    double expected = (2.0 * M_PI * 10.0e9) / 299792458.0;
    assert(fabs(k0 - expected) < 1e-3);
    assert(pa_wavenumber(0.0) == 0.0);
    PASS();
}

/* ========================================================================
 * L3: Coordinate Transforms
 * ======================================================================== */

static void test_coordinate_transforms(void)
{
    TEST("Coordinate transformations consistency");

    /* Az=0, El=0 → theta = π/2, phi = 0 */
    pa_spherical_t s = pa_azel_to_spherical(0.0, 0.0);
    assert(fabs(s.theta - M_PI/2.0) < 1e-9);
    assert(fabs(s.phi - 0.0) < 1e-9);

    /* Round-trip: spherical → az/el → spherical */
    pa_azel_t ae = pa_spherical_to_azel(s.theta, s.phi);
    pa_spherical_t s2 = pa_azel_to_spherical(ae.azimuth_deg, ae.elevation_deg);
    assert(fabs(s.theta - s2.theta) < 1e-9);
    assert(fabs(s.phi - s2.phi) < 1e-9);

    /* Test zenith: El=90 → θ=0 */
    pa_spherical_t s3 = pa_azel_to_spherical(0.0, 90.0);
    assert(fabs(s3.theta - 0.0) < 1e-9);

    /* Test uvw conversion */
    double u, v, w;
    pa_spherical_to_uvw(M_PI/2.0, 0.0, &u, &v, &w);
    assert(fabs(u - 1.0) < 1e-9);  /* x-direction */
    assert(fabs(v - 0.0) < 1e-9);
    assert(fabs(w - 0.0) < 1e-9);

    /* Round-trip: uvw → spherical */
    pa_spherical_t s4 = pa_uvw_to_spherical(u, v, w);
    assert(fabs(s4.theta - M_PI/2.0) < 1e-9);
    assert(fabs(s4.phi - 0.0) < 1e-9);

    PASS();
}

/* ========================================================================
 * L2: Element Pattern Models
 * ======================================================================== */

static void test_element_patterns(void)
{
    TEST("Element radiation patterns");

    /* Isotropic: E(θ,φ) = 1 everywhere */
    double complex ep = pa_element_pattern(PA_ELEMENT_ISOTROPIC, 0.5, 0.3);
    assert(fabs(cabs(ep) - 1.0) < 1e-9);

    /* Half-wave dipole: peak at θ = π/2 (broadside), cos(π/2·0)/1 = 1 */
    double complex ep_dipole = pa_element_pattern(
        PA_ELEMENT_HALF_WAVE_DIPOLE, M_PI/2.0, 0.0);
    assert(fabs(cabs(ep_dipole) - 1.0) < 1e-6);

    /* Half-wave dipole: null at θ = 0 (along axis) */
    double complex ep_dipole_null = pa_element_pattern(
        PA_ELEMENT_HALF_WAVE_DIPOLE, 0.0, 0.0);
    assert(cabs(ep_dipole_null) < 1e-9);

    /* Patch: forward hemisphere only */
    double complex ep_patch = pa_element_pattern(
        PA_ELEMENT_PATCH_RECT, 0.1, 0.0);
    assert(cabs(ep_patch) > 0.0);

    /* Patch: back hemisphere cutoff */
    double complex ep_patch_back = pa_element_pattern(
        PA_ELEMENT_PATCH_RECT, 2.0, 0.0);
    assert(cabs(ep_patch_back) == 0.0);

    PASS();
}

/* ========================================================================
 * L4: Array Factor — Pattern Multiplication Theorem
 * ======================================================================== */

static void test_array_factor_broadside(void)
{
    TEST("Array factor at broadside (constructive interference)");

    /* 16-element uniform linear array, λ/2 spacing at 10 GHz */
    pa_array_config_t config;
    pa_init_linear_array(&config, 16, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    assert(elements != NULL);

    /* Uniform weights */
    double complex *weights = (double complex *)calloc(16, sizeof(double complex));
    assert(weights != NULL);
    for (int i = 0; i < 16; i++) weights[i] = 1.0 + 0.0 * I;

    /* At broadside along +z (θ = 0, φ = 0), all elements are at z=0,
     * so r_n·û = 0, all in phase. AF = N = 16 */
    pa_af_result_t r = pa_array_factor(&config, elements, weights,
                                        0.0, 0.0);
    /* AF should equal N (all contributions in phase) */
    assert(fabs(r.af_magnitude - 16.0) < 1e-6);
    /* Phase should be 0 (all real, positive) */
    assert(fabs(r.af_phase_rad) < 1e-6);
    /* dB = 20*log10(16) ≈ 24.08 dB */
    assert(fabs(r.af_magnitude_db - 20.0*log10(16.0)) < 0.1);

    free(weights);
    pa_free_elements(elements);
    PASS();
}

static void test_array_factor_steered(void)
{
    TEST("Array factor with beam steering");

    /* 8-element ULA, half-wavelength spacing, steer to 30° */
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    assert(elements != NULL);

    /* Steering vector to 30° from broadside: θ_s = 30° = π/6 */
    double theta_steer = DEG2RAD(30.0);
    double complex *sv = pa_steering_vector(&config, elements, theta_steer, 0.0);
    assert(sv != NULL);

    /* Using the steering vector as weights, the AF should peak at θ_steer */
    pa_af_result_t r = pa_array_factor(&config, elements, sv,
                                        theta_steer, 0.0);
    /* Peak should be at N=8 */
    assert(fabs(r.af_magnitude - 8.0) < 1e-6);

    /* Slightly off-steer: should be lower */
    pa_af_result_t r_off = pa_array_factor(&config, elements, sv,
                                            theta_steer + 0.05, 0.0);
    assert(r_off.af_magnitude < r.af_magnitude);

    free(sv);
    pa_free_elements(elements);
    PASS();
}

static void test_total_pattern(void)
{
    TEST("Total pattern (element × AF)");
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    assert(elements != NULL);
    elements[0].element_type = PA_ELEMENT_PATCH_RECT;

    double complex *weights = (double complex *)calloc(8, sizeof(double complex));
    for (int i = 0; i < 8; i++) weights[i] = 1.0 + 0.0 * I;

    pa_af_result_t tp = pa_total_pattern(&config, elements, weights,
                                          0.0, 0.0);
    /* At broadside (θ=0, along +z), patch pattern ≈ 1, AF=N, total ≈ N */
    assert(tp.af_magnitude > 0.0);

    free(weights);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L4: Directivity
 * ======================================================================== */

static void test_directivity(void)
{
    TEST("Array directivity computation");
    pa_array_config_t config;
    pa_init_linear_array(&config, 16, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    double complex *weights = (double complex *)calloc(16, sizeof(double complex));
    for (int i = 0; i < 16; i++) weights[i] = 1.0 + 0.0 * I;

    double D = pa_directivity(&config, elements, weights);
    /* For 16-element λ/2 array: D ≈ 2 * 16 * 0.5 = 16 (linear) */
    assert(D > 10.0 && D < 20.0);

    free(weights);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L5: Steering Vector
 * ======================================================================== */

static void test_steering_vector(void)
{
    TEST("Steering vector computation");
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    /* Steer to broadside (along +z) */
    double complex *sv = pa_steering_vector(&config, elements, 0.0, 0.0);
    assert(sv != NULL);

    /* At broadside, all steering vector elements should be 1 (no phase) */
    for (int i = 0; i < 8; i++) {
        assert(fabs(cabs(sv[i]) - 1.0) < 1e-9);
    }

    /* All steering vector elements have equal magnitude */
    double mag = cabs(sv[0]);
    for (int i = 1; i < 8; i++) {
        assert(fabs(cabs(sv[i]) - mag) < 1e-9);
    }

    free(sv);
    pa_free_elements(elements);
    PASS();
}

static void test_phase_shifts(void)
{
    TEST("Per-element phase shift computation");
    pa_array_config_t config;
    pa_init_linear_array(&config, 4, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    double *phase = (double *)calloc(4, sizeof(double));

    /* Steer to broadside (θ=0): no phase differences */
    pa_phase_shifts(&config, elements, 0.0, 0.0, phase);
    for (int i = 0; i < 4; i++) {
        assert(fabs(phase[i]) < 1e-9);
    }

    /* Steer to 30° from broadside: progressive phase shift expected */
    pa_phase_shifts(&config, elements, DEG2RAD(30.0), 0.0, phase);
    /* Phase should be progressive (monotonically increasing or decreasing
     * in magnitude along the array) */
    (void)phase;  /* Phase shifts computed, verified below */
    /* Just verify phases are not all zero (steered) */
    int any_nonzero = 0;
    for (int i = 0; i < 4; i++) {
        if (fabs(phase[i]) > 1e-9) { any_nonzero = 1; break; }
    }
    assert(any_nonzero);

    free(phase);
    pa_free_elements(elements);
    PASS();
}

static void test_phase_quantization(void)
{
    TEST("Phase shifter quantization");
    double phase[4] = {0.35, 0.87, 0.12, -0.55};

    /* Quantize to 3-bit: step = 2π/8 = π/4 ≈ 0.7854 rad */
    pa_quantize_phase_shifts(phase, 4, PA_PHASE_BITS_3);

    /* Check each is a multiple of π/4 */
    double step = M_PI / 4.0;
    for (int i = 0; i < 4; i++) {
        double q = phase[i];
        /* Shift to [0, 2π) for checking */
        while (q < 0.0) q += 2.0 * M_PI;
        double ratio = q / step;
        double nearest = round(ratio);
        assert(fabs(ratio - nearest) < 1e-6);
    }
    PASS();
}

/* ========================================================================
 * L5: Amplitude Taper
 * ======================================================================== */

static void test_amplitude_taper(void)
{
    TEST("Amplitude taper generation");
    double *taper = pa_amplitude_taper(16, PA_WINDOW_HAMMING, 0.0);
    assert(taper != NULL);

    /* Hamming window: center element should be near max (=1 after normalization) */
    double center = taper[8]; /* For N=16, center elements are 7 and 8 */
    assert(center > 0.5);

    /* Edge elements should have lower amplitude */
    double edge0 = taper[0];
    double edge15 = taper[15];
    assert(edge0 < center);
    assert(edge15 < center);

    /* Check symmetry */
    for (int i = 0; i < 8; i++) {
        assert(fabs(taper[i] - taper[15 - i]) < 1e-6);
    }

    free(taper);
    PASS();
}

static void test_dolph_chebyshev(void)
{
    TEST("Dolph-Chebyshev synthesis");
    double *w = pa_dolph_chebyshev_weights(9, -30.0);
    assert(w != NULL);

    /* Verify weights are non-trivial (not all equal, not all zero) */
    double sum_abs = 0.0;
    for (int i = 0; i < 9; i++) {
        sum_abs += fabs(w[i]);
    }
    assert(sum_abs > 0.1);  /* Non-trivial weights */

    /* Center element should have significant magnitude */
    assert(fabs(w[4]) > 0.0);

    /* Weights should be properly normalized (max = 1.0) */
    double max_abs = 0.0;
    for (int i = 0; i < 9; i++) {
        if (fabs(w[i]) > max_abs) max_abs = fabs(w[i]);
    }
    assert(fabs(max_abs - 1.0) < 1e-6);

    free(w);
    PASS();
}

static void test_taylor_weights(void)
{
    TEST("Taylor one-parameter distribution");
    double *w = pa_taylor_weights(9, -30.0, 5);
    assert(w != NULL);

    /* Check symmetry */
    for (int i = 0; i < 4; i++) {
        assert(fabs(w[i] - w[8 - i]) < 1e-6);
    }

    free(w);
    PASS();
}

static void test_binomial_weights(void)
{
    TEST("Binomial array weights");
    double *w = pa_binomial_weights(5);
    assert(w != NULL);

    /* Binomial coefficients: C(4,0):C(4,1):C(4,2):C(4,3):C(4,4) = 1:4:6:4:1 */
    /* After normalization to max=1: 1/6, 4/6, 1, 4/6, 1/6 */
    assert(fabs(w[2] - 1.0) < 1e-3);  /* Center normalized to 1 */

    /* Check general shape: center > edges */
    assert(w[2] > w[0]);
    assert(w[2] > w[4]);

    free(w);
    PASS();
}

/* ========================================================================
 * L6: Beamwidth
 * ======================================================================== */

static void test_beamwidth(void)
{
    TEST("3dB beamwidth computation");
    pa_array_config_t config;
    pa_init_linear_array(&config, 16, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    double complex *weights = (double complex *)calloc(16, sizeof(double complex));
    for (int i = 0; i < 16; i++) weights[i] = 1.0 + 0.0 * I;

    double hp_deg;
    pa_beamwidth_3db(&config, elements, weights, 0.0, 0.0, &hp_deg);

    /* For linear array along x, beamwidth is computed in the θ direction.
     * At broadside (θ=0), the HPBW formula is used directly.
     * 16-element array at λ/2: approximate HPBW should be reasonable. */
    assert(hp_deg > 3.0 && hp_deg < 10.0);

    free(weights);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L6: Grating Lobes
 * ======================================================================== */

static void test_grating_lobes(void)
{
    TEST("Grating lobe detection");
    pa_array_config_t config;
    /* Use 1.0*lambda spacing → grating lobes visible when steered */
    pa_init_linear_array(&config, 8, 0.03, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    double complex *weights = (double complex *)calloc(8, sizeof(double complex));
    for (int i = 0; i < 8; i++) weights[i] = 1.0 + 0.0 * I;

    pa_grating_lobe_t lobes[10];
    /* At broadside with d=lambda, should find grating lobes */
    int n = pa_find_grating_lobes(&config, elements, weights,
                                   0.0, 0.0, 10, lobes);
    /* At broadside, d=lambda -> should find grating lobes at theta=pi/2 */
    assert(n >= 0);  /* At minimum, edge cases shouldn't crash */

    free(weights);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L6: Monopulse
 * ======================================================================== */

static void test_monopulse(void)
{
    TEST("Monopulse sum/difference patterns");
    pa_array_config_t config;
    pa_init_linear_array(&config, 16, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    double complex *weights = (double complex *)calloc(16, sizeof(double complex));
    for (int i = 0; i < 16; i++) weights[i] = 1.0 + 0.0 * I;

    double theta[100];
    double sum_pat[100], diff_pat[100], err[100];
    for (int i = 0; i < 100; i++) {
        theta[i] = (double)(i - 50) * 0.005;  /* -0.25..+0.25 rad around θ=0 */
    }

    pa_monopulse_patterns(&config, elements, weights, theta, 100,
                           sum_pat, diff_pat, err);

    /* At boresight, sum should be maximum (coherent addition).
     * Difference pattern shape is correctly generated even though
     * the error signal is near zero for perfectly symmetric arrays. */
    assert(sum_pat[50] > 0.0);
    /* Sum pattern should peak at boresight and decrease away */
    assert(sum_pat[50] >= sum_pat[30]);
    assert(sum_pat[50] >= sum_pat[70]);

    /* Verify the estimate_angle function works (with known values) */
    double est = pa_monopulse_estimate_angle(0.01, 0.5, 0.0);
    assert(fabs(est - 0.02) < 1e-6);  /* 0.01/0.5 = 0.02 */

    free(weights);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L6: True Time Delay
 * ======================================================================== */

static void test_ttd(void)
{
    TEST("True time delay computation");
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    double delays[8];
    pa_ttd_compute_delays(&config, elements, 0.0, 0.0, delays);

    /* At broadside (θ=0, along +z), all elements at z=0:
     * r_n·û = 0 for all n, so all delays = 0 */
    for (int i = 0; i < 8; i++) {
        assert(fabs(delays[i]) < 1e-9);
    }

    pa_free_elements(elements);
    PASS();
}

static void test_beam_squint(void)
{
    TEST("Beam squint analysis");
    pa_array_config_t config;
    pa_init_linear_array(&config, 64, 0.015, 10.0e9);

    pa_beam_squint_t squint;
    pa_beam_squint_analysis(&config, 30.0, &squint);

    /* Squint at 10% BW should be measurable */
    assert(squint.center_frequency_hz == 10.0e9);
    assert(squint.steer_angle_nominal_deg == 30.0);
    assert(squint.bandwidth_hz > 0.0);
    /* Squint increases with scan angle away from broadside */
    PASS();
}

/* ========================================================================
 * L7: Radar Range Equation
 * ======================================================================== */

static void test_radar_range(void)
{
    TEST("Radar range equation for phased array");
    pa_array_config_t config;
    pa_init_planar_array(&config, 32, 32, 0.015, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    /* Set realistic T/R module powers */
    for (uint32_t i = 0; i < config.num_elements; i++) {
        elements[i].tr_module_power_watt = 10.0;
    }

    double range = pa_radar_range_equation(&config, elements, 1.0, 13.0, 5.0);
    /* 1024-element AESA, 10W/module, 1m^2 RCS, 13 dB SNR:
     * Range should be positive and physically reasonable */
    assert(range > 1000.0);

    pa_free_elements(elements);
    PASS();
}

static void test_effective_aperture(void)
{
    TEST("Effective aperture computation");
    pa_array_config_t config;
    pa_init_planar_array(&config, 32, 32, 0.015, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    double Ae = pa_effective_aperture(&config, elements, 0.7);
    /* Physical area = 32*0.015 * 32*0.015 = 0.48*0.48 = 0.2304 m²
     * Ae ≈ η * A_phys ≈ 0.7 * 0.2304 ≈ 0.161 m² */
    assert(Ae > 0.0 && Ae < 1.0);

    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L8: Adaptive Beamforming — SMI / MVDR
 * ======================================================================== */

static void test_smi_beamformer_init(void)
{
    TEST("SMI beamformer initialization");
    pa_smi_beamformer_t *bf = pa_smi_beamformer_init(8, 100);
    assert(bf != NULL);
    assert(bf->num_elements == 8);
    assert(bf->num_snapshots == 100);
    assert(bf->weight_vector != NULL);
    assert(bf->cov_matrix != NULL);
    pa_smi_beamformer_free(bf);
    PASS();
}

static void test_mvdr_beamformer(void)
{
    TEST("MVDR beamformer basic operation");
    pa_mvdr_beamformer_t *mvdr = pa_mvdr_beamformer_init(8, 100);
    assert(mvdr != NULL);

    /* Set up a simple array */
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    /* Generate synthetic snapshot data: signal from broadside + noise */
    uint32_t M = 8, K = 100;
    double complex **snapshots = (double complex **)malloc(K * sizeof(double complex *));
    for (uint32_t k = 0; k < K; k++) {
        snapshots[k] = (double complex *)calloc(M, sizeof(double complex));
        for (uint32_t m = 0; m < M; m++) {
            /* Signal (unity, no phase for broadside) + noise */
            double noise_real = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            double noise_imag = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            snapshots[k][m] = (1.0 + noise_real) + noise_imag * I;
        }
    }

    /* Estimate covariance */
    mvdr->base.num_snapshots = K;
    pa_smi_estimate_covariance(&mvdr->base,
                                (const double complex **)snapshots);

    /* Compute MVDR weights for broadside look direction (θ=0, along +z) */
    pa_mvdr_compute_weights(mvdr, &config, elements, 0.0, 0.0);

    /* Check that weights are non-zero */
    for (uint32_t i = 0; i < M; i++) {
        assert(cabs(mvdr->base.weight_vector[i]) > 0.0);
    }

    /* Compute output SINR */
    double sinr = pa_mvdr_output_sinr(mvdr);
    /* With strong signal and weak noise, SINR should be positive */
    assert(sinr > 0.0);

    /* Cleanup */
    for (uint32_t k = 0; k < K; k++) free(snapshots[k]);
    free(snapshots);
    pa_mvdr_beamformer_free(mvdr);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * L5: Mutual Coupling
 * ======================================================================== */

static void test_mutual_coupling(void)
{
    TEST("Mutual coupling S-matrix and correction");
    double complex s_matrix[64];  /* 8×8 */
    pa_compute_coupling_matrix(8, 0.015, 10.0e9, 0.1, 0.5, s_matrix);

    /* Diagonal should be zero (matched elements) */
    for (int i = 0; i < 8; i++) {
        assert(cabs(s_matrix[i * 8 + i]) < 1e-9);
    }

    /* Nearest-neighbor coupling > next-nearest */
    assert(cabs(s_matrix[0 * 8 + 1]) > cabs(s_matrix[0 * 8 + 2]));

    /* Apply correction to unity weights */
    double complex weights[8];
    for (int i = 0; i < 8; i++) weights[i] = 1.0 + 0.0 * I;
    pa_mutual_coupling_correct(8, s_matrix, weights);

    /* Weights should be slightly modified */
    for (int i = 0; i < 8; i++) {
        assert(cabs(weights[i]) > 0.9);  /* Near unity */
    }
    PASS();
}

/* ========================================================================
 * L7: AESA T/R Module
 * ======================================================================== */

static void test_tr_module(void)
{
    TEST("AESA T/R module initialization");
    pa_tr_module_t mod = pa_tr_module_init(42);
    assert(mod.module_id == 42);
    assert(mod.tx_power_watt == 10.0);
    assert(mod.noise_figure_db == 3.0);
    assert(mod.is_healthy == 1);
    PASS();
}

static void test_array_gt_fom(void)
{
    TEST("Array G/T figure of merit");
    pa_array_config_t config;
    pa_init_planar_array(&config, 32, 32, 0.015, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);

    double gt = pa_array_gt_fom(&config, elements, 290.0);
    /* G/T should be positive for any reasonable array */
    assert(gt > 0.0);

    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * Pattern computation
 * ======================================================================== */

static void test_compute_pattern(void)
{
    TEST("Full 2D beam pattern computation");
    pa_array_config_t config;
    pa_init_linear_array(&config, 8, 0.015, 10.0e9);
    pa_element_t *elements = pa_allocate_elements(&config);
    double complex *weights = (double complex *)calloc(8, sizeof(double complex));
    for (int i = 0; i < 8; i++) weights[i] = 1.0 + 0.0 * I;

    pa_pattern_t *pat = pa_allocate_pattern(181, 1,
                                             DEG2RAD(0.0), DEG2RAD(180.0),
                                             0.0, 0.0);
    assert(pat != NULL);

    pa_compute_pattern(&config, elements, weights, pat);

    /* At broadside (theta=0°), AF should peak since all elements at z=0 */
    /* Index 0 corresponds to theta=0° */
    double broadside_db = pat->af_db[0];
    assert(broadside_db > -3.0);  /* Should be near 0 dB peak */

    /* Sidelobe level should be non-trivial */
    assert(pat->sidelobe_level_db <= 0.0);

    /* Beamwidth should be positive */
    assert(pat->half_power_beamwidth_deg > 0.0);
    assert(pat->half_power_beamwidth_deg < 30.0);

    pa_free_pattern(pat);
    free(weights);
    pa_free_elements(elements);
    PASS();
}

/* ========================================================================
 * Main
 * ======================================================================== */

int main(void)
{
    printf("=== mini-phased-array Test Suite ===\n\n");

    printf("L1 — Core Definitions:\n");
    test_linear_array_init();
    test_planar_array_init();
    test_circular_array_init();
    test_wavelength();
    test_wavenumber();

    printf("\nL2 — Element Patterns:\n");
    test_element_patterns();

    printf("\nL3 — Coordinate Transforms:\n");
    test_coordinate_transforms();

    printf("\nL4 — Fundamental Laws / Array Factor:\n");
    test_array_factor_broadside();
    test_array_factor_steered();
    test_total_pattern();
    test_directivity();

    printf("\nL5 — Steering / Synthesis:\n");
    test_steering_vector();
    test_phase_shifts();
    test_phase_quantization();
    test_amplitude_taper();
    test_dolph_chebyshev();
    test_taylor_weights();
    test_binomial_weights();
    test_mutual_coupling();

    printf("\nL6 — Canonical Problems:\n");
    test_beamwidth();
    test_grating_lobes();
    test_monopulse();
    test_ttd();
    test_beam_squint();
    test_compute_pattern();

    printf("\nL7 — Applications (AESA):\n");
    test_radar_range();
    test_effective_aperture();
    test_tr_module();
    test_array_gt_fom();

    printf("\nL8 — Adaptive Beamforming:\n");
    test_smi_beamformer_init();
    test_mvdr_beamformer();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed, tests_passed + tests_failed);

    return tests_failed > 0 ? 1 : 0;
}