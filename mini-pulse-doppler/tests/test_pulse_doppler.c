#include "../include/pulse_doppler.h"
#include "../include/radar_waveform.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define ASSERT_EQ(expr, expected, fmt) do {     typeof(expected) _val = (expr);     if (_val != (expected)) {         printf("FAIL: got " fmt ", expected " fmt "\n", _val, (expected));         return;     } } while(0)
#define ASSERT_DOUBLE_EQ(expr, expected, tol) do {     double _val = (expr);     if (fabs(_val - (expected)) > (tol)) {         printf("FAIL: got %g, expected %g (tol=%g)\n", _val, (expected), (tol));         return;     } } while(0)

static void test_matched_filter_init(void) {
    TEST("matched_filter_init");
    double complex ref[8];
    for (int i = 0; i < 8; i++) ref[i] = 1.0 + 0.0 * I;

    matched_filter_t mf;
    int rc = matched_filter_init(&mf, ref, 8);
    assert(rc == 0);
    assert(mf.filter_length == 8);
    assert(mf.coefficients != NULL);
    assert(mf.is_normalized == 1);

    double energy_check = 0.0;
    for (size_t i = 0; i < 8; i++) {
        double complex v = mf.coefficients[i];
        energy_check += creal(v)*creal(v) + cimag(v)*cimag(v);
    }
    ASSERT_DOUBLE_EQ(energy_check, 1.0, 1e-6);
    matched_filter_free(&mf);
    PASS();
}

static void test_matched_filter_apply(void) {
    TEST("matched_filter_apply");
    double complex ref[4] = {1+0*I, 1+0*I, 1+0*I, 1+0*I};
    double complex input[8] = {1+0*I,1+0*I,1+0*I,1+0*I,0+0*I,0+0*I,0+0*I,0+0*I};
    double complex output[16];

    matched_filter_t mf;
    matched_filter_init(&mf, ref, 4);
    int rc = matched_filter_apply(&mf, input, 8, output);
    assert(rc == 0);

    double peak = 0.0;
    for (int i = 0; i < 11; i++) {
        double mag = sqrt(creal(output[i])*creal(output[i])
                        + cimag(output[i])*cimag(output[i]));
        if (mag > peak) peak = mag;
    }
    assert(peak > 0.0);
    matched_filter_free(&mf);
    PASS();
}

static void test_albersheim_equation(void) {
    TEST("albersheim_equation");
    double snr;
    int rc = albersheim_equation(1e-6, 0.9, 10, &snr);
    assert(rc == 0);
    assert(snr > 0.0 && snr < 30.0);
    PASS();
}

static void test_shnidman_equation(void) {
    TEST("shnidman_equation");
    double snr0, snr1;
    albersheim_equation(1e-6, 0.9, 10, &snr0);
    shnidman_equation(1e-6, 0.9, 10, SWERLING_1, &snr1);
    assert(snr1 >= snr0);
    PASS();
}

static void test_detection_threshold(void) {
    TEST("detection_threshold");
    double threshold;
    int rc = detection_threshold_compute(1e-6, 1, 1, 1.0, &threshold);
    assert(rc == 0);
    assert(threshold > 0.0);
    PASS();
}

static void test_range_doppler_map(void) {
    TEST("range_doppler_map");
    range_doppler_map_t rdmap;
    int rc = range_doppler_map_allocate(&rdmap, 64, 64);
    assert(rc == 0);
    assert(rdmap.data != NULL);
    assert(rdmap.num_range_bins == 64);
    assert(rdmap.num_doppler_bins == 64);

    rdmap.data[32 + 32 * 64] = 100.0 + 0.0 * I;
    range_doppler_cell_t peaks[10];
    uint32_t found;
    rc = range_doppler_map_find_peaks(&rdmap, 1.0, 10, peaks, &found);
    assert(rc == 0);
    assert(found >= 1);
    range_doppler_map_free(&rdmap);
    PASS();
}

static void test_coherent_integration(void) {
    TEST("coherent_integration");
    double complex pulse_matrix[16];
    double complex integrated[4];
    for (int i = 0; i < 16; i++) pulse_matrix[i] = 1.0 + 0.0 * I;
    int rc = coherent_integration(pulse_matrix, 4, 4, integrated);
    assert(rc == 0);
    ASSERT_DOUBLE_EQ(creal(integrated[0]), 4.0, 1e-6);
    PASS();
}

static void test_cpi_params(void) {
    TEST("cpi_params_from_waveform");
    radar_waveform_params_t wf = {0};
    wf.pri = 0.001;
    wf.prf = 1000.0;
    wf.center_frequency = 10e9;
    cpi_params_t cpi;
    int rc = cpi_params_from_waveform(&wf, 16, &cpi);
    assert(rc == 0);
    assert(cpi.num_pulses == 16);
    ASSERT_DOUBLE_EQ(cpi.wavelength_m, 0.0299792458, 1e-6);
    ASSERT_DOUBLE_EQ(cpi.coherent_gain_db, 10.0*log10(16.0), 0.01);
    PASS();
}

int main(void) {
    printf("=== pulse_doppler tests ===\n");
    test_matched_filter_init();
    test_matched_filter_apply();
    test_albersheim_equation();
    test_shnidman_equation();
    test_detection_threshold();
    test_range_doppler_map();
    test_coherent_integration();
    test_cpi_params();
    printf("\n=== Results: %d/%d passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
