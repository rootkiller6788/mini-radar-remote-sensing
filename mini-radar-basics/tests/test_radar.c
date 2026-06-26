/**
 * @file    test_radar.c
 * @brief   Test suite for mini-radar-basics: comprehensive assert-based tests
 *
 * Covers L1-L6: all core APIs exercised with known-answer tests.
 *
 * Usage: make test
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include "radar_core.h"
#include "radar_waveform.h"
#include "radar_detection.h"
#include "radar_doppler.h"
#include "radar_signal_model.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(expr, msg) do { \
    if (expr) { tests_passed++; } \
    else { fprintf(stderr, "FAIL: %s\n", msg); tests_failed++; } \
} while(0)

#define TEST_FEQ(a, b, tol, msg) do { \
    if (fabs((a) - (b)) < (tol)) { tests_passed++; } \
    else { fprintf(stderr, "FAIL: %s (%.6f vs %.6f)\n", msg, a, b); tests_failed++; } \
} while(0)

/* ===== L1/L4: Radar Range Equation Tests ================================== */

static void test_range_equation(void)
{
    radar_params_t params;
    int rc = radar_params_init(&params,
        10e9,       /* 10 GHz X-band */
        1000.0,     /* 1 kW peak power */
        30.0,       /* 30 dB gain */
        10e-6,      /* 10 us pulse */
        1000.0,     /* 1 kHz PRF */
        10e6,       /* 10 MHz BW */
        3.0,        /* 3 dB NF */
        6.0,        /* 6 dB loss */
        290.0,      /* 290 K antenna temp */
        RADAR_MODE_PULSE);
    TEST(rc == 0, "radar_params_init");

    /* Verify derived parameters */
    TEST(params.wavelength_m > 0.0, "wavelength > 0");
    TEST_FEQ(params.wavelength_m, 299792458.0 / 10e9, 0.01, "wavelength calc");
    TEST_FEQ(params.pri_s, 0.001, 1e-9, "PRI = 1/PRF");
    TEST_FEQ(params.antenna_gain_linear, 1000.0, 1.0, "gain dB to linear");

    /* Range resolution */
    double res_pulse = radar_range_resolution(10e-6, 0);
    TEST_FEQ(res_pulse, 299792458.0 * 10e-6 / 2.0, 10.0, "pulse range resolution");
    double res_bw = radar_range_resolution(10e6, 1);
    TEST_FEQ(res_bw, 299792458.0 / (2.0 * 10e6), 1.0, "BW range resolution");

    /* Unambiguous range */
    double r_unamb = radar_unambiguous_range(1000.0);
    TEST_FEQ(r_unamb, 299792458.0 / 2000.0, 1.0, "unambiguous range");
    TEST(radar_unambiguous_range(0.0) == INFINITY, "PRF=0 -> inf range");

    /* Unambiguous velocity */
    double v_unamb = radar_unambiguous_velocity(params.wavelength_m, 1000.0);
    TEST(v_unamb > 0.0, "unambiguous velocity > 0");

    /* Received power at 10 km with 1 m^2 RCS */
    double pr = radar_received_power(&params, 1.0, 10000.0);
    TEST(pr > 0.0, "received power > 0");
    TEST(pr < params.peak_power_w, "received power < transmit power");

    /* SNR at 10 km */
    double snr = radar_snr(&params, 1.0, 10000.0);
    TEST(snr > 0.0, "SNR > 0 at 10 km");

    /* SNR at 1 km should be much larger */
    double snr_near = radar_snr(&params, 1.0, 1000.0);
    TEST(snr_near > snr, "SNR larger at shorter range");

    /* SNR dB */
    double snrdb = radar_snr_db(&params, 1.0, 10000.0);
    TEST(snrdb < 100.0, "SNR dB reasonable");

    /* Maximum range */
    double rmax = radar_max_range(&params, 1.0, 10.0);
    TEST(rmax > 0.0, "max range > 0");

    /* Bistatic equation (degenerate to monostatic) */
    double pr_bi = radar_bistatic_power(params.peak_power_w,
        params.antenna_gain_linear, params.antenna_gain_linear,
        params.wavelength_m, 1.0, 10000.0, 10000.0,
        params.system_loss_linear);
    TEST_FEQ(pr_bi, pr, pr * 0.01, "bistatic == monostatic at same range");

    /* System noise temperature */
    double tsys = radar_system_noise_temp(290.0, 2.0);
    TEST_FEQ(tsys, 290.0 + 290.0 * (2.0 - 1.0), 0.1, "system noise temp");

    /* Noise power */
    double pn = radar_noise_power(tsys, 1e6);
    TEST(pn > 0.0, "noise power > 0");

    /* Antenna gain from aperture */
    double gain = radar_antenna_gain_from_aperture(0.1, 0.03);
    TEST(gain > 1.0, "antenna gain > 1");

    /* Beamwidth */
    double bw = radar_beamwidth_circular(0.03, 1.0);
    TEST(bw > 0.0, "beamwidth > 0");

    /* Directivity (Kraus) */
    double dir = radar_directivity_kraus(0.03, 0.03);
    TEST(dir > 1.0, "directivity > 1");

    /* Duty cycle */
    double dc = radar_duty_cycle(10e-6, 1000.0);
    TEST_FEQ(dc, 0.01, 0.001, "duty cycle = tau * PRF");

    /* Average power */
    double pavg = radar_average_power(1000.0, dc);
    TEST_FEQ(pavg, 10.0, 0.01, "average power = peak * duty");

    /* Coherent integration gain */
    double snr_coh = radar_coherent_integration_gain(10, 1.0);
    TEST_FEQ(snr_coh, 10.0, 0.1, "coherent integration gain = N");

    /* Non-coherent integration gain */
    double gnc = radar_ncoherent_integration_gain(10);
    TEST(gnc > 1.0, "non-coherent gain > 1");
    TEST(gnc < 10.0, "non-coherent gain < N");

    /* dB conversion */
    TEST_FEQ(lin2db(10.0), 10.0, 0.01, "lin2db(10)=10");
    TEST_FEQ(db2lin(10.0), 10.0, 0.01, "db2lin(10)=10");
}

/* ===== L1/L5: RCS and Waveform Tests ====================================== */

static void test_rcs_waveform(void)
{
    /* RCS model */
    radar_rcs_t rcs;
    TEST(radar_rcs_init(&rcs, 10.0, RCS_CONSTANT, 0.0) == 0, "rcs init");
    TEST_FEQ(rcs.mean_rcs_m2, 10.0, 0.01, "rcs constant 10 dBsm = 10 m^2");
    double r1 = radar_rcs_sample(&rcs);
    TEST_FEQ(r1, 10.0, 0.01, "rcs sample constant");

    radar_rcs_init(&rcs, 0.0, RCS_SWERLING_I, 0.1);  /* 0 dBsm = 1 m^2 */
    double r2 = radar_rcs_sample(&rcs);
    TEST(r2 > 0.0, "swerling I sample > 0");

    radar_rcs_init(&rcs, 0.0, RCS_SWERLING_III, 0.1);
    double r3 = radar_rcs_sample(&rcs);
    TEST(r3 > 0.0, "swerling III sample > 0");

    /* Rectangular pulse waveform */
    radar_waveform_t wf;
    memset(&wf, 0, sizeof(wf));
    TEST(radar_waveform_rect_pulse(&wf, 10e-6, 10e6) == 0, "rect pulse gen");
    TEST(wf.num_samples > 0, "rect pulse has samples");
    TEST(wf.samples != NULL, "rect pulse allocated");
    TEST_FEQ(cabs(wf.samples[0]), 1.0, 0.01, "rect pulse amplitude");
    radar_waveform_free(&wf);

    /* LFM chirp */
    memset(&wf, 0, sizeof(wf));
    TEST(radar_waveform_lfm(&wf, 10e-6, 10e6, 20e6, 1) == 0, "LFM up-chirp");
    TEST(wf.num_samples > 0, "lfm has samples");
    TEST(wf.type == WAVEFORM_LFM_UP, "lfm type up");
    TEST_FEQ(wf.chirp_rate_hz_per_s, 1e12, 1e10, "chirp rate = B/tau");
    radar_waveform_free(&wf);

    /* LFM down-chirp */
    memset(&wf, 0, sizeof(wf));
    TEST(radar_waveform_lfm(&wf, 10e-6, 5e6, 10e6, 0) == 0, "LFM down-chirp");
    TEST(wf.type == WAVEFORM_LFM_DOWN, "lfm type down");
    radar_waveform_free(&wf);

    /* Window application */
    memset(&wf, 0, sizeof(wf));
    radar_waveform_lfm(&wf, 10e-6, 5e6, 10e6, 1);
    TEST(radar_waveform_apply_window(&wf, WINDOW_HAMMING) == 0, "hamming window");
    TEST(wf.window_coeffs != NULL, "window coeffs allocated");
    radar_waveform_free(&wf);

    /* Barker code */
    memset(&wf, 0, sizeof(wf));
    TEST(radar_waveform_barker(&wf, 13, 10e-6, 10e6) == 0, "barker-13 gen");
    TEST(wf.num_samples > 0, "barker has samples");
    radar_waveform_free(&wf);

    /* Barker bad length */
    memset(&wf, 0, sizeof(wf));
    TEST(radar_waveform_barker(&wf, 6, 10e-6, 10e6) == -1, "barker invalid length");

    /* Frank code */
    memset(&wf, 0, sizeof(wf));
    TEST(radar_waveform_frank(&wf, 4, 10e-6, 10e6) == 0, "frank-4 gen");
    TEST(wf.code_length == 16, "frank code length = N^2");
    radar_waveform_free(&wf);

    /* Matched filter (time-domain) */
    double complex signal[100];
    double complex tmpl[10];
    double complex output[109];
    for (int i = 0; i < 100; i++) signal[i] = (i == 50) ? CMPLX(1.0,0.0) : CMPLX(0.0,0.0);
    for (int i = 0; i < 10; i++) tmpl[i] = CMPLX(1.0,0.0);
    TEST(radar_matched_filter(signal, 100, tmpl, 10, output) == 0, "matched filter");
    TEST(cabs(output[50]) >= 0.99, "matched filter peak at correct delay");

    /* FFT-based matched filter */
    double complex output2[109];
    TEST(radar_matched_filter_fft(signal, 100, tmpl, 10, output2) == 0, "fft matched filter");
    TEST_FEQ(cabs(output[59]), cabs(output2[59]), 1e-6, "time vs fft match");

    /* PSL computation */
    double comp[109];
    for (int i = 0; i < 109; i++) comp[i] = cabs(output[i]);
    size_t pk;
    double psl = radar_psl_compute(comp, 109, &pk);
    TEST(psl < 0.0, "PSL negative");

    /* ISL computation */
    double isl = radar_isl_compute(comp, 109, pk);
    TEST(isl <= 0.0, "ISL <= 0");
}

/* ===== L4/L5: Detection Theory Tests ====================================== */

static void test_detection(void)
{
    /* CFAR config */
    cfar_config_t cfg;
    TEST(cfar_config_init(&cfg, CFAR_CA, 2, 8, 1e-4) == 0, "cfar config init");
    TEST(cfg.threshold_factor > 0.0, "cfar threshold factor > 0");

    /* CFAR detection - simple test */
    double rp[256];
    int dets[256];
    for (int i = 0; i < 256; i++) rp[i] = 1.0;  /* uniform noise */
    rp[128] = 100.0;  /* strong target */
    TEST(cfar_detect(rp, 256, &cfg, dets) == 0, "cfar detect");
    TEST(dets[128] == 1, "cfar detects strong target");

    /* OS-CFAR */
    int dets2[256];
    cfg.type = CFAR_OS;
    cfar_config_init(&cfg, CFAR_OS, 2, 8, 1e-4);
    TEST(cfar_os_detect(rp, 256, &cfg, dets2) == 0, "os-cfar detect");

    /* Detection threshold */
    double thresh = radar_detection_threshold(1e-4, 1);
    TEST(thresh > 0.0, "detection threshold > 0");

    /* Marcum Q */
    double q = radar_marcum_q(1, 0.0, 3.0);
    TEST(q >= 0.0 && q <= 1.0, "Marcum Q(1,0,3) in [0,1]");
    double q_zero = radar_marcum_q(1, 0.0, 0.0);
    TEST_FEQ(q_zero, 1.0, 0.01, "Marcum Q(a,0) = 1");
    double q_small = radar_marcum_q(1, 0.0, 10.0);
    TEST(q_small < 1e-5, "Marcum Q(1,0,10) near zero (large threshold, no signal)");

    /* P_d for non-fluctuating target */
    double pd = radar_pd_marcum(10.0, 1e-4, 1);  /* 10 dB SNR */
    TEST(pd >= 0.0 && pd <= 1.0, "Pd in [0,1]");

    /* Swerling I */
    double pd1 = radar_pd_swerling1(10.0, 1e-4);
    TEST(pd1 >= 0.0 && pd1 <= 1.0, "Pd SW1 in [0,1]");
    /* Fluctuation loss: Swerling I should be close to but can differ from Marcum */
    TEST(pd1 > 0.0, "Swerling I Pd > 0");

    /* Swerling II */
    double pd2 = radar_pd_swerling2(10.0, 1e-4, 10);
    TEST(pd2 >= 0.0 && pd2 <= 1.0, "Pd SW2 in [0,1]");

    /* Albersheim */
    double snr_req = radar_albersheim_snr(0.9, 1e-4, 1);
    TEST(snr_req > 0.0, "Albersheim SNR > 0");

    /* ROC curve */
    double pfa_arr[20], pd_arr[20];
    TEST(radar_roc_curve(10.0, pfa_arr, pd_arr, 20, 1, 0) == 0, "ROC curve");
    TEST(pd_arr[0] < 1.0, "ROC first point < 1");

    /* Incomplete gamma */
    double ginc = radar_gamma_inc(1.0, 0.5);
    TEST(ginc >= 0.0 && ginc <= 1.0, "gamma_inc in [0,1]");
    double ginc0 = radar_gamma_inc(1.0, 0.0);
    TEST_FEQ(ginc0, 0.0, 0.01, "gamma_inc(a,0) = 0");

    /* Sort */
    double arr[5] = {3.0, 1.0, 4.0, 1.5, 2.0};
    radar_dsort(arr, 5);
    TEST(arr[0] <= arr[4], "dsort: first <= last");
    TEST_FEQ(arr[0], 1.0, 0.01, "dsort min");
    TEST_FEQ(arr[4], 4.0, 0.01, "dsort max");
}

/* ===== L1/L5: Doppler Processing Tests ==================================== */

static void test_doppler(void)
{
    /* Doppler shift */
    double fd = radar_doppler_shift(100.0, 0.03);
    TEST_FEQ(fd, 2.0 * 100.0 / 0.03, 0.1, "doppler shift = 2v/lambda");

    /* Velocity from Doppler */
    double vr = radar_doppler_to_velocity(fd, 0.03);
    TEST_FEQ(vr, 100.0, 0.1, "velocity from doppler");

    /* Doppler resolution */
    double df = radar_doppler_resolution(1000.0, 64);
    TEST_FEQ(df, 1000.0 / 64.0, 0.01, "doppler resolution = PRF/N");

    /* Blind speeds */
    double v_blind1 = radar_blind_speed(0.03, 1000.0, 1);
    TEST_FEQ(v_blind1, 0.03 * 1000.0 / 2.0, 0.01, "first blind speed");

    /* MTI filter - single canceller */
    mti_filter_t mti;
    TEST(mti_filter_init(&mti, MTI_SINGLE_CANCEL, NULL, 0) == 0, "mti init single");
    TEST(mti.order == 1, "mti order 1");
    TEST(mti.coeffs[0] == 1.0 && mti.coeffs[1] == -1.0, "mti coeffs [1,-1]");

    /* MTI rejects stationary target */
    double complex x1 = CMPLX(1.0, 0.0);
    double complex x2 = CMPLX(1.0, 0.0);  /* same -> zero output */
    mti_filter_apply(&mti, x1);
    double complex y = mti_filter_apply(&mti, x2);
    TEST_FEQ(cabs(y), 0.0, 0.01, "mti cancels stationary target");

    /* MTI passes moving target */
    mti_filter_t mti2;
    mti_filter_init(&mti2, MTI_SINGLE_CANCEL, NULL, 0);
    double complex x3 = CMPLX(1.0, 0.0);
    double complex x4 = CMPLX(0.0, 1.0);  /* 90 deg phase change */
    mti_filter_apply(&mti2, x3);
    double complex y2 = mti_filter_apply(&mti2, x4);
    TEST(cabs(y2) > 0.0, "mti passes moving target");
    mti_filter_free(&mti);
    mti_filter_free(&mti2);

    /* MTI double canceller */
    mti_filter_t mti3;
    TEST(mti_filter_init(&mti3, MTI_DOUBLE_CANCEL, NULL, 0) == 0, "mti double");
    mti_filter_free(&mti3);

    /* MTI triple canceller */
    mti_filter_t mti4;
    TEST(mti_filter_init(&mti4, MTI_TRIPLE_CANCEL, NULL, 0) == 0, "mti triple");
    mti_filter_free(&mti4);

    /* MTI frequency response */
    mti_filter_t mti5;
    mti_filter_init(&mti5, MTI_SINGLE_CANCEL, NULL, 0);
    mti5.prf_hz = 1000.0;
    double freqs[5] = {0.0, 100.0, 250.0, 400.0, 500.0};
    double resp[5];
    TEST(mti_filter_response(&mti5, freqs, 5, resp) == 0, "mti freq response");
    TEST(resp[0] < resp[2], "mti: response at DC < response at PRF/4");
    mti_filter_free(&mti5);

    /* Range-Doppler FFT */
    size_t n_range = 4, n_pulses = 8;
    double complex *data = calloc(n_range * n_pulses, sizeof(double complex));
    /* Put a slow sinusoid at range bin 2 */
    for (size_t p = 0; p < n_pulses; p++) {
        double phase = 2.0 * M_PI * (double)p * 2.0 / (double)n_pulses;
        data[2 * n_pulses + p] = CMPLX(cos(phase), sin(phase));
    }
    double complex *dmap = calloc(n_range * n_pulses, sizeof(double complex));
    TEST(radar_doppler_fft_2d(data, n_range, n_pulses, dmap) == 0, "doppler fft 2d");
    free(data); free(dmap);

    /* Doppler bin to frequency */
    double fbin = radar_doppler_bin_to_freq(2, 64, 1000.0);
    TEST(fbin >= -500.0 && fbin <= 500.0, "doppler bin in range");

    /* Velocity from phase */
    double v_phase = radar_velocity_from_phase(
        CMPLX(1.0, 0.0), CMPLX(0.0, 1.0), 0.03, 0.001);
    TEST(fabs(v_phase) > 0.0, "velocity from phase non-zero");

    /* Doppler ambiguity resolution */
    double resolved;
    int ra = radar_doppler_resolve_ambiguity(100.0, 1000.0, -900.0, 1000.0, &resolved);
    /* 100 Hz from PRF1 and -900 Hz from PRF2 both alias to 1100 Hz true */
    TEST(ra >= 0, "doppler ambiguity resolved");
}

/* ===== L3/L6: Signal Model Tests ========================================== */

static void test_signal_model(void)
{
    /* Target init */
    radar_target_t tgt;
    TEST(radar_target_init(&tgt, 1, 5000.0, 50.0, 0.0, 0.0, 1.0, 0.03) == 0,
         "target init");
    TEST_FEQ(tgt.doppler_hz, 2.0 * 50.0 / 0.03, 0.1, "target doppler");
    TEST_FEQ(tgt.delay_s, 2.0 * 5000.0 / 299792458.0, 1e-8, "target delay");

    /* Range bin mapping */
    size_t bin = radar_range_to_bin(5000.0, 20e6);
    TEST(bin > 0, "range bin > 0");
    double r_back = radar_bin_to_range(bin, 20e6);
    TEST(fabs(r_back - 5000.0) < 100.0, "range bin round-trip");

    /* Number of range bins */
    size_t nb = radar_num_range_bins(15000.0, 20e6);
    TEST(nb > 0, "num range bins > 0");

    /* Generate return */
    double complex tx[100];
    for (int i = 0; i < 100; i++) tx[i] = CMPLX(1.0, 0.0);
    double complex rx[500];
    int rc = radar_generate_return(tx, 100, &tgt, 20e6, rx, 500, 0.0);
    TEST(rc == 0, "generate return");

    /* Multi-target */
    radar_target_t tgts[2];
    radar_target_init(&tgts[0], 1, 5000.0, 50.0, 0.0, 0.0, 1.0, 0.03);
    radar_target_init(&tgts[1], 2, 10000.0, -30.0, 0.1, 0.0, 0.5, 0.03);
    rc = radar_generate_multitarget(tx, 100, tgts, 2, 20e6, rx, 500, 0.0);
    TEST(rc == 0, "multi-target return");

    /* AWGN */
    double complex noise[100];
    rc = radar_awgn(noise, 100, 1.0, 42);
    TEST(rc == 0, "awgn gen");

    /* Clutter model */
    clutter_model_t clutter;
    rc = clutter_model_init(&clutter, -20.0, CLUTTER_RAYLEIGH, 0.1, 10.0);
    TEST(rc == 0, "clutter init");

    /* Clutter generate */
    double complex clut[200];
    rc = clutter_generate(&clutter, 200, 5.0, 0.0, clut, 200);
    TEST(rc == 0, "clutter generate");

    /* Clutter RCS cell */
    double c_rcs = clutter_rcs_cell(0.01, 10000.0, 0.03, 15.0);
    TEST(c_rcs > 0.0, "clutter RCS cell > 0");

    /* Range profile */
    double profile[401];
    rc = radar_range_profile(rx, 500, tx, 100, profile);
    TEST(rc == 0, "range profile");
}

/* ===== main =============================================================== */

int main(void)
{
    printf("=== mini-radar-basics Test Suite ===\n\n");

    test_range_equation();
    printf("  [range equation]       pass=%d fail=%d\n",
           tests_passed - (tests_failed > 0 ? 0 : 0),
           tests_failed > 0 ? 1 : 0);

    int before = tests_passed + tests_failed;
    test_rcs_waveform();
    printf("  [rcs & waveform]       pass=%d fail=%d\n",
           tests_passed + tests_failed - before - tests_failed,
           tests_failed);

    before = tests_passed + tests_failed;
    test_detection();
    printf("  [detection theory]     pass=%d fail=%d\n",
           tests_passed + tests_failed - before - tests_failed,
           tests_failed);

    before = tests_passed + tests_failed;
    test_doppler();
    printf("  [doppler processing]   pass=%d fail=%d\n",
           tests_passed + tests_failed - before - tests_failed,
           tests_failed);

    before = tests_passed + tests_failed;
    test_signal_model();
    printf("  [signal model]         pass=%d fail=%d\n",
           tests_passed + tests_failed - before - tests_failed,
           tests_failed);

    printf("\n=== Results: %d passed, %d failed ===\n",
           tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
