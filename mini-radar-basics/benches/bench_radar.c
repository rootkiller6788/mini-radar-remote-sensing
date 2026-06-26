/**
 * @file    bench_radar.c
 * @brief   Micro-benchmarks for radar signal processing operations
 *
 * Benchmarks: matched filter (time vs FFT), CFAR, RCS sampling,
 *             Doppler FFT, multi-target generation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "radar_core.h"
#include "radar_waveform.h"
#include "radar_detection.h"
#include "radar_doppler.h"
#include "radar_signal_model.h"

#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif

static double now_ms(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC * 1000.0;
}

static void bench_matched_filter(void)
{
    printf("=== Benchmark: Matched Filter ===\n");
    size_t sig_len = 10000, tmpl_len = 1000;
    double complex *sig = calloc(sig_len, sizeof(double complex));
    double complex *tmpl = calloc(tmpl_len, sizeof(double complex));
    double complex *out_td = calloc(sig_len + tmpl_len - 1, sizeof(double complex));
    double complex *out_fft = calloc(sig_len + tmpl_len - 1, sizeof(double complex));
    for (size_t i = 0; i < sig_len; i++) sig[i] = CMPLX((i==5000)?1.0:0.0, 0.0);
    for (size_t i = 0; i < tmpl_len; i++) tmpl[i] = CMPLX(1.0, 0.0);

    double t0 = now_ms();
    radar_matched_filter(sig, sig_len, tmpl, tmpl_len, out_td);
    double t1 = now_ms();
    printf("  Time-domain:  sig=%zu, tmpl=%zu -> %.1f ms\n", sig_len, tmpl_len, t1-t0);

    t0 = now_ms();
    radar_matched_filter_fft(sig, sig_len, tmpl, tmpl_len, out_fft);
    t1 = now_ms();
    printf("  FFT-based:    sig=%zu, tmpl=%zu -> %.1f ms\n", sig_len, tmpl_len, t1-t0);

    free(sig); free(tmpl); free(out_td); free(out_fft);
}

static void bench_cfar(void)
{
    printf("\n=== Benchmark: CFAR Detection ===\n");
    size_t nbins = 100000;
    double *rp = calloc(nbins, sizeof(double));
    int *dets = calloc(nbins, sizeof(int));
    cfar_config_t cfg;
    cfar_config_init(&cfg, CFAR_CA, 4, 16, 1e-4);
    srand(42);
    for (size_t i = 0; i < nbins; i++) {
        double u = (double)rand() / RAND_MAX;
        rp[i] = (u < 1e-12) ? 0.0 : -log(u);
    }
    rp[50000] = 50.0;

    double t0 = now_ms();
    cfar_detect(rp, nbins, &cfg, dets);
    double t1 = now_ms();
    printf("  CA-CFAR: %zu bins, %zu ref, %zu guard -> %.1f ms\n",
           nbins, cfg.num_ref_cells, cfg.num_guard_cells, t1-t0);
    free(rp); free(dets);
}

static void bench_doppler(void)
{
    printf("\n=== Benchmark: Doppler FFT ===\n");
    size_t n_range = 256, n_pulses = 128;
    double complex *data = calloc(n_range * n_pulses, sizeof(double complex));
    double complex *dmap = calloc(n_range * n_pulses, sizeof(double complex));
    for (size_t i = 0; i < n_range * n_pulses; i++)
        data[i] = CMPLX((double)rand()/RAND_MAX, (double)rand()/RAND_MAX);

    double t0 = now_ms();
    radar_doppler_fft_2d(data, n_range, n_pulses, dmap);
    double t1 = now_ms();
    printf("  2D FFT: %zu x %zu -> %.1f ms\n", n_range, n_pulses, t1-t0);
    free(data); free(dmap);
}

int main(void)
{
    printf("=== mini-radar-basics Benchmarks ===\n\n");
    bench_matched_filter();
    bench_cfar();
    bench_doppler();
    printf("\nBenchmarks complete.\n");
    return 0;
}
