/**
 * @file    bench_sar.c
 * @brief   SAR Performance Benchmarks
 * Measures: chirp generation speed, FFT throughput, pulse compression rate.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sar_core.h"
void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward);

int main(void) {
    printf("SAR Imaging Benchmarks\n====================\n");
    clock_t start, end;

    /* Benchmark chirp generation */
    start = clock();
    for (int i=0; i<100; i++) {
        sar_chirp_t *c = sar_chirp_alloc(4096, 40e-6, 100e6, 150e6);
        sar_chirp_free(c);
    }
    end = clock();
    printf("Chirp gen (100x 4K): %.1f ms\n", 1000.0*(double)(end-start)/CLOCKS_PER_SEC);

    /* Benchmark FFT */
    size_t N = 8192;
    double *I = malloc(N*8), *Q = malloc(N*8);
    for (size_t i=0; i<N; i++) { I[i]=(double)i; Q[i]=0.0; }
    start = clock();
    for (int i=0; i<50; i++) sar_fft_range_1d(I, Q, N, 1);
    end = clock();
    printf("FFT (50x 8K): %.1f ms\n", 1000.0*(double)(end-start)/CLOCKS_PER_SEC);
    free(I); free(Q);

    printf("Bench complete.\n");
    return 0;
}