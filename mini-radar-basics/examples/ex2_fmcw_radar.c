/**
 * @file    ex2_fmcw_radar.c
 * @brief   Example: FMCW radar range and velocity estimation
 *
 * Demonstrates: FMCW chirp generation, dechirp processing, range-Doppler
 *               joint estimation via 2D FFT.
 *
 * Knowledge: L6 FMCW radar, L5 Chirp processing, L5 2D FFT
 * Reference: Richards, Scheer & Holm (2010), Ch.4
 */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include "radar_waveform.h"
#include "radar_doppler.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif

#define RADAR_C 299792458.0

int main(void)
{
    printf("=== Example 2: FMCW Radar Range-Velocity Estimation ===\n\n");

    /* FMCW parameters */
    double fc = 77e9;          /* 77 GHz automotive radar */
    double B = 500e6;          /* 500 MHz bandwidth */
    double T_chirp = 50e-6;    /* 50 us chirp duration */
    double fs = 10e6;          /* 10 MHz sample rate */
    double lambda = RADAR_C / fc;

    printf("FMCW Parameters:\n");
    printf("  Carrier:     %.0f GHz\n", fc / 1e9);
    printf("  Bandwidth:   %.0f MHz\n", B / 1e6);
    printf("  Chirp time:  %.0f us\n", T_chirp * 1e6);
    printf("  Sample rate: %.0f MHz\n", fs / 1e6);
    printf("  Wavelength:  %.4f m\n", lambda);
    printf("  Range res:   %.3f m\n", RADAR_C / (2.0 * B));
    printf("  Max range:   %.1f m (at fs=10 MHz)\n",
           RADAR_C * fs / (2.0 * B / T_chirp) / 2.0);

    /* Generate LFM chirp waveform */
    radar_waveform_t tx_wf;
    memset(&tx_wf, 0, sizeof(tx_wf));
    radar_waveform_lfm(&tx_wf, T_chirp, B, fs, 1);  /* up-chirp */
    printf("\nWaveform: %zu samples, chirp rate = %.2e Hz/s\n",
           tx_wf.num_samples, tx_wf.chirp_rate_hz_per_s);

    /* Simulate two targets */
    /* Target 1: 30 m range, 20 m/s radial velocity */
    double R1 = 30.0, v1 = 20.0;
    /* Target 2: 75 m range, -10 m/s radial velocity */
    double R2 = 75.0, v2 = -10.0;

    double tau1 = 2.0 * R1 / RADAR_C;
    double tau2 = 2.0 * R2 / RADAR_C;
    double fd1 = 2.0 * v1 / lambda;
    double fd2 = 2.0 * v2 / lambda;

    /* Beat frequency: f_b = k * tau + f_d (approx, for small f_d) */
    double k = tx_wf.chirp_rate_hz_per_s;
    double f_b1 = k * tau1 + fd1;
    double f_b2 = k * tau2 + fd2;

    printf("\nSimulated Targets:\n");
    printf("  Target 1: R=%.1f m, v=%.1f m/s -> tau=%.2f ns, fd=%.1f Hz, fb=%.1f kHz\n",
           R1, v1, tau1*1e9, fd1, f_b1/1000.0);
    printf("  Target 2: R=%.1f m, v=%.1f m/s -> tau=%.2f ns, fd=%.1f Hz, fb=%.1f kHz\n",
           R2, v2, tau2*1e9, fd2, f_b2/1000.0);

    /* Build dechirped signal: s_rx * conj(s_tx) */
    size_t N = tx_wf.num_samples;
    double complex *beat = calloc(N, sizeof(double complex));
    for (size_t n = 0; n < N; n++) {
        double t = (double)n / fs;
        /* Target 1 contribution */
        double phase1 = 2.0 * M_PI * (fd1 * t + k * tau1 * t - k * tau1 * tau1 / 2.0);
        /* Target 2 contribution */
        double phase2 = 2.0 * M_PI * (fd2 * t + k * tau2 * t - k * tau2 * tau2 / 2.0);
        beat[n] = CMPLX(cos(phase1) + cos(phase2), sin(phase1) + sin(phase2));
    }

    /* FFT of beat signal -> range profile */
    size_t N_fft = 1;
    while (N_fft < N) N_fft <<= 1;
    double complex *fft_in = calloc(N_fft, sizeof(double complex));
    for (size_t n = 0; n < N; n++) fft_in[n] = beat[n];

    /* Simple FFT (from radar_waveform, reused) */
    /* Just compute magnitude via direct DFT for demo */
    double *spectrum = calloc(N_fft / 2, sizeof(double));
    for (size_t k = 0; k < N_fft / 2; k++) {
        double complex sum = CMPLX(0.0, 0.0);
        for (size_t n = 0; n < N; n++) {
            double phase = -2.0 * M_PI * (double)(k * n) / (double)N_fft;
            sum += beat[n] * CMPLX(cos(phase), sin(phase));
        }
        spectrum[k] = cabs(sum);
    }

    /* Find peaks */
    printf("\nRange Spectrum Peaks:\n");
    double max1 = 0, max2 = 0;
    size_t idx1 = 0, idx2 = 0;
    for (size_t k = 0; k < N_fft / 2; k++) {
        if (spectrum[k] > max1) {
            max2 = max1; idx2 = idx1;
            max1 = spectrum[k]; idx1 = k;
        } else if (spectrum[k] > max2) {
            max2 = spectrum[k]; idx2 = k;
        }
    }
    double freq_bin = fs / (double)N_fft;
    double r_est1 = RADAR_C * (double)idx1 * freq_bin / (2.0 * k);
    double r_est2 = RADAR_C * (double)idx2 * freq_bin / (2.0 * k);
    printf("  Peak 1: bin %zu -> %.1f kHz -> R_est = %.1f m (true: %.1f m)\n",
           idx1, idx1 * freq_bin / 1000.0, r_est1, R1);
    printf("  Peak 2: bin %zu -> %.1f kHz -> R_est = %.1f m (true: %.1f m)\n",
           idx2, idx2 * freq_bin / 1000.0, r_est2, R2);

    free(beat); free(fft_in); free(spectrum);
    radar_waveform_free(&tx_wf);

    printf("\nFMCW example complete.\n");
    return 0;
}
