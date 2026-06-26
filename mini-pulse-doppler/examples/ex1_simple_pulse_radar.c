#include "../include/radar_waveform.h"
#include "../include/pulse_doppler.h"
#include "../include/doppler_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Ex1: Simple Pulse Radar with LFM Compression ===\n\n");

    radar_waveform_params_t wf = {0};
    wf.pulse_width = 10e-6;
    wf.bandwidth = 10e6;
    wf.center_frequency = 10e9;
    wf.prf = 20000.0;
    wf.pri = 1.0 / 20000.0;
    wf.duty_cycle = wf.pulse_width * wf.prf;
    wf.sampling_rate = 50e6;
    wf.num_samples = 1024;
    wf.type = WAVEFORM_LFM_UP;

    double complex *chirp = malloc(wf.num_samples * sizeof(double complex));
    lfm_params_t lfm = {.chirp_rate = wf.bandwidth / wf.pulse_width,
                        .time_bw_product = wf.pulse_width * wf.bandwidth,
                        .is_upchirp = 1};
    lfm_chirp_generate(&wf, &lfm, chirp);

    printf("Waveform: LFM up-chirp\n");
    printf("  Pulse width:    %.1f us\n", wf.pulse_width * 1e6);
    printf("  Bandwidth:      %.1f MHz\n", wf.bandwidth / 1e6);
    printf("  T-B Product:    %.1f\n", lfm.time_bw_product);
    printf("  PRF:            %.1f kHz\n", wf.prf / 1e3);
    printf("  Range resolution: %.2f m\n", 299792458.0 / (2.0 * wf.bandwidth));
    printf("  Max unambiguous range: %.1f km\n", 299792458.0 / (2.0 * wf.prf) / 1e3);

    /* Simulate a target at 5 km (within receive window) */
    double target_range = 5000.0;
    double delay = 2.0 * target_range / 299792458.0;
    uint32_t delay_samples = (uint32_t)(delay * wf.sampling_rate);
    printf("\nTarget at %.1f km -> delay = %.1f us (%u samples)\n",
           target_range/1e3, delay*1e6, delay_samples);

    double complex *rx = calloc(2048, sizeof(double complex));
    for (uint32_t i = 0; i < wf.num_samples && (i + delay_samples) < 2048; i++)
        rx[i + delay_samples] = chirp[i] * 0.1;

    double complex *compressed = calloc(3072, sizeof(double complex));
    pulse_compress_lfm(chirp, rx, 2048, wf.num_samples, compressed);

    double peak_mag = 0.0;
    uint32_t peak_idx = 0;
    for (uint32_t i = 0; i < 3072; i++) {
        double mag = sqrt(creal(compressed[i])*creal(compressed[i])
                        + cimag(compressed[i])*cimag(compressed[i]));
        if (mag > peak_mag) { peak_mag = mag; peak_idx = i; }
    }
    printf("Compressed peak at sample %u\n", peak_idx);
    printf("Estimated range: %.1f km\n",
           (double)peak_idx / wf.sampling_rate * 299792458.0 / 2.0 / 1e3);

    free(chirp); free(rx); free(compressed);
    printf("\nDone.\n");
    return 0;
}
