#include "../include/radar_waveform.h"
#include "../include/pulse_doppler.h"
#include "../include/doppler_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Ex2: Range-Doppler Map Generation ===\n\n");

    uint32_t num_range = 128, num_pulses = 64;

    radar_waveform_params_t wf = {0};
    wf.pulse_width = 5e-6;
    wf.bandwidth = 20e6;
    wf.center_frequency = 10e9;
    wf.prf = 5000.0;
    wf.pri = 1.0 / 5000.0;
    wf.duty_cycle = wf.pulse_width * wf.prf;
    wf.sampling_rate = 40e6;
    wf.num_samples = num_range;

    double wavelength = 299792458.0 / wf.center_frequency;
    printf("Wavelength: %.4f m\n", wavelength);
    printf("Max unambiguous velocity: %.1f m/s\n", wavelength*wf.prf/4.0);
    printf("Velocity resolution: %.2f m/s\n", wavelength*wf.prf/(2.0*num_pulses));
    printf("Doppler resolution: %.2f Hz\n", wf.prf/num_pulses);

    double complex *chirp = malloc(num_range * sizeof(double complex));
    lfm_params_t lfm = {.chirp_rate = wf.bandwidth/wf.pulse_width,
                        .time_bw_product = wf.pulse_width*wf.bandwidth,
                        .is_upchirp = 1};
    lfm_chirp_generate(&wf, &lfm, chirp);

    double complex *pulse_matrix = calloc(num_range * num_pulses, sizeof(double complex));

    /* Simulate target at range bin 40 with Doppler shift */
    uint32_t target_bin = 40;
    double target_velocity = 30.0;
    double target_doppler = 2.0 * target_velocity / wavelength;
    printf("\nTarget: range bin %u, velocity %.1f m/s, Doppler %.1f Hz\n",
           target_bin, target_velocity, target_doppler);

    for (uint32_t p = 0; p < num_pulses; p++) {
        double phase = 2.0 * M_PI * target_doppler * (double)p * wf.pri;
        double complex doppler_shift = cos(phase) + sin(phase) * I;
        for (uint32_t r = 0; r < num_range && (r + target_bin) < num_range; r++) {
            if (r < num_range - target_bin)
                pulse_matrix[r + target_bin + p * num_range] = chirp[r] * doppler_shift * 1.0;
        }
    }

    range_doppler_map_t rdmap;
    range_doppler_map_allocate(&rdmap, num_range, num_pulses);

    doppler_fft_params_t dparams = {0};
    dparams.fft_size = num_pulses;
    dparams.num_pulses = num_pulses;
    dparams.prf_hz = wf.prf;
    dparams.doppler_resolution_hz = wf.prf / num_pulses;
    dparams.wavelength_m = wavelength;
    dparams.max_velocity_mps = wavelength * wf.prf / 4.0;
    dparams.velocity_resolution_mps = wavelength * wf.prf / (2.0 * num_pulses);
    dparams.window = WINDOW_HAMMING;

    doppler_fft_process(pulse_matrix, num_range, num_pulses, &dparams, &rdmap);

    range_doppler_cell_t peaks[10];
    uint32_t found;
    rdmap.range_resolution_m = 299792458.0 / (2.0 * wf.bandwidth);
    rdmap.doppler_resolution_hz = dparams.doppler_resolution_hz;
    rdmap.velocity_resolution_mps = dparams.velocity_resolution_mps;
    range_doppler_map_find_peaks(&rdmap, 1.0, 10, peaks, &found);

    printf("\nDetected %u peaks:\n", found);
    for (uint32_t i = 0; i < found && i < 5; i++) {
        printf("  Peak %u: range=%.1fm (bin %u), velocity=%.1fm/s (bin %u), mag=%.1f\n",
               i+1, peaks[i].range_meters, peaks[i].range_bin,
               peaks[i].velocity_mps, peaks[i].doppler_bin, peaks[i].magnitude);
    }

    range_doppler_map_free(&rdmap);
    free(chirp); free(pulse_matrix);
    printf("\nDone.\n");
    return 0;
}
