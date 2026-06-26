#include "ambiguity_function.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int ambiguity_surface_alloc(uint32_t num_delays, uint32_t num_dopplers,
                            ambiguity_surface_t *surface)
{
    if (surface == NULL || num_delays == 0 || num_dopplers == 0) return -1;
    size_t total = (size_t)num_delays * num_dopplers;
    surface->surface = (double complex *)calloc(total, sizeof(double complex));
    if (surface->surface == NULL) return -1;
    surface->num_delays = num_delays;
    surface->num_dopplers = num_dopplers;
    return 0;
}

void ambiguity_surface_free(ambiguity_surface_t *surface)
{
    if (surface != NULL && surface->surface != NULL) {
        free(surface->surface);
        surface->surface = NULL;
    }
}

int ambiguity_function_compute(const double complex *waveform,
                               size_t signal_length, double sampling_rate,
                               uint32_t num_delays, uint32_t num_dopplers,
                               double max_doppler_hz, ambiguity_surface_t *result)
{
    if (waveform == NULL || result == NULL || signal_length == 0) return -1;
    if (result->surface == NULL) return -1;

    double dt = 1.0 / sampling_rate;
    double delay_step = dt;
    double doppler_step = 2.0 * max_doppler_hz / (double)(num_dopplers - 1);

    result->delay_resolution_s = delay_step;
    result->doppler_resolution_hz = doppler_step;
    result->max_delay_s = (double)(num_delays - 1) * delay_step;
    result->max_doppler_hz = max_doppler_hz;

    for (uint32_t d = 0; d < num_dopplers; d++) {
        double fd = -max_doppler_hz + (double)d * doppler_step;
        for (uint32_t k = 0; k < num_delays; k++) {
            (void)((double)k * delay_step); /* tau = k * dt used implicitly in loop */
            double complex accum = 0.0 + 0.0 * I;
            for (size_t n = 0; n < signal_length; n++) {
                if (n + k >= signal_length) break;
                double phase = -2.0 * M_PI * fd * (double)n * dt;
                accum += waveform[n] * conj(waveform[n + k])
                       * (cos(phase) + sin(phase) * I);
            }
            result->surface[k + d * num_delays] = accum;
        }
    }
    return 0;
}

int ambiguity_cut_range(const ambiguity_surface_t *surface, double *range_cut)
{
    if (surface == NULL || range_cut == NULL || surface->surface == NULL) return -1;
    uint32_t mid_dop = surface->num_dopplers / 2;
    for (uint32_t k = 0; k < surface->num_delays; k++) {
        double complex val = surface->surface[k + mid_dop * surface->num_delays];
        range_cut[k] = creal(val)*creal(val) + cimag(val)*cimag(val);
    }
    return 0;
}

int ambiguity_cut_doppler(const ambiguity_surface_t *surface, double *doppler_cut)
{
    if (surface == NULL || doppler_cut == NULL || surface->surface == NULL) return -1;
    for (uint32_t d = 0; d < surface->num_dopplers; d++) {
        double complex val = surface->surface[0 + d * surface->num_delays];
        doppler_cut[d] = creal(val)*creal(val) + cimag(val)*cimag(val);
    }
    return 0;
}

int ambiguity_metrics_compute(const ambiguity_surface_t *surface,
                              double signal_energy, ambiguity_metrics_t *metrics)
{
    if (surface == NULL || metrics == NULL || surface->surface == NULL) return -1;

    uint32_t ND = surface->num_delays;
    uint32_t NF = surface->num_dopplers;
    double max_val = 0.0;
    double total_energy = 0.0;

    for (uint32_t d = 0; d < NF; d++) {
        for (uint32_t k = 0; k < ND; k++) {
            double complex v = surface->surface[k + d * ND];
            double mag_sq = creal(v)*creal(v) + cimag(v)*cimag(v);
            total_energy += mag_sq;
            if (mag_sq > max_val) max_val = mag_sq;
        }
    }

    metrics->peak_value = max_val;
    double peak_sidelobe = 0.0;
    for (uint32_t d = 0; d < NF; d++) {
        for (uint32_t k = 0; k < ND; k++) {
            if (k == 0 && d == (NF/2)) continue;
            double complex v = surface->surface[k + d * ND];
            double mag_sq = creal(v)*creal(v) + cimag(v)*cimag(v);
            if (mag_sq > peak_sidelobe) peak_sidelobe = mag_sq;
        }
    }
    metrics->peak_sidelobe_ratio_db = 10.0 * log10(peak_sidelobe / max_val);
    metrics->delay_doppler_volume = total_energy * surface->delay_resolution_s
                                    * surface->doppler_resolution_hz;
    metrics->q_function_value = signal_energy * signal_energy / (metrics->delay_doppler_volume + 1e-30);

    uint32_t mid = NF/2;
    double half_max = max_val / 2.0;
    int found_3db = 0;
    for (uint32_t k = 1; k < ND && !found_3db; k++) {
        double complex v = surface->surface[k + mid * ND];
        double mag_sq = creal(v)*creal(v) + cimag(v)*cimag(v);
        if (mag_sq <= half_max) {
            metrics->_3db_delay_width_s = (double)k * surface->delay_resolution_s;
            found_3db = 1;
        }
    }
    found_3db = 0;
    for (uint32_t d = 1; d < NF && !found_3db; d++) {
        double complex v = surface->surface[0 + d * ND];
        double mag_sq = creal(v)*creal(v) + cimag(v)*cimag(v);
        if (mag_sq <= half_max) {
            metrics->_3db_doppler_width_hz = (double)d * surface->doppler_resolution_hz;
            found_3db = 1;
        }
    }
    metrics->integrated_sidelobe_ratio_db = 10.0
        * log10((total_energy - max_val) / (max_val + 1e-30));
    return 0;
}

int ambiguity_compare(const ambiguity_surface_t *af1,
                      const ambiguity_surface_t *af2,
                      double *similarity_metric)
{
    if (af1 == NULL || af2 == NULL || similarity_metric == NULL) return -1;
    if (af1->surface == NULL || af2->surface == NULL) return -1;

    uint32_t N = af1->num_delays * af1->num_dopplers;
    uint32_t M = af2->num_delays * af2->num_dopplers;
    uint32_t min_N = (N < M) ? N : M;

    double correlation = 0.0;
    double energy1 = 0.0, energy2 = 0.0;
    for (uint32_t i = 0; i < min_N; i++) {
        double complex v1 = af1->surface[i];
        double complex v2 = af2->surface[i];
        double m1 = creal(v1)*creal(v1) + cimag(v1)*cimag(v1);
        double m2 = creal(v2)*creal(v2) + cimag(v2)*cimag(v2);
        if (m1 > 1e-30 && m2 > 1e-30)
            correlation += sqrt(m1 * m2);
        energy1 += m1;
        energy2 += m2;
    }
    *similarity_metric = (energy1+energy2 > 1e-30)
        ? correlation / sqrt(energy1 * energy2) : 0.0;
    return 0;
}

int rd_coupling_analyze(double chirp_rate, double carrier_frequency_hz,
                        rd_coupling_t *coupling)
{
    if (coupling == NULL) return -1;
    double c = 299792458.0;
    double lambda = c / carrier_frequency_hz;
    coupling->coupling_factor = chirp_rate;
    coupling->range_error_per_doppler_hz = c * lambda / (2.0 * chirp_rate);
    coupling->doppler_error_per_delay_s = chirp_rate;
    coupling->is_range_doppler_coupled = (fabs(chirp_rate) > 1e-30) ? 1 : 0;
    return 0;
}

int q_function_compute(const ambiguity_surface_t *surface, double *q_value)
{
    if (surface == NULL || q_value == NULL || surface->surface == NULL) return -1;
    double volume = 0.0;
    for (uint32_t d = 0; d < surface->num_dopplers; d++) {
        for (uint32_t k = 0; k < surface->num_delays; k++) {
            double complex v = surface->surface[k + d * surface->num_delays];
            volume += creal(v)*creal(v) + cimag(v)*cimag(v);
        }
    }
    *q_value = volume * surface->delay_resolution_s * surface->doppler_resolution_hz;
    return 0;
}
