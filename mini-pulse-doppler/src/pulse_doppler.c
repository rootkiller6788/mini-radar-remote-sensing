#include "pulse_doppler.h"
#include "radar_waveform.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int matched_filter_init(matched_filter_t *mf,
                        const double complex *reference, size_t length)
{
    if (mf == NULL || reference == NULL || length == 0) return -1;
    mf->filter_length = length;
    mf->coefficients = (double complex *)malloc(length * sizeof(double complex));
    if (mf->coefficients == NULL) return -1;
    double energy = 0.0;
    for (size_t n = 0; n < length; n++) {
        double complex val = reference[n];
        energy += creal(val)*creal(val) + cimag(val)*cimag(val);
    }
    for (size_t n = 0; n < length; n++)
        mf->coefficients[n] = conj(reference[length - 1 - n]);
    if (energy > 1e-30) {
        double scale = 1.0 / sqrt(energy);
        for (size_t n = 0; n < length; n++)
            mf->coefficients[n] *= scale;
        mf->is_normalized = 1;
    } else { mf->is_normalized = 0; }
    mf->processing_gain_db = 10.0 * log10(energy);
    return 0;
}

void matched_filter_free(matched_filter_t *mf)
{
    if (mf != NULL && mf->coefficients != NULL) {
        free(mf->coefficients);
        mf->coefficients = NULL;
    }
}

int matched_filter_apply(const matched_filter_t *mf,
                         const double complex *input, size_t input_len,
                         double complex *output)
{
    if (mf == NULL || input == NULL || output == NULL) return -1;
    if (mf->coefficients == NULL || mf->filter_length == 0) return -1;
    if (input_len == 0) return -1;
    size_t M = mf->filter_length;
    size_t out_len = input_len + M - 1;
    for (size_t n = 0; n < out_len; n++) output[n] = 0.0 + 0.0 * I;
    for (size_t n = 0; n < out_len; n++) {
        double complex accum = 0.0 + 0.0 * I;
        for (size_t k = 0; k < M; k++) {
            if (k <= n && (n - k) < input_len)
                accum += mf->coefficients[k] * input[n - k];
        }
        output[n] = accum;
    }
    return 0;
}

int matched_filter_apply_freq_domain(const double complex *reference,
                                     const double complex *input,
                                     size_t length, double complex *output)
{
    if (reference == NULL || input == NULL || output == NULL) return -1;
    if (length == 0) return -1;
    matched_filter_t mf;
    if (matched_filter_init(&mf, reference, length) != 0) return -1;
    int rc = matched_filter_apply(&mf, input, length, output);
    matched_filter_free(&mf);
    return rc;
}

int pulse_compress_lfm(const double complex *chirp_reference,
                       const double complex *received_signal,
                       size_t signal_length, size_t ref_length,
                       double complex *compressed)
{
    if (chirp_reference == NULL || received_signal == NULL || compressed == NULL)
        return -1;
    if (signal_length == 0 || ref_length == 0) return -1;
    matched_filter_t mf;
    if (matched_filter_init(&mf, chirp_reference, ref_length) != 0) return -1;
    int rc = matched_filter_apply(&mf, received_signal, signal_length, compressed);
    matched_filter_free(&mf);
    return rc;
}

int range_gate_extract(const double complex *pulse_data,
                       uint32_t num_range_bins, uint32_t num_pulses,
                       uint32_t gate_start, uint32_t gate_width,
                       double complex *gate_data)
{
    if (pulse_data == NULL || gate_data == NULL) return -1;
    if (num_range_bins == 0 || num_pulses == 0) return -1;
    if (gate_start + gate_width > num_range_bins) return -1;
    if (gate_width == 0) return -1;
    for (uint32_t p = 0; p < num_pulses; p++) {
        double complex sum = 0.0 + 0.0 * I;
        for (uint32_t g = 0; g < gate_width; g++) {
            uint32_t idx = (uint32_t)(gate_start + g) + p * num_range_bins;
            sum += pulse_data[idx];
        }
        gate_data[p] = sum;
    }
    return 0;
}

int coherent_integration(double complex *pulse_matrix,
                         uint32_t num_range_bins, uint32_t num_pulses,
                         double complex *integrated)
{
    if (pulse_matrix == NULL || integrated == NULL) return -1;
    if (num_range_bins == 0 || num_pulses == 0) return -1;
    for (uint32_t r = 0; r < num_range_bins; r++) {
        double complex sum = 0.0 + 0.0 * I;
        for (uint32_t p = 0; p < num_pulses; p++)
            sum += pulse_matrix[r + p * num_range_bins];
        integrated[r] = sum;
    }
    return 0;
}

int noncoherent_integration(const double *magnitude_matrix,
                            uint32_t num_range_bins, uint32_t num_pulses,
                            double *integrated)
{
    if (magnitude_matrix == NULL || integrated == NULL) return -1;
    if (num_range_bins == 0 || num_pulses == 0) return -1;
    for (uint32_t r = 0; r < num_range_bins; r++) {
        double sum = 0.0;
        for (uint32_t p = 0; p < num_pulses; p++)
            sum += magnitude_matrix[r + p * num_range_bins];
        integrated[r] = sum;
    }
    return 0;
}

int compute_coherent_gain(uint32_t num_pulses, double *gain_db) {
    if (gain_db == NULL) return -1;
    if (num_pulses == 0) { *gain_db = 0.0; return 0; }
    *gain_db = 10.0 * log10((double)num_pulses);
    return 0;
}

int compute_noncoherent_gain(uint32_t num_pulses, double *gain_db) {
    if (gain_db == NULL) return -1;
    if (num_pulses == 0) { *gain_db = 0.0; return 0; }
    *gain_db = 10.0 * 0.8 * log10((double)num_pulses);
    return 0;
}

int compute_integration_loss(uint32_t num_pulses,
                             double pfa, double pd, double *loss_db) {
    if (loss_db == NULL) return -1;
    if (num_pulses <= 1) { *loss_db = 0.0; return 0; }
    (void)pfa; (void)pd;
    *loss_db = 10.0*log10((double)num_pulses)-10.0*0.8*log10((double)num_pulses);
    return 0;
}

int detection_threshold_compute(double pfa, uint32_t num_samples,
                                int noise_power_known, double noise_power,
                                double *threshold)
{
    if (threshold == NULL) return -1;
    if (pfa <= 0.0 || pfa >= 1.0 || num_samples == 0) return -1;
    double sigma = (noise_power_known && noise_power > 0.0)
                   ? sqrt(noise_power) : 1.0;
    (void)num_samples;
    *threshold = sigma * sqrt(-2.0 * log(pfa));
    return 0;
}

int albersheim_equation(double pfa, double pd,
                        uint32_t num_pulses, double *required_snr_db)
{
    if (required_snr_db == NULL) return -1;
    if (pfa <= 0.0 || pfa >= 1.0 || pd <= 0.0 || pd >= 1.0) return -1;
    if (num_pulses == 0) return -1;
    double A = log(0.62 / pfa);
    double B = log(pd / (1.0 - pd));
    double N = (double)num_pulses;
    double term = 6.2 + 4.54 / sqrt(N + 0.44);
    double arg = A + 0.12 * A * B + 1.7 * B;
    *required_snr_db = -5.0 * log10(N) + term * log10(arg);
    return 0;
}

int shnidman_equation(double pfa, double pd, uint32_t num_pulses,
                      swerling_case_t swerling, double *required_snr_db)
{
    if (required_snr_db == NULL) return -1;
    if (pfa <= 0.0 || pfa >= 1.0 || pd <= 0.0 || pd >= 1.0) return -1;
    if (num_pulses == 0) return -1;
    double snr_base;
    if (albersheim_equation(pfa, pd, num_pulses, &snr_base) != 0) return -1;
    double adj = 0.0;
    switch (swerling) {
    case SWERLING_0: adj = 0.0; break;
    case SWERLING_1: adj = (pd>=0.9)?8.0:(pd>=0.7)?6.0:(pd>=0.5)?4.0:3.0; break;
    case SWERLING_2: adj = (pd>=0.9)?4.0:(pd>=0.7)?3.0:(pd>=0.5)?2.0:1.5; break;
    case SWERLING_3: adj = (pd>=0.9)?3.0:(pd>=0.7)?2.0:(pd>=0.5)?1.5:1.0; break;
    case SWERLING_4: adj = (pd>=0.9)?1.5:(pd>=0.7)?1.0:(pd>=0.5)?0.5:0.25; break;
    default: return -1;
    }
    *required_snr_db = snr_base + adj;
    return 0;
}

int target_detection_report_init(target_detection_t *report) {
    if (report == NULL) return -1;
    report->target_id = 0;
    report->range_m = 0.0; report->velocity_mps = 0.0;
    report->doppler_hz = 0.0; report->snr_db = 0.0;
    report->azimuth_deg = 0.0; report->elevation_deg = 0.0;
    report->rcs_dbsm = 0.0; report->detected = DETECTION_H0;
    return 0;
}

int cpi_params_from_waveform(const radar_waveform_params_t *wf,
                             uint32_t num_pulses, cpi_params_t *cpi) {
    if (wf == NULL || cpi == NULL || num_pulses == 0) return -1;
    cpi->num_pulses = num_pulses;
    cpi->pri_seconds = wf->pri; cpi->prf_hz = wf->prf;
    cpi->center_frequency_hz = wf->center_frequency;
    cpi->wavelength_m = 299792458.0 / wf->center_frequency;
    cpi->integration_time_s = (double)num_pulses * wf->pri;
    cpi->coherent_gain_db = 10.0 * log10((double)num_pulses);
    return 0;
}

int range_doppler_map_allocate(range_doppler_map_t *rdmap,
                               uint32_t num_range, uint32_t num_doppler) {
    if (rdmap == NULL || num_range == 0 || num_doppler == 0) return -1;
    size_t total = (size_t)num_range * num_doppler;
    rdmap->data = (double complex *)calloc(total, sizeof(double complex));
    if (rdmap->data == NULL) return -1;
    rdmap->num_range_bins = num_range;
    rdmap->num_doppler_bins = num_doppler;
    rdmap->range_resolution_m = 0.0;
    rdmap->doppler_resolution_hz = 0.0;
    rdmap->velocity_resolution_mps = 0.0;
    rdmap->max_unambiguous_range_m = 0.0;
    rdmap->max_unambiguous_velocity_mps = 0.0;
    return 0;
}

void range_doppler_map_free(range_doppler_map_t *rdmap) {
    if (rdmap != NULL && rdmap->data != NULL) {
        free(rdmap->data);
        rdmap->data = NULL;
    }
}

static int compare_peaks_desc(const void *a, const void *b) {
    const range_doppler_cell_t *pa = (const range_doppler_cell_t *)a;
    const range_doppler_cell_t *pb = (const range_doppler_cell_t *)b;
    if (pa->magnitude < pb->magnitude) return 1;
    if (pa->magnitude > pb->magnitude) return -1;
    return 0;
}

int range_doppler_map_find_peaks(const range_doppler_map_t *rdmap,
                                 double threshold, uint32_t max_peaks,
                                 range_doppler_cell_t *peaks, uint32_t *num_found) {
    if (rdmap == NULL || peaks == NULL || num_found == NULL) return -1;
    if (rdmap->data == NULL) return -1;
    uint32_t NR = rdmap->num_range_bins, ND = rdmap->num_doppler_bins;
    uint32_t count = 0;
    for (uint32_t r = 0; r < NR && count < max_peaks; r++) {
        for (uint32_t d = 0; d < ND && count < max_peaks; d++) {
            double complex val = rdmap->data[r + d * NR];
            double mag = sqrt(creal(val)*creal(val)+cimag(val)*cimag(val));
            if (mag < threshold) continue;
            int is_peak = 1;
            for (int dr = -1; dr <= 1 && is_peak; dr++) {
                for (int dd = -1; dd <= 1 && is_peak; dd++) {
                    if (dr == 0 && dd == 0) continue;
                    int nr = (int)r + dr, nd = (int)d + dd;
                    if (nr >= 0 && nr < (int)NR && nd >= 0 && nd < (int)ND) {
                        double complex nv = rdmap->data[nr + nd * NR];
                        double nm = sqrt(creal(nv)*creal(nv)+cimag(nv)*cimag(nv));
                        if (nm > mag) is_peak = 0;
                    }
                }
            }
            if (is_peak) {
                peaks[count].range_bin = r;
                peaks[count].doppler_bin = d;
                peaks[count].magnitude = mag;
                peaks[count].phase_rad = atan2(cimag(val), creal(val));
                peaks[count].range_meters = (double)r*rdmap->range_resolution_m;
                peaks[count].doppler_hz = ((double)d-(double)ND/2.0)
                                          * rdmap->doppler_resolution_hz;
                peaks[count].velocity_mps = peaks[count].doppler_hz
                    * rdmap->velocity_resolution_mps/rdmap->doppler_resolution_hz;
                count++;
            }
        }
    }
    qsort(peaks, count, sizeof(range_doppler_cell_t), compare_peaks_desc);
    *num_found = count;
    return 0;
}
