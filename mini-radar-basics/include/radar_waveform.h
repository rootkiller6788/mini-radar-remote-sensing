/**
 * @file    radar_waveform.h
 * @brief   Radar Waveform Design - L1 Definitions + L5 Algorithms
 */
#ifndef RADAR_WAVEFORM_H
#define RADAR_WAVEFORM_H
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif
typedef enum { WAVEFORM_SIMPLE_PULSE, WAVEFORM_LFM_UP, WAVEFORM_LFM_DOWN, WAVEFORM_NLFM, WAVEFORM_BARKER, WAVEFORM_FRANK, WAVEFORM_P1, WAVEFORM_P2, WAVEFORM_P3, WAVEFORM_P4, WAVEFORM_STEPPED_FREQ, WAVEFORM_COSTAS } waveform_type_t;
typedef enum { WINDOW_RECTANGULAR, WINDOW_HAMMING, WINDOW_HANN, WINDOW_BLACKMAN, WINDOW_KAISER, WINDOW_CHEBYSHEV, WINDOW_TAYLOR } window_type_t;
typedef struct { waveform_type_t type; double pulse_width_s; double bandwidth_hz; double sample_rate_hz; double center_freq_hz; double chirp_rate_hz_per_s; size_t num_samples; int code_length; double complex *samples; window_type_t window_type; double *window_coeffs; } radar_waveform_t;
int radar_waveform_rect_pulse(radar_waveform_t *wf, double tau, double fs);
int radar_waveform_lfm(radar_waveform_t *wf, double tau, double bw, double fs, int up);
int radar_waveform_apply_window(radar_waveform_t *wf, window_type_t wt);
int radar_waveform_barker(radar_waveform_t *wf, int len, double tau, double fs);
int radar_waveform_frank(radar_waveform_t *wf, int nc, double tau, double fs);
int radar_matched_filter(const double complex *sig, size_t slen, const double complex *tmpl, size_t tlen, double complex *out);
int radar_matched_filter_fft(const double complex *sig, size_t slen, const double complex *tmpl, size_t tlen, double complex *out);
double radar_psl_compute(const double *comp, size_t len, size_t *pk);
double radar_isl_compute(const double *comp, size_t len, size_t pk);
void radar_waveform_free(radar_waveform_t *wf);
#endif
