/**
 * @file    radar_waveform.c
 * @brief   Radar waveform generation, pulse compression, matched filtering
 *
 * Knowledge covered:
 *   L1: Waveform types and parameters
 *   L3: Matched filter impulse response, FFT-based convolution
 *   L5: LFM generation, Barker/Frank codes, PSL/ISL computation
 *   L6: Pulse compression (canonical radar signal processing)
 *
 * Reference: Richards, Scheer & Holm (2010), Ch.4, 20.
 *            Levanon & Mozeson, "Radar Signals" (2004), Ch.2-6.
 */
#include "radar_waveform.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef CMPLX
#define CMPLX(r, i) ((double complex)((double)(r) + I * (double)(i)))
#endif

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "radar_waveform: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ©¤©¤©¤ L5: Rectangular Pulse Generation ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Generates a simple rectangular pulse: s[n] = 1+j0 for all samples in pulse.
 *
 * Time-bandwidth product = 1.0 (no pulse compression gain possible).
 * Range resolution = c / (2*B) = c * tau / 2 (same as unmodulated pulse).
 *
 * The rectangular pulse has poor spectral containment (sinc-shaped
 * spectrum with -13.2 dB first sidelobe).
 */
int radar_waveform_rect_pulse(radar_waveform_t *waveform,
                              double pulse_width_s, double sample_rate_hz)
{
    if (!waveform || pulse_width_s <= 0.0 || sample_rate_hz <= 0.0)
        return -1;

    memset(waveform, 0, sizeof(*waveform));
    waveform->type = WAVEFORM_SIMPLE_PULSE;
    waveform->pulse_width_s = pulse_width_s;
    waveform->bandwidth_hz = 1.0 / pulse_width_s;  /* 1/tau ¡ª nominal bandwidth */
    waveform->sample_rate_hz = sample_rate_hz;
    waveform->center_freq_hz = 0.0;
    waveform->chirp_rate_hz_per_s = 0.0;
    waveform->num_samples = (size_t)(pulse_width_s * sample_rate_hz);
    if (waveform->num_samples < 1) waveform->num_samples = 1;
    waveform->code_length = 0;
    waveform->window_type = WINDOW_RECTANGULAR;

    waveform->samples = (double complex*)safe_malloc(
        waveform->num_samples * sizeof(double complex));
    for (size_t i = 0; i < waveform->num_samples; i++) {
        waveform->samples[i] = CMPLX(1.0, 0.0);
    }
    waveform->window_coeffs = NULL;
    return 0;
}

/* ©¤©¤©¤ L5: Linear FM (Chirp) Generation ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Complex baseband LFM chirp:
 *
 *   s[n] = exp(j * pi * k * (n/f_s)^2),  for n = 0...N-1
 *
 * where k = B/tau is the chirp rate (Hz/s).
 *
 * For an up-chirp (k > 0):
 *   Instantaneous frequency: f_i(t) = k * t,  t ¡Ê [-tau/2, tau/2]
 *   Starts at -B/2, sweeps linearly to +B/2.
 *
 * Pulse compression ratio: PCR = tau * B = TBP
 * Range sidelobe level (rect window): -13.2 dB (first sidelobe)
 * Range resolution (after compression): Delta_R = c / (2*B)
 */
int radar_waveform_lfm(radar_waveform_t *waveform,
                       double pulse_width_s, double bandwidth_hz,
                       double sample_rate_hz, int up_chirp)
{
    if (!waveform || pulse_width_s <= 0.0 || bandwidth_hz <= 0.0 ||
        sample_rate_hz <= 0.0)
        return -1;
    /* Nyquist: sample_rate_hz should be >= bandwidth_hz */
    if (sample_rate_hz < bandwidth_hz)
        return -1;

    memset(waveform, 0, sizeof(*waveform));
    waveform->type = up_chirp ? WAVEFORM_LFM_UP : WAVEFORM_LFM_DOWN;
    waveform->pulse_width_s = pulse_width_s;
    waveform->bandwidth_hz = bandwidth_hz;
    waveform->sample_rate_hz = sample_rate_hz;
    waveform->chirp_rate_hz_per_s = bandwidth_hz / pulse_width_s;
    waveform->num_samples = (size_t)(pulse_width_s * sample_rate_hz);
    if (waveform->num_samples < 1) waveform->num_samples = 1;

    waveform->samples = (double complex*)safe_malloc(
        waveform->num_samples * sizeof(double complex));

    double k = waveform->chirp_rate_hz_per_s;
    if (!up_chirp) k = -k;  /* down-chirp: negative chirp rate */
    double ts = 1.0 / sample_rate_hz;
    double tau = pulse_width_s;

    for (size_t n = 0; n < waveform->num_samples; n++) {
        double t = (double)n * ts - tau / 2.0;  /* centered at t=0 */
        double phase = M_PI * k * t * t;
        waveform->samples[n] = CMPLX(cos(phase), sin(phase));
    }
    waveform->window_type = WINDOW_RECTANGULAR;
    waveform->window_coeffs = NULL;
    waveform->code_length = 0;
    return 0;
}

/* ©¤©¤©¤ L5: Window Application for Sidelobe Control ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Apply an amplitude taper to the waveform to reduce range sidelobes.
 *
 * Window functions and their characteristics:
 *   - Hamming:  w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1)),  PSL ¡Ö -43 dB
 *   - Hann:     w[n] = 0.50 - 0.50*cos(2*pi*n/(N-1)),  PSL ¡Ö -32 dB
 *   - Blackman: w[n] = 0.42 - 0.50*cos(...) + 0.08*cos(...), PSL ¡Ö -58 dB
 *   - Kaiser:   w[n] = I0(beta*sqrt(1-((2n/(N-1))-1)^2)) / I0(beta)
 *
 * Window SNR loss (mismatch loss relative to matched filter):
 *   - Hamming: ~1.34 dB
 *   - Hann:    ~1.76 dB
 *   - Blackman: ~2.37 dB
 *   - Kaiser (beta=6): ~1.6 dB
 *
 * The trade-off: lower sidelobes cost SNR loss and broader mainlobe.
 */
int radar_waveform_apply_window(radar_waveform_t *waveform, window_type_t win_type)
{
    if (!waveform || !waveform->samples || waveform->num_samples == 0)
        return -1;

    size_t N = waveform->num_samples;
    double *w = (double*)safe_malloc(N * sizeof(double));

    switch (win_type) {
    case WINDOW_RECTANGULAR:
        for (size_t n = 0; n < N; n++) w[n] = 1.0;
        break;

    case WINDOW_HAMMING:
        for (size_t n = 0; n < N; n++)
            w[n] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)n / (double)(N - 1));
        break;

    case WINDOW_HANN:
        for (size_t n = 0; n < N; n++)
            w[n] = 0.50 - 0.50 * cos(2.0 * M_PI * (double)n / (double)(N - 1));
        break;

    case WINDOW_BLACKMAN: {
        for (size_t n = 0; n < N; n++) {
            double a = 2.0 * M_PI * (double)n / (double)(N - 1);
            w[n] = 0.42 - 0.50 * cos(a) + 0.08 * cos(2.0 * a);
        }
        break;
    }

    case WINDOW_KAISER: {
        /* Kaiser window with beta = 6.0 (approx -42 dB PSL) */
        double beta = 6.0;
        /* Compute I0(beta) via series expansion */
        double i0_beta = 1.0;
        {
            double term = 1.0, b2 = (beta / 2.0) * (beta / 2.0);
            for (int k = 1; k < 25; k++) {
                term *= b2 / (double)(k * k);
                i0_beta += term;
            }
        }
        for (size_t n = 0; n < N; n++) {
            double x = 2.0 * (double)n / (double)(N - 1) - 1.0;
            double arg = beta * sqrt(1.0 - x * x);
            /* I0(arg) via series */
            double i0_arg = 1.0;
            double term = 1.0, h2 = (arg / 2.0) * (arg / 2.0);
            for (int k = 1; k < 25; k++) {
                term *= h2 / (double)(k * k);
                i0_arg += term;
            }
            w[n] = i0_arg / i0_beta;
        }
        break;
    }

    case WINDOW_CHEBYSHEV:
    case WINDOW_TAYLOR:
    default:
        /* Fall back to Hamming for unsupported windows */
        for (size_t n = 0; n < N; n++)
            w[n] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)n / (double)(N - 1));
        break;
    }

    /* Apply window in-place */
    for (size_t n = 0; n < N; n++) {
        waveform->samples[n] *= w[n];
    }

    /* Store window coefficients for reference */
    waveform->window_coeffs = w;
    waveform->window_type = win_type;
    return 0;
}

/* ©¤©¤©¤ L5: Barker Code Generation ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Known Barker codes (binary phase codes with minimal autocorrelation sidelobes).
 *
 * Barker codes have the property that all aperiodic autocorrelation
 * sidelobes have magnitude <= 1. The peak-to-sidelobe ratio (in amplitude)
 * equals the code length N.
 *
 * Known lengths and PSL:
 *   N=2:  [+1, -1] or [+1, +1], PSL = -6.0 dB
 *   N=3:  [+1, +1, -1],          PSL = -9.5 dB
 *   N=4:  [+1, +1, -1, +1],      PSL = -12.0 dB
 *   N=5:  [+1, +1, +1, -1, +1],  PSL = -14.0 dB
 *   N=7:  [+1,+1,+1,-1,-1,+1,-1], PSL = -16.9 dB
 *   N=11: [+1,+1,+1,-1,-1,-1,+1,-1,-1,+1,-1], PSL = -20.8 dB
 *   N=13: [+1,+1,+1,+1,+1,-1,-1,+1,+1,-1,+1,-1,+1], PSL = -22.3 dB
 *
 * No Barker codes longer than 13 are known to exist.
 */
int radar_waveform_barker(radar_waveform_t *waveform,
                          int code_length, double pulse_width_s,
                          double sample_rate_hz)
{
    if (!waveform || code_length < 2 || code_length > 13 ||
        pulse_width_s <= 0.0 || sample_rate_hz <= 0.0)
        return -1;

    /* Known Barker codes */
    static const int barker_codes[14][14] = {
        {},
        {},
        {+1, -1},                    /* N=2 */
        {+1, +1, -1},                /* N=3 */
        {+1, +1, -1, +1},            /* N=4 */
        {+1, +1, +1, -1, +1},        /* N=5 */
        {},
        {+1, +1, +1, -1, -1, +1, -1}, /* N=7 */
        {},
        {},
        {},
        {+1,+1,+1,-1,-1,-1,+1,-1,-1,+1,-1}, /* N=11 */
        {},
        {+1,+1,+1,+1,+1,-1,-1,+1,+1,-1,+1,-1,+1} /* N=13 */
    };

    if (code_length != 2 && code_length != 3 && code_length != 4 &&
        code_length != 5 && code_length != 7 && code_length != 11 &&
        code_length != 13)
        return -1;

    memset(waveform, 0, sizeof(*waveform));
    waveform->type = WAVEFORM_BARKER;
    waveform->pulse_width_s = pulse_width_s;
    waveform->code_length = code_length;

    double chip_width_s = pulse_width_s / (double)code_length;
    size_t chips_per_sample = (size_t)(chip_width_s * sample_rate_hz);
    if (chips_per_sample < 1) chips_per_sample = 1;

    waveform->num_samples = (size_t)code_length * chips_per_sample;
    waveform->sample_rate_hz = sample_rate_hz;
    waveform->bandwidth_hz = 1.0 / chip_width_s;
    waveform->chirp_rate_hz_per_s = 0.0;

    waveform->samples = (double complex*)safe_malloc(
        waveform->num_samples * sizeof(double complex));

    for (int c = 0; c < code_length; c++) {
        double phase = (barker_codes[code_length][c] == 1) ? 0.0 : M_PI;
        double complex chip_val = CMPLX(cos(phase), sin(phase));
        for (size_t s = 0; s < chips_per_sample; s++) {
            size_t idx = (size_t)c * chips_per_sample + s;
            if (idx < waveform->num_samples)
                waveform->samples[idx] = chip_val;
        }
    }
    waveform->window_type = WINDOW_RECTANGULAR;
    waveform->window_coeffs = NULL;
    return 0;
}

/* ©¤©¤©¤ L5: Frank Code Generation ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Frank code: N_c x N_c polyphase code.
 *
 * Phase matrix: phi(i, j) = (2*pi/N_c) * i * j,  i,j = 0...N_c-1
 *
 * The code is read row-by-row to produce a sequence of N_c^2 chips.
 * Each chip phase is quantized to N_c levels.
 *
 * The Frank code approximates a stepped-chirp waveform.
 * Its periodic autocorrelation has zero sidelobes (perfect),
 * while its aperiodic autocorrelation has PSL ¡Ö 20*log10(pi/N_c^2).
 *
 * These are frequently used in pulse compression because of their
 * good Doppler tolerance compared to Barker codes.
 */
int radar_waveform_frank(radar_waveform_t *waveform,
                         int nc, double pulse_width_s, double sample_rate_hz)
{
    if (!waveform || nc < 2 || pulse_width_s <= 0.0 || sample_rate_hz <= 0.0)
        return -1;

    int total_chips = nc * nc;
    memset(waveform, 0, sizeof(*waveform));
    waveform->type = WAVEFORM_FRANK;
    waveform->pulse_width_s = pulse_width_s;
    waveform->code_length = total_chips;

    double chip_width_s = pulse_width_s / (double)total_chips;
    size_t chips_per_sample = (size_t)(chip_width_s * sample_rate_hz);
    if (chips_per_sample < 1) chips_per_sample = 1;

    waveform->num_samples = (size_t)total_chips * chips_per_sample;
    waveform->sample_rate_hz = sample_rate_hz;
    waveform->bandwidth_hz = (double)nc / pulse_width_s;
    waveform->chirp_rate_hz_per_s = 0.0;

    waveform->samples = (double complex*)safe_malloc(
        waveform->num_samples * sizeof(double complex));

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            double phase = (2.0 * M_PI / (double)nc) * (double)(i * j);
            double complex val = CMPLX(cos(phase), sin(phase));
            int chip_idx = i * nc + j;
            for (size_t s = 0; s < chips_per_sample; s++) {
                size_t idx = (size_t)chip_idx * chips_per_sample + s;
                if (idx < waveform->num_samples)
                    waveform->samples[idx] = val;
            }
        }
    }
    waveform->window_type = WINDOW_RECTANGULAR;
    waveform->window_coeffs = NULL;
    return 0;
}

/* ©¤©¤©¤ L5: Matched Filter (Time-Domain Correlation) ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Direct time-domain matched filter via cross-correlation.
 *
 * y[k] = sum_{n=0}^{M-1} signal[k+n] * conj(template[n])
 *
 * where template = conj(reversed transmitted waveform).
 *
 * O(N * M) complexity. Use radar_matched_filter_fft() for
 * large time-bandwidth products.
 *
 * Output length = signal_len + template_len - 1.
 * Peak occurs at index corresponding to target range.
 */
int radar_matched_filter(const double complex *signal, size_t signal_len,
                         const double complex *templ, size_t template_len,
                         double complex *output)
{
    if (!signal || !templ || !output || signal_len == 0 || template_len == 0)
        return -1;
    if (template_len > signal_len)
        return -1;

    size_t out_len = signal_len + template_len - 1;

    /* Zero output */
    for (size_t k = 0; k < out_len; k++) {
        output[k] = CMPLX(0.0, 0.0);
    }

    /* Cross-correlation: y[k] = sum signal[k+n] * conj(templ[n]) */
    for (size_t n = 0; n < template_len; n++) {
        double complex conj_tn = conj(templ[n]);
        for (size_t k = 0; k < signal_len; k++) {
            output[k + n] += signal[k] * conj_tn;
        }
    }

    return 0;
}

/* ©¤©¤©¤ L5: FFT-Based Matched Filter ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Compute the next power of 2 greater than or equal to n.
 */
static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

/**
 * In-place radix-2 decimation-in-time FFT.
 *
 * Implements the Cooley-Tukey algorithm.
 * Bit-reversal permutation followed by butterfly stages.
 *
 * Complexity: O(N log N)
 * Reference: Oppenheim & Schafer (2010), Ch.9.
 */
static void fft_radix2(double complex *x, size_t N, int inverse) {
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < N; i++) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j) {
            double complex tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
    }

    /* Butterfly stages */
    double sign = inverse ? 1.0 : -1.0;
    for (size_t len = 2; len <= N; len <<= 1) {
        double ang = 2.0 * M_PI / (double)len * sign;
        double complex wlen = CMPLX(cos(ang), sin(ang));
        for (size_t i = 0; i < N; i += len) {
            double complex w = CMPLX(1.0, 0.0);
            for (size_t j = 0; j < len / 2; j++) {
                double complex u = x[i + j];
                double complex v = x[i + j + len / 2] * w;
                x[i + j] = u + v;
                x[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    /* Scale for inverse FFT */
    if (inverse) {
        for (size_t i = 0; i < N; i++)
            x[i] /= (double)N;
    }
}

/**
 * FFT-based fast convolution matched filter.
 *
 * Algorithm:
 *   1. Zero-pad signal and template to N >= signal_len + template_len - 1
 *   2. N = next_pow2(signal_len + template_len - 1)
 *   3. FFT(signal_padded), FFT(template_padded_reversed)
 *   4. Multiply pointwise: Y[k] = S[k] * conj(T[k])
 *   5. IFFT(Y[k]) ¡ú time-domain matched filter output
 *
 * This is O(N log N) vs O(N*M) for time-domain implementation.
 * Efficient for large time-bandwidth products (TBP > 100).
 * Matches the result of the time-domain convolution to within
 * floating-point precision.
 */
int radar_matched_filter_fft(const double complex *signal, size_t signal_len,
                             const double complex *templ, size_t template_len,
                             double complex *output)
{
    if (!signal || !templ || !output || signal_len == 0 || template_len == 0)
        return -1;

    size_t out_len = signal_len + template_len - 1;
    size_t N = next_pow2(out_len);

    /* Allocate padded arrays */
    double complex *sig_pad = (double complex*)calloc(N, sizeof(double complex));
    double complex *tpl_pad = (double complex*)calloc(N, sizeof(double complex));
    if (!sig_pad || !tpl_pad) {
        free(sig_pad); free(tpl_pad);
        return -1;
    }

    /* Copy and zero-pad signal */
    for (size_t i = 0; i < signal_len; i++) sig_pad[i] = signal[i];

    /* Template: zero-pad with time-reversed, conjugated version */
    for (size_t i = 0; i < template_len; i++)
        tpl_pad[i] = conj(templ[template_len - 1 - i]);

    /* Forward FFTs */
    fft_radix2(sig_pad, N, 0);
    fft_radix2(tpl_pad, N, 0);

    /* Frequency-domain multiplication: Y = S * T */
    for (size_t i = 0; i < N; i++)
        sig_pad[i] *= tpl_pad[i];

    /* Inverse FFT */
    fft_radix2(sig_pad, N, 1);

    /* Copy output (first out_len samples) */
    for (size_t i = 0; i < out_len; i++)
        output[i] = sig_pad[i];

    free(sig_pad);
    free(tpl_pad);
    return 0;
}

/* ©¤©¤©¤ L5: Sidelobe Analysis ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

/**
 * Peak sidelobe level:
 *   PSL = 20 * log10(max_|sidelobe| / |mainlobe|)
 *
 * Searches for the mainlobe peak, then finds the maximum
 * value outside the mainlobe region (within 3 dB of peak).
 */
double radar_psl_compute(const double *compressed, size_t len, size_t *peak_idx)
{
    if (!compressed || len < 3 || !peak_idx) return 0.0;

    /* Find peak */
    size_t pk = 0;
    double pk_val = compressed[0];
    for (size_t i = 1; i < len; i++) {
        if (compressed[i] > pk_val) {
            pk_val = compressed[i];
            pk = i;
        }
    }
    *peak_idx = pk;

    if (pk_val <= 0.0) return 0.0;

    /* Find mainlobe width (3 dB down on each side) */
    double half_power = pk_val / sqrt(2.0);  /* 3 dB = 1/sqrt(2) in amplitude */
    size_t left = pk, right = pk;
    while (left > 0 && compressed[left - 1] > half_power) left--;
    while (right < len - 1 && compressed[right + 1] > half_power) right++;

    /* Expand by one bin on each side for safety */
    if (left > 0) left--;
    if (right < len - 1) right++;

    /* Find max sidelobe outside mainlobe */
    double max_sl = 0.0;
    for (size_t i = 0; i < len; i++) {
        if (i >= left && i <= right) continue;  /* skip mainlobe */
        if (compressed[i] > max_sl) max_sl = compressed[i];
    }

    if (max_sl <= 0.0) return -100.0;  /* essentially no sidelobes */
    return 20.0 * log10(max_sl / pk_val);
}

/**
 * Integrated sidelobe level:
 *   ISL = 10 * log10( sum_|sidelobe|^2 / |mainlobe|^2 )
 *
 * Measures total energy in sidelobes. Lower values mean
 * less energy is spread into sidelobes.
 *
 * For an ideal matched filter with rectangular LFM:
 *   ISL ¡Ö -10 dB (most energy in mainlobe plus close-in sidelobes)
 *
 * With Hamming weighting:
 *   ISL ¡Ö -20 dB
 */
double radar_isl_compute(const double *compressed, size_t len, size_t peak_idx)
{
    if (!compressed || len < 3 || peak_idx >= len) return 0.0;

    double pk_val = compressed[peak_idx];
    if (pk_val <= 0.0) return 0.0;

    /* Find mainlobe region */
    double half_power = pk_val / sqrt(2.0);
    size_t left = peak_idx, right = peak_idx;
    while (left > 0 && compressed[left - 1] > half_power) left--;
    while (right < len - 1 && compressed[right + 1] > half_power) right++;
    if (left > 0) left--;
    if (right < len - 1) right++;

    double mainlobe_pwr = 0.0;
    double sidelobe_pwr = 0.0;

    for (size_t i = 0; i < len; i++) {
        double pwr = compressed[i] * compressed[i];
        if (i >= left && i <= right)
            mainlobe_pwr += pwr;
        else
            sidelobe_pwr += pwr;
    }

    if (mainlobe_pwr <= 0.0) return 0.0;
    return 10.0 * log10(sidelobe_pwr / mainlobe_pwr);
}

/* ©¤©¤©¤ Cleanup ©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤©¤ */

void radar_waveform_free(radar_waveform_t *waveform)
{
    if (!waveform) return;
    free(waveform->samples);
    free(waveform->window_coeffs);
    waveform->samples = NULL;
    waveform->window_coeffs = NULL;
}

/* â”€â”€â”€ L8: Nonlinear FM Waveform â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */

int radar_waveform_nlfm(radar_waveform_t *waveform,
                         double pulse_width_s, double bandwidth_hz,
                         double sample_rate_hz, double alpha)
{
    if (!waveform || pulse_width_s <= 0.0 || bandwidth_hz <= 0.0 ||
        sample_rate_hz <= 0.0 || alpha <= 0.0 || alpha >= M_PI / 2.0)
        return -1;

    memset(waveform, 0, sizeof(*waveform));
    waveform->type = WAVEFORM_NLFM;
    waveform->pulse_width_s = pulse_width_s;
    waveform->bandwidth_hz = bandwidth_hz;
    waveform->sample_rate_hz = sample_rate_hz;
    waveform->num_samples = (size_t)(pulse_width_s * sample_rate_hz);
    if (waveform->num_samples < 1) waveform->num_samples = 1;

    waveform->samples = (double complex*)safe_malloc(
        waveform->num_samples * sizeof(double complex));

    double tau = pulse_width_s;
    double B = bandwidth_hz;
    double ts = 1.0 / sample_rate_hz;

    for (size_t n = 0; n < waveform->num_samples; n++) {
        double t_norm = 2.0 * (double)n * ts / tau - 1.0;
        double phi = M_PI * B * t_norm / 2.0;
        waveform->samples[n] = CMPLX(cos(phi), sin(phi));
    }
    waveform->window_type = WINDOW_RECTANGULAR;
    waveform->window_coeffs = NULL;
    return 0;
}

/* â”€â”€â”€ L7: Coherent Pulse Train â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */

size_t radar_pulse_train_generate(const radar_waveform_t *base_wf,
                                   size_t n_pulses, double prf_hz,
                                   double complex *train, size_t train_max_len)
{
    if (!base_wf || !base_wf->samples || !train || n_pulses == 0 ||
        prf_hz <= 0.0 || train_max_len == 0) return 0;

    size_t samples_per_pri = (size_t)(base_wf->sample_rate_hz / prf_hz);
    if (samples_per_pri < base_wf->num_samples)
        samples_per_pri = base_wf->num_samples + 1;

    size_t total = n_pulses * samples_per_pri;
    if (total > train_max_len) return 0;

    for (size_t i = 0; i < total; i++) train[i] = CMPLX(0.0, 0.0);
    for (size_t p = 0; p < n_pulses; p++) {
        size_t off = p * samples_per_pri;
        for (size_t s = 0; s < base_wf->num_samples; s++)
            train[off + s] = base_wf->samples[s];
    }
    return total;
}

/* â”€â”€â”€ L5: Stepped-Frequency Waveform â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */

int radar_waveform_stepped_freq(radar_waveform_t *waveform,
                                 size_t n_steps, double delta_f_hz,
                                 double pulse_width_s, double sample_rate_hz)
{
    if (!waveform || n_steps < 2 || delta_f_hz <= 0.0 ||
        pulse_width_s <= 0.0 || sample_rate_hz <= 0.0) return -1;
    memset(waveform, 0, sizeof(*waveform));
    waveform->type = WAVEFORM_STEPPED_FREQ;
    waveform->pulse_width_s = pulse_width_s;
    waveform->bandwidth_hz = (double)n_steps * delta_f_hz;
    waveform->sample_rate_hz = sample_rate_hz;
    waveform->num_samples = (size_t)(pulse_width_s * sample_rate_hz);
    if (waveform->num_samples < 1) waveform->num_samples = 1;
    waveform->code_length = (int)n_steps;
    waveform->samples = (double complex*)safe_malloc(
        waveform->num_samples * sizeof(double complex));
    for (size_t n = 0; n < waveform->num_samples; n++)
        waveform->samples[n] = CMPLX(1.0, 0.0);
    waveform->window_type = WINDOW_RECTANGULAR;
    waveform->window_coeffs = NULL;
    return 0;
}

/* â”€â”€â”€ L7: Waveform Auto-Correlation â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */

int radar_waveform_autocorr(const double complex *wf, size_t len,
                             double complex *acorr, size_t max_lag)
{
    if (!wf || !acorr || len == 0 || max_lag > len) return -1;
    for (size_t k = 0; k < max_lag; k++) {
        double complex sum = CMPLX(0.0, 0.0);
        for (size_t n = 0; n < len - k; n++)
            sum += wf[n] * conj(wf[n + k]);
        acorr[k] = sum;
    }
    return 0;
}

/* â”€â”€â”€ L8: Ambiguity Function Sample â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */

double radar_ambiguity_sample(const double complex *wf, size_t len,
                               size_t delay_bin, double doppler_hz,
                               double fs_hz)
{
    if (!wf || len == 0 || delay_bin >= len || fs_hz <= 0.0) return 0.0;
    double complex sum = CMPLX(0.0, 0.0);
    double energy = 0.0;
    for (size_t n = delay_bin; n < len; n++) {
        double phase = -2.0 * M_PI * doppler_hz * (double)n / fs_hz;
        double complex dopp = CMPLX(cos(phase), sin(phase));
        sum += wf[n] * conj(wf[n - delay_bin]) * dopp;
        energy += cabs(wf[n]) * cabs(wf[n]);
    }
    if (energy <= 0.0) return 0.0;
    return cabs(sum) / energy;
}
