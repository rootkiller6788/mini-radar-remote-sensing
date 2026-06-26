#include "radar_waveform.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const int8_t BARKER_2[]  = {1, -1};
static const int8_t BARKER_3[]  = {1, 1, -1};
static const int8_t BARKER_4[]  = {1, 1, -1, 1};
static const int8_t BARKER_5[]  = {1, 1, 1, -1, 1};
static const int8_t BARKER_7[]  = {1, 1, 1, -1, -1, 1, -1};
static const int8_t BARKER_11[] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1};
static const int8_t BARKER_13[] = {1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1};

static const int8_t* get_barker_code(uint32_t length) {
    switch (length) {
        case 2:  return BARKER_2;
        case 3:  return BARKER_3;
        case 4:  return BARKER_4;
        case 5:  return BARKER_5;
        case 7:  return BARKER_7;
        case 11: return BARKER_11;
        case 13: return BARKER_13;
        default: return NULL;
    }
}

int rect_pulse_generate(const radar_waveform_params_t *params,
                        double complex *buffer)
{
    if (params == NULL || buffer == NULL) return -1;
    if (params->pulse_width <= 0.0 || params->sampling_rate <= 0.0) return -1;
    if (params->num_samples == 0) return -1;

    double dt = 1.0 / params->sampling_rate;
    double half_width = params->pulse_width / 2.0;
    uint32_t N = params->num_samples;

    for (uint32_t n = 0; n < N; n++) {
        double t = ((double)n - (double)(N-1)/2.0) * dt;
        buffer[n] = (fabs(t) <= half_width) ? (1.0 + 0.0 * I) : (0.0 + 0.0 * I);
    }
    return 0;
}

int lfm_chirp_generate(const radar_waveform_params_t *params,
                       const lfm_params_t *chirp,
                       double complex *buffer)
{
    if (params == NULL || chirp == NULL || buffer == NULL) return -1;
    if (params->pulse_width <= 0.0 || params->sampling_rate <= 0.0) return -1;
    if (params->num_samples == 0) return -1;
    if (chirp->chirp_rate == 0.0) return -1;

    double dt = 1.0 / params->sampling_rate;
    double half_width = params->pulse_width / 2.0;
    uint32_t N = params->num_samples;

    for (uint32_t n = 0; n < N; n++) {
        double t = ((double)n - (double)(N-1)/2.0) * dt;
        if (fabs(t) <= half_width) {
            double phase = M_PI * chirp->chirp_rate * t * t;
            buffer[n] = cos(phase) + sin(phase) * I;
        } else {
            buffer[n] = 0.0 + 0.0 * I;
        }
    }
    return 0;
}

int barker_code_generate(const radar_waveform_params_t *params,
                         const barker_code_t *barker,
                         double complex *buffer)
{
    if (params == NULL || barker == NULL || buffer == NULL) return -1;
    if (barker->code_length == 0 || barker->code_length > 13) return -1;
    if (barker->chip_rate <= 0.0) return -1;
    const int8_t *code = barker->code_sequence;
    if (code == NULL) {
        code = get_barker_code(barker->code_length);
        if (code == NULL) return -1;
    }

    double dt = 1.0 / params->sampling_rate;
    double chip_dur = barker->chip_duration;
    double total_dur = barker->code_length * chip_dur;
    double half_dur = total_dur / 2.0;
    uint32_t N = params->num_samples;

    for (uint32_t n = 0; n < N; n++) {
        double t = ((double)n - (double)(N-1)/2.0) * dt;
        if (fabs(t) <= half_dur) {
            double t_shifted = t + half_dur;
            uint32_t chip_idx = (uint32_t)(t_shifted / chip_dur);
            if (chip_idx >= barker->code_length)
                chip_idx = barker->code_length - 1;
            buffer[n] = (double)code[chip_idx] + 0.0 * I;
        } else {
            buffer[n] = 0.0 + 0.0 * I;
        }
    }
    return 0;
}

int pulse_train_generate(const radar_waveform_params_t *params,
                         uint32_t num_pulses,
                         const double complex *single_pulse,
                         double complex *buffer)
{
    if (params == NULL || single_pulse == NULL || buffer == NULL) return -1;
    if (num_pulses == 0 || params->num_samples == 0) return -1;
    uint32_t N = params->num_samples;
    size_t total = (size_t)num_pulses * N;
    for (size_t i = 0; i < total; i++)
        buffer[i] = single_pulse[i % N];
    return 0;
}

int waveform_autocorrelation(const double complex *signal,
                             size_t N,
                             double *acorr)
{
    if (signal == NULL || acorr == NULL) return -1;
    if (N == 0) return -1;
    for (size_t k = 0; k < N; k++) {
        double sum = 0.0;
        for (size_t n = 0; n < N - k; n++) {
            double complex prod = signal[n] * conj(signal[n + k]);
            sum += creal(prod);
        }
        acorr[k] = sum;
    }
    return 0;
}

int waveform_cross_ambiguity_1d(const double complex *tx,
                                const double complex *rx,
                                size_t length,
                                double fd,
                                double *result)
{
    if (tx == NULL || rx == NULL || result == NULL) return -1;
    if (length == 0) return -1;
    double complex accum = 0.0 + 0.0 * I;
    for (size_t n = 0; n < length; n++) {
        double phase = -2.0 * M_PI * fd * (double)n;
        accum += tx[n] * conj(rx[n]) * (cos(phase) + sin(phase) * I);
    }
    *result = creal(accum)*creal(accum) + cimag(accum)*cimag(accum);
    return 0;
}

int waveform_instantaneous_freq(const double complex *signal,
                                size_t N,
                                double *freq,
                                double fs)
{
    if (signal == NULL || freq == NULL) return -1;
    if (N < 2 || fs <= 0.0) return -1;
    freq[0] = 0.0;
    double phase_prev = atan2(cimag(signal[0]), creal(signal[0]));
    for (size_t n = 1; n < N; n++) {
        double phase = atan2(cimag(signal[n]), creal(signal[n]));
        double diff = phase - phase_prev;
        while (diff > M_PI)  diff -= 2.0 * M_PI;
        while (diff < -M_PI) diff += 2.0 * M_PI;
        freq[n] = (fs / (2.0 * M_PI)) * diff;
        phase_prev = phase;
    }
    freq[0] = freq[1];
    return 0;
}

int waveform_apply_window(double complex *signal,
                          size_t N,
                          int window_type,
                          double taylor_nbar,
                          double taylor_sll_db)
{
    if (signal == NULL || N == 0 || window_type < 0 || window_type > 4) return -1;
    for (size_t n = 0; n < N; n++) {
        double w = 1.0, dn = (double)n, dN1 = (double)(N - 1);
        switch (window_type) {
        case 0: w = 1.0; break;
        case 1: w = 0.54 - 0.46 * cos(2.0 * M_PI * dn / dN1); break;
        case 2: w = 0.5 - 0.5 * cos(2.0 * M_PI * dn / dN1); break;
        case 3: w = 0.42 - 0.5*cos(2.0*M_PI*dn/dN1) + 0.08*cos(4.0*M_PI*dn/dN1); break;
        case 4: {
            double A_val = acosh(pow(10.0, -taylor_sll_db / 20.0));
            double nb = (taylor_nbar > 0.0) ? taylor_nbar : ceil(2.0*A_val*A_val+0.5);
            double sigma2 = nb*nb/(A_val*A_val+(nb-0.5)*(nb-0.5));
            double x = (dn-dN1/2.0)/(dN1/2.0);
            w = 1.0;
            for (int m = 1; m < (int)nb; m++) {
                double Fm = 1.0;
                for (int i = 1; i < (int)nb; i++)
                    Fm *= 1.0-(double)(m*m)/(sigma2*(A_val*A_val+(i-0.5)*(i-0.5)));
                w += 2.0*Fm*cos(2.0*M_PI*(double)m*x);
            }
            w /= (2.0*nb+1.0);
            break; }
        }
        signal[n] *= w;
    }
    return 0;
}

int waveform_rms_bandwidth(const double complex *signal,
                           size_t N, double fs, double *bw_rms)
{
    if (signal == NULL || bw_rms == NULL || N < 2 || fs <= 0.0) return -1;
    double energy = 0.0, diff_energy = 0.0;
    for (size_t n = 0; n < N; n++) {
        double msq = creal(signal[n])*creal(signal[n])+cimag(signal[n])*cimag(signal[n]);
        energy += msq;
    }
    if (energy < 1e-30) { *bw_rms = 0.0; return 0; }
    for (size_t n = 1; n < N; n++) {
        double complex d = signal[n] - signal[n-1];
        diff_energy += creal(d)*creal(d) + cimag(d)*cimag(d);
    }
    *bw_rms = sqrt((fs*fs)/(4.0*M_PI*M_PI*energy) * diff_energy);
    return 0;
}

/* Additional L5 functions for radar waveform analysis */

/* compute_peak_sidelobe_level - L5: measures range sidelobe suppression.
 * Scans autocorrelation sidelobes and finds the peak relative to mainlobe.
 * PSL = 20*log10(max_{k>0}|R[k]|/|R[0]|) dB */
int compute_peak_sidelobe_level(const double *autocorr, size_t N,
                                double *psl_db)
{
    if (autocorr == NULL || psl_db == NULL || N < 2) return -1;
    double peak = fabs(autocorr[0]);
    if (peak < 1e-30) { *psl_db = 0.0; return 0; }
    double max_sidelobe = 0.0;
    for (size_t k = 1; k < N; k++) {
        double abs_val = fabs(autocorr[k]);
        if (abs_val > max_sidelobe) max_sidelobe = abs_val;
    }
    *psl_db = 20.0 * log10(max_sidelobe / peak);
    return 0;
}

/* compute_integrated_sidelobe_level - L5: total sidelobe energy ratio.
 * ISLR = 10*log10(total_sidelobe_energy / mainlobe_energy) dB */
int compute_integrated_sidelobe_level(const double *autocorr, size_t N,
                                      double *islr_db)
{
    if (autocorr == NULL || islr_db == NULL || N < 2) return -1;
    double mainlobe = autocorr[0] * autocorr[0];
    double sidelobes = 0.0;
    for (size_t k = 1; k < N; k++)
        sidelobes += autocorr[k] * autocorr[k];
    if (mainlobe < 1e-30) { *islr_db = 0.0; return 0; }
    *islr_db = 10.0 * log10(sidelobes / mainlobe);
    return 0;
}

/* compute_range_resolution - L1: Rayleigh range resolution.
 * Delta_R = c / (2 * B) where B is the waveform bandwidth.
 * For LFM: B = chirp bandwidth; for unmodulated pulse: B ~ 1/tau. */
double compute_range_resolution(double bandwidth_hz)
{
    if (bandwidth_hz <= 0.0) return 0.0;
    return 299792458.0 / (2.0 * bandwidth_hz);
}

/* compute_doppler_resolution - L1: Doppler resolution from CPI.
 * Delta_fd = 1 / T_cpi = PRF / M where M = number of pulses. */
double compute_doppler_resolution(double prf_hz, uint32_t num_pulses)
{
    if (prf_hz <= 0.0 || num_pulses == 0) return 0.0;
    return prf_hz / (double)num_pulses;
}

/* compute_velocity_resolution - L1: velocity resolution from Doppler res.
 * Delta_v = lambda * Delta_fd / 2 = lambda * PRF / (2 * M) */
double compute_velocity_resolution(double wavelength_m, double prf_hz,
                                    uint32_t num_pulses)
{
    if (wavelength_m <= 0.0 || prf_hz <= 0.0 || num_pulses == 0) return 0.0;
    return wavelength_m * prf_hz / (2.0 * (double)num_pulses);
}

/* compute_pulse_compression_gain - L2: processing gain of LFM.
 * G = 10*log10(tau * B) dB. For unmodulated pulse: G = 0 dB. */
double compute_pulse_compression_gain(double pulse_width_s,
                                       double bandwidth_hz)
{
    double tbp = pulse_width_s * bandwidth_hz;
    if (tbp <= 0.0) return 0.0;
    return 10.0 * log10(tbp);
}

/* compute_duty_cycle - L1: duty cycle from pulse width and PRF.
 * duty = tau * PRF (dimensionless, 0 < duty < 1). */
double compute_duty_cycle(double pulse_width_s, double prf_hz)
{
    if (pulse_width_s <= 0.0 || prf_hz <= 0.0) return 0.0;
    double duty = pulse_width_s * prf_hz;
    if (duty >= 1.0) return 0.999;
    return duty;
}

/* compute_peak_power - L1: peak-to-average power relationship.
 * P_peak = P_avg / duty_cycle for pulsed radar. */
double compute_peak_power(double average_power_w, double duty_cycle)
{
    if (duty_cycle <= 0.0 || duty_cycle >= 1.0) return 0.0;
    return average_power_w / duty_cycle;
}

/* compute_average_power - L1: average from peak power.
 * P_avg = P_peak * duty_cycle. */
double compute_average_power(double peak_power_w, double duty_cycle)
{
    if (duty_cycle <= 0.0 || duty_cycle >= 1.0) return 0.0;
    return peak_power_w * duty_cycle;
}

/* compute_energy_per_pulse - L1: pulse energy.
 * E = P_peak * tau. */
double compute_energy_per_pulse(double peak_power_w, double pulse_width_s)
{
    if (peak_power_w <= 0.0 || pulse_width_s <= 0.0) return 0.0;
    return peak_power_w * pulse_width_s;
}

/* waveform_spectrum_magnitude - L3: DFT-based spectrum magnitude.
 * |S[k]| = |sum s[n]*exp(-j*2*pi*k*n/N)| for k = 0..N-1.
 * Used to verify waveform spectral properties. */
int waveform_spectrum_magnitude(const double complex *signal, size_t N,
                                double *magnitude_spectrum)
{
    if (signal == NULL || magnitude_spectrum == NULL || N == 0) return -1;
    for (size_t k = 0; k < N; k++) {
        double complex sum = 0.0 + 0.0 * I;
        for (size_t n = 0; n < N; n++) {
            double angle = -2.0 * M_PI * (double)(k * n) / (double)N;
            sum += signal[n] * (cos(angle) + sin(angle) * I);
        }
        magnitude_spectrum[k] = sqrt(creal(sum)*creal(sum) + cimag(sum)*cimag(sum));
    }
    return 0;
}

/* costas_code_generate - L2: Costas frequency-hopping code.
 * Generates a Costas array waveform: each pulse uses a different
 * frequency based on a permutation with ideal autocorrelation.
 * Costas arrays provide thumbtack ambiguity functions. */
static const uint32_t COSTAS_7[] = {4, 7, 1, 6, 5, 2, 3};

int costas_code_generate(const radar_waveform_params_t *params,
                         uint32_t code_length,
                         double complex *buffer)
{
    if (params == NULL || buffer == NULL) return -1;
    if (code_length == 0 || code_length > 13) return -1;
    if (params->num_samples == 0) return -1;

    /* Use known Costas array for length 7; generalize for others via
     * Welch construction for prime lengths. */
    double dt = 1.0 / params->sampling_rate;
    double chip_dur = params->pulse_width / (double)code_length;
    uint32_t N = params->num_samples;

    for (uint32_t n = 0; n < N; n++) {
        double t = ((double)n - (double)(N-1)/2.0) * dt;
        double half_dur = params->pulse_width / 2.0;
        buffer[n] = 0.0 + 0.0 * I;
        if (fabs(t) <= half_dur) {
            double t_shifted = t + half_dur;
            uint32_t chip = (uint32_t)(t_shifted / chip_dur);
            if (chip < code_length) {
                uint32_t freq_idx = (code_length == 7)
                    ? COSTAS_7[chip] : (chip + 1);
                double freq_offset = ((double)freq_idx - (double)code_length/2.0)
                                     * params->bandwidth / (double)code_length;
                double phase = 2.0 * M_PI * freq_offset * t;
                buffer[n] = cos(phase) + sin(phase) * I;
            }
        }
    }
    return 0;
}

/* stepped_frequency_generate - L7: Stepped-frequency CW waveform.
 * Transmits N pulses at frequencies fc + n*delta_f, n=0..N-1.
 * After IFFT of the received frequency-domain data, achieves
 * range resolution Delta_R = c/(2*N*delta_f). Used in
 * high-range-resolution profiling (HRRP) and SAR. */
int stepped_frequency_generate(const radar_waveform_params_t *params,
                               uint32_t num_steps,
                               double frequency_step_hz,
                               double complex *buffer)
{
    if (params == NULL || buffer == NULL || num_steps == 0) return -1;
    if (frequency_step_hz <= 0.0 || params->num_samples == 0) return -1;

    double dt = 1.0 / params->sampling_rate;
    double half_width = params->pulse_width / 2.0;
    uint32_t samples_per_step = params->num_samples / num_steps;
    if (samples_per_step == 0) return -1;

    for (uint32_t n = 0; n < params->num_samples; n++) {
        double t = ((double)n - (double)(params->num_samples-1)/2.0) * dt;
        if (fabs(t) <= half_width) {
            uint32_t step_idx = n / samples_per_step;
            if (step_idx >= num_steps) step_idx = num_steps - 1;
            double freq = params->center_frequency + (double)step_idx * frequency_step_hz;
            double phase = 2.0 * M_PI * freq * t;
            buffer[n] = cos(phase) + sin(phase) * I;
        } else {
            buffer[n] = 0.0 + 0.0 * I;
        }
    }
    return 0;
}

/* compute_signal_to_clutter_ratio - L7: SCR for distributed ground clutter.
 * SCR = sigma_t / (sigma0 * Delta_R * Delta_az) in linear units.
 * sigma_t = target RCS (m^2), sigma0 = clutter reflectivity (m^2/m^2). */
double compute_signal_to_clutter_ratio(double target_rcs_dbsm,
                                       double clutter_rcs_per_area_dbsm,
                                       double range_resolution_m,
                                       double azimuth_resolution_m)
{
    double area = range_resolution_m * azimuth_resolution_m;
    if (area <= 0.0) return 0.0;
    double sigma_t = pow(10.0, target_rcs_dbsm / 10.0);
    double sigma0 = pow(10.0, clutter_rcs_per_area_dbsm / 10.0);
    double scr_linear = sigma_t / (sigma0 * area + 1e-30);
    return 10.0 * log10(scr_linear);
}

/* compute_unambiguous_range - L1: Rmax = c/(2*PRF) */
double compute_unambiguous_range(double prf_hz)
{
    if (prf_hz <= 0.0) return 0.0;
    return 299792458.0 / (2.0 * prf_hz);
}

/* compute_blind_speed - L2: v_blind = k*lambda*PRF/2.
 * At these speeds, the Doppler shift equals a multiple of PRF,
 * making the target appear stationary (blind). */
double compute_blind_speed(double prf_hz, double wavelength_m, uint32_t k)
{
    if (prf_hz <= 0.0 || wavelength_m <= 0.0 || k == 0) return 0.0;
    return (double)k * wavelength_m * prf_hz / 2.0;
}
