#include "doppler_processing.h"
#include "pulse_doppler.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int window_function_create(window_type_t type, size_t length,
                           double param, window_func_t *win)
{
    if (win == NULL || length == 0) return -1;
    win->type = type;
    win->length = length;
    win->coefficients = (double *)malloc(length * sizeof(double));
    if (win->coefficients == NULL) return -1;

    double N = (double)length;
    double coherent_sum = 0.0;
    double max_val = 0.0;

    for (size_t n = 0; n < length; n++) {
        double w = 1.0;
        switch (type) {
        case WINDOW_RECTANGULAR:
            w = 1.0; break;
        case WINDOW_HAMMING:
            w = 0.54 - 0.46 * cos(2.0*M_PI*(double)n/(N-1.0)); break;
        case WINDOW_HANNING:
            w = 0.5 - 0.5 * cos(2.0*M_PI*(double)n/(N-1.0)); break;
        case WINDOW_BLACKMAN:
            w = 0.42 - 0.5*cos(2.0*M_PI*(double)n/(N-1.0))
                   + 0.08*cos(4.0*M_PI*(double)n/(N-1.0)); break;
        case WINDOW_BLACKMAN_HARRIS:
            w = 0.35875 - 0.48829*cos(2.0*M_PI*(double)n/(N-1.0))
                   + 0.14128*cos(4.0*M_PI*(double)n/(N-1.0))
                   - 0.01168*cos(6.0*M_PI*(double)n/(N-1.0)); break;
        case WINDOW_FLATTOP:
            w = 0.21557895 - 0.41663158*cos(2.0*M_PI*(double)n/(N-1.0))
                   + 0.277263158*cos(4.0*M_PI*(double)n/(N-1.0))
                   - 0.083578947*cos(6.0*M_PI*(double)n/(N-1.0))
                   + 0.006947368*cos(8.0*M_PI*(double)n/(N-1.0)); break;
        case WINDOW_CHEBYSHEV:
        {   double alpha = cosh(acosh(pow(10.0, param/20.0))/(N-1.0));
            w = 0.0;
            for (size_t k = 0; k <= length/2; k++) {
                double theta = M_PI*(double)k/(N-1.0);
                if (theta <= acos(1.0/alpha)) {
                    double t = alpha*cos(theta);
                    w += cos(2.0*(double)k*M_PI*(double)n/(N-1.0))*cosh((N-1.0)*acosh(t));
                }
            }
            w = fabs(w);
            if (w > max_val) max_val = w;
            break; }
        case WINDOW_KAISER:
        {   double beta = M_PI * param;
            double denom = 0.0;
            double x = 2.0*(double)n/(N-1.0) - 1.0;
            if (fabs(x) < 1.0) {
                double arg = beta*sqrt(1.0-x*x);
                double term = 1.0, sum = 1.0;
                for (int k = 1; k < 20; k++) {
                    term *= (arg*arg)/(4.0*(double)k*(double)k);
                    sum += term;
                }
                w = sum;
            } else { w = 1.0; }
            denom = 1.0; { double t2=1.0,s2=1.0;
            for(int k=1;k<20;k++){t2*=(beta*beta)/(4.0*k*k);s2+=t2;}denom=s2;}
            w /= denom;
            break; }
        case WINDOW_TAYLOR:
        {   double A_val = acosh(pow(10.0, param/20.0));
            double nb = ceil(2.0*A_val*A_val+0.5);
            double sigma2 = nb*nb/(A_val*A_val+(nb-0.5)*(nb-0.5));
            double x = ((double)n-(N-1.0)/2.0)/((N-1.0)/2.0);
            w = 1.0;
            for (int m = 1; m < (int)nb; m++) {
                double Fm = 1.0;
                for (int i = 1; i < (int)nb; i++)
                    Fm *= 1.0-(double)(m*m)/(sigma2*(A_val*A_val+(i-0.5)*(i-0.5)));
                w += 2.0*Fm*cos(2.0*M_PI*(double)m*x);
            }
            w /= (2.0*nb+1.0);
            break; }
        default: free(win->coefficients); return -1;
        }
        win->coefficients[n] = w;
        coherent_sum += w;
    }

    if (type == WINDOW_CHEBYSHEV && max_val > 0.0)
        for (size_t n = 0; n < length; n++)
            win->coefficients[n] /= max_val;

    win->coherent_gain = coherent_sum / N;
    win->peak_sidelobe_db = (type == WINDOW_RECTANGULAR) ? -13.3 :
        (type == WINDOW_HAMMING) ? -42.7 : (type == WINDOW_HANNING) ? -31.5 :
        (type == WINDOW_BLACKMAN) ? -58.1 : (type == WINDOW_BLACKMAN_HARRIS) ? -92.0 :
        (type == WINDOW_FLATTOP) ? -93.0 : -param;

    if (type == WINDOW_RECTANGULAR) {
        win->_3db_bandwidth_bins = 0.89; win->_6db_bandwidth_bins = 1.21;
    } else if (type == WINDOW_HAMMING) {
        win->_3db_bandwidth_bins = 1.30; win->_6db_bandwidth_bins = 1.81;
    } else if (type == WINDOW_HANNING) {
        win->_3db_bandwidth_bins = 1.44; win->_6db_bandwidth_bins = 2.00;
    } else if (type == WINDOW_BLACKMAN) {
        win->_3db_bandwidth_bins = 1.68; win->_6db_bandwidth_bins = 2.35;
    } else {
        win->_3db_bandwidth_bins = 1.5; win->_6db_bandwidth_bins = 2.0;
    }
    win->scalloping_loss = 20.0*log10(1.0/fabs(coherent_sum/N));
    win->worst_case_processing_loss = -20.0*log10(win->coherent_gain);
    return 0;
}

void window_function_free(window_func_t *win) {
    if (win != NULL && win->coefficients != NULL) {
        free(win->coefficients);
        win->coefficients = NULL;
    }
}

int window_function_apply(const window_func_t *win,
                          double complex *data, size_t length) {
    if (win == NULL || data == NULL || win->coefficients == NULL) return -1;
    if (length != win->length) return -1;
    for (size_t n = 0; n < length; n++)
        data[n] *= win->coefficients[n];
    return 0;
}

int doppler_fft_process(const double complex *pulse_matrix,
                        uint32_t num_range_bins, uint32_t num_pulses,
                        const doppler_fft_params_t *params,
                        range_doppler_map_t *rdmap)
{
    if (pulse_matrix == NULL || params == NULL || rdmap == NULL) return -1;
    if (rdmap->data == NULL) return -1;

    uint32_t fft_size = params->fft_size;
    window_func_t win;
    if (params->window != WINDOW_RECTANGULAR) {
        window_function_create(params->window, num_pulses, 60.0, &win);
    }

    double complex *slow_time = (double complex *)malloc(fft_size * sizeof(double complex));
    if (slow_time == NULL) return -1;

    for (uint32_t r = 0; r < num_range_bins; r++) {
        for (uint32_t p = 0; p < fft_size; p++) {
            if (p < num_pulses) {
                slow_time[p] = pulse_matrix[r + p * num_range_bins];
                if (params->window != WINDOW_RECTANGULAR)
                    slow_time[p] *= win.coefficients[p];
            } else {
                slow_time[p] = 0.0 + 0.0 * I;
            }
        }
        for (uint32_t k = 0; k < fft_size; k++) {
            double complex sum = 0.0 + 0.0 * I;
            for (uint32_t n = 0; n < fft_size; n++) {
                double angle = -2.0 * M_PI * (double)(k * n) / (double)fft_size;
                sum += slow_time[n] * (cos(angle) + sin(angle) * I);
            }
            rdmap->data[r + k * num_range_bins] = sum;
        }
    }

    if (params->window != WINDOW_RECTANGULAR) window_function_free(&win);
    free(slow_time);

    rdmap->num_range_bins = num_range_bins;
    rdmap->num_doppler_bins = fft_size;
    rdmap->doppler_resolution_hz = params->doppler_resolution_hz;
    rdmap->velocity_resolution_mps = params->velocity_resolution_mps;
    rdmap->max_unambiguous_velocity_mps = params->max_velocity_mps;
    return 0;
}

int doppler_spectrum_compute(const double complex *slow_time_data,
                             uint32_t num_pulses, double prf_hz,
                             double wavelength_m, window_type_t window,
                             doppler_spectrum_t *spectrum)
{
    if (slow_time_data == NULL || spectrum == NULL || num_pulses == 0) return -1;

    spectrum->num_bins = num_pulses;
    spectrum->velocity_resolution = wavelength_m * prf_hz / (2.0 * num_pulses);
    spectrum->max_velocity = wavelength_m * prf_hz / 4.0;

    spectrum->velocity_bins = (double *)malloc(num_pulses * sizeof(double));
    spectrum->magnitude_spectrum = (double *)malloc(num_pulses * sizeof(double));
    spectrum->phase_spectrum = (double *)malloc(num_pulses * sizeof(double));
    if (!spectrum->velocity_bins || !spectrum->magnitude_spectrum || !spectrum->phase_spectrum) {
        free(spectrum->velocity_bins); free(spectrum->magnitude_spectrum);
        free(spectrum->phase_spectrum); return -1;
    }

    window_func_t win;
    window_function_create(window, num_pulses, 60.0, &win);

    double complex *windowed = (double complex *)malloc(num_pulses * sizeof(double complex));
    for (uint32_t n = 0; n < num_pulses; n++) windowed[n] = slow_time_data[n] * win.coefficients[n];

    for (uint32_t k = 0; k < num_pulses; k++) {
        double complex sum = 0.0 + 0.0 * I;
        for (uint32_t n = 0; n < num_pulses; n++) {
            double angle = -2.0*M_PI*(double)(k*n)/(double)num_pulses;
            sum += windowed[n] * (cos(angle) + sin(angle)*I);
        }
        double mag = sqrt(creal(sum)*creal(sum)+cimag(sum)*cimag(sum));
        spectrum->magnitude_spectrum[k] = mag;
        spectrum->phase_spectrum[k] = atan2(cimag(sum), creal(sum));

        double freq_hz;
        if (k <= num_pulses/2) freq_hz = (double)k * prf_hz / (double)num_pulses;
        else freq_hz = ((double)k - (double)num_pulses) * prf_hz / (double)num_pulses;
        spectrum->velocity_bins[k] = freq_hz * wavelength_m / 2.0;
    }

    window_function_free(&win);
    free(windowed);
    return 0;
}

void doppler_spectrum_free(doppler_spectrum_t *spectrum) {
    if (spectrum != NULL) {
        free(spectrum->velocity_bins);
        free(spectrum->magnitude_spectrum);
        free(spectrum->phase_spectrum);
        spectrum->velocity_bins = NULL;
        spectrum->magnitude_spectrum = NULL;
        spectrum->phase_spectrum = NULL;
    }
}

int velocity_from_doppler(double doppler_hz, double wavelength_m,
                          double *velocity_mps) {
    if (velocity_mps == NULL) return -1;
    *velocity_mps = doppler_hz * wavelength_m / 2.0;
    return 0;
}

int doppler_from_velocity(double velocity_mps, double wavelength_m,
                          double *doppler_hz) {
    if (doppler_hz == NULL) return -1;
    *doppler_hz = 2.0 * velocity_mps / wavelength_m;
    return 0;
}

int max_unambiguous_velocity(double prf_hz, double wavelength_m,
                             double *vmax_mps) {
    if (vmax_mps == NULL || prf_hz <= 0.0 || wavelength_m <= 0.0) return -1;
    *vmax_mps = wavelength_m * prf_hz / 4.0;
    return 0;
}

int blind_speed_compute(double prf_hz, double wavelength_m,
                        uint32_t k, double *blind_speed_mps) {
    if (blind_speed_mps == NULL || prf_hz <= 0.0 || k == 0) return -1;
    *blind_speed_mps = (double)k * wavelength_m * prf_hz / 2.0;
    return 0;
}

int stagger_prf_design(double wavelength_m, double max_velocity_mps,
                       uint32_t num_prfs, staggered_prf_t *stagger) {
    if (stagger == NULL || num_prfs < 2 || max_velocity_mps <= 0.0) return -1;
    stagger->num_prfs = num_prfs;
    stagger->prf_values = (double *)malloc(num_prfs * sizeof(double));
    stagger->velocity_estimates = (double *)malloc(num_prfs * sizeof(double));
    if (!stagger->prf_values || !stagger->velocity_estimates) {
        free(stagger->prf_values); free(stagger->velocity_estimates); return -1;
    }
    double base_prf = 4.0 * max_velocity_mps / wavelength_m;
    for (uint32_t i = 0; i < num_prfs; i++) {
        stagger->prf_values[i] = base_prf * (1.0 + 0.15 * (double)i);
        stagger->velocity_estimates[i] = wavelength_m * stagger->prf_values[i] / 4.0;
    }
    stagger->unambiguous_velocity = stagger->velocity_estimates[num_prfs-1];
    return 0;
}

void stagger_prf_free(staggered_prf_t *stagger) {
    if (stagger != NULL) {
        free(stagger->prf_values);
        free(stagger->velocity_estimates);
    }
}

int resolve_doppler_ambiguity_crt(const double *ambiguous_frequencies,
                                  const double *prf_values,
                                  uint32_t num_prfs,
                                  double *true_frequency_hz) {
    if (ambiguous_frequencies == NULL || prf_values == NULL
        || true_frequency_hz == NULL || num_prfs < 2) return -1;

    double M = 1.0;
    for (uint32_t i = 0; i < num_prfs; i++) M *= prf_values[i];

    double result = 0.0;
    for (uint32_t i = 0; i < num_prfs; i++) {
        double Mi = M / prf_values[i];
        double inv = 1.0;
        for (double t = 1.0; t < prf_values[i]; t += 1.0) {
            if (fmod(Mi * t, prf_values[i]) < 1e-6 || fmod(Mi * t, prf_values[i]) > prf_values[i]-1e-6) {
                inv = t; break;
            }
        }
        result += ambiguous_frequencies[i] * Mi * inv;
    }
    *true_frequency_hz = fmod(result, M);
    if (*true_frequency_hz > M/2.0) *true_frequency_hz -= M;
    return 0;
}

int mti_filter_design(mti_filter_type_t type, double prf_hz,
                      double clutter_width_hz, mti_filter_t *filter) {
    if (filter == NULL || prf_hz <= 0.0) return -1;
    (void)clutter_width_hz; /* reserved for optimal MTI weight computation */
    filter->type = type;
    switch (type) {
    case MTI_SINGLE_DELAY:
        filter->filter_order = 2;
        filter->coefficients = (double complex *)malloc(2*sizeof(double complex));
        filter->coefficients[0] = 1.0 + 0.0*I;
        filter->coefficients[1] = -1.0 + 0.0*I;
        filter->clutter_attenuation_db = 20.0;
        filter->improvement_factor_db = 25.0;
        break;
    case MTI_DOUBLE_DELAY:
        filter->filter_order = 3;
        filter->coefficients = (double complex *)malloc(3*sizeof(double complex));
        filter->coefficients[0] = 1.0 + 0.0*I;
        filter->coefficients[1] = -2.0 + 0.0*I;
        filter->coefficients[2] = 1.0 + 0.0*I;
        filter->clutter_attenuation_db = 40.0;
        filter->improvement_factor_db = 45.0;
        break;
    case MTI_TRIPLE_DELAY:
        filter->filter_order = 4;
        filter->coefficients = (double complex *)malloc(4*sizeof(double complex));
        filter->coefficients[0] = 1.0 + 0.0*I;
        filter->coefficients[1] = -3.0 + 0.0*I;
        filter->coefficients[2] = 3.0 + 0.0*I;
        filter->coefficients[3] = -1.0 + 0.0*I;
        filter->clutter_attenuation_db = 60.0;
        filter->improvement_factor_db = 60.0;
        break;
    default: return -1;
    }
    filter->notch_width_hz = prf_hz / (double)filter->filter_order;
    return 0;
}

int mti_filter_apply(const mti_filter_t *filter,
                     double complex *slow_time_data, size_t length,
                     double complex *output) {
    if (filter == NULL || slow_time_data == NULL || output == NULL) return -1;
    if (filter->coefficients == NULL) return -1;
    size_t M = filter->filter_order;
    for (size_t n = 0; n < length; n++) {
        output[n] = 0.0 + 0.0*I;
        for (size_t k = 0; k < M; k++) {
            if (k <= n)
                output[n] += filter->coefficients[k] * slow_time_data[n - k];
        }
    }
    return 0;
}

void mti_filter_free(mti_filter_t *filter) {
    if (filter != NULL && filter->coefficients != NULL) {
        free(filter->coefficients);
        filter->coefficients = NULL;
    }
}

int adaptive_doppler_filter(const double complex *data_matrix,
                            uint32_t num_range_bins, uint32_t num_pulses,
                            uint32_t num_guard_cells,
                            double complex *filtered) {
    if (data_matrix == NULL || filtered == NULL) return -1;
    if (num_range_bins == 0 || num_pulses == 0) return -1;

    for (uint32_t r = 0; r < num_range_bins; r++) {
        for (uint32_t p = 0; p < num_pulses; p++) {
            double complex noise_est = 0.0 + 0.0*I;
            uint32_t count = 0;
            for (uint32_t rr = 0; rr < num_range_bins; rr++) {
                if (rr >= r - num_guard_cells && rr <= r + num_guard_cells) continue;
                if (rr < r || rr > r) {
                    noise_est += data_matrix[rr + p * num_range_bins];
                    count++;
                }
            }
            double complex signal_val = data_matrix[r + p * num_range_bins];
            if (count > 0) {
                double noise_mag = sqrt(creal(noise_est)*creal(noise_est)
                                      + cimag(noise_est)*cimag(noise_est)) / (double)count;
                double sig_mag = sqrt(creal(signal_val)*creal(signal_val)
                                    + cimag(signal_val)*cimag(signal_val));
                if (noise_mag > 1e-30)
                    filtered[r + p * num_range_bins] = signal_val * (sig_mag / noise_mag);
                else
                    filtered[r + p * num_range_bins] = signal_val;
            } else {
                filtered[r + p * num_range_bins] = signal_val;
            }
        }
    }
    return 0;
}
