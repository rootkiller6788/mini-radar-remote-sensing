/**
 * @file    radar_signal_model.c
 * @brief   Radar signal model: point targets, clutter, AWGN, range profiles
 */
#include "radar_signal_model.h"
#include "radar_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif
#define RADAR_C 299792458.0

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "signal_model: alloc(%zu) failed\n", sz); abort(); }
    return p;
}

int radar_target_init(radar_target_t *t, int id,
                      double range_m, double vr_ms,
                      double az_rad, double el_rad,
                      double rcs_m2, double lambda_m)
{
    if (!t || range_m <= 0.0 || lambda_m <= 0.0 || rcs_m2 < 0.0) return -1;
    t->target_id = id;
    t->range_m = range_m;
    t->velocity_ms = vr_ms;
    t->azimuth_rad = az_rad;
    t->elevation_rad = el_rad;
    t->rcs_m2 = rcs_m2;
    t->doppler_hz = 2.0 * vr_ms / lambda_m;
    t->delay_s = 2.0 * range_m / RADAR_C;
    t->amplitude = sqrt(rcs_m2) / (range_m * range_m);
    return 0;
}

size_t radar_range_to_bin(double range_m, double fs_hz)
{
    if (range_m <= 0.0 || fs_hz <= 0.0) return 0;
    return (size_t)(2.0 * range_m * fs_hz / RADAR_C);
}

double radar_bin_to_range(size_t bin, double fs_hz)
{
    if (fs_hz <= 0.0) return 0.0;
    return (double)bin * RADAR_C / (2.0 * fs_hz);
}

size_t radar_num_range_bins(double rmax_m, double fs_hz)
{
    if (rmax_m <= 0.0 || fs_hz <= 0.0) return 0;
    double nb = 2.0 * rmax_m * fs_hz / RADAR_C;
    return (size_t)(nb + 1.0);
}

int radar_generate_return(const double complex *tx, size_t tx_len,
                          const radar_target_t *tgt, double fs_hz,
                          double complex *rx, size_t rx_len,
                          double noise_pwr)
{
    if (!tx || !tgt || !rx || tx_len == 0 || rx_len == 0 || fs_hz <= 0.0)
        return -1;
    for (size_t i = 0; i < rx_len; i++) rx[i] = CMPLX(0.0, 0.0);
    double delay_samples = tgt->delay_s * fs_hz;
    size_t d = (size_t)delay_samples;
    if (d >= rx_len) return 0;
    double frac = delay_samples - (double)d;
    double omega = 2.0 * M_PI * tgt->doppler_hz / fs_hz;
    double A = tgt->amplitude;
    for (size_t n = 0; n < tx_len && (d + n) < rx_len; n++) {
        double complex s0 = tx[n];
        double complex s1 = (n + 1 < tx_len) ? tx[n+1] : tx[n];
        double complex s = s0 * (1.0 - frac) + s1 * frac;
        double phase = omega * (double)n;
        rx[d + n] += A * s * CMPLX(cos(phase), sin(phase));
    }
    if (noise_pwr > 0.0) {
        double std = sqrt(noise_pwr / 2.0);
        for (size_t i = 0; i < rx_len; i++) {
            double u1 = (double)rand() / RAND_MAX;
            double u2 = (double)rand() / RAND_MAX;
            if (u1 < 1e-12) u1 = 1e-12;
            double r = sqrt(-2.0 * log(u1));
            rx[i] += CMPLX(r * cos(2.0 * M_PI * u2) * std,
                           r * sin(2.0 * M_PI * u2) * std);
        }
    }
    return 0;
}

int radar_generate_multitarget(const double complex *tx, size_t tx_len,
                               const radar_target_t *tgs, size_t n_tgs,
                               double fs_hz,
                               double complex *rx, size_t rx_len,
                               double noise_pwr)
{
    if (!tx || !tgs || !rx || tx_len == 0 || rx_len == 0 || n_tgs == 0)
        return -1;
    for (size_t i = 0; i < rx_len; i++) rx[i] = CMPLX(0.0, 0.0);
    for (size_t t = 0; t < n_tgs; t++) {
        const radar_target_t *tgt = &tgs[t];
        double ds = tgt->delay_s * fs_hz;
        size_t d = (size_t)ds;
        if (d >= rx_len) continue;
        double frac = ds - (double)d;
        double omega = 2.0 * M_PI * tgt->doppler_hz / fs_hz;
        double A = tgt->amplitude;
        for (size_t n = 0; n < tx_len && (d + n) < rx_len; n++) {
            double complex s0 = tx[n];
            double complex s1 = (n + 1 < tx_len) ? tx[n+1] : tx[n];
            double complex s = s0 * (1.0 - frac) + s1 * frac;
            rx[d + n] += A * s * CMPLX(cos(omega * (double)n),
                                        sin(omega * (double)n));
        }
    }
    if (noise_pwr > 0.0) {
        double std = sqrt(noise_pwr / 2.0);
        for (size_t i = 0; i < rx_len; i++) {
            double u1 = (double)rand() / RAND_MAX;
            double u2 = (double)rand() / RAND_MAX;
            if (u1 < 1e-12) u1 = 1e-12;
            double r = sqrt(-2.0 * log(u1));
            rx[i] += CMPLX(r * cos(2.0 * M_PI * u2) * std,
                           r * sin(2.0 * M_PI * u2) * std);
        }
    }
    return 0;
}

int clutter_model_init(clutter_model_t *c, double s0_db,
                       clutter_distribution_t dist,
                       double grazing_rad, double corr_len_m)
{
    if (!c) return -1;
    c->dist = dist;
    c->sigma0_db = s0_db;
    c->sigma0_linear = pow(10.0, s0_db / 10.0);
    c->grazing_angle_rad = grazing_rad;
    c->correlation_length_m = corr_len_m;
    c->gamma_param = 1.0;
    c->mean_power = c->sigma0_linear;
    return 0;
}

int clutter_generate(const clutter_model_t *c,
                     size_t range_bins, double range_res_m,
                     double range_start_m,
                     double complex *samples, size_t n_samples)
{
    if (!c || !samples || n_samples == 0) return -1;
    double sigma = sqrt(c->mean_power / 2.0);
    for (size_t i = 0; i < n_samples; i++) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        if (u1 < 1e-12) u1 = 1e-12;
        double r = sqrt(-2.0 * log(u1));
        samples[i] = CMPLX(r * cos(2.0 * M_PI * u2) * sigma,
                           r * sin(2.0 * M_PI * u2) * sigma);
    }
    if (c->correlation_length_m > 0.0) {
        double alpha = exp(-range_res_m / c->correlation_length_m);
        for (size_t i = 1; i < n_samples; i++)
            samples[i] = alpha * samples[i-1] + (1.0 - alpha) * samples[i];
    }
    return 0;
}

double clutter_rcs_cell(double s0_lin, double range_m,
                        double bw_az_rad, double range_res_m)
{
    if (range_m <= 0.0 || bw_az_rad <= 0.0 || range_res_m <= 0.0) return 0.0;
    return s0_lin * range_m * bw_az_rad * range_res_m;
}

int radar_awgn(double complex *noise, size_t n_samples,
               double noise_pwr, unsigned int seed)
{
    if (!noise || n_samples == 0 || noise_pwr < 0.0) return -1;
    if (seed != 0) srand(seed);
    double std = sqrt(noise_pwr / 2.0);
    for (size_t i = 0; i < n_samples; i++) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        if (u1 < 1e-12) u1 = 1e-12;
        double r = sqrt(-2.0 * log(u1));
        noise[i] = CMPLX(r * cos(2.0 * M_PI * u2) * std,
                         r * sin(2.0 * M_PI * u2) * std);
    }
    return 0;
}

int radar_range_profile(const double complex *rx, size_t rx_len,
                        const double complex *tx, size_t tx_len,
                        double *profile)
{
    if (!rx || !tx || !profile || rx_len == 0 || tx_len == 0) return -1;
    if (tx_len > rx_len) return -1;
    size_t out_len = rx_len - tx_len + 1;
    for (size_t k = 0; k < out_len; k++) {
        double complex sum = CMPLX(0.0, 0.0);
        for (size_t n = 0; n < tx_len; n++)
            sum += rx[k + n] * conj(tx[n]);
        profile[k] = cabs(sum) * cabs(sum);
    }
    return 0;
}

/* L7: Ground clutter constant-gamma model */
double clutter_sigma0_gamma(double gamma_db, double grazing_angle_rad)
{
    double gamma = pow(10.0, gamma_db / 10.0);
    return gamma * sin(grazing_angle_rad);
}

/* L7: Sea clutter model (Nathanson empirical) */
double sea_clutter_sigma0_db(double wind_speed_kts, double freq_ghz,
                              double grazing_angle_rad, int sea_state)
{
    /* Nathanson model: sigma0 varies with sea state, wind, frequency */
    /* Simplified empirical formula */
    double base = -40.0 + (double)sea_state * 5.0;
    double freq_factor = 10.0 * log10(freq_ghz);
    double angle_factor = 20.0 * log10(sin(grazing_angle_rad));
    double wind_factor = 5.0 * log10(wind_speed_kts / 10.0);
    return base + freq_factor + angle_factor + wind_factor;
}

/* L7: Swerling target model selection helper */
const char* swerling_case_name(rcs_model_t model)
{
    switch (model) {
    case RCS_CONSTANT:    return "Swerling 0 (Non-fluctuating / Marcum)";
    case RCS_SWERLING_I:  return "Swerling I (Slow Rayleigh, scan-to-scan)";
    case RCS_SWERLING_II: return "Swerling II (Fast Rayleigh, pulse-to-pulse)";
    case RCS_SWERLING_III:return "Swerling III (Slow chi-square, 4 DoF)";
    case RCS_SWERLING_IV: return "Swerling IV (Fast chi-square, 4 DoF)";
    default:              return "Unknown";
    }
}

/* L3: Radar range equation in dB form */
double radar_range_equation_db(double pt_dbm, double gain_db,
                                double lambda_m, double rcs_dbsm,
                                double range_m, double loss_db)
{
    double lambda_db = 20.0 * log10(lambda_m);
    double range_db = 40.0 * log10(range_m);
    double const_db = -30.0 * log10(4.0 * M_PI);
    return pt_dbm + 2.0 * gain_db + lambda_db + rcs_dbsm +
           const_db - range_db - loss_db;
}

/* L7: Convert RCS from m^2 to dBsm */
double rcs_to_dbsm(double rcs_m2)
{
    if (rcs_m2 <= 0.0) return -INFINITY;
    return 10.0 * log10(rcs_m2);
}

/* L7: Convert RCS from dBsm to m^2 */
double rcs_from_dbsm(double rcs_dbsm)
{
    return pow(10.0, rcs_dbsm / 10.0);
}

/* L7: Compute radar sensitivity time control (STC) gain curve */
void radar_stc_curve(double max_range_m, size_t n_points,
                     double *range_vals, double *gain_db)
{
    if (!range_vals || !gain_db || n_points < 2) return;
    for (size_t i = 0; i < n_points; i++) {
        double r = max_range_m * (double)i / (double)(n_points - 1);
        range_vals[i] = r;
        if (r <= 0.0)
            gain_db[i] = 0.0;
        else
            gain_db[i] = 40.0 * log10(r);
    }
}

/* L7: Rain clutter volume reflectivity */
double rain_clutter_reflectivity_dbz(double rainfall_rate_mm_per_hr,
                                      double freq_ghz)
{
    /* Marshall-Palmer Z-R relationship: Z = 200 * R^1.6 */
    double Z = 200.0 * pow(rainfall_rate_mm_per_hr, 1.6);
    double Z_dbz = 10.0 * log10(Z);
    /* Frequency-dependent attenuation (ITU-R P.838) */
    double atten = pow(freq_ghz / 10.0, 1.5) * rainfall_rate_mm_per_hr * 0.01;
    return Z_dbz - atten;
}

/* L7: Chaff cloud RCS model */
double chaff_rcs_dbsm(size_t n_dipoles, double dipole_length_m,
                       double wavelength_m)
{
    if (n_dipoles == 0 || dipole_length_m <= 0.0 || wavelength_m <= 0.0)
        return -INFINITY;
    /* Resonant dipole RCS: sigma = 0.86 * lambda^2 */
    double single_rcs = 0.86 * wavelength_m * wavelength_m;
    double total_rcs = (double)n_dipoles * single_rcs;
    return 10.0 * log10(total_rcs);
}
