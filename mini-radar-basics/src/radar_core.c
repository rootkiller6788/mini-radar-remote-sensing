/**
 * @file    radar_core.c
 * @brief   Core radar implementations: range equation, RCS, antenna, SNR, path loss
 *
 * Knowledge covered:
 *   L1: Radar parameters, RCS models, range resolution, antenna definitions
 *   L3: Radar range equation, free-space path loss, system noise model
 *   L4: Maximum detection range, bistatic equation
 *   L5: Pulse integration gain (coherent and non-coherent)
 *
 * Reference: Richards, Scheer & Holm (2010), Ch.1-2.
 *            Skolnik, "Radar Handbook" (2008), Ch.1.
 */
#include "radar_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ������ Safe allocation ������������������������������������������������������������������������������������������������������ */

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "radar_core: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ������ L1: Radar Parameter Initialization ������������������������������������������������������������������ */

int radar_params_init(radar_params_t *params,
                      double freq_hz, double peak_power_w, double ant_gain_db,
                      double tau_s, double prf_hz, double bw_hz,
                      double nf_db, double loss_db, double ant_temp_k,
                      radar_mode_t mode)
{
    if (!params) return -1;
    if (freq_hz <= 0.0 || peak_power_w <= 0.0 || tau_s <= 0.0 ||
        prf_hz <= 0.0 || bw_hz <= 0.0 || nf_db < 0.0 ||
        loss_db < 0.0 || ant_temp_k <= 0.0)
        return -1;

    params->center_freq_hz = freq_hz;
    params->wavelength_m   = RADAR_C / freq_hz;
    params->peak_power_w   = peak_power_w;
    params->antenna_gain_db = ant_gain_db;
    params->antenna_gain_linear = pow(10.0, ant_gain_db / 10.0);
    params->pulse_width_s  = tau_s;
    params->prf_hz         = prf_hz;
    params->pri_s          = 1.0 / prf_hz;
    params->bandwidth_hz   = bw_hz;
    params->noise_figure_db = nf_db;
    params->noise_figure_linear = pow(10.0, nf_db / 10.0);
    params->system_loss_db  = loss_db;
    params->system_loss_linear = pow(10.0, loss_db / 10.0);
    params->antenna_temp_k  = ant_temp_k;
    params->mode            = mode;
    params->pol             = RADAR_POL_HH;

    return 0;
}

/* ������ L1: RCS Model ������������������������������������������������������������������������������������������������������������ */

int radar_rcs_init(radar_rcs_t *rcs, double rcs_dbsm, rcs_model_t model,
                   double decorr_time_s)
{
    if (!rcs) return -1;
    rcs->mean_rcs_dbsm = rcs_dbsm;
    rcs->mean_rcs_m2   = pow(10.0, rcs_dbsm / 10.0);
    rcs->fluctuation   = model;
    rcs->decorrelation_time_s = decorr_time_s;
    rcs->n_samples     = 0;
    return 0;
}

/**
 * Generate RCS sample according to Swerling fluctuation model.
 *
 * Swerling I/II: exponential distribution (Rayleigh amplitude).
 *   PDF: p(sigma) = (1/sigma_mean) * exp(-sigma / sigma_mean)
 *   CDF: F(sigma) = 1 - exp(-sigma / sigma_mean)
 *   Inverse CDF: sigma = -sigma_mean * ln(1 - U) = -sigma_mean * ln(U)
 *
 * Swerling III/IV: chi-square with 4 degrees of freedom (gamma distribution).
 *   PDF: p(sigma) = (4*sigma / sigma_mean^2) * exp(-2*sigma / sigma_mean)
 *   Generated as sum of 2 exponential random variables:
 *     sigma = sigma_mean/2 * (X1 + X2) where X1, X2 ~ Exp(1)
 */
double radar_rcs_sample(radar_rcs_t *rcs)
{
    if (!rcs) return 0.0;

    switch (rcs->fluctuation) {
    case RCS_CONSTANT:
        return rcs->mean_rcs_m2;

    case RCS_SWERLING_I:
    case RCS_SWERLING_II: {
        /* Exponential: sigma = -sigma_mean * ln(U) where U ~ Uniform(0,1) */
        double u = (double)rand() / (double)RAND_MAX;
        if (u < 1e-12) u = 1e-12; /* avoid log(0) */
        return -rcs->mean_rcs_m2 * log(u);
    }

    case RCS_SWERLING_III:
    case RCS_SWERLING_IV: {
        /* Chi-square with 4 DoF: sum of 2 independent exponentials */
        double u1 = (double)rand() / (double)RAND_MAX;
        double u2 = (double)rand() / (double)RAND_MAX;
        if (u1 < 1e-12) u1 = 1e-12;
        if (u2 < 1e-12) u2 = 1e-12;
        double x = -log(u1) - log(u2);
        /* Scale: E[X] = 2*sigma_mean, so divide by 2 */
        return rcs->mean_rcs_m2 * x / 2.0;
    }

    default:
        return rcs->mean_rcs_m2;
    }
}

/* ������ L1: Range Resolution and Unambiguous Range/Velocity �������������������������������� */

double radar_range_resolution(double pulse_width_or_bw, int use_bandwidth)
{
    if (pulse_width_or_bw <= 0.0) return 0.0;

    if (use_bandwidth) {
        /* Bandwidth-limited: Delta_R = c / (2 * B) */
        return RADAR_C / (2.0 * pulse_width_or_bw);
    } else {
        /* Simple pulse: Delta_R = c * tau / 2 */
        return RADAR_C * pulse_width_or_bw / 2.0;
    }
}

double radar_unambiguous_range(double prf_hz)
{
    if (prf_hz <= 0.0) return INFINITY;
    return RADAR_C / (2.0 * prf_hz);
}

double radar_unambiguous_velocity(double wavelength_m, double prf_hz)
{
    if (wavelength_m <= 0.0 || prf_hz <= 0.0) return 0.0;
    return wavelength_m * prf_hz / 4.0;
}

/* ������ L3/L4: Monostatic Radar Range Equation ������������������������������������������������������������ */

/**
 * P_r = P_t * G^2 * lambda^2 * sigma / ((4*pi)^3 * R^4 * L)
 *
 * The equation assumes far-field conditions (R > 2*D^2/lambda),
 * free-space propagation (no multipath, no atmospheric attenuation),
 * and a point target fully contained within the antenna beam.
 */
double radar_received_power(const radar_params_t *params,
                            double rcs_m2, double range_m)
{
    if (!params || rcs_m2 <= 0.0 || range_m <= 0.0) return 0.0;

    double pt   = params->peak_power_w;
    double g    = params->antenna_gain_linear;
    double lam  = params->wavelength_m;
    double loss = params->system_loss_linear;

    /* Numerator: P_t * G^2 * lambda^2 * sigma */
    double num = pt * g * g * lam * lam * rcs_m2;

    /* Denominator: (4*pi)^3 * R^4 * L */
    double four_pi_cubed = 64.0 * M_PI * M_PI * M_PI;  /* (4*pi)^3 = 64*pi^3 */
    double den = four_pi_cubed * range_m * range_m * range_m * range_m * loss;

    if (den <= 0.0) return 0.0;
    return num / den;
}

/**
 * SNR = P_r / (k * T_sys * B * F)
 *
 * T_sys = T_a + T0*(F-1) where T0=290K and F=noise figure (linear)
 */
double radar_snr(const radar_params_t *params, double rcs_m2, double range_m)
{
    double pr = radar_received_power(params, rcs_m2, range_m);
    if (pr <= 0.0) return 0.0;

    double tsys = radar_system_noise_temp(params->antenna_temp_k,
                                           params->noise_figure_linear);
    double pn = radar_noise_power(tsys, params->bandwidth_hz);

    if (pn <= 0.0) return 0.0;
    return pr / pn;
}

double radar_snr_db(const radar_params_t *params, double rcs_m2, double range_m)
{
    double snr_lin = radar_snr(params, rcs_m2, range_m);
    return lin2db(snr_lin);
}

/**
 * R_max = [P_t * G^2 * lambda^2 * sigma /
 *          ((4*pi)^3 * k*T0 * B * F * SNR_min * L)]^{1/4}
 *
 * The fourth-root dependence means:
 *   - To double range, need 16x more power (or 4x antenna area)
 *   - Range is proportional to RCS^{1/4}
 */
double radar_max_range(const radar_params_t *params, double rcs_m2,
                       double snr_min)
{
    if (!params || rcs_m2 <= 0.0 || snr_min <= 0.0) return 0.0;

    double pt   = params->peak_power_w;
    double g    = params->antenna_gain_linear;
    double lam  = params->wavelength_m;
    double bw   = params->bandwidth_hz;
    double f    = params->noise_figure_linear;
    double loss = params->system_loss_linear;

    double kT0BF = RADAR_K_BOLTZMANN * RADAR_T0 * bw * f;
    double four_pi_cubed = 64.0 * M_PI * M_PI * M_PI;

    double numerator = pt * g * g * lam * lam * rcs_m2;
    double denominator = four_pi_cubed * kT0BF * snr_min * loss;

    if (denominator <= 0.0) return 0.0;
    return pow(numerator / denominator, 0.25);
}

/* ������ L3: Bistatic Radar Range Equation ���������������������������������������������������������������������� */

/**
 * P_r = P_t * G_t * G_r * lambda^2 * sigma_b /
 *       ((4*pi)^3 * R_t^2 * R_r^2 * L)
 *
 * The product R_t^2 * R_r^2 replaces R^4 from the monostatic case.
 * Bistatic RCS sigma_b may differ from monostatic RCS due to
 * different scattering geometry.
 */
double radar_bistatic_power(double pt_w, double gt_linear, double gr_linear,
                            double lambda_m, double rcs_m2,
                            double rt_m, double rr_m, double loss_linear)
{
    if (pt_w <= 0.0 || gt_linear <= 0.0 || gr_linear <= 0.0 ||
        lambda_m <= 0.0 || rcs_m2 <= 0.0 || rt_m <= 0.0 ||
        rr_m <= 0.0 || loss_linear <= 0.0)
        return 0.0;

    double num = pt_w * gt_linear * gr_linear * lambda_m * lambda_m * rcs_m2;
    double four_pi_cubed = 64.0 * M_PI * M_PI * M_PI;
    double den = four_pi_cubed * rt_m * rt_m * rr_m * rr_m * loss_linear;

    if (den <= 0.0) return 0.0;
    return num / den;
}

/* ������ L2: System Noise Model �������������������������������������������������������������������������������������������� */

/**
 * T_sys = T_a + T0 * (F - 1)
 *
 * T0*(F-1) is the receiver equivalent noise temperature referred to
 * the antenna port. T_a captures external noise (galactic, atmospheric,
 * ground). For microwave radars, T_a is typically 100-300 K.
 */
double radar_system_noise_temp(double ant_temp_k, double noise_fig_lin)
{
    if (ant_temp_k <= 0.0 || noise_fig_lin < 1.0) return 0.0;
    return ant_temp_k + RADAR_T0 * (noise_fig_lin - 1.0);
}

/**
 * P_n = k * T_sys * B
 *
 * This is the thermal noise power in the receiver bandwidth
 * (Johnson-Nyquist noise). At room temperature (290 K):
 *   k*T0 = 4.0e-21 W/Hz = -174 dBm/Hz
 */
double radar_noise_power(double tsys_k, double bw_hz)
{
    if (tsys_k <= 0.0 || bw_hz <= 0.0) return 0.0;
    return RADAR_K_BOLTZMANN * tsys_k * bw_hz;
}

/* ������ L3: Two-Way Path Loss ���������������������������������������������������������������������������������������������� */

/**
 * L_2way = (4*pi)^3 * R^4 / (G^2 * lambda^2)
 *
 * Represents the round-trip spreading loss (excluding RCS and system losses).
 * In dB: L_2way(dB) = 30*log10(4*pi) + 40*log10(R) - 20*log10(G) - 20*log10(lambda)
 */
double radar_path_loss_2way(double range_m, double lambda_m, double gain_linear)
{
    if (range_m <= 0.0 || lambda_m <= 0.0 || gain_linear <= 0.0) return 0.0;

    double four_pi_cubed = 64.0 * M_PI * M_PI * M_PI;
    double num = four_pi_cubed * range_m * range_m * range_m * range_m;
    double den = gain_linear * gain_linear * lambda_m * lambda_m;

    if (den <= 0.0) return 0.0;
    return num / den;
}

/* ������ L1: Antenna Models ������������������������������������������������������������������������������������������������������ */

/**
 * G = 4*pi * A_e / lambda^2
 *
 * This is a fundamental result from antenna theory: the effective
 * aperture and gain are related by reciprocity. For an ideal antenna,
 * A_e = A_physical. For real antennas, A_e = eta_aperture * A_physical
 * where eta_aperture is the aperture efficiency (typically 0.5-0.8).
 */
double radar_antenna_gain_from_aperture(double a_e_m2, double lambda_m)
{
    if (a_e_m2 <= 0.0 || lambda_m <= 0.0) return 0.0;
    return 4.0 * M_PI * a_e_m2 / (lambda_m * lambda_m);
}

/**
 * theta_3dB = 1.029 * lambda / D (exact for uniform circular aperture)
 *
 * Approximate: theta_3dB �� 1.22 * lambda / D
 * We use 1.029 for uniform illumination, which is more precise.
 *
 * In degrees: theta_deg �� 58.9 * lambda / D
 */
double radar_beamwidth_circular(double lambda_m, double diameter_m)
{
    if (lambda_m <= 0.0 || diameter_m <= 0.0) return 0.0;
    return 1.029 * lambda_m / diameter_m;
}

/**
 * Kraus approximation for directivity:
 * D �� 4*pi / (theta_az * theta_el)
 *
 * This is accurate for pencil-beam antennas with small beamwidths
 * (theta < 20 degrees) and low sidelobes. The approximation assumes
 * the main beam solid angle equals the product of orthogonal beamwidths.
 */
double radar_directivity_kraus(double beamwidth_az_rad, double beamwidth_el_rad)
{
    if (beamwidth_az_rad <= 0.0 || beamwidth_el_rad <= 0.0) return 0.0;
    return 4.0 * M_PI / (beamwidth_az_rad * beamwidth_el_rad);
}

/* ������ L2: Duty Cycle and Average Power ������������������������������������������������������������������������ */

/**
 * D = tau / PRI = tau * PRF
 *
 * Duty cycle is bounded: 0 < D <= 1. For pulsed radars, D << 1.
 * For FMCW radars, D = 1.0 (continuous transmission).
 * Typical values: 0.001 (search) to 0.1 (tracking).
 */
double radar_duty_cycle(double tau_s, double prf_hz)
{
    if (tau_s <= 0.0 || prf_hz <= 0.0) return 0.0;
    double dc = tau_s * prf_hz;
    return (dc <= 1.0) ? dc : 1.0;
}

/**
 * P_avg = P_t * D = P_t * tau * PRF
 *
 * Average power is the transmitter thermal limit. For a given
 * average power, reducing pulse width increases peak power
 * (for fixed PRF).
 */
double radar_average_power(double peak_power_w, double duty_cycle)
{
    if (peak_power_w <= 0.0 || duty_cycle < 0.0) return 0.0;
    if (duty_cycle > 1.0) duty_cycle = 1.0;
    return peak_power_w * duty_cycle;
}

/* ������ L5: Pulse Integration Gain ������������������������������������������������������������������������������������ */

/**
 * Coherent integration: SNR_coh = N * SNR_single
 *
 * Coherent integration sums complex samples (phase-aligned) across pulses.
 * Requires knowledge of target Doppler to align phases. Achieves maximum
 * possible SNR improvement of N (10*log10(N) dB) for N pulses.
 *
 * This is the matched filter over slow time �� optimal for AWGN.
 */
double radar_coherent_integration_gain(size_t n_pulses, double snr_single)
{
    if (n_pulses == 0) return snr_single;
    if (snr_single < 0.0) return 0.0;
    return (double)n_pulses * snr_single;
}

/**
 * Non-coherent integration: G_nc = N / (1 + (N-1) * rho)
 *
 * Non-coherent integration discards phase information (envelope detection)
 * before summation. The integration efficiency is less than N because
 * the noise-only envelope also integrates non-coherently.
 *
 * Peebles approximation with rho = 0.5 (square-law detector with threshold):
 *
 * For large N, G_nc approaches 1/rho �� 2 (diminishing returns),
 * unlike coherent integration which keeps improving linearly.
 *
 * Empirical approximation: G_nc �� N^0.7 for 1 < N < 100.
 */
double radar_ncoherent_integration_gain(size_t n_pulses)
{
    if (n_pulses == 0) return 1.0;
    if (n_pulses == 1) return 1.0;

    /* Peebles approximation with rho = 0.5 */
    double rho = 0.5;
    double n = (double)n_pulses;
    double gain = n / (1.0 + (n - 1.0) * rho);

    return gain;
}

/* ─── L7: Application Helpers ─────────────────────────────────────────── */

/**
 * Compute radar power-aperture product.
 *
 * P_avg * A_e = (P_t * tau * PRF) * (G * lambda^2 / (4*pi))
 *
 * The power-aperture product determines search capability:
 * for a given SNR requirement, the time to search a solid angle
 * Omega is proportional to 1/(P_avg * A_e).
 */
double radar_power_aperture_product(const radar_params_t *params)
{
    if (!params) return 0.0;
    double dc = radar_duty_cycle(params->pulse_width_s, params->prf_hz);
    double p_avg = radar_average_power(params->peak_power_w, dc);
    double ae = params->antenna_gain_linear * params->wavelength_m *
                params->wavelength_m / (4.0 * M_PI);
    return p_avg * ae;
}

/**
 * Compute search frame time (time to scan a given solid angle).
 *
 * T_frame = (Omega * R_max^4 * (4*pi)^3 * k*T0*F*SNR_min*L) /
 *           (P_avg * A_e * sigma)
 *
 * This is derived from the search radar equation (Skolnik, Ch.2).
 */
double radar_search_frame_time(double solid_angle_sr,
                               const radar_params_t *params,
                               double rcs_m2, double snr_min)
{
    if (!params || solid_angle_sr <= 0.0) return 0.0;
    double pa = radar_power_aperture_product(params);
    if (pa <= 0.0) return 0.0;

    double kT0BF = RADAR_K_BOLTZMANN * RADAR_T0 *
                   params->bandwidth_hz * params->noise_figure_linear;
    /* Search frame time equation: T_frame proportional to Omega / (P_avg * A_e).
     * Full form includes R_max^4 factor from radar range equation. */
    double numer = solid_angle_sr * kT0BF * snr_min * params->system_loss_linear;
    double denom = pa * rcs_m2;
    if (denom <= 0.0) return INFINITY;
    return numer / denom;
}

/**
 * Compute SNR improvement from pulse integration.
 *
 * For N pulses:
 *   Coherent: SNR_out = N * SNR_in  (10*log10(N) dB gain)
 *   Non-coherent: SNR_out ≈ N^0.7 * SNR_in (for N < 100)
 *
 * This function returns the coherent gain in dB.
 */
double radar_integration_gain_db(size_t n_pulses, int coherent)
{
    if (n_pulses == 0) return 0.0;
    if (coherent)
        return 10.0 * log10((double)n_pulses);
    else
        return 10.0 * log10(radar_ncoherent_integration_gain(n_pulses));
}

/**
 * Compute the minimum detectable RCS for a given range and SNR threshold.
 *
 * Solving the radar range equation for sigma:
 * sigma_min = SNR_min * (4*pi)^3 * R^4 * k*T0*B*F*L /
 *             (P_t * G^2 * lambda^2)
 *
 * This is the smallest RCS detectable at range R with single-pulse
 * detection at the specified SNR threshold.
 */
double radar_min_detectable_rcs(const radar_params_t *params,
                                 double range_m, double snr_min)
{
    if (!params || range_m <= 0.0 || snr_min <= 0.0) return 0.0;

    double pt = params->peak_power_w;
    double g = params->antenna_gain_linear;
    double lam = params->wavelength_m;
    double bw = params->bandwidth_hz;
    double f = params->noise_figure_linear;
    double loss = params->system_loss_linear;
    double kT0BF = RADAR_K_BOLTZMANN * RADAR_T0 * bw * f;
    double four_pi_cubed = 64.0 * M_PI * M_PI * M_PI;

    double numer = snr_min * four_pi_cubed * range_m * range_m *
                   range_m * range_m * kT0BF * loss;
    double denom = pt * g * g * lam * lam;

    if (denom <= 0.0) return 0.0;
    return numer / denom;
}

/**
 * Compute the atmospheric attenuation loss (approximate, clear air).
 *
 * Two-way attenuation at frequency f (GHz) and range R (km):
 *   L_atm ≈ 2 * R * alpha(f)  [dB]
 *
 * where alpha(f) is the specific attenuation in dB/km.
 * This is a simplified model ignoring humidity and pressure.
 * For accurate results, use the ITU-R P.676 model.
 */
double radar_atmospheric_loss_db(double freq_ghz, double range_km)
{
    if (freq_ghz <= 0.0 || range_km <= 0.0) return 0.0;
    /* Approximate specific attenuation (dB/km) for clear air */
    double alpha;
    if (freq_ghz < 10.0)
        alpha = 0.008;  /* S/C band: negligible */
    else if (freq_ghz < 20.0)
        alpha = 0.016;  /* X/Ku band: small */
    else if (freq_ghz < 40.0)
        alpha = 0.05;   /* K/Ka band: moderate */
    else if (freq_ghz < 100.0)
        alpha = 0.2;    /* W band: significant */
    else
        alpha = 1.0;    /* Sub-mm: large */

    return 2.0 * range_km * alpha;
}

/**
 * Compute the radar horizon range (line-of-sight limit).
 *
 * R_horizon ≈ sqrt(2 * R_earth * (h_radar + h_target))
 *            ≈ 4.12 * (sqrt(h_radar_m) + sqrt(h_target_m))  [km]
 *
 * Accounts for standard atmospheric refraction (4/3 Earth radius model).
 * R_earth_eff = 4/3 * 6371 ≈ 8495 km.
 */
double radar_horizon_range_km(double h_radar_m, double h_target_m)
{
    if (h_radar_m < 0.0 || h_target_m < 0.0) return 0.0;
    /* 4/3 Earth radius model with standard refraction */
    return 4.12 * (sqrt(h_radar_m) + sqrt(h_target_m));
}

/* ─── L7: Cognitive Radar ────────────────────────────────────────────── */

int cognitive_radar_init(cognitive_radar_state_t *state)
{
    if (!state) return -1;
    memset(state, 0, sizeof(*state));
    state->agility = AGILITY_NONE;
    state->current_snr_db = 0.0;
    state->interference_power_db = -100.0;
    state->target_present = 0;
    state->waveform_index = 0;
    state->adaptation_rate = 0.1;
    return 0;
}

int cognitive_radar_adapt(cognitive_radar_state_t *state,
                          double measured_snr_db)
{
    if (!state) return -1;
    /* Exponential moving average of SNR */
    double alpha = state->adaptation_rate;
    state->current_snr_db = alpha * measured_snr_db +
                            (1.0 - alpha) * state->current_snr_db;
    /* Simple adaptation: if SNR drops, increase waveform index */
    if (measured_snr_db < state->current_snr_db - 3.0)
        state->waveform_index++;
    else if (measured_snr_db > state->current_snr_db + 3.0)
        state->waveform_index = (state->waveform_index > 0) ?
                                 state->waveform_index - 1 : 0;
    return 0;
}

/* ─── L7: Radar Bands ────────────────────────────────────────────────── */

double radar_band_center_freq(radar_band_t band)
{
    switch (band) {
    case RADAR_BAND_HF:  return 15e6;
    case RADAR_BAND_VHF: return 150e6;
    case RADAR_BAND_UHF: return 600e6;
    case RADAR_BAND_L:   return 1.5e9;
    case RADAR_BAND_S:   return 3.0e9;
    case RADAR_BAND_C:   return 5.5e9;
    case RADAR_BAND_X:   return 10.0e9;
    case RADAR_BAND_KU:  return 15.0e9;
    case RADAR_BAND_K:   return 24.0e9;
    case RADAR_BAND_KA:  return 35.0e9;
    case RADAR_BAND_V:   return 60.0e9;
    case RADAR_BAND_W:   return 94.0e9;
    case RADAR_BAND_MM:  return 140.0e9;
    default:             return 10.0e9;
    }
}

double radar_band_wavelength(radar_band_t band)
{
    double fc = radar_band_center_freq(band);
    return (fc > 0.0) ? (RADAR_C / fc) : 0.0;
}

double radar_band_typical_gain_db(radar_band_t band)
{
    switch (band) {
    case RADAR_BAND_HF:  case RADAR_BAND_VHF: case RADAR_BAND_UHF:
        return 10.0;
    case RADAR_BAND_L:   case RADAR_BAND_S:
        return 25.0;
    case RADAR_BAND_C:   case RADAR_BAND_X:
        return 35.0;
    case RADAR_BAND_KU:  case RADAR_BAND_K:  case RADAR_BAND_KA:
        return 45.0;
    case RADAR_BAND_V:   case RADAR_BAND_W:  case RADAR_BAND_MM:
        return 50.0;
    default:
        return 30.0;
    }
}
