/**
 * @file    ir_detector.c
 * @brief   Infrared Detector Models - Implementation
 *
 * Implements detector characterization: responsivity, noise, NEP, D*, NETD, MRTD.
 * Covers photon detectors (InSb, MCT, QWIP) and thermal detectors (microbolometers).
 *
 * References:
 *   Rogalski (2011) "Infrared Detectors", 2nd Ed.
 *   Vincent (2015) "Fundamentals of IR Detector Operation and Testing"
 */

#include "ir_detector.h"
#include "ir_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* =========================================================================
 * Material Database
 * ========================================================================= */

static const ir_detector_material_t g_materials[IR_DET_COUNT] = {
    /* Thermal detectors */
    [IR_DET_THERMAL_BOLOMETER] = {
        "VOx / a-Si", 0.0, 30.0, 300.0, -0.025
    },
    [IR_DET_THERMAL_THERMOPILE] = {
        "Bi-Sb / Poly-Si", 0.0, 30.0, 300.0, -0.001
    },
    [IR_DET_THERMAL_PYROELECTRIC] = {
        "LiTaO3 / BST / PZT", 0.0, 30.0, 300.0, -0.005
    },
    /* Photon detectors */
    [IR_DET_PHOTON_IN_SB] = {
        "InSb", 0.23, 5.4, 77.0, 0.0
    },
    [IR_DET_PHOTON_MCT] = {
        "Hg1-xCdxTe", 0.10, 12.4, 77.0, 0.0
    },
    [IR_DET_PHOTON_QWIP] = {
        "GaAs/AlGaAs", 0.15, 8.3, 70.0, 0.0
    },
    [IR_DET_PHOTON_T2SL] = {
        "InAs/GaSb", 0.12, 10.3, 77.0, 0.0
    },
    [IR_DET_PHOTON_PBS] = {
        "PbS", 0.42, 3.0, 300.0, 0.0
    }
};

const ir_detector_material_t* ir_detector_material_get(ir_detector_type_t type) {
    if (type < 0 || type >= IR_DET_COUNT) return NULL;
    return &g_materials[type];
}

/* =========================================================================
 * Detector Initialization
 * ========================================================================= */

int ir_detector_params_init(ir_detector_params_t *p, ir_detector_type_t type) {
    if (!p || type < 0 || type >= IR_DET_COUNT) return -1;

    memset(p, 0, sizeof(*p));
    p->type = type;
    p->material = g_materials[type];

    /* Set sensible defaults based on detector type */
    switch (type) {
    case IR_DET_THERMAL_BOLOMETER:
        p->cooling = IR_DET_MODE_UNCOOLED;
        p->pixel_pitch_um = 17.0;
        p->fill_factor = 0.80;
        p->active_area_cm2 = 2.31e-6;  /* 17um pitch, 80% FF */
        p->quantum_efficiency = 0.80;  /* absorption efficiency */
        p->spectral_response_peak_um = 10.0;
        p->spectral_bandwidth_um = 6.0;
        p->resistance_ohm = 1.0e6;
        p->bias_voltage_V = 2.0;
        p->capacitance_F = 1.0e-12;
        p->thermal_conductance_W_per_K = 1.0e-7;
        p->thermal_mass_J_per_K = 1.0e-9;
        p->thermal_time_constant_ms = 10.0;
        p->integration_time_us = 33;
        p->frame_rate_Hz = 30.0;
        break;
    case IR_DET_PHOTON_MCT:
        p->cooling = IR_DET_MODE_STIRLING;
        p->pixel_pitch_um = 15.0;
        p->fill_factor = 0.90;
        p->active_area_cm2 = 2.025e-6;
        p->quantum_efficiency = 0.70;
        p->spectral_response_peak_um = 10.0;
        p->spectral_bandwidth_um = 4.0;
        p->resistance_ohm = 100.0;
        p->bias_voltage_V = 0.5;
        p->capacitance_F = 0.5e-12;
        p->thermal_conductance_W_per_K = 0.0;
        p->thermal_mass_J_per_K = 0.0;
        p->thermal_time_constant_ms = 0.0;
        p->integration_time_us = 16;
        p->frame_rate_Hz = 60.0;
        break;
    case IR_DET_PHOTON_IN_SB:
        p->cooling = IR_DET_MODE_STIRLING;
        p->pixel_pitch_um = 15.0;
        p->fill_factor = 0.90;
        p->active_area_cm2 = 2.025e-6;
        p->quantum_efficiency = 0.80;
        p->spectral_response_peak_um = 4.5;
        p->spectral_bandwidth_um = 2.0;
        p->resistance_ohm = 500.0;
        p->bias_voltage_V = 0.3;
        p->capacitance_F = 1.0e-12;
        p->integration_time_us = 10;
        p->frame_rate_Hz = 100.0;
        break;
    default:
        /* Generic defaults */
        p->pixel_pitch_um = 20.0;
        p->fill_factor = 0.80;
        p->active_area_cm2 = 3.2e-6;
        p->quantum_efficiency = 0.50;
        p->spectral_response_peak_um = 5.0;
        p->spectral_bandwidth_um = 2.0;
        p->integration_time_us = 20;
        p->frame_rate_Hz = 50.0;
        break;
    }
    return 0;
}

/* =========================================================================
 * Responsivity Calculations
 *
 * Photon detector (voltage): R_v = eta * q * lambda * R_load / (h * c)
 * Photon detector (current): R_i = eta * q * lambda / (h * c)
 * Thermal detector:          R_v = eta * alpha * I_bias * R / (G_th * sqrt(1+(w*tau)^2))
 * ========================================================================= */

double ir_detector_responsivity_voltage(const ir_detector_params_t *params,
                                         double wavelength_um) {
    if (!params || wavelength_um <= 0.0) return -1.0;

    double lambda_m = wavelength_um * 1e-6;
    double hc = IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT;
    double q = IR_ELEMENTARY_CHARGE;

    double R_i = params->quantum_efficiency * q * lambda_m / hc;  /* A/W */
    return R_i * params->resistance_ohm;  /* V/W */
}

double ir_detector_responsivity_current(double quantum_efficiency,
                                         double wavelength_um) {
    if (quantum_efficiency < 0.0 || quantum_efficiency > 1.0) return -1.0;
    if (wavelength_um <= 0.0) return -1.0;

    double lambda_m = wavelength_um * 1e-6;
    double hc = IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT;
    return quantum_efficiency * IR_ELEMENTARY_CHARGE * lambda_m / hc;
}

double ir_detector_responsivity_thermal(const ir_detector_params_t *params,
                                         double freq_Hz) {
    if (!params || freq_Hz < 0.0) return -1.0;
    if (params->thermal_conductance_W_per_K <= 0.0) return -1.0;

    double G_th = params->thermal_conductance_W_per_K;
    double alpha = params->material.thermal_coeff_K_inv;  /* TCR */
    double I_bias = params->bias_voltage_V / params->resistance_ohm;
    double R = params->resistance_ohm;
    double tau = params->thermal_time_constant_ms * 1e-3;  /* convert ms to s */
    double omega = 2.0 * M_PI * freq_Hz;

    double denom = G_th * sqrt(1.0 + omega * omega * tau * tau);
    if (denom <= 0.0) return -1.0;

    return fabs(params->quantum_efficiency * alpha * I_bias * R / denom);
}

/* =========================================================================
 * Noise Mechanisms
 *
 * Johnson-Nyquist: V_n_rms = sqrt(4 * kB * T * R * delta_f)
 * Shot noise:      I_n_rms = sqrt(2 * q * I_dc * delta_f)
 * 1/f (flicker):   S_I(f) = K_f * I_dc^alpha / f^beta
 * ========================================================================= */

double ir_noise_johnson_voltage(double temperature_K,
                                 double resistance_ohm,
                                 double bandwidth_Hz) {
    if (temperature_K <= 0.0 || resistance_ohm <= 0.0 || bandwidth_Hz <= 0.0)
        return -1.0;
    return sqrt(4.0 * IR_BOLTZMANN_CONSTANT * temperature_K
                * resistance_ohm * bandwidth_Hz);
}

double ir_noise_shot_current(double dc_current_A, double bandwidth_Hz) {
    if (dc_current_A < 0.0 || bandwidth_Hz <= 0.0) return -1.0;
    return sqrt(2.0 * IR_ELEMENTARY_CHARGE * dc_current_A * bandwidth_Hz);
}

double ir_noise_flicker_psd(double dc_current_A, double frequency_Hz,
                             double k_factor, double alpha, double beta) {
    if (dc_current_A < 0.0 || frequency_Hz <= 0.0) return -1.0;
    if (k_factor <= 0.0 || beta <= 0.0) return -1.0;
    return k_factor * pow(dc_current_A, alpha) / pow(frequency_Hz, beta);
}

double ir_noise_total_current(const ir_detector_params_t *params,
                               double bandwidth_Hz, double frequency_Hz) {
    if (!params || bandwidth_Hz <= 0.0) return -1.0;

    /* Johnson noise current = V_johnson / R */
    double V_j = ir_noise_johnson_voltage(
        params->material.operating_temperature_K,
        params->resistance_ohm, bandwidth_Hz);
    double I_j = (params->resistance_ohm > 0.0)
                 ? V_j / params->resistance_ohm : 0.0;

    /* Shot noise from bias current */
    double I_bias = params->bias_voltage_V / params->resistance_ohm;
    double I_shot = ir_noise_shot_current(I_bias, bandwidth_Hz);

    /* 1/f noise at operating frequency */
    double I_flicker = ir_noise_flicker_psd(I_bias, frequency_Hz,
                                             1e-11, 1.0, 1.0);
    /* Convert PSD to RMS: I_rms = sqrt(S_I * delta_f) */
    I_flicker = sqrt(I_flicker * bandwidth_Hz);

    return sqrt(I_j * I_j + I_shot * I_shot + I_flicker * I_flicker);
}

/* =========================================================================
 * Figures of Merit: NEP and D*
 * ========================================================================= */

double ir_nep(const ir_detector_params_t *params, double bandwidth_Hz,
               double frequency_Hz, double wavelength_um) {
    if (!params || bandwidth_Hz <= 0.0) return -1.0;

    double I_noise = ir_noise_total_current(params, bandwidth_Hz, frequency_Hz);
    double R_i = ir_detector_responsivity_current(
        params->quantum_efficiency, wavelength_um);

    if (R_i <= 0.0) return -1.0;
    return I_noise / R_i;  /* NEP = noise_current / responsivity [W] */
}

double ir_specific_detectivity(const ir_detector_params_t *params,
                                double nep_W, double bandwidth_Hz) {
    if (!params || nep_W <= 0.0 || bandwidth_Hz <= 0.0) return -1.0;
    double A_d = params->active_area_cm2;
    if (A_d <= 0.0) return -1.0;
    return sqrt(A_d * bandwidth_Hz) / nep_W;
}

double ir_d_star_blip(double quantum_eff, double wavelength_um,
                       double photon_flux_bg) {
    if (quantum_eff <= 0.0 || quantum_eff > 1.0) return -1.0;
    if (wavelength_um <= 0.0 || photon_flux_bg <= 0.0) return -1.0;

    double lambda_m = wavelength_um * 1e-6;
    double hc = IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT;
    return (lambda_m / hc) * sqrt(quantum_eff / (2.0 * photon_flux_bg));
}

double ir_d_star_compare(double d_star1, double d_star2) {
    if (d_star1 <= 0.0 || d_star2 <= 0.0) return -1.0;
    return d_star1 / d_star2;
}

/* =========================================================================
 * NETD - Noise Equivalent Temperature Difference
 *
 * NETD = 4 * F_no^2 * V_n / (tau_opt * A_d * dL/dT * R_v)
 *
 * The NETD represents the minimum temperature difference detectable
 * above the sensor noise floor. Typical values:
 *   Uncooled microbolometer: 30-100 mK
 *   Cooled MWIR photon:      10-25 mK
 *   Cooled LWIR photon:      15-30 mK
 * ========================================================================= */

double ir_radiance_derivative_wrt_temperature(double T_K,
                                               double lambda1_um,
                                               double lambda2_um,
                                               int n_steps) {
    if (T_K <= 0.0 || lambda1_um <= 0.0 || lambda2_um <= lambda1_um)
        return -1.0;
    if (n_steps < 10) n_steps = 100;

    double lambda1_m = lambda1_um * 1e-6;
    double lambda2_m = lambda2_um * 1e-6;
    double delta_T = T_K * 0.001;  /* 0.1% perturbation */

    double L1 = ir_in_band_radiance(lambda1_m, lambda2_m, T_K, n_steps);
    double L2 = ir_in_band_radiance(lambda1_m, lambda2_m, T_K + delta_T, n_steps);

    if (L1 < 0.0 || L2 < 0.0) return -1.0;
    return (L2 - L1) / delta_T;  /* dL/dT [W/(sr*m2*K)] */
}

double ir_netd(double f_number, double tau_opt,
                double detector_area_cm2, double nep_W,
                double dL_dT) {
    if (f_number <= 0.0 || tau_opt <= 0.0 || tau_opt > 1.0) return -1.0;
    if (detector_area_cm2 <= 0.0 || nep_W <= 0.0 || dL_dT <= 0.0) return -1.0;

    double A_d_m2 = detector_area_cm2 * 1e-4;  /* cm2 to m2 */
    double F2 = f_number * f_number;

    return 4.0 * F2 * nep_W / (tau_opt * A_d_m2 * dL_dT);
}

/* =========================================================================
 * MRTD - Minimum Resolvable Temperature Difference
 *
 * MRTD combines NETD with spatial resolution (MTF) and human perception.
 * It determines whether a standard 4-bar target can be distinguished.
 *
 * Simplified MRTD (Lloyd, 1975):
 *   MRTD(f) = (pi^2/8) * NETD * (f/MTF(f)) * sqrt(alpha*beta/(t_eye*F*tau))
 *
 * where f = spatial frequency, MTF(f) = modulation transfer function,
 * alpha,beta = IFOV (angular), t_eye = eye integration (0.1-0.2 s),
 * F = frame rate, tau = atmospheric transmission
 * ========================================================================= */

double ir_mrtd(double netd_K, double spatial_freq_cyc_per_mrad,
                double mtf_value, double ifov_mrad,
                double eye_integration_s, double frame_rate_Hz) {
    if (netd_K <= 0.0 || spatial_freq_cyc_per_mrad <= 0.0) return -1.0;
    if (mtf_value <= 0.0 || mtf_value > 1.0) return -1.0;
    if (ifov_mrad <= 0.0 || eye_integration_s <= 0.0 || frame_rate_Hz <= 0.0)
        return -1.0;

    double pi_factor = M_PI * M_PI / 8.0;  /* pi^2/8 from bar pattern analysis */
    double freq_mtf_ratio = spatial_freq_cyc_per_mrad / mtf_value;
    double geo_factor = ifov_mrad * ifov_mrad;
    double temporal_factor = eye_integration_s * frame_rate_Hz;

    if (temporal_factor <= 0.0) return -1.0;

    return pi_factor * netd_K * freq_mtf_ratio
           * sqrt(geo_factor / temporal_factor);
}

/* =========================================================================
 * Microbolometer Thermal Model
 *
 * Temperature rise: Delta_T = eta * P_inc / (G_th * sqrt(1 + (omega*tau)^2))
 * Time constant:    tau_th = C_th / G_th
 *
 * Microbolometers are thermal detectors: absorbed IR heats a thin membrane,
 * changing its electrical resistance (TCR effect). The thermal isolation
 * structure determines responsivity and speed.
 *
 * Reference: Wood, R.A. (1993) "Uncooled thermal imaging with monolithic
 *   silicon focal planes", SPIE Vol. 2020
 * ========================================================================= */

double ir_microbolometer_temp_rise(const ir_detector_params_t *params,
                                    double incident_power_W, double freq_Hz) {
    if (!params || incident_power_W < 0.0 || freq_Hz < 0.0) return -1.0;
    if (params->thermal_conductance_W_per_K <= 0.0) return -1.0;

    double G_th = params->thermal_conductance_W_per_K;
    double tau = params->thermal_time_constant_ms * 1e-3;
    double omega = 2.0 * M_PI * freq_Hz;

    double denom = G_th * sqrt(1.0 + omega * omega * tau * tau);
    return params->quantum_efficiency * incident_power_W / denom;
}

double ir_microbolometer_time_constant(const ir_detector_params_t *params) {
    if (!params || params->thermal_conductance_W_per_K <= 0.0) return -1.0;
    return params->thermal_mass_J_per_K / params->thermal_conductance_W_per_K;
}

/* =========================================================================
 * Quantum Efficiency and Cutoff Wavelength
 *
 * External QE: eta_ext = (1-R) * (1-exp(-alpha*d)) * eta_int
 * Cutoff:       lambda_c = 1.24 / Eg  [um]
 *
 * The cutoff wavelength determines the longest wavelength a photon
 * detector can detect. Photons with energy below the bandgap (lambda > lambda_c)
 * pass through without generating electron-hole pairs.
 * ========================================================================= */

double ir_quantum_efficiency_external(double reflectance,
                                       double absorption_coeff_cm_inv,
                                       double thickness_cm,
                                       double internal_qe) {
    if (reflectance < 0.0 || reflectance >= 1.0) return -1.0;
    if (absorption_coeff_cm_inv <= 0.0 || thickness_cm <= 0.0) return -1.0;
    if (internal_qe < 0.0 || internal_qe > 1.0) return -1.0;

    double transmission = 1.0 - reflectance;
    double absorption = 1.0 - exp(-absorption_coeff_cm_inv * thickness_cm);
    return transmission * absorption * internal_qe;
}

double ir_cutoff_wavelength(double bandgap_eV) {
    if (bandgap_eV <= 0.0) return -1.0;
    /* lambda_c [um] = h*c/Eg = 1.2398/Eg */
    double hc_eV_um = 1.239841984;  /* h*c in eV*um */
    return hc_eV_um / bandgap_eV;
}
