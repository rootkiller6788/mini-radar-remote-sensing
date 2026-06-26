/**
 * @file    ir_radiometry.c
 * @brief   Radiative Transfer and IR Range Equation - Implementation
 *
 * Implements the complete radiometric chain: source to detector propagation,
 * etendue, solid angle calculations, SNR and maximum detection range.
 *
 * References: Wolfe (1998) "Introduction to Radiometry", SPIE Press
 */

#include "ir_radiometry.h"
#include "ir_core.h"
#include <stdio.h>
#include <math.h>

/* ========================================================================
 * Radiative Transfer: Source -> Atmosphere -> Sensor
 * L_sensor = tau_atm * L_source + L_path
 * ======================================================================== */

double ir_source_radiance_at_sensor(const ir_source_t *source,
                                     double tau_atmosphere,
                                     double path_radiance) {
    if (!source) return -1.0;
    if (tau_atmosphere < 0.0 || tau_atmosphere > 1.0) return -1.0;

    double L_source = ir_stefan_boltzmann_exitance(source->temperature_K)
                      * source->emissivity / M_PI;
    return tau_atmosphere * L_source + path_radiance;
}

double ir_point_source_irradiance(const ir_source_t *source,
                                   const ir_optical_system_t *optics,
                                   double tau_atmosphere) {
    if (!source || !optics) return -1.0;
    if (source->range_m <= 0.0) return -1.0;

    double I_src = source->emissivity
                   * ir_stefan_boltzmann_exitance(source->temperature_K)
                   * source->area_m2 / M_PI;
    double A_aperture = M_PI * optics->aperture_diameter_m
                        * optics->aperture_diameter_m / 4.0;
    double A_det = optics->detector_area_cm2 * 1e-4;
    if (A_det <= 0.0) return -1.0;

    return I_src * tau_atmosphere * optics->optical_transmission
           * A_aperture / (source->range_m * source->range_m * A_det);
}

double ir_extended_source_irradiance(const ir_source_t *source,
                                      const ir_optical_system_t *optics,
                                      double tau_atmosphere) {
    if (!source || !optics) return -1.0;
    if (optics->f_number <= 0.0) return -1.0;

    double L_src = source->emissivity
                   * ir_stefan_boltzmann_exitance(source->temperature_K) / M_PI;
    double F2 = optics->f_number * optics->f_number;
    return L_src * tau_atmosphere * optics->optical_transmission
           * M_PI / (4.0 * F2);
}

double ir_power_on_detector(const ir_source_t *source,
                             const ir_optical_system_t *optics,
                             double tau_atmosphere) {
    if (!source || !optics) return -1.0;
    double E;
    if (source->type == IR_SOURCE_POINT) {
        E = ir_point_source_irradiance(source, optics, tau_atmosphere);
    } else {
        E = ir_extended_source_irradiance(source, optics, tau_atmosphere);
    }
    if (E < 0.0) return -1.0;
    return E * optics->detector_area_cm2 * 1e-4;
}

/* ========================================================================
 * Etendue and Optical Throughput
 * G = n^2 * A * Omega   [m2*sr]
 * Phi = L * G * tau      [W]
 * ======================================================================== */

double ir_etendue_object_space(double aperture_area_m2, double solid_angle_sr) {
    if (aperture_area_m2 <= 0.0 || solid_angle_sr <= 0.0) return -1.0;
    return aperture_area_m2 * solid_angle_sr;
}

double ir_etendue_image_space(double detector_area_m2, double f_number) {
    if (detector_area_m2 <= 0.0 || f_number <= 0.0) return -1.0;
    return M_PI * detector_area_m2 / (4.0 * f_number * f_number);
}

double ir_optical_throughput(double radiance, double etendue,
                              double optical_transmission) {
    if (radiance < 0.0 || etendue <= 0.0) return -1.0;
    if (optical_transmission < 0.0 || optical_transmission > 1.0) return -1.0;
    return radiance * etendue * optical_transmission;
}

/* ========================================================================
 * Solid Angle Calculations
 * Omega = 2*pi*(1-cos(theta_half)) = pi*sin^2(theta_half) (small angle)
 * ======================================================================== */

double ir_solid_angle_circular(double diameter_m, double range_m) {
    if (diameter_m <= 0.0 || range_m <= 0.0) return -1.0;
    double half_angle = atan(diameter_m / (2.0 * range_m));
    return 2.0 * M_PI * (1.0 - cos(half_angle));
}

double ir_solid_angle_rectangular(double width_m, double height_m,
                                   double range_m) {
    if (width_m <= 0.0 || height_m <= 0.0 || range_m <= 0.0) return -1.0;
    double A = width_m * height_m;
    return A / (range_m * range_m);  /* small angle approximation */
}

double ir_ifov_solid_angle(double ifov_h_mrad, double ifov_v_mrad) {
    if (ifov_h_mrad <= 0.0 || ifov_v_mrad <= 0.0) return -1.0;
    return ifov_h_mrad * ifov_v_mrad * 1e-6;  /* mrad^2 to sr */
}

/* ========================================================================
 * IR Range Equation - SNR Models
 *
 * Point source: SNR = I_src * tau_atm * tau_opt * A_aperture * R_v / (R^2 * V_n)
 * Extended:     SNR = L_src * tau_atm * tau_opt * pi * A_det * R_v / (4 * F^2 * V_n)
 * ======================================================================== */

double ir_snr_point_source(double source_intensity_W_per_sr,
                            double tau_atm, double tau_opt,
                            double aperture_area_m2,
                            double responsivity_V_per_W,
                            double noise_voltage, double range_m) {
    if (source_intensity_W_per_sr <= 0.0 || range_m <= 0.0) return -1.0;
    if (tau_atm < 0.0 || tau_opt < 0.0 || noise_voltage <= 0.0) return -1.0;
    if (responsivity_V_per_W <= 0.0 || aperture_area_m2 <= 0.0) return -1.0;

    double irradiance = source_intensity_W_per_sr * tau_atm * tau_opt
                        * aperture_area_m2 / (range_m * range_m);
    double signal_voltage = irradiance * responsivity_V_per_W;
    return signal_voltage / noise_voltage;
}

double ir_snr_extended_source(double source_radiance,
                               double tau_atm, double tau_opt,
                               double detector_area_cm2, double f_number,
                               double responsivity_V_per_W,
                               double noise_voltage) {
    if (source_radiance <= 0.0 || f_number <= 0.0) return -1.0;
    if (tau_atm < 0.0 || tau_opt < 0.0 || noise_voltage <= 0.0) return -1.0;
    if (detector_area_cm2 <= 0.0 || responsivity_V_per_W <= 0.0) return -1.0;

    double A_det_m2 = detector_area_cm2 * 1e-4;
    double F2 = f_number * f_number;
    double power = source_radiance * tau_atm * tau_opt * M_PI
                   * A_det_m2 / (4.0 * F2);
    return power * responsivity_V_per_W / noise_voltage;
}

double ir_max_detection_range(double source_intensity_W_per_sr,
                               double tau_atm, double tau_opt,
                               double aperture_area_m2,
                               double responsivity_V_per_W,
                               double noise_voltage, double snr_min) {
    if (source_intensity_W_per_sr <= 0.0 || snr_min <= 0.0) return -1.0;
    if (tau_atm < 0.0 || tau_opt < 0.0 || noise_voltage <= 0.0) return -1.0;
    if (responsivity_V_per_W <= 0.0 || aperture_area_m2 <= 0.0) return -1.0;

    double numerator = source_intensity_W_per_sr * tau_atm * tau_opt
                       * aperture_area_m2 * responsivity_V_per_W;
    return sqrt(numerator / (snr_min * noise_voltage));
}

/* ========================================================================
 * Temperature-Radiance Conversion
 * ======================================================================== */

double ir_temperature_to_radiance(double T_K, double lambda1_um,
                                   double lambda2_um, int n_steps) {
    if (T_K <= 0.0) return -1.0;
    double l1_m = lambda1_um * 1e-6;
    double l2_m = lambda2_um * 1e-6;
    return ir_in_band_radiance(l1_m, l2_m, T_K, n_steps);
}

double ir_radiance_to_temperature(double radiance, double lambda1_um,
                                   double lambda2_um, int n_steps,
                                   double T_guess_K) {
    if (radiance <= 0.0 || T_guess_K <= 0.0) return -1.0;

    double l1_m = lambda1_um * 1e-6;
    double l2_m = lambda2_um * 1e-6;

    /* Bisection method for inverse Planck integration */
    double T_lo = T_guess_K * 0.5;
    double T_hi = T_guess_K * 2.0;
    int max_iter = 50;
    double tol = 0.01;

    for (int i = 0; i < max_iter; i++) {
        double T_mid = (T_lo + T_hi) / 2.0;
        double L_mid = ir_in_band_radiance(l1_m, l2_m, T_mid, n_steps);
        if (L_mid < 0.0) return -1.0;
        if (fabs(L_mid - radiance) < tol * radiance / 100.0) return T_mid;
        if ((T_hi - T_lo) < tol) return T_mid;
        if (L_mid < radiance) T_lo = T_mid;
        else T_hi = T_mid;
    }
    return (T_lo + T_hi) / 2.0;
}

double ir_radiance_temperature_derivative(double T_K, double lambda1_um,
                                           double lambda2_um, int n_steps) {
    if (T_K <= 0.0) return -1.0;
    double dT = T_K * 0.001;
    double L1 = ir_temperature_to_radiance(T_K, lambda1_um, lambda2_um, n_steps);
    double L2 = ir_temperature_to_radiance(T_K + dT, lambda1_um, lambda2_um, n_steps);
    if (L1 < 0.0 || L2 < 0.0) return -1.0;
    return (L2 - L1) / dT;
}
