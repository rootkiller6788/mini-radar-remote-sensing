/**
 * @file    ir_radiometry.h
 * @brief   Infrared Radiometry - Radiative Transfer and IR Range Equation
 *
 * L1 - Definitions: radiance, irradiance, radiant intensity, etendue
 * L2 - Core Concepts: optical throughput, range equation, atmospheric attenuation
 * L3 - Math Structures: radiative transfer equation, solid angle integrals
 * L4 - Fundamental Laws: conservation of radiance, Lambert-Beer law
 *
 * References:
 *   Wolfe, W.L. (1998) "Introduction to Radiometry", SPIE Press
 *   Accetta & Shumaker (1993) "The IR & EO Systems Handbook"
 *
 * Curriculum: MIT 6.630 EM, Stanford EE247 Optical, Berkeley EE117 EM,
 *             TU Munich High-Frequency, ETH 227-0455 EM
 */

#ifndef IR_RADIOMETRY_H
#define IR_RADIOMETRY_H

#include <stddef.h>
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "ir_core.h"

/* Fundamental radiometric quantities */
typedef struct {
    double radiance;            /* L [W/(sr*m2)] */
    double irradiance;          /* E [W/m2] */
    double radiant_intensity;   /* I [W/sr] */
    double radiant_power;       /* Phi [W] */
} ir_radiometric_quantities_t;

/* Source geometry classification */
typedef enum {
    IR_SOURCE_POINT = 0,        /* unresolved (angular size << IFOV) */
    IR_SOURCE_EXTENDED,         /* resolved (fills FOV) */
    IR_SOURCE_LINE              /* linear source */
} ir_source_type_t;

typedef struct {
    ir_source_type_t type;
    double temperature_K;
    double emissivity;
    double area_m2;
    double solid_angle_sr;
    double range_m;
} ir_source_t;

/* Optical system model */
typedef struct {
    double aperture_diameter_m;
    double focal_length_m;
    double f_number;
    double ifov_horizontal_mrad;
    double ifov_vertical_mrad;
    double optical_transmission;
    double detector_area_cm2;
} ir_optical_system_t;

/* Radiative transfer: source -> atmosphere -> optics -> detector */
double ir_source_radiance_at_sensor(const ir_source_t *source,
                                     double tau_atmosphere, double path_radiance);
double ir_point_source_irradiance(const ir_source_t *source,
                                   const ir_optical_system_t *optics,
                                   double tau_atmosphere);
double ir_extended_source_irradiance(const ir_source_t *source,
                                      const ir_optical_system_t *optics,
                                      double tau_atmosphere);
double ir_power_on_detector(const ir_source_t *source,
                             const ir_optical_system_t *optics,
                             double tau_atmosphere);

/* Etendue (geometric extent) and optical throughput */
double ir_etendue_object_space(double aperture_area_m2, double solid_angle_sr);
double ir_etendue_image_space(double detector_area_m2, double f_number);
double ir_optical_throughput(double radiance, double etendue, double tau);

/* Solid angle calculations */
double ir_solid_angle_circular(double diameter_m, double range_m);
double ir_solid_angle_rectangular(double width_m, double height_m, double range_m);
double ir_ifov_solid_angle(double ifov_h_mrad, double ifov_v_mrad);

/* IR Range Equation - SNR models */
double ir_snr_point_source(double source_intensity_W_per_sr,
                            double tau_atm, double tau_opt,
                            double aperture_area_m2,
                            double responsivity_V_per_W,
                            double noise_voltage, double range_m);
double ir_snr_extended_source(double source_radiance,
                               double tau_atm, double tau_opt,
                               double detector_area_cm2, double f_number,
                               double responsivity_V_per_W, double noise_voltage);
double ir_max_detection_range(double source_intensity_W_per_sr,
                               double tau_atm, double tau_opt,
                               double aperture_area_m2,
                               double responsivity_V_per_W,
                               double noise_voltage, double snr_min);

/* Temperature-radiance conversion */
double ir_temperature_to_radiance(double T_K, double lambda1_um,
                                   double lambda2_um, int n_steps);
double ir_radiance_to_temperature(double radiance, double lambda1_um,
                                   double lambda2_um, int n_steps,
                                   double T_guess_K);
double ir_radiance_temperature_derivative(double T_K, double lambda1_um,
                                           double lambda2_um, int n_steps);

#endif /* IR_RADIOMETRY_H */
