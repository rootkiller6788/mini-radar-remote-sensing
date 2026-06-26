/**
 * @file    ir_core.h
 * @brief   Infrared Thermal Imaging - Core Physical Laws and Definitions
 *
 * Nine-Level Knowledge Coverage:
 *   L1 - Definitions: blackbody, emissivity, spectral radiance, IR bands
 *   L2 - Core Concepts: gray body, thermal contrast, photon flux, BLIP
 *   L3 - Math Structures: Planck function (lambda, nu, nutilde, photon domains)
 *   L4 - Fundamental Laws: Planck, Wien displacement, Stefan-Boltzmann,
 *        Rayleigh-Jeans, Wien approximation, Kirchhoff's law
 *
 * Primary References:
 *   Planck, M. (1901) Annalen der Physik, 4(3), 553-563
 *   Wien, W. (1896) Annalen der Physik, 294(8), 662-669
 *   Stefan, J. (1879) / Boltzmann, L. (1884) T^4 law
 *   Kirchhoff, G. (1860) thermal radiation law
 *   Richards, Scheer & Holm (2010) "Principles of Modern Radar", Ch.9
 *   Rogalski, A. (2011) "Infrared Detectors", 2nd Ed.
 *
 * Nine-School Curriculum Mapping:
 *   MIT 6.630 EM Waves - blackbody radiation, Planck spectrum
 *   Stanford EE247 Optical - radiometry, IR detectors
 *   Berkeley EE117 EM - thermal radiation fundamentals
 *   Illinois ECE 451 EM - radiation from thermal sources
 *   Michigan EECS 411 Microwave - mm-wave/IR radiometry
 *   Georgia Tech ECE 6350 EM - thermal emission modeling
 *   TU Munich High-Frequency Eng - IR sensor physics
 *   ETH 227-0455 EM - thermal radiation systems
 *   Tsinghua EM - blackbody radiation theory
 */

#ifndef IR_CORE_H
#define IR_CORE_H

#include <stddef.h>
#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L1: Fundamental Physical Constants (CODATA 2018 recommended values)
 * ======================================================================== */

/** Planck constant h [J*s] - quantum of action */
#define IR_PLANCK_CONSTANT        6.62607015e-34

/** Boltzmann constant k_B [J/K] */
#define IR_BOLTZMANN_CONSTANT     1.380649e-23

/** Speed of light in vacuum c [m/s] - exact defined value */
#define IR_SPEED_OF_LIGHT         2.99792458e8

/** Stefan-Boltzmann constant sigma [W/(m2*K4)]
 *  sigma = 2*pi^5*kB^4/(15*h^3*c^2) = 5.670374419e-8 */
#define IR_STEFAN_BOLTZMANN       5.670374419e-8

/** Wien displacement constant b = lambda_max*T [m*K]
 *  b = h*c/(kB*4.965114...) = 2.897771955e-3 */
#define IR_WIEN_DISPLACEMENT      2.897771955e-3

/** First radiation constant c1 = 2*pi*h*c^2 [W*m2] */
#define IR_FIRST_RADIATION_C1     (2.0 * M_PI * IR_PLANCK_CONSTANT                                    * IR_SPEED_OF_LIGHT * IR_SPEED_OF_LIGHT)

/** Second radiation constant c2 = h*c/kB [m*K] */
#define IR_SECOND_RADIATION_C2    ((IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT)                                    / IR_BOLTZMANN_CONSTANT)

/** Elementary charge e [C] */
#define IR_ELEMENTARY_CHARGE      1.602176634e-19

/** Photons per second per watt at 1 um wavelength */
#define IR_PHOTONS_PER_SEC_PER_WATT_AT_1UM  5.034117e18

/* ========================================================================
 * L1: IR Spectral Bands (ISO 20473:2007)
 * ======================================================================== */

typedef enum {
    IR_BAND_NIR   = 0,   /**< Near-Infrared:       0.75 - 1.4 um   */
    IR_BAND_SWIR  = 1,   /**< Short-Wave Infrared:  1.4  - 3.0 um   */
    IR_BAND_MWIR  = 2,   /**< Mid-Wave Infrared:    3.0  - 5.0 um   */
    IR_BAND_LWIR  = 3,   /**< Long-Wave Infrared:   8.0  - 14.0 um  */
    IR_BAND_VLWIR = 4,   /**< Very-Long-Wave IR:    14.0 - 30.0 um  */
    IR_BAND_FIR   = 5,   /**< Far-Infrared:         30.0 - 1000 um  */
    IR_BAND_COUNT = 6
} ir_band_t;

/** Band range descriptor with detector and application information */
typedef struct {
    double      lambda_min_um;
    double      lambda_max_um;
    const char *name;
    ir_band_t   band;
    const char *common_detectors;
    const char *typical_applications;
} ir_band_range_t;

const ir_band_range_t* ir_band_get_range(ir_band_t band);
void ir_band_print_all(void);

/* ========================================================================
 * L1: Blackbody Radiator Definition
 * ======================================================================== */

typedef struct {
    double temperature_K;    /**< absolute temperature [K] */
    double emissivity;       /**< emissivity epsilon in [0, 1] */
    double area_m2;          /**< radiating surface area [m2] */
    int    is_lambertian;    /**< 1 = Lambertian (diffuse) emitter */
    double reflectivity;     /**< reflectivity rho in [0, 1] */
    double transmissivity;   /**< transmissivity tau in [0, 1] */
} ir_blackbody_t;

int ir_blackbody_init(ir_blackbody_t *bb, double T, double eps,
                      double area, int lambertian);

/* ========================================================================
 * L3 & L4: Planck's Law of Blackbody Radiation (1901)
 *
 * B_lambda(lambda,T) = (2*h*c^2/lambda^5) / (exp(h*c/(lambda*kB*T)) - 1)
 * Units: [W*sr^-1*m^-3]
 *
 * Planck's law resolved the ultraviolet catastrophe and introduced
 * energy quantization E = h*nu. All thermal radiation laws derive from it.
 *
 * Domain representations (energy conservation: B_lambda*d_lambda = B_nu*d_nu):
 *   Wavelength: B_lambda(lambda,T)
 *   Frequency:  B_nu(nu,T) = 2*h*nu^3/c^2 / (exp(h*nu/(kB*T))-1)
 *   Wavenumber: B_nutilde(nutilde,T) = 2*h*c^2*nutilde^3/(exp(h*c*nutilde/(kB*T))-1)
 *   Photon:     B_q(lambda,T) = 2*c/lambda^4 / (exp(h*c/(lambda*kB*T))-1)
 * ======================================================================== */

/** Planck spectral radiance - wavelength domain [W*sr^-1*m^-3] */
double ir_planck_spectral_radiance_wavelength(double lambda_m, double T_K);

/** Planck spectral radiance - frequency domain [W*sr^-1*m^-2*Hz^-1] */
double ir_planck_spectral_radiance_frequency(double nu_Hz, double T_K);

/** Planck spectral radiance - wavenumber domain [W*sr^-1*m^-1] */
double ir_planck_spectral_radiance_wavenumber(double nu_tilde_per_m, double T_K);

/** Photon spectral radiance [photons*s^-1*sr^-1*m^-3] */
double ir_planck_photon_radiance(double lambda_m, double T_K);

/** Peak spectral radiance at Wien peak wavelength [W*sr^-1*m^-3] */
double ir_planck_peak_radiance(double T_K);

/* ========================================================================
 * L4: Wien's Displacement Law (1893)
 *
 * lambda_max * T = b = 2897.8 um*K
 * nu_max = 5.87893e10 * T [Hz/K]
 *
 * Describes shift of peak emission to shorter wavelengths at higher T.
 * Derived from d/d_lambda of Planck function = 0.
 * ======================================================================== */

/** Wien peak wavelength [m]: lambda_max = b/T */
double ir_wien_peak_wavelength(double T_K);

/** Wien peak frequency [Hz]: nu_max = alpha*kB*T/h */
double ir_wien_peak_frequency(double T_K);

/* ========================================================================
 * L4: Stefan-Boltzmann Law (1879/1884)
 *
 * M(T) = sigma * T^4    [W/m2]
 * P = epsilon * sigma * A * T^4    [W]
 *
 * Total radiant exitance integrated over all wavelengths.
 * Derived by integrating Planck function over all lambda.
 * ======================================================================== */

/** Stefan-Boltzmann radiant exitance [W/m2] */
double ir_stefan_boltzmann_exitance(double T_K);

/** Total radiated power from a blackbody surface [W] */
double ir_blackbody_total_power(const ir_blackbody_t *bb);

/** Verify Stefan-Boltzmann integral numerically.
 *  Returns |numerical - exact| / exact. */
double ir_stefan_boltzmann_verify(double T_K, int n_steps);

/* ========================================================================
 * L4: Rayleigh-Jeans Law (1900) - classical long-wavelength limit
 *
 * B_RJ(lambda,T) = 2*c*kB*T / lambda^4
 * Valid when h*c/(lambda*kB*T) << 1.
 * Classical EM prediction; fails at short wavelengths (UV catastrophe).
 * ======================================================================== */

/** Rayleigh-Jeans approximation [W*sr^-1*m^-3] */
double ir_rayleigh_jeans_radiance(double lambda_m, double T_K);

/** RJ brightness temperature [K]: T_B = B*lambda^4/(2*c*kB) */
double ir_rayleigh_jeans_temperature(double lambda_m, double radiance);

/* ========================================================================
 * L4: Wien Approximation (1896) - short-wavelength limit
 *
 * B_Wien = (2*h*c^2/lambda^5) * exp(-h*c/(lambda*kB*T))
 * Accurate within 1%% when lambda*T < 3000 um*K.
 * ======================================================================== */

/** Wien approximation [W*sr^-1*m^-3] */
double ir_wien_approx_radiance(double lambda_m, double T_K);

/* ========================================================================
 * L2: Emissivity and Kirchhoff's Law (1860)
 *
 * At thermal equilibrium: emissivity equals absorptivity.
 * Blackbody: epsilon=1; Gray body: epsilon=const<1; Selective: epsilon(lambda)
 * Kirchhoff: epsilon = 1 - R - T
 * ======================================================================== */

/** Gray-body spectral radiance: B_gray = epsilon * B_blackbody */
double ir_graybody_radiance(double lambda_m, double T_K, double emissivity);

/** Compute emissivity from reflectance and transmittance */
double ir_kirchhoff_emissivity(double reflectance, double transmittance);

/* ========================================================================
 * L3: Band-Integrated Radiometric Quantities
 *
 * In-band radiance: L_band = integral_{lambda1}^{lambda2} B_lambda(lambda,T) d_lambda
 * Uses Simpson's rule numerical integration.
 * M_band = pi * L_band (Lambertian approximation)
 * ======================================================================== */

/** In-band radiance by Simpson integration [W*sr^-1*m^-2] */
double ir_in_band_radiance(double lambda1_m, double lambda2_m,
                            double T_K, int n_steps);

/** In-band radiant exitance [W/m2]: M_band = pi * L_band */
double ir_in_band_exitance(double lambda1_m, double lambda2_m,
                            double T_K, int n_steps);

/** In-band photon flux for photon detector analysis [ph*s^-1*sr^-1*m^-2] */
double ir_photon_flux_band(double lambda1_m, double lambda2_m,
                            double T_K, int n_steps);

/* ========================================================================
 * L2: Thermal Contrast and Radiometer Sensitivity
 * ======================================================================== */

/** Thermal contrast: C = (M(T1)-M(T2)) / M(T2) [dimensionless] */
double ir_thermal_contrast(double T_target_K, double T_bg_K);

/** Apparent temperature difference: dT_app = dT * tau_atm [K] */
double ir_apparent_temperature_difference(double delta_T, double tau_atm);

/** Radiometer noise-equivalent temperature difference [K]
 *  NE_dT = T_sys / sqrt(B * tau_int) */
double ir_radiometer_NE_dT(double T_sys_K, double bandwidth_Hz,
                            double integration_s);

/* ========================================================================
 * L2: BLIP - Background-Limited Infrared Photodetector
 *
 * Fundamental sensitivity limit set by background photon shot noise.
 * D*_BLIP = lambda/(h*c) * sqrt(eta/(2*Phi_bg))
 * ======================================================================== */

/** BLIP-limited specific detectivity [cm*sqrt(Hz)/W] (Jones) */
double ir_blip_d_star(double quantum_eff, double wavelength_um,
                       double photon_flux_bg);

/** Photon irradiance from a blackbody scene [ph/(s*m2)] */
double ir_photon_irradiance(double T_K, double lambda1_m, double lambda2_m,
                             double solid_angle_sr, int n_steps);

/* ========================================================================
 * L3: Inverse Planck - Brightness Temperature
 *
 * T_B = h*c / (lambda*kB * ln(1 + 2*h*c^2/(lambda^5 * B_lambda)))
 *
 * Converts measured radiance to equivalent blackbody temperature.
 * ======================================================================== */

/** Brightness temperature from spectral radiance [K] */
double ir_brightness_temperature(double lambda_m, double radiance);

/** Blackbody fractional function F(lambda*T) in [0,1]
 *  Fraction of total exitance emitted in [0, lambda].
 *  Uses Chang & Rhee (1984) infinite series approximation. */
double ir_blackbody_fraction(double lambda_m, double T_K);

/** Dual-band exitance ratio for ratio thermometry.
 *  R = M_band1(T) / M_band2(T).
 *  Provides emissivity-independent temperature measurement. */
double ir_dual_band_ratio(double T_K,
                           double lambda1_1_m, double lambda1_2_m,
                           double lambda2_1_m, double lambda2_2_m,
                           int n_steps);

#endif /* IR_CORE_H */
