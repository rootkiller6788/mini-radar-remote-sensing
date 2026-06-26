/**
 * @file    ir_core.c
 * @brief   Planck, Wien, Stefan-Boltzmann law implementations
 *
 * Reference: Planck (1901), Wien (1896), Stefan (1879), Boltzmann (1884)
 */

#include "ir_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Band database */
static const ir_band_range_t g_bands[] = {
    { 0.75,  1.4,  "NIR",   IR_BAND_NIR,   "Si, InGaAs",        "Night vision, telecom" },
    { 1.4,   3.0,  "SWIR",  IR_BAND_SWIR,  "InGaAs, HgCdTe",    "Surveillance, agriculture" },
    { 3.0,   5.0,  "MWIR",  IR_BAND_MWIR,  "InSb, HgCdTe, PbSe", "Missile warning, gas detection" },
    { 8.0,   14.0, "LWIR",  IR_BAND_LWIR,  "VOx, a-Si, HgCdTe",  "Thermal imaging, thermography" },
    { 14.0,  30.0, "VLWIR", IR_BAND_VLWIR, "Si:As, Ge:Ga",       "Astronomy, spectroscopy" },
    { 30.0, 1000., "FIR",   IR_BAND_FIR,   "Bolometers, Golay",   "Radio astronomy, plasma diagnostics" }
};

const ir_band_range_t* ir_band_get_range(ir_band_t band) {
    if (band < 0 || band >= IR_BAND_COUNT) return NULL;
    return &g_bands[band];
}

void ir_band_print_all(void) {
    printf("IR Spectral Bands (ISO 20473:2007):\\n");
    for (int i = 0; i < IR_BAND_COUNT; i++) {
        printf("  %s: %.2f - %.2f um | Detectors: %s | Use: %s\\n",               g_bands[i].name, g_bands[i].lambda_min_um,               g_bands[i].lambda_max_um,               g_bands[i].common_detectors,               g_bands[i].typical_applications);
    }
}

int ir_blackbody_init(ir_blackbody_t *bb, double T, double eps,
                      double area, int lambertian) {
    if (!bb) return -1;
    if (T <= 0.0) return -1;
    if (eps < 0.0 || eps > 1.0) return -1;
    if (area <= 0.0) return -1;
    bb->temperature_K = T;
    bb->emissivity = eps;
    bb->area_m2 = area;
    bb->is_lambertian = lambertian;
    bb->reflectivity = 1.0 - eps;
    bb->transmissivity = 0.0;
    return 0;
}

/* =========================================================================
 * Planck spectral radiance - wavelength domain [W*sr^-1*m^-3]
 *
 * B_lambda(lambda,T) = 2*h*c^2/lambda^5 * 1/(exp(h*c/(lambda*kB*T)) - 1)
 *
 * Computational strategy:
 *   Let x = h*c/(lambda*kB*T)
 *   If x < 1e-6: use Rayleigh-Jeans limit B ~ 2*c*kB*T/lambda^4
 *   If x > 50:   use Wien approximation B ~ 2*h*c^2/lambda^5 * exp(-x)
 *   Otherwise:   compute directly with expm1 for numerical stability
 * ========================================================================= */
double ir_planck_spectral_radiance_wavelength(double lambda_m, double T_K) {
    if (lambda_m <= 0.0 || T_K <= 0.0) return -1.0;

    double c1 = IR_FIRST_RADIATION_C1;
    double c2 = IR_SECOND_RADIATION_C2;
    double x = c2 / (lambda_m * T_K);

    /* For very small x, use Rayleigh-Jeans limit to avoid 1/0 */
    if (x < 1e-7) {
        return (2.0 * IR_SPEED_OF_LIGHT * IR_BOLTZMANN_CONSTANT * T_K)
               / (lambda_m * lambda_m * lambda_m * lambda_m);
    }

    /* For large x, use Wien approximation to avoid overflow */
    if (x > 100.0) {
        return (c1 / (lambda_m * lambda_m * lambda_m * lambda_m * lambda_m))
               * exp(-x);
    }

    /* Direct Planck computation with expm1 for stability */
    double denom = expm1(x);
    if (denom <= 0.0) return -1.0;
    return (c1 / (lambda_m * lambda_m * lambda_m * lambda_m * lambda_m))
           / denom;
}

/* =========================================================================
 * Planck spectral radiance - frequency domain [W*sr^-1*m^-2*Hz^-1]
 *
 * B_nu(nu,T) = 2*h*nu^3/c^2 * 1/(exp(h*nu/(kB*T))-1)
 * ========================================================================= */
double ir_planck_spectral_radiance_frequency(double nu_Hz, double T_K) {
    if (nu_Hz <= 0.0 || T_K <= 0.0) return -1.0;

    double h = IR_PLANCK_CONSTANT;
    double kB = IR_BOLTZMANN_CONSTANT;
    double c = IR_SPEED_OF_LIGHT;
    double x = h * nu_Hz / (kB * T_K);

    if (x < 1e-7) {
        return 2.0 * nu_Hz * nu_Hz * kB * T_K / (c * c);
    }
    if (x > 100.0) {
        return 2.0 * h * nu_Hz * nu_Hz * nu_Hz / (c * c) * exp(-x);
    }

    double denom = expm1(x);
    if (denom <= 0.0) return -1.0;
    return (2.0 * h * nu_Hz * nu_Hz * nu_Hz / (c * c)) / denom;
}

/* =========================================================================
 * Planck spectral radiance - wavenumber domain [W*sr^-1*m^-1]
 *
 * B_nutilde = 2*h*c^2*nutilde^3 / (exp(h*c*nutilde/(kB*T)) - 1)
 * ========================================================================= */
double ir_planck_spectral_radiance_wavenumber(double nu_tilde_per_m, double T_K) {
    if (nu_tilde_per_m <= 0.0 || T_K <= 0.0) return -1.0;

    double hc = IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT;
    double x = hc * nu_tilde_per_m / (IR_BOLTZMANN_CONSTANT * T_K);
    double c1 = IR_FIRST_RADIATION_C1;

    if (x < 1e-7) {
        return 2.0 * IR_SPEED_OF_LIGHT * IR_BOLTZMANN_CONSTANT * T_K
               * nu_tilde_per_m * nu_tilde_per_m;
    }
    if (x > 100.0) {
        return c1 * nu_tilde_per_m * nu_tilde_per_m * nu_tilde_per_m * exp(-x);
    }

    double denom = expm1(x);
    if (denom <= 0.0) return -1.0;
    return (c1 * nu_tilde_per_m * nu_tilde_per_m * nu_tilde_per_m) / denom;
}

/* =========================================================================
 * Photon spectral radiance [photons*s^-1*sr^-1*m^-3]
 *
 * B_q = 2*c/lambda^4 * 1/(exp(h*c/(lambda*kB*T)) - 1)
 *
 * Important for photon-counting detectors where signal is proportional
 * to photon arrival rate, not energy flux.
 * ========================================================================= */
double ir_planck_photon_radiance(double lambda_m, double T_K) {
    if (lambda_m <= 0.0 || T_K <= 0.0) return -1.0;

    double x = IR_SECOND_RADIATION_C2 / (lambda_m * T_K);
    double prefactor = 2.0 * IR_SPEED_OF_LIGHT
                       / (lambda_m * lambda_m * lambda_m * lambda_m);

    if (x < 1e-7) {
        return prefactor * IR_BOLTZMANN_CONSTANT * T_K
               / (IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT) * lambda_m;
    }
    if (x > 100.0) {
        return prefactor * exp(-x);
    }

    double denom = expm1(x);
    if (denom <= 0.0) return -1.0;
    return prefactor / denom;
}

/* =========================================================================
 * Peak spectral radiance - evaluates Planck function at Wien peak
 * ========================================================================= */
double ir_planck_peak_radiance(double T_K) {
    if (T_K <= 0.0) return -1.0;
    double lambda_peak = ir_wien_peak_wavelength(T_K);
    if (lambda_peak <= 0.0) return -1.0;
    return ir_planck_spectral_radiance_wavelength(lambda_peak, T_K);
}

/* =========================================================================
 * Wien's Displacement Law: lambda_max = b / T
 *
 * b = 2.897771955...e-3 m*K
 *
 * Derivation: set d/d_lambda of Planck function = 0
 * Result: x = h*c/(lambda_max*kB*T) = 4.9651142317...
 *         lambda_max*T = h*c/(kB*4.965114...) = b
 * ========================================================================= */
double ir_wien_peak_wavelength(double T_K) {
    if (T_K <= 0.0) return -1.0;
    return IR_WIEN_DISPLACEMENT / T_K;
}

/* Wien peak in frequency domain: nu_max = 5.87893e10 * T [Hz/K] */
double ir_wien_peak_frequency(double T_K) {
    if (T_K <= 0.0) return -1.0;
    double alpha = 2.821439372122079;  /* solution of 3*(1-exp(-x)) = x */
    return alpha * IR_BOLTZMANN_CONSTANT * T_K / IR_PLANCK_CONSTANT;
}

/* =========================================================================
 * Stefan-Boltzmann Law: M = sigma * T^4
 *
 * sigma = 2*pi^5*kB^4 / (15*h^3*c^2) = 5.670374419e-8 W/(m2*K4)
 *
 * Total radiant exitance over all wavelengths and hemisphere.
 * ========================================================================= */
double ir_stefan_boltzmann_exitance(double T_K) {
    if (T_K <= 0.0) return -1.0;
    return IR_STEFAN_BOLTZMANN * T_K * T_K * T_K * T_K;
}

double ir_blackbody_total_power(const ir_blackbody_t *bb) {
    if (!bb || bb->temperature_K <= 0.0) return -1.0;
    double M = ir_stefan_boltzmann_exitance(bb->temperature_K);
    return bb->emissivity * M * bb->area_m2;
}

/* =========================================================================
 * Numerical verification of Stefan-Boltzmann by integrating Planck function
 *
 * Uses truncated integration: 0 to lambda_max where contribution < epsilon.
 * Transformation: lambda = b/x where b = Wien constant
 * Integral: M = pi * integral_0^inf B_lambda d_lambda
 * ========================================================================= */
double ir_stefan_boltzmann_verify(double T_K, int n_steps) {
    if (T_K <= 0.0 || n_steps < 10) return -1.0;

    /* Integration range: 10% to 1000% of Wien peak */
    double lambda_peak = IR_WIEN_DISPLACEMENT / T_K;
    double lambda_min = lambda_peak * 0.01;
    double lambda_max = lambda_peak * 100.0;
    double h = (lambda_max - lambda_min) / n_steps;

    double sum = 0.0;
    for (int i = 0; i <= n_steps; i++) {
        double lambda = lambda_min + i * h;
        double B = ir_planck_spectral_radiance_wavelength(lambda, T_K);
        if (B < 0.0) continue;
        double weight = (i == 0 || i == n_steps) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        sum += weight * B;
    }
    double L_total = (h / 3.0) * sum;   /* Simpson integrated radiance */
    double M_numerical = M_PI * L_total; /* Lambertian hemispherical exitance */
    double M_exact = IR_STEFAN_BOLTZMANN * T_K * T_K * T_K * T_K;

    return fabs(M_numerical - M_exact) / M_exact;
}

/* =========================================================================
 * Rayleigh-Jeans Law: B_RJ = 2*c*kB*T / lambda^4
 *
 * Classical limit: h*c/(lambda*kB*T) << 1
 * Valid when lambda*T > 0.77 m*K (microwave/radio domain at room temp).
 * Historical note: This is the classical EM prediction that leads to
 * the ultraviolet catastrophe at short wavelengths.
 * ========================================================================= */
double ir_rayleigh_jeans_radiance(double lambda_m, double T_K) {
    if (lambda_m <= 0.0 || T_K <= 0.0) return -1.0;
    return 2.0 * IR_SPEED_OF_LIGHT * IR_BOLTZMANN_CONSTANT * T_K
           / (lambda_m * lambda_m * lambda_m * lambda_m);
}

double ir_rayleigh_jeans_temperature(double lambda_m, double radiance) {
    if (lambda_m <= 0.0 || radiance < 0.0) return -1.0;
    double l4 = lambda_m * lambda_m * lambda_m * lambda_m;
    return radiance * l4 / (2.0 * IR_SPEED_OF_LIGHT * IR_BOLTZMANN_CONSTANT);
}

/* =========================================================================
 * Wien Approximation: B_Wien = (2*h*c^2/lambda^5) * exp(-h*c/(lambda*kB*T))
 *
 * Short-wavelength limit: h*c/(lambda*kB*T) >> 1
 * Accurate within 1% when lambda*T < 3000 um*K.
 * Useful for remote temperature sensing in the visible/NIR range.
 * ========================================================================= */
double ir_wien_approx_radiance(double lambda_m, double T_K) {
    if (lambda_m <= 0.0 || T_K <= 0.0) return -1.0;
    double c1 = IR_FIRST_RADIATION_C1;
    double c2 = IR_SECOND_RADIATION_C2;
    double x = c2 / (lambda_m * T_K);
    double l5 = lambda_m * lambda_m * lambda_m * lambda_m * lambda_m;
    return (c1 / l5) * exp(-x);
}

/* =========================================================================
 * Emissivity Models
 *
 * Gray body: B_gray = epsilon * B_blackbody
 * Kirchhoff's law: epsilon = alpha = 1 - R - T
 * ========================================================================= */
double ir_graybody_radiance(double lambda_m, double T_K, double emissivity) {
    if (lambda_m <= 0.0 || T_K <= 0.0) return -1.0;
    if (emissivity < 0.0 || emissivity > 1.0) return -1.0;
    double B_bb = ir_planck_spectral_radiance_wavelength(lambda_m, T_K);
    if (B_bb < 0.0) return -1.0;
    return emissivity * B_bb;
}

double ir_kirchhoff_emissivity(double reflectance, double transmittance) {
    if (reflectance < 0.0 || reflectance > 1.0) return -1.0;
    if (transmittance < 0.0 || transmittance > 1.0) return -1.0;
    if (reflectance + transmittance > 1.0) return 0.0;
    double eps = 1.0 - reflectance - transmittance;
    if (eps < 0.0) eps = 0.0;
    if (eps > 1.0) eps = 1.0;
    return eps;
}

/* =========================================================================
 * Simpson's Rule Integration for band radiance
 *
 * Integral_{a}^{b} f(x) dx ~ (h/3) * [f0 + 4*f1 + 2*f2 + 4*f3 + ... + fn]
 * where h = (b-a)/n and n must be even.
 *
 * Complexity: O(n) evaluations of Planck function.
 * Accuracy: O(h^4) for smooth f (4th-order method).
 * ========================================================================= */

/** Adaptive Simpson integration helper */
static double simpson_integrate(double (*f)(double, double),
                                double param, double a, double b, int n) {
    if (n <= 0 || (n % 2) != 0) n += 1;
    double h = (b - a) / n;
    double sum = f(a, param) + f(b, param);
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x, param);
    }
    return (h / 3.0) * sum;
}

/* Wrapper for Planck function with fixed T */
static double planck_wrapper(double lambda, double T) {
    return ir_planck_spectral_radiance_wavelength(lambda, T);
}

/* Wrapper for photon radiance with fixed T */
static double photon_wrapper(double lambda, double T) {
    return ir_planck_photon_radiance(lambda, T);
}

double ir_in_band_radiance(double lambda1_m, double lambda2_m,
                            double T_K, int n_steps) {
    if (lambda1_m <= 0.0 || lambda2_m <= lambda1_m || T_K <= 0.0)
        return -1.0;
    if (n_steps < 10) n_steps = 100;
    return simpson_integrate(planck_wrapper, T_K, lambda1_m, lambda2_m, n_steps);
}

double ir_in_band_exitance(double lambda1_m, double lambda2_m,
                            double T_K, int n_steps) {
    double L = ir_in_band_radiance(lambda1_m, lambda2_m, T_K, n_steps);
    if (L < 0.0) return -1.0;
    return M_PI * L;
}

double ir_photon_flux_band(double lambda1_m, double lambda2_m,
                            double T_K, int n_steps) {
    if (lambda1_m <= 0.0 || lambda2_m <= lambda1_m || T_K <= 0.0)
        return -1.0;
    if (n_steps < 10) n_steps = 100;
    return simpson_integrate(photon_wrapper, T_K, lambda1_m, lambda2_m, n_steps);
}

/* =========================================================================
 * Thermal Contrast & Radiometer Sensitivity
 * ========================================================================= */
double ir_thermal_contrast(double T_target_K, double T_bg_K) {
    if (T_target_K <= 0.0 || T_bg_K <= 0.0) return -1.0;
    double M_target = ir_stefan_boltzmann_exitance(T_target_K);
    double M_bg = ir_stefan_boltzmann_exitance(T_bg_K);
    if (M_bg <= 0.0) return -1.0;
    return (M_target - M_bg) / M_bg;
}

double ir_apparent_temperature_difference(double delta_T, double tau_atm) {
    if (delta_T < 0.0 || tau_atm < 0.0 || tau_atm > 1.0) return -1.0;
    return delta_T * tau_atm;
}

double ir_radiometer_NE_dT(double T_sys_K, double bandwidth_Hz,
                            double integration_s) {
    if (T_sys_K <= 0.0 || bandwidth_Hz <= 0.0 || integration_s <= 0.0)
        return -1.0;
    return T_sys_K / sqrt(bandwidth_Hz * integration_s);
}

/* =========================================================================
 * BLIP (Background-Limited Infrared Photodetector) detectivity
 *
 * D*_BLIP = lambda/(h*c) * sqrt(eta/(2*Phi_bg))
 *
 * This is the theoretical maximum D* when the detector is limited
 * by background photon shot noise. Achieving BLIP-limited performance
 * is a key goal of cooled photon detector design.
 *
 * Reference: Rogalski (2011) "Infrared Detectors", 2nd Ed.
 * ========================================================================= */
double ir_blip_d_star(double quantum_eff, double wavelength_um,
                       double photon_flux_bg) {
    if (quantum_eff <= 0.0 || quantum_eff > 1.0) return -1.0;
    if (wavelength_um <= 0.0 || photon_flux_bg <= 0.0) return -1.0;

    double lambda_m = wavelength_um * 1e-6;  /* convert um to m */
    double hc = IR_PLANCK_CONSTANT * IR_SPEED_OF_LIGHT;

    return (lambda_m / hc) * sqrt(quantum_eff / (2.0 * photon_flux_bg));
}

double ir_photon_irradiance(double T_K, double lambda1_m, double lambda2_m,
                             double solid_angle_sr, int n_steps) {
    if (T_K <= 0.0 || lambda1_m <= 0.0 || lambda2_m <= lambda1_m
        || solid_angle_sr <= 0.0) return -1.0;
    double Q = ir_photon_flux_band(lambda1_m, lambda2_m, T_K, n_steps);
    if (Q < 0.0) return -1.0;
    return Q * solid_angle_sr;
}

/* =========================================================================
 * Brightness Temperature from spectral radiance
 *
 * T_B = h*c / (lambda*kB * ln(1 + 2*h*c^2/(lambda^5 * B_lambda)))
 *
 * Inverts Planck's law: find temperature of a blackbody that would
 * produce the measured spectral radiance at a given wavelength.
 *
 * Used in IR thermometry and remote temperature sensing.
 * ========================================================================= */
double ir_brightness_temperature(double lambda_m, double radiance) {
    if (lambda_m <= 0.0 || radiance <= 0.0) return -1.0;

    double c1 = IR_FIRST_RADIATION_C1;
    double c2 = IR_SECOND_RADIATION_C2;
    double l5 = lambda_m * lambda_m * lambda_m * lambda_m * lambda_m;

    double arg = 1.0 + c1 / (l5 * radiance);
    if (arg <= 1.0) return -1.0;

    return c2 / (lambda_m * log(arg));
}

/* =========================================================================
 * Blackbody Fractional Function F(lambda*T)
 *
 * F = integral_0^lambda M_lambda' d_lambda' / (sigma * T^4)
 *
 * Uses Chang & Rhee (1984) infinite series approximation.
 *
 * f(n) = 15/pi^4 * sum_{m=1}^{inf} [ e^{-m*x} * (x^3/m + 3x^2/m^2 + 6x/m^3 + 6/m^4) ]
 *
 * Reference: Chang, S.L. & Rhee, K.T. (1984) JOSA
 * ========================================================================= */
double ir_blackbody_fraction(double lambda_m, double T_K) {
    if (lambda_m <= 0.0 || T_K <= 0.0) return -1.0;

    double x = IR_SECOND_RADIATION_C2 / (lambda_m * T_K);
    double sum = 0.0;

    /* Converge series: 15 terms gives ~1e-6 accuracy for x > 0.5 */
    for (int m = 1; m <= 20; m++) {
        double mx = m * x;
        double term = exp(-mx) * (x*x*x/m + 3.0*x*x/(m*m)
                                  + 6.0*x/(m*m*m) + 6.0/(m*m*m*m));
        sum += term;
        if (term < 1e-15) break;
    }

    return (15.0 / (M_PI * M_PI * M_PI * M_PI)) * sum;
}

/* =========================================================================
 * Dual-Band Exitance Ratio for Ratio Thermometry
 *
 * R(T) = M_band1(T) / M_band2(T)
 *
 * Ratio thermometry provides emissivity-independent temperature
 * measurement for gray bodies. If emissivity is the same (or has
 * known ratio) in both bands, the ratio M1/M2 depends only on T.
 *
 * Applications: industrial furnaces, gas turbine blades, molten metals.
 * ========================================================================= */
double ir_dual_band_ratio(double T_K,
                           double l1_1, double l1_2,
                           double l2_1, double l2_2,
                           int n_steps) {
    if (T_K <= 0.0) return -1.0;
    if (l1_1 <= 0.0 || l1_2 <= l1_1) return -1.0;
    if (l2_1 <= 0.0 || l2_2 <= l2_1) return -1.0;

    double M1 = ir_in_band_exitance(l1_1, l1_2, T_K, n_steps);
    double M2 = ir_in_band_exitance(l2_1, l2_2, T_K, n_steps);

    if (M1 < 0.0 || M2 <= 0.0) return -1.0;
    return M1 / M2;
}
