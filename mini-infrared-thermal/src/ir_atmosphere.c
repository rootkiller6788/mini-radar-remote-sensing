/**
 * @file    ir_atmosphere.c
 * @brief   Atmospheric Transmission Modeling for IR Systems
 *
 * Implements simplified atmospheric transmission models based on molecular
 * absorption (H2O, CO2, O3, CH4), aerosol scattering (Mie, Rayleigh),
 * and path radiance calculations. Models are comparable to LOWTRAN-7.
 *
 * References:
 *   Berk et al. (2014) "MODTRAN 6: A major upgrade..."
 *   Gordon et al. (2017) "The HITRAN2016 molecular database"
 *   Smith, F.G. (1993) "Atmospheric Propagation of Radiation", SPIE
 */

#include "ir_atmosphere.h"
#include "ir_core.h"
#include <stdio.h>
#include <math.h>

/* Molecular constants */
#define IR_H2O_MOLAR_MASS   18.01528   /* g/mol */
#define IR_CO2_DEFAULT_PPM  420.0      /* 2024 level */
#define IR_CH4_DEFAULT_PPM  1.9        /* atmospheric CH4 */
#define IR_DRY_AIR_MOLAR_MASS 28.9647  /* g/mol */

/* ========================================================================
 * Beer-Lambert Transmission
 *
 * tau = exp(-sigma_ext * R)
 * sigma_ext = sigma_abs + sigma_scat
 * ======================================================================== */

int ir_atmosphere_transmission(const ir_atmosphere_params_t *params,
                                double wavelength_um,
                                ir_atmospheric_transmission_t *result) {
    if (!params || !result || wavelength_um <= 0.0) return -1;

    double sigma_abs = 0.0;

    /* H2O continuum (rotational + vibrational) */
    sigma_abs += ir_h2o_continuum_absorption(wavelength_um,
                                              params->temperature_K,
                                              params->relative_humidity,
                                              params->pressure_hPa);

    /* CO2 absorption */
    if (params->co2_ppm > 0.0)
        sigma_abs += ir_co2_absorption(wavelength_um, params->co2_ppm,
                                        params->temperature_K, params->pressure_hPa);

    /* O3 absorption */
    sigma_abs += ir_o3_absorption(wavelength_um, params->temperature_K);

    /* CH4 absorption */
    sigma_abs += ir_ch4_absorption(wavelength_um, IR_CH4_DEFAULT_PPM);

    /* Aerosol scattering */
    double sigma_scat = ir_aerosol_scattering(wavelength_um,
                                               params->visibility_km,
                                               IR_AEROSOL_RURAL);

    /* Rayleigh scattering */
    sigma_scat += ir_rayleigh_scattering_coefficient(wavelength_um,
                                                      params->pressure_hPa,
                                                      params->temperature_K);

    double sigma_ext = sigma_abs + sigma_scat;
    double R_km = params->path_length_m / 1000.0;

    result->extinction_coeff_per_km = sigma_ext;
    result->absorption_coeff_per_km = sigma_abs;
    result->scattering_coeff_per_km = sigma_scat;
    result->optical_depth = sigma_ext * R_km;
    result->transmission = exp(-result->optical_depth);

    /* Path radiance from atmospheric self-emission (simplified) */
    double L_atm = ir_stefan_boltzmann_exitance(params->temperature_K) / M_PI;
    result->path_radiance = (1.0 - result->transmission) * L_atm;
    result->sky_temperature_K = params->temperature_K
                                 * (1.0 - result->transmission);

    return 0;
}

/* ========================================================================
 * H2O Continuum Absorption
 *
 * Water vapor is the dominant absorber in the 8-14 um LWIR window.
 * Uses Roberts et al. (1976) semi-empirical model.
 * ======================================================================== */

double ir_h2o_continuum_absorption(double wavelength_um,
                                    double temperature_K,
                                    double relative_humidity,
                                    double pressure_hPa) {
    if (wavelength_um <= 0.0 || temperature_K <= 0.0) return -1.0;

    /* Saturation vapor pressure (Tetens formula) */
    double es_hPa = 6.112 * exp(17.67 * (temperature_K - 273.15)
                                 / (temperature_K - 29.65));
    double e_hPa = es_hPa * relative_humidity / 100.0;

    /* H2O self-broadening continuum coefficient [1/(km*hPa^2)] */
    double Cs = 1.0e-20;  /* approximate for 8-14 um window */

    /* Foreign-broadening */
    double Cf = 1.0e-22;

    double k_self = Cs * e_hPa * e_hPa;
    double k_foreign = Cf * e_hPa * (pressure_hPa - e_hPa);
    double k_cont = k_self + k_foreign;

    /* H2O line absorption: approximate using 6.3 um band model */
    double k_line = 0.0;
    if (wavelength_um > 5.0 && wavelength_um < 8.0)
        k_line = 0.5 * exp(-((wavelength_um - 6.3) * (wavelength_um - 6.3)) / 0.5);

    return k_cont + k_line;
}

/* ========================================================================
 * CO2 Absorption
 *
 * Major bands: 2.7 um, 4.3 um (strongest), 15 um
 * The 4.3 um band is essentially opaque; 15 um band is important for
 * climate (greenhouse effect).
 * ======================================================================== */

double ir_co2_absorption(double wavelength_um, double co2_ppm,
                          double temperature_K, double pressure_hPa) {
    if (wavelength_um <= 0.0 || temperature_K <= 0.0) return -1.0;

    (void)pressure_hPa;
    double k = 0.0;
    double co2_norm = co2_ppm / IR_CO2_DEFAULT_PPM;

    /* 4.3 um band - very strong */
    if (wavelength_um > 4.15 && wavelength_um < 4.45)
        k = co2_norm * 100.0 * exp(-(wavelength_um - 4.3)
                                    * (wavelength_um - 4.3) / 0.005);

    /* 15 um band - strong */
    if (wavelength_um > 13.0 && wavelength_um < 17.0)
        k = co2_norm * 10.0 * exp(-(wavelength_um - 15.0)
                                   * (wavelength_um - 15.0) / 2.0);

    /* 2.7 um band - moderate */
    if (wavelength_um > 2.6 && wavelength_um < 2.8)
        k = co2_norm * 1.0;

    return k;
}

/* ========================================================================
 * Ozone (O3) and Methane (CH4) Absorption
 * ======================================================================== */

double ir_o3_absorption(double wavelength_um, double temperature_K) {
    if (wavelength_um <= 0.0 || temperature_K <= 0.0) return -1.0;

    /* O3 has a strong band centered at 9.6 um */
    double k = 0.0;
    if (wavelength_um > 9.0 && wavelength_um < 10.2)
        k = 0.05 * exp(-(wavelength_um - 9.6) * (wavelength_um - 9.6) / 0.05);
    return k;
}

double ir_ch4_absorption(double wavelength_um, double ch4_ppm) {
    if (wavelength_um <= 0.0 || ch4_ppm < 0.0) return -1.0;
    double k = 0.0;
    /* CH4 band near 7.7 um */
    if (wavelength_um > 7.0 && wavelength_um < 8.4)
        k = 0.001 * ch4_ppm * exp(-(wavelength_um - 7.7)
                                   * (wavelength_um - 7.7) / 0.3);
    return k;
}

/* ========================================================================
 * Aerosol Scattering (Mie)
 * ======================================================================== */

double ir_aerosol_scattering(double wavelength_um, double visibility_km,
                              ir_aerosol_type_t aerosol_type) {
    if (wavelength_um <= 0.0 || visibility_km <= 0.0) return -1.0;

    /* Koschmieder: sigma_scat(0.55um) = 3.912 / V [1/km] */
    double sigma_550 = 3.912 / visibility_km;

    /* Angstrom exponent p depends on aerosol size distribution */
    double p;
    switch (aerosol_type) {
        case IR_AEROSOL_RURAL:    p = 1.3; break;
        case IR_AEROSOL_URBAN:    p = 1.0; break;
        case IR_AEROSOL_MARITIME: p = 0.5; break;
        case IR_AEROSOL_DESERT:   p = 0.2; break;
        default:                  p = 1.3; break;
    }

    /* Scaled to IR wavelength: sigma(lambda) = sigma(0.55) * (0.55/lambda)^p */
    return sigma_550 * pow(0.55 / wavelength_um, p);
}

double ir_rayleigh_scattering_coefficient(double wavelength_um,
                                           double pressure_hPa,
                                           double temperature_K) {
    if (wavelength_um <= 0.0 || pressure_hPa <= 0.0 || temperature_K <= 0.0)
        return -1.0;

    /* Rayleigh scattering cross section per molecule ~ 1/lambda^4 */
    double sigma_rayleigh_m2 = 5.45e-32 * pow(wavelength_um, -4.0);

    /* Number density from ideal gas law: n = P/(kB*T) */
    double P_Pa = pressure_hPa * 100.0;
    double n = P_Pa / (IR_BOLTZMANN_CONSTANT * temperature_K);  /* 1/m3 */

    return sigma_rayleigh_m2 * n * 1000.0;  /* convert m^-1 to km^-1 */
}

/* ========================================================================
 * Atmospheric Windows
 * ======================================================================== */

int ir_is_atmospheric_window(double wavelength_um) {
    if (wavelength_um <= 0.0) return 0;

    /* Major windows: 0.95-1.05 (NIR), 1.5-1.8 (SWIR), 2.0-2.4 (SWIR),
       3.5-4.1 (MWIR), 4.5-5.0 (MWIR), 8.0-9.5 (LWIR), 10.0-13.0 (LWIR) */
    if (wavelength_um >= 3.5 && wavelength_um <= 4.1) return 1;
    if (wavelength_um >= 4.5 && wavelength_um <= 5.0) return 1;
    if (wavelength_um >= 8.0 && wavelength_um <= 14.0) return 1;
    if (wavelength_um >= 0.95 && wavelength_um <= 1.05) return 1;
    if (wavelength_um >= 1.5 && wavelength_um <= 1.8) return 1;
    if (wavelength_um >= 2.0 && wavelength_um <= 2.4) return 1;

    return 0;
}

double ir_band_average_transmission(const ir_atmosphere_params_t *params,
                                     double lambda1_um, double lambda2_um,
                                     int n_samples) {
    if (!params || lambda1_um <= 0.0 || lambda2_um <= lambda1_um || n_samples < 2)
        return -1.0;

    double d_lambda = (lambda2_um - lambda1_um) / (n_samples - 1);
    double sum_tau = 0.0;

    for (int i = 0; i < n_samples; i++) {
        double lambda = lambda1_um + i * d_lambda;
        ir_atmospheric_transmission_t result;
        if (ir_atmosphere_transmission(params, lambda, &result) == 0)
            sum_tau += result.transmission;
    }

    return sum_tau / n_samples;
}

/* ========================================================================
 * Weather Attenuation: Rain and Fog
 * ======================================================================== */

double ir_rain_attenuation(double frequency_GHz, double rain_rate_mm_per_hr) {
    if (frequency_GHz <= 0.0 || rain_rate_mm_per_hr < 0.0) return -1.0;

    /* ITU-R P.838: gamma_r = k * R^alpha [dB/km] */
    double k, alpha;
    if (frequency_GHz < 100.0) {
        /* Microwave band */
        k = 0.01;
        alpha = 1.2;
    } else {
        /* IR band (convert to equivalent attenuation) */
        k = 0.1;
        alpha = 0.8;
    }

    return k * pow(rain_rate_mm_per_hr, alpha);
}

double ir_fog_attenuation(double wavelength_um, double visibility_km) {
    if (wavelength_um <= 0.0 || visibility_km <= 0.0) return -1.0;

    /* Fog attenuation: sigma ~ 3.912 / V * (0.55/lambda)^q */
    /* q ~ 0 (wavelength-independent) for dense fog (large droplets) */
    /* q ~ 1.6 for haze (small droplets) */
    double q = (visibility_km < 0.5) ? 0.0 : 1.6;
    double sigma_550 = 3.912 / visibility_km;
    return sigma_550 * pow(0.55 / wavelength_um, q);
}
