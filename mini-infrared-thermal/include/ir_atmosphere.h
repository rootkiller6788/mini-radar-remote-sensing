/**
 * @file    ir_atmosphere.h
 * @brief   Atmospheric Transmission Modeling for Infrared Systems
 *
 * Knowledge: L1 (atmospheric windows, extinction), L2 (absorption, scattering),
 *            L4 (Beer-Lambert law), L7 (applications: weather, climate)
 *
 * The atmosphere significantly affects IR system performance through:
 * 1. Molecular absorption (H2O, CO2, O3, CH4, N2O, CO)
 * 2. Aerosol scattering (Mie, Rayleigh)
 * 3. Turbulence-induced scintillation
 *
 * Primary atmospheric windows: 3-5 um (MWIR) and 8-14 um (LWIR)
 *
 * References:
 *   MODTRAN (Berk et al., 2014) - MODerate resolution atmospheric TRANsmission
 *   HITRAN database (Gordon et al., 2017)
 *   Smith, F.G. (1993) "Atmospheric Propagation of Radiation", SPIE
 */

#ifndef IR_ATMOSPHERE_H
#define IR_ATMOSPHERE_H

#include <stddef.h>
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L1: Atmospheric Model Parameters
 * ======================================================================== */

typedef enum {
    IR_ATMOSPHERE_STANDARD = 0,   /**< 1976 US Standard Atmosphere */
    IR_ATMOSPHERE_TROPICAL,       /**< Tropical (15 deg N) */
    IR_ATMOSPHERE_SUB_ARCTIC_WINTER, /**< Sub-Arctic Winter */
    IR_ATMOSPHERE_SUB_ARCTIC_SUMMER, /**< Sub-Arctic Summer */
    IR_ATMOSPHERE_MID_LATITUDE_WINTER,
    IR_ATMOSPHERE_MID_LATITUDE_SUMMER,
    IR_ATMOSPHERE_DESERT,
    IR_ATMOSPHERE_COUNT
} ir_atmosphere_model_t;

typedef struct {
    ir_atmosphere_model_t model;
    double temperature_K;       /**< surface air temperature [K] */
    double pressure_hPa;        /**< surface pressure [hPa] */
    double relative_humidity;   /**< RH [0-100] */
    double visibility_km;       /**< meteorological visibility [km] */
    double altitude_m;          /**< sensor altitude [m] */
    double zenith_angle_deg;    /**< path zenith angle [deg] */
    double co2_ppm;             /**< CO2 concentration [ppm] (default 420) */
    double path_length_m;       /**< slant path length [m] */
} ir_atmosphere_params_t;

/* ========================================================================
 * L4: Beer-Lambert Law for Atmospheric Transmission
 *
 * tau(lambda) = exp(-sigma_ext(lambda) * R)
 *
 * where sigma_ext = sigma_absorption + sigma_scattering
 * ======================================================================== */

typedef struct {
    double transmission;        /**< tau [0,1] */
    double optical_depth;       /**< tau_od = -ln(transmission) */
    double extinction_coeff_per_km; /**< sigma_ext [1/km] */
    double absorption_coeff_per_km; /**< sigma_abs [1/km] */
    double scattering_coeff_per_km; /**< sigma_scat [1/km] */
    double path_radiance;       /**< L_path [W/(sr*m2)] from atmospheric emission */
    double sky_temperature_K;   /**< apparent sky temperature [K] */
} ir_atmospheric_transmission_t;

/* ========================================================================
 * L2: Molecular Absorption Models
 * ======================================================================== */

/** Compute atmospheric transmission for a given wavelength and path.
 *  Simplified band-model approach (comparable to LOWTRAN). */
int ir_atmosphere_transmission(const ir_atmosphere_params_t *params,
                                double wavelength_um,
                                ir_atmospheric_transmission_t *result);

/** Water vapor continuum absorption [1/km].
 *  Dominant in the 8-14 um window. Uses Roberts et al. (1976) model. */
double ir_h2o_continuum_absorption(double wavelength_um,
                                    double temperature_K,
                                    double relative_humidity,
                                    double pressure_hPa);

/** CO2 absorption coefficient at 4.3 um and 15 um bands [1/km] */
double ir_co2_absorption(double wavelength_um, double co2_ppm,
                          double temperature_K, double pressure_hPa);

/** Ozone (O3) absorption in 9.6 um band [1/km] */
double ir_o3_absorption(double wavelength_um, double temperature_K);

/** Methane (CH4) absorption near 7.7 um [1/km] */
double ir_ch4_absorption(double wavelength_um, double ch4_ppm);

/* ========================================================================
 * L2: Aerosol Scattering Models
 * ======================================================================== */

typedef enum {
    IR_AEROSOL_RURAL = 0,
    IR_AEROSOL_URBAN,
    IR_AEROSOL_MARITIME,
    IR_AEROSOL_DESERT,
    IR_AEROSOL_COUNT
} ir_aerosol_type_t;

/** Mie scattering coefficient for aerosols [1/km].
 *  Uses Koschmieder formula: sigma_scat = 3.912/V for visible,
 *  scaled to IR by lambda^-p where p depends on aerosol type. */
double ir_aerosol_scattering(double wavelength_um, double visibility_km,
                              ir_aerosol_type_t aerosol_type);

/** Rayleigh scattering coefficient [1/km].
 *  sigma_Rayleigh = const / lambda^4.
 *  Negligible in thermal IR compared to aerosol scattering. */
double ir_rayleigh_scattering_coefficient(double wavelength_um,
                                           double pressure_hPa,
                                           double temperature_K);

/* ========================================================================
 * L2: Atmospheric Windows
 * ======================================================================== */

/** Check if wavelength is within a major atmospheric window.
 *  Returns 1 if in window (high transmission), 0 otherwise. */
int ir_is_atmospheric_window(double wavelength_um);

/** Integrated transmission over a spectral band.
 *  tau_band = integral tau(lambda) d_lambda / (lambda2 - lambda1) */
double ir_band_average_transmission(const ir_atmosphere_params_t *params,
                                     double lambda1_um, double lambda2_um,
                                     int n_samples);

/** Path radiance from atmospheric emission [W/(sr*m2)].
 *  L_path = integral (1 - tau(lambda)) * B_lambda(lambda, T_atm) d_lambda */
double ir_path_radiance(const ir_atmosphere_params_t *params,
                         double lambda1_um, double lambda2_um,
                         int n_samples);

/** Effective sky temperature for downwelling radiance [K] */
double ir_sky_temperature(const ir_atmosphere_params_t *params,
                           double lambda1_um, double lambda2_um,
                           int n_samples);

/** Rain attenuation at mm-wave and IR [dB/km].
 *  Uses ITU-R P.838 model. */
double ir_rain_attenuation(double frequency_GHz, double rain_rate_mm_per_hr);

/** Fog attenuation model [1/km].
 *  Uses empirical relation: sigma_fog = A / V where V is visibility [km]. */
double ir_fog_attenuation(double wavelength_um, double visibility_km);

#endif /* IR_ATMOSPHERE_H */
