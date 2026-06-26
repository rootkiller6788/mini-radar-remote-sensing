/**
 * @file    hyperspectral_radiometry.h
 * @brief   Radiometric Calibration and Atmospheric Correction — L2 + L4 + L6
 *
 * @details Implements radiometric calibration converting raw DN to radiance,
 *          top-of-atmosphere reflectance computation, atmospheric correction
 *          (empirical line, dark object subtraction, simplified radiative
 *          transfer), and BRDF modeling.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Digital Number (DN), gain/offset calibration
 *     - Top-of-atmosphere (TOA) radiance and reflectance
 *     - Surface reflectance (atmospherically corrected)
 *     - BRDF (Bidirectional Reflectance Distribution Function)
 *     - Solar irradiance spectrum (Thuillier, Kurucz, ASTM E490)
 *   L2 - Core Concepts:
 *     - Radiometric calibration chain: DN → radiance → reflectance
 *     - Atmospheric path radiance, adjacency effect
 *     - Dark object subtraction (DOS)
 *     - Empirical line calibration
 *   L3 - Mathematical Structures:
 *     - Linear calibration models
 *     - Radiative transfer equation (simplified)
 *   L4 - Fundamental Laws:
 *     - Inverse square law for solar irradiance
 *     - Cosine law for irradiance (Lambert's cosine law)
 *   L6 - Canonical Problems:
 *     - Atmospheric correction of AVIRIS data
 *     - Water vapor retrieval from 940nm / 1140nm bands
 *     - Aerosol optical depth estimation
 *
 * Reference:
 *   - Schott, "Remote Sensing: The Image Chain Approach" (2007)
 *   - Chavez, "Image-Based Atmospheric Corrections" (1996), PE&RS
 *   - Vermote et al., "Second Simulation of the Satellite Signal" (6S model) (1997)
 */

#ifndef HYPERSPECTRAL_RADIOMETRY_H
#define HYPERSPECTRAL_RADIOMETRY_H

#include "hyperspectral_core.h"
#include "hyperspectral_spectroscopy.h"

/* ─── L1: Radiometric Calibration Coefficients ──────────────────── */

/**
 * @brief Per-band radiometric calibration parameters.
 */
typedef struct {
    size_t   band_index;
    double   gain;                /**< Gain [(W·m⁻²·sr⁻¹·μm⁻¹)/DN] */
    double   offset;              /**< Offset [W·m⁻²·sr⁻¹·μm⁻¹] */
    double   calibration_uncertainty; /**< 1-σ uncertainty in calibration */
    double   nonlinearity_coeff;  /**< Nonlinearity correction coefficient */
    char     cal_source[64];      /**< "laboratory", "on-board", "vicarious" */
} hspec_cal_coefficient_t;

/* ─── L1: Solar Irradiance Model ────────────────────────────────── */

/**
 * @brief Solar spectral irradiance at top of atmosphere.
 *
 * Typically from Thuillier (2004), Kurucz, or ASTM E490-00a.
 * Key spectral lines: Fraunhofer lines, H-α (656.3nm), Ca II.
 */
typedef struct {
    size_t   nbands;
    double  *wavelengths;
    double  *solar_irradiance;    /**< E_sun(λ) [W·m⁻²·μm⁻¹] */
    double   total_solar_irradiance; /**< Integrated TSI ≈ 1361 W/m² */
    char     model_name[64];      /**< "Thuillier2004", "Kurucz", "ASTM_E490" */
} hspec_solar_irradiance_t;

/* ─── L1: Atmospheric Correction Parameters ──────────────────────── */

/**
 * @brief Parameters for atmospheric correction.
 */
typedef struct {
    double   visibility_km;       /**< Meteorological visibility [km] */
    double   aerosol_od_550;      /**< Aerosol optical depth at 550 nm */
    double   water_vapor_cm;      /**< Total column water vapor [g·cm⁻²] */
    double   ozone_cm_atm;        /**< Total column ozone [cm-atm] */
    double   co2_ppm;             /**< CO₂ concentration [ppm] */
    double   surface_pressure_mb; /**< Surface pressure [mb] */
    double   sensor_altitude_km;  /**< Sensor altitude [km] */
    double   target_altitude_km;  /**< Target surface altitude [km] */
    double   solar_zenith_deg;    /**< Solar zenith angle [degrees] */
    double   view_zenith_deg;     /**< View zenith angle [degrees] */
    double   relative_azimuth_deg;/**< Relative azimuth [degrees] */
    char     aerosol_model[32];   /**< "continental", "maritime", "urban", "desert" */
} hspec_atm_params_t;

/**
 * @brief Results of atmospheric correction.
 */
typedef struct {
    size_t   nbands;
    double  *path_radiance;       /**< Atmospheric path radiance L_p(λ) */
    double  *transmittance_down;  /**< Downwelling transmittance T↓(λ) */
    double  *transmittance_up;    /**< Upwelling transmittance T↑(λ) */
    double  *spherical_albedo;    /**< Atmospheric spherical albedo S(λ) */
    double  *direct_irradiance;   /**< Direct solar irradiance at surface */
    double  *diffuse_irradiance;  /**< Diffuse skylight irradiance at surface */
    double  *surface_reflectance; /**< Retrieved ρ_surface(λ) */
    double  *water_vapor_map;     /**< Per-pixel water vapor estimate */
    double  *aod_map;             /**< Per-pixel aerosol optical depth */
} hspec_atm_correction_result_t;

/* ─── L2: Dark Object Subtraction Parameters ───────────────────── */

typedef struct {
    double   haze_removal_start;  /**< Starting wavelength for haze model */
    double   haze_removal_end;    /**< Ending wavelength for haze model */
    double   dark_object_percentile; /**< Percentile for dark object (default 1%) */
    int      use_rayleigh_lut;    /**< Use Rayleigh LUT vs simplified */
} hspec_dos_params_t;

/* ─── L1: BRDF Model ────────────────────────────────────────────── */

typedef enum {
    HSPEC_BRDF_LAMBERTIAN    = 0,  /**< ρ = constant (ideal diffuse) */
    HSPEC_BRDF_MINNAERT      = 1,  /**< ρ' = ρ·(cosθ_i)^{k-1}·(cosθ_v)^{k-1} */
    HSPEC_BRDF_ROUJEAN       = 2,  /**< Three-parameter kernel-driven */
    HSPEC_BRDF_ROSS_LI       = 3,  /**< Ross-Thick + Li-Sparse kernel */
    HSPEC_BRDF_RAHMAN        = 4   /**< Rahman-Pinty-Verstraete model */
} hspec_brdf_model_type_t;

/**
 * @brief BRDF model parameters.
 */
typedef struct {
    hspec_brdf_model_type_t model_type;
    /** Minnaert: k exponent */
    double   minnaert_k;
    /** Ross-Li: f_iso, f_vol, f_geo weights */
    double   ross_li_f_iso;
    double   ross_li_f_vol;
    double   ross_li_f_geo;
    /** Rahman-Pinty-Verstraete */
    double   rpv_rho0;           /**< Overall amplitude */
    double   rpv_k;              /**< Minnaert-like parameter */
    double   rpv_theta;          /**< Asymmetry parameter (Henvey-Greenstein) */
    double   hot_spot_width;     /**< Hot spot width parameter */
} hspec_brdf_params_t;

/* ─── API: Radiometric Calibration ───────────────────────────────── */

/**
 * Convert raw Digital Numbers to at-sensor radiance.
 *
 * @math L_sensor(λ_i) = gain_i · DN_i + offset_i
 *
 * @param dn             Raw digital number array [nbands]
 * @param cal_coeffs     Calibration coefficients [nbands]
 * @param nbands         Number of bands
 * @param radiance_out   Output radiance [W·m⁻²·sr⁻¹·μm⁻¹] [nbands]
 * @return               0 on success
 *
 * Complexity: O(nbands)
 */
int hspec_dn_to_radiance(const uint16_t *dn, const hspec_cal_coefficient_t *cal_coeffs,
                          size_t nbands, double *radiance_out);

/**
 * Compute top-of-atmosphere reflectance from radiance.
 *
 * @math ρ_TOA(λ) = π·L_sensor(λ)·d² / (E_sun(λ)·cos(θ_s))
 *       where d = Earth-Sun distance [AU]
 *
 * @param radiance      At-sensor radiance [nbands]
 * @param solar_irrad   Solar irradiance E_sun(λ) [nbands]
 * @param nbands        Number of bands
 * @param earth_sun_au  Earth-Sun distance [AU]
 * @param solar_zenith  Solar zenith angle [radians]
 * @param toa_refl_out  Output TOA reflectance [nbands]
 * @return              0 on success
 *
 * Reference: Chander et al. (2009), Remote Sensing of Environment
 */
int hspec_radiance_to_toa_refl(const double *radiance, const double *solar_irrad,
                                size_t nbands, double earth_sun_au,
                                double solar_zenith, double *toa_refl_out);

/* ─── API: Atmospheric Correction ────────────────────────────────── */

/**
 * Dark Object Subtraction (DOS) — simplest atmospheric correction.
 *
 * Assumes dark objects (shadows, clear water) have near-zero reflectance
 * in certain bands. Subtracts the minimum (or low-percentile) DN as an
 * estimate of atmospheric path radiance.
 *
 * @param dc           Datacube (in DN or radiance)
 * @param params       DOS parameters
 * @param corrected_out Output surface-reflectance-like values [npixels][nbands]
 * @return             0 on success
 *
 * Complexity: O(n_pixels·nbands)
 * Reference: Chavez (1988, 1996)
 */
int hspec_dos_correction(const hspec_datacube_t *dc,
                          const hspec_dos_params_t *params,
                          double *corrected_out);

/**
 * Empirical Line calibration — uses known reflectance targets in scene.
 *
 * Requires at least 2 targets with known reflectance at each band.
 *
 * @math ρ = a(λ)·L(λ) + b(λ)
 *       where a,b are solved per-band from known targets.
 *
 * @param dc              Datacube
 * @param known_reflectance Known reflectance values [n_targets][nbands]
 * @param target_pixels   Pixel indices of known targets [n_targets]
 * @param n_targets       Number of calibration targets (≥ 2)
 * @param corrected_out   Output surface reflectance
 * @return                0 on success
 *
 * Complexity: O(n_pixels·nbands)
 */
int hspec_empirical_line_correction(const hspec_datacube_t *dc,
                                     const double *known_reflectance,
                                     const size_t *target_pixels,
                                     size_t n_targets,
                                     double *corrected_out);

/**
 * Simplified radiative transfer atmospheric correction.
 *
 * Uses a parameterized atmospheric model (6S-based simplifications)
 * to estimate path radiance and transmittance, then invert to
 * obtain surface reflectance.
 *
 * @math ρ_surface = (ρ_TOA/T_g - ρ_path) / (T↑·T↓ + S·(ρ_TOA/T_g - ρ_path))
 *
 * @param dc         Datacube (should contain TOA radiance)
 * @param atm_params Atmospheric parameters
 * @param solar_irr  Solar irradiance spectrum
 * @return           Atmospheric correction result
 *
 * Complexity: O(n_pixels·nbands)
 * Reference: Vermote et al. (1997), IEEE TGRS
 */
hspec_atm_correction_result_t hspec_atmospheric_correction(
    const hspec_datacube_t *dc, const hspec_atm_params_t *atm_params,
    const hspec_solar_irradiance_t *solar_irr);

void hspec_atm_correction_result_free(hspec_atm_correction_result_t *ac);

/* ─── API: Solar Irradiance ──────────────────────────────────────── */

/**
 * Compute standard solar irradiance (Thuillier 2004) at a given wavelength
 * using interpolation of tabulated values.
 *
 * @param lambda  Wavelength [μm]
 * @return        Solar spectral irradiance [W·m⁻²·μm⁻¹]
 */
double hspec_solar_irradiance_at_wavelength(double lambda);

/**
 * Initialize a solar irradiance spectrum for given bands.
 *
 * @param bands        Wavelength array [nbands] (μm)
 * @param nbands       Number of bands
 * @param sirr         Output solar irradiance structure (pre-allocated)
 * @return             0 on success
 */
int hspec_solar_irradiance_init(const double *bands, size_t nbands,
                                 hspec_solar_irradiance_t *sirr);
void hspec_solar_irradiance_free(hspec_solar_irradiance_t *sirr);

/* ─── API: BRDF ──────────────────────────────────────────────────── */

/**
 * Evaluate Minnaert BRDF-corrected reflectance.
 *
 * @math ρ_corrected = ρ_observed · (cos θ_i)^(1-k) · (cos θ_v)^(1-k)
 *
 * @param reflectance_obs  Observed reflectance per band [nbands]
 * @param nbands           Number of bands
 * @param solar_zenith     Solar zenith angle [radians]
 * @param view_zenith      View zenith angle [radians]
 * @param minnaert_k       Minnaert k exponent per band [nbands]
 * @param corrected_out    Output corrected reflectance [nbands]
 * @return                 0 on success
 *
 * Complexity: O(nbands)
 * Reference: Minnaert (1941), Tellus
 */
int hspec_brdf_minnaert_correct(const double *reflectance_obs, size_t nbands,
                                 double solar_zenith, double view_zenith,
                                 const double *minnaert_k, double *corrected_out);

/**
 * Evaluate Ross-Li kernel-driven BRDF model.
 *
 * @math ρ(θ_i,θ_v,φ) = f_iso + f_vol·K_vol(θ_i,θ_v,φ) + f_geo·K_geo(θ_i,θ_v,φ)
 *
 * @param f_iso         Isotropic weight
 * @param f_vol         Volumetric (Ross-Thick) weight
 * @param f_geo         Geometric-optical (Li-Sparse) weight
 * @param solar_zenith  θ_i [radians]
 * @param view_zenith   θ_v [radians]
 * @param rel_azimuth   φ [radians]
 * @return              BRDF reflectance value
 *
 * Reference: Roujean et al. (1992), JGR;
 *            Wanner, Li, Strahler (1995), JGR
 */
double hspec_brdf_ross_li(double f_iso, double f_vol, double f_geo,
                           double solar_zenith, double view_zenith,
                           double rel_azimuth);

/* ─── API: Water Vapor Estimation ────────────────────────────────── */

/**
 * Estimate column water vapor from the 940nm water absorption feature.
 *
 * Uses the Continuum Interpolated Band Ratio (CIBR) method.
 *
 * @math CIBR = R(measurement) / (a·R(λ_left) + b·R(λ_right))
 *
 * @param reflectance  Reflectance spectrum [nbands]
 * @param wavelengths  Wavelengths [nbands]
 * @param nbands       Number of bands
 * @return             Estimated water vapor column [g·cm⁻²]
 *
 * Reference: Gao & Goetz (1990), JGR
 */
double hspec_water_vapor_estimate(const double *reflectance,
                                   const double *wavelengths, size_t nbands);

#endif /* HYPERSPECTRAL_RADIOMETRY_H */