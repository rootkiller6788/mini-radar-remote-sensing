/**
 * @file    hyperspectral_spectroscopy.h
 * @brief   Spectroscopy Principles — L1 Definitions + L4 Fundamental Laws
 *
 * @details Covers the physical laws governing electromagnetic radiation,
 *          absorption, emission, and scattering relevant to hyperspectral
 *          remote sensing. Provides the theoretical foundation for spectral
 *          signature interpretation.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Radiance L [W·m⁻²·sr⁻¹·μm⁻¹], Irradiance E [W·m⁻²·μm⁻¹]
 *     - Reflectance ρ, Transmittance τ, Absorptance α
 *     - Emissivity ε, Optical depth
 *     - Absorption features (electronic, vibrational, rotational)
 *   L2 - Core Concepts:
 *     - Kirchhoff's law: ε(λ) = α(λ) at thermal equilibrium
 *     - Atmospheric absorption bands (H₂O, CO₂, O₃, CH₄)
 *     - Continuum removal for absorption feature analysis
 *   L3 - Mathematical Structures:
 *     - Exponential attenuation (Beer-Lambert-Bouguer)
 *     - Spectral integration
 *   L4 - Fundamental Laws:
 *     - Planck's law: spectral radiance of a blackbody
 *     - Wien's displacement law: λ_max·T = b = 2.898×10⁻³ m·K
 *     - Stefan-Boltzmann law: M = σ·T⁴
 *     - Beer-Lambert-Bouguer law: I = I₀·exp(-α·d)
 *     - Kirchhoff's law of thermal radiation
 *
 * Reference:
 *   - Planck, "On the Law of Distribution of Energy in the Normal Spectrum" (1901)
 *   - Goody & Yung, "Atmospheric Radiation: Theoretical Basis" (1989)
 *   - Clark & Roush, "Reflectance Spectroscopy" (1984), JGR
 */

#ifndef HYPERSPECTRAL_SPECTROSCOPY_H
#define HYPERSPECTRAL_SPECTROSCOPY_H

#include "hyperspectral_core.h"

/* ─── L1: Spectral Absorption Feature ────────────────────────────────── */

/**
 * @brief Characterized absorption band in a reflectance spectrum.
 *
 * Absorption features arise from electronic transitions (VIS),
 * vibrational overtones (SWIR), and fundamental vibrations (TIR).
 */
typedef struct {
    double   center_wavelength;    /**< Center of absorption [same unit as wavelengths] */
    double   depth;                /**< Depth relative to continuum (0-1) */
    double   width_fwhm;           /**< Full width at half maximum */
    double   area;                 /**< Integrated area of absorption */
    double   asymmetry;            /**< Asymmetry ratio (left/right half-widths) */
    char     feature_type[64];     /**< E.g. "Fe3+ crystal field", "Al-OH bend" */
    char     mineral_association[64]; /**< Associated mineral class */
} hspec_absorption_feature_t;

/* ─── L1: Continuum-Removed Spectrum ─────────────────────────────────── */

/**
 * @brief Results of continuum removal spectral normalization.
 *
 * Continuum removal (CR) isolates absorption features by dividing
 * the spectrum by a convex hull envelope.
 *
 * @math R_CR(λ) = R(λ) / R_continuum(λ)
 */
typedef struct {
    size_t   nbands;
    double  *original_reflectance;   /**< Input spectrum */
    double  *continuum;              /**< Convex hull continuum [nbands] */
    double  *continuum_removed;      /**< R_CR(λ) = R(λ) / continuum(λ) [nbands] */
    size_t   num_absorption_features;
    hspec_absorption_feature_t *features; /**< Detected absorption features */
} hspec_continuum_removed_t;

/* ─── L4: Planck Blackbody Result ──────────────────────────────────── */

/**
 * @brief Blackbody spectral radiance as computed from Planck's law.
 *
 * @math B_λ(λ, T) = (2·h·c²/λ⁵) · 1/(exp(h·c/(λ·k·T)) - 1) [W·m⁻²·sr⁻¹·μm⁻¹]
 */
typedef struct {
    double   temperature_k;     /**< Temperature [K] */
    size_t   npoints;           /**< Number of wavelength points */
    double  *wavelengths;       /**< Wavelength grid [m] */
    double  *spectral_radiance; /**< B_λ(λ, T) [W·m⁻²·sr⁻¹·m⁻¹] */
    double   peak_wavelength;   /**< Wien peak λ_max [m] */
    double   total_exitance;    /**< Stefan-Boltzmann: M = σ·T⁴ [W·m⁻²] */
} hspec_blackbody_t;

/* ─── L4: Beer-Lambert-Bouguer Transmission ─────────────────────────── */

/**
 * @brief Transmission through an absorbing/scattering medium.
 *
 * @math T(λ) = exp(-τ(λ)) = exp(-Σ_i σ_i(λ)·n_i·d)
 *       where τ(λ) is the optical depth.
 */
typedef struct {
    size_t   nbands;
    double  *incident_radiance;    /**< I₀(λ) entering the medium */
    double  *transmitted_radiance; /**< I(λ) after attenuation */
    double  *optical_depth;        /**< τ(λ) = -ln(T(λ)) */
    double  *absorption_coeff;     /**< α(λ) per unit path length */
    double   path_length;          /**< Geometric path length [m] */
    double   mass_column;          /**< Absorber amount [kg·m⁻²] */
} hspec_beer_lambert_result_t;

/* ─── L1: Spectral Signature Type ───────────────────────────────────── */

typedef enum {
    HSPEC_SIG_FLAT         = 0,  /**< Featureless, spectrally neutral */
    HSPEC_SIG_RED_SLOPE    = 1,  /**< Monotonic increase with wavelength */
    HSPEC_SIG_BLUE_SLOPE   = 2,  /**< Monotonic decrease with wavelength */
    HSPEC_SIG_VEGETATION   = 3,  /**< Red edge + NIR plateau */
    HSPEC_SIG_IRON_OXIDE   = 4,  /**< Fe³⁺ absorption at ~900nm */
    HSPEC_SIG_CLAY         = 5,  /**< Al-OH absorption at ~2200nm */
    HSPEC_SIG_CARBONATE    = 6,  /**< CO₃²⁻ absorption at ~2350nm */
    HSPEC_SIG_WATER        = 7,  /**< Strong absorption features */
    HSPEC_SIG_SNOW_ICE     = 8,  /**< High VIS, decreasing NIR/SWIR */
    HSPEC_SIG_MANMADE      = 9,  /**< Paints, metals, asphalt */
    HSPEC_SIG_COUNT        = 10
} hspec_signature_type_t;

/* ─── API Declarations ───────────────────────────────────────────────── */

/**
 * Planck's law: compute blackbody spectral radiance.
 *
 * @math B_λ(λ,T) = (2hc²/λ⁵) * 1/(exp(hc/(λkT)) - 1)
 *
 * @param T              Temperature [K]
 * @param wavelengths    Wavelength array [m]
 * @param npoints        Number of wavelength points
 * @param radiance_out   Output B_λ(λ,T) [W·m⁻²·sr⁻¹·m⁻¹]
 * @return               0 on success, -1 on invalid input
 *
 * Complexity: O(npoints)
 * Reference: Planck (1901), "On the Law of Distribution of Energy"
 */
int hspec_planck_radiance(double T, const double *wavelengths, size_t npoints,
                           double *radiance_out);

/**
 * Wien's displacement law: find lambda_max for a given temperature.
 *
 * @param T  Temperature [K]
 * @return   lambda_max = b/T [m] (returns -1.0 for T <= 0)
 *
 * @math λ_max = b / T,  b = 2.897771955×10⁻³ m·K
 */
double hspec_wien_peak_wavelength(double T);

/**
 * Stefan-Boltzmann law: total exitance from a blackbody.
 *
 * @param T  Temperature [K]
 * @return   M = σ·T⁴ [W·m⁻²]
 *
 * @math M = σ T⁴,  σ = 5.670374419×10⁻⁸ W·m⁻²·K⁻⁴
 */
double hspec_stefan_boltzmann_exitance(double T);

/**
 * Beer-Lambert-Bouguer attenuation.
 *
 * @param I0              Incident intensity/radiance array [nbands]
 * @param absorption_coeff Absorption coefficient α(λ) [m⁻¹]
 * @param path_length     Geometric path length [m]
 * @param nbands          Number of spectral bands
 * @param transmitted_out Output transmitted radiance [nbands]
 * @return                0 on success, -1 on null input
 *
 * @math I(λ) = I₀(λ) · exp(-α(λ) · d)
 *
 * Complexity: O(nbands)
 * Reference: Bouguer (1729), Lambert (1760), Beer (1852)
 */
int hspec_beer_lambert_transmission(const double *I0,
                                     const double *absorption_coeff,
                                     double path_length, size_t nbands,
                                     double *transmitted_out);

/**
 * Optical depth from transmission.
 *
 * @param transmittance  T(λ) ∈ (0, 1] [nbands]
 * @param nbands         Number of bands
 * @param tau_out        Output optical depth τ(λ) [nbands]
 * @return               0 on success, -1 on error
 *
 * @math τ(λ) = -ln(T(λ))
 */
int hspec_optical_depth(const double *transmittance, size_t nbands,
                         double *tau_out);

/**
 * Compute spectral emissivity from reflectance via Kirchhoff's law.
 *
 * For opaque materials at thermal equilibrium:
 *   ε(λ) = 1 - ρ(λ)  (if transmittance τ ≈ 0)
 *
 * @param reflectance  ρ(λ) ∈ [0,1] [nbands]
 * @param nbands       Number of bands
 * @param emissivity_out Output ε(λ) [nbands]
 * @return             0 on success
 *
 * Reference: Kirchhoff (1860)
 */
int hspec_emissivity_from_reflectance(const double *reflectance, size_t nbands,
                                       double *emissivity_out);

/**
 * Continuum removal (convex-hull baseline removal).
 *
 * Identifies the convex hull of a reflectance spectrum and normalizes
 * by it to isolate absorption features.
 *
 * @param reflectance   Original spectrum [nbands]
 * @param wavelengths   Corresponding wavelengths [nbands]
 * @param nbands        Number of bands
 * @param cr            Output continuum-removed result (pre-allocated)
 * @return              0 on success, -1 on error
 *
 * Complexity: O(nbands²) for hull identification
 * Reference: Clark & Roush (1984)
 */
int hspec_continuum_removal(const double *reflectance, const double *wavelengths,
                             size_t nbands, hspec_continuum_removed_t *cr);

/**
 * Free continuum-removed data structure.
 */
void hspec_continuum_removed_free(hspec_continuum_removed_t *cr);

/**
 * Classify a spectral signature into predefined types based on
 * band ratios and absorption feature patterns.
 *
 * @param reflectance   Spectrum [nbands]
 * @param wavelengths   Wavelengths [nbands]
 * @param nbands        Number of bands
 * @return              Signature type enum value
 *
 * Complexity: O(nbands)
 */
hspec_signature_type_t hspec_classify_signature(const double *reflectance,
                                                 const double *wavelengths,
                                                 size_t nbands);

/**
 * Compute the Simple Ratio vegetation index.
 *
 * @math SR = ρ_NIR / ρ_RED
 *
 * @param nir_band  NIR reflectance (~800nm)
 * @param red_band  RED reflectance (~670nm)
 * @return          SR value (clamped to reasonable range)
 */
double hspec_simple_ratio(double nir_band, double red_band);

/**
 * Detect and characterize absorption features in a continuum-removed spectrum.
 *
 * @param cr           Continuum-removed spectrum
 * @param min_depth    Minimum depth threshold (0-1)
 * @return             Number of features detected
 *
 * Complexity: O(nbands)
 */
size_t hspec_detect_absorption_features(hspec_continuum_removed_t *cr,
                                         double min_depth);

/**
 * Compute atmospheric transmission for a given wavelength through
 * a simplified MODTRAN-like parameterized model.
 *
 * @param wavelength   Wavelength [μm]
 * @param water_vapor  Water vapor column [g·cm⁻²]
 * @param aerosol_od   Aerosol optical depth at 550nm
 * @param co2_ppm      CO₂ concentration [ppm]
 * @return             Atmospheric transmission ∈ [0, 1]
 */
double hspec_atmospheric_transmission(double wavelength, double water_vapor,
                                       double aerosol_od, double co2_ppm);

#endif /* HYPERSPECTRAL_SPECTROSCOPY_H */