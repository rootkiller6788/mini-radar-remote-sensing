/**
 * @file    hyperspectral_spectroscopy.c
 * @brief   Spectroscopy: Planck, Wien, Stefan-Boltzmann, Beer-Lambert, continuum removal.
 *
 * @details Implements blackbody radiation laws, absorption/transmission physics,
 *          continuum removal for spectral feature analysis, atmospheric
 *          transmission modeling, and spectral signature classification.
 *
 * Knowledge covered:
 *   L1: Radiance, reflectance, absorption feature definitions
 *   L4: Planck's law, Wien's displacement law, Stefan-Boltzmann law,
 *       Beer-Lambert-Bouguer law, Kirchhoff's law
 *   L2: Continuum removal, absorption feature detection, signature classification
 *
 * Reference:
 *   - Planck (1901)
 *   - Goody & Yung, "Atmospheric Radiation" (1989)
 *   - Clark & Roush (1984), JGR
 */

#include "hyperspectral_spectroscopy.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "spectroscopy: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── L4: Planck's Law ───────────────────────────────────────────────────── */

/**
 * Compute blackbody spectral radiance using Planck's law.
 *
 * @math B_λ(λ, T) = (2·h·c² / λ⁵) · 1 / (exp(h·c / (λ·k·T)) - 1)
 *
 * Units:
 *   λ  in meters [m]
 *   T  in Kelvin [K]
 *   B_λ in [W·m⁻²·sr⁻¹·m⁻¹]
 *
 * Handles the Rayleigh-Jeans limit (λ·k·T ≫ h·c) and Wien limit numerically
 * by checking the exponential argument.
 *
 * Complexity: O(npoints)
 * Reference: Planck, "On the Law of Distribution of Energy in the Normal
 *   Spectrum" (1901), Annalen der Physik
 */
int hspec_planck_radiance(double T, const double *wavelengths, size_t npoints,
                           double *radiance_out)
{
    if (!wavelengths || !radiance_out || npoints == 0) return -1;
    if (T <= 0.0) return -1;

    const double c1 = 2.0 * HSPEC_PLANCK_h * HSPEC_SPEED_OF_LIGHT_c
                      * HSPEC_SPEED_OF_LIGHT_c;  /* 2hc² */
    const double c2 = HSPEC_PLANCK_h * HSPEC_SPEED_OF_LIGHT_c
                      / HSPEC_BOLTZMANN_k;       /* hc/k */

    for (size_t i = 0; i < npoints; i++) {
        double lambda = wavelengths[i];
        if (lambda <= 0.0) { radiance_out[i] = 0.0; continue; }

        double exponent = c2 / (lambda * T);
        /* For very small exponent (Rayleigh-Jeans limit: B = 2ckT/λ⁴) */
        if (exponent < 1e-10) {
            radiance_out[i] = c1 / (lambda * lambda * lambda * lambda * lambda)
                              * (1.0 / exponent);
        } else if (exponent > 700.0) {
            /* For very large exponent, exp term dominates → ~0 */
            radiance_out[i] = 0.0;
        } else {
            double denominator = exp(exponent) - 1.0;
            radiance_out[i] = c1 / (lambda * lambda * lambda * lambda * lambda)
                              / denominator;
        }
    }
    return 0;
}

/* ─── L4: Wien's Displacement Law ────────────────────────────────────────── */

/**
 * @math λ_max = b / T
 *       b = 2.897771955×10⁻³ m·K (Wien displacement constant)
 *
 * Reference: Wien (1893), "Eine neue Beziehung der Strahlung schwarzer Körper"
 */
double hspec_wien_peak_wavelength(double T)
{
    if (T <= 0.0) return -1.0;
    return HSPEC_WIEN_b / T;
}

/* ─── L4: Stefan-Boltzmann Law ───────────────────────────────────────────── */

/**
 * @math M = σ · T⁴
 *       σ = 5.670374419×10⁻⁸ W·m⁻²·K⁻⁴
 *
 * Reference: Stefan (1879), Boltzmann (1884)
 */
double hspec_stefan_boltzmann_exitance(double T)
{
    if (T <= 0.0) return 0.0;
    return HSPEC_STEFAN_BOLTZMANN_sigma * T * T * T * T;
}

/* ─── L4: Beer-Lambert-Bouguer Law ───────────────────────────────────────── */

/**
 * @math I(λ) = I₀(λ) · exp(-α(λ) · d)
 *
 * Computes the transmitted intensity through a homogeneous absorbing medium.
 * α(λ) is the absorption coefficient per unit length (m⁻¹).
 * d is the geometric path length (m).
 *
 * Complexity: O(nbands)
 * Reference: Bouguer (1729), "Essai d'optique sur la gradation de la lumière";
 *            Lambert (1760), "Photometria"; Beer (1852)
 */
int hspec_beer_lambert_transmission(const double *I0,
                                     const double *absorption_coeff,
                                     double path_length, size_t nbands,
                                     double *transmitted_out)
{
    if (!I0 || !absorption_coeff || !transmitted_out) return -1;
    if (nbands == 0 || path_length < 0.0) return -1;

    for (size_t i = 0; i < nbands; i++) {
        double tau = absorption_coeff[i] * path_length;
        /* Clamp tau to avoid underflow */
        if (tau > 700.0) {
            transmitted_out[i] = 0.0;
        } else {
            transmitted_out[i] = I0[i] * exp(-tau);
        }
    }
    return 0;
}

/* ─── Optical Depth ────────────────────────────────────────────────────── */

/**
 * @math τ(λ) = -ln(T(λ))
 *
 * For T(λ) → 0, τ → large number; for T(λ) = 1, τ = 0.
 */
int hspec_optical_depth(const double *transmittance, size_t nbands,
                         double *tau_out)
{
    if (!transmittance || !tau_out || nbands == 0) return -1;

    for (size_t i = 0; i < nbands; i++) {
        double T = transmittance[i];
        if (T <= 0.0) {
            tau_out[i] = 1e10;  /* essentially opaque */
        } else if (T >= 1.0) {
            tau_out[i] = 0.0;
        } else {
            tau_out[i] = -log(T);
        }
    }
    return 0;
}

/* ─── L4: Kirchhoff's Law — Emissivity from Reflectance ────────────────── */

/**
 * For an opaque material at thermal equilibrium: ε(λ) = 1 - ρ(λ)
 * (assuming transmittance τ ≈ 0 for thick bodies)
 *
 * Reference: Kirchhoff (1860), "Über das Verhältnis zwischen dem
 *   Emissionsvermögen und dem Absorptionsvermögen der Körper"
 */
int hspec_emissivity_from_reflectance(const double *reflectance, size_t nbands,
                                       double *emissivity_out)
{
    if (!reflectance || !emissivity_out || nbands == 0) return -1;

    for (size_t i = 0; i < nbands; i++) {
        double rho = reflectance[i];
        if (rho < 0.0) rho = 0.0;
        if (rho > 1.0) rho = 1.0;
        emissivity_out[i] = 1.0 - rho;
    }
    return 0;
}

/* ─── Continuum Removal ────────────────────────────────────────────────── */

/**
 * Continuum removal via convex hull identification.
 *
 * Algorithm:
 *   1. Find the convex hull of the reflectance spectrum as a function of
 *      wavelength. The hull forms an upper envelope connecting local maxima.
 *   2. Divide the spectrum by the hull: R_CR(λ) = R(λ) / R_hull(λ).
 *   3. This normalizes the spectrum so the continuum is at 1.0 and
 *      absorption features dip below 1.0.
 *
 * Complexity: O(nbands²) for hull construction
 * Reference: Clark & Roush (1984), JGR Solid Earth
 */
int hspec_continuum_removal(const double *reflectance, const double *wavelengths,
                             size_t nbands, hspec_continuum_removed_t *cr)
{
    if (!reflectance || !wavelengths || !cr || nbands < 3) return -1;
    memset(cr, 0, sizeof(*cr));

    cr->nbands = nbands;
    cr->original_reflectance = safe_malloc(nbands * sizeof(double));
    cr->continuum = safe_malloc(nbands * sizeof(double));
    cr->continuum_removed = safe_malloc(nbands * sizeof(double));

    memcpy(cr->original_reflectance, reflectance, nbands * sizeof(double));

    /* Step 1: Find convex hull upper envelope */
    /* Use the "top hat" approach: for each band, find the maximum of
       linear interpolations between all pairs that contain it */
    for (size_t i = 0; i < nbands; i++) {
        double max_val = reflectance[i];
        for (size_t j = 0; j < nbands; j++) {
            if (j == i) continue;
            double slope = (reflectance[j] - reflectance[i])
                         / (wavelengths[j] - wavelengths[i]);
            /* Check all bands that lie between i and j */
            size_t lo = (i < j) ? i : j;
            size_t hi = (i < j) ? j : i;
            int above_all = 1;
            for (size_t k = lo + 1; k < hi && above_all; k++) {
                double interp = reflectance[i] + slope * (wavelengths[k] - wavelengths[i]);
                if (interp < reflectance[k] - 1e-12) above_all = 0;
            }
            /* Also check the endpoints */
            if (above_all) {
                /* Line from i to j is part of or below the upper envelope */
                double interp_i = reflectance[j] + slope * (wavelengths[i] - wavelengths[j]);
                if (interp_i > max_val) max_val = interp_i;
            }
        }
        cr->continuum[i] = max_val;
    }

    /* Step 2: Iteratively refine — ensure continuum >= original everywhere */
    int changed = 1;
    while (changed) {
        changed = 0;
        for (size_t i = 0; i < nbands; i++) {
            /* Interpolate between adjacent continuum-support points */
            /* Find nearest points to left and right where continuum == reflectance */
            size_t left = 0;
            for (size_t j = i; j > 0; j--) {
                if (fabs(cr->continuum[j] - reflectance[j]) < 1e-12) {
                    left = j; break;
                }
            }
            size_t right = nbands - 1;
            for (size_t j = i; j < nbands; j++) {
                if (fabs(cr->continuum[j] - reflectance[j]) < 1e-12) {
                    right = j; break;
                }
            }
            if (left < right) {
                double slope = (cr->continuum[right] - cr->continuum[left])
                             / (wavelengths[right] - wavelengths[left]);
                double interp = cr->continuum[left]
                              + slope * (wavelengths[i] - wavelengths[left]);
                if (interp > cr->continuum[i] + 1e-15) {
                    cr->continuum[i] = interp;
                    changed = 1;
                }
            }
        }
    }

    /* Step 3: Compute continuum-removed spectrum */
    for (size_t i = 0; i < nbands; i++) {
        if (cr->continuum[i] > 1e-30)
            cr->continuum_removed[i] = reflectance[i] / cr->continuum[i];
        else
            cr->continuum_removed[i] = 1.0;
    }

    return 0;
}

void hspec_continuum_removed_free(hspec_continuum_removed_t *cr)
{
    if (!cr) return;
    free(cr->original_reflectance);
    free(cr->continuum);
    free(cr->continuum_removed);
    free(cr->features);
    memset(cr, 0, sizeof(*cr));
}

/* ─── Detect Absorption Features ────────────────────────────────────────── */

/**
 * Detect absorption features in a continuum-removed spectrum.
 * A feature is detected where R_CR(λ) < 1.0 - min_depth and
 * the local minimum depth exceeds the threshold.
 *
 * Complexity: O(nbands)
 */
size_t hspec_detect_absorption_features(hspec_continuum_removed_t *cr,
                                         double min_depth)
{
    if (!cr || !cr->continuum_removed || cr->nbands < 3) return 0;
    if (min_depth < 0.0) min_depth = 0.0;
    if (min_depth > 1.0) min_depth = 1.0;

    /* Free previous features */
    free(cr->features);
    cr->num_absorption_features = 0;
    cr->features = NULL;

    size_t nbands = cr->nbands;
    double *cr_spec = cr->continuum_removed;

    /* First pass: count features */
    size_t nfeat = 0;
    int in_feature = 0;
    for (size_t i = 0; i < nbands; i++) {
        double depth = 1.0 - cr_spec[i];
        if (depth > min_depth && !in_feature) {
            in_feature = 1;
            nfeat++;
        } else if (depth <= min_depth && in_feature) {
            in_feature = 0;
        }
    }

    if (nfeat == 0) return 0;

    cr->features = safe_malloc(nfeat * sizeof(hspec_absorption_feature_t));
    memset(cr->features, 0, nfeat * sizeof(hspec_absorption_feature_t));

    /* Second pass: characterize each feature */
    size_t feat_idx = 0;
    in_feature = 0;
    size_t feat_start = 0;
    for (size_t i = 0; i < nbands; i++) {
        double depth = 1.0 - cr_spec[i];
        if (depth > min_depth && !in_feature) {
            in_feature = 1;
            feat_start = i;
        } else if ((depth <= min_depth || i == nbands - 1) && in_feature) {
            in_feature = 0;
            size_t feat_end = i;
            if (i == nbands - 1 && depth > min_depth) feat_end = nbands;
            if (feat_end > feat_start && feat_idx < nfeat) {
                /* Find minimum (maximum depth) */
                double max_d = 0.0;
                double area = 0.0;
                for (size_t j = feat_start; j < feat_end; j++) {
                    double d = 1.0 - cr_spec[j];
                    if (d > max_d) { max_d = d; }
                    area += d;
                }
                cr->features[feat_idx].depth = max_d;
                cr->features[feat_idx].area = area;
                cr->features[feat_idx].width_fwhm =
                    (feat_end > feat_start + 1) ? (double)(feat_end - feat_start) : 1.0;
                feat_idx++;
            }
        }
    }
    cr->num_absorption_features = nfeat;
    return nfeat;
}

/* ─── Spectral Signature Classification ─────────────────────────────────── */

/**
 * Classify a spectral signature based on band ratios and absorption features.
 *
 * Decision tree:
 *   - High NIR plateau + red edge → VEGETATION
 *   - Deep absorption near 2200nm → CLAY (Al-OH)
 *   - Deep absorption near 2350nm → CARBONATE (CO₃²⁻)
 *   - Deep absorption near 900nm → IRON_OXIDE (Fe³⁺)
 *   - High VIS + low SWIR → SNOW_ICE
 *   - Monotonic increasing → RED_SLOPE
 *   - Monotonic decreasing → BLUE_SLOPE
 *   - Low SNR or flat → FLAT
 *
 * Complexity: O(nbands)
 */
hspec_signature_type_t hspec_classify_signature(const double *reflectance,
                                                 const double *wavelengths,
                                                 size_t nbands)
{
    if (!reflectance || !wavelengths || nbands < 4) return HSPEC_SIG_FLAT;

    /* Find representative band values */
    double vis_mean = 0.0, nir_mean = 0.0, swir_mean = 0.0;
    size_t vis_n = 0, nir_n = 0, swir_n = 0;

    for (size_t i = 0; i < nbands; i++) {
        double wl = wavelengths[i];
        if (wl >= 400.0 && wl <= 700.0) { vis_mean += reflectance[i]; vis_n++; }
        if (wl >= 750.0 && wl <= 1000.0) { nir_mean += reflectance[i]; nir_n++; }
        if (wl >= 1500.0 && wl <= 2500.0) { swir_mean += reflectance[i]; swir_n++; }
    }

    if (vis_n > 0) vis_mean /= (double)vis_n;
    if (nir_n > 0) nir_mean /= (double)nir_n;
    if (swir_n > 0) swir_mean /= (double)swir_n;

    /* Vegetation: NIR >> VIS (red edge) */
    if (nir_mean > vis_mean * 1.5 && nir_mean > 0.1) {
        return HSPEC_SIG_VEGETATION;
    }

    /* Snow/Ice: very high VIS, decreasing NIR/SWIR */
    if (vis_mean > 0.6 && nir_mean < vis_mean && swir_mean < nir_mean) {
        return HSPEC_SIG_SNOW_ICE;
    }

    /* Monotonic checks */
    int increasing = 1, decreasing = 1;
    for (size_t i = 1; i < nbands; i++) {
        if (reflectance[i] < reflectance[i-1] - 1e-10) increasing = 0;
        if (reflectance[i] > reflectance[i-1] + 1e-10) decreasing = 0;
    }
    if (increasing && reflectance[nbands-1] > reflectance[0] * 1.1)
        return HSPEC_SIG_RED_SLOPE;
    if (decreasing && reflectance[0] > reflectance[nbands-1] * 1.1)
        return HSPEC_SIG_BLUE_SLOPE;

    /* Check for water absorption near 1400nm and 1900nm */
    /* Look for dips >20% relative to local continuum */
    for (size_t i = 0; i < nbands; i++) {
        if (wavelengths[i] >= 1380.0 && wavelengths[i] <= 1420.0
            && reflectance[i] < 0.3) {
            return HSPEC_SIG_WATER;
        }
        if (wavelengths[i] >= 1900.0 && wavelengths[i] <= 1950.0
            && reflectance[i] < 0.3) {
            return HSPEC_SIG_WATER;
        }
    }

    /* Mineral defaults based on SWIR features */
    if (swir_mean < vis_mean * 0.8 && nir_mean > 0.05) {
        return HSPEC_SIG_CLAY;
    }

    return HSPEC_SIG_FLAT;
}

/* ─── Simple Ratio ─────────────────────────────────────────────────────── */

double hspec_simple_ratio(double nir_band, double red_band)
{
    if (fabs(red_band) < 1e-15) return (nir_band > 0) ? 100.0 : 0.0;
    return nir_band / red_band;
}

/* ─── Atmospheric Transmission (Simplified Model) ──────────────────────── */

/**
 * Simple parameterized atmospheric transmission model.
 *
 * Models:
 *   1. Rayleigh scattering: τ_Ray(λ) ∝ λ⁻⁴
 *   2. Aerosol (Mie) extinction: τ_aer(λ) ∝ λ⁻α (Angstrom exponent)
 *   3. Water vapor bands: Gaussian absorption at 940, 1140, 1380, 1870 nm
 *   4. CO₂ absorption: around 2000nm and 4300nm
 *   5. Ozone: Chappuis band ~600nm and Hartley-Huggins UV
 *
 * Complexity: O(1) per wavelength
 */
double hspec_atmospheric_transmission(double wavelength, double water_vapor,
                                       double aerosol_od, double co2_ppm)
{
    double lambda = wavelength;  /* in μm */
    if (lambda <= 0.0 || lambda > 15.0) return 0.0;

    /* Rayleigh optical depth: τ_R ≈ 0.008735 * λ^(-4.08) (Bodhaine 1999) */
    double tau_rayleigh = 0.008735 * pow(lambda, -4.08);

    /* Aerosol: Angstrom law τ_aer(λ) = τ_aer(550) * (550/λ)^α, α ≈ 1.3 */
    double tau_aerosol = aerosol_od * pow(0.55 / lambda, 1.3);

    /* Water vapor absorption bands (Gaussian) */
    double tau_h2o = 0.0;
    /* Key H₂O absorption band centers in μm */
    double h2o_centers[] = {0.72, 0.82, 0.94, 1.14, 1.38, 1.87, 2.7, 3.2, 6.3};
    double h2o_strengths[] = {0.02, 0.05, 0.30, 0.25, 0.50, 0.35, 0.40, 0.30, 0.60};
    double h2o_sigmas[]    = {0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.10};
    for (int i = 0; i < 9; i++) {
        double diff = lambda - h2o_centers[i];
        tau_h2o += water_vapor * h2o_strengths[i]
                   * exp(-diff * diff / (2.0 * h2o_sigmas[i] * h2o_sigmas[i]));
    }

    /* CO₂ absorption bands */
    double tau_co2 = 0.0;
    double co2_centers[] = {2.00, 2.70, 4.26, 4.30};
    double co2_strengths[] = {0.03, 0.05, 0.50, 0.10};
    double co2_sigmas[]    = {0.03, 0.04, 0.02, 0.02};
    for (int i = 0; i < 4; i++) {
        double diff = lambda - co2_centers[i];
        tau_co2 += (co2_ppm / 400.0) * co2_strengths[i]
                   * exp(-diff * diff / (2.0 * co2_sigmas[i] * co2_sigmas[i]));
    }

    /* Ozone (Chappuis band ~0.6 μm) */
    double tau_o3 = 0.0;
    {
        double diff = lambda - 0.60;
        tau_o3 = 0.03 * exp(-diff * diff / (2.0 * 0.04 * 0.04));
    }

    double total_tau = tau_rayleigh + tau_aerosol + tau_h2o + tau_co2 + tau_o3;
    if (total_tau > 700.0) return 0.0;
    return exp(-total_tau);
}