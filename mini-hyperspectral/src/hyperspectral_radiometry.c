/**
 * @file    hyperspectral_radiometry.c
 * @brief   Radiometric calibration and atmospheric correction.
 *
 * @details Implements DN-to-radiance conversion, TOA reflectance computation,
 *          atmospheric correction (DOS, empirical line, simplified radiative
 *          transfer), BRDF modeling (Minnaert, Ross-Li), solar irradiance,
 *          and water vapor estimation.
 *
 * Knowledge covered:
 *   L1: DN, radiance, reflectance, BRDF definitions
 *   L2: Radiometric calibration chain, atmospheric correction concepts
 *   L4: Inverse-square law, Lambert cosine law
 *   L6: Atmospheric correction of AVIRIS data, water vapor retrieval
 *   L7: Application — NASA AVIRIS data processing, solar irradiance modeling
 *
 * Reference:
 *   - Schott (2007), "Remote Sensing: The Image Chain Approach"
 *   - Chavez (1996), PE&RS
 *   - Vermote et al. (1997), IEEE TGRS (6S model)
 */

#include "hyperspectral_radiometry.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "radiometry: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── DN to Radiance ──────────────────────────────────────────────────── */

/**
 * @math L_sensor(λ_i) = gain_i · DN_i + offset_i
 *
 * Complexity: O(nbands)
 */
int hspec_dn_to_radiance(const uint16_t *dn, const hspec_cal_coefficient_t *cal_coeffs,
                          size_t nbands, double *radiance_out)
{
    if (!dn || !cal_coeffs || !radiance_out || nbands == 0) return -1;

    for (size_t i = 0; i < nbands; i++) {
        radiance_out[i] = cal_coeffs[i].gain * (double)dn[i] + cal_coeffs[i].offset;
        /* Apply nonlinearity correction if significant */
        if (fabs(cal_coeffs[i].nonlinearity_coeff) > 1e-15) {
            double L = radiance_out[i];
            radiance_out[i] += cal_coeffs[i].nonlinearity_coeff * L * L;
        }
    }
    return 0;
}

/* ─── Radiance to TOA Reflectance ────────────────────────────────────── */

/**
 * @math ρ_TOA(λ) = π·L(λ)·d² / (E_sun(λ)·cos(θ_s))
 *
 * where d = Earth-Sun distance in AU (≈ 1 for average),
 *       θ_s = solar zenith angle (radians from zenith).
 *
 * Reference: Chander et al. (2009), RSE
 */
int hspec_radiance_to_toa_refl(const double *radiance, const double *solar_irrad,
                                size_t nbands, double earth_sun_au,
                                double solar_zenith, double *toa_refl_out)
{
    if (!radiance || !solar_irrad || !toa_refl_out || nbands == 0) return -1;

    double cos_theta = cos(solar_zenith);
    if (cos_theta <= 0.0) return -1;  /* sun below horizon */

    double d2 = earth_sun_au * earth_sun_au;

    for (size_t i = 0; i < nbands; i++) {
        double esun = solar_irrad[i];
        if (esun < 1e-30) {
            toa_refl_out[i] = 0.0;
            continue;
        }
        toa_refl_out[i] = M_PI * radiance[i] * d2 / (esun * cos_theta);
        /* Clamp to reasonable range [0, 2] for bright targets like clouds */
        if (toa_refl_out[i] < 0.0) toa_refl_out[i] = 0.0;
        if (toa_refl_out[i] > 2.0) toa_refl_out[i] = 2.0;
    }
    return 0;
}

/* ─── Dark Object Subtraction (DOS) ───────────────────────────────────── */

/**
 * Simple DOS atmospheric correction.
 *
 * For each band, find the dark-object DN (typically 1st percentile),
 * subtract it as an estimate of path radiance, then convert to
 * approximate surface reflectance.
 *
 * @math ρ_surface(λ) ≈ π·(L(λ) - L_path(λ))·d² / (E_sun(λ)·cos(θ_s))
 *
 * Complexity: O(n_pixels·nbands)
 * Reference: Chavez (1988, 1996)
 */
int hspec_dos_correction(const hspec_datacube_t *dc,
                          const hspec_dos_params_t *params,
                          double *corrected_out)
{
    if (!dc || !dc->data || !corrected_out) return -1;

    size_t N = dc->npixels, B = dc->nbands;
    double percentile = (params && params->dark_object_percentile > 0.0)
                        ? params->dark_object_percentile : 0.01;

    /* Compute dark object value per band */
    for (size_t b = 0; b < B; b++) {
        /* Extract and sort this band */
        double *band_data = safe_malloc(N * sizeof(double));
        for (size_t p = 0; p < N; p++)
            band_data[p] = dc->data[p * B + b];

        /* Simple sort (bubble - OK for small N, suboptimal for large) */
        for (size_t i = 0; i < N - 1; i++)
            for (size_t j = i + 1; j < N; j++)
                if (band_data[j] < band_data[i]) {
                    double t = band_data[i]; band_data[i] = band_data[j]; band_data[j] = t;
                }

        /* Get dark object value at given percentile */
        size_t idx = (size_t)(percentile * (double)(N - 1));
        if (idx >= N) idx = N - 1;
        double dark_val = band_data[idx];

        /* Subtract dark value from all pixels */
        for (size_t p = 0; p < N; p++) {
            double val = dc->data[p * B + b] - dark_val;
            corrected_out[p * B + b] = (val > 0.0) ? val : 0.0;
        }

        free(band_data);
    }
    return 0;
}

/* ─── Empirical Line Calibration ─────────────────────────────────────── */

/**
 * Uses known reflectance targets to solve per-band linear calibration.
 *
 * For each band b, solve:
 *   ρ_target = a_b · L_target + b_b
 * using least squares on n_targets known points.
 *
 * Complexity: O(n_pixels·nbands + n_targets·nbands)
 */
int hspec_empirical_line_correction(const hspec_datacube_t *dc,
                                     const double *known_reflectance,
                                     const size_t *target_pixels,
                                     size_t n_targets,
                                     double *corrected_out)
{
    if (!dc || !dc->data || !known_reflectance || !target_pixels
        || !corrected_out || n_targets < 2) return -1;

    size_t N = dc->npixels, B = dc->nbands;

    for (size_t b = 0; b < B; b++) {
        /* Solve linear regression: ρ = a·L + b */
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
        for (size_t t = 0; t < n_targets; t++) {
            double L = dc->data[target_pixels[t] * B + b];
            double rho = known_reflectance[t * B + b];
            sum_x += L;
            sum_y += rho;
            sum_xy += L * rho;
            sum_xx += L * L;
        }
        double denom = (double)n_targets * sum_xx - sum_x * sum_x;
        double slope = 0.0, intercept = 0.0;
        if (fabs(denom) > 1e-30) {
            slope = ((double)n_targets * sum_xy - sum_x * sum_y) / denom;
            intercept = (sum_y - slope * sum_x) / (double)n_targets;
        }

        /* Apply per-band */
        for (size_t p = 0; p < N; p++) {
            double val = slope * dc->data[p * B + b] + intercept;
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
            corrected_out[p * B + b] = val;
        }
    }
    return 0;
}

/* ─── Simplified Radiative Transfer Atmospheric Correction ──────────── */

hspec_atm_correction_result_t hspec_atmospheric_correction(
    const hspec_datacube_t *dc, const hspec_atm_params_t *atm_params,
    const hspec_solar_irradiance_t *solar_irr)
{
    hspec_atm_correction_result_t ac;
    memset(&ac, 0, sizeof(ac));
    if (!dc || !dc->data || !atm_params) return ac;

    size_t N = dc->npixels, B = dc->nbands;
    ac.nbands = B;
    ac.path_radiance = safe_malloc(B * sizeof(double));
    ac.transmittance_down = safe_malloc(B * sizeof(double));
    ac.transmittance_up = safe_malloc(B * sizeof(double));
    ac.spherical_albedo = safe_malloc(B * sizeof(double));
    ac.surface_reflectance = safe_malloc(N * B * sizeof(double));

    /* Per-band atmospheric parameters using simplified 6S model */
    for (size_t b = 0; b < B; b++) {
        double wl = dc->wavelengths[b];  /* in μm */

        /* Rayleigh optical depth */
        double tau_ray = 0.008735 * pow(wl, -4.08);
        double pres_ratio = atm_params->surface_pressure_mb / 1013.25;
        tau_ray *= pres_ratio;

        /* Aerosol optical depth (Angstrom law) */
        double tau_aer = atm_params->aerosol_od_550
                        * pow(0.55 / wl, 1.3);

        double tau_total = tau_ray + tau_aer;

        /* Transmittance (Beer-Lambert through atmosphere) */
        double airmass_down = 1.0 / cos(atm_params->solar_zenith_deg * M_PI / 180.0);
        double airmass_up = 1.0 / cos(atm_params->view_zenith_deg * M_PI / 180.0);
        if (airmass_down > 10.0) airmass_down = 10.0;
        if (airmass_up > 10.0) airmass_up = 10.0;

        ac.transmittance_down[b] = exp(-tau_total * airmass_down);
        ac.transmittance_up[b] = exp(-tau_total * airmass_up);

        /* Path radiance: simplified as fraction of TOA radiance */
        ac.path_radiance[b] = 0.02 + 0.05 * tau_total;  /* heuristic */

        /* Spherical albedo */
        ac.spherical_albedo[b] = 0.1 * tau_aer;

        /* Direct/diffuse irradiance at surface */
        ac.direct_irradiance[b] = solar_irr->solar_irradiance[b]
                                  * ac.transmittance_down[b];
        ac.diffuse_irradiance[b] = solar_irr->solar_irradiance[b]
                                   * (0.1 * tau_total);
    }

    /* Surface reflectance retrieval using 6S formula */
    for (size_t p = 0; p < N; p++) {
        for (size_t b = 0; b < B; b++) {
            double L_sensor = dc->data[p * B + b];
            double L_path = ac.path_radiance[b];

            /* Simplification: ρ = (L_sensor - L_path) / (T_up * E_surface/π) */
            double Esurface = ac.direct_irradiance[b] + ac.diffuse_irradiance[b];
            double numerator = M_PI * (L_sensor - L_path);
            double denominator = ac.transmittance_up[b] * Esurface
                               + ac.spherical_albedo[b] * numerator;

            double rho = 0.0;
            if (fabs(denominator) > 1e-30)
                rho = numerator / denominator;

            if (rho < 0.0) rho = 0.0;
            if (rho > 2.0) rho = 2.0;

            ac.surface_reflectance[p * B + b] = rho;
        }
    }

    return ac;
}

void hspec_atm_correction_result_free(hspec_atm_correction_result_t *ac)
{
    if (!ac) return;
    free(ac->path_radiance);
    free(ac->transmittance_down);
    free(ac->transmittance_up);
    free(ac->spherical_albedo);
    free(ac->direct_irradiance);
    free(ac->diffuse_irradiance);
    free(ac->surface_reflectance);
    free(ac->water_vapor_map);
    free(ac->aod_map);
    memset(ac, 0, sizeof(*ac));
}

/* ─── Solar Irradiance ────────────────────────────────────────────────── */

/**
 * Thuillier (2004) solar irradiance at 1 AU.
 * Tabulated key reference values (simplified model).
 *
 * @param lambda  Wavelength [μm]
 * @return        Solar spectral irradiance [W·m⁻²·μm⁻¹]
 */
double hspec_solar_irradiance_at_wavelength(double lambda)
{
    /* Piecewise model based on Thuillier 2004 / ASTM E490 */
    if (lambda < 0.2) return 0.0;
    if (lambda < 0.3) return 100.0 + (lambda - 0.2) * 500.0;
    if (lambda < 0.4) return 300.0 + (lambda - 0.3) * 1000.0;
    if (lambda < 0.5) return 1850.0;
    if (lambda < 0.6) return 1950.0 - (lambda - 0.5) * 100.0;
    if (lambda < 0.7) return 1550.0;
    if (lambda < 0.8) return 1200.0 - (lambda - 0.7) * 500.0;
    if (lambda < 1.0) return 1000.0;
    if (lambda < 1.5) return 500.0 - (lambda - 1.0) * 400.0;
    if (lambda < 2.0) return 200.0;
    if (lambda < 2.5) return 80.0;
    if (lambda < 3.0) return 30.0 - (lambda - 2.5) * 20.0;
    if (lambda < 4.0) return 10.0;
    if (lambda < 10.0) return 5.0 - (lambda - 4.0) * 0.5;
    return 1.0;
}

int hspec_solar_irradiance_init(const double *bands, size_t nbands,
                                 hspec_solar_irradiance_t *sirr)
{
    if (!bands || !sirr || nbands == 0) return -1;
    memset(sirr, 0, sizeof(*sirr));
    sirr->nbands = nbands;
    sirr->wavelengths = safe_malloc(nbands * sizeof(double));
    sirr->solar_irradiance = safe_malloc(nbands * sizeof(double));
    memcpy(sirr->wavelengths, bands, nbands * sizeof(double));

    sirr->total_solar_irradiance = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        sirr->solar_irradiance[i] = hspec_solar_irradiance_at_wavelength(bands[i]);
        sirr->total_solar_irradiance += sirr->solar_irradiance[i];
    }
    strcpy(sirr->model_name, "Thuillier2004_simplified");
    return 0;
}

void hspec_solar_irradiance_free(hspec_solar_irradiance_t *sirr)
{
    if (!sirr) return;
    free(sirr->wavelengths);
    free(sirr->solar_irradiance);
    memset(sirr, 0, sizeof(*sirr));
}

/* ─── Minnaert BRDF Correction ──────────────────────────────────────── */

/**
 * @math ρ_corrected = ρ_observed · (cos θ_i)^(1-k) · (cos θ_v)^(1-k)
 *
 * Reference: Minnaert (1941)
 */
int hspec_brdf_minnaert_correct(const double *reflectance_obs, size_t nbands,
                                 double solar_zenith, double view_zenith,
                                 const double *minnaert_k, double *corrected_out)
{
    if (!reflectance_obs || !corrected_out || nbands == 0) return -1;

    double cos_i = cos(solar_zenith);
    double cos_v = cos(view_zenith);
    if (cos_i <= 0.0 || cos_v <= 0.0) return -1;

    for (size_t b = 0; b < nbands; b++) {
        double k = minnaert_k ? minnaert_k[b] : 0.7;
        double factor = pow(cos_i, 1.0 - k) * pow(cos_v, 1.0 - k);
        corrected_out[b] = reflectance_obs[b] * factor;
        if (corrected_out[b] < 0.0) corrected_out[b] = 0.0;
        if (corrected_out[b] > 1.0) corrected_out[b] = 1.0;
    }
    return 0;
}

/* ─── Ross-Li BRDF ─────────────────────────────────────────────────── */

/**
 * Ross-Thick kernel (volumetric scattering):
 *   K_vol = ((π/2 - ξ)·cos ξ + sin ξ) / (cos θ_i + cos θ_v) - π/4
 *
 * Li-Sparse kernel (geometric-optical):
 *   K_geo = O − sec θ_i' − sec θ_v' + 0.5·(1 + cos ξ')·sec θ_v'
 *
 * @math ρ = f_iso + f_vol·K_vol + f_geo·K_geo
 *
 * Reference: Roujean et al. (1992), JGR
 */
double hspec_brdf_ross_li(double f_iso, double f_vol, double f_geo,
                           double solar_zenith, double view_zenith,
                           double rel_azimuth)
{
    double cos_i = cos(solar_zenith);
    double cos_v = cos(view_zenith);
    if (cos_i <= 0.0 || cos_v <= 0.0) return f_iso;

    double sin_i = sin(solar_zenith);
    double sin_v = sin(view_zenith);
    double cos_phi = cos(rel_azimuth);

    /* Phase angle ξ */
    double cos_xi = cos_i * cos_v + sin_i * sin_v * cos_phi;
    if (cos_xi > 1.0) cos_xi = 1.0;
    if (cos_xi < -1.0) cos_xi = -1.0;
    double xi = acos(cos_xi);
    double sin_xi = sin(xi);

    /* Ross-Thick kernel */
    double K_vol = 0.0;
    if (cos_i + cos_v > 1e-30) {
        K_vol = ((M_PI / 2.0 - xi) * cos_xi + sin_xi)
              / (cos_i + cos_v) - M_PI / 4.0;
    }

    /* Li-Sparse kernel */
    /* Compute secant of transformed angles */
    double tan_i = sin_i / cos_i;
    double tan_v = sin_v / cos_v;
    double sec_i_prime = sqrt(1.0 + tan_i * tan_i);
    double sec_v_prime = sqrt(1.0 + tan_v * tan_v);

    /* Distance between sun and view directions */
    double D = sqrt(tan_i * tan_i + tan_v * tan_v
                   - 2.0 * tan_i * tan_v * cos_phi);

    /* Simplified phase angle */
    double t = (D > 0.01) ? acos(cos_i * cos_v + sin_i * sin_v * cos_phi) : 0.0;
    double O = (1.0 / M_PI) * (t - sin(t) * cos(t)) * (sec_i_prime + sec_v_prime);

    double cos_xi_prime = cos_i * cos_v + sin_i * sin_v * cos_phi;
    double K_geo = O - sec_i_prime - sec_v_prime
                   + 0.5 * (1.0 + cos_xi_prime) * sec_v_prime;

    double rho = f_iso + f_vol * K_vol + f_geo * K_geo;
    if (rho < 0.0) rho = 0.0;
    return rho;
}

/* ─── Water Vapor Estimation ────────────────────────────────────────── */

/**
 * Estimate column water vapor using CIBR method at 940nm.
 *
 * @math CIBR = R_measurement / (a·R_left + b·R_right)
 *
 * where the continuum is interpolated between shoulder wavelengths
 * (typically 860nm and 1040nm for the 940nm absorption feature).
 *
 * Reference: Gao & Goetz (1990), JGR
 */
double hspec_water_vapor_estimate(const double *reflectance,
                                   const double *wavelengths, size_t nbands)
{
    if (!reflectance || !wavelengths || nbands < 5) return 0.0;

    /* Find band closest to 940nm (water absorption center) */
    size_t idx_center = 0, idx_left = 0, idx_right = 0;
    double d_center = 1e9, d_left = 1e9, d_right = 1e9;

    for (size_t i = 0; i < nbands; i++) {
        double wl = wavelengths[i];
        if (fabs(wl - 940.0) < d_center) { d_center = fabs(wl - 940.0); idx_center = i; }
        if (fabs(wl - 865.0) < d_left) { d_left = fabs(wl - 865.0); idx_left = i; }
        if (fabs(wl - 1040.0) < d_right) { d_right = fabs(wl - 1040.0); idx_right = i; }
    }

    if (idx_left == idx_right) return 0.0;

    /* Linear continuum interpolation */
    double wl1 = wavelengths[idx_left],  R1 = reflectance[idx_left];
    double wl2 = wavelengths[idx_right], R2 = reflectance[idx_right];
    double wlc = wavelengths[idx_center];

    double slope = (R2 - R1) / (wl2 - wl1);
    double continuum = R1 + slope * (wlc - wl1);

    double R_center = reflectance[idx_center];

    /* CIBR ratio */
    if (continuum < 1e-10) return 0.0;
    double cibr = R_center / continuum;

    /* Empirical conversion to water vapor column [g/cm²] */
    /* Based on Gao & Goetz 1990 calibration */
    double wv = (1.0 - cibr) * 10.0;  /* rough scaling */

    if (wv < 0.0) wv = 0.0;
    if (wv > 6.0) wv = 6.0;  /* max plausible column for AVIRIS */
    return wv;
}