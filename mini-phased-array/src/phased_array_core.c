/* ============================================================================
 * phased_array_core.c — Phased Array Core Implementations
 *
 * Implements fundamental phased array operations:
 *   - Array geometry initialization (linear, planar, circular)
 *   - Coordinate transforms (spherical <-> azimuth/elevation <-> unit vectors)
 *   - Element pattern models (analytical)
 *   - Array factor computation (Pattern Multiplication Theorem)
 *   - Wavelength / wavenumber utilities
 *   - Memory management for array structures
 *
 * Reference Textbooks:
 *   Balanis, C.A. (2016) "Antenna Theory: Analysis and Design", 4th ed. Ch. 6
 *   Mailloux, R.J. (2005) "Phased Array Antenna Handbook", 2nd ed.
 *   Richards, M.A., Scheer, J.A., Holm, W.A. (2010) "Principles of Modern Radar"
 *   Skolnik, M.I. (2008) "Radar Handbook", 3rd ed. Ch. 13
 * ============================================================================ */

#include "phased_array.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#define C0 299792458.0
#define DEG2RAD(d) ((d) * M_PI / 180.0)
#define RAD2DEG(r) ((r) * 180.0 / M_PI)

/* === L1: Wavelength and Wavenumber === */

double pa_wavelength(double frequency_hz)
{
    if (frequency_hz <= 0.0) return 0.0;
    return C0 / frequency_hz;
}

double pa_wavenumber(double frequency_hz)
{
    if (frequency_hz <= 0.0) return 0.0;
    return (2.0 * M_PI * frequency_hz) / C0;
}

/* === L3: Coordinate Transformations === */

pa_spherical_t pa_azel_to_spherical(double azimuth_deg, double elevation_deg)
{
    pa_spherical_t s;
    s.theta = M_PI / 2.0 - DEG2RAD(elevation_deg);
    if (s.theta < 0.0) s.theta = 0.0;
    if (s.theta > M_PI) s.theta = M_PI;
    s.phi = DEG2RAD(azimuth_deg);
    while (s.phi < 0.0) s.phi += 2.0 * M_PI;
    while (s.phi >= 2.0 * M_PI) s.phi -= 2.0 * M_PI;
    return s;
}

pa_azel_t pa_spherical_to_azel(double theta_rad, double phi_rad)
{
    pa_azel_t ae;
    ae.elevation_deg = 90.0 - RAD2DEG(theta_rad);
    if (ae.elevation_deg > 90.0) ae.elevation_deg = 90.0;
    if (ae.elevation_deg < -90.0) ae.elevation_deg = -90.0;
    ae.azimuth_deg = RAD2DEG(phi_rad);
    while (ae.azimuth_deg > 180.0) ae.azimuth_deg -= 360.0;
    while (ae.azimuth_deg < -180.0) ae.azimuth_deg += 360.0;
    return ae;
}

void pa_spherical_to_uvw(double theta, double phi,
                         double *u, double *v, double *w)
{
    double st = sin(theta);
    *u = st * cos(phi);
    *v = st * sin(phi);
    *w = cos(theta);
}

pa_spherical_t pa_uvw_to_spherical(double u, double v, double w)
{
    pa_spherical_t s;
    double norm = sqrt(u*u + v*v + w*w);
    if (norm < 1e-15) { s.theta = 0.0; s.phi = 0.0; return s; }
    double wn = w / norm;
    if (wn > 1.0) wn = 1.0; if (wn < -1.0) wn = -1.0;
    s.theta = acos(wn);
    double un = u / norm, vn = v / norm;
    if (fabs(un) < 1e-12 && fabs(vn) < 1e-12) { s.phi = 0.0; }
    else { s.phi = atan2(vn, un); if (s.phi < 0.0) s.phi += 2.0 * M_PI; }
    return s;
}

/* === L2: Element Radiation Pattern Models === */

double complex pa_element_pattern(pa_element_type_t element_type,
                                   double theta, double phi)
{
    double st;
    (void) phi;
    switch (element_type) {
    case PA_ELEMENT_ISOTROPIC:
        return 1.0 + 0.0 * I;
    case PA_ELEMENT_HALF_WAVE_DIPOLE:
        st = sin(theta);
        if (fabs(st) < 1e-10) return 0.0 + 0.0 * I;
        return (cos(M_PI / 2.0 * cos(theta)) / st) + 0.0 * I;
    case PA_ELEMENT_PATCH_RECT:
        if (theta > M_PI / 2.0) return 0.0 + 0.0 * I;
        return (cos(theta) > 0.0 ? cos(theta) : 0.0) + 0.0 * I;
    case PA_ELEMENT_PATCH_CIRCULAR:
        if (theta > M_PI / 2.0) return 0.0 + 0.0 * I;
        { double j0a = 1.0 - 0.25 * theta * theta;
          if (j0a < 0.0) j0a = 0.0;
          return (cos(theta) * j0a) + 0.0 * I; }
    case PA_ELEMENT_SLOT:
        if (theta > M_PI / 2.0) return 0.0 + 0.0 * I;
        return 1.0 + 0.0 * I;
    case PA_ELEMENT_VIVALDI:
        if (theta > M_PI / 2.0) return 0.0 + 0.0 * I;
        return cos(theta) + 0.0 * I;
    case PA_ELEMENT_DIPOLE_CROSSED:
    case PA_ELEMENT_PATCH_DUAL_FEED:
        if (theta > M_PI / 2.0) return 0.0 + 0.0 * I;
        return cos(theta) + 0.0 * I;
    default:
        return 1.0 + 0.0 * I;
    }
}

double pa_element_pattern_dbi(pa_element_type_t element_type,
                               double theta, double phi)
{
    double mag = cabs(pa_element_pattern(element_type, theta, phi));
    return (mag < 1e-15) ? -300.0 : 20.0 * log10(mag);
}

/* === L4: Array Factor — Pattern Multiplication Theorem ===
 *
 * AF(theta,phi) = sum_{n=0}^{N-1} w_n * exp(j k0 (r_n · u_hat))
 *
 * Total far-field pattern:
 *   E_total(theta,phi) = E_element(theta,phi) × AF(theta,phi)
 *
 * This is the cornerstone of array theory (Balanis §6.3).
 * The element pattern acts as a spatial envelope that limits
 * the useful scan range of the array.
 */

pa_af_result_t pa_array_factor(const pa_array_config_t *config,
                                const pa_element_t *elements,
                                const double complex *weights,
                                double theta, double phi)
{
    pa_af_result_t r;
    memset(&r, 0, sizeof(r));
    r.theta_rad = theta; r.phi_rad = phi;

    if (!config || !elements || !weights || config->num_elements == 0) {
        r.af_magnitude_db = -300.0; return r;
    }

    uint32_t N = config->num_elements;
    double u, v, w;
    pa_spherical_to_uvw(theta, phi, &u, &v, &w);

    double k0 = pa_wavenumber(config->frequency_hz);
    if (k0 <= 0.0) { r.af_magnitude_db = -300.0; return r; }

    double complex sum = 0.0 + 0.0 * I;
    for (uint32_t n = 0; n < N; n++) {
        double dot = elements[n].position.x * u
                   + elements[n].position.y * v
                   + elements[n].position.z * w;
        double phase = k0 * dot;
        sum += weights[n] * (cos(phase) + sin(phase) * I);
    }

    r.af_complex   = sum;
    r.af_magnitude = cabs(sum);
    r.af_phase_rad = carg(sum);
    r.af_magnitude_db = (r.af_magnitude < 1e-15) ? -300.0
                        : 20.0 * log10(r.af_magnitude);
    r.directivity_dbi = r.af_magnitude_db;
    return r;
}

pa_af_result_t pa_total_pattern(const pa_array_config_t *config,
                                 const pa_element_t *elements,
                                 const double complex *weights,
                                 double theta, double phi)
{
    pa_af_result_t r = pa_array_factor(config, elements, weights, theta, phi);
    pa_element_type_t et = (elements && config && config->num_elements > 0)
                           ? elements[0].element_type : PA_ELEMENT_ISOTROPIC;
    double ep_mag = cabs(pa_element_pattern(et, theta, phi));
    r.af_complex   *= ep_mag;
    r.af_magnitude *= ep_mag;
    r.af_magnitude_db = (r.af_magnitude < 1e-15) ? -300.0
                        : 20.0 * log10(r.af_magnitude);
    return r;
}

void pa_af_normalize_db(pa_af_result_t *result, double max_af_magnitude)
{
    if (!result) return;
    if (max_af_magnitude < 1e-15 || result->af_magnitude < 1e-15) {
        result->af_magnitude_db = -300.0; return;
    }
    result->af_magnitude_db = 20.0 * log10(result->af_magnitude / max_af_magnitude);
}

/* === L5: Array Geometry Initialization === */

void pa_init_linear_array(pa_array_config_t *config, uint32_t num_elements,
                          double element_spacing, double frequency_hz)
{
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->geometry          = PA_GEOMETRY_LINEAR;
    config->num_elements      = num_elements;
    config->frequency_hz      = frequency_hz;
    config->element_spacing_x = element_spacing;
    config->rows = 1; config->cols = num_elements;
    config->scan_limit_az_deg = 60.0;
    config->scan_limit_el_deg = 60.0;
}

void pa_init_planar_array(pa_array_config_t *config,
                          uint32_t rows, uint32_t cols,
                          double dx, double dy, double frequency_hz)
{
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->geometry          = PA_GEOMETRY_PLANAR_RECT;
    config->num_elements      = rows * cols;
    config->frequency_hz      = frequency_hz;
    config->element_spacing_x = dx;
    config->element_spacing_y = dy;
    config->rows = rows; config->cols = cols;
    config->scan_limit_az_deg = 60.0;
    config->scan_limit_el_deg = 60.0;
}

void pa_init_circular_array(pa_array_config_t *config, uint32_t num_elements,
                            double radius, double frequency_hz)
{
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->geometry            = PA_GEOMETRY_CIRCULAR;
    config->num_elements        = num_elements;
    config->frequency_hz        = frequency_hz;
    config->element_spacing_radius = radius;
    config->rows = 1; config->cols = num_elements;
    config->scan_limit_az_deg   = 180.0;
    config->scan_limit_el_deg   = 60.0;
}

/* === L5: Element Array Allocation === */

pa_element_t *pa_allocate_elements(const pa_array_config_t *config)
{
    if (!config || config->num_elements == 0) return NULL;
    uint32_t N = config->num_elements;
    pa_element_t *el = (pa_element_t *)calloc(N, sizeof(pa_element_t));
    if (!el) return NULL;

    switch (config->geometry) {
    case PA_GEOMETRY_LINEAR: {
        double off = -((double)(N - 1)) / 2.0 * config->element_spacing_x;
        for (uint32_t i = 0; i < N; i++) {
            el[i].position.x = off + (double)i * config->element_spacing_x;
            el[i].excitation = 1.0 + 0.0 * I;
            el[i].element_type = PA_ELEMENT_ISOTROPIC;
            el[i].tr_module_id = -1;
        }
        break;
    }
    case PA_GEOMETRY_PLANAR_RECT: {
        uint32_t R = config->rows, C = config->cols;
        double dx = config->element_spacing_x;
        double dy = config->element_spacing_y;
        double x0 = -((double)(C - 1)) / 2.0 * dx;
        double y0 = -((double)(R - 1)) / 2.0 * dy;
        for (uint32_t r = 0; r < R; r++) {
            for (uint32_t c = 0; c < C; c++) {
                uint32_t idx = r * C + c;
                el[idx].position.x = x0 + (double)c * dx;
                el[idx].position.y = y0 + (double)r * dy;
                el[idx].excitation = 1.0 + 0.0 * I;
                el[idx].element_type = PA_ELEMENT_ISOTROPIC;
                el[idx].tr_module_id = -1;
            }
        }
        break;
    }
    case PA_GEOMETRY_CIRCULAR: {
        double Rc = config->element_spacing_radius;
        for (uint32_t i = 0; i < N; i++) {
            double ang = 2.0 * M_PI * (double)i / (double)N;
            el[i].position.x = Rc * cos(ang);
            el[i].position.y = Rc * sin(ang);
            el[i].excitation = 1.0 + 0.0 * I;
            el[i].element_type = PA_ELEMENT_ISOTROPIC;
            el[i].tr_module_id = -1;
        }
        break;
    }
    default:
        break;
    }
    return el;
}

void pa_free_elements(pa_element_t *elements) { free(elements); }

/* === L5: Pattern Memory Management === */

pa_pattern_t *pa_allocate_pattern(uint32_t num_theta, uint32_t num_phi,
                                   double th_s, double th_e,
                                   double ph_s, double ph_e)
{
    if (num_theta == 0 || num_phi == 0) return NULL;
    pa_pattern_t *p = (pa_pattern_t *)calloc(1, sizeof(pa_pattern_t));
    if (!p) return NULL;
    p->af_db = (double *)calloc((size_t)num_theta * num_phi, sizeof(double));
    if (!p->af_db) { free(p); return NULL; }
    p->num_theta = num_theta; p->num_phi = num_phi;
    p->theta_start = th_s; p->theta_end = th_e;
    p->phi_start = ph_s;   p->phi_end = ph_e;
    return p;
}

void pa_free_pattern(pa_pattern_t *pattern)
{
    if (!pattern) return;
    free(pattern->af_db);
    free(pattern);
}

/* === L4: Array Directivity ===
 *
 * D_max = 4*pi*U_max / P_rad
 * For uniform linear array: D ≈ 2*N*d/lambda (d <= lambda/2)
 * For planar: D ≈ D_x * D_y
 * For circular: D ≈ N
 */

double pa_directivity(const pa_array_config_t *config,
                      const pa_element_t *elements,
                      const double complex *weights)
{
    if (!config || !elements || config->num_elements == 0) return 0.0;
    (void) weights;
    double lambda = pa_wavelength(config->frequency_hz);
    if (lambda <= 0.0) return 0.0;
    double N = (double)config->num_elements;

    switch (config->geometry) {
    case PA_GEOMETRY_LINEAR: {
        double d_l = config->element_spacing_x / lambda;
        if (d_l > 0.5) d_l = 0.5;
        return 2.0 * N * d_l;
    }
    case PA_GEOMETRY_PLANAR_RECT: {
        double dx_l = config->element_spacing_x / lambda;
        double dy_l = config->element_spacing_y / lambda;
        if (dx_l > 0.5) dx_l = 0.5;
        if (dy_l > 0.5) dy_l = 0.5;
        return (2.0 * config->cols * dx_l) * (2.0 * config->rows * dy_l);
    }
    case PA_GEOMETRY_CIRCULAR:
        return N;
    default:
        return N;
    }
}

/* === L6: 3dB Beamwidth ===
 *
 * ULA broadside: HPBW ≈ 0.886*lambda/(N*d) rad
 * Scanned:        HPBW(theta_s) ≈ HPBW(0) / cos(theta_s)
 * UCA:            HPBW ≈ 29.1*lambda/(2*R) deg
 */

void pa_beamwidth_3db(const pa_array_config_t *config,
                      const pa_element_t *elements,
                      const double complex *weights,
                      double theta_steer, double phi_steer,
                      double *hp_beamwidth_deg)
{
    (void) elements; (void) weights; (void) phi_steer;
    if (!config || config->num_elements == 0) {
        *hp_beamwidth_deg = 0.0; return;
    }
    double lambda = pa_wavelength(config->frequency_hz);
    if (lambda <= 0.0) { *hp_beamwidth_deg = 0.0; return; }

    switch (config->geometry) {
    case PA_GEOMETRY_LINEAR: {
        double hpbw_r = 0.886 * lambda / (config->num_elements * config->element_spacing_x);
        double cs = cos(theta_steer);
        if (fabs(cs) < 0.02) cs = 0.02;
        *hp_beamwidth_deg = RAD2DEG(hpbw_r / cs);
        break;
    }
    case PA_GEOMETRY_PLANAR_RECT: {
        double hpbw_r = 0.886 * lambda / (config->cols * config->element_spacing_x);
        double cs = cos(theta_steer);
        if (fabs(cs) < 0.02) cs = 0.02;
        *hp_beamwidth_deg = RAD2DEG(hpbw_r / cs);
        break;
    }
    case PA_GEOMETRY_CIRCULAR:
        *hp_beamwidth_deg = (config->element_spacing_radius > 0.0)
            ? 29.1 * lambda / (2.0 * config->element_spacing_radius) : 0.0;
        break;
    default:
        *hp_beamwidth_deg = 0.0;
    }
}

/* === L7: Radar Range Equation ===
 *
 * R^4 = Pt*Gt*Gr*lambda^2*sigma / ((4*pi)^3*k*T0*B*F*SNR_min*L)
 *
 * Typical AESA: 1000+ T/R modules, each ~10W, detect 1m^2 RCS at 100 km+
 */

double pa_radar_range_equation(const pa_array_config_t *config,
                                const pa_element_t *elements,
                                double target_rcs_m2,
                                double snr_min_db, double system_loss_db)
{
    if (!config || !elements || target_rcs_m2 <= 0.0) return 0.0;
    double lambda = pa_wavelength(config->frequency_hz);
    if (lambda <= 0.0) return 0.0;

    const double kb = 1.380649e-23;
    const double T0 = 290.0;
    const double B  = 1.0e6;
    const double F  = 2.0;

    double SNR_min = pow(10.0, snr_min_db / 10.0);
    double L       = pow(10.0, system_loss_db / 10.0);

    double Pt = 0.0;
    for (uint32_t i = 0; i < config->num_elements; i++)
        Pt += elements[i].tr_module_power_watt;
    if (Pt <= 0.0) Pt = 2.5 * (double)config->num_elements;

    double G = pa_directivity(config, elements, NULL);
    if (G <= 0.0) G = (double)config->num_elements;

    double num = Pt * G * G * lambda * lambda * target_rcs_m2;
    double den = pow(4.0 * M_PI, 3.0) * kb * T0 * B * F * SNR_min * L;
    if (den <= 0.0) return 0.0;
    double R4 = num / den;
    return (R4 > 0.0) ? pow(R4, 0.25) : 0.0;
}

double pa_effective_aperture(const pa_array_config_t *config,
                             const pa_element_t *elements,
                             double aperture_efficiency)
{
    if (!config || !elements || config->num_elements == 0) return 0.0;
    double lambda = pa_wavelength(config->frequency_hz);
    if (lambda <= 0.0) return 0.0;
    double N = (double)config->num_elements;
    double Ae = 0.0;

    switch (config->geometry) {
    case PA_GEOMETRY_LINEAR:
        Ae = lambda * config->element_spacing_x / (4.0 * M_PI); break;
    case PA_GEOMETRY_PLANAR_RECT:
        Ae = config->element_spacing_x * config->element_spacing_y; break;
    case PA_GEOMETRY_CIRCULAR:
        Ae = (2.0 * M_PI * config->element_spacing_radius / N) * (lambda / 2.0); break;
    default:
        Ae = lambda * lambda / (4.0 * M_PI); break;
    }

    if (aperture_efficiency < 0.0) aperture_efficiency = 0.0;
    if (aperture_efficiency > 1.0) aperture_efficiency = 1.0;
    return aperture_efficiency * N * Ae;
}