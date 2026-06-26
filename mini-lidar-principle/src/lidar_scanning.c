/**
 * @file    lidar_scanning.c
 * @brief   Laser beam propagation and scanning pattern implementation
 *
 * Knowledge covered:
 *   L1: Gaussian beam parameters, beam divergence effects
 *   L3: Gaussian beam waist, Rayleigh range, beam radius
 *   L4: Gaussian beam propagation law, diffraction limit
 *   L5: Scan pattern generation (raster, sinusoidal, Lissajous, spiral,
 *       rosette, circular), point spacing, ground range computation
 *
 * Reference:
 *   - Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., Wiley, 2019.
 *   - Siegman, A.E., *Lasers*, University Science Books, 1986.
 */

#include "lidar_scanning.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Gaussian beam initialization
 * ═══════════════════════════════════════════════════════════════════════════ */

void lidar_beam_init(lidar_gaussian_beam_t *beam,
                      double wavelength, double w0, double M2) {
    if (!beam) return;
    beam->wavelength = wavelength;
    beam->w0 = w0;
    beam->M2 = (M2 >= 1.0) ? M2 : 1.0;

    /* Rayleigh range: z_R = pi * w0^2 / (lambda * M2) */
    beam->z_R = M_PI * w0 * w0 / (wavelength * beam->M2);

    /* Far-field divergence half-angle: theta = lambda * M2 / (pi * w0) */
    beam->theta_div = wavelength * beam->M2 / (M_PI * w0);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L4: Gaussian beam propagation
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_beam_radius(const lidar_gaussian_beam_t *beam, double z) {
    if (!beam || beam->z_R <= 0.0) return 0.0;
    /* w(z) = w0 * sqrt(1 + (z/z_R)^2) * sqrt(M2)
       Actually M^2 just scales the divergence, so w(z)^2 = w0^2 * (1 + (z*lambda*M2/(pi*w0^2))^2)
       = w0^2 + (z * lambda * M2 / (pi * w0))^2 */
    double z_eff = z;
    double wz_sq = beam->w0 * beam->w0
                   + (z_eff * beam->wavelength * beam->M2 / (M_PI * beam->w0))
                     * (z_eff * beam->wavelength * beam->M2 / (M_PI * beam->w0));
    return sqrt(wz_sq);
}

double lidar_beam_intensity(const lidar_gaussian_beam_t *beam, double z) {
    if (!beam || beam->w0 <= 0.0) return 0.0;
    double wz = lidar_beam_radius(beam, z);
    if (wz <= 0.0) return 0.0;
    /* I(z)/I(0) = (w0/w(z))^2 */
    double ratio = beam->w0 / wz;
    return ratio * ratio;
}

double lidar_beam_enclosed_power(const lidar_gaussian_beam_t *beam,
                                   double z, double r) {
    if (!beam) return 0.0;
    double wz = lidar_beam_radius(beam, z);
    if (wz <= 0.0) return (r > 0.0) ? 1.0 : 0.0;
    /* P(r)/P_total = 1 - exp(-2 * r^2 / w(z)^2) */
    double arg = -2.0 * r * r / (wz * wz);
    return 1.0 - exp(arg);
}

double lidar_beam_footprint(const lidar_gaussian_beam_t *beam, double range) {
    if (!beam) return 0.0;
    return 2.0 * lidar_beam_radius(beam, range);
}

double lidar_aperture_coupling(const lidar_gaussian_beam_t *beam,
                                 double range, double aperture) {
    if (!beam || range <= 0.0 || aperture <= 0.0) return 0.0;

    /* Lambertian scattering assumption:
       Received fraction = (aperture area) / (pi * range^2)
       This is the fraction of scattered power collected by the receiver
       for a Lambertian target. */
    double A_r = M_PI * (aperture / 2.0) * (aperture / 2.0);
    double hemisphere_area = 2.0 * M_PI * range * range;
    double coupling = A_r / hemisphere_area;

    /* Cap at 1.0 */
    return (coupling < 1.0) ? coupling : 1.0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Scan pattern generation
 * ═══════════════════════════════════════════════════════════════════════════ */

void lidar_scan_pattern_angles(const lidar_scan_pattern_config_t *config,
                                double t, double *azimuth, double *elev) {
    if (!config || !azimuth || !elev) return;

    double fov_h = config->fov_h;
    double fov_v = config->fov_v;
    double fx = config->f_x;
    double fy = config->f_y;
    double omega_x = 2.0 * M_PI * fx;
    double omega_y = 2.0 * M_PI * fy;

    switch (config->pattern) {
        case SCAN_PATTERN_RASTER: {
            /* Raster: horizontal fast axis, vertical slow axis */
            double period_h = 1.0 / fx;
            double period_v = period_h * config->num_lines;
            double phase_h = fmod(t, period_h) / period_h; /* 0 to 1 */
            double phase_v = fmod(t, period_v) / period_v;
            *azimuth = fov_h * (phase_h - 0.5);
            *elev    = fov_v * (0.5 - phase_v); /* top to bottom */
            break;
        }
        case SCAN_PATTERN_SINUSOIDAL: {
            /* Resonant sinusoidal scan (Lissajous with integer ratio) */
            *azimuth = (fov_h / 2.0) * sin(omega_x * t);
            *elev    = (fov_v / 2.0) * sin(omega_y * t);
            break;
        }
        case SCAN_PATTERN_LISSAJOUS: {
            /* General Lissajous pattern */
            *azimuth = (fov_h / 2.0) * sin(omega_x * t);
            *elev    = (fov_v / 2.0) * sin(omega_y * t + M_PI / 4.0);
            break;
        }
        case SCAN_PATTERN_SPIRAL: {
            /* Archimedean spiral: r = a * theta
               theta increases linearly with time */
            double theta = omega_x * t;
            double r_max = fmin(fov_h, fov_v) / 2.0;
            double max_theta = 2.0 * M_PI * 5.0; /* 5 turns */
            double r = r_max * fmod(theta, max_theta) / max_theta;
            *azimuth = r * cos(theta);
            *elev    = r * sin(theta);
            break;
        }
        case SCAN_PATTERN_ROSETTE: {
            /* Rosette: (a*cos(w1*t) + b*cos(w2*t), a*sin(w1*t) + b*sin(w2*t)) */
            double w1 = omega_x;
            double w2 = omega_y;
            double a = fov_h / 4.0;
            double b = fov_h / 4.0;
            *azimuth = a * cos(w1 * t) + b * cos(w2 * t);
            *elev    = a * sin(w1 * t) + b * sin(w2 * t);
            break;
        }
        case SCAN_PATTERN_CIRCULAR: {
            /* Conical scan: circular pattern at fixed offset from center */
            double offset_angle = fov_h / 4.0;
            *azimuth = offset_angle * cos(omega_x * t);
            *elev    = offset_angle * sin(omega_x * t);
            break;
        }
        default:
            *azimuth = 0.0;
            *elev    = 0.0;
    }
}

size_t lidar_scan_pattern_generate(const lidar_scan_pattern_config_t *config,
                                     double *azimuths, double *elevations,
                                     size_t max_points, double dt) {
    if (!config || !azimuths || !elevations || max_points == 0) return 0;

    for (size_t i = 0; i < max_points; i++) {
        double t = (double)i * dt;
        lidar_scan_pattern_angles(config, t, &azimuths[i], &elevations[i]);
    }
    return max_points;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Scan geometry computations
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_ground_range(double altitude, double elev) {
    if (altitude <= 0.0) return 0.0;
    /* x_g = H * tan(elevation_from_nadir)
       Here elev is elevation above horizontal (0 = level).
       Convert to angle from nadir: theta_nadir = pi/2 - elev */
    double theta_nadir = M_PI / 2.0 - elev;
    /* Ground range from nadir point */
    return altitude * tan(theta_nadir);
}

void lidar_point_spacing(const lidar_scan_pattern_config_t *config,
                          double altitude, double speed,
                          double *dx, double *dy) {
    if (!config || !dx || !dy) return;

    /* Along-track spacing: dx = v_platform / f_scan */
    double f_scan = config->f_y;
    if (f_scan < 1e-9) f_scan = config->f_x;
    if (f_scan < 1e-9) f_scan = 1.0;
    *dx = speed / f_scan;

    /* Cross-track spacing at nadir: dy = 2*H*tan(FOV_h/2) / N_lines */
    double half_fov = config->fov_h / 2.0;
    double swath_width = 2.0 * altitude * tan(half_fov);
    double n_lines = config->num_lines;
    if (n_lines < 1.0) n_lines = config->points_per_line;
    if (n_lines < 1.0) n_lines = 100.0;
    *dy = swath_width / n_lines;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L5: Diffraction limit
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_diffraction_limit(double wavelength, double aperture) {
    if (wavelength <= 0.0 || aperture <= 0.0) return 0.0;
    /* Rayleigh criterion: theta_min = 1.22 * lambda / D */
    return 1.22 * wavelength / aperture;
}