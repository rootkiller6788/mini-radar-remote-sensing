/**
 * @file    lidar_scanning.h
 * @brief   Laser beam propagation and scanning pattern generation
 *
 * Knowledge covered:
 *   L1: Beam divergence, Gaussian beam optics, scanning patterns
 *   L3: Gaussian beam waist, Rayleigh range, M² beam quality factor
 *   L4: Gaussian beam propagation law, diffraction limit
 *   L5: Scan pattern generation (raster, sinusoidal, Lissajous, spiral)
 *
 * Reference:
 *   - Saleh, B.E.A. & Teich, M.C., *Fundamentals of Photonics*, 3rd ed.,
 *     Wiley, 2019, Ch.3 (Gaussian beams).
 *   - Siegman, A.E., *Lasers*, University Science Books, 1986, Ch.17.
 */

#ifndef LIDAR_SCANNING_H
#define LIDAR_SCANNING_H

#include "lidar_core.h"
#include <stddef.h>

/* ─── L1: Gaussian beam parameters ───────────────────────────────────────── */

/**
 * @brief Gaussian beam propagation parameters
 *
 * A Gaussian beam is characterized by:
 *   w(z)   = w₀ · sqrt(1 + (z/z_R)²)   — beam radius at distance z
 *   R(z)   = z · (1 + (z_R/z)²)        — wavefront radius of curvature
 *   Θ(z)   = θ(z)                      — local divergence
 *   z_R    = π · w₀² / λ               — Rayleigh range
 *   θ_div  = λ / (π · w₀)              — far-field half-angle divergence
 *   M²     = beam quality factor (M² ≥ 1, = 1 for ideal TEM₀₀)
 */
typedef struct {
    double wavelength;    /**< Laser wavelength [m] */
    double w0;            /**< Beam waist radius (1/e² intensity) [m] */
    double z_R;           /**< Rayleigh range [m] */
    double theta_div;     /**< Far-field half-angle divergence [rad] */
    double M2;            /**< Beam quality factor M² [dimensionless] */
} lidar_gaussian_beam_t;

/**
 * @brief Initialize Gaussian beam parameters from waist size and wavelength
 *
 * Computes derived parameters:
 *   z_R = π · w₀² / λ
 *   θ_div = λ / (π · w₀)
 *
 * @param beam       Beam struct to initialize
 * @param wavelength Laser wavelength [m]
 * @param w0         Beam waist radius [m]
 * @param M2         Beam quality factor (1.0 = diffraction-limited)
 */
void lidar_beam_init(lidar_gaussian_beam_t *beam,
                      double wavelength, double w0, double M2);

/**
 * @brief Beam radius at propagation distance z
 *
 * w(z) = w₀ · sqrt(1 + (z / z_R)²), scaled by sqrt(M²) for non-ideal beams.
 *
 * @param beam  Gaussian beam parameters
 * @param z     Propagation distance from waist [m]
 * @return      Beam radius (1/e² intensity) [m]
 */
double lidar_beam_radius(const lidar_gaussian_beam_t *beam, double z);

/**
 * @brief On-axis intensity at distance z, normalized to waist intensity
 *
 * I(z) / I(0) = (w₀ / w(z))²
 *
 * @param beam  Gaussian beam parameters
 * @param z     Propagation distance from waist [m]
 * @return      Normalized on-axis intensity [0-1]
 */
double lidar_beam_intensity(const lidar_gaussian_beam_t *beam, double z);

/**
 * @brief Fraction of power within a radius r at distance z
 *
 * P(r) / P_total = 1 - exp(-2 · r² / w(z)²)
 *
 * Applications: compute aperture coupling efficiency (power through
 * receiver aperture vs. total backscattered power).
 *
 * @param beam  Gaussian beam parameters
 * @param z     Propagation distance [m]
 * @param r     Radial distance [m]
 * @return      Enclosed power fraction [0-1]
 */
double lidar_beam_enclosed_power(const lidar_gaussian_beam_t *beam,
                                   double z, double r);

/**
 * @brief Compute beam footprint diameter on target at range R
 *
 * For a Gaussian beam focused at infinity:
 *   D_footprint ≈ θ_div · R        (far-field approximation)
 * or more precisely:
 *   D_footprint = 2 · w(R)
 *
 * @param beam   Gaussian beam parameters
 * @param range  Target range [m]
 * @return       Beam footprint diameter [m]
 */
double lidar_beam_footprint(const lidar_gaussian_beam_t *beam, double range);

/**
 * @brief Receiver aperture coupling efficiency
 *
 * Fraction of reflected/scattered power captured by the receiver aperture
 * at distance R.  Accounts for:
 *   - Beam divergence (expanding spot on target)
 *   - Lambertian scattering into hemisphere
 *   - Receiver aperture solid angle
 *
 * η_coupling = A_r / (π · R² · Ω_scatter)
 *
 * For a Lambertian scatterer: Ω_scatter = π (hemisphere)
 *   η_coupling = A_r / (π · R²)        [simplified]
 *
 * @param beam      Gaussian beam parameters
 * @param range     Target range [m]
 * @param aperture  Receiver aperture diameter [m]
 * @return          Coupling efficiency [0-1]
 */
double lidar_aperture_coupling(const lidar_gaussian_beam_t *beam,
                                 double range, double aperture);

/* ─── L5: Scan pattern generation ────────────────────────────────────────── */

/**
 * @brief Scan pattern type
 */
typedef enum {
    SCAN_PATTERN_RASTER,      /**< Raster (TV-like) scan */
    SCAN_PATTERN_SINUSOIDAL,  /**< Sinusoidal (resonant mirror) */
    SCAN_PATTERN_LISSAJOUS,   /**< Lissajous (dual resonant) */
    SCAN_PATTERN_SPIRAL,      /**< Archimedean spiral */
    SCAN_PATTERN_ROSETTE,     /**< Rosette (rotating offset) */
    SCAN_PATTERN_CIRCULAR     /**< Circular scan (conical mirror) */
} lidar_scan_pattern_t;

/**
 * @brief Scan pattern configuration
 */
typedef struct {
    lidar_scan_pattern_t pattern;  /**< Pattern type */
    double fov_h;                  /**< Horizontal FOV [rad] */
    double fov_v;                  /**< Vertical FOV [rad] */
    double f_x;                    /**< Horizontal scan frequency [Hz] */
    double f_y;                    /**< Vertical scan frequency [Hz] */
    double num_lines;              /**< Number of scan lines (raster) */
    double points_per_line;        /**< Points per line (raster) */
} lidar_scan_pattern_config_t;

/**
 * @brief Generate scan pattern — returns beam pointing angles at time t
 *
 * @param config   Scan pattern configuration
 * @param t        Time [s]
 * @param azimuth  Output: azimuth angle [rad]
 * @param elev     Output: elevation angle [rad]
 */
void lidar_scan_pattern_angles(const lidar_scan_pattern_config_t *config,
                                double t, double *azimuth, double *elev);

/**
 * @brief Generate complete scan lines (array of pointing angles)
 *
 * Creates a full scan frame — useful for simulating LiDAR point cloud
 * acquisition.
 *
 * @param config       Scan pattern config
 * @param azimuths     Output: array of azimuth angles [rad]
 * @param elevations   Output: array of elevation angles [rad]
 * @param max_points   Maximum points to generate
 * @param dt           Time between consecutive measurements [s]
 * @return             Actual number of points generated
 */
size_t lidar_scan_pattern_generate(const lidar_scan_pattern_config_t *config,
                                     double *azimuths, double *elevations,
                                     size_t max_points, double dt);

/**
 * @brief Compute the scan line index on the target surface
 *
 * Given LiDAR altitude H and scan angle θ, the ground range is:
 *   x_g = H · tan(θ)     (flat Earth approximation)
 *
 * @param altitude  Sensor altitude above ground [m]
 * @param elev      Scan elevation angle from nadir [rad]
 * @return          Ground range from nadir point [m]
 */
double lidar_ground_range(double altitude, double elev);

/**
 * @brief Compute point spacing on ground for given scan pattern
 *
 * Along-track spacing:  Δx = v_platform / f_scan
 * Cross-track spacing:  Δy = 2·H·tan(FOV_h/2) / N_lines
 *
 * @param config     Scan pattern config
 * @param altitude   Platform altitude [m]
 * @param speed      Platform speed [m/s]
 * @param dx         Output: along-track spacing [m]
 * @param dy         Output: cross-track spacing [m]
 */
void lidar_point_spacing(const lidar_scan_pattern_config_t *config,
                          double altitude, double speed,
                          double *dx, double *dy);

/* ─── L5: Beam divergence effects ────────────────────────────────────────── */

/**
 * @brief Compute angular resolution from aperture size (diffraction limit)
 *
 * Rayleigh criterion for diffraction-limited optics:
 *   θ_min = 1.22 · λ / D
 *
 * where D is the aperture diameter.
 *
 * @param wavelength  Laser wavelength [m]
 * @param aperture    Aperture diameter [m]
 * @return            Diffraction-limited angular resolution [rad]
 */
double lidar_diffraction_limit(double wavelength, double aperture);

#endif /* LIDAR_SCANNING_H */