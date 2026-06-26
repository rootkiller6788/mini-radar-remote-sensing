/**
 * @file    lidar_core.h
 * @brief   Core LiDAR definitions, data structures, and system parameters
 *
 * @details Defines the fundamental types and constants for LiDAR (Light
 *          Detection and Ranging) systems.  Covers the LiDAR range equation,
 *          point representations, system configuration, and performance
 *          metrics.
 *
 * Knowledge covered:
 *   L1: LiDAR range equation, time-of-flight, range resolution,
 *       angular resolution, beam divergence, PRF, laser pulse energy,
 *       point cloud representation, cross-section/reflectance
 *   L3: Range equation mathematical form, photon statistics structures
 *   L4: LiDAR range equation (derived from radar equation)
 *
 * Reference:
 *   - Richards, Scheer & Holm, *Principles of Modern Radar* (2010), Ch.13.
 *   - Wehr & Lohr, "Airborne Laser Scanning", *ISPRS JPRS* 54(2-3), 1999.
 */

#ifndef LIDAR_CORE_H
#define LIDAR_CORE_H

#include <stddef.h>
#include <stdint.h>


/* Physical constants for optical / LiDAR domain */
#define LIDAR_C                          299792458.0
#define LIDAR_PLANCK_H                   6.62607015e-34
#define LIDAR_ELECTRON_Q                 1.602176634e-19
#define LIDAR_BOLTZMANN_K                1.380649e-23

/* L1: LiDAR core enumerations */
typedef enum {
    LIDAR_TYPE_PULSED       = 0,
    LIDAR_TYPE_FMCW         = 1,
    LIDAR_TYPE_PHASE_BASED  = 2,
    LIDAR_TYPE_GEIGER       = 3,
    LIDAR_TYPE_FLASH        = 4,
    LIDAR_TYPE_COHERENT     = 5
} lidar_type_t;

typedef enum {
    LIDAR_WAVELENGTH_355NM  = 355,
    LIDAR_WAVELENGTH_532NM  = 532,
    LIDAR_WAVELENGTH_905NM  = 905,
    LIDAR_WAVELENGTH_1064NM = 1064,
    LIDAR_WAVELENGTH_1550NM = 1550,
    LIDAR_WAVELENGTH_2000NM = 2000
} lidar_wavelength_t;

typedef enum {
    LIDAR_SCAN_NONE         = 0,
    LIDAR_SCAN_OPTO_MECH    = 1,
    LIDAR_SCAN_GALVO        = 2,
    LIDAR_SCAN_RISLEY       = 3,
    LIDAR_SCAN_MEMS         = 4,
    LIDAR_SCAN_OPA          = 5,
    LIDAR_SCAN_FIBER_SWITCH = 6
} lidar_scan_type_t;

typedef enum {
    LIDAR_DETECTOR_PIN      = 0,
    LIDAR_DETECTOR_APD      = 1,
    LIDAR_DETECTOR_SPAD     = 2,
    LIDAR_DETECTOR_SIPM     = 3,
    LIDAR_DETECTOR_PMT      = 4,
    LIDAR_DETECTOR_BALANCED = 5
} lidar_detector_type_t;

/* L1: 3D point with intensity — fundamental LiDAR data element */
typedef struct {
    double x;            /* X coordinate in local frame [m] */
    double y;            /* Y coordinate in local frame [m] */
    double z;            /* Z coordinate (elevation) [m] */
    double intensity;    /* Return intensity [0, 1] or raw ADC counts */
    uint8_t return_num;  /* Return number (1=first, 2=second, etc.) */
    uint8_t num_returns; /* Total returns for this pulse */
    uint16_t scan_angle; /* Scan angle [0.001 deg units] */
    double gps_time;     /* GPS time of acquisition [s] */
    uint8_t class_label; /* Classification label (0=unclassified) */
    double pulse_width;  /* Measured pulse width at detector [ns] */
    double peak_ampl;    /* Peak amplitude of return [V or ADC counts] */
} lidar_point_t;

/* L1: LiDAR system configuration parameters */
typedef struct {
    lidar_type_t          type;
    lidar_wavelength_t    wavelength;
    double                pulse_energy;
    double                pulse_width;
    double                prf;
    double                beam_divergence;
    double                beam_diameter;
    lidar_detector_type_t detector;
    double                aperture_diam;
    double                opt_transmit;
    double                opt_receive;
    double                detector_neb;
    double                detector_bw;
    double                excess_noise_F;
    double                responsivity;
    lidar_scan_type_t     scan_type;
    double                scan_fov_h;
    double                scan_fov_v;
    double                scan_rate;
    double                ang_res_h;
    double                ang_res_v;
    double                sample_rate;
    double                range_min;
    double                range_max;
} lidar_config_t;

/* L1: LiDAR performance metrics */
typedef struct {
    double snr;
    double received_power;
    double nefd;
    double range_resolution;
    double ang_resolution;
    double max_range;
    double prob_detection;
    double prob_false_alarm;
    double range_accuracy;
    double cross_range_res;
} lidar_performance_t;

/* L1: Target optical properties */
typedef struct {
    double reflectivity;
    double backscatter_xsec;
    double area;
    double normal_angle;
    int    is_distributed;
    double roughness_rms;
} lidar_target_t;

/* L1: Atmospheric transmission parameters */
typedef struct {
    double extinction_coeff;
    double visibility;
    double aerosol_conc;
    double two_way_transmit;
    double background_rad;
    double solar_zenith_angle;
} lidar_atmosphere_t;

/* L1: Complete LiDAR scan (point cloud frame) */
typedef struct {
    lidar_point_t *points;
    size_t         num_points;
    size_t         capacity;
    double         start_time;
    double         end_time;
    double         origin_lat;
    double         origin_lon;
    double         origin_alt;
    double         roll;
    double         pitch;
    double         yaw;
} lidar_scan_t;

/* ─── L4: LiDAR Range Equation ─────────────────────────────────────────── */

/**
 * Compute received power using the LiDAR range equation.
 * Lambertian point target (standard form):
 *   P_r = P_t * rho * A_r * eta_sys * T^2 / (pi * R^2)
 * where P_t = pulse_energy / pulse_width (peak power).
 */
double lidar_range_equation_received_power(const lidar_config_t *config,
                                            const lidar_target_t *target,
                                            const lidar_atmosphere_t *atm,
                                            double range);

/**
 * Maximum detection range for given SNR threshold.
 * Solves range equation via Newton-Raphson iteration.
 */
double lidar_max_range(const lidar_config_t *config,
                       const lidar_target_t *target,
                       const lidar_atmosphere_t *atm,
                       double snr_min_db);

/**
 * Compute SNR for LiDAR measurement (linear-mode detector).
 * Includes signal shot, background shot, dark current, and thermal noise.
 */
double lidar_snr(const lidar_config_t *config,
                 const lidar_target_t *target,
                 const lidar_atmosphere_t *atm,
                 double range);

/**
 * Range resolution: delta_R = c * tau_pulse / 2  (pulsed)
 *                   delta_R = c / (2 * B_chirp)   (FMCW)
 */
double lidar_range_resolution(const lidar_config_t *config);

/**
 * Cross-range resolution: delta_x = R * theta_div  (beam limited)
 */
double lidar_cross_range_resolution(const lidar_config_t *config, double range);

/* ─── L1: Time-of-Flight core conversions ──────────────────────────────── */

double lidar_tof_to_range(double tof);
double lidar_range_to_tof(double range);
double lidar_unambiguous_range(double prf);

/* ─── L1: Point cloud management ──────────────────────────────────────── */

void   lidar_scan_init(lidar_scan_t *scan, size_t capacity);
int    lidar_scan_add_point(lidar_scan_t *scan, const lidar_point_t *point);
void   lidar_scan_free(lidar_scan_t *scan);
int    lidar_scan_bounding_box(const lidar_scan_t *scan, double min_xyz[3], double max_xyz[3]);
size_t lidar_scan_size(const lidar_scan_t *scan);

/* ─── L1: Configuration initialization ────────────────────────────────── */

void lidar_config_init_automotive(lidar_config_t *config);
void lidar_config_init_airborne(lidar_config_t *config);
void lidar_config_init_terrestrial(lidar_config_t *config);
int  lidar_config_validate(const lidar_config_t *config);

/* ─── L1: Atmosphere computation ──────────────────────────────────────── */

double lidar_atmosphere_extinction(double visibility, double wavelength);
double lidar_atmosphere_transmission(double extinction_coeff, double range);

/* ─── L4: Poisson photon statistics ───────────────────────────────────── */

double lidar_expected_signal_pe(const lidar_config_t *config,
                                 double range,
                                 const lidar_target_t *target,
                                 const lidar_atmosphere_t *atm);
double lidar_poisson_prob(int k, double lambda);
int    lidar_compute_performance(const lidar_config_t *config,
                                  const lidar_target_t *target,
                                  const lidar_atmosphere_t *atm,
                                  double range,
                                  lidar_performance_t *perf);

#endif /* LIDAR_CORE_H */
