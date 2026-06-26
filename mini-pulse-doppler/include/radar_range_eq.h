#ifndef RADAR_RANGE_EQ_H
#define RADAR_RANGE_EQ_H

#include <stddef.h>
#include <stdint.h>
#include "pulse_doppler.h"

/* ---------------------------------------------------------------------------
 * L1/L4: Radar Range Equation Definitions
 *
 * The radar range equation is the fundamental law relating transmitted power,
 * antenna gains, target RCS, range, and received SNR.
 *
 * Monostatic form (L4):
 *   Pr = (Pt * Gt * Gr * λ² * σ) / ((4π)³ * R⁴ * L)
 *
 * SNR form:
 *   SNR = (Pt * Gt * Gr * λ² * σ * τ * Npulses) /
 *         ((4π)³ * R⁴ * k * T0 * F * L)
 *
 * Reference: Skolnik "Radar Handbook" (2008), Ch.2
 *            Richards et al. "Principles of Modern Radar" (2010), Ch.2
 * --------------------------------------------------------------------------- */

/** Radar system parameters for range equation */
typedef struct {
    double tx_power_watts;
    double tx_antenna_gain_db;
    double rx_antenna_gain_db;
    double wavelength_m;
    double center_frequency_hz;
    double noise_figure_db;
    double system_losses_db;
    double bandwidth_hz;
    double pulse_width_s;
    double num_coherent_pulses;
    double noise_temperature_k;
} radar_system_params_t;

/** Target Radar Cross Section (RCS) model */
typedef enum {
    RCS_CONSTANT = 0,
    RCS_SPHERE,
    RCS_FLAT_PLATE,
    RCS_CORNER_REFLECTOR,
    RCS_CYLINDER,
    RCS_DIPOLE,
    RCS_SWERLING_MODEL,
    RCS_USER_DEFINED
} rcs_model_t;

/** RCS descriptor */
typedef struct {
    rcs_model_t model;
    double rcs_value_dbsm;
    double characteristic_dimension_m;
    double aspect_angle_deg;
    double frequency_hz;
    swerling_case_t swerling_case;
} rcs_target_t;

/** Propagation model */
typedef enum {
    PROP_FREE_SPACE = 0,
    PROP_TWO_RAY,
    PROP_KNIFE_EDGE,
    PROP_ROUND_EARTH,
    PROP_ITU_R_P528,
    PROP_LONGLEY_RICE
} propagation_model_t;

/** Atmospheric attenuation parameters */
typedef struct {
    double oxygen_atten_db_per_km;
    double water_vapor_atten_db_per_km;
    double rain_atten_db_per_km;
    double fog_atten_db_per_km;
    double total_atmospheric_atten_db;
    double frequency_ghz;
    double rain_rate_mm_per_hr;
} atmospheric_atten_t;

/** Range equation result */
typedef struct {
    double snr_db;
    double received_power_dbm;
    double max_detection_range_m;
    double snr_at_range_db;
    double propagation_loss_db;
    double one_way_loss_db;
} range_eq_result_t;

/** Blake chart parameters (comprehensive range calculation) */
typedef struct {
    double tx_power_dbm;
    double tx_line_loss_db;
    double tx_antenna_gain_db;
    double tx_beamwidth_deg;
    double radome_loss_db;
    double propagation_factor_db;
    double rx_antenna_gain_db;
    double rx_line_loss_db;
    double rx_noise_figure_db;
    double integration_gain_db;
    double signal_processing_gain_db;
    double required_snr_db;
} blake_chart_params_t;

/* ---------------------------------------------------------------------------
 * L5: Range Equation Computation Algorithms
 * --------------------------------------------------------------------------- */

/** L4: Monostatic radar range equation — SNR at given range */
int radar_range_snr(const radar_system_params_t *sys,
                    double range_m,
                    double target_rcs_dbsm,
                    range_eq_result_t *result);

/** L4: Maximum detection range for given minimum SNR */
int radar_max_range(const radar_system_params_t *sys,
                    double target_rcs_dbsm,
                    double min_snr_db,
                    double *max_range_m);

/** L4: Bistatic radar range equation */
int bistatic_range_snr(const radar_system_params_t *sys,
                       double tx_range_m,
                       double rx_range_m,
                       double bistatic_rcs_dbsm,
                       range_eq_result_t *result);

/** L4: Radar horizon range (4/3 Earth radius model) */
int radar_horizon_range(double antenna_height_m,
                        double target_height_m,
                        double *horizon_range_m);

/** L4: Friis transmission equation (one-way link budget) */
int friis_transmission(double tx_power_dbm,
                       double tx_gain_db,
                       double rx_gain_db,
                       double frequency_hz,
                       double range_m,
                       double *rx_power_dbm);

/** L3: Free-space path loss */
int free_space_path_loss(double frequency_hz,
                         double range_m,
                         double *loss_db);

/** L3: Two-ray ground reflection propagation loss */
int two_ray_path_loss(double tx_height_m,
                      double rx_height_m,
                      double range_m,
                      double frequency_hz,
                      double reflection_coefficient,
                      double *loss_db);

/** L3: Knife-edge diffraction loss */
int knife_edge_diffraction_loss(double obstacle_height_m,
                                double tx_height_m,
                                double rx_height_m,
                                double d1_m,
                                double d2_m,
                                double frequency_hz,
                                double *loss_db);

/** L1: RCS of a conducting sphere (Mie/Rayleigh/Optical regions) */
int rcs_sphere_compute(double radius_m,
                       double frequency_hz,
                       double *rcs_dbsm);

/** L1: RCS of a flat rectangular plate */
int rcs_flat_plate_compute(double width_m,
                           double height_m,
                           double frequency_hz,
                           double angle_deg,
                           double *rcs_dbsm);

/** L1: RCS of a trihedral corner reflector */
int rcs_corner_reflector_compute(double edge_length_m,
                                 double frequency_hz,
                                 double *rcs_dbsm);

/** L1: RCS of a cylinder */
int rcs_cylinder_compute(double radius_m,
                         double length_m,
                         double frequency_hz,
                         double angle_deg,
                         double *rcs_dbsm);

/** L3: Atmospheric attenuation (ITU-R P.676 model) */
int atmospheric_attenuation_compute(double frequency_ghz,
                                    double range_km,
                                    double temperature_c,
                                    double humidity_percent,
                                    double rain_rate_mm_per_hr,
                                    atmospheric_atten_t *atten);

/** L6: Link budget (Blake chart) for complete system analysis */
int blake_chart_compute(const blake_chart_params_t *params,
                        double target_rcs_dbsm,
                        double range_m,
                        range_eq_result_t *result);

/** L7: SNR vs range curve (range profile) */
int snr_vs_range_profile(const radar_system_params_t *sys,
                         double target_rcs_dbsm,
                         double r_min_m,
                         double r_max_m,
                         uint32_t num_points,
                         double *ranges_m,
                         double *snr_values_db);

#endif /* RADAR_RANGE_EQ_H */
