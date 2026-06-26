/**
 * @file    ir_detector.h
 * @brief   Infrared Detector Models and Figures of Merit
 *
 * Knowledge Coverage:
 *   L1 - Definitions: responsivity, NEP, D*, NETD, MRTD, quantum efficiency
 *   L2 - Core Concepts: photon vs thermal detectors, cooling methods
 *   L4 - Fundamental Laws: Johnson-Nyquist noise, shot noise, 1/f noise
 *   L5 - Algorithms: detector characterization, noise analysis
 *
 * Detector types covered:
 *   Thermal: microbolometer (VOx, a-Si), thermopile, pyroelectric
 *   Photon: InSb, MCT (HgCdTe), QWIP, Type-II superlattice, PbS/PbSe
 *
 * References:
 *   Rogalski, A. (2011) "Infrared Detectors", 2nd Ed., CRC Press
 *   Vincent, J.D. (2015) "Fundamentals of IR Detector Operation & Testing"
 *
 * Curriculum Mapping:
 *   MIT 6.450 Digital Comm (detection theory)
 *   Stanford EE247 Optical (photodetectors)
 *   Berkeley EE105 Analog (noise analysis)
 *   TU Munich High-Frequency (detector physics)
 */

#ifndef IR_DETECTOR_H
#define IR_DETECTOR_H

#include <stddef.h>
#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Detector Types */
typedef enum {
    IR_DET_THERMAL_BOLOMETER = 0,
    IR_DET_THERMAL_THERMOPILE,
    IR_DET_THERMAL_PYROELECTRIC,
    IR_DET_PHOTON_IN_SB,
    IR_DET_PHOTON_MCT,
    IR_DET_PHOTON_QWIP,
    IR_DET_PHOTON_T2SL,
    IR_DET_PHOTON_PBS,
    IR_DET_COUNT
} ir_detector_type_t;

/* Cooling modes */
typedef enum {
    IR_DET_MODE_UNCOOLED = 0,
    IR_DET_MODE_THERMOELECTRIC,
    IR_DET_MODE_STIRLING,
    IR_DET_MODE_JOULE_THOMSON
} ir_detector_cooling_t;

/* Material properties */
typedef struct {
    const char *material;
    double bandgap_eV;
    double cutoff_wavelength_um;
    double operating_temperature_K;
    double thermal_coeff_K_inv;
} ir_detector_material_t;

/* Complete detector characterization */
typedef struct {
    ir_detector_type_t    type;
    ir_detector_cooling_t cooling;
    ir_detector_material_t material;
    double pixel_pitch_um;
    double fill_factor;
    double active_area_cm2;
    double quantum_efficiency;
    double spectral_response_peak_um;
    double spectral_bandwidth_um;
    double resistance_ohm;
    double bias_voltage_V;
    double capacitance_F;
    double thermal_conductance_W_per_K;
    double thermal_mass_J_per_K;
    double thermal_time_constant_ms;
    double integration_time_us;
    double frame_rate_Hz;
} ir_detector_params_t;

/* Initialization */
int ir_detector_params_init(ir_detector_params_t *p, ir_detector_type_t type);
const ir_detector_material_t* ir_detector_material_get(ir_detector_type_t type);

/* Responsivity */
double ir_detector_responsivity_voltage(const ir_detector_params_t *params, double wavelength_um);
double ir_detector_responsivity_current(double quantum_efficiency, double wavelength_um);
double ir_detector_responsivity_thermal(const ir_detector_params_t *params, double freq_Hz);

/* Noise mechanisms */
double ir_noise_johnson_voltage(double temperature_K, double resistance_ohm, double bandwidth_Hz);
double ir_noise_shot_current(double dc_current_A, double bandwidth_Hz);
double ir_noise_flicker_psd(double dc_current_A, double frequency_Hz, double k_factor, double alpha, double beta);
double ir_noise_total_current(const ir_detector_params_t *params, double bandwidth_Hz, double frequency_Hz);

/* Figures of Merit */
double ir_nep(const ir_detector_params_t *params, double bandwidth_Hz, double frequency_Hz, double wavelength_um);
double ir_specific_detectivity(const ir_detector_params_t *params, double nep_W, double bandwidth_Hz);
double ir_d_star_blip(double quantum_eff, double wavelength_um, double photon_flux_bg);
double ir_netd(double f_number, double tau_opt, double detector_area_cm2, double nep_W, double dL_dT);
double ir_radiance_derivative_wrt_temperature(double T_K, double lambda1_um, double lambda2_um, int n_steps);
double ir_mrtd(double netd_K, double spatial_freq_cyc_per_mrad, double mtf_value, double ifov_mrad, double eye_integration_s, double frame_rate_Hz);

/* Microbolometer thermal model */
double ir_microbolometer_temp_rise(const ir_detector_params_t *params, double incident_power_W, double freq_Hz);
double ir_microbolometer_time_constant(const ir_detector_params_t *params);

/* Quantum efficiency */
double ir_quantum_efficiency_external(double reflectance, double absorption_coeff_cm_inv, double thickness_cm, double internal_qe);
double ir_cutoff_wavelength(double bandgap_eV);
double ir_d_star_compare(double d_star1, double d_star2);

#endif /* IR_DETECTOR_H */
