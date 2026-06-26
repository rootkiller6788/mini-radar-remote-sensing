/**
 * @file    lidar_detection.h
 * @brief   LiDAR detection theory — SNR, probability of detection, false alarms
 *
 * Knowledge covered:
 *   L1: SNR, NEP, NEFD, probability of detection/false alarm
 *   L2: Detection threshold, CFAR, ROC curves
 *   L3: Poisson statistics, Gaussian approximation, chi-squared distributions
 *   L4: Neyman-Pearson detection criterion
 *   L5: CFAR (Constant False Alarm Rate) detection, matched filter
 *
 * Reference:
 *   - Kingston, R.H., *Optical Sources, Detectors, and Systems*, 1995.
 *   - Kay, S.M., *Fundamentals of Statistical Signal Processing:
 *     Detection Theory*, Prentice Hall, 1998.
 *   - Goodman, J.W., *Statistical Optics*, 2nd ed., Wiley, 2015.
 */

#ifndef LIDAR_DETECTION_H
#define LIDAR_DETECTION_H

#include "lidar_core.h"
#include <stddef.h>

/* ─── L1: Detection result ──────────────────────────────────────────────── */

typedef struct {
    int    detected;         /**< 1 if target detected at this range bin */
    double range;            /**< Range of detection [m] */
    double amplitude;        /**< Measured signal amplitude */
    double snr;              /**< Measured SNR [dB] */
    double threshold;        /**< Detection threshold used */
    double prob_detection;   /**< Probability of detection for this measurement */
    double prob_false_alarm; /**< Single-look PFA at this threshold */
} lidar_detection_result_t;

/* ─── L4: Neyman-Pearson detection ──────────────────────────────────────── */

/**
 * @brief Compute detection threshold for desired PFA (Gaussian noise)
 *
 * Under H0 (noise only):   y ~ N(mu_n, sigma_n²)
 * Under H1 (signal+noise): y ~ N(A + mu_n, sigma_n²)
 *
 * For a desired PFA, the threshold is:
 *   T = mu_n + sigma_n · Q^{-1}(PFA)
 *
 * where Q^{-1} is the inverse Q-function (complementary CDF of standard normal).
 *
 * Probability of detection:
 *   P_d = Q( Q^{-1}(PFA) - sqrt(SNR_linear) )
 *
 * This is the standard radar detection model applied to LiDAR.
 *
 * Reference: Kay (1998), Ch.3; Richards et al. (2010), Ch.15.
 *
 * @param pfa        Desired probability of false alarm
 * @param noise_mean Noise mean
 * @param noise_std  Noise standard deviation
 * @return           Detection threshold
 */
double lidar_detection_threshold_gaussian(double pfa,
                                            double noise_mean,
                                            double noise_std);

/**
 * @brief Compute probability of detection for a given SNR and PFA
 *
 * P_d = Q( Q^{-1}(PFA) - sqrt(SNR_linear) )
 *
 * where SNR_linear = 10^(SNR_dB / 10).
 *
 * @param snr_db  Signal-to-noise ratio [dB]
 * @param pfa     Probability of false alarm
 * @return        Probability of detection [0-1]
 */
double lidar_prob_detection(double snr_db, double pfa);

/**
 * @brief Compute minimum SNR required for given (P_d, PFA) pair
 *
 * SNR_min = ( Q^{-1}(PFA) - Q^{-1}(P_d) )²
 *
 * This is the detectability factor — the SNR needed at the detector
 * input to achieve the required detection performance.
 *
 * @param pd   Desired probability of detection
 * @param pfa  Maximum acceptable probability of false alarm
 * @return     Minimum SNR [linear, not dB]
 */
double lidar_detectability_factor(double pd, double pfa);

/* ─── L2: Constant False Alarm Rate (CFAR) ──────────────────────────────── */

/**
 * @brief Cell-Averaging CFAR detection
 *
 * CA-CFAR estimates the local noise level from neighboring range bins
 * (guard cells excluded) and sets an adaptive threshold:
 *
 *   T_i = alpha · (1/N_ref) · Σ x_j
 *
 * where j indexes reference cells around the cell under test,
 * and alpha is the threshold multiplier computed from the desired PFA:
 *
 *   alpha = N_ref · (PFA^{-1/N_ref} - 1)
 *
 * Reference: Finn & Johnson, "Adaptive Detection Mode with Threshold
 *            Control as a Function of Spatially Sampled Clutter-Level
 *            Estimates", *RCA Review* 29, pp.414-464, 1968.
 *
 * @param signal           Input signal (range bins)
 * @param num_bins         Number of range bins
 * @param pfa              Desired probability of false alarm
 * @param guard_cells      Number of guard cells on each side
 * @param ref_cells        Number of reference cells on each side
 * @param detections       Output: detection results (1 per bin)
 * @param max_detections   Max detections to store
 * @return                 Number of detections found
 */
size_t lidar_cfar_detection(const double *signal, size_t num_bins,
                              double pfa,
                              size_t guard_cells, size_t ref_cells,
                              lidar_detection_result_t *detections,
                              size_t max_detections);

/* ─── L5: Matched filter for LiDAR pulse detection ──────────────────────── */

/**
 * @brief Apply matched filter to LiDAR waveform
 *
 * The matched filter maximizes SNR for a known pulse shape in white noise.
 * For a transmitted pulse p(t), the matched filter impulse response is:
 *
 *   h(t) = p*(-t)    (time-reversed complex conjugate)
 *
 * Output:  y(t) = x(t) ∗ h(t) = ∫ x(τ) · p(τ - t) dτ
 *
 * For a Gaussian pulse (typical LiDAR), the matched filter is a Gaussian
 * of the same width, yielding a Gaussian output with sqrt(2) wider.
 *
 * Note: This is equivalent to cross-correlation with the pulse template.
 *
 * @param waveform    Input waveform samples
 * @param num_samples Number of samples
 * @param pulse       Pulse template samples (should be same dt)
 * @param pulse_len   Length of pulse template
 * @param output      Output: matched filter response (must be pre-allocated,
 *                     length = num_samples + pulse_len - 1, or use same buffer)
 * @param out_len     Output length
 * @return            0 on success
 */
int lidar_matched_filter(const double *waveform, size_t num_samples,
                          const double *pulse, size_t pulse_len,
                          double *output, size_t *out_len);

/**
 * @brief Generate a Gaussian pulse template for matched filtering
 *
 * p(t) = exp( -t² / (2·σ²) )
 * where σ = fwhm / (2·sqrt(2·ln(2))) ≈ fwhm / 2.35482
 *
 * @param pulse    Output pulse array
 * @param length   Number of samples
 * @param dt       Sample interval [ns]
 * @param fwhm     Pulse FWHM [ns]
 */
void lidar_gaussian_pulse_template(double *pulse, size_t length,
                                     double dt, double fwhm);

/* ─── L2: Optical noise sources ─────────────────────────────────────────── */

/**
 * @brief Compute shot noise variance (Poisson → Gaussian approximation)
 *
 * For a photocurrent I, the shot noise variance is:
 *   sigma²_shot = 2 · q · I · B
 *
 * where q = electron charge, B = bandwidth.
 *
 * For APD with gain M and excess noise factor F:
 *   sigma²_shot = 2 · q · I_primary · M² · F · B
 *
 * @param current    Mean photocurrent [A]
 * @param bandwidth  Detector bandwidth [Hz]
 * @param M          APD gain (1 for PIN)
 * @param F          Excess noise factor (1 for PIN)
 * @return           Shot noise variance [A²]
 */
double lidar_shot_noise_variance(double current, double bandwidth,
                                   double M, double F);

/**
 * @brief Compute thermal (Johnson-Nyquist) noise variance
 *
 * sigma²_thermal = 4 · k · T · B / R_fb
 *
 * where k = Boltzmann constant, T = temperature [K], B = bandwidth,
 * R_fb = transimpedance feedback resistance.
 *
 * @param temperature  Temperature [K]
 * @param bandwidth    Detector bandwidth [Hz]
 * @param resistance   Feedback resistance [Ω]
 * @return             Thermal noise variance [A²]
 */
double lidar_thermal_noise_variance(double temperature,
                                      double bandwidth,
                                      double resistance);

/**
 * @brief Compute background (solar) photocurrent
 *
 * I_bg = L_solar · Δλ · A_r · Ω_r · η_opt · R_det
 *
 * where:
 *   L_solar   = solar spectral radiance [W/(m²·sr·μm)]
 *   Δλ       = optical filter bandwidth [μm]
 *   A_r       = receiver aperture area [m²]
 *   Ω_r       = receiver FOV solid angle [sr]
 *   η_opt     = optical efficiency
 *   R_det     = detector responsivity [A/W]
 *
 * @param config  LiDAR configuration
 * @param atm     Atmospheric parameters
 * @return        Background photocurrent [A]
 */
double lidar_background_current(const lidar_config_t *config,
                                  const lidar_atmosphere_t *atm);

/* ─── L1: Receiver Operating Characteristic (ROC) ───────────────────────── */

/**
 * @brief Compute ROC curve (P_d vs PFA) for Gaussian detection model
 *
 * For each PFA value in [pfa_min, pfa_max], compute the corresponding P_d.
 * This generates the standard ROC curve used to assess detector performance.
 *
 * @param snr_db    Signal-to-noise ratio [dB]
 * @param pfa       Output: array of PFA values (log-spaced)
 * @param pd        Output: array of P_d values
 * @param num_pts   Number of points on the curve
 * @param pfa_min   Minimum PFA (e.g., 1e-6)
 * @param pfa_max   Maximum PFA (e.g., 0.5)
 * @return          0 on success
 */
int lidar_roc_curve(double snr_db,
                      double *pfa, double *pd,
                      size_t num_pts,
                      double pfa_min, double pfa_max);

/* ─── L1: Range accuracy ────────────────────────────────────────────────── */

/**
 * @brief Cramer-Rao Lower Bound (CRLB) for range estimation
 *
 * For a Gaussian pulse in white noise, the CRLB for range (time-delay)
 * estimation is:
 *
 *   σ_R² ≥ c² / (8 · SNR · B_rms²)
 *
 * where B_rms is the RMS bandwidth of the pulse:
 *
 *   B_rms² = ∫ f²·|P(f)|² df / ∫ |P(f)|² df
 *
 * For a Gaussian pulse with FWHM τ:
 *   B_rms ≈ 0.187 / τ     (for Gaussian envelope)
 *   σ_R_min ≈ c·τ / (2 · sqrt(2·SNR))
 *
 * Reference: Van Trees, H.L., *Detection, Estimation, and Modulation
 *            Theory*, Part I, Wiley, 1968.
 *
 * @param snr_linear  SNR (linear, not dB)
 * @param b_rms       RMS bandwidth of pulse [Hz]
 * @return            CRLB on range standard deviation [m]
 */
double lidar_range_crlb(double snr_linear, double b_rms);

#endif /* LIDAR_DETECTION_H */