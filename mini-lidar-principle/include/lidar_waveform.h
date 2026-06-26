/**
 * @file    lidar_waveform.h
 * @brief   LiDAR full-waveform processing — pulse analysis, Gaussian decomposition
 *
 * Knowledge covered:
 *   L1: Waveform digitization, pulse characteristics, full-width half-maximum
 *   L3: Gaussian function, superposition, non-linear least squares
 *   L5: Peak detection, Gaussian decomposition, pulse fitting algorithms
 *   L6: Full-waveform LiDAR — decomposing overlapping returns
 *
 * Reference:
 *   - Wagner et al., "Gaussian Decomposition and Calibration of a Novel
 *     Small-Footprint Full-Waveform Digitising Airborne Laser Scanner",
 *     *ISPRS JPRS* 60(2), pp.100-112, 2006.
 *   - Mallet & Bretar, "Full-Waveform Topographic LiDAR: State-of-the-Art",
 *     *ISPRS JPRS* 64(1), pp.1-16, 2009.
 */

#ifndef LIDAR_WAVEFORM_H
#define LIDAR_WAVEFORM_H

#include "lidar_core.h"
#include <stddef.h>

/* ─── L1: Waveform representation ────────────────────────────────────────── */

/** Maximum number of returns detectable in one waveform */
#define LIDAR_WAVEFORM_MAX_RETURNS  16

/** Maximum waveform sample length */
#define LIDAR_WAVEFORM_MAX_SAMPLES  4096

/**
 * @brief Single Gaussian component representing one laser return
 *
 * Model:  f(t) = amplitude * exp( -(t - center)^2 / (2 * sigma^2) )
 *
 *   amplitude  = peak amplitude (same units as raw waveform)
 *   center     = temporal position of peak [sample index or ns]
 *   sigma      = standard deviation [samples or ns]
 *   fwhm       = 2 * sqrt(2*ln(2)) * sigma ≈ 2.35482 * sigma
 */
typedef struct {
    double amplitude;   /**< Peak amplitude */
    double center;      /**< Peak center position [ns] */
    double sigma;       /**< Standard deviation [ns] */
    double fwhm;        /**< Full Width at Half Maximum [ns] */
    double range;       /**< Range derived from center: R = c*center/(2*1000) [m] */
    double energy;      /**< Integral under Gaussian = ampl*sigma*sqrt(2*pi) */
    double r_squared;   /**< R² fit quality for this component */
} lidar_gaussian_component_t;

/**
 * @brief Digitized LiDAR waveform
 *
 * Represents the time-resolved intensity of a single laser shot.
 * Range bins are: R[i] = (t_offset + i * dt) * c / 2
 */
typedef struct {
    double  *samples;       /**< Waveform amplitude samples (raw ADC or normalized) */
    size_t   num_samples;   /**< Number of samples */
    double   dt;            /**< Sample interval (bin width) [ns] */
    double   t_offset;      /**< Time offset to first sample [ns] */
    double   noise_floor;   /**< Estimated noise floor [same units as samples] */
    double   noise_std;     /**< Estimated noise standard deviation */
    double   max_amplitude; /**< Maximum amplitude in waveform */
    size_t   max_index;     /**< Index of maximum sample */
    double   pulse_width;   /**< Expected pulse FWHM [ns], for initial estimates */
} lidar_waveform_t;

/**
 * @brief Result of Gaussian decomposition of a waveform
 */
typedef struct {
    lidar_gaussian_component_t components[LIDAR_WAVEFORM_MAX_RETURNS];
    size_t                     num_components;
    double                     residual_rms;  /**< RMS of fitting residual */
    double                     r_squared;     /**< Overall R² of fit */
    int                        converged;     /**< 1 if iterative fitting converged */
    int                        num_iterations;/**< Number of iterations used */
} lidar_waveform_decomp_t;

/* ─── L1: Waveform creation and management ─────────────────────────────── */

/**
 * @brief Initialize a waveform with given sample count
 *
 * @param wf           Waveform struct to initialize
 * @param num_samples  Number of time bins
 * @param dt           Sample interval [ns]
 * @param t_offset     Time of first sample [ns]
 * @return             0 on success, -1 on error
 */
int lidar_waveform_init(lidar_waveform_t *wf, size_t num_samples,
                         double dt, double t_offset, double pulse_width);

/**
 * @brief Free waveform memory
 */
void lidar_waveform_free(lidar_waveform_t *wf);

/**
 * @brief Generate a synthetic waveform with known Gaussian returns
 *
 * Creates a waveform as sum of N Gaussian pulses plus noise.
 * Useful for algorithm validation and testing.
 *
 * @param wf           Output waveform (must be initialized)
 * @param components   Array of Gaussian component parameters
 * @param num_comp     Number of components
 * @param noise_std    Standard deviation of additive Gaussian noise
 * @return             0 on success
 */
int lidar_waveform_synthesize(lidar_waveform_t *wf,
                               const lidar_gaussian_component_t *components,
                               size_t num_comp,
                               double noise_std);

/* ─── L5: Peak detection algorithms ─────────────────────────────────────── */

/**
 * @brief Constant Fraction Discrimination (CFD) for pulse timing
 *
 * CFD finds the time when the pulse reaches a fixed fraction (typically 50%)
 * of its peak amplitude.  More robust than leading-edge detection because
 * it is independent of pulse amplitude variations (walk error compensation).
 *
 * Algorithm:
 *   threshold = noise_floor + fraction * (peak_amplitude - noise_floor)
 *   Find crossing point via linear interpolation between samples.
 *
 * Reference: Knoll, G.F., *Radiation Detection and Measurement*, 4th ed.,
 *            Wiley, 2010, Ch.17.
 *
 * @param wf        Input waveform
 * @param fraction  Fraction of peak for threshold (e.g., 0.5 for 50%)
 * @return          Timing of crossing [ns], or -1 if not found
 */
double lidar_cfd_timing(const lidar_waveform_t *wf, double fraction);

/**
 * @brief Leading-edge detection (simple threshold crossing)
 *
 * Finds the first sample exceeding noise_floor + k * noise_std.
 * Fast but suffers from walk error (timing shifts with amplitude).
 *
 * @param wf        Input waveform
 * @param k_sigma   Threshold in units of noise standard deviation (typical: 3-5)
 * @return          Timing of first threshold crossing [ns], or -1 if not found
 */
double lidar_leading_edge_timing(const lidar_waveform_t *wf, double k_sigma);

/**
 * @brief Peak detection via 1st-derivative zero-crossing
 *
 * Finds local maxima by detecting sign changes in the first derivative.
 * Smooths with a moving average to reduce false peaks from noise.
 *
 * @param wf          Input waveform
 * @param peak_times  Output: array of detected peak times [ns]
 * @param max_peaks   Maximum number of peaks to find
 * @param min_sep     Minimum separation between peaks [ns]
 * @param smooth_win  Smoothing window half-width [samples]
 * @return            Number of peaks found
 */
int lidar_detect_peaks_derivative(const lidar_waveform_t *wf,
                                   double *peak_times,
                                   size_t max_peaks,
                                   double min_sep,
                                   size_t smooth_win);

/**
 * @brief Compute noise statistics from waveform tail (noise-only region)
 *
 * Estimates noise floor and standard deviation from the last N samples
 * of the waveform, assumed to contain no signal returns.
 *
 * @param wf          Input waveform
 * @param tail_ratio  Fraction of samples at end used for noise estimation (e.g., 0.1)
 * @return            0 on success
 */
int lidar_waveform_noise_estimate(lidar_waveform_t *wf, double tail_ratio);

/* ─── L5: Gaussian decomposition ─────────────────────────────────────────── */

/**
 * @brief Decompose waveform into Gaussian components via Levenberg-Marquardt
 *
 * Full-waveform LiDAR decomposes the digitized return signal into a sum of
 * Gaussian pulses, each corresponding to a distinct reflecting surface.
 *
 * Model:  y(t) = offset + Σ A_i · exp( -(t - μ_i)² / (2·σ_i²) )
 *
 * Algorithm:
 *   1. Initial peak detection (derivative zero-crossing)
 *   2. Initial parameter estimation for each peak
 *   3. Levenberg-Marquardt non-linear least squares fitting
 *       min Σ (y_measured - y_model)²
 *   4. Removal of spurious components (amplitude < threshold)
 *
 * Complexity: O(N · M · I) where N = samples, M = components, I = iterations
 *
 * Reference: Wagner et al. (2006); Marquardt, D.W., "An Algorithm for
 *            Least-Squares Estimation of Nonlinear Parameters",
 *            *J. SIAM* 11(2), pp.431-441, 1963.
 *
 * @param wf      Input waveform
 * @param result  Output decomposition result
 * @param max_iter Maximum LM iterations
 * @param tol     Convergence tolerance on parameter change
 * @return        0 on success, -1 on error
 */
int lidar_gaussian_decompose(const lidar_waveform_t *wf,
                              lidar_waveform_decomp_t *result,
                              int max_iter,
                              double tol);

/**
 * @brief Levenberg-Marquardt step for Gaussian model fitting
 *
 * Core nonlinear optimization for the waveform decomposition.
 *
 * @param wf        Waveform data
 * @param params    Parameter vector: {A_1, mu_1, sigma_1, ..., A_n, mu_n, sigma_n, offset}
 * @param num_comp  Number of Gaussian components
 * @param max_iter  Maximum iterations
 * @param lambda    Initial damping parameter
 * @param tol       Convergence tolerance
 * @return          Number of iterations performed, -1 on failure
 */
int lidar_levenberg_marquardt_gaussian(const lidar_waveform_t *wf,
                                        double *params,
                                        size_t num_comp,
                                        int max_iter,
                                        double lambda,
                                        double tol);

/* ─── L5: Pulse characterization ────────────────────────────────────────── */

/**
 * @brief Compute Full Width at Half Maximum (FWHM) of a pulse
 *
 * FWHM is measured as the width of the pulse at 50% of its peak amplitude
 * above the noise floor.
 *
 * @param wf      Input waveform
 * @param center  Center index around which to measure FWHM
 * @return        FWHM [ns], or -1 if cannot be determined
 */
double lidar_pulse_fwhm(const lidar_waveform_t *wf, size_t center);

/**
 * @brief Compute pulse energy (integral) over a window
 *
 * Energy = Σ samples[i] * dt  over [center - half_win, center + half_win]
 * with noise floor subtracted.
 *
 * @param wf        Input waveform
 * @param center    Center sample index
 * @param half_win  Half-window width [samples]
 * @return          Integrated energy (noise-subtracted)
 */
double lidar_pulse_energy(const lidar_waveform_t *wf,
                           size_t center, size_t half_win);

/**
 * @brief Compute pulse skewness (asymmetry measure)
 *
 * Skewness > 0: tail extends to right (longer fall time)
 * Skewness < 0: tail extends to left  (longer rise time)
 * Skewness = 0: symmetric pulse
 *
 * @param wf        Input waveform
 * @param center    Center sample index
 * @param half_win  Half-window width
 * @return          Skewness (dimensionless)
 */
double lidar_pulse_skewness(const lidar_waveform_t *wf,
                             size_t center, size_t half_win);

/**
 * @brief Evaluate Gaussian model at time t
 *
 * f(t) = amplitude * exp( -(t - center)^2 / (2 * sigma^2) )
 *
 * @param comp   Gaussian component parameters
 * @param t      Evaluation time [ns]
 * @return       Model value at time t
 */
double lidar_gaussian_eval(const lidar_gaussian_component_t *comp, double t);

/**
 * @brief Evaluate multi-Gaussian model (sum of components + offset)
 *
 * @param components  Array of Gaussian components
 * @param num_comp    Number of components
 * @param offset      Baseline offset
 * @param t           Evaluation time [ns]
 * @return            Sum of all components at time t
 */
double lidar_multigaussian_eval(const lidar_gaussian_component_t *components,
                                 size_t num_comp, double offset, double t);

#endif /* LIDAR_WAVEFORM_H */