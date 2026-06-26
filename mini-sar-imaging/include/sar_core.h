/**
 * @file    sar_core.h
 * @brief   SAR Core Definitions -- L1 Definitions + L3 Mathematical Structures
 *
 * @details Fundamental SAR imaging data structures: radar waveform parameters,
 *          chirp signals, SAR complex image, raw echo data model, and core
 *          mathematical structures for SAR signal processing in 2D Fourier domain.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Range resolution: rho_r = c/(2B) -- governed by bandwidth B
 *     - Azimuth resolution: rho_a = L_a/2 (focused SAR) -- antenna length dependent
 *     - Chirp waveform / Linear FM: s(t) = exp(j * pi * K_r * t^2)
 *     - Pulse Repetition Frequency (PRF)
 *     - Range cell migration (RCM)
 *     - Complex SAR image: magnitude (backscatter) + phase (range information)
 *     - Backscatter coefficient sigma_0, Normalized Radar Cross Section (NRCS)
 *     - Integration time / synthetic aperture length
 *     - Slant range vs ground range geometry
 *     - Doppler centroid f_Dc, Doppler rate f_R
 *   L2 - Core Concepts:
 *     - Pulse compression via matched filtering
 *     - Synthetic aperture formation
 *     - Coherent integration principle
 *     - Range-Doppler coupling in chirp signals
 *   L3 - Mathematical Structures:
 *     - 2D Fourier transform of raw SAR data
 *     - Range-Doppler domain representation
 *     - Wavenumber (omega-k / k_x-k_r) domain representation
 *     - Stationary phase principle for SAR signal analysis
 *
 * Reference:
 *   - Cumming & Wong, "Digital Processing of SAR Data" (2005)
 *   - Curlander & McDonough, "Synthetic Aperture Radar: Systems & Signal Processing" (1991)
 *   - Richards, Scheer & Holm, "Principles of Modern Radar" (2010)
 *   - Chen & Martorella, "Inverse Synthetic Aperture Radar Imaging" (2014)
 */

#ifndef SAR_CORE_H
#define SAR_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --- Physical Constants --- */

/** Speed of light in vacuum [m/s] */
#define SAR_C         2.99792458e8
/** Speed of light / 2 -- half-light factor appears throughout SAR equations */
#define SAR_C_OVER_2  1.49896229e8

/* --- L1: SAR Operating Mode Enumeration --- */

/**
 * @brief SAR imaging mode classification.
 *
 * Stripmap: fixed antenna pointing, continuous along-track strip
 * Spotlight: steering antenna to dwell on a spot, finer azimuth resolution
 * ScanSAR:   scanning across multiple sub-swaths, wider coverage
 * TOPS:      Terrain Observation by Progressive Scan, rotating antenna
 * Inverse:   ISAR, radar stationary, target moving
 */
typedef enum {
    SAR_MODE_STRIPMAP = 0,
    SAR_MODE_SPOTLIGHT = 1,
    SAR_MODE_SCANSAR   = 2,
    SAR_MODE_TOPS      = 3,
    SAR_MODE_INVERSE   = 4,
    SAR_MODE_BISTATIC  = 5
} sar_mode_t;

/**
 * @brief Radar band designation -- IEEE standard letter bands.
 */
typedef enum {
    SAR_BAND_P  = 0,   /**< ~0.3 GHz (foliage penetration) */
    SAR_BAND_L  = 1,   /**< ~1.25 GHz */
    SAR_BAND_S  = 2,   /**< ~3 GHz */
    SAR_BAND_C  = 3,   /**< ~5.3 GHz (ERS, Radarsat) */
    SAR_BAND_X  = 4,   /**< ~9.65 GHz (TerraSAR-X) */
    SAR_BAND_Ku = 5,   /**< ~15 GHz */
    SAR_BAND_K  = 6,   /**< ~24 GHz */
    SAR_BAND_Ka = 7,   /**< ~35 GHz */
    SAR_BAND_W  = 8    /**< ~94 GHz */
} sar_band_t;

/* --- L1: Radar System Parameters --- */

/**
 * @brief SAR system parameter structure.
 *
 * Core radar parameters that define the imaging system configuration.
 *
 * Range resolution (slant):  rho_r = c / (2 * B_r)
 * Azimuth resolution (focused stripmap): rho_a = L_a / 2
 * PRF constraint (Nyquist in azimuth):  PRF >= 2 * v_plat / L_a
 * Maximum unambiguous range:           R_max = c / (2 * PRF)
 */
typedef struct {
    double   carrier_freq_hz;     /**< f0, RF carrier frequency [Hz] */
    double   wavelength_m;        /**< lambda = c/f0 [m] */
    double   bandwidth_hz;        /**< B_r, Range bandwidth [Hz] */
    double   pulse_width_s;       /**< tau_p, Transmitted pulse duration [s] */
    double   prf_hz;              /**< PRF, Pulse Repetition Frequency [Hz] */
    double   chirp_rate;          /**< K_r = B_r/tau_p, Chirp rate [Hz/s] */
    double   sample_rate_hz;      /**< f_s, range sampling rate [Hz] */
    double   range_resolution_m;  /**< rho_r = c/(2*B_r), slant range resolution [m] */
    double   azimuth_resolution_m;/**< rho_a, azimuth resolution [m] */
    double   antenna_length_m;    /**< L_a, antenna length along-track [m] */
    double   platform_velocity_ms;/**< v, platform velocity [m/s] */
    double   squint_angle_rad;    /**< theta_sq, squint angle [rad] (0 = broadside) */
    double   look_angle_rad;      /**< theta_l, look angle from nadir [rad] */
    double   incidence_angle_rad; /**< theta_i, incidence angle at surface [rad] */
    double   platform_altitude_m; /**< H, platform altitude above surface [m] */
    double   near_range_m;        /**< R_n, slant range to near edge of swath [m] */
    double   far_range_m;         /**< R_f, slant range to far edge of swath [m] */
    double   swath_width_m;       /**< Ground range swath width [m] */
    sar_mode_t mode;              /**< Imaging mode */
    sar_band_t band;              /**< Frequency band */
    char     sensor_name[64];     /**< Sensor/platform identifier string */
} sar_params_t;

/**
 * @brief Set SAR parameters from frequency and compute derived quantities.
 *
 * Given carrier frequency f0, bandwidth B_r, pulse width tau_p, PRF, antenna
 * length L_a, platform velocity v, altitude H, and geometry angles, this
 * function fills in all derived fields of sar_params_t:
 *   wavelength = c/f0, chirp_rate = B_r/tau_p,
 *   range_resolution = c/(2*B_r), azimuth_resolution = L_a/2
 *
 * @param p     Pointer to sar_params_t (modified in-place)
 * @param freq  Carrier frequency [Hz]
 * @param bw    Bandwidth [Hz]
 * @param pw    Pulse width [s]
 * @param prf   PRF [Hz]
 * @param ant   Antenna length [m]
 * @param vel   Platform velocity [m/s]
 * @param alt   Platform altitude [m]
 * @param squint Squint angle [rad]
 * @param look  Look angle [rad]
 * @param near  Near slant range [m]
 * @param far   Far slant range [m]
 */
void sar_params_init(sar_params_t *p,
                     double freq, double bw, double pw, double prf,
                     double ant, double vel, double alt,
                     double squint, double look,
                     double near, double far);

/* --- L1: Chirp / Linear FM Waveform --- */

/**
 * @brief Complex baseband chirp (Linear FM) signal definition.
 *
 * @math  s(t) = exp( j * pi * K_r * t^2 )  for  -tau/2 <= t <= tau/2
 *        K_r = B / tau   (chirp rate in Hz/s)
 *
 * The chirp is the fundamental transmitted waveform in SAR.
 * Time-bandwidth product: TBP = B * tau = K_r * tau^2
 */
typedef struct {
    size_t   num_samples;         /**< N, number of samples */
    double   pulse_width_s;       /**< tau, pulse duration [s] */
    double   bandwidth_hz;        /**< B, bandwidth [Hz] */
    double   chirp_rate;          /**< K = B/tau, chirp rate [Hz/s] */
    double   sample_rate_hz;      /**< f_s, sampling rate [Hz] */
    double   dt;                  /**< Delta_t = 1/f_s, sampling interval [s] */
    double   t_min;               /**< -tau/2, start time [s] */
    double   t_max;               /**< +tau/2, end time [s] */
    double   time_bandwidth_product; /**< TBP = B*tau */
    double  *I;                   /**< In-phase samples [N] */
    double  *Q;                   /**< Quadrature samples [N] */
} sar_chirp_t;

/**
 * @brief Allocate and generate a complex baseband chirp.
 *
 * s[n] = exp( j * pi * K * (n*dt)^2 )  for n = 0, ..., N-1
 * I[n] = cos( pi * K * (n*dt)^2 )
 * Q[n] = sin( pi * K * (n*dt)^2 )
 *
 * @param N        Number of samples
 * @param tau      Pulse width [s]
 * @param B        Bandwidth [Hz]
 * @param fs       Sample rate [Hz]
 * @return         Allocated sar_chirp_t (caller must free with sar_chirp_free)
 */
sar_chirp_t *sar_chirp_alloc(size_t N, double tau, double B, double fs);

/**
 * @brief Free a sar_chirp_t structure.
 */
void sar_chirp_free(sar_chirp_t *c);

/**
 * @brief Compute the autocorrelation (matched filter response) of the chirp.
 *
 * Autocorrelation: R(tau) = integral s(t) * conj(s(t-tau)) dt
 * For a chirp: peak width ~ 1/B, peak-to-sidelobe ratio ~ 13.3 dB
 * (rectangular window gives 13.26 dB theoretical PSLR)
 *
 * @param c      Chirp signal
 * @param R_out  Output autocorrelation [2*N-1], caller pre-allocates
 */
void sar_chirp_autocorrelation(const sar_chirp_t *c, double *R_out);

/* --- L1: Matched Filter for Pulse Compression --- */

/**
 * @brief Generate the matched filter impulse response for a chirp.
 *
 * Matched filter maximizes SNR at output: h(t) = conj(s(-t))
 * For the baseband chirp s(t) = exp(j*pi*K*t^2)
 * The matched filter: h(t) = exp(-j*pi*K*t^2)
 *
 * @param c          Input chirp
 * @param h_I_out    Output I component of matched filter [N]
 * @param h_Q_out    Output Q component of matched filter [N]
 */
void sar_matched_filter_coeff(const sar_chirp_t *c,
                              double *h_I_out, double *h_Q_out);

/**
 * @brief Apply pulse compression (matched filtering) to raw range data.
 *
 * Time-domain convolution of input with matched filter:
 *   y[n] = sum_k x[k] * h[n-k]  (for all n)
 *
 * R, x, h, y: all double arrays
 *
 * @param x_I, x_Q  Input complex signal (I/Q) [Nx]
 * @param Nx        Input length
 * @param h_I, h_Q  Matched filter coeffs [Nh]
 * @param Nh        Filter length
 * @param y_I, y_Q  Output compressed signal [Nx+Nh-1], caller pre-allocates
 */
void sar_pulse_compression(const double *x_I, const double *x_Q, size_t Nx,
                           const double *h_I, const double *h_Q, size_t Nh,
                           double *y_I, double *y_Q);

/**
 * @brief Frequency-domain pulse compression via FFT (faster O(N log N)).
 *
 * Steps:
 *   1. FFT of x[n] -> X[k]
 *   2. Multiply by conj(FFT of h[n]) = H*[k]
 *   3. IFFT of product
 *
 * @param x_I, x_Q  Input complex signal [N]
 * @param N         Signal length (power of 2 recommended)
 * @param h_I, h_Q  Matched filter coeffs [M]
 * @param M         Filter length
 * @param y_I, y_Q  Output compressed signal [N], caller pre-allocates
 */
void sar_pulse_compression_fft(const double *x_I, const double *x_Q, size_t N,
                               const double *h_I, const double *h_Q, size_t M,
                               double *y_I, double *y_Q);

/* --- L1: SAR Raw Echo Data Model --- */

/**
 * @brief SAR raw echo data (before processing).
 *
 * The 2D raw data represents the output of the radar receiver,
 * organized as range lines x range bins:
 *   Range direction (fast time)  -> index r (columns)
 *   Azimuth direction (slow time) -> index a (rows)
 *
 * For a point target at (R0, eta_0):
 *   s_raw(tau, eta) = A0 * w_r(tau - 2R(eta)/c) * w_a(eta - eta_0)
 *                   * exp(-j * 4*pi * f0 * R(eta)/c)
 *                   * exp(j * pi * K_r * (tau - 2R(eta)/c)^2)
 *
 * where R(eta) = sqrt(R0^2 + v^2*(eta-eta_0)^2)
 *             ~= R0 + v^2*(eta-eta_0)^2/(2*R0)  (parabolic approx)
 */
typedef struct {
    double   **data_I;           /**< In-phase raw data [naz][nrng] */
    double   **data_Q;           /**< Quadrature raw data [naz][nrng] */
    size_t     naz;              /**< Number of azimuth lines */
    size_t     nrng;             /**< Number of range bins */
    sar_params_t params;         /**< Radar parameters */
    double     range_sampling_interval; /**< delta_r = c/(2*f_s), range bin spacing [m] */
    double     azimuth_sampling_interval; /**< delta_eta = 1/PRF, pulse spacing [s] */
    double     doppler_centroid_hz;  /**< f_Dc, Doppler centroid [Hz] */
    double     doppler_rate_hz_per_s;/**< f_R, Doppler rate (FM rate) [Hz/s] */
} sar_raw_data_t;

/**
 * @brief Allocate a sar_raw_data_t structure for given dimensions.
 *
 * @param naz   Number of azimuth lines
 * @param nrng  Number of range bins
 * @return      Allocated structure (caller frees with sar_raw_data_free)
 */
sar_raw_data_t *sar_raw_data_alloc(size_t naz, size_t nrng);

/**
 * @brief Free a sar_raw_data_t structure.
 */
void sar_raw_data_free(sar_raw_data_t *raw);

/**
 * @brief Fill raw data with simulated point target echo using exact range equation.
 *
 * Simulates a single point target at (range_idx, az_idx) using
 * the exact range equation R(eta) = sqrt(R0^2 + v^2*(eta-eta_0)^2).
 *
 * @param raw       Raw data to fill (I/Q accumulation)
 * @param rng_idx   Range index of target
 * @param az_idx    Azimuth index of target
 * @param amplitude Backscatter amplitude
 */
void sar_raw_data_point_target(sar_raw_data_t *raw,
                               size_t rng_idx, size_t az_idx,
                               double amplitude);

/**
 * @brief Add additive white Gaussian noise to the raw data.
 *
 * @param raw    Raw data (modified in-place)
 * @param sigma  Noise standard deviation
 */
void sar_raw_data_add_noise(sar_raw_data_t *raw, double sigma);

/* --- L1: SAR Complex Image --- */

/**
 * @brief Single-Look Complex (SLC) SAR image.
 *
 * An SLC image is the focused output of SAR processing.
 * Each pixel is a complex number (I + jQ):
 *   - Magnitude: |z| = sqrt(I^2 + Q^2) -> backscatter intensity
 *   - Phase:     phi = atan2(Q, I)     -> range + scatterer phase
 *
 * The complex nature enables interferometry (InSAR) where phase
 * differences between two acquisitions give range change.
 */
typedef struct {
    double  **data_I;            /**< Real part (in-phase) [nrows][ncols] */
    double  **data_Q;            /**< Imaginary part (quadrature) [nrows][ncols] */
    size_t    nrows;             /**< Number of rows (azimuth pixels) */
    size_t    ncols;             /**< Number of columns (range pixels) */
    double    range_pixel_spacing_m;  /**< Range pixel spacing [m] */
    double    azimuth_pixel_spacing_m;/**< Azimuth pixel spacing [m] */
    char      sensor_name[64];        /**< Sensor/platform identifier */
} sar_image_t;

/**
 * @brief Allocate a sar_image_t with given dimensions.
 */
sar_image_t *sar_image_alloc(size_t nrows, size_t ncols);

/**
 * @brief Free a sar_image_t.
 */
void sar_image_free(sar_image_t *img);

/**
 * @brief Compute magnitude image from complex data.
 * @math  M[i][j] = sqrt(I[i][j]^2 + Q[i][j]^2)
 * @param img     Input complex SLC image
 * @param mag_out Output magnitude image [nrows*ncols], caller pre-allocates
 */
void sar_image_magnitude(const sar_image_t *img, double *mag_out);

/**
 * @brief Compute phase image from complex data.
 * @math phi[i][j] = atan2(Q[i][j], I[i][j])  in [-pi, pi]
 * @param img       Input complex SLC image
 * @param phase_out Output phase image [nrows*ncols], caller pre-allocates
 */
void sar_image_phase(const sar_image_t *img, double *phase_out);

/**
 * @brief Multi-look averaging to reduce speckle noise.
 *
 * Multi-looking divides the azimuth bandwidth into N_look sub-bands,
 * processes each independently, then incoherently averages:
 *   M_ML = (1/N_l) * sum_k |SLC_k|^2
 *
 * Equivalent to spatial averaging with N_l degrees of freedom.
 * Speckle reduction factor = 1/sqrt(N_l)
 *
 * @param img      Input SLC
 * @param nlooks   Number of looks
 * @param ml_out   Output multi-looked magnitude [nrows*ncols]
 */
void sar_multilook(const sar_image_t *img, int nlooks, double *ml_out);

/**
 * @brief Radiometric calibration: convert DN to sigma_0 (NRCS).
 *
 * sigma_0[dB] = 10*log10( DN^2 * K_cal )
 *
 * where K_cal is the calibration constant from corner reflectors.
 *
 * @param mag       Input magnitude (Digital Number)
 * @param n         Number of pixels
 * @param cal_const Calibration constant K_cal
 * @param sigma0_db Output sigma_0 (Normalized Radar Cross Section) in dB [n]
 */
void sar_calibrate_sigma0(const double *mag, size_t n,
                          double cal_const, double *sigma0_db);

/* --- L3: 2D Fourier Domain Representations --- */

/**
 * @brief Compute the 2D FFT of the raw data for frequency-domain algorithms.
 *
 * The range-Doppler algorithm operates on data transformed into
 * range-frequency / azimuth-frequency domain (f_tau, f_eta):
 *
 *   S_2df(f_tau, f_eta) = double_integral s(tau, eta)
 *                        * exp(-j 2*pi (f_tau*tau + f_eta*eta)) dtau deta
 *
 * @param data_I, data_Q  Complex 2D input [nrows][ncols]
 * @param nrows           Azimuth dimension (slow time)
 * @param ncols           Range dimension (fast time)
 * @param spec_I, spec_Q  Complex 2D output [nrows][ncols], caller pre-allocates
 * @param forward         Non-zero = forward FFT, zero = inverse FFT
 */
void sar_fft2d(const double **data_I, const double **data_Q,
               size_t nrows, size_t ncols,
               double **spec_I, double **spec_Q, int forward);

/**
 * @brief Compute a 1D FFT along range (rows stay, FFT on columns).
 *
 * Transforms each azimuth line from time-domain to range-frequency domain.
 * Used in the first step of RDA: range compression in frequency domain.
 *
 * @param data_I, data_Q  Complex 2D input [nrows][ncols]
 * @param nrows           Number of rows
 * @param ncols           Number of columns (FFT by)
 * @param forward         Non-zero = forward, zero = inverse
 * @param out_I, out_Q    Complex 2D output [nrows][ncols], caller pre-allocates
 */
void sar_fft_range(const double **data_I, const double **data_Q,
                   size_t nrows, size_t ncols,
                   double **out_I, double **out_Q, int forward);

/**
 * @brief Compute a 1D FFT along azimuth (columns stay, FFT on rows).
 *
 * Used after range compression to bring data into range-Doppler domain.
 *
 * @param data_I, data_Q  Complex 2D input [nrows][ncols]
 * @param nrows           Number of rows (FFT by)
 * @param ncols           Number of columns
 * @param forward         Non-zero = forward, zero = inverse
 * @param out_I, out_Q    Complex 2D output [nrows][ncols], caller pre-allocates
 */
void sar_fft_azimuth(const double **data_I, const double **data_Q,
                     size_t nrows, size_t ncols,
                     double **out_I, double **out_Q, int forward);

/* --- L3: Wavenumber Domain (omega-k) Mapping --- */

/**
 * @brief Wavenumber domain parameters.
 *
 * In omega-k algorithm, data is processed in the 2D frequency
 * domain (f_tau, f_eta). The fundamental relationship:
 *
 *   k_r = sqrt( (4*pi*(f0+f_tau)/c)^2 - k_x^2 )
 *
 * where k_x = 2*pi*f_eta/v  is the azimuth wavenumber.
 *
 * Stolt interpolation maps from (f_tau, f_eta) -> (k_r, k_x).
 */
typedef struct {
    double   f_tau_min;          /**< Minimum range frequency [Hz] */
    double   f_tau_max;          /**< Maximum range frequency [Hz] */
    double   f_eta_min;          /**< Minimum azimuth frequency [Hz] */
    double   f_eta_max;          /**< Maximum azimuth frequency [Hz] */
    double   k_r_min;            /**< Minimum range wavenumber [rad/m] */
    double   k_r_max;            /**< Maximum range wavenumber [rad/m] */
    double   k_x_min;            /**< Minimum azimuth wavenumber [rad/m] */
    double   k_x_max;            /**< Maximum azimuth wavenumber [rad/m] */
    double   f0;                 /**< Carrier frequency [Hz] */
    double   v;                  /**< Platform velocity [m/s] */
    double   c;                  /**< Speed of light [m/s] */
} sar_wavenumber_params_t;

/**
 * @brief Initialize wavenumber parameters from SAR system parameters.
 *
 * Computes the mapping between temporal frequencies and wavenumbers:
 *   k_r(f_tau) = 4*pi*(f0 + f_tau)/c
 *   k_x(f_eta) = 2*pi*f_eta/v
 *
 * @param wp     Output wavenumber params
 * @param sp     Input SAR system params
 */
void sar_wavenumber_init(sar_wavenumber_params_t *wp, const sar_params_t *sp);

/**
 * @brief Stolt interpolation: map (f_tau, f_eta) -> k_r for omega-k algorithm.
 *
 * For given azimuth wavenumber k_x, the range wavenumber:
 *   k_r(k_x, f_tau) = sqrt( (4*pi*(f0+f_tau)/c)^2 - k_x^2 )
 *
 * Stolt interpolation resamples data from uniform f_tau grid to uniform
 * k_r grid, correcting range-dependent range cell migration.
 *
 * @param wp          Wavenumber parameters
 * @param f_eta       Azimuth frequency for this column [Hz]
 * @param f_tau_in    Input range frequency grid [Nr]
 * @param Nr          Number of range frequency bins
 * @param k_r_out     Output range wavenumber grid [Nr] (uniform spacing)
 */
void sar_stolt_interpolation(const sar_wavenumber_params_t *wp,
                             double f_eta,
                             const double *f_tau_in, size_t Nr,
                             double *k_r_out);

#endif /* SAR_CORE_H */