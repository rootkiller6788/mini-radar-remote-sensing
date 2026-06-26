#ifndef PHASED_ARRAY_H
#define PHASED_ARRAY_H

#include <stddef.h>
#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions — Phased Array Geometry, Element, and Parameter Types
 *
 * Reference: Balanis (2016) Chapter 6 — Array Theory
 *            Richards, Scheer, Holm (2010) Principles of Modern Radar Ch. 9-10
 *            Mailloux (2005) Phased Array Antenna Handbook Ch. 1-2
 * ============================================================================ */

/* --- Array geometry topology --- */
typedef enum {
    PA_GEOMETRY_LINEAR           = 0,   /* 1D uniform/non-uniform linear array */
    PA_GEOMETRY_PLANAR_RECT      = 1,   /* 2D rectangular grid planar array */
    PA_GEOMETRY_PLANAR_TRIANG    = 2,   /* 2D triangular lattice planar array */
    PA_GEOMETRY_CIRCULAR         = 3,   /* Circular ring array */
    PA_GEOMETRY_CYLINDRICAL      = 4,   /* Cylindrical conformal array */
    PA_GEOMETRY_SPHERICAL        = 5,   /* Spherical conformal array */
    PA_GEOMETRY_HEXAGONAL        = 6    /* Hexagonal grid planar array */
} pa_geometry_t;

/* --- Element type — antenna element used in the array --- */
typedef enum {
    PA_ELEMENT_ISOTROPIC         = 0,   /* Ideal isotropic radiator (reference) */
    PA_ELEMENT_HALF_WAVE_DIPOLE  = 1,   /* λ/2 dipole */
    PA_ELEMENT_PATCH_RECT        = 2,   /* Rectangular microstrip patch */
    PA_ELEMENT_PATCH_CIRCULAR    = 3,   /* Circular microstrip patch */
    PA_ELEMENT_SLOT              = 4,   /* Slot antenna */
    PA_ELEMENT_VIVALDI           = 5,   /* Tapered slot / Vivaldi */
    PA_ELEMENT_DIPOLE_CROSSED    = 6,   /* Crossed dipole (dual-pol) */
    PA_ELEMENT_PATCH_DUAL_FEED   = 7    /* Dual-feed patch (circular pol) */
} pa_element_type_t;

/* --- Amplitude taper window functions --- */
typedef enum {
    PA_WINDOW_UNIFORM            = 0,   /* Uniform weighting (no taper) */
    PA_WINDOW_HAMMING            = 1,   /* Hamming window */
    PA_WINDOW_HANN               = 2,   /* Hann (Hanning) window */
    PA_WINDOW_BLACKMAN           = 3,   /* Blackman window */
    PA_WINDOW_DOLPH_CHEBYSHEV    = 4,   /* Dolph-Chebyshev (constant SLL) */
    PA_WINDOW_TAYLOR             = 5,   /* Taylor distribution */
    PA_WINDOW_BINOMIAL           = 6,   /* Binomial distribution (no sidelobes) */
    PA_WINDOW_KAISER             = 7,   /* Kaiser-Bessel window */
    PA_WINDOW_GAUSSIAN           = 8,   /* Gaussian taper */
    PA_WINDOW_RAISED_COSINE      = 9,   /* Raised cosine taper */
    PA_WINDOW_CUSTOM             = 100  /* User-supplied coefficients */
} pa_window_t;

/* --- Beamforming architecture --- */
typedef enum {
    PA_BF_ANALOG_RF              = 0,   /* Phase shifters at RF */
    PA_BF_ANALOG_IF              = 1,   /* Phase shifters at IF */
    PA_BF_DIGITAL_ELEMENT        = 2,   /* Digital beamforming per element */
    PA_BF_DIGITAL_SUBARRAY       = 3,   /* Hybrid beamforming */
    PA_BF_OPTICAL_TTD            = 4    /* True-time-delay beamforming */
} pa_beamforming_arch_t;

/* --- Steering mechanism --- */
typedef enum {
    PA_STEERING_PHASE            = 0,   /* Phase shift steering */
    PA_STEERING_TRUE_TIME_DELAY  = 1,   /* True-time-delay steering */
    PA_STEERING_FREQUENCY        = 2    /* Frequency scanned array */
} pa_steering_t;

/* --- Bit resolution of phase shifter --- */
typedef enum {
    PA_PHASE_BITS_1              = 1,
    PA_PHASE_BITS_2              = 2,
    PA_PHASE_BITS_3              = 3,
    PA_PHASE_BITS_4              = 4,
    PA_PHASE_BITS_5              = 5,
    PA_PHASE_BITS_6              = 6,
    PA_PHASE_BITS_8              = 8,
    PA_PHASE_BITS_CONTINUOUS     = 64   /* Ideal continuous phase */
} pa_phase_bits_t;

/* --- Monopulse mode --- */
typedef enum {
    PA_MONOPULSE_AMPLITUDE       = 0,   /* Amplitude comparison monopulse */
    PA_MONOPULSE_PHASE           = 1,   /* Phase comparison monopulse */
    PA_MONOPULSE_SUM_DIFF        = 2    /* Full Σ/Δ monopulse */
} pa_monopulse_mode_t;

/* ============================================================================
 * L1: Coordinate and Position Structures
 * ============================================================================ */

/** 3D position vector (meters, element placement) */
typedef struct {
    double x, y, z;
} pa_position3d_t;

/** Spherical coordinate (radians) — theta: polar angle from z-axis,
 *  phi: azimuthal angle from x-axis in xy-plane */
typedef struct {
    double theta;    /* 0 .. π */
    double phi;      /* 0 .. 2π */
} pa_spherical_t;

/** Azimuth/Elevation pair for display convention (degrees, radar coords) */
typedef struct {
    double azimuth_deg;    /* Horizontal scan angle, ± range */
    double elevation_deg;  /* Vertical scan angle, ± range */
} pa_azel_t;

/** Single array element description */
typedef struct {
    pa_position3d_t position;                    /* (x,y,z) location (meters) */
    pa_element_type_t element_type;              /* Element type */
    double complex excitation;                   /* Complex weight: A*exp(j*phase) */
    double orientation_theta;                    /* Element rotation θ (rad) */
    double orientation_phi;                      /* Element rotation φ (rad) */
    int tr_module_id;                            /* T/R module ID (-1 if passive) */
    double tr_module_gain_db;                    /* T/R module gain (dB) */
    double tr_module_noise_figure_db;            /* T/R module noise figure (dB) */
    double tr_module_power_watt;                 /* T/R module output power (W) */
} pa_element_t;

/** Array configuration — full geometry + feed parameters */
typedef struct {
    /* Geometry */
    pa_geometry_t geometry;                      /* Array geometry type */
    uint32_t num_elements;                       /* Total number of elements */
    double frequency_hz;                         /* Operating frequency (Hz) */
    double element_spacing_x;                    /* dx spacing (m) */
    double element_spacing_y;                    /* dy spacing (m) */
    double element_spacing_radius;               /* For circular: ring radius (m) */
    /* Wiring topology */
    uint32_t rows;                               /* Number of rows (planar) */
    uint32_t cols;                               /* Number of columns (planar) */
    /* Scan limits */
    double scan_limit_az_deg;                    /* Max azimuth scan angle */
    double scan_limit_el_deg;                    /* Max elevation scan angle */
} pa_array_config_t;

/* ============================================================================
 * L1: Array State — Runtime beamforming state
 * ============================================================================ */

/** Beam steering state for a single beam */
typedef struct {
    double theta_steer_rad;                      /* Steered polar angle */
    double phi_steer_rad;                        /* Steered azimuth angle */
    double complex *weights;                     /* Per-element complex weights */
    double *phase_shifts_rad;                    /* Per-element phase shifts */
    double *amplitude_taper;                     /* Per-element amplitude taper */
    pa_window_t window_type;                     /* Amplitude window type */
    double window_param;                         /* Window param (e.g., SLL for Chebyshev) */
    pa_steering_t steering_type;                 /* Phase vs TTD steering */
    pa_phase_bits_t phase_quant_bits;            /* Phase shifter resolution */
} pa_steering_state_t;

/** Array factor evaluation result at a single angle */
typedef struct {
    double theta_rad;                            /* Angle evaluated */
    double phi_rad;                              /* Azimuth evaluated */
    double af_magnitude;                         /* |AF(θ,φ)| */
    double af_magnitude_db;                      /* |AF(θ,φ)| in dB */
    double af_phase_rad;                         /* ∠AF(θ,φ) */
    double complex af_complex;                   /* AF(θ,φ) as complex number */
    double directivity_dbi;                      /* Directivity (dBi) at this angle */
} pa_af_result_t;

/** Full beam pattern scan over angle range */
typedef struct {
    uint32_t num_theta;                          /* Number of θ samples */
    uint32_t num_phi;                            /* Number of φ samples */
    double theta_start;                          /* Start of θ range (rad) */
    double theta_end;                            /* End of θ range (rad) */
    double phi_start;                            /* Start of φ range (rad) */
    double phi_end;                              /* End of φ range (rad) */
    double *af_db;                               /* 2D pattern: af_db[θ_idx][φ_idx] */
    double max_directivity_dbi;                  /* Peak directivity */
    double sidelobe_level_db;                    /* Worst sidelobe level (dB relative to peak) */
    double half_power_beamwidth_deg;             /* 3 dB beamwidth (degrees) */
    double first_null_beamwidth_deg;             /* FNBW (degrees) */
} pa_pattern_t;

/** Grating lobe location descriptor */
typedef struct {
    double theta_grating_rad;                    /* Grating lobe polar angle */
    double phi_grating_rad;                      /* Grating lobe azimuth */
    double relative_level_db;                    /* Level relative to main lobe */
    int order_m;                                 /* Grating lobe order index */
} pa_grating_lobe_t;

/* ============================================================================
 * L1: T/R Module (AESA Building Block)
 * ============================================================================ */

/** Transmit/Receive module for AESA */
typedef struct {
    int module_id;                               /* Unique T/R module ID */
    double tx_power_watt;                        /* Peak transmit power (W) */
    double tx_gain_db;                           /* Transmit gain (dB) */
    double rx_gain_db;                           /* Receive gain (dB) */
    double noise_figure_db;                      /* Noise figure (dB) */
    double phase_resolution_deg;                 /* Phase shifter step (deg) */
    pa_phase_bits_t phase_bits;                  /* Number of phase bits */
    double attenuator_range_db;                  /* Attenuator dynamic range (dB) */
    double attenuator_step_db;                   /* Attenuator LSB (dB) */
    double switching_time_ns;                    /* T/R switch time (ns) */
    double bandwidth_mhz;                        /* RF bandwidth (MHz) */
    int is_healthy;                              /* Built-in test flag */
    double temperature_degc;                     /* Junction temperature (°C) */
} pa_tr_module_t;

/* ============================================================================
 * L1: Adaptive Beamforming Structures
 * ============================================================================ */

/** Sample Matrix Inversion (SMI) adaptive beamformer state */
typedef struct {
    uint32_t num_elements;                       /* Number of array elements */
    uint32_t num_snapshots;                      /* Number of temporal snapshots */
    double forgetting_factor;                    /* λ for RLS/SMI (0 < λ ≤ 1) */
    double diagonal_loading;                     /* Regularization for ill-conditioned R */
    double complex *weight_vector;               /* Current adaptive weight vector */
    double complex *steering_vector;             /* Look-direction steering vector */
    double complex *cov_matrix;                  /* [M×M] sample covariance matrix */
    double complex *cov_matrix_inv;              /* [M×M] inverse covariance */
    double output_sinr_db;                       /* Beamformer output SINR (dB) */
} pa_smi_beamformer_t;

/** MVDR (Capon) beamformer state */
typedef struct {
    pa_smi_beamformer_t base;                    /* Inherits SMI structure */
    double output_power;                         /* MVDR output power */
    double *angular_spectrum;                    /* Capon spatial spectrum */
    uint32_t num_angles;                         /* Spectrum angular bins */
} pa_mvdr_beamformer_t;

/** LCMV (Linearly Constrained Minimum Variance) state */
typedef struct {
    pa_smi_beamformer_t base;                    /* Inherits SMI structure */
    double complex *constraint_matrix;           /* [M×K] constraint matrix */
    double complex *constraint_response;         /* [K×1] desired response vector */
    uint32_t num_constraints;                    /* K: number of constraints */
} pa_lcmv_beamformer_t;

/* ============================================================================
 * L1: Wideband / Time-Delay Beamforming
 * ============================================================================ */

/** True-time-delay unit descriptor */
typedef struct {
    uint32_t element_index;                      /* Which element */
    double delay_ps;                             /* Delay in picoseconds */
    double delay_fractional_lambda;              /* Delay in fractions of λ */
    double squint_error_deg;                     /* Beam squint angle at band edge */
    int delay_line_taps;                         /* Number of delay line taps */
    double tap_spacing_ps;                       /* Tap spacing (ps) */
} pa_ttd_element_t;

/** Wideband beam squint analysis result */
typedef struct {
    double center_frequency_hz;
    double bandwidth_hz;
    double steer_angle_nominal_deg;
    double squint_at_band_edge_deg;              /* Beam squint at f_max */
    double max_allowable_bandwidth_hz;           /* BW before unacceptable squint */
    double fractional_bandwidth_limit;           /* % BW where squint = BW/2 */
} pa_beam_squint_t;

/* ============================================================================
 * L2: Core Concept Functions — Element Pattern Models
 *
 * Element pattern × Array factor = Total radiation pattern
 * (Pattern Multiplication Theorem)
 * ============================================================================ */

/**
 * @brief Compute the element radiation pattern for a given element type.
 *
 * Implements the individual antenna element pattern E(θ,φ).
 * For isotropic: E(θ,φ) = 1.
 * For half-wave dipole: E(θ) = cos(½π cos θ) / sin θ.
 * For patch: approximate E(θ) = cos(θ) for elevation, broad in azimuth.
 *
 * @param element_type Type of antenna element.
 * @param theta Polar angle (rad).
 * @param phi Azimuthal angle (rad).
 * @return double complex Element pattern value.
 *
 * Reference: Balanis (2016) §4.5–4.7, §14.2.
 */
double complex pa_element_pattern(pa_element_type_t element_type,
                                   double theta, double phi);

/**
 * @brief Compute the element radiation pattern magnitude in dBi.
 */
double pa_element_pattern_dbi(pa_element_type_t element_type,
                               double theta, double phi);

/* ============================================================================
 * L4: Pattern Multiplication Theorem — Core Array Factor
 * ============================================================================ */

/**
 * @brief Compute the array factor AF(θ,φ) for an arbitrary array geometry.
 *
 * Array Factor:
 *   AF(θ,φ) = Σ_{n=0}^{N-1} w_n exp(j k₀ (r_n · û))
 *   where û = (sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)) is the unit direction vector,
 *   r_n is the n-th element position, and w_n is the complex excitation weight.
 *   k₀ = 2π/λ is the free-space wavenumber.
 *
 * Theoretical basis:
 *   The total far-field pattern E_total(θ,φ) = E_element(θ,φ) × AF(θ,φ)
 *   This is the Pattern Multiplication Theorem (Balanis §6.3).
 *
 * @param config Array configuration.
 * @param elements Array of elements (size = config.num_elements).
 * @param weights Complex excitation weights (size = config.num_elements).
 * @param theta Polar scan angle (rad).
 * @param phi Azimuthal scan angle (rad).
 * @return pa_af_result_t Array factor magnitude/phase.
 *
 * O(N) per evaluation point.
 */
pa_af_result_t pa_array_factor(const pa_array_config_t *config,
                                const pa_element_t *elements,
                                const double complex *weights,
                                double theta, double phi);

/**
 * @brief Compute total radiation pattern (element pattern × array factor).
 *
 * E_total(θ,φ) = AF(θ,φ) × E_element(θ,φ)  [Pattern Multiplication Theorem]
 *
 * @param config Array configuration.
 * @param elements Array elements.
 * @param weights Complex excitation weights.
 * @param theta Polar angle (rad).
 * @param phi Azimuthal angle (rad).
 * @return pa_af_result_t Total pattern magnitude/phase.
 */
pa_af_result_t pa_total_pattern(const pa_array_config_t *config,
                                 const pa_element_t *elements,
                                 const double complex *weights,
                                 double theta, double phi);

/**
 * @brief Compute the array factor magnitude in dB (normalized to peak).
 *
 * AF_dB(θ,φ) = 20 log₁₀(|AF(θ,φ)| / max|AF|)
 *
 * @param result Pointer to result, magnitude_db field will be set.
 */
void pa_af_normalize_db(pa_af_result_t *result, double max_af_magnitude);

/**
 * @brief Compute the wavelength λ = c / f.
 */
double pa_wavelength(double frequency_hz);

/**
 * @brief Compute the free-space wavenumber k₀ = 2π / λ.
 */
double pa_wavenumber(double frequency_hz);

/* ============================================================================
 * L5: Steering Vector Computation
 * ============================================================================ */

/**
 * @brief Compute the steering vector for a given scan direction.
 *
 * The steering vector a(θ,φ) has elements:
 *   a_n = exp(-j k₀ (r_n · û_scan))
 * where û_scan is the unit vector in the desired beam direction.
 *
 * Physically: applying weights w_n = a_n^* produces constructive interference
 * in direction (θ_scan, φ_scan), i.e., a phased-array beam.
 *
 * Length: config->num_elements. Caller must free the returned array.
 *
 * @param config Array geometry config.
 * @param elements Array of elements.
 * @param theta_steer Target polar angle (rad).
 * @param phi_steer Target azimuth (rad).
 * @return double complex* Steering vector (caller frees).
 *
 * Reference: Van Trees (2002) Optimum Array Processing Ch. 2.
 */
double complex *pa_steering_vector(const pa_array_config_t *config,
                                    const pa_element_t *elements,
                                    double theta_steer, double phi_steer);

/**
 * @brief Compute per-element phase shifts (radians) for beam steering.
 *
 * Phase shift for element n:
 *   φ_n = -k₀ (r_n · û_scan)  [mod 2π]
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Target θ (rad).
 * @param phi_steer Target φ (rad).
 * @param phase_shifts Output array (caller allocates, size = num_elements).
 */
void pa_phase_shifts(const pa_array_config_t *config,
                     const pa_element_t *elements,
                     double theta_steer, double phi_steer,
                     double *phase_shifts);

/**
 * @brief Quantize phase shifts to finite-bit resolution.
 *
 * Each phase shift is rounded to the nearest 2π/2^N step.
 * Phase quantization causes beam pointing error and raised sidelobes.
 *
 * @param phase_shifts Input/output array of phase shifts (rad).
 * @param num_elements Number of elements.
 * @param bits Number of phase shifter bits.
 *
 * Reference: Mailloux (2005) §3.4 — Quantization effects.
 */
void pa_quantize_phase_shifts(double *phase_shifts, uint32_t num_elements,
                              pa_phase_bits_t bits);

/* ============================================================================
 * L5: Array Synthesis — Amplitude Taper Windows
 * ============================================================================ */

/**
 * @brief Generate amplitude taper coefficients for sidelobe control.
 *
 * Returns pointer to taper array (length = num_elements). Caller frees.
 *
 * @param num_elements Number of elements.
 * @param window_type Window/taper type.
 * @param window_param Parameter:
 *   Chebyshev: desired SLL (dB, negative, e.g., -30)
 *   Taylor: desired SLL (dB, negative) + nbar (embedded in integer part:
 *           abs(param) = SLL, fraction = nbar/10, e.g. -30.5 → SLL=-30dB, nbar=5)
 *   Kaiser: β parameter
 *   Gaussian: σ parameter
 *   Others: ignored.
 * @return double* Taper coefficients (amplitude weights, callers frees).
 */
double *pa_amplitude_taper(uint32_t num_elements, pa_window_t window_type,
                           double window_param);

/**
 * @brief Compute Dolph-Chebyshev amplitude coefficients for a linear array.
 *
 * Yields constant sidelobe level at specified dB below peak.
 * Based on Chebyshev polynomials T_n(x).
 * For symmetric linear arrays only; returns N-element taper.
 *
 * Derivation:
 *   Given desired SLL = -R_dB, define R = 10^{-R_dB/20}.
 *   x₀ = cosh(acosh(R)/(N-1)).
 *   Weights determined from roots of shifted Chebyshev polynomial.
 *
 * @param num_elements Number of elements.
 * @param sll_db Desired sidelobe level (dB, negative, e.g., -30.0).
 * @return double* Dolph-Chebyshev amplitude weights.
 *
 * Reference: Balanis (2016) §6.8.3 — Dolph-Chebyshev Array.
 */
double *pa_dolph_chebyshev_weights(uint32_t num_elements, double sll_db);

/**
 * @brief Compute Taylor (one-parameter) amplitude distribution for line source.
 *
 * Produces a tapered distribution with specified near-in sidelobe level
 * that then decays monotonically. More practical than Chebyshev for
 * large arrays because far-out sidelobes decay.
 *
 * @param num_elements Number of array points along the aperture.
 * @param sll_db Desired peak sidelobe level (dB, negative).
 * @param nbar Number of equal-level sidelobes before decay (≥ 2).
 * @return double* Taylor amplitude weights.
 *
 * Reference: Balanis (2016) §6.8.4 — Taylor Line-Source (Tschebyscheff Error).
 *            Taylor (1955) "Design of Line-Source Antennas for Narrow Beamwidth
 *            and Low Side Lobes", IRE Trans. AP-3.
 */
double *pa_taylor_weights(uint32_t num_elements, double sll_db, uint32_t nbar);

/**
 * @brief Compute binomial array coefficients (zero sidelobes, maximum BW).
 *
 * Coefficients follow the binomial distribution: C(N-1, n) normalized.
 * Produces the widest beam but has exactly zero sidelobes.
 *
 * @param num_elements Number of elements.
 * @return double* Binomial weights.
 */
double *pa_binomial_weights(uint32_t num_elements);

/* ============================================================================
 * L3: Coordinate Transforms
 * ============================================================================ */

/**
 * @brief Convert azimuth/elevation (degrees) to spherical (radians).
 *
 * Standard radar convention:
 *   azimuth: angle from x-axis in the xy-plane
 *   elevation: angle from xy-plane toward z-axis
 *
 * Transformation:
 *   θ = π/2 - elevation     (polar angle from z-axis)
 *   φ = azimuth              (azimuthal angle from x-axis)
 */
pa_spherical_t pa_azel_to_spherical(double azimuth_deg, double elevation_deg);

/**
 * @brief Convert spherical (radians) to azimuth/elevation (degrees).
 */
pa_azel_t pa_spherical_to_azel(double theta_rad, double phi_rad);

/**
 * @brief Convert spherical direction to unit vector (x,y,z).
 */
void pa_spherical_to_uvw(double theta, double phi, double *u, double *v, double *w);

/**
 * @brief Convert unit vector to spherical coordinates.
 */
pa_spherical_t pa_uvw_to_spherical(double u, double v, double w);

/* ============================================================================
 * L5: Array Initialization and Memory Management
 * ============================================================================ */

/**
 * @brief Initialize a uniform linear array (ULA) along x-axis.
 *
 * Creates N equally spaced elements along the x-axis at y = z = 0.
 * Spacing d = config->element_spacing_x.
 * All elements initialized as isotropic with unit excitation.
 *
 * @param config Array configuration to fill.
 * @param num_elements Number of elements along x.
 * @param element_spacing Spacing between elements (meters).
 * @param frequency_hz Operating frequency (Hz).
 */
void pa_init_linear_array(pa_array_config_t *config, uint32_t num_elements,
                          double element_spacing, double frequency_hz);

/**
 * @brief Initialize a uniform rectangular planar array (URPA) in xy-plane.
 *
 * Creates rows × cols elements in a rectangular grid.
 * Spacing dx, dy. Elements at z = 0.
 *
 * @param config Array config to fill.
 * @param rows Number of rows (y-direction).
 * @param cols Number of columns (x-direction).
 * @param dx Column spacing (x, meters).
 * @param dy Row spacing (y, meters).
 * @param frequency_hz Operating frequency.
 */
void pa_init_planar_array(pa_array_config_t *config,
                          uint32_t rows, uint32_t cols,
                          double dx, double dy, double frequency_hz);

/**
 * @brief Initialize a uniform circular array (UCA) in xy-plane.
 *
 * N elements equally spaced around a circle of given radius.
 *
 * @param config Array config to fill.
 * @param num_elements Number of elements around the ring.
 * @param radius Radius of the circular array (meters).
 * @param frequency_hz Operating frequency.
 */
void pa_init_circular_array(pa_array_config_t *config, uint32_t num_elements,
                            double radius, double frequency_hz);

/**
 * @brief Allocate and initialize element array from geometry config.
 *
 * Elements are placed according to the geometry type.
 * All elements default to isotropic with unit excitation.
 *
 * @param config Array geometry config (must be initialized).
 * @return pa_element_t* Array of elements (caller frees).
 */
pa_element_t *pa_allocate_elements(const pa_array_config_t *config);

/**
 * @brief Free an element array.
 */
void pa_free_elements(pa_element_t *elements);

/**
 * @brief Allocate and populate a pa_pattern_t for full scan.
 *
 * @param num_theta Number of θ samples.
 * @param num_phi Number of φ samples.
 * @param theta_start Start angle (rad).
 * @param theta_end End angle (rad).
 * @param phi_start Start angle (rad).
 * @param phi_end End angle (rad).
 * @return pa_pattern_t* Fully allocated pattern (caller frees via pa_free_pattern).
 */
pa_pattern_t *pa_allocate_pattern(uint32_t num_theta, uint32_t num_phi,
                                   double theta_start, double theta_end,
                                   double phi_start, double phi_end);

/**
 * @brief Free a full pattern structure.
 */
void pa_free_pattern(pa_pattern_t *pattern);

/* ============================================================================
 * L6: Beam Pattern Computation and Analysis
 * ============================================================================ */

/**
 * @brief Compute the full 2D beam pattern over a specified angular range.
 *
 * For each (θ,φ) in the scan grid, computes the array factor and fills
 * pattern->af_db. Also computes key metrics: peak directivity,
 * sidelobe level, 3dB beamwidth, first-null beamwidth.
 *
 * @param config Array configuration.
 * @param elements Array elements.
 * @param weights Complex excitation weights.
 * @param pattern Pre-allocated pattern structure (from pa_allocate_pattern).
 *
 * O(N · num_theta · num_phi) complexity.
 */
void pa_compute_pattern(const pa_array_config_t *config,
                        const pa_element_t *elements,
                        const double complex *weights,
                        pa_pattern_t *pattern);

/**
 * @brief Compute and return the peak directivity from the pattern.
 *
 * Directivity D_max = 4π / Ω_A (Ω_A = total radiated solid angle).
 * For uniform arrays: D_max ≈ 4π A_eff / λ².
 *
 * @param config Array configuration.
 * @param elements Array elements.
 * @param weights Complex weights.
 * @return double Maximum directivity (linear, not dB).
 */
double pa_directivity(const pa_array_config_t *config,
                      const pa_element_t *elements,
                      const double complex *weights);

/**
 * @brief Compute 3dB beamwidth for a specified scan angle.
 *
 * Finds the angular width between -3dB points of the main lobe
 * when the beam is steered to (θ_steer,φ_steer).
 *
 * Approximation for uniform linear broadside array (no scan):
 *   HPBW ≈ 0.886 λ/(N d)  [radians]  ≈ 50.8° λ/(N d)
 *
 * For steered beams:
 *   HPBW(θ_s) ≈ HPBW_broadside / cos(θ_s)     (1D scan)
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param weights Complex weights.
 * @param theta_steer Steered direction θ (rad).
 * @param phi_steer Steered direction φ (rad).
 * @param hp_beamwidth_deg Output: 3dB beamwidth (degrees).
 *
 * Reference: Balanis §6.5.2.
 */
void pa_beamwidth_3db(const pa_array_config_t *config,
                      const pa_element_t *elements,
                      const double complex *weights,
                      double theta_steer, double phi_steer,
                      double *hp_beamwidth_deg);

/**
 * @brief Find all grating lobes visible in visible space.
 *
 * Grating lobes occur at angles where:
 *   sin(θ_g) = sin(θ_s) ± m λ/d   (linear array)
 * Condition for visible grating lobe: |sin(θ_g)| ≤ 1.
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param weights Complex weights.
 * @param theta_steer Steering angle θ (rad).
 * @param phi_steer Steering angle φ (rad).
 * @param max_lobes Max number of lobes to return.
 * @param lobes Output array (caller allocates size ≥ max_lobes).
 * @return int Number of grating lobes found.
 *
 * Reference: Balanis §6.9, Mailloux §2.1.
 */
int pa_find_grating_lobes(const pa_array_config_t *config,
                          const pa_element_t *elements,
                          const double complex *weights,
                          double theta_steer, double phi_steer,
                          uint32_t max_lobes, pa_grating_lobe_t *lobes);

/* ============================================================================
 * L7: AESA Radar Equation — Radar Performance Metrics
 * ============================================================================ */

/**
 * @brief Compute radar range equation for a phased array.
 *
 * Radar Range Equation (monostatic, single pulse):
 *   R⁴ = (P_t · G_t · G_r · λ² · σ) / ((4π)³ · k · T₀ · B · F_n · SNR_min · L)
 *
 * For phased array with N elements:
 *   G_t = D_t ≈ 4π A_eff / λ²  (effective aperture from N elements)
 *   A_eff = N · A_element · η   (η = aperture efficiency)
 *
 * @param config Array configuration.
 * @param elements Array elements.
 * @param target_rcs_m2 Target radar cross section (m²).
 * @param snr_min_db Minimum detectable SNR (dB).
 * @param system_loss_db Total system losses (dB, positive).
 * @return double Maximum detection range (meters).
 *
 * Reference: Richards et al. (2010) §2.2–2.5,
 *            Skolnik (2008) Radar Handbook §1.5.
 */
double pa_radar_range_equation(const pa_array_config_t *config,
                                const pa_element_t *elements,
                                double target_rcs_m2,
                                double snr_min_db, double system_loss_db);

/**
 * @brief Compute the effective aperture of the phased array.
 *
 * A_eff = N · A_element · η
 * A_element = effective area per element
 * η = aperture efficiency (0.5–0.8 typical for patches)
 *
 * For half-wave dipole: A_element ≈ λ² · G_element / (4π)
 * For patch: A_element ≈ dx · dy · η_element
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param aperture_efficiency η (0.0 to 1.0).
 * @return double Effective aperture area (m²).
 */
double pa_effective_aperture(const pa_array_config_t *config,
                             const pa_element_t *elements,
                             double aperture_efficiency);

/* ============================================================================
 * L8: Adaptive Beamforming — MVDR and LCMV
 * ============================================================================ */

/**
 * @brief Initialize an SMI-based beamformer.
 *
 * @param num_elements Number of array elements.
 * @param num_snapshots Number of temporal snapshots for covariance estimation.
 * @return pa_smi_beamformer_t* Initialized beamformer (caller frees).
 */
pa_smi_beamformer_t *pa_smi_beamformer_init(uint32_t num_elements,
                                             uint32_t num_snapshots);

/**
 * @brief Free an SMI beamformer.
 */
void pa_smi_beamformer_free(pa_smi_beamformer_t *bf);

/**
 * @brief Estimate the sample covariance matrix from snapshot data.
 *
 * R̂ = (1/K) Σ_{k=0}^{K-1} x[k] x[k]^H
 * where x[k] is the M-element snapshot vector at time k.
 *
 * Diagonal loading is applied for regularization:
 * R̂_loaded = R̂ + γ I  (γ = diagonal_loading * trace(R̂)/M)
 *
 * @param bf Beamformer state.
 * @param snapshots 2D array: snapshots[snap_idx][elem_idx], size K×M.
 */
void pa_smi_estimate_covariance(pa_smi_beamformer_t *bf,
                                 const double complex **snapshots);

/**
 * @brief Compute the MVDR weight vector.
 *
 * MVDR (Minimum Variance Distortionless Response, aka Capon beamformer):
 *   minimize w^H R w   subject to   w^H a(θ_s) = 1
 *
 * Solution:
 *   w_MVDR = R^{-1} a(θ_s) / (a(θ_s)^H R^{-1} a(θ_s))
 *
 * This preserves the signal from direction (θ_s, φ_s) while minimizing
 * the total output power from all other directions (interference + noise).
 *
 * @param bf MVDR beamformer state.
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Look direction θ (rad).
 * @param phi_steer Look direction φ (rad).
 *
 * Reference: Capon (1969) "High-Resolution Frequency-Wavenumber Spectrum
 *            Analysis", Proc. IEEE. Van Trees (2002) Optimum Array Processing.
 */
void pa_mvdr_compute_weights(pa_mvdr_beamformer_t *bf,
                              const pa_array_config_t *config,
                              const pa_element_t *elements,
                              double theta_steer, double phi_steer);

/**
 * @brief Compute the Capon spatial power spectrum.
 *
 * P_Capon(θ,φ) = 1 / (a(θ,φ)^H R^{-1} a(θ,φ))
 *
 * Evaluated over a grid of angles. Peaks indicate source directions.
 *
 * @param bf MVDR beamformer.
 * @param config Array config.
 * @param elements Array elements.
 * @param spectrum Output array (size = bf->num_angles).
 * @param theta_vals Theta angles to evaluate.
 * @param phi_vals Phi angles to evaluate.
 * @param num_angles Number of angular bins.
 */
void pa_mvdr_capon_spectrum(pa_mvdr_beamformer_t *bf,
                             const pa_array_config_t *config,
                             const pa_element_t *elements,
                             double *spectrum,
                             const double *theta_vals,
                             const double *phi_vals,
                             uint32_t num_angles);

/**
 * @brief Compute the output SINR of the adaptive beamformer.
 *
 * SINR = (w^H R_s w) / (w^H R_{i+n} w)
 * where R_s = signal covariance, R_{i+n} = interference+noise covariance.
 *
 * For MVDR at optimality:
 *   SINR_opt = a^H R^{-1} a  (if signal perfectly known)
 *
 * @param bf MVDR beamformer.
 * @return double Output SINR (dB).
 */
double pa_mvdr_output_sinr(const pa_mvdr_beamformer_t *bf);

/**
 * @brief Initialize an MVDR beamformer.
 */
pa_mvdr_beamformer_t *pa_mvdr_beamformer_init(uint32_t num_elements,
                                               uint32_t num_snapshots);

/**
 * @brief Free an MVDR beamformer.
 */
void pa_mvdr_beamformer_free(pa_mvdr_beamformer_t *bf);

/* ============================================================================
 * L6: Monopulse Angle Estimation
 * ============================================================================ */

/**
 * @brief Generate sum (Σ) and difference (Δ) patterns for monopulse tracking.
 *
 * For a linear array split into two halves:
 *   Σ(θ) = AF_left(θ) + AF_right(θ)   (used for detection)
 *   Δ(θ) = AF_left(θ) - AF_right(θ)   (used for angle error)
 *
 * The error signal ε = Re{Δ/Σ} is proportional to the angle offset
 * near boresight (zero at boresight, linear in small neighborhood).
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param weights Complex weights (for Σ pattern — Δ uses modified weights).
 * @param theta_vals Array of theta values to evaluate.
 * @param num_theta Number of theta values.
 * @param sum_pattern Output Σ pattern magnitudes.
 * @param diff_pattern Output Δ pattern magnitudes.
 * @param error_signal Output error signal (Re{Δ/Σ}).
 *
 * Reference: Skolnik (2008) §8.7, Sherman & Barton (2011)
 *            "Monopulse Principles and Techniques".
 */
void pa_monopulse_patterns(const pa_array_config_t *config,
                           const pa_element_t *elements,
                           const double complex *weights,
                           const double *theta_vals, uint32_t num_theta,
                           double *sum_pattern, double *diff_pattern,
                           double *error_signal);

/**
 * @brief Estimate DOA from monopulse error signal using slope calibration.
 *
 * θ_estimated = θ_boresight + ε / k_m
 * where k_m = dε/dθ|_{θ=0} is the monopulse slope factor.
 *
 * @param error_signal Measured ε = Re{Δ/Σ}.
 * @param slope_factor Monopulse slope k_m (rad⁻¹).
 * @param boresight_rad Boresight angle (rad).
 * @return double Estimated angle (rad).
 */
double pa_monopulse_estimate_angle(double error_signal,
                                    double slope_factor,
                                    double boresight_rad);

/* ============================================================================
 * L6: True Time Delay and Beam Squint
 * ============================================================================ */

/**
 * @brief Compute per-element true-time-delay values for wideband beamforming.
 *
 * For element at position r_n, delay to steer to direction û_s:
 *   τ_n = (r_n · û_s) / c
 *
 * TTD eliminates beam squint because the delay is independent of frequency.
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param theta_steer Steered θ (rad).
 * @param phi_steer Steered φ (rad).
 * @param delays_ps Output delays in picoseconds (caller allocates).
 */
void pa_ttd_compute_delays(const pa_array_config_t *config,
                            const pa_element_t *elements,
                            double theta_steer, double phi_steer,
                            double *delays_ps);

/**
 * @brief Analyze beam squint for phase-steered arrays.
 *
 * For a phase-steered array at center frequency f₀, at frequency f the beam
 * squints because the phase shift φ = -2π(f₀/c)(r_n·û_s) no longer produces
 * the correct spatial phase at frequency f.
 *
 * Actual steer angle at frequency f:
 *   sin(θ(f)) = (f₀/f) sin(θ_s)
 *
 * Squint = θ(f) - θ(f₀).
 *
 * @param config Array config.
 * @param theta_steer_deg Nominal steer angle (deg).
 * @param squint Output beam squint analysis.
 *
 * Reference: Mailloux (2005) §5.2 — Bandwidth and Squint.
 */
void pa_beam_squint_analysis(const pa_array_config_t *config,
                              double theta_steer_deg,
                              pa_beam_squint_t *squint);

/* ============================================================================
 * L5: Mutual Coupling
 * ============================================================================ */

/**
 * @brief Apply mutual coupling correction using an S-parameter matrix.
 *
 * For M elements, the mutual coupling matrix C is M×M.
 * The isolated-element excitation vector w is modified to:
 *   w_corrected = C^{-1} w
 *
 * This compensates for element-to-element coupling effects.
 *
 * @param num_elements Number of elements.
 * @param s_matrix M×M complex S-parameter matrix (row-major).
 * @param weights Input/output: excitation weights (corrected in place).
 *
 * Reference: Gupta & Ksienski (1983) "Effect of Mutual Coupling on the
 *            Performance of Adaptive Arrays", IEEE Trans. AP-31.
 */
void pa_mutual_coupling_correct(uint32_t num_elements,
                                 const double complex *s_matrix,
                                 double complex *weights);

/**
 * @brief Generate an S-parameter coupling matrix for uniform linear array.
 *
 * Simplified exponential decay model for identical elements:
 *   S_{m,n} = s0 · exp(-|m-n| · α) · exp(-j k₀ d |m-n|)
 * where s0 is the nearest-neighbor coupling coefficient,
 * α is the decay rate, d is element spacing.
 *
 * Self-coupling (diagonal): S_{n,n} = 0 (reflectionless, matched).
 *
 * @param num_elements Number of elements.
 * @param spacing_m Element spacing (m).
 * @param frequency_hz Frequency (Hz).
 * @param coupling_coefficient Nearest-neighbor S-parameter |S_{n,n+1}|.
 * @param decay_rate Exponential decay rate α.
 * @param s_matrix Output: M×M complex S-matrix (row-major, caller allocates).
 */
void pa_compute_coupling_matrix(uint32_t num_elements,
                                 double spacing_m, double frequency_hz,
                                 double coupling_coefficient,
                                 double decay_rate,
                                 double complex *s_matrix);

/* ============================================================================
 * L5: AESA T/R Module Functions
 * ============================================================================ */

/**
 * @brief Initialize a T/R module with realistic AESA parameters.
 *
 * Models typical X-band T/R module (F-35 AN/APG-81 style):
 *   GaAs/GaN MMIC, ~10W peak TX power, ~4-bit phase shifter, etc.
 *
 * @param module_id T/R module identifier.
 * @return pa_tr_module_t Initialized module.
 */
pa_tr_module_t pa_tr_module_init(int module_id);

/**
 * @brief Compute noise-equivalent delta temperature (NEDT) for the array.
 *
 * @param tr_array Array of T/R modules.
 * @param num_modules Number of modules.
 * @return double NEDT (Kelvin).
 */
double pa_array_nedt(const pa_tr_module_t *tr_array, uint32_t num_modules);

/**
 * @brief Compute array-level G/T (gain over noise temperature) figure of merit.
 *
 * G/T = G_rx (dB) - 10 log₁₀(T_sys)
 * where G_rx is the receive gain and T_sys is the system noise temperature.
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param system_noise_temp_k System noise temperature (K).
 * @return double G/T (dB/K).
 */
double pa_array_gt_fom(const pa_array_config_t *config,
                        const pa_element_t *elements,
                        double system_noise_temp_k);

#ifdef __cplusplus
}
#endif

#endif /* PHASED_ARRAY_H */
