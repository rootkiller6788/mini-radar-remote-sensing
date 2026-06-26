/**
 * @file    hyperspectral_core.h
 * @brief   Core Hyperspectral Remote Sensing — L1 Definitions + L3 Mathematical Structures
 *
 * @details Fundamental data structures: spectral datacube, band definitions,
 *          noise models, radiometric quantities, spectral sampling relations.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Hyperspectral datacube (x, y, lambda) — 3D spatial-spectral data
 *     - Spectral band / channel — narrow contiguous wavelength interval
 *     - Spectral resolution (FWHM), Spatial resolution (IFOV, GSD)
 *     - Radiometric resolution (bit depth, NEDL)
 *     - Signal-to-Noise Ratio (SNR), NEdT, NESR
 *     - Atmospheric windows (VIS, SWIR, MWIR, LWIR)
 *   L2 - Core Concepts:
 *     - Imaging spectroscopy — spatial and spectral data simultaneously
 *     - Pushbroom vs whiskbroom vs staring sensor geometries
 *   L3 - Mathematical Structures:
 *     - Tensor/3D array representation of datacube
 *     - Linear algebra for spectra (R^B vector spaces)
 *     - Convex geometry of spectral mixtures
 *
 * Reference:
 *   - Manolakis, Lockwood, Cooley, "Hyperspectral Imaging Remote Sensing" (2016)
 *   - Borengasser, Hungate, Watkins, "Hyperspectral Remote Sensing" (2008)
 *   - Shaw & Burke, "Spectral Imaging for Remote Sensing" (2003)
 */

#ifndef HYPERSPECTRAL_CORE_H
#define HYPERSPECTRAL_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ─── Fundamental Physical Constants (L1) ────────────────────────────── */

/** Planck constant [J*s] */
#define HSPEC_PLANCK_h              6.62607015e-34
/** Speed of light in vacuum [m/s] */
#define HSPEC_SPEED_OF_LIGHT_c      2.99792458e8
/** Boltzmann constant [J/K] */
#define HSPEC_BOLTZMANN_k           1.380649e-23
/** Stefan-Boltzmann constant [W*m^-2*K^-4] */
#define HSPEC_STEFAN_BOLTZMANN_sigma 5.670374419e-8
/** Wien displacement constant [m*K] */
#define HSPEC_WIEN_b                2.897771955e-3

/* ─── Radiometric bit depth (L1) ─────────────────────────────────────── */

typedef enum {
    HSPEC_BITS_8  = 8,
    HSPEC_BITS_10 = 10,
    HSPEC_BITS_12 = 12,
    HSPEC_BITS_14 = 14,
    HSPEC_BITS_16 = 16
} hspec_bit_depth_t;

typedef enum {
    HSPEC_UNIT_NANOMETER  = 0,
    HSPEC_UNIT_MICROMETER = 1,
    HSPEC_UNIT_WAVENUMBER = 2
} hspec_wavelength_unit_t;

/* ─── L1: Spectral Band ──────────────────────────────────────────────── */
typedef struct {
    size_t   band_index;
    double   center_wavelength;
    double   fwhm;
    double   lower_bound;
    double   upper_bound;
    double   spectral_sampling;
    char     band_name[32];
} hspec_band_t;

/* ─── L1: Spectral Response Function ─────────────────────────────────── */
typedef struct {
    double   center_wavelength;
    double   fwhm;
    double   sigma;
    double   efficiency;
} hspec_srf_t;

/* ─── L1: Hyperspectral Datacube (BIP layout) ────────────────────────── */
typedef struct {
    double  *data;
    size_t   nrows;
    size_t   ncols;
    size_t   nbands;
    size_t   npixels;
    size_t   total_elements;
    double  *wavelengths;
    double  *fwhm;
    char     sensor_name[64];
    double   spatial_resolution;
    double   temporal_acq_time;
} hspec_datacube_t;

/* ─── L1: Single-Pixel Spectrum ──────────────────────────────────────── */
typedef struct {
    size_t   nbands;
    double  *reflectance;
    size_t   row;
    size_t   col;
    double   latitude;
    double   longitude;
    char     label[64];
} hspec_pixel_t;

/* ─── L1: Spectral Library Entry ─────────────────────────────────────── */
typedef struct {
    size_t   nbands;
    double  *reflectance;
    double  *wavelengths;
    char     material_name[128];
    char     category[64];
    char     source[64];
} hspec_spectral_lib_entry_t;

/* ─── L1: Sensor Noise Model ─────────────────────────────────────────── */
typedef struct {
    double   photon_noise_coeff;
    double   dark_current_e;
    double   readout_noise_e;
    double   quantization_noise;
    double   total_noise_model;
    double   saturation_level;
    double   operating_temp_k;
} hspec_noise_model_t;

/* ─── L1: Radiometric Quantities ─────────────────────────────────────── */
typedef struct {
    size_t   nbands;
    double  *radiance;
    double  *gain;
    double  *offset;
    double  *noise_equivalent_delta_radiance;
    double   solar_zenith_angle;
    double   solar_azimuth;
    double   earth_sun_distance;
} hspec_radiance_t;

/* ─── L1: Reflectance ────────────────────────────────────────────────── */
typedef struct {
    size_t   nbands;
    double  *toa_reflectance;
    double  *surface_reflectance;
    double  *exo_atmospheric_irradiance;
} hspec_reflectance_t;

/* ─── L1: Detection Figures of Merit ─────────────────────────────────── */
typedef struct {
    size_t   nbands;
    double  *snr_per_band;
    double  *nedt_per_band;
    double  *nesr_per_band;
    double  *contrast_per_band;
    double   integrated_snr;
    double   peak_snr;
    double   peak_snr_band;
    double  *wavelengths;          /**< Optional: band center wavelengths for peak SNR band ID */
} hspec_detection_metrics_t;

/* ─── L1: Spectral Indices ───────────────────────────────────────────── */
typedef struct {
    double   ndvi;
    double   ndwi;
    double   evi;
    double   savi;
    double   ndbi;
    double   pri;
    double   red_edge_position;
    double   red_edge_slope;
    double   water_absorption_depth;
} hspec_spectral_indices_t;

/* ─── L3: Band Statistics ────────────────────────────────────────────── */
typedef struct {
    double   mean;
    double   variance;
    double   skewness;
    double   kurtosis;
    double   min_val;
    double   max_val;
    double   median;
} hspec_band_statistics_t;

/* ─── L3: Band Correlation ───────────────────────────────────────────── */
typedef struct {
    double   covariance;
    double   correlation;
    double   mutual_information;
} hspec_band_correlation_t;

/* ─── L3: Spectral Similarity Metrics ────────────────────────────────── */
typedef struct {
    double   euclidean_distance;
    double   spectral_angle_rad;
    double   spectral_angle_deg;
    double   sid;
    double   cross_correlation;
    double   chebyshev_distance;
} hspec_spectral_similarity_t;

/* ─── L1: Atmospheric Window ─────────────────────────────────────────── */
typedef struct {
    double   window_start;
    double   window_end;
    double   peak_transmission;
    char     window_name[32];
} hspec_atmospheric_window_t;

/* Standard atmospheric windows (MODTRAN/LOWTRAN) */
#define HSPEC_WINDOW_VIS_START   0.38
#define HSPEC_WINDOW_VIS_END     0.75
#define HSPEC_WINDOW_SWIR1_START 1.00
#define HSPEC_WINDOW_SWIR1_END   1.30
#define HSPEC_WINDOW_SWIR2_START 1.50
#define HSPEC_WINDOW_SWIR2_END   1.75
#define HSPEC_WINDOW_SWIR3_START 2.00
#define HSPEC_WINDOW_SWIR3_END   2.50
#define HSPEC_WINDOW_MWIR_START  3.00
#define HSPEC_WINDOW_MWIR_END    5.00
#define HSPEC_WINDOW_LWIR_START  8.00
#define HSPEC_WINDOW_LWIR_END   14.00

/* Spectral region enumeration */
typedef enum {
    HSPEC_REGION_UV    = 0,
    HSPEC_REGION_VIS   = 1,
    HSPEC_REGION_VNIR  = 2,
    HSPEC_REGION_NIR   = 3,
    HSPEC_REGION_SWIR  = 4,
    HSPEC_REGION_MWIR  = 5,
    HSPEC_REGION_LWIR  = 6,
    HSPEC_REGION_COUNT = 7
} hspec_spectral_region_t;

typedef enum {
    HSPEC_GEOM_PUSHBROOM  = 0,
    HSPEC_GEOM_WHISKBROOM = 1,
    HSPEC_GEOM_STARING    = 2,
    HSPEC_GEOM_FOURIER    = 3
} hspec_sensor_geometry_t;

/* ─── API Declarations ───────────────────────────────────────────────── */

/* Datacube: O(nrows*ncols*nbands) */
hspec_datacube_t hspec_datacube_alloc(size_t nrows, size_t ncols, size_t nbands);
void hspec_datacube_free(hspec_datacube_t *dc);

/* Pixel/band extraction */
int hspec_datacube_extract_pixel(const hspec_datacube_t *dc, size_t row,
                                  size_t col, hspec_pixel_t *pixel);
int hspec_datacube_extract_band(const hspec_datacube_t *dc, size_t band_idx,
                                 double *image_out);

/* Statistics */
int hspec_datacube_band_statistics(const hspec_datacube_t *dc,
                                    hspec_band_statistics_t *stats);
int hspec_datacube_covariance(const hspec_datacube_t *dc, double *cov_matrix);
int hspec_datacube_correlation(const hspec_datacube_t *dc, double *corr_matrix);
int hspec_datacube_standardize(hspec_datacube_t *dc);

/* Spectral library */
int hspec_lib_entry_init(hspec_spectral_lib_entry_t *entry, size_t nbands,
                          const char *name);
void hspec_lib_entry_free(hspec_spectral_lib_entry_t *entry);

/* Spectral similarity (SAM, SID, Euclidean, etc.) */
int hspec_spectral_similarity(const double *a, const double *b, size_t nbands,
                               hspec_spectral_similarity_t *sim);

/* Spectral indices (NDVI, NDWI, EVI, SAVI, etc.) */
int hspec_compute_spectral_indices(const hspec_pixel_t *pixel,
                                    const double *wavelengths,
                                    hspec_spectral_indices_t *indices);

/* Band configuration */
int hspec_create_uniform_bands(double start_wavelength, double end_wavelength,
                                size_t nbands, double *wavelengths_out,
                                double *fwhm_out);

/* SRF evaluation: SRF(lambda) = exp(-(lambda-lambda_c)^2/(2*sigma^2)) */
double hspec_evaluate_srf(const hspec_srf_t *srf, double lambda);

/* Spectral resampling via SRF convolution */
int hspec_spectral_resample(const double *highres_lambda,
                             const double *highres_radiance,
                             size_t M, const hspec_srf_t *srf_array,
                             size_t N, double *band_radiance_out);

/* Detection metrics (SNR, NEdT, NESR per band) */
int hspec_compute_detection_metrics(const double *signal, const double *noise_std,
                                     size_t nbands,
                                     hspec_detection_metrics_t *metrics_out);

/* Virtual dimensionality via HFC method (Harsanyi-Farrand-Chang) */
size_t hspec_virtual_dimensionality(const double *cov_matrix, size_t nbands,
                                     size_t npixels, double p_fa);

#endif /* HYPERSPECTRAL_CORE_H */