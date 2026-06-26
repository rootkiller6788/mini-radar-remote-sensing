/**
 * @file    hyperspectral_classification.h
 * @brief   Hyperspectral Classification and Detection — L5 Algorithms + L6 Canonical Problems
 *
 * @details Implements spectral classification algorithms: Spectral Angle Mapper
 *          (SAM), Matched Filter, Adaptive Coherence Estimator (ACE), Reed-Xiaoli
 *          (RX) anomaly detection, k-means clustering, and maximum likelihood.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Confusion matrix, Overall Accuracy (OA), Kappa coefficient
 *     - Endmember, class prototype, training/testing split
 *   L2 - Core Concepts:
 *     - Supervised vs unsupervised classification
 *     - Hard vs soft classification
 *     - Spectral variability, intra-class vs inter-class
 *   L3 - Mathematical Structures:
 *     - Mahalanobis distance, Bayes decision rule
 *     - Quadratic/linear discriminant functions
 *   L5 - Algorithms:
 *     - Spectral Angle Mapper (SAM)
 *     - Matched Filter (MF), Adaptive Coherence Estimator (ACE)
 *     - Reed-Xiaoli (RX) anomaly detector
 *     - K-means / ISODATA clustering
 *     - Maximum Likelihood (ML) with Gaussian class models
 *   L6 - Canonical Problems:
 *     - Mineral mapping from AVIRIS data
 *     - Target detection (subpixel)
 *     - Anomaly detection in cluttered backgrounds
 *
 * Reference:
 *   - Richards & Jia, "Remote Sensing Digital Image Analysis" (2006)
 *   - Manolakis & Shaw, "Detection Algorithms for Hyperspectral Imaging" (2002)
 *   - Reed & Yu, "Adaptive Multiple-Band CFAR Detection" (1990), IEEE TASSP
 */

#ifndef HYPERSPECTRAL_CLASSIFICATION_H
#define HYPERSPECTRAL_CLASSIFICATION_H

#include "hyperspectral_core.h"

/* ─── L1: Classification Label Map ──────────────────────────────────── */

typedef struct {
    size_t   nrows;
    size_t   ncols;
    size_t  *labels;             /**< Class label per pixel [nrows*ncols] */
    size_t   nclasses;           /**< Number of distinct classes */
    char   (*class_names)[64];   /**< Class name strings [nclasses] */
} hspec_label_map_t;

/* ─── L1: Confusion Matrix ──────────────────────────────────────────── */

/**
 * @brief Confusion matrix for classification accuracy assessment.
 *
 * M_{ij} = number of pixels of true class i classified as class j.
 * OA = Σ_i M_{ii} / Σ_{i,j} M_{ij}   (overall accuracy)
 * κ  = (OA - p_e) / (1 - p_e)        (Cohen's kappa)
 *    where p_e = Σ_i (row_i / N)·(col_i / N)
 */
typedef struct {
    size_t   nclasses;
    size_t **matrix;             /**< Confusion matrix [nclasses][nclasses] */
    size_t  *row_sums;           /**< Per-class reference counts */
    size_t  *col_sums;           /**< Per-class predicted counts */
    size_t   total_samples;
    double   overall_accuracy;   /**< OA ∈ [0, 1] */
    double   kappa;              /**< Cohen's kappa ∈ [-1, 1] */
    double  *producer_accuracy;  /**< Per-class: diagonal/row_sum */
    double  *user_accuracy;      /**< Per-class: diagonal/col_sum */
} hspec_confusion_matrix_t;

/* ─── L5: Spectral Angle Mapper (SAM) ───────────────────────────────── */

/**
 * @brief SAM classification result for a single pixel.
 *
 * SAM measures the angle between the pixel spectrum vector and the
 * reference spectrum vector in B-dimensional space.
 *
 * @math θ = arccos( ⟨x, r⟩ / (‖x‖·‖r‖) )
 */
typedef struct {
    size_t   nclasses;
    double  *angles_rad;         /**< SAM angle to each class [nclasses] */
    size_t   best_class;         /**< Index of closest class (smallest angle) */
    double   best_angle_rad;     /**< Smallest SAM angle */
    double   best_angle_deg;     /**< Smallest SAM angle in degrees */
    double   confidence;         /**< Normalized confidence ∈ [0, 1] */
} hspec_sam_result_t;

/* ─── L5: Matched Filter ────────────────────────────────────────────── */

/**
 * @brief Matched Filter detection result.
 *
 * MF(x) = (sᵀ·Σ⁻¹·x) / (sᵀ·Σ⁻¹·s)
 * where s = target spectrum, Σ = background covariance, x = test pixel.
 *
 * The MF maximizes SNR under Gaussian background assumption.
 */
typedef struct {
    double   mf_score;           /**< Matched filter output statistic */
    double   threshold;          /**< Detection threshold at given P_FA */
    int      detection;          /**< 1 if mf_score >= threshold, 0 otherwise */
    double   snr_estimate;       /**< Estimated SNR of the match */
} hspec_matched_filter_t;

/* ─── L5: Adaptive Coherence Estimator (ACE) ────────────────────────── */

/**
 * @brief ACE detection result.
 *
 * ACE(x) = (sᵀ·Σ⁻¹·x)² / [(sᵀ·Σ⁻¹·s)·(xᵀ·Σ⁻¹·x)]
 * ACE is invariant to scaling of x (angle-based detection).
 */
typedef struct {
    double   ace_score;          /**< ACE statistic ∈ [0, 1] */
    double   threshold;          /**< Detection threshold */
    int      detection;          /**< 1 = target present */
} hspec_ace_result_t;

/* ─── L5: Reed-Xiaoli (RX) Anomaly Detector ─────────────────────────── */

/**
 * @brief RX anomaly detection — finds pixels that differ from background.
 *
 * RX(x) = (x - μ)ᵀ·Σ⁻¹·(x - μ)  (Mahalanobis distance squared)
 *
 * The RX detector is the GLRT for a known background, unknown target.
 * It detects pixels that are spectrally anomalous relative to the
 * local or global background.
 *
 * Reference: Reed & Yu (1990), IEEE TASSP
 */
typedef struct {
    size_t   nrows;
    size_t   ncols;
    double  *rx_scores;          /**< RX anomaly score per pixel [nrows*ncols] */
    double   threshold;          /**< Detection threshold */
    size_t  *anomaly_mask;       /**< 1 = anomalous, 0 = normal [nrows*ncols] */
    size_t   nanomalies;         /**< Number of detected anomalies */
    double   mean_score;         /**< Mean RX score (should follow χ²_B) */
    double   std_score;          /**< Std of RX scores */
} hspec_rx_result_t;

/* ─── L5: K-Means Clustering ────────────────────────────────────────── */

typedef struct {
    size_t   k;                  /**< Number of clusters */
    size_t   nbands;
    double  *centroids;          /**< Cluster centroids [k][nbands] row-major */
    size_t  *assignments;        /**< Cluster assignment per pixel [npixels] */
    size_t  *cluster_sizes;      /**< Number of pixels per cluster [k] */
    size_t   niterations;        /**< Iterations performed */
    double   total_within_ss;    /**< Total within-cluster sum of squares */
    int      converged;          /**< 1 if converged */
} hspec_kmeans_result_t;

/* ─── API: Label Map ───────────────────────────────────────────────── */

hspec_label_map_t hspec_label_map_alloc(size_t nrows, size_t ncols, size_t nclasses);
void hspec_label_map_free(hspec_label_map_t *lm);

/* ─── API: Confusion Matrix ─────────────────────────────────────────── */

hspec_confusion_matrix_t hspec_confusion_matrix_alloc(size_t nclasses);
void hspec_confusion_matrix_free(hspec_confusion_matrix_t *cm);
int hspec_confusion_matrix_compute(hspec_confusion_matrix_t *cm,
                                    const size_t *truth, const size_t *predicted,
                                    size_t npixels);

/* ─── API: Spectral Angle Mapper ────────────────────────────────────── */

/**
 * Compute SAM angles between a test pixel and all reference spectra.
 *
 * @param pixel         Test pixel spectrum
 * @param ref_spectra   Reference spectra [nclasses][nbands] row-major
 * @param nclasses      Number of classes
 * @param nbands        Number of bands (must match pixel->nbands)
 * @return             SAM result
 *
 * @math θ_i = arccos(⟨x, r_i⟩ / (‖x‖·‖r_i‖))
 *
 * Complexity: O(nclasses·nbands)
 * Reference: Kruse et al. (1993), Remote Sensing of Environment
 */
hspec_sam_result_t hspec_sam_classify(const hspec_pixel_t *pixel,
                                       const double *ref_spectra,
                                       size_t nclasses, size_t nbands);

/**
 * Full-scene SAM classification: classify every pixel against reference spectra.
 *
 * @param dc           Datacube
 * @param ref_spectra  Reference spectra [nclasses][nbands] row-major
 * @param nclasses     Number of classes
 * @return             Label map with class assignments
 *
 * Complexity: O(n_pixels·nclasses·nbands)
 */
hspec_label_map_t hspec_sam_classify_scene(const hspec_datacube_t *dc,
                                            const double *ref_spectra,
                                            size_t nclasses);

/* ─── API: Matched Filter ───────────────────────────────────────────── */

/**
 * Apply matched filter to a pixel given target and background statistics.
 *
 * @math MF(x) = (tᵀ·Σ⁻¹·(x-μ)) / √(tᵀ·Σ⁻¹·t)
 *
 * @param pixel        Test pixel
 * @param target       Target spectrum [nbands]
 * @param bg_mean      Background mean [nbands]
 * @param bg_inv_cov   Background inverse covariance [nbands][nbands] row-major
 * @param nbands       Number of bands
 * @param p_fa         Desired false-alarm probability
 * @return             MF detection result
 *
 * Complexity: O(nbands²)
 * Reference: Manolakis & Shaw (2002), Lincoln Lab Journal
 */
hspec_matched_filter_t hspec_matched_filter_detect(const hspec_pixel_t *pixel,
                                                     const double *target,
                                                     const double *bg_mean,
                                                     const double *bg_inv_cov,
                                                     size_t nbands, double p_fa);

/* ─── API: ACE ──────────────────────────────────────────────────────── */

/**
 * Adaptive Coherence Estimator for a single pixel.
 *
 * @math ACE(x) = (tᵀ·Σ⁻¹·(x-μ))² / [(tᵀ·Σ⁻¹·t)·((x-μ)ᵀ·Σ⁻¹·(x-μ))]
 *
 * @param pixel        Test pixel
 * @param target       Target spectrum [nbands]
 * @param bg_mean      Background mean [nbands]
 * @param bg_inv_cov   Background inverse covariance [nbands][nbands]
 * @param nbands       Number of bands
 * @param p_fa         False-alarm probability
 * @return             ACE detection result
 *
 * Complexity: O(nbands²)
 * Reference: Kraut & Scharf (1999), IEEE TSP
 */
hspec_ace_result_t hspec_ace_detect(const hspec_pixel_t *pixel,
                                     const double *target,
                                     const double *bg_mean,
                                     const double *bg_inv_cov,
                                     size_t nbands, double p_fa);

/* ─── API: RX Anomaly Detection ──────────────────────────────────── */

/**
 * Global RX anomaly detector over the entire datacube.
 *
 * @math RX(x_p) = (x_p - μ)ᵀ·Σ⁻¹·(x_p - μ)
 *
 * @param dc     Datacube
 * @param p_fa   False-alarm probability
 * @return       RX detection result (caller must hspec_rx_result_free)
 *
 * Complexity: O(n_pixels·nbands²)
 * Reference: Reed & Yu (1990)
 */
hspec_rx_result_t hspec_rx_detect_global(const hspec_datacube_t *dc, double p_fa);

void hspec_rx_result_free(hspec_rx_result_t *rx);

/* ─── API: K-Means Clustering ─────────────────────────────────────── */

/**
 * K-means clustering of hyperspectral pixels.
 *
 * @param dc          Datacube
 * @param k           Number of clusters
 * @param max_iters   Maximum iterations
 * @param tol         Convergence tolerance (relative centroid shift)
 * @return            K-means result
 *
 * Complexity: O(n_pixels·k·nbands·iterations)
 * Reference: MacQueen (1967), "Some methods for classification..."
 */
hspec_kmeans_result_t hspec_kmeans_cluster(const hspec_datacube_t *dc, size_t k,
                                            size_t max_iters, double tol);

void hspec_kmeans_result_free(hspec_kmeans_result_t *km);

/* ─── Utility Functions ────────────────────────────────────────────── */

/**
 * Estimate a detection threshold from Mahalanobis distance distribution
 * using chi-squared approximation (B degrees of freedom).
 *
 * @param nbands  Spectral dimensionality
 * @param p_fa    Desired false-alarm probability
 * @return        Threshold value
 */
double hspec_chi2_threshold(size_t nbands, double p_fa);

/**
 * Compute Mahalanobis distance between a pixel and background.
 *
 * @math d²(x, μ, Σ) = (x - μ)ᵀ·Σ⁻¹·(x - μ)
 *
 * @param x        Pixel spectrum [nbands]
 * @param mu       Background mean [nbands]
 * @param inv_cov  Inverse covariance [nbands][nbands] row-major
 * @param nbands   Number of bands
 * @return         Squared Mahalanobis distance
 */
double hspec_mahalanobis_distance(const double *x, const double *mu,
                                   const double *inv_cov, size_t nbands);

#endif /* HYPERSPECTRAL_CLASSIFICATION_H */