/**
 * @file    hyperspectral_dimensionality.h
 * @brief   Dimensionality Reduction — L3 Mathematical Structures + L5 Algorithms
 *
 * @details Implements PCA (Principal Component Analysis), MNF (Minimum Noise
 *          Fraction), ICA (Independent Component Analysis via FastICA), and
 *          NMF (Non-negative Matrix Factorization) for hyperspectral data.
 *          These methods address the "curse of dimensionality" (Hughes phenomenon)
 *          inherent to high-dimensional hyperspectral data.
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Intrinsic dimensionality, Virtual Dimensionality (VD)
 *     - Eigenvalue spectrum, explained variance ratio
 *     - Principal components (PCs), Independent components (ICs)
 *   L2 - Core Concepts:
 *     - Hughes phenomenon: classification accuracy degrades with too many bands
 *     - Data redundancy: adjacent bands are highly correlated
 *     - Dimensionality reduction for classification and unmixing
 *   L3 - Mathematical Structures:
 *     - Covariance matrix, eigendecomposition
 *     - Singular Value Decomposition (SVD)
 *     - Non-negative matrix factorization
 *     - Mutual information, negentropy
 *   L5 - Algorithms:
 *     - PCA via eigendecomposition (covariance) or SVD
 *     - Minimum Noise Fraction (MNF) transform
 *     - FastICA (fixed-point iteration, kurtosis/negentropy)
 *     - NMF via multiplicative updates (Lee-Seung)
 *     - Band selection via mutual information maximization
 *
 * Reference:
 *   - Jolliffe, "Principal Component Analysis" (2002)
 *   - Green et al., "A Transformation for Ordering Multispectral Data" (1988)
 *   - Hyvarinen & Oja, "Independent Component Analysis" (2001)
 *   - Lee & Seung, "Learning the Parts of Objects by NMF" (1999), Nature
 */

#ifndef HYPERSPECTRAL_DIMENSIONALITY_H
#define HYPERSPECTRAL_DIMENSIONALITY_H

#include "hyperspectral_core.h"

/* ─── L3/L5: PCA Result ───────────────────────────────────────────── */

/**
 * @brief Principal Component Analysis result.
 *
 * PCA finds orthogonal directions (principal components) that maximize
 * variance. PCs are eigenvectors of the covariance matrix sorted by
 * decreasing eigenvalue.
 *
 * @math X_c = X - μ   (centered)
 * @math Σ = (1/(N-1))·X_cᵀ·X_c
 * @math Σ = V·Λ·Vᵀ  (eigendecomposition)
 * @math Y = X_c·V    (projection onto PCs)
 */
typedef struct {
    size_t   nbands;              /**< Original dimensionality */
    size_t   ncomponents;         /**< Number of retained PCs */
    double  *eigenvalues;         /**< Eigenvalues λ_i descending [ncomponents] */
    double  *eigenvectors;        /**< Eigenvectors v_i [ncomponents][nbands] row-major */
    double  *explained_variance;  /**< λ_i / Σλ_j per component */
    double  *cumulative_variance; /**< Cumulative explained variance ratio */
    double  *projected_data;      /**< Projected data [npixels][ncomponents] row-major */
    double  *mean_spectrum;       /**< Mean spectrum μ [nbands] */
    size_t   npixels;
    double   total_variance;      /**< Sum of all eigenvalues */
} hspec_pca_result_t;

/* ─── L5: MNF Result ─────────────────────────────────────────────── */

/**
 * @brief Minimum Noise Fraction (MNF) transform result.
 *
 * MNF is a noise-adjusted PCA that orders components by SNR rather
 * than variance. It requires noise covariance estimation (e.g., from
 * dark-current measurements or spatial differencing).
 *
 * Steps:
 *   1. Estimate noise covariance Σ_n
 *   2. Whiten noise: Z = Σ_n^(-1/2)·X
 *   3. PCA on Z yields MNF components
 *   4. Inverse transform for physical interpretation
 *
 * MNF components are ordered by decreasing signal-to-noise ratio.
 */
typedef struct {
    size_t   nbands;
    size_t   ncomponents;
    double  *eigenvalues;         /**< MNF eigenvalues (λ_signal/λ_noise) */
    double  *eigenvectors;        /**< MNF eigenvectors [ncomponents][nbands] */
    double  *snr_per_component;   /**< SNR estimate per MNF component */
    double  *projected_data;      /**< MNF-transformed data [npixels][ncomponents] */
    double  *noise_fraction;      /**< 1/(1+SNR) per component */
    size_t   npixels;
    int      noise_estimated;     /**< 1 = noise cov estimated, 0 = identity assumed */
} hspec_mnf_result_t;

/* ─── L5: ICA Result ──────────────────────────────────────────────── */

/**
 * @brief Independent Component Analysis (FastICA) result.
 *
 * ICA finds statistically independent components (as opposed to merely
 * uncorrelated in PCA). Uses negentropy or kurtosis as the independence
 * measure via fixed-point iteration.
 *
 * @math X = A·S  (mixing model)
 * @math S = W·X  (demixing, W = A⁻¹)
 */
typedef struct {
    size_t   ncomponents;         /**< Number of independent components */
    size_t   nbands;              /**< Original bands */
    double  *mixing_matrix;       /**< A [nbands][ncomponents] */
    double  *unmixing_matrix;     /**< W = pinv(A) [ncomponents][nbands] */
    double  *components;          /**< S = W·X [npixels][ncomponents] */
    size_t   npixels;
    size_t   niterations;         /**< Number of FastICA iterations */
    int      converged;           /**< 1 = converged */
    double   negentropy;          /**< Total estimated negentropy */
} hspec_ica_result_t;

/* ─── L5: NMF Result ──────────────────────────────────────────────── */

/**
 * @brief Non-negative Matrix Factorization result.
 *
 * NMF decomposes X ≈ W·H where X ≥ 0, W ≥ 0, H ≥ 0.
 * In hyperspectral context:
 *   - W are spectral basis vectors (endmember-like)
 *   - H are abundance coefficients
 *
 * Multiplicative update rules:
 *   H ← H ⊙ (Wᵀ·X) ⊘ (Wᵀ·W·H + ε)
 *   W ← W ⊙ (X·Hᵀ) ⊘ (W·H·Hᵀ + ε)
 */
typedef struct {
    size_t   nbands;
    size_t   ncomponents;         /**< Rank of factorization (≥ 1) */
    size_t   npixels;
    double  *W;                   /**< Basis matrix [nbands][ncomponents] */
    double  *H;                   /**< Coefficient matrix [ncomponents][npixels] */
    double   reconstruction_error; /**< Frobenius norm ||X - W·H||_F */
    size_t   niterations;
    int      converged;
} hspec_nmf_result_t;

/* ─── L3: Eigenvalue Spectrum for VD Estimation ──────────────────── */

/**
 * @brief Complete eigenvalue analysis of the covariance matrix.
 */
typedef struct {
    size_t   nbands;
    double  *eigenvalues;         /**< All eigenvalues [nbands] descending */
    double  *eigenvectors;        /**< Full eigendecomposition [nbands][nbands] */
    double   trace;               /**< Sum of eigenvalues (total variance) */
    double   condition_number;    /**< λ_max / λ_min */
    size_t   effective_rank;      /**< Number of eigenvalues > tolerance */
} hspec_eigenvalue_spectrum_t;

/* ─── API: PCA ────────────────────────────────────────────────────── */

/**
 * PCA via covariance eigendecomposition using QR iteration.
 *
 * @param dc          Datacube
 * @param ncomponents Number of PCs to retain (0 = use 99% variance threshold)
 * @return            PCA result (caller must hspec_pca_result_free)
 *
 * Complexity: O(n_pixels·nbands² + nbands³)
 * Reference: Jolliffe (2002)
 */
hspec_pca_result_t hspec_pca_compute(const hspec_datacube_t *dc, size_t ncomponents);
void hspec_pca_result_free(hspec_pca_result_t *pca);

/**
 * Reconstruct data from PCA projection.
 *
 * @param pca         PCA result
 * @param ncomponents Number of components to use
 * @param data_out    Reconstructed data [npixels][nbands] row-major
 * @return            0 on success
 */
int hspec_pca_reconstruct(const hspec_pca_result_t *pca, size_t ncomponents,
                           double *data_out);

/* ─── API: MNF ────────────────────────────────────────────────────── */

/**
 * Minimum Noise Fraction transform.
 *
 * @param dc         Datacube
 * @param noise_cov  Noise covariance estimate [nbands][nbands] (NULL = identity)
 * @param ncomponents Number of MNF components to retain
 * @return           MNF result
 *
 * Complexity: O(n_pixels·nbands² + nbands³)
 * Reference: Green et al. (1988), Remote Sensing of Environment
 */
hspec_mnf_result_t hspec_mnf_compute(const hspec_datacube_t *dc,
                                      const double *noise_cov,
                                      size_t ncomponents);
void hspec_mnf_result_free(hspec_mnf_result_t *mnf);

/* ─── API: ICA ────────────────────────────────────────────────────── */

/**
 * FastICA algorithm with deflationary orthogonalization.
 *
 * Uses negentropy approximation via tanh or kurtosis contrast function.
 *
 * @param dc           Datacube (should be centered and whitened)
 * @param ncomponents  Number of ICs to extract
 * @param max_iters    Maximum iterations per component
 * @param tol          Convergence tolerance
 * @return             ICA result
 *
 * Complexity: O(n_pixels·nbands·ncomponents·n_iters)
 * Reference: Hyvarinen (1999), IEEE TNN
 */
hspec_ica_result_t hspec_fastica(const hspec_datacube_t *dc, size_t ncomponents,
                                  size_t max_iters, double tol);
void hspec_ica_result_free(hspec_ica_result_t *ica);

/* ─── API: NMF ─────────────────────────────────────────────────── */

/**
 * NMF via multiplicative updates (Lee-Seung algorithm).
 *
 * @param dc             Datacube (all values must be ≥ 0)
 * @param ncomponents    Factorization rank
 * @param max_iterations Maximum iterations
 * @param tol            Convergence tolerance
 * @return               NMF result
 *
 * Complexity: O(n_pixels·nbands·ncomponents·n_iters)
 * Reference: Lee & Seung (2001), NIPS
 */
hspec_nmf_result_t hspec_nmf_factorize(const hspec_datacube_t *dc,
                                        size_t ncomponents,
                                        size_t max_iterations, double tol);
void hspec_nmf_result_free(hspec_nmf_result_t *nmf);

/* ─── API: Eigendecomposition ─────────────────────────────────⟧ */

/**
 * Symmetric matrix eigendecomposition via QR iteration.
 *
 * @param A            Symmetric matrix [n][n] row-major
 * @param n            Dimension
 * @param max_iters    Maximum QR iterations
 * @param tol          Convergence tolerance
 * @param eigenvalues  Output eigenvalues [n] descending
 * @param eigenvectors Output eigenvectors [n][n] row-major
 * @return             0 on success
 *
 * Complexity: O(n³·n_iters)
 */
int hspec_symmetric_eigen(const double *A, size_t n, size_t max_iters,
                           double tol, double *eigenvalues,
                           double *eigenvectors);

/**
 * Compute the full eigenvalue spectrum from a datacube covariance.
 *
 * @param dc       Datacube
 * @param max_iters Max QR iterations
 * @param tol      Convergence tolerance
 * @return         Eigenvalue spectrum
 */
hspec_eigenvalue_spectrum_t hspec_eigenvalue_spectrum_compute(
    const hspec_datacube_t *dc, size_t max_iters, double tol);
void hspec_eigenvalue_spectrum_free(hspec_eigenvalue_spectrum_t *es);

/* ─── Band Selection (L3) ─────────────────────────────────────── */

/**
 * Select the k most informative bands via mutual information maximization.
 *
 * Uses forward greedy search: starting from empty set, iteratively add
 * the band that maximizes mutual information with class labels.
 *
 * @param dc           Datacube
 * @param labels       Ground truth labels [npixels] (can pass NULL for unsupervised)
 * @param npixels      Number of pixels
 * @param k            Number of bands to select
 * @param selected_out Output selected band indices [k]
 * @return             0 on success
 *
 * Complexity: O(k·nbands·npixels)
 */
int hspec_band_selection_mi(const hspec_datacube_t *dc, const size_t *labels,
                             size_t npixels, size_t k, size_t *selected_out);

#endif /* HYPERSPECTRAL_DIMENSIONALITY_H */