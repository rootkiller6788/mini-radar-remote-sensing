/**
 * @file    hyperspectral_unmixing.h
 * @brief   Spectral Unmixing — L2 Core Concepts + L5 Algorithms + L6 Canonical Problems
 *
 * @details Spectral unmixing decomposes each pixel's spectrum into a set of
 *          pure material spectra (endmembers) and their fractional abundances.
 *          Covers linear mixing model (LMM), bilinear models, endmember
 *          extraction (PPI, N-FINDR, VCA), and abundance estimation (FCLS, NCLS).
 *
 * Knowledge Mapping:
 *   L1 - Definitions:
 *     - Endmember: pure spectral signature of a single material
 *     - Abundance: fractional area coverage of each endmember in a pixel
 *     - Mixing model: mathematical relationship between endmembers and pixel
 *   L2 - Core Concepts:
 *     - Linear mixing: photons interact with single material before reaching sensor
 *     - Intimate (nonlinear) mixing: multiple scattering, grain mixtures
 *     - Additivity constraint: abundances sum to 1
 *     - Non-negativity constraint: abundances >= 0
 *   L3 - Mathematical Structures:
 *     - Convex simplex geometry: pixel spectra lie within endmember convex hull
 *     - Matrix factorizations for unmixing
 *     - Constrained least squares optimization
 *   L5 - Algorithms:
 *     - Vertex Component Analysis (VCA) — endmember extraction
 *     - N-FINDR — maximum simplex volume
 *     - Pixel Purity Index (PPI) — convex geometry
 *     - Fully Constrained Least Squares (FCLS) — abundance estimation
 *     - Non-negatively Constrained Least Squares (NCLS)
 *   L6 - Canonical Problems:
 *     - Subpixel mineral mapping
 *     - Land cover fraction estimation
 *     - Urban impervious surface mapping
 *
 * Reference:
 *   - Keshava & Mustard, "Spectral Unmixing" (2002), IEEE Signal Processing Mag
 *   - Nascimento & Dias, "Vertex Component Analysis" (2005), IEEE TGRS
 *   - Heinz & Chang, "Fully Constrained Least Squares" (2001), IEEE TGRS
 */

#ifndef HYPERSPECTRAL_UNMIXING_H
#define HYPERSPECTRAL_UNMIXING_H

#include "hyperspectral_core.h"

/* ─── L1: Endmember ────────────────────────────────────────────────── */

/**
 * @brief An endmember — pure material spectral signature.
 *
 * Endmembers represent the spectral reflectance of a homogeneous
 * surface material (e.g., kaolinite, green vegetation, water).
 * They form the vertices of the data simplex in B-dimensional
 * spectral space.
 */
typedef struct {
    size_t   nbands;
    double  *spectrum;           /**< Reflectance [nbands] */
    char     material_name[128];
    size_t   source_pixel_row;   /**< Originating pixel (if from image) */
    size_t   source_pixel_col;
    double   purity_index;       /**< Estimated purity score */
} hspec_endmember_t;

/* ─── L1: Abundance Map ────────────────────────────────────────────── */

/**
 * @brief Fractional abundance map for one endmember.
 *
 * For each pixel, abundance ∈ [0, 1] and Σ a_i = 1 (sum-to-one).
 */
typedef struct {
    size_t   nrows;
    size_t   ncols;
    size_t   endmember_index;
    double  *abundance;          /**< Fractional abundance per pixel [nrows*ncols] */
    double   mean_abundance;     /**< Mean abundance over the scene */
    double   std_abundance;      /**< Standard deviation of abundance */
} hspec_abundance_map_t;

/* ─── L2: Mixing Model Result ──────────────────────────────────────── */

/**
 * @brief Complete unmixing result for a hyperspectral scene.
 */
typedef struct {
    size_t   n_endmembers;
    hspec_endmember_t *endmembers; /**< Extracted endmembers [n_endmembers] */
    size_t   nbands;
    size_t   nrows;
    size_t   ncols;
    double  *abundances;         /**< Abundance maps [n_endmembers][nrows*ncols] row-major */
    double  *reconstruction;     /**< Reconstructed pixel spectra [nrows*ncols][nbands] */
    double  *residual;           /**< Residual error per pixel per band */
    double   rmse;               /**< Root-mean-square reconstruction error */
    double   mean_sad;           /**< Mean spectral angle between original and reconstructed */
    double   sum_to_one_error;   /**< Mean |Σa - 1| constraint violation */
    double   nonneg_error;       /**< Mean |max(-a, 0)| constraint violation */
    int      fully_constrained;  /**< 1 = FCLS, 0 = unconstrained/partially */
} hspec_unmixing_result_t;

/* ─── L2: Mixing Model Types ──────────────────────────────────────── */

typedef enum {
    HSPEC_MIXING_LINEAR       = 0,  /**< x = M·a + n (standard LMM) */
    HSPEC_MIXING_FAN          = 1,  /**< Bilinear: x = M·a + Σ_{i<j} a_i·a_j·m_i⊙m_j */
    HSPEC_MIXING_NASCIMENTO   = 2,  /**< Bilinear: Nascimento model */
    HSPEC_MIXING_GBM          = 3,  /**< Generalized Bilinear Model */
    HSPEC_MIXING_PPNM         = 4   /**< Polynomial Post-Nonlinear Model */
} hspec_mixing_model_t;

/* ─── L5: VCA (Vertex Component Analysis) Parameters ───────────────── */

typedef struct {
    size_t   n_endmembers;       /**< Number of endmembers to extract */
    size_t   max_iterations;     /**< Maximum iterations */
    double   convergence_tol;    /**< Convergence tolerance */
    int      use_snr_proj;       /**< 1 = project using SNR, 0 = skip */
} hspec_vca_params_t;

/* ─── L5: N-FINDR Parameters ─────────────────────────────────────── */

typedef struct {
    size_t   n_endmembers;       /**< Number of endmembers */
    size_t   max_iterations;     /**< Maximum iterations */
    size_t   n_random_starts;    /**< Number of random initializations */
} hspec_nfindr_params_t;

/* ─── L5: FCLS Parameters ─────────────────────────────────────────── */

typedef struct {
    size_t   max_iterations;     /**< Maximum iterations for NNLS solver */
    double   tolerance;          /**< Convergence tolerance */
    double   regularization;     /**< Tikhonov regularization parameter λ */
    int      enforce_sum_to_one; /**< 1 = FCLS, 0 = NCLS only */
} hspec_fcls_params_t;

/* ─── API: Endmember Extraction ────────────────────────────────────── */

/**
 * Vertex Component Analysis — extracts endmembers by projecting data
 * onto directions orthogonal to previously extracted endmembers.
 *
 * @param dc       Datacube
 * @param params   VCA parameters
 * @return         Unmixing result with extracted endmembers (caller must free)
 *
 * @math VCA iteratively finds extreme points in the spectral simplex
 *       by projecting onto directions orthogonal to the already-identified
 *       endmember subspace.
 *
 * Complexity: O(n_pixels·nbands·n_endmembers)
 * Reference: Nascimento & Dias (2005), IEEE TGRS
 */
hspec_unmixing_result_t hspec_vca_extract(const hspec_datacube_t *dc,
                                           hspec_vca_params_t params);

/**
 * N-FINDR — finds endmembers by maximizing the simplex volume.
 *
 * @param dc       Datacube
 * @param params   N-FINDR parameters
 * @return         Unmixing result with extracted endmembers
 *
 * Complexity: O(n_pixels·n_endmembers²·nbands·n_iterations)
 * Reference: Winter (1999), Proc. SPIE
 */
hspec_unmixing_result_t hspec_nfindr_extract(const hspec_datacube_t *dc,
                                              hspec_nfindr_params_t params);

/**
 * Pixel Purity Index — ranks pixels by "extremeness" via random projections.
 *
 * @param dc           Datacube
 * @param n_projections Number of random projections
 * @param purity_scores_out Output purity score per pixel [npixels]
 * @return              0 on success
 *
 * Complexity: O(n_projections·n_pixels·nbands)
 * Reference: Boardman et al. (1995), JPL
 */
int hspec_ppi_compute(const hspec_datacube_t *dc, size_t n_projections,
                       double *purity_scores_out);

/* ─── API: Abundance Estimation ──────────────────────────────────── */

/**
 * Fully Constrained Least Squares (FCLS) abundance estimation.
 *
 * Solves: min ||M·a - x||²  subject to  a ≥ 0,  Σ a_i = 1
 *
 * @param x          Pixel spectrum [nbands]
 * @param E          Endmember matrix [n_endmembers][nbands] (each row is an EM)
 * @param n_endmembers Number of endmembers
 * @param nbands     Number of spectral bands
 * @param params     FCLS parameters
 * @param abundance_out Output abundance vector [n_endmembers]
 * @return            0 on success, -1 on error
 *
 * Complexity: O(n_endmembers²·nbands·n_iterations)
 * Reference: Heinz & Chang (2001), IEEE TGRS
 */
int hspec_fcls_estimate(const double *x, const double *E,
                         size_t n_endmembers, size_t nbands,
                         hspec_fcls_params_t params, double *abundance_out);

/**
 * Non-negatively Constrained Least Squares (NCLS) — no sum-to-one.
 *
 * Solves: min ||M·a - x||²  subject to  a ≥ 0
 *
 * @param x          Pixel spectrum [nbands]
 * @param E          Endmember matrix [n_endmembers][nbands] row-major
 * @param n_endmembers Number of endmembers
 * @param nbands     Number of spectral bands
 * @param max_iters  Maximum iterations for active-set solver
 * @param tol        Convergence tolerance
 * @param abundance_out Output abundance vector [n_endmembers]
 * @return            0 on success
 *
 * Complexity: O(n_endmembers²·nbands·n_max_iters)
 * Reference: Lawson & Hanson (1974), "Solving Least Squares Problems"
 */
int hspec_ncls_estimate(const double *x, const double *E,
                         size_t n_endmembers, size_t nbands,
                         size_t max_iters, double tol, double *abundance_out);

/* ─── API: Full Scene Unmixing ────────────────────────────────────── */

/**
 * Unmix the entire datacube using extracted endmembers and FCLS.
 *
 * @param dc           Datacube
 * @param endmembers   Endmember array [n_endmembers]
 * @param n_endmembers Number of endmembers
 * @param params       FCLS parameters
 * @return             Complete unmixing result
 *
 * Complexity: O(n_pixels·n_endmembers²·nbands·n_iterations)
 */
hspec_unmixing_result_t hspec_unmix_scene(const hspec_datacube_t *dc,
                                           const hspec_endmember_t *endmembers,
                                           size_t n_endmembers,
                                           hspec_fcls_params_t params);

/* ─── Memory Management ───────────────────────────────────────────── */

void hspec_unmixing_result_free(hspec_unmixing_result_t *ur);
void hspec_abundance_map_free(hspec_abundance_map_t *am);

/* ─── Bilinear Mixing Models (L8 Advanced) ─────────────────────────── */

/**
 * Fan bilinear mixing model reconstruction.
 *
 * @math x = M·a + Σ_{i<j} γ_{ij}·a_i·a_j·(m_i ⊙ m_j)
 *
 * @param E            Endmember matrix [n_endmembers][nbands]
 * @param a            Abundances [n_endmembers]
 * @param gamma        Bilinear coefficients (set to 1.0 for standard Fan)
 * @param n_endmembers Number of endmembers
 * @param nbands       Number of bands
 * @param x_out        Output mixed spectrum [nbands]
 * @return             0 on success
 *
 * Reference: Fan et al. (2009), IEEE TGRS
 */
int hspec_fan_mixing_model(const double *E, const double *a,
                            const double *gamma, size_t n_endmembers,
                            size_t nbands, double *x_out);

/**
 * Hapke intimate mixing model for particulate surfaces.
 *
 * Models multiple scattering in granular media where particles
 * are in intimate contact (nonlinear mixing).
 *
 * @math R = w/(4) · μ₀/(μ₀+μ) · [(1+B(g))·p(g) + H(μ₀)·H(μ) - 1]
 *
 * @param ssa          Single-scattering albedo per band [nbands]
 * @param nbands       Number of bands
 * @param mu0_cos      cos(solar zenith)
 * @param mu_cos       cos(viewing zenith)
 * @param refl_out     Output reflectance [nbands]
 * @return             0 on success
 *
 * Complexity: O(nbands)
 * Reference: Hapke (2012), "Theory of Reflectance and Emittance Spectroscopy"
 */
int hspec_hapke_reflectance(const double *ssa, size_t nbands,
                             double mu0_cos, double mu_cos,
                             double *refl_out);

#endif /* HYPERSPECTRAL_UNMIXING_H */