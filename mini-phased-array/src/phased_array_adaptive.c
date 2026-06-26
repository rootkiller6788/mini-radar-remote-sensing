/* ============================================================================
 * phased_array_adaptive.c — Adaptive Beamforming Algorithms
 *
 * Implements adaptive beamforming for interference rejection and
 * spatial filtering:
 *   - Sample Matrix Inversion (SMI) for covariance estimation
 *   - MVDR (Minimum Variance Distortionless Response) / Capon beamformer
 *   - LCMV (Linearly Constrained Minimum Variance) beamformer
 *   - Capon spatial power spectrum estimation
 *   - Hermitian matrix utilities (complex conjugate transpose, inverse)
 *
 * Adaptive beamforming automatically adjusts weights based on the
 * received data statistics to suppress interference while preserving
 * the desired signal from the look direction.
 *
 * Mathematical foundation (Van Trees §6):
 *
 * MVDR optimization problem:
 *   minimize  w^H R w     (minimize total output power)
 *   subject to  w^H a(theta_s) = 1   (unity gain in look direction)
 *
 * Solution: w_MVDR = R^{-1} a / (a^H R^{-1} a)
 *
 * where R is the spatial covariance matrix of the received data:
 *   R = E{x x^H} = sigma_s^2 a a^H + sum_{j} sigma_j^2 v_j v_j^H + sigma_n^2 I
 *
 * The MVDR beamformer places nulls in the directions of interferers
 * while maintaining a fixed response in the look direction.
 * Null depth depends on interferer power: stronger interferers
 * get deeper nulls (since they contribute more to the cost function).
 *
 * References:
 *   Van Trees, H.L. (2002) "Optimum Array Processing", Part IV.
 *   Capon, J. (1969) "High-Resolution Frequency-Wavenumber Spectrum
 *     Analysis", Proc. IEEE, vol. 57.
 *   Monzingo & Miller (2004) "Introduction to Adaptive Arrays".
 * ============================================================================ */

#include "phased_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* ============================================================================
 * L8: Complex Matrix Utilities — Hermitian Operations
 *
 * For an M×M complex matrix:
 *   - Conjugate transpose (Hermitian): A^H_{ij} = conj(A_{ji})
 *   - Matrix-vector multiply: y = A x
 *   - Inner product: a^H b = Σ conj(a_i) b_i
 *   - Outer product: a b^H = matrix with (i,j) = a_i conj(b_j)
 * ============================================================================ */

/**
 * Conjugate transpose of an M×M matrix.
 * B = A^H where B[i][j] = conj(A[j][i]).
 *
 * @param A Input matrix (M×M, row-major).
 * @param B Output matrix (M×M, row-major). May equal A (in-place not supported).
 * @param M Matrix dimension.
 */
/**
 * Matrix-vector multiply: y = A x.
 *
 * @param A M×M complex matrix (row-major).
 * @param x M-element complex vector.
 * @param y Output: M-element complex vector.
 * @param M Dimension.
 */
static void mat_vec_mul(const double complex *A,
                         const double complex *x,
                         double complex *y, uint32_t M)
{
    for (uint32_t i = 0; i < M; i++) {
        y[i] = 0.0 + 0.0 * I;
        for (uint32_t j = 0; j < M; j++) {
            y[i] += A[i * M + j] * x[j];
        }
    }
}

/**
 * Inner (dot) product of two complex vectors: a^H b = Σ conj(a_i) * b_i.
 */
static double complex inner_product(const double complex *a,
                                     const double complex *b,
                                     uint32_t M)
{
    double complex result = 0.0 + 0.0 * I;
    for (uint32_t i = 0; i < M; i++) {
        result += conj(a[i]) * b[i];
    }
    return result;
}

/* ============================================================================
 * L8: Gauss-Jordan Matrix Inversion for Complex Matrices
 *
 * Inverts an M×M complex matrix using Gauss-Jordan elimination with
 * partial pivoting. This is numerically stable for moderate condition
 * numbers (κ < 10⁶). For ill-conditioned matrices, diagonal loading
 * is applied externally.
 *
 * Algorithm complexity: O(M³).
 *
 * For M > 100, consider using Cholesky decomposition (O(M³/3)) for
 * Hermitian positive-definite matrices (covariance matrices are HPD).
 * ============================================================================ */

/**
 * Invert an M×M complex matrix using Gauss-Jordan with partial pivoting.
 *
 * @param A Input matrix (M×M, row-major). Content is destroyed.
 * @param A_inv Output: inverse matrix (M×M, row-major).
 * @param M Matrix dimension.
 * @return 0 on success, -1 if singular.
 */
static int matrix_invert(double complex *A, double complex *A_inv, uint32_t M)
{
    if (M == 0 || !A || !A_inv) return -1;
    if (M == 1) {
        if (cabs(A[0]) < 1e-15) return -1;
        A_inv[0] = 1.0 / A[0];
        return 0;
    }

    /* Initialize A_inv as identity matrix */
    for (uint32_t i = 0; i < M; i++) {
        for (uint32_t j = 0; j < M; j++) {
            A_inv[i * M + j] = (i == j) ? (1.0 + 0.0 * I) : (0.0 + 0.0 * I);
        }
    }

    /* Gauss-Jordan elimination */
    for (uint32_t col = 0; col < M; col++) {
        /* Find pivot (largest magnitude in column col, rows col..M-1) */
        uint32_t pivot_row = col;
        double pivot_mag = cabs(A[col * M + col]);
        for (uint32_t row = col + 1; row < M; row++) {
            double mag = cabs(A[row * M + col]);
            if (mag > pivot_mag) {
                pivot_mag = mag;
                pivot_row = row;
            }
        }

        if (pivot_mag < 1e-15) {
            return -1;  /* Singular matrix */
        }

        /* Swap rows if needed */
        if (pivot_row != col) {
            for (uint32_t j = 0; j < M; j++) {
                double complex tmp = A[col * M + j];
                A[col * M + j] = A[pivot_row * M + j];
                A[pivot_row * M + j] = tmp;

                tmp = A_inv[col * M + j];
                A_inv[col * M + j] = A_inv[pivot_row * M + j];
                A_inv[pivot_row * M + j] = tmp;
            }
        }

        /* Normalize pivot row */
        double complex pivot = A[col * M + col];
        for (uint32_t j = 0; j < M; j++) {
            A[col * M + j] /= pivot;
            A_inv[col * M + j] /= pivot;
        }

        /* Eliminate all other rows */
        for (uint32_t row = 0; row < M; row++) {
            if (row == col) continue;
            double complex factor = A[row * M + col];
            for (uint32_t j = 0; j < M; j++) {
                A[row * M + j] -= factor * A[col * M + j];
                A_inv[row * M + j] -= factor * A_inv[col * M + j];
            }
        }
    }

    return 0;
}

/* ============================================================================
 * L8: SMI Beamformer Initialization and Covariance Estimation
 * ============================================================================ */

pa_smi_beamformer_t *pa_smi_beamformer_init(uint32_t num_elements,
                                             uint32_t num_snapshots)
{
    if (num_elements == 0) return NULL;

    pa_smi_beamformer_t *bf = (pa_smi_beamformer_t *)calloc(1, sizeof(*bf));
    if (!bf) return NULL;

    bf->num_elements   = num_elements;
    bf->num_snapshots  = num_snapshots;
    bf->forgetting_factor = 1.0;    /* No forgetting (block processing) */
    bf->diagonal_loading  = 0.01;   /* 1% of trace for regularization */

    bf->weight_vector    = (double complex *)calloc(num_elements, sizeof(double complex));
    bf->steering_vector  = (double complex *)calloc(num_elements, sizeof(double complex));
    bf->cov_matrix       = (double complex *)calloc((size_t)num_elements * num_elements,
                                                      sizeof(double complex));
    bf->cov_matrix_inv   = (double complex *)calloc((size_t)num_elements * num_elements,
                                                      sizeof(double complex));

    if (!bf->weight_vector || !bf->steering_vector
        || !bf->cov_matrix || !bf->cov_matrix_inv) {
        pa_smi_beamformer_free(bf);
        return NULL;
    }

    /* Initialize weight vector to uniform (delay-and-sum for broadside) */
    for (uint32_t i = 0; i < num_elements; i++) {
        bf->weight_vector[i] = 1.0 / sqrt((double)num_elements);
    }

    return bf;
}

void pa_smi_beamformer_free(pa_smi_beamformer_t *bf)
{
    if (!bf) return;
    free(bf->weight_vector);
    free(bf->steering_vector);
    free(bf->cov_matrix);
    free(bf->cov_matrix_inv);
    free(bf);
}

/**
 * Estimate the sample covariance matrix from snapshot data.
 *
 * R̂ = (1/K) Σ_{k=0}^{K-1} x[k] x[k]^H
 *
 * This is the Maximum Likelihood estimate of the covariance matrix
 * when the snapshots are i.i.d. zero-mean complex Gaussian.
 *
 * For convergence, the number of snapshots K should satisfy:
 *   K ≥ 2M  (to ensure invertibility)
 * where M is the number of array elements.
 *
 * Diagonal loading:
 *   R̂_loaded = R̂ + γ · (tr(R̂)/M) · I
 *
 * This regularizes the inverse for poorly conditioned matrices
 * (e.g., when K < 2M or when sources are closely spaced).
 * γ is the diagonal loading factor (typical: 0.001–0.1).
 *
 * @param bf Beamformer state (updated in place).
 * @param snapshots 2D array: snapshots[k][m] where k = 0..K-1,
 *                   m = 0..M-1 (element index).
 */
void pa_smi_estimate_covariance(pa_smi_beamformer_t *bf,
                                 const double complex **snapshots)
{
    if (!bf || !snapshots || bf->num_elements == 0
        || bf->num_snapshots == 0) return;

    uint32_t M = bf->num_elements;
    uint32_t K = bf->num_snapshots;

    /* Zero out covariance matrix */
    memset(bf->cov_matrix, 0,
           (size_t)M * M * sizeof(double complex));

    /* Accumulate outer products: R = Σ x[k] x[k]^H */
    for (uint32_t k = 0; k < K; k++) {
        for (uint32_t i = 0; i < M; i++) {
            for (uint32_t j = 0; j < M; j++) {
                bf->cov_matrix[i * M + j] +=
                    snapshots[k][i] * conj(snapshots[k][j]);
            }
        }
    }

    /* Normalize by number of snapshots */
    double inv_K = 1.0 / (double)K;
    for (size_t idx = 0; idx < (size_t)M * M; idx++) {
        bf->cov_matrix[idx] *= inv_K;
    }

    /* Diagonal loading for regularization */
    if (bf->diagonal_loading > 0.0) {
        /* Compute trace: tr(R) = Σ R_{ii} */
        double trace = 0.0;
        for (uint32_t i = 0; i < M; i++) {
            trace += creal(bf->cov_matrix[i * M + i]);
        }

        double loading = bf->diagonal_loading * trace / (double)M;

        /* Add gamma * I to the diagonal */
        for (uint32_t i = 0; i < M; i++) {
            bf->cov_matrix[i * M + i] += loading;
        }
    }
}

/* ============================================================================
 * L8: MVDR (Capon) Beamformer
 *
 * The MVDR beamformer solves:
 *   minimize  w^H R w
 *   s.t.      w^H a = 1
 *
 * Solution (closed form):
 *   w_MVDR = R^{-1} a / (a^H R^{-1} a)
 *
 * The denominator a^H R^{-1} a is a real scalar (since R is Hermitian).
 *
 * Output power:
 *   P_MVDR = w^H R w = 1 / (a^H R^{-1} a)
 *
 * Capon spatial spectrum:
 *   P_Capon(theta) = 1 / (a(theta)^H R^{-1} a(theta))
 *
 * This is not a true power spectrum (it does not integrate to total power),
 * but it gives high-resolution DOA estimates: peaks correspond to source
 * directions, and the peak height is proportional to source power.
 *
 * SINR optimality:
 *   The MVDR beamformer maximizes the output SINR among all beamformers
 *   that preserve the look-direction response.
 *
 *   SINR_opt = sigma_s^2 · a^H R_{i+n}^{-1} a
 *   where R_{i+n} is the interference+noise covariance.
 *
 * The number of interferers that can be nulled is at most M-1
 * (M elements gives M-1 degrees of freedom beyond the look-direction
 * constraint).
 * ============================================================================ */

pa_mvdr_beamformer_t *pa_mvdr_beamformer_init(uint32_t num_elements,
                                               uint32_t num_snapshots)
{
    pa_mvdr_beamformer_t *mvdr = (pa_mvdr_beamformer_t *)calloc(1, sizeof(*mvdr));
    if (!mvdr) return NULL;

    /* Initialize the SMI base */
    pa_smi_beamformer_t *base = pa_smi_beamformer_init(num_elements, num_snapshots);
    if (!base) {
        free(mvdr);
        return NULL;
    }
    mvdr->base = *base;
    free(base);  /* Free the temp, we've copied the contents */

    return mvdr;
}

void pa_mvdr_beamformer_free(pa_mvdr_beamformer_t *bf)
{
    if (!bf) return;
    /* Free the internal arrays (allocated via pa_smi_beamformer_init) */
    free(bf->base.weight_vector);
    free(bf->base.steering_vector);
    free(bf->base.cov_matrix);
    free(bf->base.cov_matrix_inv);
    free(bf->angular_spectrum);
    free(bf);
}

/**
 * Compute MVDR weight vector given the array geometry and look direction.
 *
 * Steps:
 *   1. Compute steering vector a(theta_s, phi_s)
 *   2. Compute R^{-1} (invert covariance matrix)
 *   3. Compute a^H R^{-1} a (denominator)
 *   4. w = R^{-1} a / (a^H R^{-1} a)
 *
 * @param bf MVDR beamformer state (updated in place).
 * @param config Array geometry.
 * @param elements Array elements.
 * @param theta_steer Look direction theta (rad).
 * @param phi_steer Look direction phi (rad).
 */
void pa_mvdr_compute_weights(pa_mvdr_beamformer_t *bf,
                              const pa_array_config_t *config,
                              const pa_element_t *elements,
                              double theta_steer, double phi_steer)
{
    if (!bf || !config || !elements) return;

    uint32_t M = bf->base.num_elements;
    if (M == 0) return;

    /* Compute steering vector for the look direction */
    double complex *sv = pa_steering_vector(config, elements,
                                             theta_steer, phi_steer);
    if (!sv) return;

    /* Store steering vector */
    memcpy(bf->base.steering_vector, sv, M * sizeof(double complex));
    free(sv);

    /* Invert covariance matrix: R^{-1} */
    double complex *R_copy = (double complex *)malloc(
        (size_t)M * M * sizeof(double complex));
    if (!R_copy) return;

    memcpy(R_copy, bf->base.cov_matrix,
           (size_t)M * M * sizeof(double complex));

    int ret = matrix_invert(R_copy, bf->base.cov_matrix_inv, M);
    free(R_copy);

    if (ret != 0) {
        /* Inverse failed (singular matrix). Use identity as fallback. */
        for (uint32_t i = 0; i < M; i++) {
            for (uint32_t j = 0; j < M; j++) {
                bf->base.cov_matrix_inv[i * M + j] =
                    (i == j) ? (1.0 + 0.0 * I) : (0.0 + 0.0 * I);
            }
        }
    }

    /* Compute denominator: d = a^H R^{-1} a
     * First compute temp = R^{-1} a */
    double complex *temp = (double complex *)calloc(M, sizeof(double complex));
    if (!temp) return;

    mat_vec_mul(bf->base.cov_matrix_inv, bf->base.steering_vector, temp, M);

    /* d = a^H · temp = Σ conj(a_i) * temp_i */
    double complex d = inner_product(bf->base.steering_vector, temp, M);

    /* Compute weight vector: w = R^{-1} a / d = temp / d */
    if (cabs(d) < 1e-15) {
        /* Division by zero: fall back to steering vector */
        memcpy(bf->base.weight_vector, bf->base.steering_vector,
               M * sizeof(double complex));
    } else {
        for (uint32_t i = 0; i < M; i++) {
            bf->base.weight_vector[i] = temp[i] / d;
        }
    }

    /* Compute output power: P = 1 / (a^H R^{-1} a) = 1 / d */
    bf->output_power = (cabs(d) > 1e-15) ? (1.0 / creal(d)) : 0.0;

    free(temp);
}

/**
 * Compute the Capon spatial power spectrum.
 *
 * P_Capon(theta, phi) = 1 / (a(theta, phi)^H R^{-1} a(theta, phi))
 *
 * Evaluated on a grid of angles. The peaks indicate source DOAs.
 * The Capon spectrum is not a power spectral density in the strict
 * sense, but it provides higher angular resolution than conventional
 * (delay-and-sum) beamforming.
 *
 * @param bf MVDR beamformer (must have computed R^{-1}).
 * @param config Array geometry.
 * @param elements Array elements.
 * @param spectrum Output: P_Capon values (size = num_angles).
 * @param theta_vals Array of theta angles.
 * @param phi_vals Array of phi angles.
 * @param num_angles Number of angular bins.
 */
void pa_mvdr_capon_spectrum(pa_mvdr_beamformer_t *bf,
                             const pa_array_config_t *config,
                             const pa_element_t *elements,
                             double *spectrum,
                             const double *theta_vals,
                             const double *phi_vals,
                             uint32_t num_angles)
{
    if (!bf || !config || !elements || !spectrum
        || !theta_vals || !phi_vals || num_angles == 0)
        return;

    uint32_t M = bf->base.num_elements;
    double complex *temp = (double complex *)calloc(M, sizeof(double complex));
    if (!temp) return;

    for (uint32_t k = 0; k < num_angles; k++) {
        /* Compute steering vector for this angle */
        double complex *sv = pa_steering_vector(config, elements,
                                                 theta_vals[k], phi_vals[k]);
        if (!sv) {
            spectrum[k] = 0.0;
            continue;
        }

        /* temp = R^{-1} a */
        mat_vec_mul(bf->base.cov_matrix_inv, sv, temp, M);

        /* P_Capon = 1 / (a^H R^{-1} a) = 1 / (a^H · temp) */
        double complex d = inner_product(sv, temp, M);
        double denom = creal(d);  /* Should be real for Hermitian R */

        if (denom > 1e-15) {
            spectrum[k] = 1.0 / denom;
        } else {
            spectrum[k] = 0.0;
        }

        free(sv);
    }

    free(temp);
}

/**
 * Compute MVDR output SINR.
 *
 * For ideal case (perfect covariance knowledge):
 *   SINR_opt = sigma_s^2 · a^H R^{-1} a
 *
 * In practice, the signal is included in the covariance estimate,
 * so the SINR is approximated from the weight vector:
 *   SINR ≈ |w^H a|^2 / (w^H R w - |w^H a|^2 · sigma_s^2)
 *
 * This returns the worst-case bound based on the current weight vector.
 *
 * @param bf MVDR beamformer.
 * @return SINR in dB.
 */
double pa_mvdr_output_sinr(const pa_mvdr_beamformer_t *bf)
{
    if (!bf || bf->base.num_elements == 0) return -300.0;

    uint32_t M = bf->base.num_elements;

    /* w^H a (should be ≈ 1 for MVDR) */
    double complex wHa = inner_product(bf->base.weight_vector,
                                        bf->base.steering_vector, M);

    /* w^H R w */
    double complex *temp = (double complex *)calloc(M, sizeof(double complex));
    if (!temp) return -300.0;

    mat_vec_mul(bf->base.cov_matrix, bf->base.weight_vector, temp, M);
    double complex wHRw = inner_product(bf->base.weight_vector, temp, M);

    free(temp);

    double signal_power = cabs(wHa) * cabs(wHa);
    double total_power  = creal(wHRw);

    if (total_power < signal_power + 1e-15) {
        return 100.0;  /* Nearly perfect, no interference */
    }

    double interference_plus_noise = total_power - signal_power;
    if (interference_plus_noise < 1e-15) {
        return 100.0;
    }

    double sinr_linear = signal_power / interference_plus_noise;
    return 10.0 * log10(sinr_linear);
}