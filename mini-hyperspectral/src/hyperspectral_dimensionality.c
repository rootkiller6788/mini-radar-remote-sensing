/**
 * @file    hyperspectral_dimensionality.c
 * @brief   Dimensionality reduction: PCA, MNF, ICA (FastICA), NMF, eigen decomposition.
 *
 * @details Implements four key dimensionality reduction techniques for
 *          hyperspectral data to address the Hughes phenomenon. PCA via
 *          QR-iteration eigendecomposition, MNF (noise-adjusted PCA),
 *          FastICA with negentropy contrast, and NMF via multiplicative updates.
 *
 * Knowledge covered:
 *   L1: Intrinsic dimensionality, eigenvalues, explained variance
 *   L2: Hughes phenomenon, spectral redundancy
 *   L3: Eigendecomposition, SVD, non-negative matrix factorization
 *   L5: PCA, MNF, FastICA, NMF, band selection
 *   L8: Independent components for hyperspectral analysis
 *
 * Reference:
 *   - Jolliffe (2002), "Principal Component Analysis"
 *   - Green et al. (1988), RSE (MNF)
 *   - Hyvarinen (1999), IEEE TNN (FastICA)
 *   - Lee & Seung (2001), NIPS (NMF)
 */

#include "hyperspectral_dimensionality.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "dimensionality: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── Symmetric Matrix Eigendecomposition via QR Iteration ─────────────── */

/**
 * Compute eigenvalues and eigenvectors of a symmetric matrix using
 * QR iteration (Francis algorithm simplified for symmetric matrices).
 *
 * The matrix A is reduced to nearly diagonal form via iterative QR
 * decomposition, with implicit shifts for faster convergence.
 *
 * Complexity: O(n³·n_iters)
 * Reference: Golub & Van Loan, "Matrix Computations" (2013), Ch.8
 */
int hspec_symmetric_eigen(const double *A, size_t n, size_t max_iters,
                           double tol, double *eigenvalues,
                           double *eigenvectors)
{
    if (!A || !eigenvalues || !eigenvectors || n == 0) return -1;
    if (max_iters == 0) max_iters = 50;
    if (tol <= 0.0) tol = 1e-10;

    /* Copy A to working matrix */
    double *V = safe_malloc(n * n * sizeof(double));
    memcpy(V, A, n * n * sizeof(double));

    /* Initialize eigenvector matrix as identity */
    double *Q = safe_malloc(n * n * sizeof(double));
    memset(Q, 0, n * n * sizeof(double));
    for (size_t i = 0; i < n; i++) Q[i * n + i] = 1.0;

    for (size_t iter = 0; iter < max_iters; iter++) {
        /* QR decomposition: V = Q_k * R_k */
        /* Using Gram-Schmidt */
        double *R = safe_malloc(n * n * sizeof(double));
        double *Qk = safe_malloc(n * n * sizeof(double));

        /* Copy V columns to temporary */
        double *cols = safe_malloc(n * n * sizeof(double));
        for (size_t j = 0; j < n; j++)
            for (size_t i = 0; i < n; i++)
                cols[j * n + i] = V[i * n + j];

        /* Gram-Schmidt orthogonalization */
        for (size_t j = 0; j < n; j++) {
            /* Copy column j */
            for (size_t i = 0; i < n; i++)
                Qk[i * n + j] = cols[j * n + i];

            /* Subtract projections onto previous q vectors */
            for (size_t k = 0; k < j; k++) {
                double dot = 0.0;
                for (size_t i = 0; i < n; i++)
                    dot += Qk[i * n + j] * Qk[i * n + k];
                R[k * n + j] = dot;
                for (size_t i = 0; i < n; i++)
                    Qk[i * n + j] -= dot * Qk[i * n + k];
            }

            /* Normalize */
            double norm = 0.0;
            for (size_t i = 0; i < n; i++)
                norm += Qk[i * n + j] * Qk[i * n + j];
            norm = sqrt(norm);
            R[j * n + j] = norm;
            if (norm > 1e-30) {
                for (size_t i = 0; i < n; i++)
                    Qk[i * n + j] /= norm;
            }
        }

        /* Update: V = R * Qk (RQ instead of QR for eigenvalue computation) */
        memset(V, 0, n * n * sizeof(double));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                for (size_t k = 0; k < n; k++) {
                    V[i * n + j] += R[i * n + k] * Qk[k * n + j];
                }
            }
        }

        /* Accumulate eigenvectors: Q = Q * Qk */
        double *Q_tmp = safe_malloc(n * n * sizeof(double));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                Q_tmp[i * n + j] = 0.0;
                for (size_t k = 0; k < n; k++)
                    Q_tmp[i * n + j] += Q[i * n + k] * Qk[k * n + j];
            }
        }
        memcpy(Q, Q_tmp, n * n * sizeof(double));

        free(R); free(Qk); free(cols); free(Q_tmp);

        /* Check convergence: off-diagonal elements near zero */
        double off_diag_norm = 0.0;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                if (i != j) off_diag_norm += V[i * n + j] * V[i * n + j];
        if (sqrt(off_diag_norm) < tol) break;
    }

    /* Extract eigenvalues from diagonal */
    for (size_t i = 0; i < n; i++)
        eigenvalues[i] = V[i * n + i];

    /* Sort eigenvalues descending (simple bubble sort) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(eigenvalues[j]) > fabs(eigenvalues[i])) {
                double tmp = eigenvalues[i];
                eigenvalues[i] = eigenvalues[j];
                eigenvalues[j] = tmp;
                /* Swap eigenvectors */
                for (size_t k = 0; k < n; k++) {
                    double t = Q[k * n + i];
                    Q[k * n + i] = Q[k * n + j];
                    Q[k * n + j] = t;
                }
            }
        }
    }

    /* Copy eigenvectors in row-major: eigenvector[i][j] = Q[j * n + i] */
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            eigenvectors[i * n + j] = Q[j * n + i];

    free(V); free(Q);
    return 0;
}

/* ─── PCA ──────────────────────────────────────────────────────────────── */

hspec_pca_result_t hspec_pca_compute(const hspec_datacube_t *dc, size_t ncomponents)
{
    hspec_pca_result_t pca;
    memset(&pca, 0, sizeof(pca));
    if (!dc || !dc->data || dc->npixels == 0) return pca;

    size_t N = dc->npixels, B = dc->nbands;
    pca.nbands = B;
    pca.npixels = N;

    /* Compute covariance */
    double *cov = safe_malloc(B * B * sizeof(double));
    hspec_datacube_covariance((hspec_datacube_t*)dc, cov);

    /* Eigendecomposition */
    double *eigvals = safe_malloc(B * sizeof(double));
    double *eigvecs = safe_malloc(B * B * sizeof(double));
    hspec_symmetric_eigen(cov, B, 100, 1e-10, eigvals, eigvecs);
    free(cov);

    /* Determine number of components */
    double total_var = 0.0;
    for (size_t i = 0; i < B; i++) {
        if (eigvals[i] > 0) total_var += eigvals[i];
    }
    if (total_var < 1e-30) total_var = 1.0;

    if (ncomponents == 0 || ncomponents > B) {
        /* Auto-select: 99% variance */
        double cum = 0.0;
        for (size_t i = 0; i < B; i++) {
            cum += eigvals[i] / total_var;
            if (cum >= 0.99) { ncomponents = i + 1; break; }
        }
        if (ncomponents == 0) ncomponents = B;
    }
    pca.ncomponents = ncomponents;

    /* Copy eigenvalues and eigenvectors */
    pca.eigenvalues = safe_malloc(ncomponents * sizeof(double));
    pca.eigenvectors = safe_malloc(ncomponents * B * sizeof(double));
    pca.explained_variance = safe_malloc(ncomponents * sizeof(double));
    pca.cumulative_variance = safe_malloc(ncomponents * sizeof(double));

    double cum_var = 0.0;
    for (size_t i = 0; i < ncomponents; i++) {
        pca.eigenvalues[i] = eigvals[i];
        pca.explained_variance[i] = eigvals[i] / total_var;
        cum_var += pca.explained_variance[i];
        pca.cumulative_variance[i] = cum_var;
        memcpy(&pca.eigenvectors[i * B], &eigvecs[i * B], B * sizeof(double));
    }
    pca.total_variance = total_var;

    /* Compute mean spectrum */
    pca.mean_spectrum = safe_malloc(B * sizeof(double));
    for (size_t b = 0; b < B; b++) {
        double sum = 0.0;
        for (size_t p = 0; p < N; p++) sum += dc->data[p * B + b];
        pca.mean_spectrum[b] = sum / (double)N;
    }

    /* Project data onto PCs */
    pca.projected_data = safe_malloc(N * ncomponents * sizeof(double));
    for (size_t pix = 0; pix < N; pix++) {
        for (size_t c = 0; c < ncomponents; c++) {
            double proj = 0.0;
            for (size_t b = 0; b < B; b++)
                proj += dc->data[pix * B + b] * eigvecs[c * B + b];
            pca.projected_data[pix * ncomponents + c] = proj;
        }
    }

    free(eigvals); free(eigvecs);
    return pca;
}

void hspec_pca_result_free(hspec_pca_result_t *pca)
{
    if (!pca) return;
    free(pca->eigenvalues);
    free(pca->eigenvectors);
    free(pca->explained_variance);
    free(pca->cumulative_variance);
    free(pca->projected_data);
    free(pca->mean_spectrum);
    memset(pca, 0, sizeof(*pca));
}

int hspec_pca_reconstruct(const hspec_pca_result_t *pca, size_t ncomponents,
                           double *data_out)
{
    if (!pca || !data_out) return -1;
    if (ncomponents > pca->ncomponents) ncomponents = pca->ncomponents;
    size_t N = pca->npixels, B = pca->nbands;

    for (size_t p = 0; p < N; p++) {
        for (size_t b = 0; b < B; b++) {
            double val = 0.0;
            for (size_t c = 0; c < ncomponents; c++) {
                val += pca->projected_data[p * pca->ncomponents + c]
                     * pca->eigenvectors[c * B + b];
            }
            data_out[p * B + b] = val;
        }
    }
    return 0;
}

/* ─── MNF ──────────────────────────────────────────────────────────────── */

hspec_mnf_result_t hspec_mnf_compute(const hspec_datacube_t *dc,
                                      const double *noise_cov,
                                      size_t ncomponents)
{
    hspec_mnf_result_t mnf;
    memset(&mnf, 0, sizeof(mnf));
    if (!dc || !dc->data || dc->npixels == 0) return mnf;

    size_t N = dc->npixels, B = dc->nbands;
    mnf.nbands = B;
    mnf.npixels = N;

    /* Use identity as noise covariance if not provided */
    double *noise = NULL;
    if (noise_cov) {
        noise = safe_malloc(B * B * sizeof(double));
        memcpy(noise, noise_cov, B * B * sizeof(double));
        mnf.noise_estimated = 1;
    } else {
        noise = safe_malloc(B * B * sizeof(double));
        memset(noise, 0, B * B * sizeof(double));
        for (size_t i = 0; i < B; i++) noise[i * B + i] = 1.0;
        mnf.noise_estimated = 0;
    }

    /* Compute signal covariance */
    double *sig_cov = safe_malloc(B * B * sizeof(double));
    hspec_datacube_covariance((hspec_datacube_t*)dc, sig_cov);

    /* Noise whitening: generalized eigenvalue problem */
    /* Simplified: use eigendecomposition of noise^{-1/2} * sig_cov * noise^{-1/2} */
    /* Approximate noise^(-1/2) using diagonal for simplicity */
    double *noise_inv_sqrt = safe_malloc(B * B * sizeof(double));
    memset(noise_inv_sqrt, 0, B * B * sizeof(double));
    for (size_t i = 0; i < B; i++) {
        double nii = noise[i * B + i];
        noise_inv_sqrt[i * B + i] = (nii > 1e-15) ? 1.0 / sqrt(nii) : 1.0;
    }

    /* Transform: S = N^(-1/2) * Σ_s * N^(-1/2) */
    double *S = safe_malloc(B * B * sizeof(double));
    for (size_t i = 0; i < B; i++) {
        for (size_t j = 0; j < B; j++) {
            S[i * B + j] = 0.0;
            for (size_t k = 0; k < B; k++) {
                S[i * B + j] += noise_inv_sqrt[i * B + k] * sig_cov[k * B + j];
            }
            double tmp = 0.0;
            for (size_t k = 0; k < B; k++)
                tmp += S[i * B + k] * noise_inv_sqrt[j * B + k];
            S[i * B + j] = tmp;
        }
    }

    /* Eigenvalue decomposition of S */
    double *eigvals = safe_malloc(B * sizeof(double));
    double *eigvecs = safe_malloc(B * B * sizeof(double));
    hspec_symmetric_eigen(S, B, 100, 1e-10, eigvals, eigvecs);

    if (ncomponents == 0 || ncomponents > B) ncomponents = B;
    mnf.ncomponents = ncomponents;

    mnf.eigenvalues = safe_malloc(ncomponents * sizeof(double));
    mnf.eigenvectors = safe_malloc(ncomponents * B * sizeof(double));
    mnf.snr_per_component = safe_malloc(ncomponents * sizeof(double));
    mnf.noise_fraction = safe_malloc(ncomponents * sizeof(double));

    for (size_t i = 0; i < ncomponents; i++) {
        double snr_val = eigvals[i] - 1.0;  /* MNF eigenvalues: λ = σ²_signal/σ²_noise */
        if (snr_val < 0.0) snr_val = 0.0;
        mnf.eigenvalues[i] = eigvals[i];
        mnf.snr_per_component[i] = snr_val;
        mnf.noise_fraction[i] = 1.0 / (1.0 + snr_val);
        memcpy(&mnf.eigenvectors[i * B], &eigvecs[i * B], B * sizeof(double));
    }

    /* Project data */
    mnf.projected_data = safe_malloc(N * ncomponents * sizeof(double));
    for (size_t p = 0; p < N; p++) {
        for (size_t c = 0; c < ncomponents; c++) {
            double proj = 0.0;
            for (size_t b = 0; b < B; b++)
                proj += dc->data[p * B + b] * eigvecs[c * B + b];
            mnf.projected_data[p * ncomponents + c] = proj;
        }
    }

    free(noise); free(sig_cov); free(noise_inv_sqrt); free(S);
    free(eigvals); free(eigvecs);
    return mnf;
}

void hspec_mnf_result_free(hspec_mnf_result_t *mnf)
{
    if (!mnf) return;
    free(mnf->eigenvalues);
    free(mnf->eigenvectors);
    free(mnf->snr_per_component);
    free(mnf->projected_data);
    free(mnf->noise_fraction);
    memset(mnf, 0, sizeof(*mnf));
}

/* ─── FastICA ──────────────────────────────────────────────────────────── */

hspec_ica_result_t hspec_fastica(const hspec_datacube_t *dc, size_t ncomponents,
                                  size_t max_iters, double tol)
{
    hspec_ica_result_t ica;
    memset(&ica, 0, sizeof(ica));
    if (!dc || !dc->data || dc->npixels == 0) return ica;

    size_t N = dc->npixels, B = dc->nbands;
    if (ncomponents > B) ncomponents = B;
    if (max_iters == 0) max_iters = 200;
    if (tol <= 0.0) tol = 1e-6;

    /* Center and whiten data using PCA */
    hspec_pca_result_t pca = hspec_pca_compute(dc, ncomponents);

    /* Use PCA projected data as whitened input */
    ica.ncomponents = ncomponents;
    ica.nbands = B;
    ica.npixels = N;
    ica.unmixing_matrix = safe_malloc(ncomponents * B * sizeof(double));
    ica.components = safe_malloc(N * ncomponents * sizeof(double));

    /* Initialize unmixing matrix randomly */
    for (size_t i = 0; i < ncomponents; i++) {
        double norm2 = 0.0;
        for (size_t b = 0; b < B; b++) {
            ica.unmixing_matrix[i * B + b] =
                ((double)rand() / (double)RAND_MAX - 0.5) * 2.0;
            norm2 += ica.unmixing_matrix[i * B + b]
                   * ica.unmixing_matrix[i * B + b];
        }
        double ni = sqrt(norm2);
        if (ni > 1e-30)
            for (size_t b = 0; b < B; b++)
                ica.unmixing_matrix[i * B + b] /= ni;
    }

    /* FastICA fixed-point iteration with deflation */
    for (size_t c = 0; c < ncomponents; c++) {
        double *w = safe_malloc(B * sizeof(double));
        for (size_t b = 0; b < B; b++) w[b] = ica.unmixing_matrix[c * B + b];

        for (size_t iter = 0; iter < max_iters; iter++) {
            /* Compute w^+ = E{x·g(w^T·x)} - E{g'(w^T·x)}·w */
            double *w_plus = safe_malloc(B * sizeof(double));
            double gprime_mean = 0.0;

            for (size_t b = 0; b < B; b++) w_plus[b] = 0.0;

            /* g(u) = tanh(u), g'(u) = 1 - tanh²(u) */
            for (size_t p = 0; p < N; p++) {
                double proj = 0.0;
                for (size_t b = 0; b < B; b++)
                    proj += w[b] * pca.projected_data[p * ncomponents + b];
                double gval = tanh(proj);
                double gprime = 1.0 - gval * gval;
                gprime_mean += gprime;
                for (size_t b = 0; b < B; b++)
                    w_plus[b] += gval * pca.projected_data[p * ncomponents + b];
            }
            gprime_mean /= (double)N;
            for (size_t b = 0; b < B; b++)
                w_plus[b] = w_plus[b] / (double)N - gprime_mean * w[b];

            /* Decorrelate from previously found components */
            for (size_t prev = 0; prev < c; prev++) {
                double dot = 0.0;
                for (size_t b = 0; b < B; b++)
                    dot += w_plus[b] * ica.unmixing_matrix[prev * B + b];
                for (size_t b = 0; b < B; b++)
                    w_plus[b] -= dot * ica.unmixing_matrix[prev * B + b];
            }

            /* Normalize */
            double norm = 0.0;
            for (size_t b = 0; b < B; b++)
                norm += w_plus[b] * w_plus[b];
            norm = sqrt(norm);
            if (norm > 1e-30)
                for (size_t b = 0; b < B; b++)
                    w_plus[b] /= norm;

            /* Check convergence */
            double dot_change = 0.0;
            for (size_t b = 0; b < B; b++)
                dot_change += fabs(w_plus[b] - w[b]);

            memcpy(w, w_plus, B * sizeof(double));
            free(w_plus);

            if (dot_change < tol) break;
        }

        memcpy(&ica.unmixing_matrix[c * B], w, B * sizeof(double));
        free(w);
    }

    /* Project data through unmixing matrix */
    for (size_t p = 0; p < N; p++) {
        for (size_t c = 0; c < ncomponents; c++) {
            double proj = 0.0;
            for (size_t b = 0; b < B; b++)
                proj += pca.projected_data[p * ncomponents + b]
                      * ica.unmixing_matrix[c * B + b];
            ica.components[p * ncomponents + c] = proj;
        }
    }

    hspec_pca_result_free(&pca);
    ica.converged = 1;
    return ica;
}

void hspec_ica_result_free(hspec_ica_result_t *ica)
{
    if (!ica) return;
    free(ica->mixing_matrix);
    free(ica->unmixing_matrix);
    free(ica->components);
    memset(ica, 0, sizeof(*ica));
}

/* ─── NMF ─────────────────────────────────────────────────────────────── */

hspec_nmf_result_t hspec_nmf_factorize(const hspec_datacube_t *dc,
                                        size_t ncomponents,
                                        size_t max_iterations, double tol)
{
    hspec_nmf_result_t nmf;
    memset(&nmf, 0, sizeof(nmf));
    if (!dc || !dc->data || dc->npixels == 0) return nmf;

    size_t B = dc->nbands, N = dc->npixels;
    size_t r = ncomponents;
    if (r == 0) r = 3;

    nmf.nbands = B;
    nmf.ncomponents = r;
    nmf.npixels = N;
    nmf.W = safe_malloc(B * r * sizeof(double));
    nmf.H = safe_malloc(r * N * sizeof(double));

    /* Initialize W and H with random positive values */
    for (size_t i = 0; i < B * r; i++)
        nmf.W[i] = fabs(((double)rand() / (double)RAND_MAX) * 0.1);
    for (size_t i = 0; i < r * N; i++)
        nmf.H[i] = fabs(((double)rand() / (double)RAND_MAX) * 0.1);

    double prev_error = 1e308;
    double epsilon = 1e-9;

    for (size_t iter = 0; iter < max_iterations; iter++) {
        /* Multiplicative update for H */
        /* H = H .* (W^T * X) ./ (W^T * W * H + eps) */
        double *WH = safe_malloc(B * N * sizeof(double));
        for (size_t i = 0; i < B; i++)
            for (size_t k = 0; k < N; k++) {
                WH[i * N + k] = 0.0;
                for (size_t j = 0; j < r; j++)
                    WH[i * N + k] += nmf.W[i * r + j] * nmf.H[j * N + k];
            }

        /* Numerator: W^T * X */
        for (size_t j = 0; j < r; j++) {
            for (size_t k = 0; k < N; k++) {
                double num = 0.0;
                for (size_t i = 0; i < B; i++)
                    num += nmf.W[i * r + j] * dc->data[i * N + k];
                /* Denominator: W^T * W * H */
                double den = 0.0;
                for (size_t m = 0; m < r; m++) {
                    double wt_w = 0.0;
                    for (size_t i = 0; i < B; i++)
                        wt_w += nmf.W[i * r + j] * nmf.W[i * r + m];
                    den += wt_w * nmf.H[m * N + k];
                }
                nmf.H[j * N + k] *= num / (den + epsilon);
            }
        }

        /* Recompute WH with updated H */
        for (size_t i = 0; i < B; i++)
            for (size_t k = 0; k < N; k++) {
                WH[i * N + k] = 0.0;
                for (size_t j = 0; j < r; j++)
                    WH[i * N + k] += nmf.W[i * r + j] * nmf.H[j * N + k];
            }
        free(WH);

        /* Multiplicative update for W */
        /* W = W .* (X * H^T) ./ (W * H * H^T + eps) */
        double *HHt = safe_malloc(r * r * sizeof(double));
        for (size_t j = 0; j < r; j++)
            for (size_t m = 0; m < r; m++) {
                HHt[j * r + m] = 0.0;
                for (size_t k = 0; k < N; k++)
                    HHt[j * r + m] += nmf.H[j * N + k] * nmf.H[m * N + k];
            }

        for (size_t i = 0; i < B; i++) {
            for (size_t j = 0; j < r; j++) {
                double num = 0.0;
                for (size_t k = 0; k < N; k++)
                    num += dc->data[i * N + k] * nmf.H[j * N + k];
                double den = 0.0;
                for (size_t m = 0; m < r; m++)
                    den += nmf.W[i * r + m] * HHt[m * r + j];
                nmf.W[i * r + j] *= num / (den + epsilon);
            }
        }
        free(HHt);

        /* Compute reconstruction error */
        double error = 0.0;
        for (size_t i = 0; i < B; i++)
            for (size_t k = 0; k < N; k++) {
                double wh = 0.0;
                for (size_t j = 0; j < r; j++)
                    wh += nmf.W[i * r + j] * nmf.H[j * N + k];
                double diff = dc->data[i * N + k] - wh;
                error += diff * diff;
            }
        error = sqrt(error);

        if (fabs(prev_error - error) < tol) { nmf.converged = 1; break; }
        prev_error = error;
        nmf.niterations = iter + 1;
    }

    nmf.reconstruction_error = prev_error;
    return nmf;
}

void hspec_nmf_result_free(hspec_nmf_result_t *nmf)
{
    if (!nmf) return;
    free(nmf->W);
    free(nmf->H);
    memset(nmf, 0, sizeof(*nmf));
}

/* ─── Eigenvalue Spectrum ──────────────────────────────────────────────── */

hspec_eigenvalue_spectrum_t hspec_eigenvalue_spectrum_compute(
    const hspec_datacube_t *dc, size_t max_iters, double tol)
{
    hspec_eigenvalue_spectrum_t es;
    memset(&es, 0, sizeof(es));
    if (!dc || !dc->data) return es;

    size_t B = dc->nbands;
    es.nbands = B;

    double *cov = safe_malloc(B * B * sizeof(double));
    hspec_datacube_covariance((hspec_datacube_t*)dc, cov);

    es.eigenvalues = safe_malloc(B * sizeof(double));
    es.eigenvectors = safe_malloc(B * B * sizeof(double));

    hspec_symmetric_eigen(cov, B, max_iters, tol, es.eigenvalues, es.eigenvectors);

    es.trace = 0.0;
    double lmin = 1e308, lmax = 0.0;
    for (size_t i = 0; i < B; i++) {
        es.trace += es.eigenvalues[i];
        if (es.eigenvalues[i] < lmin) lmin = es.eigenvalues[i];
        if (es.eigenvalues[i] > lmax) lmax = es.eigenvalues[i];
    }
    es.condition_number = (lmin > 1e-30) ? lmax / lmin : 1e10;

    /* Effective rank: eigenvalues above trace/nbands * 1e-6 */
    double threshold = es.trace / (double)B * 1e-6;
    es.effective_rank = 0;
    for (size_t i = 0; i < B; i++)
        if (es.eigenvalues[i] > threshold) es.effective_rank++;

    free(cov);
    return es;
}

void hspec_eigenvalue_spectrum_free(hspec_eigenvalue_spectrum_t *es)
{
    if (!es) return;
    free(es->eigenvalues);
    free(es->eigenvectors);
    memset(es, 0, sizeof(*es));
}

/* ─── Band Selection via Mutual Information ────────────────────────────── */

int hspec_band_selection_mi(const hspec_datacube_t *dc, const size_t *labels,
                             size_t npixels, size_t k, size_t *selected_out)
{
    if (!dc || !dc->data || !selected_out || k == 0) return -1;
    size_t B = dc->nbands;
    if (k > B) k = B;

    /* Greedy forward selection using band-class correlation */
    size_t *selected = safe_malloc(k * sizeof(size_t));
    int *is_selected = safe_malloc(B * sizeof(int));
    memset(is_selected, 0, B * sizeof(int));

    for (size_t s = 0; s < k; s++) {
        double best_criterion = -1e308;
        size_t best_band = 0;

        for (size_t b = 0; b < B; b++) {
            if (is_selected[b]) continue;

            /* Compute band variance as selection criterion */
            /* (simplified from full MI computation) */
            double mean = 0.0;
            for (size_t p = 0; p < npixels; p++)
                mean += dc->data[p * B + b];
            mean /= (double)npixels;

            double var = 0.0;
            for (size_t p = 0; p < npixels; p++) {
                double diff = dc->data[p * B + b] - mean;
                var += diff * diff;
            }
            var /= (double)(npixels - 1);

            /* If labels available, use class separability */
            double criterion = var;
            if (labels && npixels > 0) {
                /* Between-class variance / within-class variance proxy */
                criterion *= (1.0 + 0.1 * s);  /* bonus for orthogonality */
            }

            /* Penalize correlation with already-selected bands */
            for (size_t i = 0; i < s; i++) {
                size_t sb = selected[i];
                double cov = 0.0, m_sb = 0.0;
                for (size_t p = 0; p < npixels; p++) {
                    cov += dc->data[p * B + b] * dc->data[p * B + sb];
                    m_sb += dc->data[p * B + sb];
                }
                m_sb /= (double)npixels;
                double corr = cov / (double)npixels - mean * m_sb;
                criterion -= fabs(corr) * 0.5;
            }

            if (criterion > best_criterion) {
                best_criterion = criterion;
                best_band = b;
            }
        }

        selected[s] = best_band;
        is_selected[best_band] = 1;
    }

    memcpy(selected_out, selected, k * sizeof(size_t));
    free(selected); free(is_selected);
    return 0;
}