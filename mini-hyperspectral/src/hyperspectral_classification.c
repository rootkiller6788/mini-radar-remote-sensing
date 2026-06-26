/**
 * @file    hyperspectral_classification.c
 * @brief   Classification algorithms: SAM, Matched Filter, ACE, RX, K-means.
 *
 * @details Implements spectral classification (SAM, ML), target detection
 *          (MF, ACE), anomaly detection (RX), and unsupervised clustering
 *          (k-means) for hyperspectral imagery. Includes confusion matrix
 *          computation for accuracy assessment.
 *
 * Knowledge covered:
 *   L1: Confusion matrix, overall accuracy, kappa coefficient
 *   L3: Mahalanobis distance, Bayes decision rule, discriminant functions
 *   L5: SAM, MF, ACE, RX, k-means algorithms with complete implementations
 *   L6: Subpixel target detection, anomaly detection in cluttered backgrounds
 *
 * Reference:
 *   - Richards & Jia, "Remote Sensing Digital Image Analysis" (2006)
 *   - Manolakis & Shaw, "Detection Algorithms for HS Imaging" (2002)
 *   - Reed & Yu, "Adaptive Multiple-Band CFAR Detection" (1990)
 *   - MacQueen (1967), "Some methods for classification..."
 */

#include "hyperspectral_classification.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "classification: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ─── Label Map ──────────────────────────────────────────────────────────── */

hspec_label_map_t hspec_label_map_alloc(size_t nrows, size_t ncols, size_t nclasses)
{
    hspec_label_map_t lm;
    memset(&lm, 0, sizeof(lm));
    if (nrows == 0 || ncols == 0) return lm;
    lm.nrows = nrows;
    lm.ncols = ncols;
    lm.nclasses = nclasses;
    lm.labels = calloc(nrows * ncols, sizeof(size_t));
    if (nclasses > 0)
        lm.class_names = calloc(nclasses, sizeof(char[64]));
    return lm;
}

void hspec_label_map_free(hspec_label_map_t *lm)
{
    if (!lm) return;
    free(lm->labels);
    free(lm->class_names);
    memset(lm, 0, sizeof(*lm));
}

/* ─── Confusion Matrix ─────────────────────────────────────────────────── */

hspec_confusion_matrix_t hspec_confusion_matrix_alloc(size_t nclasses)
{
    hspec_confusion_matrix_t cm;
    memset(&cm, 0, sizeof(cm));
    if (nclasses == 0) return cm;
    cm.nclasses = nclasses;
    cm.matrix = safe_malloc(nclasses * sizeof(size_t*));
    for (size_t i = 0; i < nclasses; i++) {
        cm.matrix[i] = calloc(nclasses, sizeof(size_t));
    }
    cm.row_sums = calloc(nclasses, sizeof(size_t));
    cm.col_sums = calloc(nclasses, sizeof(size_t));
    cm.producer_accuracy = safe_malloc(nclasses * sizeof(double));
    cm.user_accuracy = safe_malloc(nclasses * sizeof(double));
    return cm;
}

void hspec_confusion_matrix_free(hspec_confusion_matrix_t *cm)
{
    if (!cm) return;
    for (size_t i = 0; i < cm->nclasses; i++) free(cm->matrix[i]);
    free(cm->matrix);
    free(cm->row_sums);
    free(cm->col_sums);
    free(cm->producer_accuracy);
    free(cm->user_accuracy);
    memset(cm, 0, sizeof(*cm));
}

/**
 * Compute confusion matrix from truth and predicted labels.
 *
 * @math M_{ij} = count(pixels: truth = i, predicted = j)
 *       OA = Σ_i M_{ii} / N
 *       κ  = (OA - p_e) / (1 - p_e)
 *       p_e = Σ_i (row_i/N)·(col_i/N)
 *
 * Complexity: O(npixels)
 * Reference: Cohen (1960), "A Coefficient of Agreement for Nominal Scales"
 */
int hspec_confusion_matrix_compute(hspec_confusion_matrix_t *cm,
                                    const size_t *truth, const size_t *predicted,
                                    size_t npixels)
{
    if (!cm || !truth || !predicted || npixels == 0) return -1;

    size_t K = cm->nclasses;
    for (size_t i = 0; i < K; i++) {
        memset(cm->matrix[i], 0, K * sizeof(size_t));
        cm->row_sums[i] = 0;
        cm->col_sums[i] = 0;
    }
    cm->total_samples = npixels;

    for (size_t p = 0; p < npixels; p++) {
        size_t t = truth[p];
        size_t pr = predicted[p];
        if (t < K && pr < K) {
            cm->matrix[t][pr]++;
            cm->row_sums[t]++;
            cm->col_sums[pr]++;
        }
    }

    /* Overall accuracy */
    size_t correct = 0;
    for (size_t i = 0; i < K; i++) correct += cm->matrix[i][i];
    cm->overall_accuracy = (double)correct / (double)npixels;

    /* Producer accuracy (recall): PA_i = M_{ii} / row_i */
    for (size_t i = 0; i < K; i++) {
        cm->producer_accuracy[i] = (cm->row_sums[i] > 0)
            ? (double)cm->matrix[i][i] / (double)cm->row_sums[i] : 0.0;
    }

    /* User accuracy (precision): UA_i = M_{ii} / col_i */
    for (size_t i = 0; i < K; i++) {
        cm->user_accuracy[i] = (cm->col_sums[i] > 0)
            ? (double)cm->matrix[i][i] / (double)cm->col_sums[i] : 0.0;
    }

    /* Cohen's kappa */
    double p_e = 0.0;
    for (size_t i = 0; i < K; i++) {
        p_e += ((double)cm->row_sums[i] / (double)npixels)
             * ((double)cm->col_sums[i] / (double)npixels);
    }
    if (fabs(1.0 - p_e) > 1e-15)
        cm->kappa = (cm->overall_accuracy - p_e) / (1.0 - p_e);
    else
        cm->kappa = 1.0;

    return 0;
}

/* ─── Spectral Angle Mapper (SAM) ─────────────────────────────────────── */

/**
 * Classify a single pixel using Spectral Angle Mapper.
 *
 * @math θ_i = arccos(⟨x, r_i⟩ / (‖x‖·‖r_i‖))
 *
 * The class with the smallest spectral angle is assigned.
 *
 * Complexity: O(nclasses·nbands)
 * Reference: Kruse et al. (1993), Remote Sensing of Environment
 */
hspec_sam_result_t hspec_sam_classify(const hspec_pixel_t *pixel,
                                       const double *ref_spectra,
                                       size_t nclasses, size_t nbands)
{
    hspec_sam_result_t result;
    memset(&result, 0, sizeof(result));
    result.nclasses = nclasses;
    result.angles_rad = safe_malloc(nclasses * sizeof(double));
    result.best_class = 0;
    result.best_angle_rad = M_PI;  /* worst possible */

    if (!pixel || !pixel->reflectance || !ref_spectra || nclasses == 0)
        return result;

    /* Norm of test pixel */
    double norm_x = 0.0;
    for (size_t b = 0; b < nbands; b++) {
        norm_x += pixel->reflectance[b] * pixel->reflectance[b];
    }
    norm_x = sqrt(norm_x);
    if (norm_x < 1e-30) return result;

    for (size_t c = 0; c < nclasses; c++) {
        double dot = 0.0, norm_r = 0.0;
        const double *ref = &ref_spectra[c * nbands];
        for (size_t b = 0; b < nbands; b++) {
            dot += pixel->reflectance[b] * ref[b];
            norm_r += ref[b] * ref[b];
        }
        norm_r = sqrt(norm_r);
        double cos_theta;
        if (norm_r < 1e-30) {
            cos_theta = 0.0;
        } else {
            cos_theta = dot / (norm_x * norm_r);
        }
        if (cos_theta > 1.0) cos_theta = 1.0;
        if (cos_theta < -1.0) cos_theta = -1.0;
        double theta = acos(cos_theta);
        result.angles_rad[c] = theta;
        if (theta < result.best_angle_rad) {
            result.best_angle_rad = theta;
            result.best_class = c;
        }
    }
    result.best_angle_deg = result.best_angle_rad * 180.0 / M_PI;

    /* Confidence: normalized where angle=0 is max confidence */
    if (result.best_angle_rad < M_PI / 2.0) {
        result.confidence = 1.0 - result.best_angle_rad / (M_PI / 2.0);
    } else {
        result.confidence = 0.0;
    }

    /* Clean up if no good match found */
    if (result.best_angle_rad >= M_PI / 2.0) result.confidence = 0.0;

    return result;
}

/**
 * Full-scene SAM classification.
 *
 * Complexity: O(n_pixels·nclasses·nbands)
 */
hspec_label_map_t hspec_sam_classify_scene(const hspec_datacube_t *dc,
                                            const double *ref_spectra,
                                            size_t nclasses)
{
    hspec_label_map_t lm = hspec_label_map_alloc(dc->nrows, dc->ncols, nclasses);
    if (!lm.labels || nclasses == 0) return lm;

    hspec_pixel_t pixel;
    pixel.nbands = dc->nbands;
    pixel.reflectance = safe_malloc(dc->nbands * sizeof(double));

    size_t B = dc->nbands;
    for (size_t p = 0; p < dc->npixels; p++) {
        /* Extract pixel inline for speed */
        for (size_t b = 0; b < B; b++) {
            pixel.reflectance[b] = dc->data[p * B + b];
        }
        hspec_sam_result_t r = hspec_sam_classify(&pixel, ref_spectra,
                                                   nclasses, B);
        lm.labels[p] = r.best_class;
        free(r.angles_rad);
    }
    free(pixel.reflectance);
    return lm;
}

/* ─── Matched Filter ─────────────────────────────────────────────────── */

/**
 * Matched filter detection statistic.
 *
 * @math MF(x) = (tᵀ·Σ⁻¹·(x-μ)) / √(tᵀ·Σ⁻¹·t)
 *
 * Under H₀ (no target), MF ~ N(0,1) if background is Gaussian.
 * Threshold set via inverse normal CDF at p_fa.
 *
 * Complexity: O(nbands²)
 * Reference: Manolakis & Shaw (2002)
 */
hspec_matched_filter_t hspec_matched_filter_detect(const hspec_pixel_t *pixel,
                                                     const double *target,
                                                     const double *bg_mean,
                                                     const double *bg_inv_cov,
                                                     size_t nbands, double p_fa)
{
    hspec_matched_filter_t result;
    memset(&result, 0, sizeof(result));

    if (!pixel || !pixel->reflectance || !target || !bg_mean || !bg_inv_cov
        || nbands == 0)
        return result;

    /* Compute t^T * Σ^{-1} * t  (denominator normalization) */
    double denom = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < nbands; j++) {
            row_sum += bg_inv_cov[i * nbands + j] * target[j];
        }
        denom += target[i] * row_sum;
    }
    if (denom < 1e-30) return result;
    double denom_sqrt = sqrt(denom);

    /* Compute numerator: t^T * Σ^{-1} * (x - μ) */
    double num = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        double centered = pixel->reflectance[i] - bg_mean[i];
        double row_sum = 0.0;
        for (size_t j = 0; j < nbands; j++) {
            row_sum += bg_inv_cov[i * nbands + j] * centered;
        }
        num += target[i] * row_sum;
    }

    result.mf_score = num / denom_sqrt;
    result.snr_estimate = num * num / denom;  /* (t^T Σ^{-1} x)^2 / (t^T Σ^{-1} t) */

    /* Threshold from normal distribution approximation */
    /* For p_fa = 1e-4, threshold ≈ 3.719 (complement of N(0,1) CDF) */
    double z = 0.0;
    if (p_fa > 0.0 && p_fa < 1.0) {
        /* Abramowitz & Stegun approximation for inverse normal CDF */
        double q = p_fa;
        if (q > 0.5) q = 1.0 - q;
        double t = sqrt(-2.0 * log(q));
        z = t - (2.515517 + 0.802853 * t + 0.010328 * t * t)
            / (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t);
    }
    result.threshold = z;
    result.detection = (result.mf_score >= result.threshold) ? 1 : 0;

    return result;
}

/* ─── ACE (Adaptive Coherence Estimator) ─────────────────────────────── */

/**
 * @math ACE(x) = (tᵀ·Σ⁻¹·(x-μ))² / [(tᵀ·Σ⁻¹·t)·((x-μ)ᵀ·Σ⁻¹·(x-μ))]
 *
 * ACE is invariant to scaling of x — it measures the cosine squared
 * of the angle between t and (x-μ) in the whitened space.
 *
 * Complexity: O(nbands²)
 * Reference: Kraut & Scharf (1999), IEEE Trans. Signal Processing
 */
hspec_ace_result_t hspec_ace_detect(const hspec_pixel_t *pixel,
                                     const double *target,
                                     const double *bg_mean,
                                     const double *bg_inv_cov,
                                     size_t nbands, double p_fa)
{
    hspec_ace_result_t result;
    memset(&result, 0, sizeof(result));

    if (!pixel || !pixel->reflectance || !target || !bg_mean || !bg_inv_cov
        || nbands == 0)
        return result;

    /* Compute t^T Σ^{-1} t */
    double sig_inv_s = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        double rs = 0.0;
        for (size_t j = 0; j < nbands; j++)
            rs += bg_inv_cov[i * nbands + j] * target[j];
        sig_inv_s += target[i] * rs;
    }

    /* Compute (x-μ)^T Σ^{-1} (x-μ) */
    double x_sig_x = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        double centered = pixel->reflectance[i] - bg_mean[i];
        double rs = 0.0;
        for (size_t j = 0; j < nbands; j++)
            rs += bg_inv_cov[i * nbands + j] * centered;
        x_sig_x += centered * rs;
    }

    /* Compute numerator */
    double num = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        double centered = pixel->reflectance[i] - bg_mean[i];
        double rs = 0.0;
        for (size_t j = 0; j < nbands; j++)
            rs += bg_inv_cov[i * nbands + j] * centered;
        num += target[i] * rs;
    }
    num = num * num;

    double denom = sig_inv_s * x_sig_x;
    if (denom > 1e-30)
        result.ace_score = num / denom;
    else
        result.ace_score = 0.0;

    if (result.ace_score > 1.0) result.ace_score = 1.0;
    if (result.ace_score < 0.0) result.ace_score = 0.0;

    /* Threshold: ACE follows a beta distribution under H_0 */
    /* Simplified: use fixed threshold based on p_fa */
    double beta_threshold = 1.0 - pow(p_fa, 1.0 / (double)nbands);
    result.threshold = beta_threshold;
    result.detection = (result.ace_score >= result.threshold) ? 1 : 0;

    return result;
}

/* ─── RX Anomaly Detector ────────────────────────────────────────────── */

/**
 * Global Reed-Xiaoli anomaly detector.
 *
 * @math RX(x_p) = (x_p - μ)ᵀ·Σ⁻¹·(x_p - μ)
 *
 * Under H₀ (no anomalies), RX follows χ²_B (chi-squared with B dof).
 * The threshold is set via chi-squared quantile at (1-p_fa).
 *
 * Complexity: O(n_pixels·nbands²)
 * Reference: Reed & Yu (1990), IEEE Trans. ASSP
 */
hspec_rx_result_t hspec_rx_detect_global(const hspec_datacube_t *dc, double p_fa)
{
    hspec_rx_result_t rx;
    memset(&rx, 0, sizeof(rx));
    if (!dc || !dc->data || dc->npixels == 0) return rx;

    size_t N = dc->npixels, B = dc->nbands;
    rx.nrows = dc->nrows;
    rx.ncols = dc->ncols;
    rx.rx_scores = safe_malloc(N * sizeof(double));
    rx.anomaly_mask = safe_malloc(N * sizeof(size_t));

    /* Compute mean vector */
    double *mu = safe_malloc(B * sizeof(double));
    for (size_t b = 0; b < B; b++) {
        double sum = 0.0;
        for (size_t p = 0; p < N; p++) sum += dc->data[p * B + b];
        mu[b] = sum / (double)N;
    }

    /* Compute covariance and its inverse */
    double *cov = safe_malloc(B * B * sizeof(double));
    for (size_t i = 0; i < B; i++) {
        for (size_t j = i; j < B; j++) {
            double s = 0.0;
            for (size_t p = 0; p < N; p++) {
                s += (dc->data[p * B + i] - mu[i])
                   * (dc->data[p * B + j] - mu[j]);
            }
            cov[i * B + j] = s / (double)(N - 1);
            cov[j * B + i] = cov[i * B + j];
        }
    }

    /* Regularized inverse (add small diagonal) */
    double reg = 1e-6;
    for (size_t i = 0; i < B; i++) cov[i * B + i] += reg;

    /* Compute inverse via Gauss-Jordan elimination */
    double *inv_cov = safe_malloc(B * B * sizeof(double));
    /* Augmented matrix [cov | I] */
    double *aug = safe_malloc(B * (2 * B) * sizeof(double));
    for (size_t i = 0; i < B; i++) {
        for (size_t j = 0; j < B; j++) {
            aug[i * (2 * B) + j] = cov[i * B + j];
            aug[i * (2 * B) + B + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    /* Gauss-Jordan */
    for (size_t k = 0; k < B; k++) {
        double pivot = aug[k * (2 * B) + k];
        if (fabs(pivot) < 1e-30) continue;
        for (size_t j = 0; j < 2 * B; j++) aug[k * (2 * B) + j] /= pivot;
        for (size_t i = 0; i < B; i++) {
            if (i == k) continue;
            double factor = aug[i * (2 * B) + k];
            for (size_t j = 0; j < 2 * B; j++)
                aug[i * (2 * B) + j] -= factor * aug[k * (2 * B) + j];
        }
    }
    for (size_t i = 0; i < B; i++)
        for (size_t j = 0; j < B; j++)
            inv_cov[i * B + j] = aug[i * (2 * B) + B + j];
    free(aug);

    /* Compute RX scores */
    double sum_score = 0.0, sum2_score = 0.0;
    for (size_t p = 0; p < N; p++) {
        double score = 0.0;
        for (size_t i = 0; i < B; i++) {
            double centered = dc->data[p * B + i] - mu[i];
            double rs = 0.0;
            for (size_t j = 0; j < B; j++)
                rs += inv_cov[i * B + j] * centered;
            score += centered * rs;
        }
        rx.rx_scores[p] = score;
        sum_score += score;
        sum2_score += score * score;
    }

    rx.mean_score = sum_score / (double)N;
    rx.std_score = sqrt(sum2_score / (double)N - rx.mean_score * rx.mean_score);

    /* Chi-squared threshold */
    rx.threshold = hspec_chi2_threshold(B, p_fa);

    /* Mark anomalies */
    for (size_t p = 0; p < N; p++) {
        if (rx.rx_scores[p] >= rx.threshold) {
            rx.anomaly_mask[p] = 1;
            rx.nanomalies++;
        }
    }

    free(mu); free(cov); free(inv_cov);
    return rx;
}

void hspec_rx_result_free(hspec_rx_result_t *rx)
{
    if (!rx) return;
    free(rx->rx_scores);
    free(rx->anomaly_mask);
    memset(rx, 0, sizeof(*rx));
}

/* ─── K-Means Clustering ────────────────────────────────────────────── */

/**
 * K-means clustering using Lloyd's algorithm.
 *
 * Complexity: O(n_pixels·k·nbands·n_iters)
 * Reference: Lloyd (1982), "Least squares quantization in PCM";
 *            MacQueen (1967)
 */
hspec_kmeans_result_t hspec_kmeans_cluster(const hspec_datacube_t *dc, size_t k,
                                            size_t max_iters, double tol)
{
    hspec_kmeans_result_t km;
    memset(&km, 0, sizeof(km));
    if (!dc || !dc->data || dc->npixels == 0 || k == 0) return km;

    size_t N = dc->npixels, B = dc->nbands;
    km.k = k;
    km.nbands = B;
    km.centroids = safe_malloc(k * B * sizeof(double));
    km.assignments = safe_malloc(N * sizeof(size_t));
    km.cluster_sizes = safe_malloc(k * sizeof(size_t));

    /* Initialize centroids: random pixels */
    for (size_t c = 0; c < k; c++) {
        size_t idx = (c * (N / k)) % N;
        for (size_t b = 0; b < B; b++)
            km.centroids[c * B + b] = dc->data[idx * B + b];
    }

    double prev_twss = 1e308;
    for (size_t iter = 0; iter < max_iters; iter++) {
        /* Assignment step */
        for (size_t c = 0; c < k; c++) km.cluster_sizes[c] = 0;
        double twss = 0.0;

        for (size_t p = 0; p < N; p++) {
            size_t best_c = 0;
            double best_dist = 1e308;
            for (size_t c = 0; c < k; c++) {
                double dist = 0.0;
                for (size_t b = 0; b < B; b++) {
                    double diff = dc->data[p * B + b]
                                - km.centroids[c * B + b];
                    dist += diff * diff;
                }
                if (dist < best_dist) { best_dist = dist; best_c = c; }
            }
            km.assignments[p] = best_c;
            km.cluster_sizes[best_c]++;
            twss += best_dist;
        }

        /* Update step */
        double *new_centroids = safe_malloc(k * B * sizeof(double));
        memset(new_centroids, 0, k * B * sizeof(double));

        for (size_t p = 0; p < N; p++) {
            size_t c = km.assignments[p];
            for (size_t b = 0; b < B; b++)
                new_centroids[c * B + b] += dc->data[p * B + b];
        }
        for (size_t c = 0; c < k; c++) {
            if (km.cluster_sizes[c] > 0) {
                for (size_t b = 0; b < B; b++)
                    new_centroids[c * B + b] /= (double)km.cluster_sizes[c];
            }
        }

        /* Compute centroid shift */
        double max_shift = 0.0;
        for (size_t c = 0; c < k; c++) {
            double shift = 0.0;
            for (size_t b = 0; b < B; b++) {
                double diff = new_centroids[c * B + b]
                            - km.centroids[c * B + b];
                shift += diff * diff;
            }
            shift = sqrt(shift);
            if (shift > max_shift) max_shift = shift;
        }

        memcpy(km.centroids, new_centroids, k * B * sizeof(double));
        free(new_centroids);

        km.niterations = iter + 1;
        km.total_within_ss = twss;

        if (fabs(prev_twss - twss) < tol * (fabs(twss) + 1e-10)) {
            km.converged = 1;
            break;
        }
        prev_twss = twss;
    }

    return km;
}

void hspec_kmeans_result_free(hspec_kmeans_result_t *km)
{
    if (!km) return;
    free(km->centroids);
    free(km->assignments);
    free(km->cluster_sizes);
    memset(km, 0, sizeof(*km));
}

/* ─── Chi-squared Threshold ──────────────────────────────────────────── */

/**
 * Approximate chi-squared quantile using Wilson-Hilferty transform.
 *
 * @math χ²_B(1-p) ≈ B·(1 - 2/(9B) + z·√(2/(9B)))³
 *       where z = Φ⁻¹(1-p) is the normal quantile.
 *
 * Reference: Wilson & Hilferty (1931), PNAS
 */
double hspec_chi2_threshold(size_t nbands, double p_fa)
{
    if (nbands == 0) return 0.0;
    if (p_fa <= 0.0) return 1e10;
    if (p_fa >= 1.0) return 0.0;

    double B = (double)nbands;
    /* Normal quantile via Abramowitz-Stegun */
    double q = p_fa;
    if (q > 0.5) q = 1.0 - q;
    double t_val = sqrt(-2.0 * log(q));
    double z = t_val - (2.515517 + 0.802853 * t_val + 0.010328 * t_val * t_val)
        / (1.0 + 1.432788 * t_val + 0.189269 * t_val * t_val
           + 0.001308 * t_val * t_val * t_val);

    /* Wilson-Hilferty approximation */
    double term = 1.0 - 2.0 / (9.0 * B) + z * sqrt(2.0 / (9.0 * B));
    double chi2 = B * term * term * term;
    if (chi2 < 0.0) chi2 = 0.0;
    return chi2;
}

/* ─── Mahalanobis Distance ──────────────────────────────────────────── */

/**
 * Squared Mahalanobis distance.
 *
 * @math d² = (x - μ)ᵀ·Σ⁻¹·(x - μ)
 *
 * This is a multivariate z-score measuring how many standard deviations
 * x is from μ accounting for correlations.
 *
 * Complexity: O(nbands²)
 */
double hspec_mahalanobis_distance(const double *x, const double *mu,
                                   const double *inv_cov, size_t nbands)
{
    if (!x || !mu || !inv_cov || nbands == 0) return 0.0;

    double dist = 0.0;
    for (size_t i = 0; i < nbands; i++) {
        double cent = x[i] - mu[i];
        double rs = 0.0;
        for (size_t j = 0; j < nbands; j++)
            rs += inv_cov[i * nbands + j] * cent;
        dist += cent * rs;
    }
    return dist;
}