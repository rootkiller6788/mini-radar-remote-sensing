/* ============================================================================
 * phased_array_aesa.c — Active Electronically Scanned Array (AESA)
 *
 * Implements AESA-specific functions:
 *   - T/R module initialization with realistic parameters
 *   - Array-level G/T (gain over noise temperature) figure of merit
 *   - NEDT (Noise-Equivalent Delta Temperature) for radiometry
 *   - Phase shifter effects on beamforming
 *   - LCMV (Linearly Constrained Minimum Variance) beamformer
 *
 * AESA is the modern standard for fighter radar (F-35 AN/APG-81,
 * F-22 AN/APG-77), ground-based air defense (Patriot MPQ-65),
 * naval radar (SPY-6), and satellite communications.
 *
 * Key AESA characteristics:
 *   - Each element has its own T/R module (GaAs/GaN MMIC)
 *   - ~500–2000+ T/R modules per array
 *   - ~2–10W peak power per module (X-band)
 *   - 4–8 bit digital phase shifters
 *   - Digital beamforming (DBF) in newer systems
 *
 * Reference: Brookner, E. (2008) "Phased Arrays Around the World",
 *            Microwave Journal.
 *            Herd, J.S. & Conway, M.D. (2016) "The Evolution to Modern
 *            Phased Array Architectures", Proc. IEEE.
 * ============================================================================ */

#include "phased_array.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* Forward declarations */
static int matrix_invert_gj(double complex *A, double complex *A_inv, uint32_t M);

/* ============================================================================
 * L7: AESA T/R Module Initialization
 *
 * Models a typical X-band T/R module (8–12 GHz).
 *
 * Modern GaN (Gallium Nitride) T/R modules offer:
 *   - Higher power density than GaAs (5–10×)
 *   - Higher efficiency (~40-50% PAE)
 *   - Better thermal conductivity
 *   - Higher breakdown voltage
 *
 * F-35 AN/APG-81 (Northrop Grumman):
 *   - ~1,676 T/R modules
 *   - X-band (8–12 GHz)
 *   - ~10W peak power per module
 *   - ~4-bit phase shifters
 *   - Detect 1 m² RCS at >100 km
 *   - Simultaneous air-to-air and air-to-ground modes
 *
 * Patriot MPQ-65 (Raytheon):
 *   - ~5,000+ elements (C-band)
 *   - GaAs MMIC T/R modules
 *   - Track 100+ targets simultaneously
 *   - Range ~170 km vs. 1 m² RCS
 * ============================================================================ */

pa_tr_module_t pa_tr_module_init(int module_id)
{
    pa_tr_module_t mod;
    memset(&mod, 0, sizeof(mod));

    mod.module_id          = module_id;
    mod.tx_power_watt      = 10.0;     /* 10W peak (X-band GaAs typical) */
    mod.tx_gain_db         = 25.0;     /* MMIC PA gain */
    mod.rx_gain_db         = 20.0;     /* MMIC LNA gain */
    mod.noise_figure_db    = 3.0;      /* LNA noise figure */
    mod.phase_resolution_deg = 22.5;   /* 4-bit: 360/16 = 22.5° */
    mod.phase_bits          = PA_PHASE_BITS_4;
    mod.attenuator_range_db = 31.5;    /* 5-bit attenuator, 0.5 dB step */
    mod.attenuator_step_db  = 0.5;
    mod.switching_time_ns   = 50.0;    /* T/R switch time */
    mod.bandwidth_mhz       = 4000.0;  /* 4 GHz instantaneous BW (X-band) */
    mod.is_healthy          = 1;       /* Built-in test passed */
    mod.temperature_degc    = 65.0;    /* Typical operating temperature */

    return mod;
}

/**
 * Compute array-level G/T ratio.
 *
 * G/T = G_rx(dB) − 10 log₁₀(T_sys)
 *
 * where:
 *   G_rx = receive gain of the array (dB)
 *   T_sys = system noise temperature (K) = T_antenna + T_receiver
 *
 * T_sys accounts for:
 *   T_antenna: sky noise, ground noise (10–100 K typical for clear sky)
 *   T_receiver: (F − 1) · T₀  where T₀ = 290 K
 *
 * G/T is the figure of merit for receive sensitivity. Higher G/T
 * means better SNR for a given signal level. Typical values:
 *   - Small array (16 elements): ~0–5 dB/K
 *   - Medium array (100 elements): ~15–20 dB/K
 *   - Large AESA (1000+ elements): ~30–35 dB/K
 *
 * The G/T determines the communication link margin or radar detection range.
 * For satellite ground stations (e.g., Iridium, Starlink), G/T > 15 dB/K
 * is typical.
 *
 * @param config Array config.
 * @param elements Array elements.
 * @param system_noise_temp_k System noise temperature (K).
 * @return G/T in dB/K.
 */
double pa_array_gt_fom(const pa_array_config_t *config,
                        const pa_element_t *elements,
                        double system_noise_temp_k)
{
    if (!config || !elements || config->num_elements == 0) return 0.0;

    double G_rx = pa_directivity(config, elements, NULL);
    if (G_rx <= 0.0) return 0.0;

    double G_rx_db = 10.0 * log10(G_rx);

    if (system_noise_temp_k < 1.0) system_noise_temp_k = 290.0;

    double G_over_T = G_rx_db - 10.0 * log10(system_noise_temp_k);
    return G_over_T;
}

/**
 * Compute array NEDT (Noise-Equivalent Delta Temperature).
 *
 * NEDT is the temperature difference that produces an SNR of 1
 * at the detector output. It is the fundamental sensitivity metric
 * for passive radiometric imaging (e.g., weather radar, radio astronomy).
 *
 * NEDT = T_sys / sqrt(B · τ)
 *
 * where:
 *   T_sys = system noise temperature (K)
 *   B     = receiver bandwidth (Hz)
 *   τ     = integration time (s)
 *
 * For an AESA used for radiometry, the large number of receivers
 * allows imaging with reduced integration time per pixel.
 *
 * Typical values:
 *   - Weather radar: NEDT ~ 0.5 K
 *   - Radio astronomy: NEDT ~ 0.01 K (with long integration)
 *   - Security screening: NEDT ~ 0.1 K (needed for concealed object detection)
 *
 * @param tr_array Array of T/R modules.
 * @param num_modules Number of T/R modules.
 * @return NEDT (Kelvin).
 */
double pa_array_nedt(const pa_tr_module_t *tr_array, uint32_t num_modules)
{
    if (!tr_array || num_modules == 0) return 0.0;

    /* Average noise figure across modules */
    double avg_nf_db = 0.0;
    for (uint32_t i = 0; i < num_modules; i++) {
        avg_nf_db += tr_array[i].noise_figure_db;
    }
    avg_nf_db /= (double)num_modules;

    /* System noise temperature: T_sys = T₀ · (F − 1) + T_antenna
     * Assume T_antenna ≈ 50 K (clear sky at X-band) */
    double F = pow(10.0, avg_nf_db / 10.0);
    double T0 = 290.0;
    double T_antenna = 50.0;
    double T_sys = T0 * (F - 1.0) + T_antenna;

    /* Typical radiometer bandwidth and integration time */
    double B = 100.0e6;  /* 100 MHz */
    double tau = 0.001;   /* 1 ms integration */

    /* NEDT = T_sys / sqrt(B * tau) */
    double BT = B * tau;
    if (BT <= 0.0) return 0.0;

    return T_sys / sqrt(BT);
}

/* ============================================================================
 * L8: LCMV (Linearly Constrained Minimum Variance) Beamformer
 *
 * The LCMV beamformer generalizes MVDR by allowing multiple linear
 * constraints. This enables:
 *
 *   - Derivative constraints (robustness against pointing errors)
 *   - Multiple look directions (multi-beam)
 *   - Null constraints (forced nulls in known interferer directions)
 *   - Eigenvector constraints (maximize SINR)
 *
 * LCMV optimization:
 *   minimize  w^H R w
 *   subject to  C^H w = f
 *
 * where C is an M×K constraint matrix (K = number of constraints)
 * and f is a K×1 desired response vector.
 *
 * Solution (Van Trees §6.3.2):
 *   w_LCMV = R^{-1} C (C^H R^{-1} C)^{-1} f
 *
 * For K = 1 with C = a (steering vector) and f = 1:
 *   w_LCMV = w_MVDR  (MVDR is a special case of LCMV)
 *
 * For point-null constraints:
 *   C = [a(theta_s), a(theta_null_1), ..., a(theta_null_{K-1})]
 *   f = [1, 0, 0, ..., 0]^T
 *   (unity gain at look direction, nulls at the other directions)
 *
 * Derivative constraints (robust beamforming):
 *   C = [a, ∂a/∂theta, ∂a/∂phi, ...]
 *   f = [1, 0, 0, ...]^T
 *   This forces the beamformer to be flat near the look direction,
 *   which reduces sensitivity to pointing errors.
 *
 * ============================================================================ */

/**
 * Static helper: constrain matrix multiplication C^H R^{-1} C.
 *
 * Computes: M = C^H · (R_inv · C)
 *
 * @param R_inv M×M inverse covariance (row-major).
 * @param C M×K constraint matrix (row-major).
 * @param M Number of array elements.
 * @param K Number of constraints.
 * @param result Output: K×K complex matrix (row-major).
 */
static void compute_quadratic_form(const double complex *R_inv,
                                    const double complex *C,
                                    uint32_t M, uint32_t K,
                                    double complex *result)
{
    /* First: temp = R_inv · C   [M×M] × [M×K] = M×K */
    double complex *temp = (double complex *)calloc((size_t)M * K,
                                                      sizeof(double complex));
    if (!temp) return;

    for (uint32_t k = 0; k < K; k++) {
        for (uint32_t i = 0; i < M; i++) {
            double complex sum = 0.0 + 0.0 * I;
            for (uint32_t j = 0; j < M; j++) {
                sum += R_inv[i * M + j] * C[j * K + k];
            }
            temp[i * K + k] = sum;
        }
    }

    /* Then: result = C^H · temp   [K×M] × [M×K] = K×K
     * result[k1][k2] = Σ_i conj(C[i][k1]) * temp[i][k2] */
    for (uint32_t k1 = 0; k1 < K; k1++) {
        for (uint32_t k2 = 0; k2 < K; k2++) {
            double complex sum = 0.0 + 0.0 * I;
            for (uint32_t i = 0; i < M; i++) {
                sum += conj(C[i * K + k1]) * temp[i * K + k2];
            }
            result[k1 * K + k2] = sum;
        }
    }

    free(temp);
}

/**
 * Compute LCMV weight vector.
 *
 * w_LCMV = R^{-1} C (C^H R^{-1} C)^{-1} f
 *
 * This function updates bf->base.weight_vector with the LCMV solution.
 *
 * @param bf LCMV beamformer (must have cov_matrix_inv computed).
 * @param config Array config.
 * @param elements Array elements.
 */
static void pa_lcmv_compute_weights_internal(pa_lcmv_beamformer_t *bf,
                                              const pa_array_config_t *config,
                                              const pa_element_t *elements)
{
    if (!bf || !config || !elements) return;

    uint32_t M = bf->base.num_elements;
    uint32_t K = bf->num_constraints;

    if (M == 0 || K == 0) return;

    /* Step 1: Compute quadratic form Q = C^H R^{-1} C (K×K) */
    double complex *Q = (double complex *)calloc((size_t)K * K,
                                                   sizeof(double complex));
    if (!Q) return;

    compute_quadratic_form(bf->base.cov_matrix_inv, bf->constraint_matrix,
                           M, K, Q);

    /* Step 2: Invert Q (K×K) */
    double complex *Q_copy = (double complex *)malloc((size_t)K * K
                                                       * sizeof(double complex));
    double complex *Q_inv  = (double complex *)malloc((size_t)K * K
                                                       * sizeof(double complex));
    if (!Q_copy || !Q_inv) {
        free(Q); free(Q_copy); free(Q_inv);
        return;
    }

    memcpy(Q_copy, Q, (size_t)K * K * sizeof(double complex));

    int ret = 0;
    if (K == 1) {
        /* Scalar inverse for K=1 */
        if (cabs(Q_copy[0]) > 1e-15) {
            Q_inv[0] = 1.0 / Q_copy[0];
        } else {
            ret = -1;
        }
    } else {
        ret = (int)matrix_invert_gj(Q_copy, Q_inv, K);
    }

    free(Q_copy);

    if (ret != 0) {
        /* Use identity for failed inverse */
        for (uint32_t i = 0; i < K; i++)
            for (uint32_t j = 0; j < K; j++)
                Q_inv[i * K + j] = (i == j) ? (1.0 + 0.0 * I) : 0.0;
    }

    /* Step 3: Compute Q_inv · f */
    double complex *qf = (double complex *)calloc(K, sizeof(double complex));
    if (!qf) { free(Q); free(Q_inv); return; }

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < K; j++) {
            qf[i] += Q_inv[i * K + j] * bf->constraint_response[j];
        }
    }

    /* Step 4: w = R^{-1} · C · (Q^{-1} f)
     * First: temp = C · (Q^{-1} f) = C · qf   [M×K] × [K×1] = M×1 */
    double complex *temp2 = (double complex *)calloc(M, sizeof(double complex));
    if (!temp2) { free(Q); free(Q_inv); free(qf); return; }

    for (uint32_t i = 0; i < M; i++) {
        for (uint32_t k = 0; k < K; k++) {
            temp2[i] += bf->constraint_matrix[i * K + k] * qf[k];
        }
    }

    /* Then: w = R^{-1} · temp2 */
    for (uint32_t i = 0; i < M; i++) {
        bf->base.weight_vector[i] = 0.0 + 0.0 * I;
        for (uint32_t j = 0; j < M; j++) {
            bf->base.weight_vector[i] +=
                bf->base.cov_matrix_inv[i * M + j] * temp2[j];
        }
    }

    free(Q); free(Q_inv); free(qf); free(temp2);
}

/* Gauss-Jordan matrix inversion wrapper (external linkage):
 * Same algorithm from phased_array_adaptive.c, duplicated here
 * for self-contained LCMV computation.
 * In a production system, this would be in a shared utility. */
static int matrix_invert_gj(double complex *A, double complex *A_inv, uint32_t M)
{
    if (M == 0 || !A || !A_inv) return -1;
    if (M == 1) {
        if (cabs(A[0]) < 1e-15) return -1;
        A_inv[0] = 1.0 / A[0];
        return 0;
    }

    /* Initialize identity */
    for (uint32_t i = 0; i < M; i++)
        for (uint32_t j = 0; j < M; j++)
            A_inv[i * M + j] = (i == j) ? (1.0 + 0.0 * I) : (0.0 + 0.0 * I);

    for (uint32_t col = 0; col < M; col++) {
        uint32_t pivot_row = col;
        double pivot_mag = cabs(A[col * M + col]);
        for (uint32_t row = col + 1; row < M; row++) {
            double mag = cabs(A[row * M + col]);
            if (mag > pivot_mag) { pivot_mag = mag; pivot_row = row; }
        }
        if (pivot_mag < 1e-15) return -1;
        if (pivot_row != col) {
            for (uint32_t j = 0; j < M; j++) {
                double complex t = A[col * M + j];
                A[col * M + j] = A[pivot_row * M + j];
                A[pivot_row * M + j] = t;
                t = A_inv[col * M + j];
                A_inv[col * M + j] = A_inv[pivot_row * M + j];
                A_inv[pivot_row * M + j] = t;
            }
        }
        double complex pivot = A[col * M + col];
        for (uint32_t j = 0; j < M; j++) {
            A[col * M + j] /= pivot;
            A_inv[col * M + j] /= pivot;
        }
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