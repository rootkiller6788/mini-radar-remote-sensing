
#ifndef SAR_ADVANCED_H
#define SAR_ADVANCED_H
#include "sar_core.h"

/* --- L8: Compressive Sensing SAR --- */
typedef struct {
    double undersampling_ratio;
    double regularization_lambda;
    int max_iterations;
    double convergence_tol;
    int algorithm_type;
} sar_cs_params_t;

int sar_cs_reconstruct(const double **measurements, size_t n_meas, size_t n_pixels, const sar_cs_params_t *csp, double *reconstructed);

/* --- L8: MIMO-SAR --- */
typedef struct {
    size_t n_tx;
    size_t n_rx;
    double *tx_positions_x;
    double *rx_positions_x;
    double *tx_waveforms;
    size_t waveform_length;
} sar_mimo_config_t;

int sar_mimo_image(const sar_mimo_config_t *cfg, double **image, size_t nrows, size_t ncols);

/* --- L8: Bistatic SAR --- */
typedef struct {
    double tx_x, tx_y, tx_z;
    double rx_x, rx_y, rx_z;
    double tx_vel_x, tx_vel_y, tx_vel_z;
    double rx_vel_x, rx_vel_y, rx_vel_z;
} sar_bistatic_config_t;

double sar_bistatic_range(const sar_bistatic_config_t *bc, double eta, double x_tgt, double y_tgt, double z_tgt);
int sar_bistatic_processor(const sar_raw_data_t *raw, const sar_bistatic_config_t *bc, sar_image_t *image);

/* --- L8: Polarimetric SAR --- */
typedef struct {
    double **HH_I, **HH_Q;
    double **HV_I, **HV_Q;
    double **VH_I, **VH_Q;
    double **VV_I, **VV_Q;
    size_t nrows, ncols;
} sar_polarimetric_t;

typedef struct {
    double T11, T22, T33;
    double T12_real, T12_imag;
    double T13_real, T13_imag;
    double T23_real, T23_imag;
} sar_coherency_matrix_t;

/* Freeman-Durden 3-component decomposition */
void sar_freeman_durden(const sar_coherency_matrix_t *T, size_t n, double *Ps, double *Pd, double *Pv);

/* H-Alpha decomposition */
void sar_h_alpha_decomp(const sar_coherency_matrix_t *T, size_t n, double *entropy, double *alpha, double *anisotropy);

#endif
