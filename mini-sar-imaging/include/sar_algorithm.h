/**
 * @file    sar_algorithm.h
 * @brief   SAR Focusing Algorithms -- L5 + L6
 *
 * @details RDA, CSA, omega-k, Backprojection, SPECAN, PGA/MapDrift autofocus,
 *          and point target impulse response analysis.
 *
 * L5: RDA (range-Doppler), CSA (chirp scaling), omega-k (Stolt),
 *     Backprojection, SPECAN, PGA, MapDrift
 * L6: Point target simulation/focusing, stripmap/spotlight processing, impulse response
 *
 * Reference: Cumming & Wong (2005), Raney et al. (1994), Wahl et al. (1994)
 */
#ifndef SAR_ALGORITHM_H
#define SAR_ALGORITHM_H
#include "sar_core.h"
#include "sar_geometry.h"

/* RDA parameters */
typedef struct {
    int use_secondary_range_compression;
    int rcmc_interpolation_method;
    double range_oversample_factor;
    size_t rcmc_kernel_length;
} sar_rda_params_t;

int sar_rda_process(const sar_raw_data_t *raw, const sar_rda_params_t *rda_p, sar_image_t *image);
void sar_rda_range_compression(const sar_raw_data_t *raw, double **compressed_I, double **compressed_Q);
void sar_rda_rcmc(double **rda_domain_I, double **rda_domain_Q, size_t naz, size_t nrng, const sar_params_t *sp, const sar_rda_params_t *rda_p);
void sar_rda_azimuth_compression(double **rda_domain_I, double **rda_domain_Q, size_t naz, size_t nrng, const sar_params_t *sp);

/* CSA parameters */
typedef struct { double scaling_factor; int use_extended_csa; } sar_csa_params_t;
int sar_csa_process(const sar_raw_data_t *raw, const sar_csa_params_t *csa_p, sar_image_t *image);
void sar_csa_chirp_scaling(double **data_I, double **data_Q, size_t naz, size_t nrng, const sar_params_t *sp, const sar_csa_params_t *csa_p);

/* omega-k */
int sar_omegak_process(const sar_raw_data_t *raw, double R_ref, sar_image_t *image);
void sar_omegak_ref_multiply(double **data_I, double **data_Q, size_t naz, size_t nrng, const sar_params_t *sp, double R_ref);
void sar_omegak_stolt(double **data_I, double **data_Q, size_t naz, size_t nrng, const sar_params_t *sp);

/* BP parameters */
typedef struct {
    double grid_x_min, grid_x_max, grid_y_min, grid_y_max;
    size_t grid_nx, grid_ny;
    double pixel_size_x, pixel_size_y;
} sar_bp_params_t;
int sar_bp_process(const sar_raw_data_t *raw, const sar_bp_params_t *bp_p, const double *traj_x, const double *traj_y, const double *traj_z, sar_image_t *image);

/* SPECAN */
int sar_specan_process(const sar_raw_data_t *raw, size_t burst_len, sar_image_t *image);

/* Autofocus */
int sar_pga_autofocus(sar_image_t *image, int niter);
int sar_mapdrift_autofocus(sar_image_t *image);

/* Impulse response metrics */
typedef struct {
    double range_resolution_3dB_m, azimuth_resolution_3dB_m;
    double range_pslr_db, azimuth_pslr_db;
    double range_islr_db, azimuth_islr_db;
    double peak_range_idx, peak_azimuth_idx;
    double peak_magnitude, peak_phase_rad;
} sar_impulse_response_t;
void sar_analyze_impulse_response(const sar_image_t *image, size_t r0, size_t c0, size_t win, sar_impulse_response_t *ir);

#endif
