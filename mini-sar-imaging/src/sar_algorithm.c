/**
 * @file    sar_algorithm.c
 * @brief   SAR Focusing Algorithms -- L5 Algorithms + L6 Canonical Problems
 * @details RDA, CSA, omega-k, Backprojection, SPECAN, Autofocus (PGA, MapDrift)
 * Reference: Cumming & Wong (2005); Raney et al. (1994); Wahl et al. (1994)
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "sar_core.h"
#include "sar_geometry.h"
#include "sar_algorithm.h"
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

/* ==================================================================
 * L5: Range-Doppler Algorithm (RDA)
 * ================================================================== */

/* Forward decl: 1D FFT helper from sar_core.c */
void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward);

/**
 * RDA processing chain (L5):
 * 1. Range compression (frequency-domain matched filter)
 * 2. Azimuth FFT -> range-Doppler domain
 * 3. RCMC (sinc interpolation)
 * 4. Azimuth compression (range-dependent matched filter)
 * 5. Azimuth IFFT -> focused SLC
 *
 * The RDA is the most widely used SAR processing algorithm because
 * it separates range and azimuth processing via the range-Doppler
 * domain, achieving O(N^2 log N) efficiency vs O(N^3) for time-domain.
 *
 * Limitation: Parabolic range equation assumption breaks down for
 * wide swaths / high squint angles (handled by CSA / omega-k).
 */
int sar_rda_process(const sar_raw_data_t *raw,
                    const sar_rda_params_t *rda_p,
                    sar_image_t *image)
{
    if (!raw || !rda_p || !image) return -1;
    size_t naz = raw->naz, nrng = raw->nrng;
    if (image->nrows != naz || image->ncols != nrng) return -2;

    /* Allocate working buffers */
    double **work_I = (double **)malloc(naz * sizeof(double *));
    double **work_Q = (double **)malloc(naz * sizeof(double *));
    if (!work_I || !work_Q) { free(work_I); free(work_Q); return -3; }
    for (size_t i=0; i<naz; i++) {
        work_I[i] = (double *)malloc(nrng * sizeof(double));
        work_Q[i] = (double *)malloc(nrng * sizeof(double));
        if (!work_I[i] || !work_Q[i]) {
            for (size_t j=0; j<i; j++) { free(work_I[j]); free(work_Q[j]); }
            free(work_I); free(work_Q); return -3;
        }
    }

    /* Step 1: Range compression (frequency domain) */
    sar_rda_range_compression(raw, work_I, work_Q);

    /* Step 2: Azimuth FFT -> range-Doppler domain */
    double **rd_I = (double **)malloc(naz * sizeof(double *));
    double **rd_Q = (double **)malloc(naz * sizeof(double *));
    for (size_t i=0; i<naz; i++) {
        rd_I[i] = (double *)malloc(nrng * sizeof(double));
        rd_Q[i] = (double *)malloc(nrng * sizeof(double));
    }
    /* Transpose + azimuth FFT */
    for (size_t c=0; c<nrng; c++) {
        double *col_I = (double *)malloc(naz * sizeof(double));
        double *col_Q = (double *)malloc(naz * sizeof(double));
        for (size_t r=0; r<naz; r++) {
            col_I[r] = work_I[r][c];
            col_Q[r] = work_Q[r][c];
        }
        sar_fft_range_1d(col_I, col_Q, naz, 1);
        for (size_t r=0; r<naz; r++) {
            rd_I[r][c] = col_I[r];
            rd_Q[r][c] = col_Q[r];
        }
        free(col_I); free(col_Q);
    }

    /* Step 3: RCMC */
    sar_rda_rcmc(rd_I, rd_Q, naz, nrng, &raw->params, rda_p);

    /* Step 4: Azimuth compression */
    sar_rda_azimuth_compression(rd_I, rd_Q, naz, nrng, &raw->params);

    /* Step 5: Azimuth IFFT -> back to image */
    for (size_t c=0; c<nrng; c++) {
        double *col_I = (double *)malloc(naz * sizeof(double));
        double *col_Q = (double *)malloc(naz * sizeof(double));
        for (size_t r=0; r<naz; r++) {
            col_I[r] = rd_I[r][c];
            col_Q[r] = rd_Q[r][c];
        }
        sar_fft_range_1d(col_I, col_Q, naz, 0);
        for (size_t r=0; r<naz; r++) {
            image->data_I[r][c] = col_I[r];
            image->data_Q[r][c] = col_Q[r];
        }
        free(col_I); free(col_Q);
    }

    /* Cleanup */
    for (size_t i=0; i<naz; i++) {
        free(work_I[i]); free(work_Q[i]);
        free(rd_I[i]); free(rd_Q[i]);
    }
    free(work_I); free(work_Q); free(rd_I); free(rd_Q);
    return 0;
}

/* RDA Step 1: Range compression in frequency domain
 * H(f_tau) = conj(FFT(chirp)) = exp(j*pi*f_tau^2/K_r) */
void sar_rda_range_compression(const sar_raw_data_t *raw,
                               double **compressed_I, double **compressed_Q)
{
    if (!raw || !compressed_I || !compressed_Q) return;
    size_t naz = raw->naz, nrng = raw->nrng;
    double K_r = raw->params.chirp_rate;
    if (K_r == 0.0) return;

    for (size_t a=0; a<naz; a++) {
        double *row_I = (double *)malloc(nrng * sizeof(double));
        double *row_Q = (double *)malloc(nrng * sizeof(double));
        memcpy(row_I, raw->data_I[a], nrng*sizeof(double));
        memcpy(row_Q, raw->data_Q[a], nrng*sizeof(double));

        /* FFT to range frequency */
        sar_fft_range_1d(row_I, row_Q, nrng, 1);

        /* Multiply by matched filter H*(f_tau) */
        double fs = raw->params.sample_rate_hz;
        for (size_t r=0; r<nrng; r++) {
            double f_tau = ((double)((int64_t)r - (int64_t)(nrng/2)) * fs) / (double)nrng;
            double phase = M_PI * f_tau * f_tau / K_r;
            double H_I = cos(phase), H_Q = -sin(phase); /* conj(H) */
            double tmp_I = row_I[r]*H_I - row_Q[r]*H_Q;
            double tmp_Q = row_I[r]*H_Q + row_Q[r]*H_I;
            row_I[r] = tmp_I; row_Q[r] = tmp_Q;
        }

        /* IFFT back to time domain */
        sar_fft_range_1d(row_I, row_Q, nrng, 0);

        memcpy(compressed_I[a], row_I, nrng*sizeof(double));
        memcpy(compressed_Q[a], row_Q, nrng*sizeof(double));
        free(row_I); free(row_Q);
    }
}


/* RDA Step 3: Range Cell Migration Correction (sinc interpolation)
 * For each azimuth frequency f_eta, shift range by delta_R(f_eta) */
void sar_rda_rcmc(double **rda_domain_I, double **rda_domain_Q,
                  size_t naz, size_t nrng,
                  const sar_params_t *sp, const sar_rda_params_t *rda_p)
{
    if (!rda_domain_I || !rda_domain_Q || !sp || !rda_p) return;
    double lambda = sp->wavelength_m, v = sp->platform_velocity_ms;
    double R0 = sp->near_range_m;
    double PRF = sp->prf_hz;
    size_t K = rda_p->rcmc_kernel_length;
    if (K < 4) K = 4;

    double *shifted_I = (double *)malloc(nrng * sizeof(double));
    double *shifted_Q = (double *)malloc(nrng * sizeof(double));

    for (size_t a=0; a<naz; a++) {
        double f_eta = ((double)((int64_t)a - (int64_t)(naz/2)) * PRF) / (double)naz;
        /* RCM shift in range bins */
        double delta_R = R0 * (1.0/sqrt(1.0 - lambda*lambda*f_eta*f_eta/(4.0*v*v)) - 1.0);
        double delta_bins = delta_R / sp->range_resolution_m;
        int shift_int = (int)floor(delta_bins);
        double frac = delta_bins - (double)shift_int;

        /* Windowed sinc interpolation (Hamming window) */
        for (size_t r=0; r<nrng; r++) {
            double sum_I=0.0, sum_Q=0.0;
            for (size_t k=0; k<=2*K; k++) {
                int idx = (int)r - (int)K + (int)k + shift_int;
                if (idx < 0 || (size_t)idx >= nrng) continue;
                double t = (double)((int)k - (int)K) + frac;
                double s = (fabs(t)<1e-12) ? 1.0 : sin(M_PI*t)/(M_PI*t);
                /* Hamming window: 0.54+0.46*cos(pi*t/K) */
                double w = 0.54 + 0.46*cos(M_PI*t/(double)K);
                if (fabs(t) >= (double)K) w = 0.0;
                sum_I += rda_domain_I[a][idx] * s * w;
                sum_Q += rda_domain_Q[a][idx] * s * w;
            }
            shifted_I[r] = sum_I; shifted_Q[r] = sum_Q;
        }
        memcpy(rda_domain_I[a], shifted_I, nrng*sizeof(double));
        memcpy(rda_domain_Q[a], shifted_Q, nrng*sizeof(double));
    }
    free(shifted_I); free(shifted_Q);
}

/* RDA Step 4: Azimuth compression with range-dependent matched filter
 * H_az(f_eta; R0) = exp(-j*pi*f_eta^2/f_R(R0)) */
void sar_rda_azimuth_compression(double **rda_domain_I, double **rda_domain_Q,
                                 size_t naz, size_t nrng,
                                 const sar_params_t *sp)
{
    if (!rda_domain_I || !rda_domain_Q || !sp) return;
    double lambda = sp->wavelength_m, v = sp->platform_velocity_ms;
    double R_near = sp->near_range_m;
    double dr = sp->range_resolution_m;
    double PRF = sp->prf_hz;

    for (size_t c=0; c<nrng; c++) {
        double R0 = R_near + (double)c * dr;
        double f_R = sar_doppler_rate(lambda, v, R0, 0.0);

        for (size_t a=0; a<naz; a++) {
            double f_eta = ((double)((int64_t)a - (int64_t)(naz/2)) * PRF) / (double)naz;
            double phase = -M_PI * f_eta * f_eta / f_R;
            double H_I = cos(phase), H_Q = sin(phase);
            double tmp_I = rda_domain_I[a][c]*H_I - rda_domain_Q[a][c]*H_Q;
            double tmp_Q = rda_domain_I[a][c]*H_Q + rda_domain_Q[a][c]*H_I;
            rda_domain_I[a][c] = tmp_I; rda_domain_Q[a][c] = tmp_Q;
        }
    }
}


/* ==================================================================
 * L5: Chirp Scaling Algorithm (CSA)
 * ================================================================== */

/* CSA main processing chain -- avoids range interpolation */
int sar_csa_process(const sar_raw_data_t *raw,
                    const sar_csa_params_t *csa_p,
                    sar_image_t *image)
{
    if (!raw || !csa_p || !image) return -1;
    size_t naz = raw->naz, nrng = raw->nrng;
    if (image->nrows != naz || image->ncols != nrng) return -2;

    /* Allocate working arrays */
    double **data_I = (double **)malloc(naz * sizeof(double *));
    double **data_Q = (double **)malloc(naz * sizeof(double *));
    for (size_t i=0; i<naz; i++) {
        data_I[i] = (double *)malloc(nrng * sizeof(double));
        data_Q[i] = (double *)malloc(nrng * sizeof(double));
        memcpy(data_I[i], raw->data_I[i], nrng*sizeof(double));
        memcpy(data_Q[i], raw->data_Q[i], nrng*sizeof(double));
    }

    /* Step 1: Azimuth FFT -> range-Doppler domain */
    for (size_t c=0; c<nrng; c++) {
        double *col_I = (double *)malloc(naz * sizeof(double));
        double *col_Q = (double *)malloc(naz * sizeof(double));
        for (size_t r=0; r<naz; r++) { col_I[r]=data_I[r][c]; col_Q[r]=data_Q[r][c]; }
        sar_fft_range_1d(col_I, col_Q, naz, 1);
        for (size_t r=0; r<naz; r++) { data_I[r][c]=col_I[r]; data_Q[r][c]=col_Q[r]; }
        free(col_I); free(col_Q);
    }

    /* Step 2: Chirp scaling phase multiply */
    sar_csa_chirp_scaling(data_I, data_Q, naz, nrng, &raw->params, csa_p);

    /* Step 3: Range FFT -> 2D frequency domain */
    for (size_t r=0; r<naz; r++) {
        sar_fft_range_1d(data_I[r], data_Q[r], nrng, 1);
    }

    /* Step 4: Range compression + bulk RCMC + SRC (combined phase multiply)
     * H_rc(f_tau, f_eta) = exp(-j*pi*f_tau^2/(K_r*(1+Cs))) * exp(-j*4pi*Rref*...)
     * Simplified: apply range compression phase */
    double K_r = raw->params.chirp_rate, fs = raw->params.sample_rate_hz;
    for (size_t r=0; r<naz; r++) {
        for (size_t c=0; c<nrng; c++) {
            double f_tau = ((double)((int64_t)c - (int64_t)(nrng/2)) * fs) / (double)nrng;
            double phase = -M_PI*f_tau*f_tau/K_r;
            double H_I=cos(phase), H_Q=sin(phase);
            double tI=data_I[r][c]*H_I - data_Q[r][c]*H_Q;
            double tQ=data_I[r][c]*H_Q + data_Q[r][c]*H_I;
            data_I[r][c]=tI; data_Q[r][c]=tQ;
        }
    }

    /* Step 5: Range IFFT */
    for (size_t r=0; r<naz; r++) {
        sar_fft_range_1d(data_I[r], data_Q[r], nrng, 0);
    }

    /* Step 6: Azimuth compression + residual phase */
    double lambda = raw->params.wavelength_m, v = raw->params.platform_velocity_ms;
    double PRF = raw->params.prf_hz;
    for (size_t c=0; c<nrng; c++) {
        double R0 = raw->params.near_range_m + (double)c*raw->params.range_resolution_m;
        double f_R = sar_doppler_rate(lambda, v, R0, 0.0);
        for (size_t r=0; r<naz; r++) {
            double f_eta = ((double)((int64_t)r - (int64_t)(naz/2)) * PRF) / (double)naz;
            double phase = -M_PI*f_eta*f_eta/f_R;
            double H_I=cos(phase), H_Q=sin(phase);
            double tI=data_I[r][c]*H_I - data_Q[r][c]*H_Q;
            double tQ=data_I[r][c]*H_Q + data_Q[r][c]*H_I;
            data_I[r][c]=tI; data_Q[r][c]=tQ;
        }
    }

    /* Step 7: Azimuth IFFT */
    for (size_t c=0; c<nrng; c++) {
        double *col_I=(double*)malloc(naz*sizeof(double)), *col_Q=(double*)malloc(naz*sizeof(double));
        for (size_t r=0; r<naz; r++) { col_I[r]=data_I[r][c]; col_Q[r]=data_Q[r][c]; }
        sar_fft_range_1d(col_I, col_Q, naz, 0);
        for (size_t r=0; r<naz; r++) { image->data_I[r][c]=col_I[r]; image->data_Q[r][c]=col_Q[r]; }
        free(col_I); free(col_Q);
    }

    for (size_t i=0; i<naz; i++) { free(data_I[i]); free(data_Q[i]); }
    free(data_I); free(data_Q);
    return 0;
}

/* Chirp scaling phase: equalizes RCM for all ranges
 * Phi_cs = exp(j*pi*K_s*Cs*(tau-tau_ref)^2) */
void sar_csa_chirp_scaling(double **data_I, double **data_Q,
                           size_t naz, size_t nrng,
                           const sar_params_t *sp, const sar_csa_params_t *csa_p)
{
    if (!data_I || !data_Q || !sp || !csa_p) return;
    /* Chirp scaling factor -- simplified for broadside case */
    double K_r = sp->chirp_rate, lambda = sp->wavelength_m, v = sp->platform_velocity_ms;
    double PRF = sp->prf_hz, fs = sp->sample_rate_hz;
    double dt = 1.0/fs;
    double R_ref = sp->near_range_m;

    for (size_t r=0; r<naz; r++) {
        double f_eta = ((double)((int64_t)r - (int64_t)(naz/2)) * PRF) / (double)naz;
        double D_f = sqrt(1.0 - lambda*lambda*f_eta*f_eta/(4.0*v*v));
        double Cs = 1.0/D_f - 1.0; /* curvature factor */
        double K_s = K_r / (1.0 - K_r*2.0*lambda*R_ref*Cs*Cs*Cs/(SAR_C*D_f*D_f*D_f));

        for (size_t c=0; c<nrng; c++) {
            double t = ((double)((int64_t)c - (int64_t)(nrng/2)) * dt);
            double tau_ref = 2.0*R_ref/(SAR_C*D_f);
            double dt2 = (t-tau_ref)*(t-tau_ref);
            double phase = M_PI*K_s*Cs*dt2;
            double phI=cos(phase), phQ=sin(phase);
            double tI=data_I[r][c]*phI - data_Q[r][c]*phQ;
            double tQ=data_I[r][c]*phQ + data_Q[r][c]*phI;
            data_I[r][c]=tI; data_Q[r][c]=tQ;
        }
    }
}


/* ==================================================================
 * L5: omega-k (Wavenumber) Algorithm
 * ================================================================== */

/* omega-k is the most accurate 2D frequency-domain SAR algorithm.
 * It uses Stolt interpolation to achieve exact focusing at all ranges.
 * Steps: 2D FFT -> ref multiply -> Stolt interp -> 2D IFFT */
int sar_omegak_process(const sar_raw_data_t *raw, double R_ref, sar_image_t *image)
{
    if (!raw || !image) return -1;
    size_t naz=raw->naz, nrng=raw->nrng;
    if (image->nrows!=naz || image->ncols!=nrng) return -2;

    double **data_I=(double**)malloc(naz*sizeof(double*));
    double **data_Q=(double**)malloc(naz*sizeof(double*));
    for (size_t i=0; i<naz; i++){
        data_I[i]=(double*)malloc(nrng*sizeof(double));
        data_Q[i]=(double*)malloc(nrng*sizeof(double));
        memcpy(data_I[i],raw->data_I[i],nrng*sizeof(double));
        memcpy(data_Q[i],raw->data_Q[i],nrng*sizeof(double));
    }

    sar_fft2d((const double**)data_I,(const double**)data_Q,naz,nrng,data_I,data_Q,1);
    sar_omegak_ref_multiply(data_I,data_Q,naz,nrng,&raw->params,R_ref);
    sar_omegak_stolt(data_I,data_Q,naz,nrng,&raw->params);
    sar_fft2d((const double**)data_I,(const double**)data_Q,naz,nrng,data_I,data_Q,0);

    for (size_t i=0; i<naz; i++){
        memcpy(image->data_I[i],data_I[i],nrng*sizeof(double));
        memcpy(image->data_Q[i],data_Q[i],nrng*sizeof(double));
        free(data_I[i]); free(data_Q[i]);
    }
    free(data_I); free(data_Q);
    return 0;
}

/* omega-k reference function multiply:
 * H_ref(f_tau,f_eta)=exp(j*R_ref*sqrt((4pi(f0+f_tau)/c)^2-(2pi*f_eta/v)^2)) */
void sar_omegak_ref_multiply(double **data_I, double **data_Q,
                             size_t naz, size_t nrng,
                             const sar_params_t *sp, double R_ref)
{
    if (!data_I || !data_Q || !sp) return;
    double f0=sp->carrier_freq_hz, v=sp->platform_velocity_ms;
    double PRF=sp->prf_hz, fs=sp->sample_rate_hz;
    double four_pi_over_c = 4.0*M_PI/SAR_C;

    for (size_t r=0; r<naz; r++){
        double f_eta = ((double)((int64_t)r-(int64_t)(naz/2))*PRF)/(double)naz;
        double kx = 2.0*M_PI*f_eta/v;
        double kx2 = kx*kx;
        for (size_t c=0; c<nrng; c++){
            double f_tau = ((double)((int64_t)c-(int64_t)(nrng/2))*fs)/(double)nrng;
            double kr = four_pi_over_c*(f0+f_tau);
            double krc2 = kr*kr - kx2;
            if (krc2 < 0.0) krc2 = 0.0;
            double phase = R_ref * sqrt(krc2);
            double H_I=cos(phase), H_Q=sin(phase);
            double tI=data_I[r][c]*H_I - data_Q[r][c]*H_Q;
            double tQ=data_I[r][c]*H_Q + data_Q[r][c]*H_I;
            data_I[r][c]=tI; data_Q[r][c]=tQ;
        }
    }
}

/* Stolt interpolation for omega-k: maps (f_tau,f_eta) -> (k_r,k_x) */
void sar_omegak_stolt(double **data_I, double **data_Q,
                      size_t naz, size_t nrng,
                      const sar_params_t *sp)
{
    if (!data_I || !data_Q || !sp) return;
    double f0=sp->carrier_freq_hz, v=sp->platform_velocity_ms;
    double PRF=sp->prf_hz, fs=sp->sample_rate_hz;
    double four_pi_over_c = 4.0*M_PI/SAR_C;

    double *new_I=(double*)malloc(nrng*sizeof(double));
    double *new_Q=(double*)malloc(nrng*sizeof(double));

    for (size_t r=0; r<naz; r++){
        double f_eta = ((double)((int64_t)r-(int64_t)(naz/2))*PRF)/(double)naz;
        double kx = 2.0*M_PI*f_eta/v;
        double kx2 = kx*kx;

        for (size_t c=0; c<nrng; c++){
            double f_tau = ((double)((int64_t)c-(int64_t)(nrng/2))*fs)/(double)nrng;
            double kr = four_pi_over_c*(f0+f_tau);
            double krc = sqrt(fmax(0.0, kr*kr - kx2));
            /* Map krc back to new f_tau grid: linear interpolation */
            double f_tau_new = krc/four_pi_over_c - f0;
            double idx_f = (f_tau_new/(fs/(double)nrng)) + (double)(nrng/2);
            int i0 = (int)floor(idx_f);
            double frac = idx_f - (double)i0;
            if (i0>=0 && (size_t)(i0+1)<nrng){
                new_I[c] = data_I[r][i0]*(1.0-frac) + data_I[r][i0+1]*frac;
                new_Q[c] = data_Q[r][i0]*(1.0-frac) + data_Q[r][i0+1]*frac;
            } else if (i0>=0 && (size_t)i0<nrng){
                new_I[c]=data_I[r][i0]; new_Q[c]=data_Q[r][i0];
            } else {
                new_I[c]=0.0; new_Q[c]=0.0;
            }
        }
        memcpy(data_I[r],new_I,nrng*sizeof(double));
        memcpy(data_Q[r],new_Q,nrng*sizeof(double));
    }
    free(new_I); free(new_Q);
}


/* ==================================================================
 * L5: Backprojection Algorithm (BP)
 * ================================================================== */

/* BP is a time-domain algorithm ideal for non-linear trajectories.
 * For each pixel (xi,yj), coherently sum range-compressed echoes
 * from all pulses with appropriate phase compensation.
 * Complexity: O(N_az*N_rng*N_pixels) -- slow but exact. */
int sar_bp_process(const sar_raw_data_t *raw,
                   const sar_bp_params_t *bp_p,
                   const double *traj_x, const double *traj_y, const double *traj_z,
                   sar_image_t *image)
{
    if (!raw || !bp_p || !traj_x || !traj_y || !traj_z || !image) return -1;

    double dx = bp_p->pixel_size_x, dy = bp_p->pixel_size_y;
    double lambda = raw->params.wavelength_m;
    double dr = raw->range_sampling_interval;
    double r0 = raw->params.near_range_m;

    for (size_t py=0; py<bp_p->grid_ny; py++){
        double y = bp_p->grid_y_min + (double)py*dy;
        for (size_t px=0; px<bp_p->grid_nx; px++){
            double x = bp_p->grid_x_min + (double)px*dx;
            double sum_I=0.0, sum_Q=0.0;

            for (size_t a=0; a<raw->naz; a++){
                double dx2=(traj_x[a]-x), dy2=(traj_y[a]-y), dz2=traj_z[a];
                double R = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
                double range_bin = (R - r0)/dr;
                int ib = (int)floor(range_bin);
                double frac = range_bin - (double)ib;

                /* Linear interpolation of range-compressed data */
                double sI=0.0, sQ=0.0;
                if (ib>=0 && (size_t)(ib+1)<raw->nrng){
                    sI = raw->data_I[a][ib]*(1.0-frac) + raw->data_I[a][ib+1]*frac;
                    sQ = raw->data_Q[a][ib]*(1.0-frac) + raw->data_Q[a][ib+1]*frac;
                } else if (ib>=0 && (size_t)ib<raw->nrng){
                    sI = raw->data_I[a][ib]; sQ = raw->data_Q[a][ib];
                }

                /* Phase compensation: exp(+j*4pi*R/lambda) */
                double phase = 4.0*M_PI*R/lambda;
                double cph=cos(phase), sph=sin(phase);
                sum_I += sI*cph - sQ*sph;
                sum_Q += sI*sph + sQ*cph;
            }
            image->data_I[py][px] = sum_I;
            image->data_Q[py][px] = sum_Q;
        }
    }
    return 0;
}

/* ==================================================================
 * L5: SPECAN (Spectral Analysis) for ScanSAR/Burst
 * ================================================================== */

/* SPECAN: azimuth deramp + FFT, trading resolution for speed.
 * Suitable for burst-mode (ScanSAR, TOPS) where azimuth extent is limited. */
int sar_specan_process(const sar_raw_data_t *raw, size_t burst_len, sar_image_t *image)
{
    if (!raw || !image || burst_len==0 || burst_len>raw->naz) return -1;

    double lambda = raw->params.wavelength_m, v = raw->params.platform_velocity_ms;
    double R0 = raw->params.near_range_m;
    double f_R = sar_doppler_rate(lambda, v, R0, 0.0);
    double PRF = raw->params.prf_hz;
    double d_eta = 1.0/PRF;

    for (size_t c=0; c<raw->nrng; c++){
        double R0c = R0 + (double)c*raw->params.range_resolution_m;
        double f_Rc = sar_doppler_rate(lambda, v, R0c, 0.0);

        /* Deramp: multiply by exp(-j*pi*f_R*eta^2) */
        double *az_I = (double*)malloc(burst_len*sizeof(double));
        double *az_Q = (double*)malloc(burst_len*sizeof(double));
        for (size_t a=0; a<burst_len; a++){
            double eta = ((double)((int64_t)a-(int64_t)(burst_len/2)))*d_eta;
            double deramp_phase = -M_PI*f_Rc*eta*eta;
            double dI=cos(deramp_phase), dQ=sin(deramp_phase);
            az_I[a] = raw->data_I[a][c]*dI - raw->data_Q[a][c]*dQ;
            az_Q[a] = raw->data_I[a][c]*dQ + raw->data_Q[a][c]*dI;
        }

        /* FFT for azimuth focusing */
        sar_fft_range_1d(az_I, az_Q, burst_len, 1);

        for (size_t a=0; a<burst_len && a<image->nrows; a++){
            image->data_I[a][c] = az_I[a];
            image->data_Q[a][c] = az_Q[a];
        }
        free(az_I); free(az_Q);
    }
    return 0;
}


/* ==================================================================
 * L5: Phase Gradient Autofocus (PGA) -- non-parametric, iterative
 * Reference: Wahl et al., IEEE TASSP 1994
 * ================================================================== */
int sar_pga_autofocus(sar_image_t *image, int niter)
{
    if (!image || niter<1) return -1;
    size_t nrows=image->nrows, ncols=image->ncols;

    for (int iter=0; iter<niter; iter++){
        for (size_t c=0; c<ncols; c++){
            double max_mag=0.0; size_t max_r=0;
            for (size_t r=0; r<nrows; r++){
                double mag = image->data_I[r][c]*image->data_I[r][c]
                           + image->data_Q[r][c]*image->data_Q[r][c];
                if (mag>max_mag){ max_mag=mag; max_r=r; }
            }
            if (max_mag <= 0.0) continue;

            double *shift_I=(double*)malloc(nrows*sizeof(double));
            double *shift_Q=(double*)malloc(nrows*sizeof(double));
            for (size_t r=0; r<nrows; r++){
                size_t sr = (r + max_r) % nrows;
                shift_I[r]=image->data_I[sr][c];
                shift_Q[r]=image->data_Q[sr][c];
            }

            /* Window around center */
            size_t wh=nrows/4, ws=nrows/2-wh, we=nrows/2+wh;
            for (size_t r=0; r<nrows; r++)
                if (r<ws || r>we){ shift_I[r]=0.0; shift_Q[r]=0.0; }

            sar_fft_range_1d(shift_I, shift_Q, nrows, 1);

            /* Phase gradient: Delta_phi(k)=arg(X[k]*conj(X[k-1])) */
            double *pg=(double*)malloc(nrows*sizeof(double));
            pg[0]=0.0;
            for (size_t r=1; r<nrows; r++){
                double pI=shift_I[r]*shift_I[r-1]+shift_Q[r]*shift_Q[r-1];
                double pQ=shift_Q[r]*shift_I[r-1]-shift_I[r]*shift_Q[r-1];
                pg[r]=atan2(pQ,pI);
            }
            double *pe=(double*)malloc(nrows*sizeof(double));
            pe[0]=0.0;
            for (size_t r=1; r<nrows; r++) pe[r]=pe[r-1]+pg[r];

            for (size_t r=0; r<nrows; r++){
                double cph=cos(-pe[r]), sph=sin(-pe[r]);
                double tI=image->data_I[r][c]*cph-image->data_Q[r][c]*sph;
                double tQ=image->data_I[r][c]*sph+image->data_Q[r][c]*cph;
                image->data_I[r][c]=tI; image->data_Q[r][c]=tQ;
            }
            free(shift_I); free(shift_Q); free(pg); free(pe);
        }
    }
    return 0;
}


/* ==================================================================
 * L5: Map Drift (MD) Autofocus
 * Split azimuth spectrum into sub-apertures, cross-correlate
 * resulting images to estimate quadratic phase error (defocus).
 * ================================================================== */
int sar_mapdrift_autofocus(sar_image_t *image)
{
    if (!image) return -1;
    size_t nrows=image->nrows, ncols=image->ncols;
    if (nrows<8) return -1;

    size_t half=nrows/2;
    double *s1_I=(double*)malloc(half*sizeof(double));
    double *s1_Q=(double*)malloc(half*sizeof(double));
    double *s2_I=(double*)malloc(half*sizeof(double));
    double *s2_Q=(double*)malloc(half*sizeof(double));

    for (size_t c=0; c<ncols; c++){
        for (size_t r=0; r<half; r++){
            s1_I[r]=image->data_I[r][c]; s1_Q[r]=image->data_Q[r][c];
            s2_I[r]=image->data_I[r+half][c]; s2_Q[r]=image->data_Q[r+half][c];
        }
        sar_fft_range_1d(s1_I,s1_Q,half,0); sar_fft_range_1d(s2_I,s2_Q,half,0);
        sar_fft_range_1d(s1_I,s1_Q,half,1); sar_fft_range_1d(s2_I,s2_Q,half,1);

        for (size_t r=0; r<half; r++){
            double pI=s1_I[r]*s2_I[r]+s1_Q[r]*s2_Q[r];
            double pQ=s1_Q[r]*s2_I[r]-s1_I[r]*s2_Q[r];
            s1_I[r]=pI; s1_Q[r]=pQ;
        }
        sar_fft_range_1d(s1_I,s1_Q,half,0);

        double max_corr=0.0; size_t peak=0;
        for (size_t r=0; r<half; r++){
            double mag=s1_I[r]*s1_I[r]+s1_Q[r]*s1_Q[r];
            if (mag>max_corr){ max_corr=mag; peak=r; }
        }
        int shift=(int)peak;
        if (shift>(int)half/2) shift-=(int)half;
        double df=(double)shift/(double)half;

        for (size_t r=0; r<nrows; r++){
            double t=((double)((int64_t)r-(int64_t)(nrows/2)))/(double)nrows;
            double phase=-M_PI*df*t*t;
            double cp=cos(phase), sp=sin(phase);
            double tI=image->data_I[r][c]*cp-image->data_Q[r][c]*sp;
            double tQ=image->data_I[r][c]*sp+image->data_Q[r][c]*cp;
            image->data_I[r][c]=tI; image->data_Q[r][c]=tQ;
        }
    }
    free(s1_I); free(s1_Q); free(s2_I); free(s2_Q);
    return 0;
}


/* ==================================================================
 * L6: Point Target Impulse Response Analysis
 * Computes -3dB resolution, PSLR, ISLR from focused point target.
 * ================================================================== */
void sar_analyze_impulse_response(const sar_image_t *image,
                                  size_t r0, size_t c0, size_t win,
                                  sar_impulse_response_t *ir)
{
    if (!image || !ir) return;
    memset(ir, 0, sizeof(*ir));

    /* Find actual peak within window */
    double peak_mag=0.0; size_t pr=r0, pc=c0;
    for (size_t r=(r0>win?r0-win:0); r<r0+win && r<image->nrows; r++){
        for (size_t c=(c0>win?c0-win:0); c<c0+win && c<image->ncols; c++){
            double mag=sqrt(image->data_I[r][c]*image->data_I[r][c]
                          + image->data_Q[r][c]*image->data_Q[r][c]);
            if (mag>peak_mag){ peak_mag=mag; pr=r; pc=c; }
        }
    }
    ir->peak_range_idx=(double)pc; ir->peak_azimuth_idx=(double)pr;
    ir->peak_magnitude=peak_mag;
    if (pr<image->nrows && pc<image->ncols)
        ir->peak_phase_rad=atan2(image->data_Q[pr][pc],image->data_I[pr][pc]);

    /* Range profile through peak */
    if (pr<image->nrows){
        double *prof=(double*)malloc(image->ncols*sizeof(double));
        for (size_t c=0; c<image->ncols; c++)
            prof[c]=sqrt(image->data_I[pr][c]*image->data_I[pr][c]
                       + image->data_Q[pr][c]*image->data_Q[pr][c]);

        double thresh=peak_mag/M_SQRT2;
        size_t left=pc, right=pc;
        while (left>0 && prof[left]>thresh) left--;
        while (right<image->ncols-1 && prof[right]>thresh) right++;
        ir->range_resolution_3dB_m=(double)(right-left)*image->range_pixel_spacing_m;

        double max_sl=0.0;
        size_t ml=(pc>20?pc-20:0), mr=(pc+20<image->ncols?pc+20:image->ncols);
        for (size_t c=0; c<image->ncols; c++){
            if (c>=ml && c<=mr) continue;
            if (prof[c]>max_sl) max_sl=prof[c];
        }
        ir->range_pslr_db=(max_sl>0.0&&peak_mag>0.0)?20.0*log10(max_sl/peak_mag):-999.0;

        double Em=0.0, Es=0.0;
        for (size_t c=0; c<image->ncols; c++){
            double p=prof[c]*prof[c];
            if (c>=ml && c<=mr) Em+=p; else Es+=p;
        }
        ir->range_islr_db=(Em>0.0&&Es>0.0)?10.0*log10(Es/Em):-999.0;
        free(prof);
    }

    /* Azimuth profile through peak */
    if (pc<image->ncols){
        double *prof=(double*)malloc(image->nrows*sizeof(double));
        for (size_t r=0; r<image->nrows; r++)
            prof[r]=sqrt(image->data_I[r][pc]*image->data_I[r][pc]
                       + image->data_Q[r][pc]*image->data_Q[r][pc]);

        double thresh=peak_mag/M_SQRT2;
        size_t top=pr, bot=pr;
        while (top>0 && prof[top]>thresh) top--;
        while (bot<image->nrows-1 && prof[bot]>thresh) bot++;
        ir->azimuth_resolution_3dB_m=(double)(bot-top)*image->azimuth_pixel_spacing_m;

        double max_sl=0.0;
        size_t mt=(pr>20?pr-20:0), mb=(pr+20<image->nrows?pr+20:image->nrows);
        for (size_t r=0; r<image->nrows; r++){
            if (r>=mt && r<=mb) continue;
            if (prof[r]>max_sl) max_sl=prof[r];
        }
        ir->azimuth_pslr_db=(max_sl>0.0&&peak_mag>0.0)?20.0*log10(max_sl/peak_mag):-999.0;

        double Em=0.0, Es=0.0;
        for (size_t r=0; r<image->nrows; r++){
            double p=prof[r]*prof[r];
            if (r>=mt && r<=mb) Em+=p; else Es+=p;
        }
        ir->azimuth_islr_db=(Em>0.0&&Es>0.0)?10.0*log10(Es/Em):-999.0;
        free(prof);
    }
}
