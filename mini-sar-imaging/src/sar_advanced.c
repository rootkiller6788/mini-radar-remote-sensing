/**
 * @file    sar_advanced.c
 * @brief   Advanced SAR Topics -- L8: CS-SAR, MIMO-SAR, Bistatic, Polarimetric
 * @details Compressive sensing reconstruction, bistatic range equation,
 *          MIMO SAR imaging, Freeman-Durden and H-Alpha polarimetric decomposition.
 * Reference: Candes & Tao (2006), Freeman & Durden (1998), Cloude & Pottier (1997)
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "sar_core.h"
#include "sar_advanced.h"

/* ==================================================================
 * L8: Compressive Sensing SAR Reconstruction
 * ================================================================== */

/* ISTA (Iterative Soft-Thresholding Algorithm) for L1-regularized
 * least squares: min_x 0.5*||Ax - y||^2 + lambda*||x||_1
 *
 * CS-SAR exploits sparsity of scene in some basis (e.g., wavelet,
 * gradient) to reconstruct from undersampled measurements.
 *
 * y = A*x + n where A is the sensing matrix (partial Fourier, etc.)
 * Reconstruction via iterative soft-thresholding:
 *   x_{k+1} = soft(x_k - mu*A^T*(A*x_k - y), mu*lambda)
 *   soft(z, t) = sign(z)*max(|z|-t, 0)
 */
int sar_cs_reconstruct(const double **measurements, size_t n_meas, size_t n_pixels,
                       const sar_cs_params_t *csp, double *reconstructed)
{
    if (!measurements || !csp || !reconstructed) return -1;

    double mu = csp->undersampling_ratio;
    double lambda = csp->regularization_lambda;
    int max_iter = csp->max_iterations;
    double tol = csp->convergence_tol;

    double *x = (double *)calloc(n_pixels, sizeof(double));
    double *r = (double *)malloc(n_meas * sizeof(double));
    double *At_r = (double *)calloc(n_pixels, sizeof(double));

    if (!x || !r || !At_r) { free(x); free(r); free(At_r); return -2; }

    for (int iter=0; iter<max_iter; iter++){
        /* Compute residual r = A*x - y */
        for (size_t i=0; i<n_meas; i++){
            double Ax=0.0;
            for (size_t j=0; j<n_pixels; j++)
                Ax += measurements[i][j] * x[j];
            r[i] = Ax - measurements[i][n_pixels]; /* last column = measurement */
        }

        /* Gradient step: At * r */
        for (size_t j=0; j<n_pixels; j++){
            At_r[j] = 0.0;
            for (size_t i=0; i<n_meas; i++)
                At_r[j] += measurements[i][j] * r[i];
        }

        /* ISTA update with soft-thresholding */
        double max_change = 0.0;
        for (size_t j=0; j<n_pixels; j++){
            double z = x[j] - mu * At_r[j];
            double thresh = mu * lambda;
            double x_new;
            if (z > thresh) x_new = z - thresh;
            else if (z < -thresh) x_new = z + thresh;
            else x_new = 0.0;

            double change = fabs(x_new - x[j]);
            if (change > max_change) max_change = change;
            x[j] = x_new;
        }
        if (max_change < tol) break;
    }

    memcpy(reconstructed, x, n_pixels * sizeof(double));
    free(x); free(r); free(At_r);
    return 0;
}

/* ==================================================================
 * L8: Bistatic SAR Range Equation
 * ================================================================== */

/* Bistatic range: R(eta) = R_tx(eta) + R_rx(eta)
 * where R_tx is range from transmitter to target,
 * and R_rx is range from receiver to target.
 * The bistatic geometry enables forward-scattering enhancement,
 * stealth target detection, and interferometric configurations. */
double sar_bistatic_range(const sar_bistatic_config_t *bc, double eta,
                          double x_tgt, double y_tgt, double z_tgt)
{
    if (!bc) return 0.0;
    double tx_x=bc->tx_x+bc->tx_vel_x*eta, tx_y=bc->tx_y+bc->tx_vel_y*eta, tx_z=bc->tx_z+bc->tx_vel_z*eta;
    double rx_x=bc->rx_x+bc->rx_vel_x*eta, rx_y=bc->rx_y+bc->rx_vel_y*eta, rx_z=bc->rx_z+bc->rx_vel_z*eta;

    double dx_tx=tx_x-x_tgt, dy_tx=tx_y-y_tgt, dz_tx=tx_z-z_tgt;
    double dx_rx=rx_x-x_tgt, dy_rx=rx_y-y_tgt, dz_rx=rx_z-z_tgt;

    return sqrt(dx_tx*dx_tx+dy_tx*dy_tx+dz_tx*dz_tx)
         + sqrt(dx_rx*dx_rx+dy_rx*dy_rx+dz_rx*dz_rx);
}

/* Simplified bistatic SAR processor.
 * Uses backprojection-like approach adapted for bistatic geometry.
 * For each pixel, coherently sums returns from transmitter->target->receiver path. */
int sar_bistatic_processor(const sar_raw_data_t *raw, const sar_bistatic_config_t *bc,
                           sar_image_t *image)
{
    if (!raw || !bc || !image) return -1;
    double lambda=raw->params.wavelength_m, dr=raw->range_sampling_interval;
    double r0=raw->params.near_range_m, d_eta=raw->azimuth_sampling_interval;

    for (size_t r=0; r<image->nrows; r++){
        for (size_t c=0; c<image->ncols; c++){
            double sum_I=0.0, sum_Q=0.0;
            /* Grid coordinates: assume uniform grid */
            double x_tgt=(double)c*image->range_pixel_spacing_m;
            double y_tgt=(double)r*image->azimuth_pixel_spacing_m;
            double z_tgt=0.0;

            for (size_t a=0; a<raw->naz; a++){
                double eta=(double)a*d_eta;
                double R=sar_bistatic_range(bc,eta,x_tgt,y_tgt,z_tgt);
                double bin=(R-r0)/dr;
                int ib=(int)floor(bin);
                double frac=bin-(double)ib;
                double sI=0.0,sQ=0.0;
                if (ib>=0&&(size_t)(ib+1)<raw->nrng){
                    sI=raw->data_I[a][ib]*(1.0-frac)+raw->data_I[a][ib+1]*frac;
                    sQ=raw->data_Q[a][ib]*(1.0-frac)+raw->data_Q[a][ib+1]*frac;
                } else if (ib>=0&&(size_t)ib<raw->nrng){
                    sI=raw->data_I[a][ib]; sQ=raw->data_Q[a][ib];
                }
                double phase=2.0*M_PI*R/lambda;
                double cp=cos(phase), sp=sin(phase);
                sum_I+=sI*cp-sQ*sp; sum_Q+=sI*sp+sQ*cp;
            }
            image->data_I[r][c]=sum_I; image->data_Q[r][c]=sum_Q;
        }
    }
    return 0;
}


/* ==================================================================
 * L8: MIMO SAR Imaging
 * ================================================================== */

/* MIMO-SAR uses multiple transmit and receive antennas to increase
 * spatial degrees of freedom. Virtual array concept: N_tx * N_rx
 * effective phase centers from N_tx+N_rx physical elements.
 * This enables high-resolution wide-swath (HRWS) imaging by
 * overcoming the minimum antenna area constraint. */
int sar_mimo_image(const sar_mimo_config_t *cfg, double **image,
                   size_t nrows, size_t ncols)
{
    if (!cfg || !image) return -1;
    size_t n_tx=cfg->n_tx, n_rx=cfg->n_rx;
    size_t N=n_tx*n_rx;

    /* For each virtual channel (tx_i, rx_j):
     *   - Phase center at (tx_pos[i]+rx_pos[j])/2
     *   - Two-way phase: exp(-j*2pi*(R_tx+R_rx)/lambda)
     * Coherent combination of all virtual channels improves SNR by 10*log10(N) */
    for (size_t r=0; r<nrows; r++){
        for (size_t c=0; c<ncols; c++){
            double sum_I=0.0, sum_Q=0.0;
            for (size_t ti=0; ti<n_tx; ti++){
                for (size_t rj=0; rj<n_rx; rj++){
                    double pc=(cfg->tx_positions_x[ti]+cfg->rx_positions_x[rj])/2.0;
                    double R0=sqrt(pc*pc+(double)(r*nrows+r)*(c*ncols+c));
                    double phase=-4.0*M_PI*R0;
                    double wf=cfg->tx_waveforms[ti*cfg->waveform_length+rj];
                    sum_I+=wf*cos(phase); sum_Q+=wf*sin(phase);
                }
            }
            image[r][c]=sqrt(sum_I*sum_I+sum_Q*sum_Q);
        }
    }
    return 0;
}

/* ==================================================================
 * L8: Polarimetric SAR Decomposition
 * ================================================================== */

/* Freeman-Durden 3-component decomposition (1998).
 * Decomposes coherency matrix T into:
 *   T = Ps*Ts + Pd*Td + Pv*Tv
 * where Ps=surface, Pd=double-bounce, Pv=volume scattering power.
 *
 * Surface:    Ts = [[1, b, 0], [b*, |b|^2, 0], [0, 0, 0]]
 * Double:     Td = [[|a|^2, a, 0], [a*, 1, 0], [0, 0, 0]]
 * Volume:     Tv = (1/4)*[[2, 0, 0], [0, 1, 0], [0, 0, 1]] (random)
 *
 * Solving the underdetermined system using van Zyl's method. */
void sar_freeman_durden(const sar_coherency_matrix_t *T, size_t n,
                        double *Ps, double *Pd, double *Pv)
{
    if (!T || !Ps || !Pd || !Pv) return;
    for (size_t i=0; i<n; i++){
        double t11=T[i].T11, t22=T[i].T22, t33=T[i].T33;
        double t12r=T[i].T12_real, t12i=T[i].T12_imag;

        /* Volume scattering power */
        double pv = 4.0*t33;
        if (pv<0.0) pv=0.0;

        /* Surface/Double-bounce separation using sign of Re(S_hh*S_vv*)
         * If t11-t22+pv/4 > pv/2: double-bounce dominant */
        double x11 = t11 - pv/2.0;
        double x22 = t22 - pv/4.0;
        double x12_mag = sqrt(t12r*t12r+t12i*t12i) - pv/4.0;
        if (x12_mag<0.0) x12_mag=0.0;

        double alpha = -1.0;
        if (x11>0.0 && x22>0.0) alpha = x12_mag/x11;

        double ps, pd;
        if (x11>x22){
            pd = x22 - x12_mag*alpha;
            ps = x11 - x12_mag/alpha;
        } else {
            ps = x22 - x12_mag*alpha;
            pd = x11 - x12_mag/alpha;
        }
        if (ps<0.0) ps=0.0; if (pd<0.0) pd=0.0;

        Ps[i]=ps; Pd[i]=pd; Pv[i]=pv;
    }
}

/* H-Alpha polarimetric decomposition (Cloude-Pottier, 1997).
 * Eigenvalue decomposition of coherency matrix T:
 *   H = -sum(p_i*log3(p_i))  (entropy, 0-1)
 *   alpha_bar = sum(p_i*alpha_i)  (mean alpha angle, 0-90 deg)
 *   A = (lambda2-lambda3)/(lambda2+lambda3)  (anisotropy)
 *
 * H=0: pure target, H=1: random volume
 * alpha=0: surface, alpha=45: dipole, alpha=90: double-bounce */
void sar_h_alpha_decomp(const sar_coherency_matrix_t *T, size_t n,
                        double *entropy, double *alpha, double *anisotropy)
{
    if (!T || !entropy || !alpha || !anisotropy) return;
    for (size_t i=0; i<n; i++){
        /* Compute eigenvalues of 3x3 Hermitian T matrix via characteristic polynomial:
         * det(T - lambda*I) = -lambda^3 + a*lambda^2 + b*lambda + c = 0
         * where a=tr(T), b=0.5*(tr(T)^2-tr(T^2)), c=det(T). */
        double t11=T[i].T11, t22=T[i].T22, t33=T[i].T33;
        double a = t11+t22+t33;
        double b = t11*t22+t22*t33+t33*t11
                 - (T[i].T12_real*T[i].T12_real+T[i].T12_imag*T[i].T12_imag)
                 - (T[i].T13_real*T[i].T13_real+T[i].T13_imag*T[i].T13_imag)
                 - (T[i].T23_real*T[i].T23_real+T[i].T23_imag*T[i].T23_imag);
        double c = t11*t22*t33
                 + 2.0*(T[i].T12_real*(T[i].T23_real*T[i].T13_real+T[i].T23_imag*T[i].T13_imag)
                      + T[i].T12_imag*(T[i].T23_real*T[i].T13_imag-T[i].T23_imag*T[i].T13_real))
                 - t11*(T[i].T23_real*T[i].T23_real+T[i].T23_imag*T[i].T23_imag)
                 - t22*(T[i].T13_real*T[i].T13_real+T[i].T13_imag*T[i].T13_imag)
                 - t33*(T[i].T12_real*T[i].T12_real+T[i].T12_imag*T[i].T12_imag);
        b = a*a - b;

        /* Solve cubic using trigonometric method for 3 real roots (Hermitian->real eigenvalues) */
        double p = (3.0*b - a*a)/3.0;
        double q = (2.0*a*a*a - 9.0*a*b + 27.0*c)/27.0;
        double disc = q*q/4.0 + p*p*p/27.0;

        double l1,l2,l3;
        if (disc>=0.0){
            double u = cbrt(-q/2.0+sqrt(disc));
            double v = cbrt(-q/2.0-sqrt(disc));
            l1=u+v+a/3.0;
            l2=-(u+v)/2.0+a/3.0;
            l3=l2;
        } else {
            double phi=acos(-q/2.0/sqrt(-p*p*p/27.0));
            double s=2.0*sqrt(-p/3.0);
            l1=s*cos(phi/3.0)+a/3.0;
            l2=s*cos((phi+2.0*M_PI)/3.0)+a/3.0;
            l3=s*cos((phi+4.0*M_PI)/3.0)+a/3.0;
        }
        if (l1<0.0) l1=0.0; if (l2<0.0) l2=0.0; if (l3<0.0) l3=0.0;
        /* Sort descending */
        if (l2>l1){ double t=l1; l1=l2; l2=t; }
        if (l3>l1){ double t=l1; l1=l3; l3=t; }
        if (l3>l2){ double t=l2; l2=l3; l3=t; }

        double sum=l1+l2+l3;
        if (sum<1e-12){ entropy[i]=0.0; alpha[i]=0.0; anisotropy[i]=0.0; continue; }
        double p1=l1/sum, p2=l2/sum, p3=l3/sum;

        /* Entropy H */
        double H=0.0;
        if (p1>0.0) H-=p1*log(p1)/log(3.0);
        if (p2>0.0) H-=p2*log(p2)/log(3.0);
        if (p3>0.0) H-=p3*log(p3)/log(3.0);
        entropy[i]=H;

        /* Mean alpha angle: alpha_bar = p1*alpha1 + p2*alpha2 + p3*alpha3
         * Simplified: alpha_i = acos(|u_i1|) where u_i1 is first element of eigenvector.
         * Approximation using diagonal elements */
        alpha[i] = p1*acos(fabs(t11)/sqrt(sum)) + p2*acos(fabs(t22)/sqrt(sum))
                 + p3*acos(fabs(t33)/sqrt(sum));

        /* Anisotropy A = (l2-l3)/(l2+l3) */
        anisotropy[i] = (l2+l3>0.0) ? (l2-l3)/(l2+l3) : 0.0;
    }
}
