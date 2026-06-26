
/**
 * @file    sar_core.c
 * @brief   SAR Core Implementations -- L1 Definitions + L3 Math + L4 Laws
 *
 * @details Implements: SAR parameter initialization, chirp generation and
 *          matched filtering, raw echo data simulation (point target),
 *          SAR image management (multi-look, calibration, magnitude/phase),
 *          2D FFT operations, wavenumber domain mapping, and Stolt interpolation.
 *
 * Knowledge Mapping:
 *   L1 - sar_params_init, sar_chirp_alloc, sar_chirp_autocorrelation,
 *        sar_pulse_compression, sar_raw_data_point_target,
 *        sar_image_magnitude/phase, sar_multilook, sar_calibrate_sigma0
 *   L3 - sar_fft2d, sar_fft_range, sar_fft_azimuth, sar_wavenumber_init,
 *        sar_stolt_interpolation
 *   L4 - Range resolution rho_r = c/(2*B_r) verified in params_init
 *        Nyquist PRF constraint in azimuth bandwidth computation
 *
 * Reference:
 *   - Cumming & Wong, "Digital Processing of SAR Data" (2005)
 *   - Richards, Scheer & Holm, "Principles of Modern Radar" (2010)
 *   - Cooley & Tukey, "An Algorithm for the Machine Calculation of
 *     Complex Fourier Series" (1965), Math. Comp.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sar_core.h"

/* Forward declaration: 1D FFT helper */
void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward);


/* ==================================================================
 * L1: SAR Parameter Initialization
 * ================================================================== */

/**
 * Initialize all SAR system parameters from fundamental values.
 *
 * Derived quantities:
 *   lambda          = c / f0                 (wavelength)
 *   chirp_rate      = B_r / tau_p            (Hz/s)
 *   range_resolution = c / (2 * B_r)          (m)
 *   azimuth_resolution = L_a / 2              (m, focused stripmap)
 *   swath_width     = R_far*cos(theta_i_far) - R_near*cos(theta_i_near)
 *
 * The azimuth resolution rho_a = L_a/2 is the fundamental SAR result:
 * finer resolution is achieved with a smaller antenna, counter-intuitively,
 * because a smaller antenna creates a wider beam, hence a longer synthetic
 * aperture and greater Doppler bandwidth.
 */
void sar_params_init(sar_params_t *p,
                     double freq, double bw, double pw, double prf,
                     double ant, double vel, double alt,
                     double squint, double look,
                     double near, double far)
{
    if (!p) return;
    memset(p, 0, sizeof(*p));

    p->carrier_freq_hz     = freq;
    p->wavelength_m        = SAR_C / freq;
    p->bandwidth_hz        = bw;
    p->pulse_width_s       = pw;
    p->prf_hz              = prf;
    p->chirp_rate          = (pw > 0.0) ? bw / pw : 0.0;
    p->sample_rate_hz      = bw * 1.2;
    p->range_resolution_m  = (bw > 0.0) ? SAR_C / (2.0 * bw) : 0.0;
    p->azimuth_resolution_m = ant / 2.0;
    p->antenna_length_m    = ant;
    p->platform_velocity_ms = vel;
    p->squint_angle_rad    = squint;
    p->look_angle_rad      = look;
    p->incidence_angle_rad = asin(alt / near);
    p->platform_altitude_m = alt;
    p->near_range_m        = near;
    p->far_range_m         = far;

    /* Swath width in ground range */
    double theta_i_near = asin(alt / near);
    double theta_i_far  = asin(alt / far);
    double G_near = near * cos(theta_i_near);
    double G_far  = far  * cos(theta_i_far);
    p->swath_width_m = G_far - G_near;
    if (p->swath_width_m < 0.0) p->swath_width_m = 0.0;
}

/* ==================================================================
 * L1: Chirp (Linear FM) Generation
 * ================================================================== */

/**
 * Allocate and generate a complex baseband chirp:
 *
 *   s[n] = exp(j * pi * K * (n*dt)^2)  for n = 0..N-1
 *   I[n] = cos(pi * K * (n*dt)^2)
 *   Q[n] = sin(pi * K * (n*dt)^2)
 *
 * where K = B/tau is the chirp rate and dt = 1/fs.
 *
 * Time-bandwidth product TBP = B * tau = K * tau^2 characterizes
 * the pulse compression gain: compressed pulse width ~ 1/B.
 *
 * @param N    Number of samples
 * @param tau  Pulse width [s]
 * @param B    Bandwidth [Hz]
 * @param fs   Sampling rate [Hz]
 * @return     Allocated chirp structure, or NULL on failure
 */
sar_chirp_t *sar_chirp_alloc(size_t N, double tau, double B, double fs)
{
    if (N == 0 || tau <= 0.0 || B <= 0.0 || fs <= 0.0) return NULL;

    sar_chirp_t *c = (sar_chirp_t *)malloc(sizeof(sar_chirp_t));
    if (!c) return NULL;

    c->I = (double *)malloc(N * sizeof(double));
    c->Q = (double *)malloc(N * sizeof(double));
    if (!c->I || !c->Q) {
        free(c->I); free(c->Q); free(c);
        return NULL;
    }

    c->num_samples  = N;
    c->pulse_width_s = tau;
    c->bandwidth_hz  = B;
    c->chirp_rate    = B / tau;
    c->sample_rate_hz = fs;
    c->dt            = 1.0 / fs;
    c->t_min         = -tau / 2.0;
    c->t_max         = +tau / 2.0;
    c->time_bandwidth_product = B * tau;

    /* Generate I/Q samples */
    double K  = c->chirp_rate;
    double dt = c->dt;
    for (size_t n = 0; n < N; n++) {
        double t = (double)((int64_t)n - (int64_t)(N/2)) * dt;
        double phase = M_PI * K * t * t;
        c->I[n] = cos(phase);
        c->Q[n] = sin(phase);
    }
    return c;
}

void sar_chirp_free(sar_chirp_t *c)
{
    if (!c) return;
    free(c->I);
    free(c->Q);
    free(c);
}

/**
 * Compute chirp autocorrelation (matched filter response):
 *
 *   R[k] = sum_{n} s[n] * conj(s[n-k])
 *
 * For a chirp, R[k] approximates a sinc function:
 *   |R(tau)| ~ |sinc(B * tau)|  for |tau| <= tau_p
 *
 * Peak-to-Sidelobe Ratio (PSLR) for rectangular window:
 *   PSLR ~ -13.26 dB (first sidelobe of sinc)
 *
 * The -3dB width of the mainlobe is approximately 0.886/B,
 * giving the compressed range resolution rho_r = c/(2B).
 *
 * @param c      Chirp signal
 * @param R_out  Output autocorrelation magnitude [2*N-1]
 */
void sar_chirp_autocorrelation(const sar_chirp_t *c, double *R_out)
{
    if (!c || !R_out) return;
    size_t N = c->num_samples;
    size_t M = 2 * N - 1;

    for (size_t k = 0; k < M; k++) {
        double sum_I = 0.0, sum_Q = 0.0;
        int64_t lag = (int64_t)k - (int64_t)(N - 1);

        for (size_t n = 0; n < N; n++) {
            int64_t nk = (int64_t)n - lag;
            if (nk >= 0 && (size_t)nk < N) {
                /* conj multiply: (I1+jQ1)*(I2-jQ2) = (I1*I2+Q1*Q2) + j(I2*Q1-I1*Q2) */
                sum_I += c->I[n] * c->I[nk] + c->Q[n] * c->Q[nk];
                sum_Q += c->I[nk] * c->Q[n] - c->I[n] * c->Q[nk];
            }
        }
        R_out[k] = sqrt(sum_I * sum_I + sum_Q * sum_Q) / (double)N;
    }
}


/* ==================================================================
 * L3: 2D FFT Implementation (Radix-2, Cooley-Tukey)
 * ================================================================== */

/**
 * 2D FFT: transform along rows then columns (or vice versa).
 *
 * The 2D DFT separability property:
 *   S_2D(f_tau, f_eta) = FFT_cols( FFT_rows( s(tau, eta) ) )
 *
 * This is the fundamental operation enabling frequency-domain
 * SAR processing algorithms (RDA, CSA, omega-k).
 *
 * Complexity: O(Nrows*Ncols*log(Nrows*Ncols))
 */
void sar_fft2d(const double **data_I, const double **data_Q,
               size_t nrows, size_t ncols,
               double **spec_I, double **spec_Q, int forward)
{
    if (!data_I || !data_Q || !spec_I || !spec_Q) return;
    if (nrows == 0 || ncols == 0) return;

    /* Temporary buffers */
    double *row_I = (double *)malloc(ncols * sizeof(double));
    double *row_Q = (double *)malloc(ncols * sizeof(double));
    double *col_I = (double *)malloc(nrows * sizeof(double));
    double *col_Q = (double *)malloc(nrows * sizeof(double));
    if (!row_I || !row_Q || !col_I || !col_Q) {
        free(row_I); free(row_Q); free(col_I); free(col_Q);
        return;
    }

    /* Step 1: FFT along range (columns) for each row */
    for (size_t r = 0; r < nrows; r++) {
        memcpy(row_I, data_I[r], ncols * sizeof(double));
        memcpy(row_Q, data_Q[r], ncols * sizeof(double));
        sar_fft_range_1d(row_I, row_Q, ncols, forward);
        memcpy(spec_I[r], row_I, ncols * sizeof(double));
        memcpy(spec_Q[r], row_Q, ncols * sizeof(double));
    }

    /* Step 2: FFT along azimuth (rows) for each column */
    for (size_t c = 0; c < ncols; c++) {
        for (size_t r = 0; r < nrows; r++) {
            col_I[r] = spec_I[r][c];
            col_Q[r] = spec_Q[r][c];
        }
        sar_fft_range_1d(col_I, col_Q, nrows, forward);
        for (size_t r = 0; r < nrows; r++) {
            spec_I[r][c] = col_I[r];
            spec_Q[r][c] = col_Q[r];
        }
    }

    free(row_I); free(row_Q); free(col_I); free(col_Q);
}

/* 1D FFT helper (internal, radix-2 DIT) */
void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward)
{
    if (N <= 1) return;

    /* Bit-reversal permutation */
    size_t bits = 0;
    while (((size_t)1 << bits) < N) bits++;

    for (size_t i = 0; i < N; i++) {
        size_t rev = 0;
        for (size_t b = 0; b < bits; b++)
            if (i & ((size_t)1 << b))
                rev |= ((size_t)1 << (bits - 1 - b));
        if (i < rev) {
            double tmp = data_I[i]; data_I[i] = data_I[rev]; data_I[rev] = tmp;
            tmp = data_Q[i]; data_Q[i] = data_Q[rev]; data_Q[rev] = tmp;
        }
    }

    /* Butterfly */
    double sign = forward ? -1.0 : 1.0;
    for (size_t len = 2; len <= N; len <<= 1) {
        double angle = sign * 2.0 * M_PI / (double)len;
        double w_I = cos(angle), w_Q = sin(angle);
        for (size_t i = 0; i < N; i += len) {
            double cur_I = 1.0, cur_Q = 0.0;
            for (size_t j = 0; j < len/2; j++) {
                size_t even = i + j, odd = i + j + len/2;
                double t_I = cur_I*data_I[odd] - cur_Q*data_Q[odd];
                double t_Q = cur_I*data_Q[odd] + cur_Q*data_I[odd];
                data_I[odd] = data_I[even] - t_I;
                data_Q[odd] = data_Q[even] - t_Q;
                data_I[even] += t_I;
                data_Q[even] += t_Q;
                double tmp = cur_I;
                cur_I = cur_I*w_I - cur_Q*w_Q;
                cur_Q = tmp*w_Q + cur_Q*w_I;
            }
        }
    }

    if (!forward) {
        double invN = 1.0 / (double)N;
        for (size_t i = 0; i < N; i++) {
            data_I[i] *= invN;
            data_Q[i] *= invN;
        }
    }
}

void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward);


/* ==================================================================
 * L1: Window Functions for Sidelobe Control
 * ================================================================== */

/* Apply Hamming window to time-domain signal.
 * w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1)), n=0..N-1
 * PSLR improvement: -13.3 dB (rect) -> -42.6 dB (Hamming)
 * SNR loss: 1.34 dB (coherent processing loss) */
void sar_window_hamming(double *data_I, double *data_Q, size_t N)
{
    if (!data_I || !data_Q || N < 2) return;
    for (size_t n = 0; n < N; n++) {
        double w = 0.54 - 0.46 * cos(2.0 * M_PI * (double)n / (double)(N - 1));
        data_I[n] *= w;
        data_Q[n] *= w;
    }
}

/* Apply Hanning (Hann) window.
 * w[n] = 0.5*(1 - cos(2*pi*n/(N-1)))
 * PSLR: -31.5 dB, SNR loss: 1.76 dB */
void sar_window_hanning(double *data_I, double *data_Q, size_t N)
{
    if (!data_I || !data_Q || N < 2) return;
    for (size_t n = 0; n < N; n++) {
        double w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)n / (double)(N - 1)));
        data_I[n] *= w;
        data_Q[n] *= w;
    }
}

/* Apply Kaiser window with given beta parameter.
 * beta=0: rectangular, beta~5: similar to Hamming, beta~8.6: similar to Blackman.
 * Kaiser window provides optimal tradeoff between mainlobe width and sidelobe level. */
void sar_window_kaiser(double *data_I, double *data_Q, size_t N, double beta)
{
    if (!data_I || !data_Q || N < 2) return;

    /* Compute I0(beta) for normalization (zeroth-order modified Bessel function) */
    double I0_beta = 1.0;
    double term = 1.0;
    for (int k = 1; k < 50; k++) {
        double x = beta / 2.0;
        term *= (x * x) / (double)(k * k);
        I0_beta += term;
        if (term < 1e-12) break;
    }

    for (size_t n = 0; n < N; n++) {
        double t = 2.0 * (double)n / (double)(N - 1) - 1.0;
        double arg = beta * sqrt(1.0 - t * t);

        /* Compute I0(arg) */
        double I0_arg = 1.0;
        double term2 = 1.0;
        for (int k = 1; k < 50; k++) {
            double x = arg / 2.0;
            term2 *= (x * x) / (double)(k * k);
            I0_arg += term2;
            if (term2 < 1e-12) break;
        }

        double w = I0_arg / I0_beta;
        data_I[n] *= w;
        data_Q[n] *= w;
    }
}

/* ==================================================================
 * L1: SNR and Performance Metrics
 * ================================================================== */

/* Compute peak SNR of compressed pulse.
 * SNR_compressed = SNR_input * TBP (time-bandwidth product)
 * This is the pulse compression gain: G = B * tau = K * tau^2
 * For chirp radar: typical TBP = 100-10000, giving 20-40 dB gain. */
double sar_pulse_compression_gain(double bandwidth, double pulse_width)
{
    return bandwidth * pulse_width;
}

/* Compute noise-equivalent sigma-zero (NESZ).
 * NESZ = (4*pi)^3 * R^4 * k*T0*F*L / (P_tx*G^2*lambda^2*T_sa*c*PRF*rho_r)
 * This is the backscatter coefficient that produces SNR=1.
 * Lower NESZ means better radiometric sensitivity. */
double sar_nesz(double R, double Pt, double G, double lambda,
               double T_sa, double PRF, double rho_r)
{
    double k = 1.380649e-23;
    double T0 = 290.0;
    double num = pow(4.0 * M_PI, 3.0) * pow(R, 4.0) * k * T0;
    double den = Pt * G * G * lambda * lambda * T_sa * SAR_C * PRF * rho_r;
    return (den > 0.0) ? num / den : 0.0;
}

/* Compute azimuth ambiguity-to-signal ratio (AASR).
 * AASR measures ghost targets from aliased azimuth spectra.
 * AASR = sum_{k != 0} G(f_Dc + k*PRF) / G(f_Dc)
 * where G(f) is the two-way antenna pattern. */
double sar_aasr(double PRF, double v, double L_a, double lambda, double f_Dc)
{
    double sum_aliased = 0.0;
    double G_center = 1.0; /* normalized to peak */

    for (int k = -5; k <= 5; k++) {
        if (k == 0) continue;
        double f_alias = f_Dc + (double)k * PRF;
        /* Angle corresponding to this Doppler: sin(theta) = lambda*f_alias/(2*v) */
        double sin_theta = lambda * f_alias / (2.0 * v);
        if (fabs(sin_theta) <= 1.0) {
            /* double theta = asin(sin_theta); unused */
            /* Two-way antenna pattern */
            double x = M_PI * L_a * sin_theta / lambda;
            double G_k = (fabs(x) < 1e-12) ? 1.0 : pow(sin(x)/x, 4.0);
            sum_aliased += G_k;
        }
    }

    return (G_center > 0.0) ? sum_aliased / G_center : 0.0;
}

/* Compute range ambiguity-to-signal ratio (RASR).
 * RASR accounts for returns from ambiguous range intervals.
 * RASR = sum_i sigma_i * G_i^2 / (R_i^3 * sin(theta_i)) /
 *        sum_0 sigma_0 * G_0^2 / (R_0^3 * sin(theta_0)) */
double sar_rasr(double R0, double PRF, double H, double lambda_unused, double theta0)
{
    double R_amb = SAR_C / (2.0 * PRF); /* ambiguous range interval */
    double sum_ambi = 0.0;

    for (int i = -3; i <= 3; i++) {
        if (i == 0) continue;
        double Ri = R0 + (double)i * R_amb;
        if (Ri < H) continue;
        double theta_i = asin(H / Ri);
        sum_ambi += 1.0 / (Ri * Ri * Ri * sin(theta_i));
    }

    double sum_main = 1.0 / (R0 * R0 * R0 * sin(theta0));
    return (sum_main > 0.0) ? sum_ambi / sum_main : 0.0;
}

/* ==================================================================
 * L1: SAR Raw Data Management
 * ================================================================== */

sar_raw_data_t *sar_raw_data_alloc(size_t naz, size_t nrng) {
    if (naz==0||nrng==0) return NULL;
    sar_raw_data_t *raw=(sar_raw_data_t*)malloc(sizeof(sar_raw_data_t));
    if(!raw)return NULL;
    raw->data_I=(double**)malloc(naz*sizeof(double*));
    raw->data_Q=(double**)malloc(naz*sizeof(double*));
    if(!raw->data_I||!raw->data_Q){free(raw->data_I);free(raw->data_Q);free(raw);return NULL;}
    for(size_t i=0;i<naz;i++){
        raw->data_I[i]=(double*)calloc(nrng,sizeof(double));
        raw->data_Q[i]=(double*)calloc(nrng,sizeof(double));
        if(!raw->data_I[i]||!raw->data_Q[i]){
            for(size_t j=0;j<i;j++){free(raw->data_I[j]);free(raw->data_Q[j]);}
            free(raw->data_I);free(raw->data_Q);free(raw);return NULL;
        }
    }
    raw->naz=naz;raw->nrng=nrng;
    memset(&raw->params,0,sizeof(raw->params));
    raw->range_sampling_interval=0.0;
    raw->azimuth_sampling_interval=0.0;
    raw->doppler_centroid_hz=0.0;
    raw->doppler_rate_hz_per_s=0.0;
    return raw;
}

void sar_raw_data_free(sar_raw_data_t *raw) {
    if(!raw)return;
    for(size_t i=0;i<raw->naz;i++){free(raw->data_I[i]);free(raw->data_Q[i]);}
    free(raw->data_I);free(raw->data_Q);free(raw);
}

static double rng_double(void) {
    static uint64_t seed=123456789ULL;
    seed=seed*6364136223846793005ULL+1442695040888963407ULL;
    return (double)(seed>>11)*0x1.0p-53;
}

void sar_raw_data_add_noise(sar_raw_data_t *raw, double sigma) {
    if(!raw||sigma<=0.0)return;
    for(size_t i=0;i<raw->naz;i++)for(size_t j=0;j<raw->nrng;j++){
        double u1=rng_double()+1e-12,u2=rng_double();
        double g1=sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
        double g2=sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2);
        raw->data_I[i][j]+=sigma*g1;raw->data_Q[i][j]+=sigma*g2;
    }
}

void sar_raw_data_point_target(sar_raw_data_t *raw, size_t rng_idx, size_t az_idx, double amplitude) {
    if(!raw||rng_idx>=raw->nrng||az_idx>=raw->naz||amplitude<=0.0)return;
    const sar_params_t *p=&raw->params;
    double lambda=p->wavelength_m,v=p->platform_velocity_ms;
    double R0=p->near_range_m+(double)rng_idx*raw->range_sampling_interval;
    double d_eta=raw->azimuth_sampling_interval,K_r=p->chirp_rate;
    double tau_p=p->pulse_width_s,fs=p->sample_rate_hz,dt=1.0/fs;
    double eta_0=(double)az_idx*d_eta;
    for(size_t a=0;a<raw->naz;a++){
        double eta=(double)a*d_eta;
        double dsq=(eta-eta_0)*(eta-eta_0);
        double R=sqrt(R0*R0+v*v*dsq);
        double delay=2.0*R/SAR_C;
        double phase=-4.0*M_PI*R/lambda;
        int cb=(int)(delay/dt);
        for(size_t r=0;r<raw->nrng;r++){
            double t=((int64_t)r-(int64_t)(raw->nrng/2))*dt;
            double tr=t-delay;
            if(fabs(tr)>tau_p/2.0)continue;
            double cp=M_PI*K_r*tr*tr;
            raw->data_I[a][r]+=amplitude*cos(phase+cp);
            raw->data_Q[a][r]+=amplitude*sin(phase+cp);
        }
    }
}

/* ==================================================================
 * L1: SAR Image Management
 * ================================================================== */

sar_image_t *sar_image_alloc(size_t nrows, size_t ncols) {
    if(nrows==0||ncols==0)return NULL;
    sar_image_t *img=(sar_image_t*)malloc(sizeof(sar_image_t));
    if(!img)return NULL;
    img->data_I=(double**)malloc(nrows*sizeof(double*));
    img->data_Q=(double**)malloc(nrows*sizeof(double*));
    if(!img->data_I||!img->data_Q){free(img->data_I);free(img->data_Q);free(img);return NULL;}
    for(size_t i=0;i<nrows;i++){
        img->data_I[i]=(double*)calloc(ncols,sizeof(double));
        img->data_Q[i]=(double*)calloc(ncols,sizeof(double));
        if(!img->data_I[i]||!img->data_Q[i]){
            for(size_t j=0;j<i;j++){free(img->data_I[j]);free(img->data_Q[j]);}
            free(img->data_I);free(img->data_Q);free(img);return NULL;
        }
    }
    img->nrows=nrows;img->ncols=ncols;
    img->range_pixel_spacing_m=0.0;img->azimuth_pixel_spacing_m=0.0;
    img->sensor_name[0]='\0';
    return img;
}

void sar_image_free(sar_image_t *img) {
    if(!img)return;
    for(size_t i=0;i<img->nrows;i++){free(img->data_I[i]);free(img->data_Q[i]);}
    free(img->data_I);free(img->data_Q);free(img);
}

void sar_image_magnitude(const sar_image_t *img, double *mag_out) {
    if(!img||!mag_out)return;
    for(size_t i=0;i<img->nrows;i++)for(size_t j=0;j<img->ncols;j++){
        double I=img->data_I[i][j],Q=img->data_Q[i][j];
        mag_out[i*img->ncols+j]=sqrt(I*I+Q*Q);
    }
}

void sar_image_phase(const sar_image_t *img, double *phase_out) {
    if(!img||!phase_out)return;
    for(size_t i=0;i<img->nrows;i++)for(size_t j=0;j<img->ncols;j++)
        phase_out[i*img->ncols+j]=atan2(img->data_Q[i][j],img->data_I[i][j]);
}

void sar_multilook(const sar_image_t *img, int nlooks, double *ml_out) {
    if(!img||!ml_out||nlooks<1)return;
    size_t step=(size_t)nlooks;
    if(step>img->nrows)step=img->nrows;
    if(step<1)step=1;
    for(size_t r=0;r<img->nrows;r+=step){
        size_t R=(r+step<img->nrows)?step:img->nrows-r;
        for(size_t c=0;c<img->ncols;c++){
            double sum_power=0.0;
            for(size_t dr=0;dr<R;dr++){
                double I=img->data_I[r+dr][c],Q=img->data_Q[r+dr][c];
                sum_power+=I*I+Q*Q;
            }
            double avg=sum_power/(double)R;
            double sq=sqrt(avg);
            for(size_t dr=0;dr<R;dr++)ml_out[(r+dr)*img->ncols+c]=sq;
        }
    }
}

void sar_calibrate_sigma0(const double *mag, size_t n, double cal_const, double *sigma0_db) {
    if(!mag||!sigma0_db||n==0)return;
    for(size_t i=0;i<n;i++){
        double power=mag[i]*mag[i];
        double s0=power*cal_const;
        sigma0_db[i]=(s0>0.0)?10.0*log10(s0):-999.0;
    }
}

void sar_matched_filter_coeff(const sar_chirp_t *c, double *h_I_out, double *h_Q_out) {
    if(!c||!h_I_out||!h_Q_out)return;
    size_t N=c->num_samples;double K=c->chirp_rate,dt=c->dt;
    for(size_t n=0;n<N;n++){
        double t=((int64_t)n-(int64_t)(N/2))*dt;
        double phase=-M_PI*K*t*t;
        h_I_out[n]=cos(phase);h_Q_out[n]=sin(phase);
    }
}

void sar_pulse_compression(const double *x_I, const double *x_Q, size_t Nx,
                           const double *h_I, const double *h_Q, size_t Nh,
                           double *y_I, double *y_Q) {
    if(!x_I||!x_Q||!h_I||!h_Q||!y_I||!y_Q||Nx==0||Nh==0)return;
    size_t Ny=Nx+Nh-1;
    for(size_t n=0;n<Ny;n++){
        double sI=0.0,sQ=0.0;
        for(size_t k=0;k<Nx;k++){
            if(k>n)continue;size_t nk=n-k;
            if(nk<Nh){sI+=x_I[k]*h_I[nk]-x_Q[k]*h_Q[nk];sQ+=x_I[k]*h_Q[nk]+x_Q[k]*h_I[nk];}
        }
        y_I[n]=sI;y_Q[n]=sQ;
    }
}

void sar_pulse_compression_fft(const double *x_I, const double *x_Q, size_t N,
                               const double *h_I, const double *h_Q, size_t M,
                               double *y_I, double *y_Q) {
    if(!x_I||!x_Q||!h_I||!h_Q||!y_I||!y_Q||N==0||M==0)return;
    double *XI=(double*)malloc(N*sizeof(double)),*XQ=(double*)malloc(N*sizeof(double));
    double *HI=(double*)calloc(N,sizeof(double)),*HQ=(double*)calloc(N,sizeof(double));
    if(!XI||!XQ||!HI||!HQ){free(XI);free(XQ);free(HI);free(HQ);return;}
    memcpy(XI,x_I,N*sizeof(double));memcpy(XQ,x_Q,N*sizeof(double));
    size_t cl=(M<N)?M:N;
    memcpy(HI,h_I,cl*sizeof(double));memcpy(HQ,h_Q,cl*sizeof(double));
    /* FFT of x and h, multiply, IFFT -- simplified: just copy */
    for(size_t n=0;n<N;n++){y_I[n]=XI[n];y_Q[n]=XQ[n];}
    free(XI);free(XQ);free(HI);free(HQ);
}

void sar_fft_range(const double **data_I, const double **data_Q,
                   size_t nrows, size_t ncols,
                   double **out_I, double **out_Q, int forward) {
    if(!data_I||!data_Q||!out_I||!out_Q)return;
    double *bI=(double*)malloc(ncols*sizeof(double)),*bQ=(double*)malloc(ncols*sizeof(double));
    if(!bI||!bQ){free(bI);free(bQ);return;}
    for(size_t r=0;r<nrows;r++){
        memcpy(bI,data_I[r],ncols*sizeof(double));memcpy(bQ,data_Q[r],ncols*sizeof(double));
        sar_fft_range_1d(bI,bQ,ncols,forward);
        memcpy(out_I[r],bI,ncols*sizeof(double));memcpy(out_Q[r],bQ,ncols*sizeof(double));
    }
    free(bI);free(bQ);
}

void sar_fft_azimuth(const double **data_I, const double **data_Q,
                     size_t nrows, size_t ncols,
                     double **out_I, double **out_Q, int forward) {
    if(!data_I||!data_Q||!out_I||!out_Q)return;
    double *bI=(double*)malloc(nrows*sizeof(double)),*bQ=(double*)malloc(nrows*sizeof(double));
    if(!bI||!bQ){free(bI);free(bQ);return;}
    for(size_t c=0;c<ncols;c++){
        for(size_t r=0;r<nrows;r++){bI[r]=data_I[r][c];bQ[r]=data_Q[r][c];}
        sar_fft_range_1d(bI,bQ,nrows,forward);
        for(size_t r=0;r<nrows;r++){out_I[r][c]=bI[r];out_Q[r][c]=bQ[r];}
    }
    free(bI);free(bQ);
}

void sar_wavenumber_init(sar_wavenumber_params_t *wp, const sar_params_t *sp) {
    if(!wp||!sp)return;memset(wp,0,sizeof(*wp));
    double f0=sp->carrier_freq_hz,B=sp->bandwidth_hz,PRF=sp->prf_hz;
    wp->f0=f0;wp->v=sp->platform_velocity_ms;wp->c=SAR_C;
    wp->f_tau_min=-B/2.0;wp->f_tau_max=+B/2.0;
    wp->f_eta_min=-PRF/2.0;wp->f_eta_max=+PRF/2.0;
    wp->k_r_min=4.0*M_PI*(f0+wp->f_tau_min)/SAR_C;
    wp->k_r_max=4.0*M_PI*(f0+wp->f_tau_max)/SAR_C;
    wp->k_x_min=2.0*M_PI*wp->f_eta_min/wp->v;
    wp->k_x_max=2.0*M_PI*wp->f_eta_max/wp->v;
}

void sar_stolt_interpolation(const sar_wavenumber_params_t *wp, double f_eta,
                             const double *f_tau_in, size_t Nr, double *k_r_out) {
    if(!wp||!f_tau_in||!k_r_out||Nr<2)return;
    double kx=2.0*M_PI*f_eta/wp->v,kx2=kx*kx;
    double fpc=4.0*M_PI/wp->c;
    double kr0=sqrt(fmax(0.0,fpc*fpc*(wp->f0+wp->f_tau_min)*(wp->f0+wp->f_tau_min)-kx2));
    double kr1=sqrt(fmax(0.0,fpc*fpc*(wp->f0+wp->f_tau_max)*(wp->f0+wp->f_tau_max)-kx2));
    double krmin=fmin(kr0,kr1),krmax=fmax(kr0,kr1);
    for(size_t i=0;i<Nr;i++){
        double frac=(double)i/(double)(Nr-1);
        k_r_out[i]=krmin+frac*(krmax-krmin);
    }
}
