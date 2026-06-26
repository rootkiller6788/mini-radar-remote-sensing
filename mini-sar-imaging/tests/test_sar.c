#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "sar_core.h"
#include "sar_geometry.h"
#include "sar_algorithm.h"
#include "sar_interferometry.h"
#include "sar_advanced.h"

#define TOL 1e-9
#define TP() printf("  PASS: %s\n", __func__)
static int tr=0, tp=0;
void sar_fft_range_1d(double *data_I, double *data_Q, size_t N, int forward);

static int t_params(void) {
    sar_params_t p;
    sar_params_init(&p, 5.3e9, 100e6, 40e-6, 1700.0, 4.8, 7500.0, 514e3, 0.0, 0.35, 830e3, 870e3);
    assert(fabs(p.carrier_freq_hz-5.3e9)<TOL);
    assert(p.range_resolution_m>0.0);
    assert(fabs(p.azimuth_resolution_m-2.4)<0.1);
    TP(); return 0;
}
static int t_chirp(void) {
    sar_chirp_t *c=sar_chirp_alloc(512,40e-6,100e6,150e6);
    assert(c!=NULL); assert(c->time_bandwidth_product>0.0);
    double sI=0,sQ=0;
    for(size_t i=0;i<c->num_samples;i++){sI+=fabs(c->I[i]);sQ+=fabs(c->Q[i]);}
    assert(sI>0.0&&sQ>0.0);
    sar_chirp_free(c); TP(); return 0;
}
static int t_mf(void) {
    sar_chirp_t *c=sar_chirp_alloc(128,10e-6,50e6,100e6);
    double *hI=malloc(128*8),*hQ=malloc(128*8);
    sar_matched_filter_coeff(c,hI,hQ);
    double s=0; for(int i=0;i<128;i++)s+=hI[i]*hI[i]+hQ[i]*hQ[i];
    assert(s>0.0); free(hI);free(hQ);sar_chirp_free(c);TP();return 0;
}
static int t_pc(void) {
    size_t N=128,M=64;
    double *xI=calloc(N,8),*xQ=calloc(N,8),*hI=calloc(M,8),*hQ=calloc(M,8);
    double *yI=calloc(N+M-1,8),*yQ=calloc(N+M-1,8);
    xI[N/2]=1.0;hI[M/2]=1.0;
    sar_pulse_compression(xI,xQ,N,hI,hQ,M,yI,yQ);
    double mx=0; for(size_t i=0;i<N+M-1;i++){double m=yI[i]*yI[i]+yQ[i]*yQ[i];if(m>mx)mx=m;}
    assert(mx>0.9); free(xI);free(xQ);free(hI);free(hQ);free(yI);free(yQ);TP();return 0;
}
static int t_ac(void) {
    sar_chirp_t *c=sar_chirp_alloc(200,20e-6,30e6,60e6);
    size_t M=2*c->num_samples-1; double *R=malloc(M*8);
    sar_chirp_autocorrelation(c,R);
    double pk=0;size_t pi=0;
    for(size_t i=0;i<M;i++){if(R[i]>pk){pk=R[i];pi=i;}}
    assert(fabs(pk-1.0)<0.1); free(R);sar_chirp_free(c);TP();return 0;
}
static int t_raw(void) {
    sar_raw_data_t *raw=sar_raw_data_alloc(64,128);assert(raw!=NULL);
    sar_raw_data_free(raw); TP();return 0;
}
static int t_pt(void) {
    sar_raw_data_t *raw=sar_raw_data_alloc(64,128);
    sar_params_init(&raw->params,5.3e9,100e6,40e-6,1700.0,4.8,7500.0,0.0,0.0,0.35,0.0,1000.0);
    raw->range_sampling_interval=SAR_C/(2.0*120e6);
    raw->azimuth_sampling_interval=1.0/1700.0;
    sar_raw_data_point_target(raw,50,30,1.0);
    double s=0;for(size_t i=0;i<raw->naz;i++)for(size_t j=0;j<raw->nrng;j++)s+=fabs(raw->data_I[i][j]);
    assert(s>0.0); sar_raw_data_free(raw);TP();return 0;
}
static int t_img(void) {
    sar_image_t *img=sar_image_alloc(32,32);assert(img!=NULL);
    for(size_t i=0;i<img->nrows;i++)for(size_t j=0;j<img->ncols;j++){img->data_I[i][j]=1.0;img->data_Q[i][j]=0.0;}
    double *mag=malloc(1024*8); sar_image_magnitude(img,mag);
    for(size_t i=0;i<1024;i++)assert(mag[i]>=0.0);
    free(mag);sar_image_free(img);TP();return 0;
}
static int t_ml(void) {
    sar_image_t *img=sar_image_alloc(32,16);
    for(size_t i=0;i<img->nrows;i++)for(size_t j=0;j<img->ncols;j++){img->data_I[i][j]=1.0;img->data_Q[i][j]=0.0;}
    double *ml=malloc(512*8); sar_multilook(img,4,ml);
    double sum=0;for(size_t i=0;i<512;i++)sum+=ml[i];
    assert(fabs(sum/512.0-1.0)<0.1); free(ml);sar_image_free(img);TP();return 0;
}
static int t_cal(void) {
    double mag[]={10.0,100.0,1000.0},s0[3];
    sar_calibrate_sigma0(mag,3,1.0,s0);
    assert(fabs(s0[0]-20.0)<0.1); TP();return 0;
}
static int t_range(void) {
    double R=sar_range_hyperbolic(1000.0,100.0,0.0,0.0);
    assert(fabs(R-1000.0)<TOL); TP();return 0;
}
static int t_dop(void) {
    double fdc=sar_doppler_centroid(0.0566,7500.0,0.0);assert(fabs(fdc)<TOL);
    double fR=sar_doppler_rate(0.0566,7500.0,850e3,0.0);assert(fR<0.0);
    TP();return 0;
}
static int t_rcm(void) {
    sar_rcm_t rcm;sar_rcm_init(&rcm,850e3,7500.0,0.0566,0.0,1.5);
    assert(rcm.max_rcm_m>=0.0); TP();return 0;
}
static int t_coord(void) {
    sar_geodetic_t geo={34.0,-118.0,500.0};sar_ecef_t ecef;
    sar_geodetic_to_ecef(&geo,&ecef);
    sar_geodetic_t geo2;sar_ecef_to_geodetic(&ecef,&geo2);
    assert(fabs(geo2.lat_deg-34.0)<0.01); TP();return 0;
}
static int t_ant(void) {
    sar_antenna_t ant;sar_antenna_init(&ant,4.8,0.0566,850e3);
    double g0=sar_antenna_gain(&ant,0.0);assert(fabs(g0-1.0)<0.01);
    TP();return 0;
}
static int t_ev(void) {
    double ve=sar_effective_velocity(7500.0,514e3,6371e3);
    assert(ve>0.0&&ve<7500.0); TP();return 0;
}
static int t_res(void) {
    sar_params_t p;sar_params_init(&p,5.3e9,100e6,40e-6,1700.0,4.8,7500.0,514e3,0.0,0.35,830e3,870e3);
    assert(sar_resolution_check(&p,0.01)==1); TP();return 0;
}
static int t_prf(void) {
    sar_params_t p;sar_params_init(&p,5.3e9,100e6,40e-6,1700.0,4.8,7500.0,514e3,0.0,0.35,830e3,870e3);
    int r=sar_prf_check(&p); assert(r==0||r==1); TP();return 0;
}
static int t_insar(void) {
    sar_insar_pair_t *pair=sar_insar_pair_alloc();
    sar_insar_set_baseline(pair,200.0,0.3,514e3,850e3,0.0566);
    assert(pair->baseline_m>0.0); sar_insar_pair_free(pair); TP();return 0;
}
static int t_coh(void) {
    sar_image_t *m=sar_image_alloc(16,16),*s=sar_image_alloc(16,16);
    for(size_t i=0;i<16;i++)for(size_t j=0;j<16;j++){m->data_I[i][j]=1.0;s->data_I[i][j]=1.0;}
    sar_coherence_map_t *cm=sar_coherence_alloc(16,16);
    sar_coherence_estimate(m,s,cm,3,3); assert(cm->mean_coherence>0.9);
    sar_coherence_free(cm);sar_image_free(m);sar_image_free(s);TP();return 0;
}
static int t_bis(void) {
    sar_bistatic_config_t bc;memset(&bc,0,sizeof(bc));
    bc.tx_x=0;bc.tx_z=500e3;bc.rx_x=1000;bc.rx_z=500e3;
    double R=sar_bistatic_range(&bc,0.0,500.0,0.0,0.0);
    assert(R>0.0); TP();return 0;
}
static int t_fft(void) {
    size_t N=64;double *I=calloc(N,8),*Q=calloc(N,8);I[0]=1.0;
    sar_fft_range_1d(I,Q,N,1);sar_fft_range_1d(I,Q,N,0);
    assert(fabs(I[0]-1.0)<1e-9);free(I);free(Q);TP();return 0;
}
static int t_fd(void) {
    sar_coherency_matrix_t T={1.0,0.5,0.25,0.3,0,0,0,0,0};
    double Ps,Pd,Pv;sar_freeman_durden(&T,1,&Ps,&Pd,&Pv);
    assert(Ps>=0&&Pd>=0&&Pv>=0); TP();return 0;
}
static int t_ha(void) {
    sar_coherency_matrix_t T={1.0,0.5,0.25,0,0,0,0,0,0};
    double H,a,A;sar_h_alpha_decomp(&T,1,&H,&a,&A);
    assert(H>=0&&H<=1); TP();return 0;
}
int main(void) {
    printf("SAR Imaging Test Suite\n====================\n");
    struct{const char*n;int(*f)(void);}tests[]={
        {"params",t_params},{"chirp",t_chirp},{"mf",t_mf},{"pc",t_pc},
        {"autocorr",t_ac},{"raw",t_raw},{"ptarget",t_pt},{"image",t_img},
        {"multilook",t_ml},{"calibrate",t_cal},{"range",t_range},
        {"doppler",t_dop},{"rcm",t_rcm},{"coord",t_coord},{"antenna",t_ant},
        {"effvel",t_ev},{"rescheck",t_res},{"prf",t_prf},
        {"insar",t_insar},{"coherence",t_coh},{"bistatic",t_bis},
        {"fft",t_fft},{"freeman",t_fd},{"halpha",t_ha},
    };
    int nt=sizeof(tests)/sizeof(tests[0]);
    for(int i=0;i<nt;i++){printf("Test %d: %s\n",i+1,tests[i].n);if(tests[i].f()==0)tp++;else break;tr++;}
    printf("\nResults: %d/%d passed\n",tp,nt);
    return(tp==nt)?0:1;
}