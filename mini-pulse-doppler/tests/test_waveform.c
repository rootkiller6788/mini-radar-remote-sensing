#include "../include/radar_waveform.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int test_failed_wf = 0;
#define ASSERT_DOUBLE_EQ(v, e, tol) do { \
    double _v=(v); if(fabs(_v-(e))>(tol)){ \
    printf("FAIL: got %g exp %g\n",_v,(e));test_failed_wf=1;} \
} while(0)

static int tr=0, tp=0;
#define T(n) do{tr++;printf("  %s... ",n);}while(0)
#define P() do{if(!test_failed_wf){tp++;printf("PASS\n");}test_failed_wf=0;}while(0)

int main(void) {
    printf("=== waveform tests ===\n");
    radar_waveform_params_t wf = {.pulse_width=10e-6, .bandwidth=10e6,
        .center_frequency=10e9, .prf=10000, .pri=1e-4,
        .duty_cycle=0.1, .sampling_rate=50e6, .num_samples=1024,
        .type=WAVEFORM_LFM_UP};

    double complex *buf = malloc(1024*sizeof(double complex));
    assert(buf);

    T("rect_pulse_generate");
    wf.type = WAVEFORM_RECT_PULSE;
    assert(rect_pulse_generate(&wf, buf) == 0);
    int nonzero = 0;
    for(int i=0;i<1024;i++) if(cabs(buf[i])>0.5) nonzero++;
    assert(nonzero > 0 && nonzero < 1024);
    P();

    T("lfm_chirp_generate");
    lfm_params_t chirp = {.chirp_rate=1e12, .time_bw_product=100.0, .is_upchirp=1};
    wf.type = WAVEFORM_LFM_UP;
    assert(lfm_chirp_generate(&wf, &chirp, buf) == 0);
    double max_mag = 0;
    for(int i=0;i<1024;i++){double m=cabs(buf[i]);if(m>max_mag)max_mag=m;}
    assert(max_mag > 0.9);
    P();

    T("barker_code_generate");
    barker_code_t bc = {.code_length=13, .chip_rate=10e6,
        .chip_duration=1e-7, .code_sequence=(int8_t*)(int8_t[]){1,1,1,1,1,-1,-1,1,1,-1,1,-1,1}};
    assert(barker_code_generate(&wf, &bc, buf) == 0);
    P();

    T("waveform_autocorrelation");
    double acorr[1024];
    for(int i=0;i<512;i++) buf[i]=1.0+0.0*I;
    for(int i=512;i<1024;i++) buf[i]=0.0+0.0*I;
    assert(waveform_autocorrelation(buf, 1024, acorr) == 0);
    assert(acorr[0] > 500.0);
    assert(acorr[100] < acorr[0]);
    P();

    T("waveform_rms_bandwidth");
    double bw;
    assert(waveform_rms_bandwidth(buf, 1024, 50e6, &bw) == 0);
    assert(bw > 0.0);
    P();

    T("waveform_apply_window");
    double complex *buf2 = malloc(128*sizeof(double complex));
    for(int i=0;i<128;i++) buf2[i]=1.0+0.0*I;
    assert(waveform_apply_window(buf2, 128, 1, 0.0, 0.0) == 0);
    ASSERT_DOUBLE_EQ(creal(buf2[0]), 0.08, 0.02);
    ASSERT_DOUBLE_EQ(creal(buf2[64]), 1.0, 0.02);
    free(buf2);
    P();

    T("waveform_cross_ambiguity");
    double result;
    double complex *tx = malloc(64*sizeof(double complex));
    for(int i=0;i<64;i++) tx[i]=1.0+0.0*I;
    assert(waveform_cross_ambiguity_1d(tx,tx,64,0.0,&result)==0);
    assert(result > 1000.0);
    free(tx);
    P();

    T("instantaneous_freq");
    double complex *chirp_buf = malloc(256*sizeof(double complex));
    for(int i=0;i<256;i++) {
        double t=((double)i-127.5)/50e6;
        chirp_buf[i]=cos(M_PI*1e12*t*t)+sin(M_PI*1e12*t*t)*I;
    }
    double freq[256];
    assert(waveform_instantaneous_freq(chirp_buf,256,freq,50e6)==0);
    free(chirp_buf);
    P();

    free(buf);
    printf("\n=== %d/%d passed ===\n", tp, tr);
    return (tp==tr)?0:1;
}
