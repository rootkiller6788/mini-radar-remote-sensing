#include "../include/doppler_processing.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tr=0,tp=0;
static int test_failed = 0;
#define T(n) do{tr++;printf("  %s... ",n);}while(0)
#define P() do{if(!test_failed){tp++;printf("PASS\n");}test_failed=0;}while(0)
#define AD(v,e,t) do{double _v=(v);if(fabs(_v-(e))>(t)){ \
    printf("FAIL: g=%g e=%g\n",_v,(e));test_failed=1;} \
}while(0)

int main(void){
    printf("=== doppler tests ===\n");

    T("window_function_create");
    window_func_t win;
    assert(window_function_create(WINDOW_HAMMING,64,60.0,&win)==0);
    assert(win.length==64);
    assert(fabs(win.coherent_gain-0.54)<0.1);
    window_function_free(&win);
    P();

    T("velocity_from_doppler");
    double v;
    assert(velocity_from_doppler(1000.0,0.03,&v)==0);
    AD(v,15.0,0.1);
    P();

    T("doppler_from_velocity");
    double fd;
    assert(doppler_from_velocity(15.0,0.03,&fd)==0);
    AD(fd,1000.0,1.0);
    P();

    T("max_unambiguous_velocity");
    double vmax;
    assert(max_unambiguous_velocity(5000.0,0.03,&vmax)==0);
    AD(vmax,37.5,0.1);
    P();

    T("blind_speed_compute");
    double blind;
    assert(blind_speed_compute(5000.0,0.03,1,&blind)==0);
    AD(blind,75.0,0.1);
    P();

    T("stagger_prf_design");
    staggered_prf_t stag;
    assert(stagger_prf_design(0.03,50.0,3,&stag)==0);
    assert(stag.num_prfs==3);
    stagger_prf_free(&stag);
    P();

    T("doppler_spectrum_compute");
    double complex *slow = malloc(64*sizeof(double complex));
    for(int i=0;i<64;i++){
        double phase=2.0*M_PI*0.1*(double)i;
        slow[i]=cos(phase)+sin(phase)*I;
    }
    doppler_spectrum_t spec;
    assert(doppler_spectrum_compute(slow,64,5000.0,0.03,WINDOW_HAMMING,&spec)==0);
    assert(spec.num_bins==64);
    doppler_spectrum_free(&spec);
    free(slow);
    P();

    T("mti_filter_design");
    mti_filter_t mti;
    assert(mti_filter_design(MTI_DOUBLE_DELAY,5000.0,10.0,&mti)==0);
    assert(mti.filter_order==3);
    mti_filter_free(&mti);
    P();

    printf("\n=== %d/%d passed ===\n",tp,tr);
    return (tp==tr)?0:1;
}
