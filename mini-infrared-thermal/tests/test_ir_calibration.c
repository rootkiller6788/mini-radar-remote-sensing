/** @file test_ir_calibration.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ir_core.h"
#include "ir_calibration.h"

#define RR 4
#define CC 4
#define NN (RR*CC)

static int run=0,pass=0;
#define T(n) do{run++;printf("  TEST %s... ",n);}while(0)
#define P() do{pass++;printf("PASS\n");}while(0)
#define C(c) do{if(!(c)){printf("FAIL\n");return 1;}}while(0)

int main(void) {
  printf("=== Test: ir_calibration ===\n");
  double lo[NN]={100,110,105,95,102,108,107,98,101,109,103,97,104,106,100,99};
  double hi[NN]={200,210,205,195,202,208,207,198,201,209,203,197,204,206,200,199};
  ir_cal_pixel_coeff_t coeffs[NN];
  T("2pt comp");
  C(ir_cal_two_point_compute(lo,hi,300,400,RR,CC,coeffs)==0); P();
  T("2pt apply");
  double corr[NN],tst[NN]={150,160,155,145,152,158,157,148,151,159,153,147,154,156,150,149};
  C(ir_cal_two_point_apply(coeffs,tst,RR,CC,corr)==0); P();
  T("res NUC"); C(ir_nuc_residual_nuc(corr,RR,CC)>=0); P();
  T("raw2temp"); C(ir_cal_raw_to_temperature(150,1,100,10)>0); P();
  printf("  %d/%d tests passed\n\n",pass,run);
  return (pass==run)?0:1;
}
