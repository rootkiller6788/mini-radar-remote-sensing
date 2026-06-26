/** @file test_ir_detector.c */
#include <stdio.h>
#include <math.h>
#include "ir_core.h"
#include "ir_detector.h"

static int run=0,pass=0;
#define T(n) do{run++;printf("  TEST %s... ",n);}while(0)
#define P() do{pass++;printf("PASS\n");}while(0)
#define C(c) do{if(!(c)){printf("FAIL\n");return 1;}}while(0)

int main(void) {
  printf("=== Test: ir_detector ===\n");
  ir_detector_params_t p;
  T("init MCT"); C(ir_detector_params_init(&p,IR_DET_PHOTON_MCT)==0); P();
  T("init bolo"); C(ir_detector_params_init(&p,IR_DET_THERMAL_BOLOMETER)==0); P();
  T("responsivity"); C(ir_detector_responsivity_current(0.7,10)>0); P();
  T("Johnson"); C(ir_noise_johnson_voltage(300,1e6,1000)>0); P();
  T("shot"); C(ir_noise_shot_current(1e-6,1000)>0); P();
  T("NEP"); C(ir_nep(&p,1000,60,10)>0); P();
  T("D*"); C(ir_specific_detectivity(&p,1e-12,1000)>0); P();
  T("BLIP"); C(ir_d_star_blip(0.7,10,1e16)>0); P();
  printf("  %d/%d tests passed\n\n",pass,run);
  return (pass==run)?0:1;
}
