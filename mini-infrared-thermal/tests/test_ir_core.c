/** @file test_ir_core.c */
#include <stdio.h>
#include <math.h>
#include "ir_core.h"

static int run=0, pass=0;
#define T(n) do{run++;printf("  TEST %s... ",n);}while(0)
#define P() do{pass++;printf("PASS\n");}while(0)
#define C(c) do{if(!(c)){printf("FAIL at %d\n",__LINE__);return 1;}}while(0)

int main(void) {
  printf("=== Test: ir_core ===\n");
  T("Planck"); double B=ir_planck_spectral_radiance_wavelength(10e-6,300); C(B>1e6&&B<1e9); P();
  T("Wien"); C(ir_wien_peak_wavelength(300)>0); P();
  T("Stefan-Boltzmann"); double M=ir_stefan_boltzmann_exitance(300); C(M>400&&M<500); P();
  T("BB init"); ir_blackbody_t bb; C(ir_blackbody_init(&bb,300,0.95,0.01,1)==0);
    C(ir_blackbody_init(NULL,300,0.5,0.01,1)==-1); P();
  T("Brightness"); B=ir_planck_spectral_radiance_wavelength(10e-6,310);
    C(fabs(ir_brightness_temperature(10e-6,B)-310)<1.0); P();
  T("Contrast"); C(ir_thermal_contrast(310,300)>0); P();
  T("Invalid inputs"); C(ir_planck_spectral_radiance_wavelength(-1,300)<0); P();
  printf("  %d/%d tests passed\n\n",pass,run);
  return (pass==run)?0:1;
}
