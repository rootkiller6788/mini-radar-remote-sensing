/** @file test_ir_radiometry.c */
#include <stdio.h>
#include <math.h>
#include "ir_core.h"
#include "ir_radiometry.h"

static int run=0,pass=0;
#define T(n) do{run++;printf("  TEST %s... ",n);}while(0)
#define P() do{pass++;printf("PASS\n");}while(0)
#define C(c) do{if(!(c)){printf("FAIL\n");return 1;}}while(0)

int main(void) {
  printf("=== Test: ir_radiometry ===\n");
  ir_source_t src={IR_SOURCE_EXTENDED,300,0.95,1,0.01,1000};
  ir_optical_system_t opt={0.05,0.05,1,0.5,0.5,0.9,2.025e-6};
  T("source rad"); C(ir_source_radiance_at_sensor(&src,0.8,0)>0); P();
  T("ext irrad"); C(ir_extended_source_irradiance(&src,&opt,0.8)>0); P();
  T("power det"); C(ir_power_on_detector(&src,&opt,0.8)>0); P();
  T("etendue"); C(ir_etendue_image_space(1e-4,1)>0); P();
  T("solid ang"); C(ir_solid_angle_circular(0.1,100)>0); P();
  T("SNR pt"); C(ir_snr_point_source(1000,0.8,0.9,0.002,1e6,1e-9,1000)>0); P();
  T("SNR ext"); C(ir_snr_extended_source(10,0.8,0.9,2.025e-6,1,1e6,1e-9)>0); P();
  T("T2L"); C(ir_temperature_to_radiance(300,8,14,200)>0); P();
  printf("  %d/%d tests passed\n\n",pass,run);
  return (pass==run)?0:1;
}
