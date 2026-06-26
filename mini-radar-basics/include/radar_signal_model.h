#ifndef RADAR_SIGNAL_MODEL_H
#define RADAR_SIGNAL_MODEL_H
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef CMPLX
#define CMPLX(r,i) ((double complex)((double)(r) + I * (double)(i)))
#endif
typedef struct { double range_m; double velocity_ms; double azimuth_rad; double elevation_rad; double rcs_m2; double doppler_hz; double delay_s; double amplitude; int target_id; } radar_target_t;
int radar_target_init(radar_target_t *t, int id, double r, double v, double az, double el, double rcs, double lam);
size_t radar_range_to_bin(double r, double fs);
double radar_bin_to_range(size_t b, double fs);
size_t radar_num_range_bins(double rmax, double fs);
int radar_generate_return(const double complex *tx, size_t txl, const radar_target_t *t, double fs, double complex *rx, size_t rxl, double n);
int radar_generate_multitarget(const double complex *tx, size_t txl, const radar_target_t *ts, size_t nt, double fs, double complex *rx, size_t rxl, double n);
typedef enum { CLUTTER_RAYLEIGH, CLUTTER_WEIBULL, CLUTTER_LOG_NORMAL, CLUTTER_K_DISTRIBUTION } clutter_distribution_t;
typedef struct { clutter_distribution_t dist; double sigma0_db; double sigma0_linear; double gamma_param; double mean_power; double correlation_length_m; double grazing_angle_rad; } clutter_model_t;
int clutter_model_init(clutter_model_t *c, double s0, clutter_distribution_t d, double gr, double cl);
int clutter_generate(const clutter_model_t *c, size_t rb, double rr, double rs, double complex *s, size_t n);
double clutter_rcs_cell(double s0, double r, double bw_az, double rr);
int radar_awgn(double complex *noise, size_t n, double pwr, unsigned int seed);
int radar_range_profile(const double complex *rx, size_t rxl, const double complex *tx, size_t txl, double *prof);
#endif
