#ifndef RADAR_DOPPLER_H
#define RADAR_DOPPLER_H
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
double radar_doppler_shift(double vr, double lam);
double radar_doppler_to_velocity(double fd, double lam);
double radar_doppler_resolution(double prf, size_t np);
double radar_blind_speed(double lam, double prf, int n);
typedef enum { MTI_SINGLE_CANCEL, MTI_DOUBLE_CANCEL, MTI_TRIPLE_CANCEL, MTI_CUSTOM } mti_filter_type_t;
typedef struct { mti_filter_type_t type; size_t order; double *coeffs; double complex *delay_line; size_t delay_index; double prf_hz; int initialized; } mti_filter_t;
int mti_filter_init(mti_filter_t *f, mti_filter_type_t t, const double *coeffs, size_t order);
double complex mti_filter_apply(mti_filter_t *f, double complex xn);
int mti_filter_response(const mti_filter_t *f, const double *freq, size_t n, double *resp);
void mti_filter_free(mti_filter_t *f);
int radar_doppler_fft_2d(const double complex *data, size_t nr, size_t np, double complex *map);
double radar_doppler_bin_to_freq(size_t bin, size_t nfft, double prf);
double radar_velocity_from_phase(double complex s0, double complex s1, double lam, double pri);
int radar_doppler_resolve_ambiguity(double fd1, double prf1, double fd2, double prf2, double *res);
#endif
