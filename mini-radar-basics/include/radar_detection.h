#ifndef RADAR_DETECTION_H
#define RADAR_DETECTION_H
#include <stddef.h>
#include <stdint.h>
#include <math.h>
typedef enum { DETECT_RESULT_TRUE_POSITIVE, DETECT_RESULT_FALSE_POSITIVE, DETECT_RESULT_TRUE_NEGATIVE, DETECT_RESULT_FALSE_NEGATIVE } detection_result_t;
typedef struct { double prob_detection; double prob_false_alarm; double prob_miss; double snr_db; size_t n_pulses; double threshold; size_t n_trials; } detection_metrics_t;
typedef enum { CFAR_CA, CFAR_GO, CFAR_SO, CFAR_OS, CFAR_CENSORED } cfar_type_t;
typedef struct { cfar_type_t type; size_t num_guard_cells; size_t num_ref_cells; double pfa_desired; double threshold_factor; size_t os_rank; int use_log_detector; } cfar_config_t;
int cfar_config_init(cfar_config_t *cfg, cfar_type_t t, size_t ng, size_t nr, double pfa);
int cfar_detect(const double *rp, size_t n, const cfar_config_t *cfg, int *dets);
int cfar_os_detect(const double *rp, size_t n, const cfar_config_t *cfg, int *dets);
double radar_detection_threshold(double pfa, size_t np);
double radar_pd_marcum(double snr, double pfa, size_t np);
double radar_pd_swerling1(double snr, double pfa);
double radar_pd_swerling2(double snr, double pfa, size_t np);
double radar_albersheim_snr(double pd, double pfa, size_t np);
double radar_marcum_q(int m, double a, double b);
int radar_roc_curve(double snr, double *pfa, double *pd, size_t n, size_t np, int sw);
double radar_gamma_inc(double a, double x);
void radar_dsort(double *arr, size_t n);
#endif
