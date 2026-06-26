/**
 * @file    lidar_core.c
 * @brief   Core LiDAR implementation — range equation, TOF, configs, atmosphere
 *
 * Knowledge covered:
 *   L1: Time-of-flight, range equation, point cloud management,
 *       LiDAR system configurations (automotive, airborne, terrestrial)
 *   L3: Range equation mathematical structure, photon statistics
 *   L4: LiDAR range equation, Koschmieder's law, Beer-Lambert law,
 *       Poisson photon statistics
 *
 * Reference:
 *   - Jelalian, A.V., *Laser Radar Systems*, Artech House, 1992.
 *   - Wehr & Lohr, *ISPRS JPRS* 54(2-3), pp.68-82, 1999.
 *   - Kruse et al., *Elements of Infrared Technology*, Wiley, 1962.
 */

#include "lidar_core.h"
#include "lidar_detection.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ─── Internal helpers ──────────────────────────────────────────────────── */

static void* safe_malloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "lidar_core: malloc(%zu) failed\n", sz); abort(); }
    return p;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Time-of-Flight conversions
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_tof_to_range(double tof) {
    if (tof <= 0.0) return 0.0;
    return LIDAR_C * tof / 2.0;
}

double lidar_range_to_tof(double range) {
    if (range <= 0.0) return 0.0;
    return 2.0 * range / LIDAR_C;
}

double lidar_unambiguous_range(double prf) {
    if (prf <= 0.0) return 0.0;
    return LIDAR_C / (2.0 * prf);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Point cloud management
 * ═══════════════════════════════════════════════════════════════════════════ */

void lidar_scan_init(lidar_scan_t *scan, size_t capacity) {
    if (!scan) return;
    memset(scan, 0, sizeof(*scan));
    if (capacity == 0) capacity = 1024;
    scan->points = (lidar_point_t*)safe_malloc(capacity * sizeof(lidar_point_t));
    scan->capacity = capacity;
    scan->num_points = 0;
}

int lidar_scan_add_point(lidar_scan_t *scan, const lidar_point_t *point) {
    if (!scan || !point) return -1;
    if (scan->num_points >= scan->capacity) {
        size_t new_cap = scan->capacity * 2;
        if (new_cap == 0) new_cap = 64;
        lidar_point_t *new_pts = (lidar_point_t*)realloc(scan->points,
            new_cap * sizeof(lidar_point_t));
        if (!new_pts) return -1;
        scan->points = new_pts;
        scan->capacity = new_cap;
    }
    scan->points[scan->num_points++] = *point;
    return 0;
}

void lidar_scan_free(lidar_scan_t *scan) {
    if (!scan) return;
    free(scan->points);
    memset(scan, 0, sizeof(*scan));
}

int lidar_scan_bounding_box(const lidar_scan_t *scan,
                             double min_xyz[3], double max_xyz[3]) {
    if (!scan || scan->num_points == 0 || !min_xyz || !max_xyz) return -1;
    double xmin = scan->points[0].x, xmax = scan->points[0].x;
    double ymin = scan->points[0].y, ymax = scan->points[0].y;
    double zmin = scan->points[0].z, zmax = scan->points[0].z;
    for (size_t i = 1; i < scan->num_points; i++) {
        double x = scan->points[i].x;
        double y = scan->points[i].y;
        double z = scan->points[i].z;
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
    }
    min_xyz[0] = xmin; min_xyz[1] = ymin; min_xyz[2] = zmin;
    max_xyz[0] = xmax; max_xyz[1] = ymax; max_xyz[2] = zmax;
    return 0;
}

size_t lidar_scan_size(const lidar_scan_t *scan) {
    return scan ? scan->num_points : 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Configuration initialization
 * ═══════════════════════════════════════════════════════════════════════════ */

void lidar_config_init_automotive(lidar_config_t *config) {
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->type            = LIDAR_TYPE_PULSED;
    config->wavelength      = LIDAR_WAVELENGTH_905NM;
    config->pulse_energy    = 1.0e-6;       /* 1 uJ per pulse */
    config->pulse_width     = 5.0e-9;       /* 5 ns FWHM */
    config->prf             = 1.0e6;        /* 1 MHz */
    config->beam_divergence = 0.003;        /* 3 mrad */
    config->beam_diameter   = 0.005;        /* 5 mm aperture */
    config->detector        = LIDAR_DETECTOR_APD;
    config->aperture_diam   = 0.030;        /* 30 mm receiver */
    config->opt_transmit    = 0.90;
    config->opt_receive     = 0.85;
    config->detector_neb    = 1.0e-14;      /* 10 fW/rtHz */
    config->detector_bw     = 200.0e6;      /* 200 MHz */
    config->excess_noise_F  = 3.0;          /* APD typical */
    config->responsivity    = 30.0;         /* 30 A/W at 905 nm with APD gain */
    config->scan_type       = LIDAR_SCAN_OPTO_MECH;
    config->scan_fov_h      = 6.283185307;  /* 360 deg */
    config->scan_fov_v      = 0.523599;     /* 30 deg */
    config->scan_rate       = 20.0;         /* 20 Hz */
    config->ang_res_h       = 0.00349;      /* 0.2 deg */
    config->ang_res_v       = 0.007;        /* 0.4 deg */
    config->sample_rate     = 2.0e9;        /* 2 GS/s */
    config->range_min       = 0.5;
    config->range_max       = 200.0;        /* 200 m max */
}

void lidar_config_init_airborne(lidar_config_t *config) {
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->type            = LIDAR_TYPE_PULSED;
    config->wavelength      = LIDAR_WAVELENGTH_1064NM;
    config->pulse_energy    = 1.0e-3;       /* 1 mJ */
    config->pulse_width     = 5.0e-9;       /* 5 ns */
    config->prf             = 100.0e3;      /* 100 kHz */
    config->beam_divergence = 0.0005;       /* 0.5 mrad */
    config->beam_diameter   = 0.010;        /* 10 mm */
    config->detector        = LIDAR_DETECTOR_APD;
    config->aperture_diam   = 0.150;        /* 15 cm receiver */
    config->opt_transmit    = 0.80;
    config->opt_receive     = 0.75;
    config->detector_neb    = 5.0e-15;
    config->detector_bw     = 100.0e6;      /* 100 MHz */
    config->excess_noise_F  = 2.5;
    config->responsivity    = 40.0;
    config->scan_type       = LIDAR_SCAN_GALVO;
    config->scan_fov_h      = 1.047198;     /* 60 deg */
    config->scan_fov_v      = 1.047198;     /* 60 deg */
    config->scan_rate       = 50.0;         /* 50 lines/s */
    config->ang_res_h       = 0.0005;
    config->ang_res_v       = 0.0005;
    config->sample_rate     = 1.0e9;        /* 1 GS/s */
    config->range_min       = 50.0;
    config->range_max       = 6000.0;       /* 6 km */
}

void lidar_config_init_terrestrial(lidar_config_t *config) {
    if (!config) return;
    memset(config, 0, sizeof(*config));
    config->type            = LIDAR_TYPE_PHASE_BASED;
    config->wavelength      = LIDAR_WAVELENGTH_1550NM;
    config->pulse_energy    = 1.0e-7;       /* 100 nJ */
    config->pulse_width     = 1.0e-9;       /* 1 ns */
    config->prf             = 500.0e3;      /* 500 kHz */
    config->beam_divergence = 0.0002;       /* 0.2 mrad */
    config->beam_diameter   = 0.008;        /* 8 mm */
    config->detector        = LIDAR_DETECTOR_PIN;
    config->aperture_diam   = 0.040;        /* 40 mm */
    config->opt_transmit    = 0.90;
    config->opt_receive     = 0.85;
    config->detector_neb    = 1.0e-13;
    config->detector_bw     = 500.0e6;      /* 500 MHz */
    config->excess_noise_F  = 1.0;          /* PIN — no excess noise */
    config->responsivity    = 1.0;          /* 1 A/W at 1550 nm */
    config->scan_type       = LIDAR_SCAN_GALVO;
    config->scan_fov_h      = 6.283185307;  /* 360 deg */
    config->scan_fov_v      = 5.235988;     /* 300 deg (nearly full sphere) */
    config->scan_rate       = 100.0;        /* 100 Hz */
    config->ang_res_h       = 0.0001;
    config->ang_res_v       = 0.0001;
    config->sample_rate     = 10.0e9;       /* 10 GS/s */
    config->range_min       = 0.3;
    config->range_max       = 1000.0;       /* 1 km */
}

int lidar_config_validate(const lidar_config_t *config) {
    if (!config) return -1;
    if (config->pulse_width <= 0.0) return -2;
    if (config->prf <= 0.0) return -3;
    if (config->beam_divergence <= 0.0) return -4;
    if (config->aperture_diam <= 0.0) return -5;
    if (config->pulse_energy <= 0.0) return -6;
    if (config->detector_bw <= 0.0) return -7;
    if (config->range_max <= config->range_min) return -8;
    if (config->range_min < 0.0) return -9;
    /* Nyquist: sample_rate >= 2 / pulse_width */
    if (config->sample_rate < 2.0 / config->pulse_width) return -10;
    if (config->opt_transmit <= 0.0 || config->opt_transmit > 1.0) return -11;
    if (config->opt_receive <= 0.0 || config->opt_receive > 1.0) return -12;
    if (config->responsivity <= 0.0) return -13;
    if (config->detector_neb < 0.0) return -14;
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L4: LiDAR Range Equation
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_range_equation_received_power(const lidar_config_t *config,
                                            const lidar_target_t *target,
                                            const lidar_atmosphere_t *atm,
                                            double range) {
    if (!config || !target || !atm || range <= 0.0) return 0.0;

    /* Peak transmitted power */
    double P_t = config->pulse_energy / config->pulse_width;

    /* Receiver aperture area */
    double A_r = M_PI * (config->aperture_diam / 2.0) * (config->aperture_diam / 2.0);

    /* System optical efficiency */
    double eta_sys = config->opt_transmit * config->opt_receive;

    /* Two-way atmospheric transmission (Beer-Lambert) */
    double T2 = exp(-2.0 * atm->extinction_coeff * range);

    /* For point target with Lambertian reflectance:
       sigma = 4 * rho * A * cos(theta_normal) */
    double sigma;
    if (target->is_distributed) {
        /* Distributed target: illuminated area */
        double A_illum = M_PI * config->beam_divergence * config->beam_divergence
                         * range * range / 4.0;
        sigma = 4.0 * target->reflectivity * A_illum;
    } else {
        double cos_theta = cos(target->normal_angle);
        if (cos_theta < 0.1) cos_theta = 0.1; /* grazing incidence floor */
        sigma = 4.0 * target->reflectivity * target->area * cos_theta;
    }

    /* LiDAR range equation:
       P_r = P_t * sigma * A_r * eta_sys * T2 / (4 * pi * R^2 * omega_t * R^2)
       where beam solid angle omega_t = pi * (theta_div/2)^2 */
    double omega_t = M_PI * (config->beam_divergence / 2.0) * (config->beam_divergence / 2.0);
    if (omega_t < 1e-20) omega_t = 1e-20;

    double denominator = 4.0 * M_PI * range * range * omega_t * range * range;
    double P_r = P_t * sigma * A_r * eta_sys * T2 / denominator;

    return (P_r > 0.0) ? P_r : 0.0;
}

double lidar_max_range(const lidar_config_t *config,
                       const lidar_target_t *target,
                       const lidar_atmosphere_t *atm,
                       double snr_min_db) {
    if (!config || !target || !atm || snr_min_db <= 0.0) return 0.0;

    double snr_min_linear = pow(10.0, snr_min_db / 10.0);
    double P_t = config->pulse_energy / config->pulse_width;
    double A_r = M_PI * (config->aperture_diam / 2.0) * (config->aperture_diam / 2.0);
    double eta_sys = config->opt_transmit * config->opt_receive;

    /* Approximate minimum detectable power:
       P_min = k * T * B / (R_det * M^2) * SNR_min  (simplified)
       Using detector NEP: P_min = NEP * sqrt(B) * SNR_min */
    double NEP = config->detector_neb * sqrt(config->detector_bw);
    double P_min = NEP * snr_min_linear;

    double sigma;
    if (target->is_distributed) {
        sigma = 4.0 * target->reflectivity * target->area;
    } else {
        double cos_theta = cos(target->normal_angle);
        if (cos_theta < 0.1) cos_theta = 0.1;
        sigma = 4.0 * target->reflectivity * target->area * cos_theta;
    }

    double omega_t = M_PI * (config->beam_divergence / 2.0) * (config->beam_divergence / 2.0);
    if (omega_t < 1e-20) omega_t = 1e-20;

    /* Newton-Raphson for R⁴ when extinction is present.
       The equation: P_r(R) = P_min
       where P_r(R) = K * exp(-2*beta*R) / R^4 */
    double beta = atm->extinction_coeff;
    double K = P_t * sigma * A_r * eta_sys / (4.0 * M_PI * omega_t);

    /* Initial guess: ignore extinction */
    double R = pow(K / (P_min * 4.0 * M_PI), 0.25);
    if (R < 1.0 || !isfinite(R)) R = 100.0;

    /* Newton iterations: f(R) = P_r(R) - P_min = 0 */
    for (int iter = 0; iter < 30; iter++) {
        double exp_term = exp(-2.0 * beta * R);
        if (exp_term < 1e-100) { R = 1.0 / (2.0 * beta); break; }

        double f = K * exp_term / (R * R * R * R) - P_min;
        /* f'(R) = -K * exp(-2beta*R) * (4/R^5 + 2*beta/R^4) */
        double fp = -K * exp_term * (4.0 / pow(R, 5.0) + 2.0 * beta / pow(R, 4.0));

        if (fabs(fp) < 1e-50) break;
        double dR = f / fp;
        R -= dR;
        if (R < 0.1) R = 0.1;
        if (fabs(dR) < 1e-3) break;
    }

    return (R > 0.0 && isfinite(R)) ? R : 0.0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: SNR computation
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_snr(const lidar_config_t *config,
                 const lidar_target_t *target,
                 const lidar_atmosphere_t *atm,
                 double range) {
    if (!config || !target || !atm || range <= 0.0) return -999.0;

    double P_r = lidar_range_equation_received_power(config, target, atm, range);
    if (P_r <= 0.0) return -999.0;

    /* Signal current */
    double I_sig = P_r * config->responsivity;

    /* Background current */
    /* Solar spectral radiance at wavelength — approximate */
    double L_solar;
    if (config->wavelength == 905)   L_solar = 500.0;  /* W/(m^2*sr*um) */
    else if (config->wavelength == 1064) L_solar = 300.0;
    else if (config->wavelength == 1550) L_solar = 150.0;
    else                              L_solar = 400.0;

    double delta_lambda = 0.01e-6; /* 10 nm optical filter (in meters as um) */
    /* Actually using nm for filter: 10 nm = 0.01 um */
    double A_r = M_PI * (config->aperture_diam / 2.0) * (config->aperture_diam / 2.0);
    double Omega_r = M_PI * (config->scan_fov_v / 2.0) * (config->scan_fov_h / 2.0);
    if (Omega_r > M_PI) Omega_r = M_PI;  /* cap at hemisphere */
    if (Omega_r < 1e-10) Omega_r = 1e-10;

    double P_bg = L_solar * delta_lambda * A_r * Omega_r * config->opt_receive;
    double I_bg = P_bg * config->responsivity;

    /* Dark current (typical values) */
    double I_dark = 1.0e-9; /* 1 nA typical */

    /* APD gain M: treat as 1 for PIN, ~50 for APD */
    double M = (config->detector == LIDAR_DETECTOR_PIN) ? 1.0 : 50.0;
    double F = config->excess_noise_F;

    /* Shot noise variance */
    double sigma_shot2 = 2.0 * LIDAR_ELECTRON_Q * (I_sig + I_bg + I_dark)
                         * M * M * F * config->detector_bw;

    /* Thermal noise: T=300K, R_fb=1kOhm typical */
    double T = 300.0;
    double R_fb = 1000.0;
    double sigma_thermal2 = 4.0 * LIDAR_BOLTZMANN_K * T * config->detector_bw / R_fb;

    double noise_power = sigma_shot2 + sigma_thermal2;
    if (noise_power <= 0.0) return -999.0;

    double signal_power = (I_sig * M) * (I_sig * M);
    double snr_linear = signal_power / noise_power;

    if (snr_linear <= 0.0) return -999.0;

    return 10.0 * log10(snr_linear);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Range resolution
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_range_resolution(const lidar_config_t *config) {
    if (!config) return 0.0;
    if (config->type == LIDAR_TYPE_FMCW) {
        /* For FMCW: delta_R = c / (2 * B_chirp), but we don't have chirp BW
           in the base config.  Use detector BW as proxy for chirp bandwidth. */
        double B_chirp = config->detector_bw;
        if (B_chirp <= 0.0) return 0.0;
        return LIDAR_C / (2.0 * B_chirp);
    }
    /* Pulsed LiDAR */
    return LIDAR_C * config->pulse_width / 2.0;
}

double lidar_cross_range_resolution(const lidar_config_t *config, double range) {
    if (!config || range <= 0.0) return 0.0;
    /* Beam-footprint limited */
    double beam_fp = range * config->beam_divergence;
    /* Angular sampling limited */
    double ang_samp = range * config->ang_res_h;
    /* Take the larger of the two */
    return (beam_fp > ang_samp) ? beam_fp : ang_samp;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L1: Atmosphere parameters
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_atmosphere_extinction(double visibility, double wavelength) {
    if (visibility <= 0.0 || wavelength <= 0.0) return 0.0;

    /* Convert wavelength from nm to um */
    double lambda_um = wavelength / 1000.0;

    /* Koschmieder: beta(550nm) = 3.912 / V */
    double beta_550 = 3.912 / visibility;

    /* Angstrom exponent correction q */
    double q;
    if (visibility > 50.0e3) {
        q = 1.6;
    } else if (visibility >= 6.0e3) {
        q = 1.3;
    } else {
        q = 0.585 * pow(visibility / 1000.0, 1.0/3.0);
    }

    /* Wavelength scaling: beta(lambda) = beta(0.55) * (0.55/lambda)^q */
    double beta = beta_550 * pow(0.55 / lambda_um, q);

    return beta;
}

double lidar_atmosphere_transmission(double extinction_coeff, double range) {
    if (extinction_coeff < 0.0 || range < 0.0) return 0.0;
    double T_single = exp(-extinction_coeff * range);
    double T_two_way = T_single * T_single;
    return (T_two_way >= 0.0 && T_two_way <= 1.0) ? T_two_way : 0.0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L4: Poisson photon statistics
 * ═══════════════════════════════════════════════════════════════════════════ */

double lidar_expected_signal_pe(const lidar_config_t *config,
                                 double range,
                                 const lidar_target_t *target,
                                 const lidar_atmosphere_t *atm) {
    if (!config || !target || !atm || range <= 0.0) return 0.0;

    double P_r = lidar_range_equation_received_power(config, target, atm, range);
    if (P_r <= 0.0) return 0.0;

    /* Photon energy: E_photon = h * c / lambda */
    double lambda_m = config->wavelength * 1e-9; /* nm -> m */
    double E_photon = LIDAR_PLANCK_H * LIDAR_C / lambda_m;

    /* Number of photons in the pulse: N_photons = P_r * tau_pulse / E_photon */
    double N_photons = P_r * config->pulse_width / E_photon;

    /* Quantum efficiency: approximate
       Si @ 905nm ~ 30%, InGaAs @ 1550nm ~ 70% */
    double qe;
    if (config->wavelength <= 1064) {
        qe = 0.30; /* Silicon */
    } else {
        qe = 0.70; /* InGaAs */
    }

    return N_photons * qe;
}

double lidar_poisson_prob(int k, double lambda) {
    if (k < 0 || lambda < 0.0) return 0.0;
    if (lambda == 0.0) return (k == 0) ? 1.0 : 0.0;

    /* Compute P(k) = lambda^k * exp(-lambda) / k!
       Use log-space for numerical stability */
    double log_p = k * log(lambda) - lambda - lgamma(k + 1.0);
    return exp(log_p);
}

int lidar_compute_performance(const lidar_config_t *config,
                               const lidar_target_t *target,
                               const lidar_atmosphere_t *atm,
                               double range,
                               lidar_performance_t *perf) {
    if (!config || !target || !atm || !perf || range <= 0.0) return -1;

    memset(perf, 0, sizeof(*perf));

    perf->received_power   = lidar_range_equation_received_power(config, target, atm, range);
    perf->snr              = lidar_snr(config, target, atm, range);
    perf->range_resolution = lidar_range_resolution(config);
    perf->ang_resolution   = config->ang_res_h;
    perf->max_range        = lidar_max_range(config, target, atm, 10.0);
    perf->cross_range_res  = lidar_cross_range_resolution(config, range);

    /* NEFD: Noise Equivalent Flux Density [W/m^2] */
    double A_r = M_PI * (config->aperture_diam / 2.0) * (config->aperture_diam / 2.0);
    perf->nefd = config->detector_neb * sqrt(config->detector_bw) / A_r;

    /* Range accuracy from CRLB approximation */
    double snr_linear = pow(10.0, perf->snr / 10.0);
    double tau = config->pulse_width;
    double b_rms = 0.187 / tau;
    perf->range_accuracy = LIDAR_C / (2.0 * sqrt(2.0) * sqrt(snr_linear) * b_rms);

    /* Probabilities for single-photon detection */
    double N_sig = lidar_expected_signal_pe(config, range, target, atm);
    if (config->detector == LIDAR_DETECTOR_SPAD || config->detector == LIDAR_DETECTOR_SIPM) {
        double N_dark = 1e-6; /* typical dark count rate per range bin */
        perf->prob_false_alarm = 1.0 - exp(-N_dark);
        perf->prob_detection   = 1.0 - exp(-(N_sig + N_dark + 1e-6));
    } else {
        /* Linear mode: Gaussian approximation */
        if (perf->snr > -900.0) {
            perf->prob_false_alarm = 0.01; /* default PFA */
            perf->prob_detection   = lidar_prob_detection(perf->snr, perf->prob_false_alarm);
        }
    }

    return 0;
}