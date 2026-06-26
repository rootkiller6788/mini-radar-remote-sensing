#include "radar_range_eq.h"
#include "pulse_doppler.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int radar_range_snr(const radar_system_params_t *sys, double range_m,
                    double target_rcs_dbsm, range_eq_result_t *result)
{
    if (sys == NULL || result == NULL || range_m <= 0.0) return -1;
    double lambda = sys->wavelength_m;
    double Pt = sys->tx_power_watts;
    double Gt = pow(10.0, sys->tx_antenna_gain_db / 10.0);
    double Gr = pow(10.0, sys->rx_antenna_gain_db / 10.0);
    double rcs = pow(10.0, target_rcs_dbsm / 10.0);
    double F = pow(10.0, sys->noise_figure_db / 10.0);
    double Ls = pow(10.0, sys->system_losses_db / 10.0);
    double k = 1.380649e-23;
    double T0 = sys->noise_temperature_k;
    double B = sys->bandwidth_hz;
    double tau = sys->pulse_width_s;
    double Np = sys->num_coherent_pulses;

    double denom = pow(4.0*M_PI, 3.0) * pow(range_m, 4.0) * k * T0 * B * F * Ls;
    if (denom < 1e-300) return -1;
    double Pr = (Pt * Gt * Gr * lambda*lambda * rcs * tau * Np) / denom;
    double N0 = k * T0 * B * F;
    double snr_linear = (Pr * tau * Np) / (N0 + 1e-300);

    result->snr_db = 10.0 * log10(snr_linear + 1e-300);
    result->received_power_dbm = 10.0 * log10(Pr * 1000.0);
    result->snr_at_range_db = result->snr_db;
    result->propagation_loss_db = 20.0 * log10(4.0*M_PI*range_m/lambda);
    result->one_way_loss_db = 10.0 * log10(4.0*M_PI*range_m*range_m/(Gt*lambda*lambda/(4.0*M_PI)));
    result->max_detection_range_m = range_m;
    return 0;
}

int radar_max_range(const radar_system_params_t *sys, double target_rcs_dbsm,
                    double min_snr_db, double *max_range_m)
{
    if (sys == NULL || max_range_m == NULL || min_snr_db <= 0.0) return -1;

    double lambda = sys->wavelength_m;
    double Pt = sys->tx_power_watts;
    double Gt = pow(10.0, sys->tx_antenna_gain_db/10.0);
    double Gr = pow(10.0, sys->rx_antenna_gain_db/10.0);
    double rcs = pow(10.0, target_rcs_dbsm/10.0);
    double F = pow(10.0, sys->noise_figure_db/10.0);
    double Ls = pow(10.0, sys->system_losses_db/10.0);
    double k = 1.380649e-23;
    double T0 = sys->noise_temperature_k;
    double B = sys->bandwidth_hz;
    double tau = sys->pulse_width_s;
    double Np = sys->num_coherent_pulses;
    double min_snr_linear = pow(10.0, min_snr_db/10.0);

    double numerator = Pt * Gt * Gr * lambda*lambda * rcs * tau * Np;
    double denominator = pow(4.0*M_PI, 3.0) * k * T0 * B * F * Ls * min_snr_linear;

    if (denominator < 1e-300 || numerator < 1e-300) return -1;
    *max_range_m = pow(numerator / denominator, 0.25);
    return 0;
}

int bistatic_range_snr(const radar_system_params_t *sys,
                       double tx_range_m, double rx_range_m,
                       double bistatic_rcs_dbsm, range_eq_result_t *result)
{
    if (sys == NULL || result == NULL || tx_range_m <= 0.0 || rx_range_m <= 0.0)
        return -1;

    double lambda = sys->wavelength_m;
    double Pt = sys->tx_power_watts;
    double Gt = pow(10.0, sys->tx_antenna_gain_db/10.0);
    double Gr = pow(10.0, sys->rx_antenna_gain_db/10.0);
    double rcs = pow(10.0, bistatic_rcs_dbsm/10.0);
    double F = pow(10.0, sys->noise_figure_db/10.0);
    double Ls = pow(10.0, sys->system_losses_db/10.0);
    double k = 1.380649e-23;
    double T0 = sys->noise_temperature_k;
    double B = sys->bandwidth_hz;
    double tau = sys->pulse_width_s;
    double Np = sys->num_coherent_pulses;

    double denom = pow(4.0*M_PI, 3.0) * tx_range_m*tx_range_m * rx_range_m*rx_range_m
                   * k * T0 * B * F * Ls;
    if (denom < 1e-300) return -1;
    double Pr = (Pt * Gt * Gr * lambda*lambda * rcs * tau * Np) / denom;

    result->snr_db = 10.0*log10(Pr/(k*T0*B*F) + 1e-300);
    result->received_power_dbm = 10.0*log10(Pr*1000.0);
    result->max_detection_range_m = sqrt(tx_range_m * rx_range_m);
    return 0;
}

int radar_horizon_range(double antenna_height_m, double target_height_m,
                        double *horizon_range_m)
{
    if (horizon_range_m == NULL || antenna_height_m < 0.0 || target_height_m < 0.0)
        return -1;
    double Re = 4.0/3.0 * 6371000.0;
    *horizon_range_m = sqrt(2.0*Re*antenna_height_m) + sqrt(2.0*Re*target_height_m);
    return 0;
}

int friis_transmission(double tx_power_dbm, double tx_gain_db,
                       double rx_gain_db, double frequency_hz,
                       double range_m, double *rx_power_dbm)
{
    if (rx_power_dbm == NULL || range_m <= 0.0 || frequency_hz <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double fspl = pow(4.0*M_PI*range_m/lambda, 2.0);
    double fspl_db = 10.0 * log10(fspl);
    *rx_power_dbm = tx_power_dbm + tx_gain_db + rx_gain_db - fspl_db;
    return 0;
}

int free_space_path_loss(double frequency_hz, double range_m, double *loss_db)
{
    if (loss_db == NULL || range_m <= 0.0 || frequency_hz <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    *loss_db = 20.0 * log10(4.0 * M_PI * range_m / lambda);
    return 0;
}

int two_ray_path_loss(double tx_height_m, double rx_height_m,
                      double range_m, double frequency_hz,
                      double reflection_coefficient, double *loss_db)
{
    if (loss_db == NULL || range_m <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double direct = range_m;
    double reflected = sqrt(range_m*range_m + 4.0*tx_height_m*rx_height_m);
    double phase_diff = 2.0*M_PI*(reflected - direct)/lambda;
    double complex E = 1.0 + reflection_coefficient*(cos(phase_diff)-sin(phase_diff)*I);
    double mag = sqrt(creal(E)*creal(E)+cimag(E)*cimag(E));
    double fspl = 20.0*log10(4.0*M_PI*range_m/lambda);
    *loss_db = fspl - 20.0*log10(mag + 1e-30);
    return 0;
}

int knife_edge_diffraction_loss(double obstacle_height_m,
                                double tx_height_m, double rx_height_m,
                                double d1_m, double d2_m,
                                double frequency_hz, double *loss_db)
{
    if (loss_db == NULL || d1_m <= 0.0 || d2_m <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double h = obstacle_height_m;
    double h1 = tx_height_m, h2 = rx_height_m;
    double h_clear = h - (h1*d2_m + h2*d1_m)/(d1_m + d2_m);
    double v = h_clear * sqrt(2.0*(d1_m+d2_m)/(lambda*d1_m*d2_m));

    if (v <= -1.0) *loss_db = 0.0;
    else if (v <= 0.0) *loss_db = 20.0*log10(0.5 - 0.62*v);
    else if (v <= 1.0) *loss_db = 20.0*log10(0.5*exp(-0.95*v));
    else if (v <= 2.4) *loss_db = 20.0*log10(0.4-sqrt(0.1184-(0.38-0.1*v)*(0.38-0.1*v)));
    else *loss_db = 20.0*log10(0.225/v);
    return 0;
}

int rcs_sphere_compute(double radius_m, double frequency_hz, double *rcs_dbsm)
{
    if (rcs_dbsm == NULL || radius_m <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double k = 2.0*M_PI/lambda;
    double ka = k * radius_m;

    double rcs;
    if (ka < 0.1)
        rcs = M_PI*radius_m*radius_m * (64.0/9.0)*pow(ka,4.0);
    else if (ka < 10.0)
        rcs = M_PI*radius_m*radius_m * (1.0 + 0.5/(ka*ka));
    else
        rcs = M_PI * radius_m * radius_m;

    *rcs_dbsm = 10.0 * log10(rcs);
    return 0;
}

int rcs_flat_plate_compute(double width_m, double height_m,
                           double frequency_hz, double angle_deg,
                           double *rcs_dbsm)
{
    if (rcs_dbsm == NULL || width_m <= 0.0 || height_m <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double area = width_m * height_m;
    double rcs_max = 4.0*M_PI*area*area/(lambda*lambda);
    double angle_rad = angle_deg * M_PI / 180.0;
    double x = (M_PI*width_m/lambda)*sin(angle_rad);
    double pattern;
    if (fabs(x) < 1e-6) pattern = 1.0;
    else pattern = sin(x)/x;
    double rcs = rcs_max * pattern*pattern * cos(angle_rad)*cos(angle_rad);
    *rcs_dbsm = 10.0*log10(rcs + 1e-30);
    return 0;
}

int rcs_corner_reflector_compute(double edge_length_m,
                                 double frequency_hz, double *rcs_dbsm)
{
    if (rcs_dbsm == NULL || edge_length_m <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double rcs = 12.0*M_PI*pow(edge_length_m,4.0)/(lambda*lambda);
    *rcs_dbsm = 10.0*log10(rcs);
    return 0;
}

int rcs_cylinder_compute(double radius_m, double length_m,
                         double frequency_hz, double angle_deg,
                         double *rcs_dbsm)
{
    if (rcs_dbsm == NULL || radius_m <= 0.0 || length_m <= 0.0) return -1;
    double c = 299792458.0;
    double lambda = c / frequency_hz;
    double angle_rad = angle_deg * M_PI / 180.0;
    double x = (2.0*M_PI*length_m/lambda)*sin(angle_rad);
    double pattern;
    if (fabs(x) < 1e-6) pattern = 1.0;
    else pattern = sin(x)/x;
    double rcs = (2.0*M_PI*radius_m*length_m*length_m/lambda)*pattern*pattern;
    *rcs_dbsm = 10.0*log10(rcs + 1e-30);
    return 0;
}

int atmospheric_attenuation_compute(double frequency_ghz,
                                    double range_km, double temperature_c,
                                    double humidity_percent,
                                    double rain_rate_mm_per_hr,
                                    atmospheric_atten_t *atten)
{
    if (atten == NULL || range_km < 0.0 || frequency_ghz <= 0.0) return -1;
    double f = frequency_ghz;

    double oxy = 0.0;
    if (f > 1.0 && f < 350.0)
        oxy = 0.00719 * (f*f / (1.0 + 0.0001*f*f)) * (1.0 + 2.0/((f-60.0)*(f-60.0)+1.0));

    double wv = 0.0;
    if (f > 1.0 && f < 350.0) {
        double rho = humidity_percent * 7.5 / (temperature_c + 237.3);
        wv = 0.000035 * rho * f*f * (1.0/((f-22.2)*(f-22.2)+2.0) + 8.5/((f-183.3)*(f-183.3)+8.0));
    }

    double rain = 0.0;
    if (rain_rate_mm_per_hr > 0.0)
        rain = 0.001 * pow(rain_rate_mm_per_hr, 1.1) * pow(f, 0.7);

    atten->oxygen_atten_db_per_km = oxy;
    atten->water_vapor_atten_db_per_km = wv;
    atten->rain_atten_db_per_km = rain;
    atten->fog_atten_db_per_km = 0.0;
    atten->total_atmospheric_atten_db = (oxy + wv + rain) * range_km;
    atten->frequency_ghz = f;
    atten->rain_rate_mm_per_hr = rain_rate_mm_per_hr;
    return 0;
}

int blake_chart_compute(const blake_chart_params_t *params,
                        double target_rcs_dbsm, double range_m,
                        range_eq_result_t *result)
{
    if (params == NULL || result == NULL || range_m <= 0.0) return -1;

    double tx_erp = params->tx_power_dbm + params->tx_antenna_gain_db
                    - params->tx_line_loss_db - params->radome_loss_db;
    double rx_sensitivity = -174.0 + params->rx_noise_figure_db
                            + 10.0*log10(1.0);
    double total_gain = params->tx_antenna_gain_db + params->rx_antenna_gain_db
                       + params->integration_gain_db + params->signal_processing_gain_db
                       + target_rcs_dbsm;
    double total_loss = params->tx_line_loss_db + params->rx_line_loss_db
                       + params->radome_loss_db + params->propagation_factor_db;

    result->snr_db = tx_erp + total_gain - total_loss - 30.0*log10(4.0*M_PI)
                     - 40.0*log10(range_m) + 20.0*log10(3e8/1e9) - rx_sensitivity;
    result->received_power_dbm = params->tx_power_dbm + total_gain - total_loss
                                 - 20.0*log10(4.0*M_PI*range_m*1e9/3e8);
    result->max_detection_range_m = range_m;
    return 0;
}

int snr_vs_range_profile(const radar_system_params_t *sys,
                         double target_rcs_dbsm,
                         double r_min_m, double r_max_m,
                         uint32_t num_points,
                         double *ranges_m, double *snr_values_db)
{
    if (sys == NULL || ranges_m == NULL || snr_values_db == NULL) return -1;
    if (r_min_m <= 0.0 || r_max_m <= r_min_m || num_points < 2) return -1;

    double step = (r_max_m - r_min_m) / (double)(num_points - 1);
    for (uint32_t i = 0; i < num_points; i++) {
        ranges_m[i] = r_min_m + step * (double)i;
        range_eq_result_t res;
        if (radar_range_snr(sys, ranges_m[i], target_rcs_dbsm, &res) == 0)
            snr_values_db[i] = res.snr_db;
        else
            snr_values_db[i] = -999.0;
    }
    return 0;
}
