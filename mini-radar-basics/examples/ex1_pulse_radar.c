/**
 * @file    ex1_pulse_radar.c
 * @brief   Example: Pulse radar range estimation and SNR analysis
 *
 * Demonstrates: Radar range equation, SNR vs range, maximum detection range,
 *               range resolution, unambiguous range/velocity.
 *
 * Knowledge: L4 Radar range equation, L6 Range estimation
 * Reference: Richards, Scheer & Holm (2010), Ch.2
 */
#include <stdio.h>
#include <math.h>
#include "radar_core.h"

int main(void)
{
    printf("=== Example 1: Pulse Radar Range Estimation ===\n\n");

    /* X-band (10 GHz) surveillance radar */
    radar_params_t radar;
    radar_params_init(&radar,
        10e9,        /* 10 GHz */
        50000.0,     /* 50 kW peak */
        35.0,        /* 35 dB gain */
        1e-6,        /* 1 us pulse */
        2000.0,      /* 2 kHz PRF */
        2e6,         /* 2 MHz BW */
        3.0,         /* 3 dB NF */
        4.0,         /* 4 dB losses */
        200.0,       /* 200 K antenna temp */
        RADAR_MODE_PULSE);

    printf("System Parameters:\n");
    printf("  Frequency:      %.1f GHz\n", radar.center_freq_hz / 1e9);
    printf("  Wavelength:     %.3f m\n", radar.wavelength_m);
    printf("  Peak Power:     %.1f kW\n", radar.peak_power_w / 1000.0);
    printf("  Antenna Gain:   %.1f dB (%.1f linear)\n",
           radar.antenna_gain_db, radar.antenna_gain_linear);
    printf("  Pulse Width:    %.1f us\n", radar.pulse_width_s * 1e6);
    printf("  PRF:            %.0f Hz (PRI = %.3f ms)\n",
           radar.prf_hz, radar.pri_s * 1000.0);
    printf("  Bandwidth:      %.1f MHz\n", radar.bandwidth_hz / 1e6);
    printf("  Noise Figure:   %.1f dB\n", radar.noise_figure_db);

    /* Range resolution */
    double delta_r = radar_range_resolution(radar.pulse_width_s, 0);
    double delta_r_bw = radar_range_resolution(radar.bandwidth_hz, 1);
    printf("\nRange Performance:\n");
    printf("  Resolution (pulse):    %.1f m\n", delta_r);
    printf("  Resolution (bandwidth): %.1f m\n", delta_r_bw);
    printf("  Unambiguous Range:     %.1f km\n",
           radar_unambiguous_range(radar.prf_hz) / 1000.0);

    /* SNR vs range for a 1 m^2 target */
    printf("\nSNR vs Range (RCS = 1 m^2, Swerling 0):\n");
    printf("  Range (km)    SNR (dB)    Received Power (dBm)\n");
    printf("  ----------    --------    ---------------------\n");

    double ranges_km[] = {1, 5, 10, 20, 50, 75, 100};
    for (int i = 0; i < 7; i++) {
        double r = ranges_km[i] * 1000.0;
        double snr_db = radar_snr_db(&radar, 1.0, r);
        double pr_dbm = lin2db(radar_received_power(&radar, 1.0, r) * 1000.0);
        printf("  %8.1f     %8.1f     %8.1f\n", ranges_km[i], snr_db, pr_dbm);
    }

    /* Maximum detection range for different SNR thresholds */
    printf("\nMaximum Detection Range (RCS = 1 m^2):\n");
    printf("  SNR_min (dB)    Max Range (km)\n");
    printf("  ------------    --------------\n");
    double snr_min_db[] = {13.0, 10.0, 5.0, 0.0};
    for (int i = 0; i < 4; i++) {
        double snr_min = db2lin(snr_min_db[i]);
        double rmax = radar_max_range(&radar, 1.0, snr_min);
        printf("  %8.1f       %9.1f\n", snr_min_db[i], rmax / 1000.0);
    }

    /* Effect of RCS on max range */
    printf("\nMax Range vs RCS (SNR_min = 13 dB):\n");
    double snr_min_13 = db2lin(13.0);
    double rcs_vals[] = {0.01, 0.1, 1.0, 10.0, 100.0};
    for (int i = 0; i < 5; i++) {
        double rmax = radar_max_range(&radar, rcs_vals[i], snr_min_13);
        printf("  RCS = %6.1f m^2 (%5.1f dBsm) -> R_max = %.1f km\n",
               rcs_vals[i], lin2db(rcs_vals[i]), rmax / 1000.0);
    }

    /* Unambiguous velocity */
    double v_unamb = radar_unambiguous_velocity(radar.wavelength_m, radar.prf_hz);
    printf("\nUnambiguous Velocity: %.1f m/s (%.1f km/h)\n",
           v_unamb, v_unamb * 3.6);

    /* Duty cycle and average power */
    double dc = radar_duty_cycle(radar.pulse_width_s, radar.prf_hz);
    double pavg = radar_average_power(radar.peak_power_w, dc);
    printf("Duty Cycle: %.4f (%.2f%%), Average Power: %.1f W\n",
           dc, dc * 100.0, pavg);

    printf("\nExample complete.\n");
    return 0;
}
