/**
 * @file    demo_radar.c
 * @brief   Demo: Interactive radar range equation calculator
 *
 * Computes and displays SNR, received power, and maximum range
 * for user-specified radar parameters.
 */
#include <stdio.h>
#include <math.h>
#include "radar_core.h"

int main(void)
{
    printf("=== Radar Range Equation Demo ===\n\n");

    radar_params_t radar;
    radar_params_init(&radar,
        10e9, 100000.0, 40.0,
        10e-6, 1000.0, 10e6,
        2.0, 5.0, 150.0,
        RADAR_MODE_PULSE);

    printf("Radar system configured:\n");
    printf("  Frequency: 10 GHz (X-band), Peak power: 100 kW\n");
    printf("  Antenna gain: 40 dB, Pulse width: 10 us\n");
    printf("  PRF: 1 kHz, BW: 10 MHz, NF: 2 dB\n");

    double rcs_m2 = 1.0;
    printf("\nRange sweep (RCS = %.1f m^2):\n", rcs_m2);
    printf("  Range      P_rx (dBm)    SNR (dB)    Status\n");
    printf("  -----      ----------    --------    ------\n");

    double ranges[] = {1e3, 5e3, 10e3, 20e3, 50e3, 100e3, 200e3, 500e3};
    for (int i = 0; i < 8; i++) {
        double r = ranges[i];
        double pr = radar_received_power(&radar, rcs_m2, r);
        double snr_db = radar_snr_db(&radar, rcs_m2, r);
        double pr_dbm = lin2db(pr * 1000.0);
        const char *status = snr_db > 13.0 ? "DETECTABLE" :
                             snr_db > 0.0  ? "MARGINAL"   : "BELOW NOISE";
        printf("  %6.1f km   %8.1f      %7.1f     %s\n",
               r/1000.0, pr_dbm, snr_db, status);
    }

    printf("\nMax range for SNR_min=13 dB: %.1f km\n",
           radar_max_range(&radar, rcs_m2, db2lin(13.0)) / 1000.0);

    double r_res = radar_range_resolution(radar.bandwidth_hz, 1);
    double r_unamb = radar_unambiguous_range(radar.prf_hz);
    double v_unamb = radar_unambiguous_velocity(radar.wavelength_m, radar.prf_hz);

    printf("\nConstraints:\n");
    printf("  Range resolution:    %.1f m\n", r_res);
    printf("  Unambiguous range:   %.1f km\n", r_unamb/1000.0);
    printf("  Unambiguous velocity: %.1f m/s\n", v_unamb);

    printf("\nDemo complete.\n");
    return 0;
}
