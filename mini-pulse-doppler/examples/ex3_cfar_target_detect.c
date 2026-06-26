#include "../include/cfar_detector.h"
#include "../include/radar_range_eq.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Ex3: CFAR Target Detection with Radar Range Equation ===\n\n");

    /* Radar system: X-band airborne radar */
    radar_system_params_t sys = {0};
    sys.tx_power_watts = 1000.0;
    sys.tx_antenna_gain_db = 30.0;
    sys.rx_antenna_gain_db = 30.0;
    sys.wavelength_m = 0.03;
    sys.center_frequency_hz = 10e9;
    sys.noise_figure_db = 4.0;
    sys.system_losses_db = 6.0;
    sys.bandwidth_hz = 10e6;
    sys.pulse_width_s = 10e-6;
    sys.num_coherent_pulses = 16;
    sys.noise_temperature_k = 290.0;

    printf("X-band Radar Parameters:\n");
    printf("  Frequency:  %.1f GHz\n", sys.center_frequency_hz/1e9);
    printf("  Tx Power:   %.0f W (%.0f dBm)\n", sys.tx_power_watts,
           10.0*log10(sys.tx_power_watts*1000.0));
    printf("  Antenna Gain: %.0f dBi\n", sys.tx_antenna_gain_db);
    printf("  Bandwidth:  %.1f MHz\n", sys.bandwidth_hz/1e6);
    printf("  Pulses/CPI: %.0f\n", sys.num_coherent_pulses);
    printf("  Coherent Gain: %.1f dB\n", 10.0*log10(sys.num_coherent_pulses));
    printf("  Range Resolution: %.1f m\n", 299792458.0/(2.0*sys.bandwidth_hz));

    double rcs = 10.0;
    range_eq_result_t res;
    radar_range_snr(&sys, 50000.0, rcs, &res);
    printf("\nAt 50 km, RCS=10 dBsm: SNR = %.1f dB\n", res.snr_db);

    double max_range;
    radar_max_range(&sys, rcs, 13.0, &max_range);
    printf("Max detection range (SNR>=13dB): %.1f km\n", max_range/1e3);

    double loss_db;
    free_space_path_loss(sys.center_frequency_hz, 50000.0, &loss_db);
    printf("Free-space path loss at 50 km: %.1f dB\n", loss_db);

    /* RCS computation */
    double sphere_rcs;
    rcs_sphere_compute(1.0, sys.center_frequency_hz, &sphere_rcs);
    printf("RCS of 1m sphere: %.1f dBsm\n", sphere_rcs);

    double corner_rcs;
    rcs_corner_reflector_compute(0.5, sys.center_frequency_hz, &corner_rcs);
    printf("RCS of 0.5m corner reflector: %.1f dBsm\n", corner_rcs);

    /* CFAR detection simulation */
    printf("\n--- CFAR Detection Simulation ---\n");
    double rd_data[200];
    for (int i = 0; i < 200; i++) rd_data[i] = 1.0 + 0.5 * ((double)rand()/RAND_MAX - 0.5);
    rd_data[50] = 25.0;
    rd_data[120] = 18.0;

    cfar_config_t cfg;
    cfar_config_init(CFAR_CA, 16, 4, 1e-4, &cfg);
    printf("CA-CFAR: %u ref cells, %u guard cells, Pfa=1e-4\n",
           cfg.num_reference_cells, cfg.num_guard_cells);
    printf("Threshold factor alpha = %.3f\n", cfg.threshold_factor);

    cfar_detection_t dets[200];
    size_t nd;
    ca_cfar_detect(rd_data, 200, &cfg, dets, &nd);
    int detected_count = 0;
    for (size_t i = 0; i < nd; i++)
        if (dets[i].detected) detected_count++;
    printf("Detections: %d out of %zu cells tested\n", detected_count, nd);
    for (size_t i = 0; i < nd; i++) {
        if (dets[i].detected)
            printf("  Target at cell %u: val=%.1f, threshold=%.2f, noise=%.2f\n",
                   dets[i].cell_index, dets[i].cell_value,
                   dets[i].threshold_value, dets[i].estimated_noise_power);
    }

    double pd_est;
    cfar_pd_compute(CFAR_CA, res.snr_db, 1e-4, 16, &pd_est);
    printf("Estimated Pd at SNR=%.1f dB: %.3f\n", res.snr_db, pd_est);

    printf("\nDone.\n");
    return 0;
}
