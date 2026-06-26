# Course Alignment — mini-pulse-doppler

## Nine-School Radar Curriculum Mapping

### MIT
- **6.003 Signal Processing**: LFM chirp, matched filter, Fourier transform (waveform_autocorrelation, doppler_fft_process)
- **6.450 Digital Communications**: Barker codes, detection theory (barker_code_generate, albersheim_equation)
- **6.630 EM Waves**: Radar range equation, RCS (radar_range_snr, rcs_sphere_compute)

### Stanford
- **EE359 Wireless**: Doppler shift, multipath propagation (two_ray_path_loss, velocity_from_doppler)
- **EE264 DSP**: Window functions, FFT (window_function_create, doppler_fft_process)
- **EE252 Antennas**: RCS, Friis equation (rcs_flat_plate_compute, friis_transmission)

### Berkeley
- **EE123 DSP**: DFT/FFT, pulse compression, matched filter (doppler_fft_process, pulse_compress_lfm)
- **EE117 EM**: Radar range equation, scattering (radar_max_range, rcs_corner_reflector_compute)

### Illinois
- **ECE 310 DSP**: Signal processing for radar, window design (waveform_apply_window)
- **ECE 459 Comm**: Detection theory, matched filter (shnidman_equation, detection_threshold_compute)

### Michigan
- **EECS 351 DSP**: Doppler processing, FFT (doppler_spectrum_compute)
- **EECS 455 Comm**: Radar waveforms, ambiguity function (ambiguity_function_compute)
- **EECS 411 Microwave**: Propagation, link budget (blake_chart_compute, atmospheric_attenuation_compute)

### Georgia Tech
- **ECE 4270 DSP**: MTI, adaptive filtering (mti_filter_design, adaptive_doppler_filter)
- **ECE 6601 Comm**: CFAR detection theory (ca_cfar_detect, os_cfar_detect, vi_cfar_detect)
- **ECE 6350 EM**: Scattering, RCS prediction (rcs_sphere_compute, bistatic_range_snr)

### TU Munich
- **Signal Processing**: Waveform design, pulse compression (lfm_chirp_generate)
- **Communications**: Radar detection, receiver design (matched_filter_init)
- **High-Frequency Engineering**: Radar systems, antenna (radar_range_snr, radar_horizon_range)

### ETH Zurich
- **227-0427 Signal Processing**: Coherent processing, integration (coherent_integration, compute_coherent_gain)
- **227-0436 Comm**: Detection and estimation (albersheim_equation, range_doppler_map_find_peaks)
- **227-0455 EM**: Propagation models (free_space_path_loss, knife_edge_diffraction_loss)

### Tsinghua (Tsinghua)
- **Signal and Systems**: LFM chirp, matched filter theorem (matched_filter_apply)
- **Communication Principles**: Detection theory, CFAR (cfar_config_init, cfar_2d_scan)
- **EM Fields**: Radar equation, scattering (radar_range_snr, rcs_cylinder_compute)
- **DSP**: Doppler FFT, window design (doppler_fft_process, window_function_create)

## Reference Textbooks

| Textbook | Chapters | Implementation Mapping |
|----------|----------|----------------------|
| Richards, Scheer, Holm "Principles of Modern Radar" (2010) | Ch.2 Radar Equation, Ch.4 Waveforms, Ch.5 Doppler, Ch.15 CFAR | `radar_range_eq.c`, `radar_waveform.c`, `doppler_processing.c`, `cfar_detector.c` |
| Skolnik "Radar Handbook" (2008) | Ch.2 Range Eq, Ch.8 Waveforms, Ch.15 Detection | `radar_range_snr`, `detection_threshold_compute` |
| Levanon & Mozeson "Radar Signals" (2004) | Ch.2-6 Waveforms, Ch.8 Ambiguity | `ambiguity_function.c`, `radar_waveform.c` |
| Barton "Radar System Analysis" (1988) | Ch.4 Detection, Ch.8 CFAR | `albersheim_equation`, `shnidman_equation` |
| Mahafza "Radar Signal Analysis" (2005) | Full text | All modules |
