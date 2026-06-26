# Knowledge Graph — mini-pulse-doppler

## L1: Definitions (Complete)

| # | Concept | C Type/Define | Lean Definition |
|---|---------|---------------|-----------------|
| 1 | Pulse width (tau) | `radar_waveform_params_t.pulse_width` | `WaveformParams.pulseWidth` |
| 2 | Bandwidth (B) | `radar_waveform_params_t.bandwidth` | `WaveformParams.bandwidth` |
| 3 | PRF / PRI | `radar_waveform_params_t.prf/pri` | `WaveformParams.prf` |
| 4 | Duty cycle | `radar_waveform_params_t.duty_cycle` | — |
| 5 | Barker code | `barker_code_t` | — |
| 6 | LFM chirp rate | `lfm_params_t.chirp_rate` | — |
| 7 | Time-bandwidth product | `lfm_params_t.time_bw_product` | `PulseCompression.timeBandwidthProduct` |
| 8 | SNR | `radar_detection_params_t.snr_db` | `DetectionParams.snrDb` |
| 9 | Pfa / Pd | `radar_detection_params_t.pfa/pd` | `DetectionParams.pfa/pd` |
| 10 | RCS | `rcs_target_t.rcs_value_dbsm` | — |
| 11 | Range resolution | `range_doppler_map_t.range_resolution_m` | `rangeResolution` |
| 12 | Doppler resolution | `range_doppler_map_t.doppler_resolution_hz` | — |
| 13 | Unambiguous range | `range_doppler_map_t.max_unambiguous_range_m` | `maxUnambiguousRange` |
| 14 | Unambiguous velocity | `range_doppler_map_t.max_unambiguous_velocity_mps` | `maxUnambiguousVelocity` |
| 15 | Blind speed | `blind_speed_compute()` | — |
| 16 | Coherent Processing Interval | `cpi_params_t` | — |
| 17 | Swerling cases | `swerling_case_t` | `SwerlingCase` |
| 18 | CFAR threshold factor | `cfar_config_t.threshold_factor` | — |
| 19 | Range-Doppler cell | `range_doppler_cell_t` | — |
| 20 | Radar system parameters | `radar_system_params_t` | — |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Matched filter | `matched_filter_init/apply` |
| 2 | Pulse compression | `pulse_compress_lfm` |
| 3 | LFM chirp | `lfm_chirp_generate` |
| 4 | Barker phase coding | `barker_code_generate` |
| 5 | Coherent integration | `coherent_integration` |
| 6 | Noncoherent integration | `noncoherent_integration` |
| 7 | Range gating | `range_gate_extract` |
| 8 | Doppler processing | `doppler_fft_process` |
| 9 | MTI filtering | `mti_filter_design/apply` |
| 10 | CA-CFAR | `ca_cfar_detect` |
| 11 | GO/SO-CFAR | `go_cfar_detect` / `so_cfar_detect` |
| 12 | OS-CFAR | `os_cfar_detect` |
| 13 | Doppler shift | `velocity_from_doppler` / `doppler_from_velocity` |
| 14 | Range-Doppler coupling | `rd_coupling_analyze` |
| 15 | Window functions | `waveform_apply_window` / `window_function_create` |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Autocorrelation function | `waveform_autocorrelation` |
| 2 | Ambiguity function (2-D) | `ambiguity_function_compute` |
| 3 | Cross-ambiguity (1-D) | `waveform_cross_ambiguity_1d` |
| 4 | Instantaneous frequency | `waveform_instantaneous_freq` |
| 5 | FFT (DFT) for Doppler | `doppler_fft_process` (DFT implementation) |
| 6 | Complex baseband signals | `double complex` throughout |
| 7 | Order statistics | `os_cfar_detect` (sorting) |
| 8 | Variability Index | `vi_cfar_detect` |
| 9 | Free-space path loss | `free_space_path_loss` |
| 10 | Two-ray propagation model | `two_ray_path_loss` |
| 11 | Knife-edge diffraction | `knife_edge_diffraction_loss` |
| 12 | Q-function (ambiguity) | `q_function_compute` |

## L4: Fundamental Laws (Complete)

| # | Law/Theorem | C Implementation | Lean Statement |
|---|-------------|-----------------|----------------|
| 1 | Radar range equation | `radar_range_snr` | `radarRangeEq` |
| 2 | Matched filter theorem | `matched_filter_init` | `matchedFilterOptimality` |
| 3 | Doppler shift formula | `velocity_from_doppler` | `dopplerShift` |
| 4 | Unambiguous velocity limit | `max_unambiguous_velocity` | `maxUnambiguousVelocity` |
| 5 | Unambiguous range limit | — | `maxUnambiguousRange` |
| 6 | Range resolution | — | `rangeResolution` |
| 7 | Albersheim equation | `albersheim_equation` | — |
| 8 | Shnidman equation | `shnidman_equation` | — |
| 9 | Neyman-Pearson criterion | `detection_threshold_compute` | — |
| 10 | Ambiguity product constant | — | `ambiguityProductConstant` |
| 11 | CFAR threshold scaling | — | `cfarThresholdScalesWithNoise` |
| 12 | Friis transmission equation | `friis_transmission` | — |
| 13 | Rayleigh resolution criterion | `waveform_rms_bandwidth` | — |
| 14 | Bistatic range equation | `bistatic_range_snr` | — |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Matched filtering (time domain) | `matched_filter_apply` |
| 2 | Matched filtering (freq domain) | `matched_filter_apply_freq_domain` |
| 3 | Pulse compression | `pulse_compress_lfm` |
| 4 | DFT-based Doppler FFT | `doppler_fft_process` |
| 5 | CA-CFAR algorithm | `ca_cfar_detect` |
| 6 | GO-CFAR algorithm | `go_cfar_detect` |
| 7 | SO-CFAR algorithm | `so_cfar_detect` |
| 8 | OS-CFAR algorithm | `os_cfar_detect` |
| 9 | VI-CFAR algorithm | `vi_cfar_detect` |
| 10 | 2-D CFAR scan | `cfar_2d_scan` |
| 11 | Window functions (8 types) | `window_function_create` |
| 12 | Ambiguity function computation | `ambiguity_function_compute` |
| 13 | RCS computation (4 shapes) | `rcs_sphere/plate/corner/cylinder_compute` |
| 14 | SNR vs range profile | `snr_vs_range_profile` |
| 15 | Blake chart link budget | `blake_chart_compute` |
| 16 | Doppler ambiguity resolution (CRT) | `resolve_doppler_ambiguity_crt` |
| 17 | Staggered PRF design | `stagger_prf_design` |
| 18 | MTI filter bank | `mti_filter_design` |

## L6: Canonical Problems (Complete)

| # | Problem | Example |
|---|---------|---------|
| 1 | LFM pulse compression + range estimation | `ex1_simple_pulse_radar.c` |
| 2 | Range-Doppler map generation + peak detection | `ex2_range_doppler_map.c` |
| 3 | CFAR detection + radar range equation + RCS | `ex3_cfar_target_detect.c` |

## L7: Applications (Partial — 5 applications)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | Airborne radar link budget | `radar_range_eq.c` (SNR vs range) |
| 2 | Weather radar attenuation | `atmospheric_attenuation_compute` |
| 3 | GMTI (Ground Moving Target Indication) | `mti_filter_apply` + `doppler_fft_process` |
| 4 | Automotive radar (77 GHz) | Configurable via `radar_system_params_t` |
| 5 | Staggered PRF for Doppler ambiguity | `stagger_prf_design` + `resolve_doppler_ambiguity_crt` |

## L8: Advanced Topics (Partial — 3 topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Adaptive Doppler filtering (simplified STAP) | `adaptive_doppler_filter` |
| 2 | Variability-Index CFAR | `vi_cfar_detect` |
| 3 | Non-linear FM waveform design | `nlfm_params_t` (structure) |
| 4 | Cognitive radar parameters | `CognitiveRadarParams` (Lean) |

## L9: Research Frontiers (Partial — documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Cognitive radar | Structure defined in Lean; waveform selection heuristic |
| 2 | MIMO radar | Documented, not implemented |
| 3 | Compressive sensing radar | Documented, not implemented |
| 4 | Quantum radar | Mentioned, not implemented |
| 5 | AI-based target classification | Mentioned, not implemented |
