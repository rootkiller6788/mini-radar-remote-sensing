# mini-pulse-doppler

Pulse Doppler Radar Signal Processing — complete implementation covering
waveform design, matched filtering, Doppler processing, CFAR detection,
ambiguity function analysis, and the radar range equation.

## Module Status: COMPLETE

- **L1-L6**: Complete
- **L7**: Partial+ (5 applications: airborne radar, weather radar, GMTI, automotive, staggered PRF)
- **L8**: Partial (3/7 advanced topics: adaptive Doppler filter, VI-CFAR, NLFM structure)
- **L9**: Partial (documented, cognitive radar structure in Lean)

### Line Count Verification

include/ + src/ total: see `make lines` for exact count (target >= 3000).

## Core Definitions (L1)

| Definition | Symbol | C Type |
|-----------|--------|--------|
| Pulse width | tau | `radar_waveform_params_t.pulse_width` |
| Bandwidth | B | `radar_waveform_params_t.bandwidth` |
| PRF | PRF | `radar_waveform_params_t.prf` |
| Duty cycle | tau*PRF | `radar_waveform_params_t.duty_cycle` |
| Time-bandwidth product | tau*B | `lfm_params_t.time_bw_product` |
| SNR | SNR | `radar_detection_params_t.snr_db` |
| Pfa | Pfa | `radar_detection_params_t.pfa` |
| RCS | sigma | `rcs_target_t.rcs_value_dbsm` |
| Range resolution | Delta_R | `range_doppler_map_t.range_resolution_m` |
| Doppler resolution | Delta_fd | `range_doppler_map_t.doppler_resolution_hz` |
| Unambiguous range | Rmax | `maxUnambiguousRange` |
| Unambiguous velocity | vmax | `maxUnambiguousVelocity` |
| Blind speed | v_blind | `blind_speed_compute` |

## Core Theorems (L4)

| Theorem | Formula | Code |
|---------|---------|------|
| Radar Range Equation | Pr = Pt*Gt*Gr*lambda^2*sigma / ((4pi)^3*R^4) | `radar_range_snr` |
| Matched Filter Theorem | h(t) = s*(-t) maximizes SNR | `matched_filter_init` |
| Doppler Shift | fd = 2*v/lambda | `velocity_from_doppler` |
| Unambiguous Velocity | vmax = lambda*PRF/4 | `max_unambiguous_velocity` |
| Unambiguous Range | Rmax = c/(2*PRF) | `maxUnambiguousRange` |
| Range Resolution | Delta_R = c/(2*B) | `rangeResolution` |
| Albersheim Equation | SNR = f(Pfa, Pd, N) | `albersheim_equation` |
| Friis Transmission | Pr = Pt*Gt*Gr*(lambda/(4pi*R))^2 | `friis_transmission` |

## Core Algorithms (L5)

| Algorithm | Function | Complexity |
|-----------|----------|------------|
| Matched Filter (time) | `matched_filter_apply` | O(N*M) |
| Matched Filter (freq) | `matched_filter_apply_freq_domain` | O(N log N) |
| CA-CFAR | `ca_cfar_detect` | O(N*R) |
| GO-CFAR | `go_cfar_detect` | O(N*R) |
| SO-CFAR | `so_cfar_detect` | O(N*R) |
| OS-CFAR | `os_cfar_detect` | O(N*R log R) |
| VI-CFAR | `vi_cfar_detect` | O(N*R) |
| 2-D CFAR | `cfar_2d_scan` | O(NR*ND*R^2) |
| Doppler DFT | `doppler_fft_process` | O(NR*ND^2) |
| Ambiguity Function | `ambiguity_function_compute` | O(ND*NF*Nsig) |
| RCS Sphere | `rcs_sphere_compute` | O(1) |
| RCS Flat Plate | `rcs_flat_plate_compute` | O(1) |
| RCS Corner Reflector | `rcs_corner_reflector_compute` | O(1) |
| Window Functions (8 types) | `window_function_create` | O(N) |

## Canonical Problems (L6)

1. **LFM Pulse Compression + Range Estimation** (`examples/ex1_simple_pulse_radar.c`)
   - Generate LFM chirp -> simulate target echo -> matched filter -> range estimate
2. **Range-Doppler Map Generation** (`examples/ex2_range_doppler_map.c`)
   - Generate CPI -> Doppler FFT -> RD map -> peak detection
3. **CFAR Detection + Link Budget** (`examples/ex3_cfar_target_detect.c`)
   - Radar range equation -> SNR calculation -> CFAR detection -> RCS models

## Course Alignment

| School | Courses | Coverage |
|--------|---------|----------|
| MIT | 6.003, 6.450, 6.630 | LFM, matched filter, radar range eq |
| Stanford | EE359, EE264, EE252 | Doppler, DSP, RCS |
| Berkeley | EE123, EE117 | DFT, pulse compression |
| Illinois | ECE 310, ECE 459 | Window design, detection theory |
| Michigan | EECS 351, 455, 411 | Doppler, ambiguity, propagation |
| Georgia Tech | ECE 4270, 6601, 6350 | MTI, CFAR, RCS |
| TU Munich | SP, Comm, HF | Waveform, detection, systems |
| ETH | 227-0427, 0436, 0455 | Coherent proc, estimation, EM |
| Tsinghua | SigSys, Comm, EM, DSP | Full radar signal processing chain |

## Building and Testing

```bash
make        # build and run tests + examples
make test   # run tests only
make lines  # count lines in include/ and src/
make clean  # remove binaries
```

## File Structure

```
mini-pulse-doppler/
  include/
    radar_waveform.h       - Waveform types, LFM, Barker, pulse train
    pulse_doppler.h        - Matched filter, integration, detection
    doppler_processing.h   - DFT, MTI, window functions, ambiguity resolution
    cfar_detector.h        - CA/GO/SO/OS/VI CFAR algorithms
    ambiguity_function.h   - 2-D ambiguity function computation
    radar_range_eq.h       - Range equation, RCS, propagation models
  src/
    radar_waveform.c       - Waveform generation + analysis
    pulse_doppler.c        - Core pulse Doppler processing
    doppler_processing.c   - Doppler FFT + MTI + ambiguity resolution
    cfar_detector.c        - 5 CFAR variants + 2-D scan
    ambiguity_function.c   - Ambiguity function + metrics
    radar_range_eq.c       - Range equation + RCS + propagation
    pulse_doppler.lean     - Lean 4 formalization (theorems + structures)
  tests/
    test_pulse_doppler.c   - Matched filter, detection, RD map tests
    test_waveform.c        - Waveform generation + correlation tests
    test_doppler.c         - Doppler, MTI, window function tests
    test_cfar.c            - CFAR detection algorithm tests
  examples/
    ex1_simple_pulse_radar.c      - LFM compression + range estimation
    ex2_range_doppler_map.c       - RD map + peak detection
    ex3_cfar_target_detect.c      - CFAR + link budget + RCS
  docs/
    knowledge-graph.md     - L1-L9 knowledge coverage table
    coverage-report.md     - Coverage assessment
    gap-report.md          - Missing items + priorities
    course-alignment.md    - Nine-school curriculum mapping
    course-tree.md         - Prerequisites dependency tree
  Makefile
  README.md
```
