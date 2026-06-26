# mini-radar-basics — Radar Fundamentals

**Module 13 — Radar Remote Sensing · Submodule: Radar Basics**

Implementation of core radar principles: range equation, waveform design, detection theory, Doppler processing, and signal modeling.

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (10+ applications: cognitive radar, radar bands, clutter models, STAP, MIMO)
- **L8**: Partial (NLFM, MIMO virtual array, MVDR, MUSIC, STAP sampling)
- **L9**: Partial (cognitive radar infrastructure present)

**Line count**: 3114 lines in include/ + src/ (>= 3000)
**Tests**: 108 passed, 0 failed
**Safety**: No filler patterns, no TODO/FIXME/stub/placeholder

## Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| L1 | Definitions | Complete | Radar range equation, RCS, PRF/PRI, range resolution, Doppler, antenna gain, SNR, noise figure, waveform types, CFAR window |
| L2 | Core Concepts | Complete | Time-of-flight, pulse/CW/FMCW, matched filter, pulse compression, CFAR principle, Doppler effect, system noise model |
| L3 | Math Structures | Complete | Range equation = Friis power budget, complex baseband, matched filter = correlation, chi-square/gamma, phase progression, ambiguity function |
| L4 | Fundamental Laws | Complete | Radar range equation, Neyman-Pearson lemma, Albersheim equation, Marcum Q, unambiguous range/Doppler constraints |
| L5 | Algorithms | Complete | CA/OS/GO/SO-CFAR, LFM generation, Barker/Frank codes, matched filter (time + FFT), MTI filters, Doppler FFT, beamforming, PSL/ISL |
| L6 | Canonical Problems | Complete | Pulse radar range, FMCW range-velocity, CFAR detection, multi-target returns, range-Doppler map, clutter rejection |
| L7 | Applications | Complete | Automotive radar (77 GHz), weather radar, ATC radar, cognitive radar, radar bands, chaff/sea/rain clutter, STC, micro-Doppler |
| L8 | Advanced Topics | Partial | NLFM, MIMO virtual array, MVDR/Capon, MUSIC DOA, STAP |
| L9 | Research Frontiers | Partial | Cognitive radar perception-action cycle, 6G joint communication-sensing |

## Core Definitions (L1)

- Radar Range Equation: P_r = P_t * G^2 * lambda^2 * sigma / ((4*pi)^3 * R^4 * L)
- Range Resolution: Delta_R = c*tau/2 (pulse) or Delta_R = c/(2*B) (modulated)
- Unambiguous Range: R_unamb = c/(2*PRF)
- Unambiguous Velocity: v_unamb = lambda * PRF / 4
- Doppler Shift: f_d = 2 * v_r / lambda
- Antenna Gain: G = 4*pi * A_e / lambda^2
- Duty Cycle: D = tau * PRF
- Pulse Compression Ratio: PCR = tau * B (TBP)

## Core Theorems (L4)

| Theorem | Formula | Function |
|---------|---------|----------|
| Radar Range Equation | P_r = P_t*G^2*lambda^2*sigma/((4*pi)^3*R^4*L) | radar_received_power() |
| Maximum Detection Range | R_max = [...]^{1/4} | radar_max_range() |
| Bistatic Range Equation | P_r = P_t*G_t*G_r*lambda^2*sigma/((4*pi)^3*R_t^2*R_r^2*L) | radar_bistatic_power() |
| Albersheim's Equation | SNR_dB = f(P_d, P_fa, N) | radar_albersheim_snr() |
| Marcum Q-Function | Q_M(a,b) | radar_marcum_q() |
| Neyman-Pearson Criterion | Lambda(x) threshold test | radar_detection_threshold() |

## Core Algorithms (L5)

- CA-CFAR, OS-CFAR, GO/SO-CFAR constant false alarm rate detection
- LFM chirp generation: s(t) = exp(j*pi*k*t^2)
- Barker codes (length 2,3,4,5,7,11,13) and Frank polyphase codes
- Matched filter: time-domain correlation + FFT-based fast convolution
- MTI filter: single/double/triple delay-line canceller
- Doppler FFT: 2D range-Doppler map processing
- ULA beamforming: array factor, steering vectors, beam pattern

## Nine-School Curriculum Mapping

| School | Key Courses | Topics |
|--------|------------|--------|
| MIT | 6.003, 6.450 | Signal model, matched filter, detection theory |
| Stanford | EE359, EE264 | Waveform design, pulse compression, arrays |
| Berkeley | EE117, EE123 | EM propagation, RCS, DSP for radar |
| Michigan | EECS 411, 455 | Automotive radar, FMCW, MIMO |
| Georgia Tech | ECE 6601, 6350 | Detection theory, CFAR, EM scattering |
| Illinois | ECE 310, 459 | DSP, communication-radar convergence |
| TU Munich | HF Engineering | mmWave radar, antenna arrays |
| ETH | 227-0455 | EM theory, RCS, wave propagation |
| Tsinghua | 信号与系统, 雷达原理 | Radar fundamentals, signal processing |

## Build & Test

```bash
make          # Build library, tests, examples
make test     # Run 108 tests (108 passed, 0 failed)
make examples # Build all examples
make clean    # Remove build artifacts
```

## Reference Textbooks

- Richards, Scheer & Holm, *"Principles of Modern Radar"* (2010)
- Skolnik, *"Radar Handbook"* (2008)
- Levanon & Mozeson, *"Radar Signals"* (2004)
- Kay, *"Fundamentals of Statistical Signal Processing: Detection Theory"* (1998)
- Peebles, *"Radar Principles"* (1998)

---

## Module Status: COMPLETE ✅

Total: 3114 lines include/ + src/ · 108 tests passed · No filler · No TODO/FIXME
