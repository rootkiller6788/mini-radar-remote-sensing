# Course Dependency Tree -- mini-sar-imaging

## Prerequisites (what you need to know before this module)

```
Digital Signal Processing
├── Fourier Transform (continuous & discrete) ─────────┐
├── FFT algorithms (Cooley-Tukey) ─────────────────────┤
├── Convolution & correlation ─────────────────────────┤
├── Filter design (FIR/IIR) ───────────────────────────┤
└── Sampling theory (Nyquist-Shannon) ─────────────────┤
                                                         │
Electromagnetic Theory                                   │
├── Maxwell's equations ────────────────────────────────┤
├── Plane wave propagation ─────────────────────────────┤
├── Antenna theory (beam pattern, gain) ────────────────┤
└── Radar cross section (RCS) ──────────────────────────┤
                                                         │
Radar Basics                                            │
├── Radar range equation ───────────────────────────────┤
├── Pulsed radar principle ─────────────────────────────┤
├── Matched filter detection ───────────────────────────┤
└── Doppler effect ─────────────────────────────────────┤
                                                         │
                                                         ▼
                                          mini-sar-imaging (THIS MODULE)
                                                         │
                                                         ▼
                                          ┌──────────────┼──────────────┐
                                          ▼              ▼              ▼
                                    InSAR/PolSAR    GMTI/MTI      SAR Tomography
                                    (interferometry) (moving targets) (3D imaging)
```

## Internal Dependencies

```
sar_core.h ──────┬──> sar_geometry.h ──┬──> sar_algorithm.h
                 │                     │
                 │                     ├──> sar_interferometry.h
                 │                     │
                 │                     └──> sar_advanced.h
                 │
                 └──> sar_core.c ──────> sar_geometry.c ──> sar_algorithm.c
                                                           ──> sar_interferometry.c
                                                           ──> sar_advanced.c
```

## Learning Path

1. **Start with**: sar_core.h/c -- basic SAR definitions, chirp, matched filter
2. **Then**: sar_geometry.h/c -- range equations, Doppler, RCM, antenna
3. **Core algorithms**: sar_algorithm.h/c -- RDA, CSA, omega-k, BP, autofocus
4. **Applications**: sar_interferometry.h/c -- InSAR, DInSAR
5. **Advanced**: sar_advanced.h/c -- CS-SAR, MIMO, bistatic, polarimetric