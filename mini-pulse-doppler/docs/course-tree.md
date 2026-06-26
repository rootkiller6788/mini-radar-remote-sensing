# Course Tree — mini-pulse-doppler

## Prerequisites Dependency Tree

```
mini-pulse-doppler
    |
    +-- [0] mini-signal-system-theory
    |       - Fourier Transform (DFT/FFT for Doppler processing)
    |       - Convolution (matched filter implementation)
    |       - Laplace Transform (filter design background)
    |       - Correlation (autocorrelation for range sidelobes)
    |
    +-- [5] mini-communication-principle
    |       - Modulation theory (BPSK for Barker codes)
    |       - Matched filter detection
    |       - Detection theory (Neyman-Pearson)
    |       - Probability of error (Pfa/Pd relationship)
    |
    +-- [6] mini-digital-signal-process
    |       - DFT and FFT algorithms
    |       - FIR/IIR filter design (MTI filters)
    |       - Window functions (Hamming, Hanning, Blackman, etc.)
    |       - Spectral analysis
    |
    +-- [7] mini-electromagnetic-wave
    |       - Radar cross section (RCS)
    |       - Free-space propagation
    |       - Scattering theory (Mie, Rayleigh, optical regions)
    |       - Knife-edge diffraction
    |       - Two-ray ground reflection
    |
    +-- [8] mini-sensor-measurement
    |       - SNR and noise figure
    |       - Detection threshold
    |       - Measurement accuracy
    |
    +-- [11] mini-wireless-mobile-comm
            - Doppler shift in mobile channels
            - Link budget analysis
            - Propagation models
            - Multipath effects

Upward dependencies:
    |
    +-- [14] mini-navigation-positioning
    |       - Radar-based positioning
    |       - Doppler navigation
    |
    +-- [15] mini-iot-edge-computing
            - Radar sensor fusion
            - Real-time radar processing
```

## Knowledge Flow Within Module

```
L1 Definitions
    |
    v
L2 Core Concepts (matched filter, pulse compression, Doppler)
    |
    v
L3 Math Structures (autocorrelation, ambiguity, FFT)
    |
    v
L4 Fundamental Laws (range eq, matched filter theorem, detection theory)
    |
    v
L5 Algorithms (CFAR, windowing, RCS computation, ambiguity function)
    |
    v
L6 Canonical Problems (pulse compression, RD map, CFAR detection)
    |
    v
L7 Applications (GMTI, weather radar, automotive, airborne)
    |
    v
L8 Advanced Topics (STAP, VI-CFAR, NLFM, cognitive radar)
    |
    v
L9 Research (cognitive radar, MIMO, compressive sensing)
```
