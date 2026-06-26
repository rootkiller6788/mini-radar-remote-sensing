# Course Tree — Prerequisite Dependencies

## Dependency Graph

```
                       mini-phased-array
                       /        |        \
                      /         |         \
           mini-radar-basics   mini-signal-basis  mini-fourier-analysis
                     |                |                  |
                     v                v                  v
              Radar range eq.    Complex numbers    DFT/FFT (spatial)
              SNR, RCS           Vector spaces      Spectral analysis
                                         |
                                         v
                                  mini-filter-theory
                                         |
                                         v
                                  FIR filter = beamforming weight vector
```

## Required Prerequisites

### From mini-signal-system-theory:
1. **mini-signal-basis**: Complex exponential representation, vector spaces
2. **mini-fourier-analysis**: Fourier transform relationship between aperture distribution and radiation pattern
3. **mini-filter-theory**: FIR filter design → weight vector design, window functions
4. **mini-adaptive-filter**: LMS/RLS algorithms → SMI/MVDR adaptive beamforming

### From mini-electromagnetic-wave:
1. Wave propagation (k₀ = 2π/λ)
2. Antenna theory (radiation patterns, directivity, effective aperture)
3. Far-field conditions

### From mini-radar-remote-sensing (mini-radar-basics):
1. Radar range equation (extended to phased array)
2. SNR, RCS concepts
3. Detection theory basics

### From mini-communication-principle:
1. Complex baseband representation
2. Signal-to-interference-plus-noise ratio (SINR)
3. Array processing for spatial filtering

## Knowledge Flow

```
1. Geometry & Coordinates (L1-L3)
   ↓
2. Array Factor & Pattern Multiplication (L4)
   ↓
3. Beam Steering & Synthesis (L5)
   ↓
4. Pattern Analysis & Metrics (L6)
   ↓
5. AESA System Performance (L7)
   ↓
6. Adaptive Beamforming (L8)
   ↓
7. Research Frontiers (L9, documented)
```

## Cross-Module Dependencies

| This Module | Depends On | Relationship |
|-------------|------------|-------------|
| Array factor | Fourier analysis | Array factor is a spatial DFT |
| Steering vector | Complex numbers | exp(j k₀ r·û) calculation |
| Window functions | Filter theory | Same windows (Hamming, Hann, etc.) |
| MVDR beamformer | Adaptive filter | Same optimization: minimize w^H R w |
| Radar range eq. | Radar basics | Extends single-antenna to array |
| Directivity | Antenna theory | Fundamental antenna metric |
| T/R module | Analog electronics | MMIC PA/LNA modeling |
