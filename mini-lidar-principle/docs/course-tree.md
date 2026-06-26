# Course Tree — mini-lidar-principle

## Prerequisites (what this module depends on)

```
mini-signal-system-theory
  ├── mini-fourier-analysis        → Matched filtering, convolution
  └── mini-laplace-z-transform     → Transfer functions, stability

mini-electromagnetic-wave         → Beam propagation, Maxwell equations
mini-sensor-measurement           → SNR, noise analysis
mini-communication-principle      → Detection theory, Neyman-Pearson
mini-digital-signal-process       → Filtering, signal processing algorithms
```

## This Module

```
mini-lidar-principle
  │
  ├── L1: Core definitions
  │     ├── LiDAR types, wavelengths, detectors, scanners
  │     ├── Point cloud, TOF, range equation structures
  │     └── Performance metrics, atmosphere models
  │
  ├── L2: Core concepts
  │     ├── TOF measurement principle
  │     ├── Detection threshold, CFAR, matched filter
  │     └── Optical noise sources
  │
  ├── L3: Mathematical structures
  │     ├── Vector algebra, matrix operations
  │     ├── SVD, Jacobi eigenvalue decomposition
  │     ├── Gaussian beam propagation
  │     └── Poisson statistics, Q-function
  │
  ├── L4: Fundamental laws
  │     ├── LiDAR range equation
  │     ├── Beer-Lambert, Koschmieder laws
  │     ├── Poisson photon statistics
  │     ├── Neyman-Pearson detection
  │     └── CRLB for range estimation
  │
  ├── L5: Algorithms
  │     ├── Peak detection (CFD, derivative, leading-edge)
  │     ├── Gaussian decomposition (Levenberg-Marquardt)
  │     ├── Normal estimation (PCA / Jacobi)
  │     ├── RANSAC plane fitting
  │     ├── ICP point cloud registration
  │     ├── CA-CFAR detection
  │     ├── Matched filter
  │     └── DEM/DSM/CHM generation
  │
  ├── L6: Canonical problems
  │     ├── Range estimation
  │     ├── Point cloud alignment
  │     ├── Ground filtering
  │     ├── Building extraction
  │     └── Waveform decomposition
  │
  ├── L7: Applications
  │     ├── Autonomous vehicle perception
  │     ├── Terrain mapping (DEM)
  │     ├── Forestry metrics
  │     ├── Urban 3D modeling
  │     └── Atmospheric LiDAR
  │
  ├── L8: Advanced topics
  │     ├── Geiger-mode / single-photon LiDAR
  │     ├── FMCW LiDAR
  │     └── Multi-spectral LiDAR
  │
  └── L9: Research frontiers
        ├── Quantum LiDAR (documented)
        ├── Chip-scale LiDAR (documented)
        └── 4D LiDAR (partially implemented via FMCW Doppler)
```

## Dependents (what depends on this module)

```
mini-radar-remote-sensing
  ├── mini-pulse-doppler     → Shares range equation, detection theory
  ├── mini-sar-imaging       → 2D similar to LiDAR raster scan
  └── mini-target-tracking   → Shares point cloud registration (ICP)

mini-navigation-positioning  → GPS/INS + LiDAR fusion

mini-automotive-electronics  → Autonomous vehicle LiDAR processing
