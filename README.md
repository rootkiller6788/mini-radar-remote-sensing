# Mini Radar Remote Sensing

A collection of **from-scratch, zero-dependency C implementations** of radar systems, remote sensing, and signal processing theory. Each module maps to MIT, Stanford, and other top-tier university courses, translating textbook equations into runnable C code covering the full radar signal chain — from waveform generation and antenna beamforming to detection, tracking, and imaging.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-radar-basics](mini-radar-basics/) | Radar signal model, detection theory (CFAR, Marcum, Swerling), Doppler/MTI processing, waveform design (LFM, Barker, Costas) | MIT 6.630, Stanford EE359, Skolnik |
| [mini-pulse-doppler](mini-pulse-doppler/) | Radar range equation, waveform generation, Doppler processing, CFAR detection (CA/GO/SO/OS), ambiguity function analysis | Stanford EE359, Georgia Tech ECE 6350, Richards |
| [mini-phased-array](mini-phased-array/) | Array geometry, phase/time-delay beamforming, adaptive beamforming (SMI/MVDR/LCMV), pattern synthesis (Dolph-Chebyshev, Taylor) | Stanford EE252, Van Trees, Balanis |
| [mini-sar-imaging](mini-sar-imaging/) | SAR geometry, focusing algorithms (RDA/CSA/ωK), interferometry (InSAR/coherence), compressive sensing SAR, MIMO-SAR | MIT 16.851, Stanford EE355, Georgia Tech ECE 6350 |
| [mini-infrared-thermal](mini-infrared-thermal/) | IR core physics/laws, atmospheric transmission (MODTRAN), detector models/NUC calibration, radiometry, target detection, image processing | MIT 22.071, Stanford EE290, ETH 227-0436 |
| [mini-hyperspectral](mini-hyperspectral/) | Spectroscopy principles, radiometric calibration, dimensionality reduction (PCA/MNF), classification (SAM/SVM), spectral unmixing (LMM) | MIT 12.710, Stanford EE369, Purdue ECE 637 |
| [mini-lidar-principle](mini-lidar-principle/) | LiDAR detection theory (SNR/Pd/Pfa), point cloud geometry, full-waveform processing, scanning patterns, ICP registration, DEM/forestry apps | MIT 16.851, Stanford EE368, TU Munich HF Eng |
| [mini-target-tracking](mini-target-tracking/) | Kalman filter variants (KF/EKF/UKF/SR-KF), motion models (CV/CA/CT/IMM), data association (NN/GNN/PDA/JPDA/MHT), track fusion, mmWave tracker | Stanford EE363, Georgia Tech ECE 6601, Bar-Shalom |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes referencing standard textbooks (Skolnik, Balanis, Van Trees, Bar-Shalom)
- **Practical demos** — pulse-Doppler processing chains, adaptive beamformer simulators, SAR image formation, multi-target trackers, and more

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-radar-basics
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-radar-remote-sensing/
├── mini-radar-basics/          # Radar fundamentals, detection, Doppler, waveform design
├── mini-pulse-doppler/         # Pulse-Doppler processing, range equation, ambiguity
├── mini-phased-array/          # Phased array beamforming, adaptive arrays, synthesis
├── mini-sar-imaging/           # Synthetic Aperture Radar imaging and interferometry
├── mini-infrared-thermal/      # Infrared physics, detectors, radiometry, image processing
├── mini-hyperspectral/         # Hyperspectral spectroscopy, unmixing, classification
├── mini-lidar-principle/       # LiDAR detection, geometry, waveform, registration
└── mini-target-tracking/       # Kalman filters, motion models, data association, fusion
```

## License

MIT
