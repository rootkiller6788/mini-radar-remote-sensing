# mini-hyperspectral

Hyperspectral Remote Sensing submodule — complete implementation of hyperspectral imaging principles: spectral datacube manipulation, spectroscopy (Planck/Wien/Stefan-Boltzmann/Beer-Lambert/Kirchhoff), classification (SAM/MF/ACE/RX/K-means), spectral unmixing (VCA/N-FINDR/FCLS), dimensionality reduction (PCA/MNF/ICA/NMF), and radiometric calibration/atmospheric correction.

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (19 core data structures: datacube, band, pixel, endmember, reflectance, noise model, spectral indices, absorption feature, etc.)
- **L2 Core Concepts**: Complete (10 concepts: imaging spectroscopy, linear mixing, atmospheric correction chain, Hughes phenomenon, spectral library matching)
- **L3 Mathematical Structures**: Complete (covariance/correlation matrices, eigendecomposition, convex hull, Mahalanobis distance, NMF, constrained least squares)
- **L4 Fundamental Laws**: Complete (Planck's law, Wien's displacement, Stefan-Boltzmann, Beer-Lambert-Bouguer, Kirchhoff's law — verified in C + 12 theorem statements in Lean 4)
- **L5 Algorithms**: Complete (16 algorithms: SAM, MF, ACE, RX, K-means, VCA, N-FINDR, PPI, FCLS, PCA, MNF, FastICA, NMF, Continuum Removal, Band Selection, Absorption Feature Detection)
- **L6 Canonical Problems**: Complete (3 end-to-end examples: mineral detection, vegetation analysis, spectral unmixing)
- **L7 Applications**: Complete (5: NASA AVIRIS mineral exploration, precision agriculture, urban land cover, water vapor retrieval, atmospheric parameter retrieval)
- **L8 Advanced Topics**: Partial (2/5: Fan bilinear mixing model, Hapke intimate mixing model)
- **L9 Research Frontiers**: Partial (documented)

## Code Metrics

| Metric | Value |
|--------|-------|
| Header files (.h) | 6 files |
| Source files (.c) | 6 files |
| Lean 4 formalization | 1 file |
| **include/ + src/ total** | **≥ 3000 lines** ✅ |
| Test coverage | 24 test functions |
| Examples | 3 end-to-end examples |
| Docs | 5 knowledge documents |

## Core Definitions (L1)

| Definition | C Type | Description |
|-----------|--------|-------------|
| Hyperspectral datacube | `hspec_datacube_t` | 3D spatial-spectral array (BIP layout) |
| Spectral band | `hspec_band_t` | Center wavelength, FWHM, bounds |
| Spectral response function | `hspec_srf_t` | Gaussian SRF model |
| Single-pixel spectrum | `hspec_pixel_t` | B-dimensional reflectance vector |
| Spectral library entry | `hspec_spectral_lib_entry_t` | USGS/JPL reference spectra |
| Sensor noise model | `hspec_noise_model_t` | Photon/thermal/readout noise |
| At-sensor radiance | `hspec_radiance_t` | L [W·m⁻²·sr⁻¹·μm⁻¹] |
| TOA/Surface reflectance | `hspec_reflectance_t` | ρ_TOA, ρ_surface |
| Detection metrics | `hspec_detection_metrics_t` | SNR, NEdT, NESR per band |
| Spectral indices | `hspec_spectral_indices_t` | NDVI, NDWI, EVI, SAVI, PRI |
| Endmember | `hspec_endmember_t` | Pure material spectral signature |
| Absorption feature | `hspec_absorption_feature_t` | Depth, FWHM, area, asymmetry |
| Confusion matrix | `hspec_confusion_matrix_t` | OA, κ, producer/user accuracy |
| BRDF model | `hspec_brdf_params_t` | Minnaert, Ross-Li kernel-driven |

## Core Theorems (L4)

| Theorem | Formula | Verification |
|---------|---------|-------------|
| **Planck's Law** | B_λ = (2hc²/λ⁵)/(e^{hc/λkT}−1) | `hspec_planck_radiance()` — 5800K solar test |
| **Wien's Displacement** | λ_max = b/T, b = 2.898·10⁻³ m·K | `hspec_wien_peak_wavelength()` — 500nm for Sun |
| **Stefan-Boltzmann** | M = σT⁴ | `hspec_stefan_boltzmann_exitance()` — 459 W/m² at 300K |
| **Beer-Lambert-Bouguer** | I = I₀·e^{-αd} | `hspec_beer_lambert_transmission()` |
| **Kirchhoff's Law** | ε = 1−ρ (opaque) | `hspec_emissivity_from_reflectance()` |
| **NDVI Bounds** | −1 ≤ NDVI ≤ 1 | Formal proof in Lean (`ndvi_bounds`) |
| **Convex Combination** | Closure under convex combination | Formal proof in Lean (`convex_combination_of_abundances`) |

## Core Algorithms (L5)

| Algorithm | Complexity | Implementation |
|-----------|-----------|---------------|
| Spectral Angle Mapper (SAM) | O(C·B) | `hspec_sam_classify()` |
| Matched Filter | O(B²) | `hspec_matched_filter_detect()` |
| Adaptive Coherence Estimator | O(B²) | `hspec_ace_detect()` |
| Reed-Xiaoli Anomaly | O(N·B²+B³) | `hspec_rx_detect_global()` |
| K-means Clustering | O(N·K·B·iters) | `hspec_kmeans_cluster()` |
| VCA Endmember Extraction | O(N·B·P) | `hspec_vca_extract()` |
| N-FINDR | O(N·P²·B·iters) | `hspec_nfindr_extract()` |
| FCLS (active-set NNLS) | O(P²·B·iters) | `hspec_fcls_estimate()` |
| PCA (QR eigendecomposition) | O(N·B²+B³) | `hspec_pca_compute()` |
| MNF Transform | O(N·B²+B³) | `hspec_mnf_compute()` |
| FastICA | O(N·B·P·iters) | `hspec_fastica()` |
| NMF (Lee-Seung) | O(N·B·R·iters) | `hspec_nmf_factorize()` |
| Continuum Removal | O(B²) | `hspec_continuum_removal()` |
| Band Selection (Greedy MI) | O(K·B·N) | `hspec_band_selection_mi()` |

## Nine-School Course Mapping

| School | Key Courses | Mapped Topics |
|--------|-------------|--------------|
| **MIT** | 6.003, 6.450, 6.630 | Signal processing, detection theory, EM waves |
| **Stanford** | EE102A, EE359, EE247 | DSP, wireless detection, optical spectroscopy |
| **Berkeley** | EE16A/B, EE117, EE123 | Sensor physics, EM, DSP |
| **Illinois** | ECE 310, ECE 459, ECE 451 | DSP, communications, EM |
| **Michigan** | EECS 351, 455, 411 | DSP, comm, microwave sensing |
| **Georgia Tech** | ECE 4270, 6601, 6350 | DSP, detection, EM |
| **TU Munich** | Signal Processing, Comm, HF Eng | DSP, detection, atmospheric sensing |
| **ETH** | 227-0427, 0436, 0455 | Signal processing, comm, EM |
| **清华** | 信号与系统, 通信原理, 电磁场, 数字信号处理 | Spectral transforms, detection, radiometry, DSP |

## File Structure

```
mini-hyperspectral/
├── Makefile              # make test → build + run tests
├── README.md             # This file
├── include/              # 6 header files
│   ├── hyperspectral_core.h           # Core types, datacube, statistics
│   ├── hyperspectral_spectroscopy.h   # Planck, Wien, Beer-Lambert, continuum
│   ├── hyperspectral_classification.h # SAM, MF, ACE, RX, K-means
│   ├── hyperspectral_unmixing.h       # VCA, N-FINDR, FCLS, NCLS, bilinear
│   ├── hyperspectral_dimensionality.h # PCA, MNF, ICA, NMF, eigendecomposition
│   └── hyperspectral_radiometry.h     # Calibration, atm correction, BRDF
├── src/                  # 6 C + 1 Lean implementation files
│   ├── hyperspectral_core.c
│   ├── hyperspectral_spectroscopy.c
│   ├── hyperspectral_classification.c
│   ├── hyperspectral_unmixing.c
│   ├── hyperspectral_dimensionality.c
│   ├── hyperspectral_radiometry.c
│   └── hyperspectral_lean.lean        # 12 formal theorems
├── tests/
│   └── test_hyperspectral.c           # 24 comprehensive tests
├── examples/
│   ├── example_mineral_detection.c    # AVIRIS Cuprite mineral mapping
│   ├── example_vegetation_analysis.c  # Precision agriculture / crop health
│   └── example_spectral_unmixing.c    # Urban land cover subpixel mapping
└── docs/
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-alignment.md
    └── course-tree.md
```

## Build & Run

```bash
make          # build everything (test + examples)
make test     # run test suite
make examples # build examples only
make clean    # clean build artifacts
```
