# Knowledge Graph — mini-hyperspectral

## L1 — Definitions: COMPLETE ✅

| # | Item | C Type / Constant | Description |
|---|------|-------------------|-------------|
| 1 | Hyperspectral datacube | `hspec_datacube_t` | 3D spatial-spectral array (x, y, λ) |
| 2 | Spectral band | `hspec_band_t` | Narrow contiguous wavelength interval with FWHM |
| 3 | Spectral Response Function (SRF) | `hspec_srf_t` | Gaussian detector sensitivity model |
| 4 | Single-pixel spectrum | `hspec_pixel_t` | B-dimensional reflectance vector |
| 5 | Spectral library entry | `hspec_spectral_lib_entry_t` | Reference spectrum (USGS, JPL) |
| 6 | Sensor noise model | `hspec_noise_model_t` | Photon/thermal/readout noise characterization |
| 7 | Radiometric quantities | `hspec_radiance_t` | At-sensor radiance L [W·m⁻²·sr⁻¹·μm⁻¹] |
| 8 | TOA/Surface reflectance | `hspec_reflectance_t` | Top-of-atmosphere and surface reflectance |
| 9 | Detection metrics | `hspec_detection_metrics_t` | SNR, NEdT, NESR, contrast |
| 10 | Spectral indices | `hspec_spectral_indices_t` | NDVI, NDWI, EVI, SAVI, PRI, Red Edge |
| 11 | Atmospheric windows | Constants | VIS (0.38-0.75), SWIR (1.0-2.5), MWIR (3-5), LWIR (8-14) |
| 12 | Spectral absorption feature | `hspec_absorption_feature_t` | Depth, FWHM, area, asymmetry |
| 13 | Continuum-removed spectrum | `hspec_continuum_removed_t` | R_CR(λ) = R(λ) / continuum(λ) |
| 14 | Blackbody radiance | `hspec_blackbody_t` | Planck function B_λ(λ, T) |
| 15 | Endmember | `hspec_endmember_t` | Pure material spectral signature |
| 16 | Abundance map | `hspec_abundance_map_t` | Fractional coverage per pixel |
| 17 | Confusion matrix | `hspec_confusion_matrix_t` | Classification accuracy (OA, κ) |
| 18 | Solar irradiance | `hspec_solar_irradiance_t` | Thuillier 2004 model E_sun(λ) |
| 19 | BRDF model | `hspec_brdf_params_t` | Minnaert, Ross-Li kernel-driven |

## L2 — Core Concepts: COMPLETE ✅

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Imaging spectroscopy | Datacube manipulation (extract pixel/band) |
| 2 | Sensor geometries | `hspec_sensor_geometry_t` (pushbroom/whiskbroom/staring/Fourier) |
| 3 | Linear mixing model (LMM) | `hspec_fcls_estimate`, `hspec_ncls_estimate` |
| 4 | Additivity + non-negativity constraints | FCLS solver with active-set method |
| 5 | Atmospheric correction chain | DN→radiance→TOA reflectance→surface reflectance |
| 6 | Spectral redundancy / Hughes phenomenon | PCA, MNF dimensionality reduction |
| 7 | Atmospheric absorption bands | H₂O, CO₂, O₃ parameterized modeling |
| 8 | Photon/thermal/readout noise | `hspec_noise_model_t` characterization |
| 9 | Dark object subtraction (DOS) | `hspec_dos_correction` |
| 10 | Spectral library matching | SAM classification against reference spectra |

## L3 — Mathematical Structures: COMPLETE ✅

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Tensor/3D array | BIP-layout datacube: data[band + row·B + col·rows·B] |
| 2 | Covariance matrix Σ | `hspec_datacube_covariance` — O(N·B²) |
| 3 | Correlation matrix ρ | `hspec_datacube_correlation` — Σ_ij/(σ_i·σ_j) |
| 4 | Eigendecomposition | `hspec_symmetric_eigen` — QR iteration |
| 5 | Convex hull geometry | VCA, N-FINDR (simplex volume maximization) |
| 6 | Constrained least squares | FCLS via augmented NNLS |
| 7 | Non-negative Matrix Factorization | NMF via multiplicative updates (Lee-Seung) |
| 8 | Mahalanobis distance | `hspec_mahalanobis_distance` — (x-μ)ᵀΣ⁻¹(x-μ) |
| 9 | Gaussian mixture / Bayes decision | ML classification framework |
| 10 | Spectral dissimilarity metrics | Euclidean, SAM, SID, Correlation, Chebyshev |

## L4 — Fundamental Laws: COMPLETE ✅

| # | Law | Formula | C Implementation | Lean Statement |
|---|-----|---------|------------------|----------------|
| 1 | Planck's law | B_λ = (2hc²/λ⁵)/(e^{hc/λkT}-1) | `hspec_planck_radiance()` | — |
| 2 | Wien's displacement | λ_max = b/T | `hspec_wien_peak_wavelength()` | `wien_peak_inverse_proportional` |
| 3 | Stefan-Boltzmann | M = σT⁴ | `hspec_stefan_boltzmann_exitance()` | `stefan_boltzmann_monotone` |
| 4 | Beer-Lambert-Bouguer | I = I₀·e^{-αd} | `hspec_beer_lambert_transmission()` | `optical_depth_additive` |
| 5 | Kirchhoff's law | ε = 1-ρ (opaque) | `hspec_emissivity_from_reflectance()` | `kirchhoff_opaque` |
| 6 | Inverse square law (solar) | E ∝ 1/d² | `hspec_radiance_to_toa_refl()` | — |
| 7 | Lambert cosine law | E ∝ cos θ | TOA reflectance computation | — |
| 8 | Rayleigh quotient bound | wᵀΣw ≤ λ_max·‖w‖² | PCA derivation | `rayleigh_quotient_bound` |
| 9 | Convex combination closure | λa+(1-λ)b valid | Abundance space | `convex_combination_of_abundances` |
| 10 | NDVI bounds | -1 ≤ NDVI ≤ 1 | Spectral indices | `ndvi_bounds` |

## L5 — Algorithms: COMPLETE ✅

| # | Algorithm | Complexity | Implementation |
|---|-----------|-----------|---------------|
| 1 | Spectral Angle Mapper (SAM) | O(nclasses·nbands) | `hspec_sam_classify()` |
| 2 | Matched Filter (MF) | O(nbands²) | `hspec_matched_filter_detect()` |
| 3 | Adaptive Coherence Estimator (ACE) | O(nbands²) | `hspec_ace_detect()` |
| 4 | Reed-Xiaoli (RX) anomaly | O(N·nbands²+nbands³) | `hspec_rx_detect_global()` |
| 5 | K-means clustering (Lloyd) | O(N·k·B·iters) | `hspec_kmeans_cluster()` |
| 6 | VCA endmember extraction | O(N·B·p) | `hspec_vca_extract()` |
| 7 | N-FINDR endmember extraction | O(N·p²·B·iters) | `hspec_nfindr_extract()` |
| 8 | PPI (Pixel Purity Index) | O(n_proj·N·B) | `hspec_ppi_compute()` |
| 9 | FCLS (active-set NNLS) | O(p²·(B+1)·iters) | `hspec_fcls_estimate()` |
| 10 | PCA (QR eigendecomposition) | O(N·B² + B³) | `hspec_pca_compute()` |
| 11 | MNF transform | O(N·B² + B³) | `hspec_mnf_compute()` |
| 12 | FastICA (fixed-point) | O(N·B·p·iters) | `hspec_fastica()` |
| 13 | NMF (multiplicative updates) | O(N·B·r·iters) | `hspec_nmf_factorize()` |
| 14 | Continum removal (convex hull) | O(nbands²) | `hspec_continuum_removal()` |
| 15 | Absorption feature detection | O(nbands) | `hspec_detect_absorption_features()` |
| 16 | Band selection (greedy MI) | O(k·B·N) | `hspec_band_selection_mi()` |

## L6 — Canonical Problems: COMPLETE ✅

| # | Problem | Implementation / Example |
|---|---------|-------------------------|
| 1 | Mineral detection & mapping | `example_mineral_detection.c` — AVIRIS Cuprite simulation |
| 2 | Vegetation health monitoring | `example_vegetation_analysis.c` — NDVI/EVI/SAVI/PRI |
| 3 | Spectral unmixing / subpixel mapping | `example_spectral_unmixing.c` — N-FINDR + FCLS |
| 4 | Target detection (subpixel) | MF and ACE detectors with CFAR thresholds |
| 5 | Anomaly detection | RX global anomaly detector with chi² threshold |
| 6 | Atmospheric correction | DOS + simplified 6S radiative transfer |
| 7 | Water vapor retrieval | CIBR method at 940nm absorption |
| 8 | Spectral resampling | SRF convolution for sensor simulation |
| 9 | Classification accuracy assessment | Confusion matrix with OA and Cohen's κ |
| 10 | Dimensionality reduction | PCA/MNF for Hughes phenomenon mitigation |

## L7 — Applications: COMPLETE ✅ (3+)

| # | Application | Key References |
|---|-------------|---------------|
| 1 | NASA AVIRIS mineral exploration | Cuprite, Nevada (Clark et al. 1993) |
| 2 | Precision agriculture / crop monitoring | NDVI/PRI for vegetation stress detection |
| 3 | Urban land cover / impervious surface mapping | Subpixel spectral unmixing |
| 4 | Water quality / water vapor estimation | CIBR water vapor retrieval (Gao & Goetz 1990) |
| 5 | Atmospheric parameter retrieval | Aerosol optical depth, water vapor column |

## L8 — Advanced Topics: PARTIAL ✅ (2+/5)

| # | Topic | Status |
|---|-------|--------|
| 1 | Fan bilinear mixing model | ✅ Implemented: `hspec_fan_mixing_model()` |
| 2 | Hapke intimate mixing model | ✅ Implemented: `hspec_hapke_reflectance()` |
| 3 | Kernel methods (K-SVM, KPCA) | 📋 Documented |
| 4 | Deep learning for HS classification | 📋 Documented |
| 5 | Sparse unmixing (OMP-based) | 📋 Documented |

## L9 — Research Frontiers: PARTIAL 📋

| # | Topic | Description |
|---|-------|-------------|
| 1 | Real-time onboard HS processing | CubeSat/Drone edge computing |
| 2 | Snapshot hyperspectral imaging | CASSI, CTIS technologies |
| 3 | Quantum illumination HS | Entangled photon spectroscopy |
| 4 | 6G RIS-integrated spectral sensing | Reconfigurable intelligent surfaces |
| 5 | Semantic spectral communication | Joint sensing and communications |

## Summary

| Level | Status | Items |
|-------|--------|-------|
| L1 Definitions | COMPLETE | 19 definitions |
| L2 Core Concepts | COMPLETE | 10 concepts |
| L3 Math Structures | COMPLETE | 10 structures |
| L4 Fundamental Laws | COMPLETE | 10 laws (C + Lean) |
| L5 Algorithms | COMPLETE | 16 algorithms |
| L6 Canonical Problems | COMPLETE | 10 problems |
| L7 Applications | COMPLETE | 5 applications |
| L8 Advanced Topics | PARTIAL (2/5) | Fan + Hapke implemented |
| L9 Research Frontiers | PARTIAL | Documented |