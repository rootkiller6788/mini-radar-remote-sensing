# Knowledge Graph — mini-lidar-principle

## L1: Definitions (Complete ✅)

| # | Concept | C Representation | Lean Definition |
|---|---------|-----------------|-----------------|
| 1 | LiDAR range equation | lidar_range_equation_received_power() | eceived_power_simple |
| 2 | Time-of-flight | lidar_tof_to_range(), lidar_range_to_tof() | 	of_to_range, ange_to_tof |
| 3 | Range resolution | lidar_range_resolution() | ange_resolution_pulsed |
| 4 | Angular resolution | ng_res_h, ng_res_v in lidar_config_t | — |
| 5 | Beam divergence | eam_divergence in lidar_config_t | — |
| 6 | Pulse repetition frequency | prf in lidar_config_t | — |
| 7 | Laser pulse energy | pulse_energy in lidar_config_t | peak_power |
| 8 | Point cloud | lidar_point_t, lidar_scan_t | — |
| 9 | LiDAR types | lidar_type_t enum | LidarType inductive |
| 10 | Wavelength bands | lidar_wavelength_t enum | LidarWavelength inductive |
| 11 | Detector types | lidar_detector_type_t enum | DetectorType inductive |
| 12 | Scan types | lidar_scan_type_t enum | ScanType inductive |
| 13 | SNR | lidar_snr() | — |
| 14 | Probability of detection | lidar_prob_detection() | geiger_detection_prob |
| 15 | Probability of false alarm | lidar_prob_false_alarm field | geiger_false_alarm_prob |
| 16 | NEP / NEFD | detector_neb, perf.nefd | — |
| 17 | Cross-range resolution | lidar_cross_range_resolution() | — |
| 18 | Gaussian beam parameters | lidar_gaussian_beam_t | GaussianComponent |
| 19 | Waveform representation | lidar_waveform_t | — |
| 20 | DEM/DSM/CHM | lidar_dem_t | — |

## L2: Core Concepts (Complete ✅)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Time-of-flight measurement | lidar_tof_to_range(), lidar_tof_to_range() |
| 2 | Pulsed vs FMCW LiDAR | LIDAR_TYPE_PULSED, LIDAR_TYPE_FMCW, lidar_fmcw_beat_frequency() |
| 3 | Scanning mechanisms | 6 scan types in lidar_scan_type_t + scan pattern generation |
| 4 | Single-photon vs linear detection | SPAD, SiPM, APD, PIN in lidar_detector_type_t |
| 5 | Multiple returns | eturn_num, 
um_returns in lidar_point_t |
| 6 | Full waveform digitization | lidar_waveform_t + Gaussian decomposition |
| 7 | Detection threshold | lidar_detection_threshold_gaussian() |
| 8 | CFAR | lidar_cfar_detection() |
| 9 | Matched filter | lidar_matched_filter() |
| 10 | Optical noise sources | lidar_shot_noise_variance(), lidar_thermal_noise_variance() |

## L3: Mathematical Structures (Complete ✅)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | LiDAR range equation | lidar_range_equation_received_power() — full form with beam solid angle |
| 2 | Gaussian function | lidar_gaussian_eval(), lidar_multigaussian_eval() |
| 3 | 3D vector algebra | lidar_vec3_add/sub/dot/cross/norm/normalize |
| 4 | 3x3/4x4 matrix operations | lidar_mat3_mul/det/inverse/transpose, lidar_mat4_transform |
| 5 | Rotation representations | Rodrigues formula (lidar_rotation_axis_angle), Euler angles (lidar_rotation_euler) |
| 6 | Singular Value Decomposition | Jacobi SVD in svd_3x3() for Procrustes problem |
| 7 | Covariance matrix (PCA) | lidar_estimate_normal_knn() — 3x3 covariance + eigenvalue decomposition |
| 8 | Poisson statistics | lidar_poisson_prob(), poisson_pmf in Lean |
| 9 | Gaussian beam propagation | lidar_beam_radius(), lidar_beam_intensity() |
| 10 | Beer-Lambert law | lidar_atmosphere_transmission(), eer_lambert_transmission |
| 11 | Angstrom exponent (wavelength scaling) | lidar_atmosphere_extinction(), ngstrom_exponent |
| 12 | Q-function (Gaussian tail) | Q_function(), Q_inverse() — rational approximations |
| 13 | Non-linear least squares | Levenberg-Marquardt in lidar_levenberg_marquardt_gaussian() |

## L4: Fundamental Laws (Complete ✅)

| # | Law / Theorem | C Verification | Lean Formalization |
|---|--------------|---------------|-------------------|
| 1 | LiDAR range equation | 	est_range_equation — received power computed and validated | eceived_power_simple |
| 2 | TOF ↔ Range identity | 	est_tof_range — round-trip verified | 	of_range_inverse (by rfl) |
| 3 | Range resolution law: ΔR = c·τ/2 | 	est_range_equation — ΔR validated | ange_resolution_pulsed |
| 4 | Unambiguous range: R_unamb = c/(2·PRF) | Implicit in config validation | unambiguous_range |
| 5 | Beer-Lambert transmission law | 	est_atmosphere — extinction + transmission | eer_lambert_transmission |
| 6 | Koschmieder visibility law | 	est_atmosphere — extinction computed | koschmieder_beta |
| 7 | Poisson photon statistics | 	est_poisson — P(k\|λ) verified | poisson_pmf |
| 8 | Geiger detection probability | 	est_performance — Pd computed | geiger_detection_prob |
| 9 | Diffraction-limited resolution: θ = 1.22λ/D | 	est_diffraction | — |
| 10 | Neyman-Pearson detection criterion | lidar_detection_threshold_gaussian(), 	est_detection_theory | — |
| 11 | Cramer-Rao Lower Bound for range | lidar_range_crlb() | — |

## L5: Algorithms / Methods (Complete ✅)

| # | Algorithm | Implementation | Complexity |
|---|----------|---------------|------------|
| 1 | Constant Fraction Discrimination | lidar_cfd_timing() | O(1) |
| 2 | Leading-edge timing | lidar_leading_edge_timing() | O(N) |
| 3 | Derivative peak detection | lidar_detect_peaks_derivative() | O(N) |
| 4 | Levenberg-Marquardt nonlinear optimization | lidar_levenberg_marquardt_gaussian() | O(N·M·I) |
| 5 | Gaussian waveform decomposition | lidar_gaussian_decompose() | O(N·M·I) |
| 6 | PCA normal estimation (Jacobi) | lidar_estimate_normal_knn() | O(K) + 20 Jacobi iterations |
| 7 | RANSAC plane fitting | lidar_ransac_plane() | O(max_iter · N) |
| 8 | Voxel grid downsampling | lidar_voxel_downsample() | O(N) |
| 9 | SVD rigid transform (Arun et al.) | lidar_estimate_rigid_transform() | O(N) |
| 10 | Iterative Closest Point (ICP) | lidar_icp_point_to_point() | O(iter · N_src · N_tgt) |
| 11 | CA-CFAR detection | lidar_cfar_detection() | O(N_bins · N_ref) |
| 12 | Matched filter (time-domain) | lidar_matched_filter() | O(N·M) |
| 13 | DEM generation (rasterization) | lidar_dem_generate() | O(N_points) |
| 14 | DSM generation (first-return max) | lidar_dsm_generate() | O(N_points) |
| 15 | Hillshade computation | lidar_dem_hillshade() | O(rows · cols) |
| 16 | Building extraction (clustering) | lidar_extract_buildings() | O(N²) BFS |
| 17 | Forestry metrics (percentiles, LAI, biomass) | lidar_forestry_metrics() | O(N log N) |
| 18 | Object detection (Euclidean clustering) | lidar_detect_objects_euclidean() | O(N²) BFS |
| 19 | Aerosol backscatter (Klett inversion) | lidar_aerosol_backscatter() | O(N_bins) |
| 20 | Boundary layer height detection | lidar_boundary_layer_height() | O(N_bins) |

## L6: Canonical Problems (Complete ✅)

| # | Problem | Solution |
|---|---------|---------|
| 1 | Range estimation from TOF | lidar_tof_to_range() + test |
| 2 | Point cloud alignment (ICP) | lidar_icp_point_to_point() + test |
| 3 | Ground/non-ground classification | lidar_ground_filter_height(), lidar_ground_filter_pmf() |
| 4 | Digital Elevation Model generation | lidar_dem_generate(), lidar_dsm_generate() |
| 5 | Building extraction from airborne LiDAR | lidar_extract_buildings() |
| 6 | Full-waveform decomposition | lidar_gaussian_decompose() |
| 7 | Pulse detection with CFAR | lidar_cfar_detection() + matched filter |
| 8 | Surface normal estimation | lidar_estimate_normal_knn() |
| 9 | Plane detection via RANSAC | lidar_ransac_plane() |

## L7: Applications (Complete — 6 applications)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | 🚗 Autonomous vehicle perception | lidar_detect_objects_euclidean() + example_object_detection.c |
| 2 | 🏔️ Terrain mapping / DEM generation | lidar_dem_generate() + example_dem_generation.c |
| 3 | 🌲 Forestry (canopy, LAI, biomass) | lidar_forestry_metrics() |
| 4 | 🏢 Urban 3D modeling / building extraction | lidar_extract_buildings() |
| 5 | 🌫️ Atmospheric sensing (aerosol, PBL) | lidar_aerosol_backscatter(), lidar_boundary_layer_height() |
| 6 | 📡 Full-waveform LiDAR analysis | example_waveform_decomp.c |

## L8: Advanced Topics (Complete — 3 topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Single-photon (Geiger-mode) LiDAR | geiger_detection_prob, geiger_false_alarm_prob, SPAD/SiPM detector types |
| 2 | FMCW LiDAR with chirp processing | lidar_fmcw_beat_frequency(), lidar_fmcw_range_from_beat() |
| 3 | Multi-spectral LiDAR (wavelength-dependent extinction) | lidar_atmosphere_extinction() with Angstrom exponent, all wavelength bands |

## L9: Research Frontiers (Partial — documented)

| # | Frontier | Status |
|---|---------|--------|
| 1 | Quantum LiDAR (entangled photons for super-resolution) | Documented — no implementation |
| 2 | Chip-scale LiDAR (silicon photonics OPA) | Documented — OPA scan type defined |
| 3 | 4D LiDAR (range + velocity per pixel via FMCW) | Partially implemented — lidar_fmcw_beat_frequency() with Doppler term |
