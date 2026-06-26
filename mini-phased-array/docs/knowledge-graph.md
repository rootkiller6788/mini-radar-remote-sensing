# Knowledge Graph — mini-phased-array

## L1: Definitions (Complete)

| ID | Concept | C Representation | Lean Representation |
|----|---------|-----------------|---------------------|
| L1.1 | Array Geometry Types | `pa_geometry_t` enum (7 values) | `ArrayConfig.geometryTag : Nat` |
| L1.2 | Antenna Element Types | `pa_element_type_t` enum (8 values) | — |
| L1.3 | Array Configuration | `pa_array_config_t` struct | `ArrayConfig` structure |
| L1.4 | Element Descriptor | `pa_element_t` struct | `ArrayElement` structure |
| L1.5 | 3D Position Vector | `pa_position3d_t` struct | `Position3D` structure |
| L1.6 | Spherical Coordinates | `pa_spherical_t` struct | `SteeringDirection` structure |
| L1.7 | Az/El Coordinates | `pa_azel_t` struct | — |
| L1.8 | Array Factor Result | `pa_af_result_t` struct | `BeamResult` structure |
| L1.9 | Beam Pattern | `pa_pattern_t` struct | — |
| L1.10 | Steering State | `pa_steering_state_t` struct | — |
| L1.11 | T/R Module | `pa_tr_module_t` struct | `TRModule` structure |
| L1.12 | Grating Lobe Descriptor | `pa_grating_lobe_t` struct | — |
| L1.13 | SMI Beamformer | `pa_smi_beamformer_t` struct | — |
| L1.14 | MVDR Beamformer | `pa_mvdr_beamformer_t` struct | `AdaptiveBeamformer` structure |
| L1.15 | Beam Squint Result | `pa_beam_squint_t` struct | — |
| L1.16 | TTD Element | `pa_ttd_element_t` struct | — |
| L1.17 | Amplitude Window Types | `pa_window_t` enum (10 values) | — |
| L1.18 | Beamforming Architecture | `pa_beamforming_arch_t` enum (5 values) | — |
| L1.19 | Steering Mechanism | `pa_steering_t` enum (3 values) | — |
| L1.20 | Phase Quantization Bits | `pa_phase_bits_t` enum | — |

## L2: Core Concepts (Complete)

| ID | Concept | Implementation |
|----|---------|---------------|
| L2.1 | Pattern Multiplication | `pa_total_pattern()` × `pa_array_factor()` |
| L2.2 | Element Pattern × AF Separation | `pa_element_pattern()` + `pa_array_factor()` |
| L2.3 | Phase Steering | `pa_phase_shifts()`, `pa_steering_vector()` |
| L2.4 | True Time Delay | `pa_ttd_compute_delays()` |
| L2.5 | Beam Squint | `pa_beam_squint_analysis()` |
| L2.6 | Grating Lobes (Spatial Aliasing) | `pa_find_grating_lobes()` |
| L2.7 | Mutual Coupling | `pa_compute_coupling_matrix()`, `pa_mutual_coupling_correct()` |
| L2.8 | Monopulse Angle Estimation | `pa_monopulse_patterns()`, `pa_monopulse_estimate_angle()` |
| L2.9 | Adaptive Nulling | `pa_mvdr_compute_weights()` |
| L2.10 | Digital Beamforming | `pa_beamforming_arch_t` — architectural enumeration |

## L3: Mathematical Structures (Complete)

| ID | Structure | Implementation |
|----|-----------|---------------|
| L3.1 | Coordinate Transforms | `pa_azel_to_spherical()`, `pa_spherical_to_azel()`, `pa_spherical_to_uvw()`, `pa_uvw_to_spherical()` |
| L3.2 | Complex Vector Spaces | `double complex *` — steering vectors, weight vectors |
| L3.3 | Hermitian Matrices | `double complex *` — covariance matrices (M×M row-major) |
| L3.4 | Wavenumber (Spatial Frequency) | `pa_wavenumber()` |
| L3.5 | Array Factor Summation | Σ w_n exp(j k₀ r_n · û) — in `pa_array_factor()` |
| L3.6 | Matrix-Vector Multiplication | `mat_vec_mul()` (adaptive.c) |
| L3.7 | Complex Inner Product | `inner_product()` (adaptive.c) |
| L3.8 | Matrix Inversion (Gauss-Jordan) | `matrix_invert()` (adaptive.c), `matrix_invert_gj()` (aesa.c) |

## L4: Fundamental Laws (Complete)

| ID | Law/Theorem | Implementation |
|----|------------|---------------|
| L4.1 | Pattern Multiplication Theorem | `pa_total_pattern()` = element × AF |
| L4.2 | Nyquist Spatial Sampling | d < λ condition in `pa_find_grating_lobes()` |
| L4.3 | Radar Range Equation | `pa_radar_range_equation()` — full 4th-power law |
| L4.4 | Aperture-Gain Relation | `pa_effective_aperture()`, `pa_directivity()` |
| L4.5 | Dolph-Chebyshev Optimality | `pa_dolph_chebyshev_weights()` — Balanis §6.8.3 |
| L4.6 | Taylor Distribution | `pa_taylor_weights()` — Taylor (1955) |
| L4.7 | MVDR Optimality | `pa_mvdr_compute_weights()` — Capon (1969) |

## L5: Algorithms/Methods (Complete)

| ID | Algorithm | Implementation |
|----|-----------|---------------|
| L5.1 | Steering Vector Computation | `pa_steering_vector()` |
| L5.2 | Phase Shift Calculation | `pa_phase_shifts()` |
| L5.3 | Phase Quantization | `pa_quantize_phase_shifts()` |
| L5.4 | Dolph-Chebyshev Synthesis | `pa_dolph_chebyshev_weights()` |
| L5.5 | Taylor Distribution | `pa_taylor_weights()` |
| L5.6 | Binomial Array | `pa_binomial_weights()` |
| L5.7 | Window Functions (9 types) | `pa_amplitude_taper()` |
| L5.8 | TTD Delay Computation | `pa_ttd_compute_delays()` |
| L5.9 | Coupling Matrix Generation | `pa_compute_coupling_matrix()` |
| L5.10 | SMI Covariance Estimation | `pa_smi_estimate_covariance()` |
| L5.11 | Gauss-Jordan Matrix Inversion | `matrix_invert()` — O(M³) |
| L5.12 | Diagonal Loading | Applied in `pa_smi_estimate_covariance()` |

## L6: Canonical Problems (Complete)

| ID | Problem | Implementation |
|----|---------|---------------|
| L6.1 | Linear Array Steering | `example_linear_array.c` |
| L6.2 | Planar Array 2D Steering | `example_planar_array.c` |
| L6.3 | Adaptive Interference Rejection | `example_adaptive_beamforming.c` |
| L6.4 | Grating Lobe Detection | `pa_find_grating_lobes()` |
| L6.5 | Beam Pattern Computation | `pa_compute_pattern()` |
| L6.6 | Sidelobe Analysis | Pattern metrics in `pa_compute_pattern()` |
| L6.7 | Monopulse Σ/Δ | `pa_monopulse_patterns()` |
| L6.8 | Beam Squint Analysis | `pa_beam_squint_analysis()` |

## L7: Applications (Complete — 3+ applications)

| ID | Application | Implementation |
|----|------------|---------------|
| L7.1 | AESA T/R Module (F-35 AN/APG-81 class) | `pa_tr_module_init()` |
| L7.2 | Radar Range Equation | `pa_radar_range_equation()` |
| L7.3 | G/T Figure of Merit | `pa_array_gt_fom()` |
| L7.4 | NEDT Radiometry | `pa_array_nedt()` |
| L7.5 | Effective Aperture | `pa_effective_aperture()` |

## L8: Advanced Topics (Partial — 2 implementations)

| ID | Topic | Implementation |
|----|-------|---------------|
| L8.1 | MVDR/Capon Beamforming | `pa_mvdr_compute_weights()`, `pa_mvdr_capon_spectrum()` |
| L8.2 | SMI Covariance Estimation | `pa_smi_estimate_covariance()` |
| L8.3 | LCMV Beamformer (structural) | `pa_lcmv_beamformer_t` + internal computation (aesa.c) |
| L8.4 | Diagonal Loading | Parameterized regularization |
| L8.5 | SINR Optimization | `pa_mvdr_output_sinr()` |

## L9: Research Frontiers (Partial — documented only)

| ID | Topic | Status |
|----|-------|--------|
| L9.1 | Metasurface Phased Arrays | Documented in Lean file |
| L9.2 | Photonic TTD Beamforming | Documented in Lean file |
| L9.3 | AI-Based Beamforming | Documented in Lean file |
