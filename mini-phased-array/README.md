# mini-phased-array — Phased Array Beamforming and Radar Library

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Complete (3+ applications: AESA T/R module, Radar Range Equation, G/T FOM)
- L8: Partial (MVDR/Capon beamforming, LCMV structural definition, SMI covariance)
- L9: Partial (documented: metasurface arrays, photonic TTD, AI beamforming)

**Line Count**: 3,719 lines (include/ + src/) ≥ 3,000 ✅  
**make test**: 31/31 tests passed ✅  
**make examples**: All 3 examples build ✅

---

## Knowledge Coverage Summary

| Level | Name | Coverage | Key Items |
|-------|------|----------|-----------|
| **L1** | Definitions | **Complete** | Array geometries (6 types), element types (8), T/R module, steering state, beam pattern result, grating lobe descriptor, SMI/MVDR/LCMV beamformer states |
| **L2** | Core Concepts | **Complete** | Element pattern × AF = total pattern, beamforming architectures (analog/digital/hybrid), TTD vs phase steering, beam squint, mutual coupling, monopulse angle estimation |
| **L3** | Math Structures | **Complete** | Spherical ↔ Az/El ↔ unit vector transforms, complex weight vectors, covariance matrices (Hermitian), steering vectors, wavenumber/spatial frequency |
| **L4** | Fundamental Laws | **Complete** | Pattern Multiplication Theorem, Nyquist spatial sampling (d < λ), Radar Range Equation (phased array form), Friis transmission, aperture directivity relation |
| **L5** | Algorithms/Methods | **Complete** | Steering vector computation, Dolph-Chebyshev synthesis, Taylor one-parameter distribution, Binomial array, amplitude taper windows (9 types), phase quantization, TTD delay computation, mutual coupling S-matrix, Gauss-Jordan matrix inversion |
| **L6** | Canonical Problems | **Complete** | Linear array steering, planar array 2D steering, grating lobe detection, sidelobe suppression, monopulse Σ/Δ, beam squint analysis, 3dB beamwidth |
| **L7** | Applications | **Complete** | AESA T/R module (F-35 AN/APG-81 class), Radar range equation, G/T figure of merit, NEDT radiometry, effective aperture |
| **L8** | Advanced Topics | **Partial** | MVDR/Capon adaptive beamforming, SMI covariance estimation, LCMV beamformer structure, diagonal loading regularization |
| **L9** | Research Frontiers | **Partial** | Metasurface phased arrays, Photonic TTD beamforming, AI-based beamforming (documented only) |

---

## Core Definitions (L1)

- `pa_geometry_t`: 6 array geometries (linear, planar rect/triang, circular, cylindrical, spherical, hexagonal)
- `pa_element_type_t`: 8 antenna element types (isotropic, dipole, patch rect/circular, slot, Vivaldi, crossed dipole, dual-feed)
- `pa_window_t`: 10 amplitude taper windows (uniform through custom, incl. Dolph-Chebyshev, Taylor, Binomial)
- `pa_beamforming_arch_t`: 5 beamforming architectures (analog RF/IF, digital element/subarray, optical TTD)
- `pa_array_config_t`: Complete array geometry descriptor
- `pa_element_t`: Single array element with position, excitation, T/R module data
- `pa_steering_state_t`: Runtime beam steering state
- `pa_pattern_t`: Full 2D beam pattern with metrics
- `pa_tr_module_t`: AESA T/R module descriptor
- `pa_smi_beamformer_t`, `pa_mvdr_beamformer_t`, `pa_lcmv_beamformer_t`: Adaptive beamformer states

---

## Core Theorems (L4)

### Pattern Multiplication Theorem
```
E_total(θ,φ) = E_element(θ,φ) × AF(θ,φ)
```
The total far-field radiation pattern of an array is the product of the single-element pattern and the array factor. (Balanis §6.3)

### Nyquist Spatial Sampling (Grating Lobe Avoidance)
```
d < λ / (1 + |sin θ_s_max|)
```
To avoid grating lobes in visible space for maximum scan angle θ_s_max. For broadside: d < λ. (Balanis §6.9)

### Radar Range Equation (Phased Array)
```
R⁴ = P_t · G² · λ² · σ / ((4π)³ · k · T₀ · B · F · SNR_min · L)
```
Detection range for monostatic AESA radar. (Richards et al. §2.2)

### Dolph-Chebyshev Mapping
```
AF(ψ) = T_{N-1}(x₀ · cos(ψ/2)),  x₀ = cosh(acosh(R₀)/(N-1))
```
Optimal beamwidth for a given sidelobe level using Chebyshev polynomials. (Dolph, 1946)

---

## Core Algorithms (L5)

- **Dolph-Chebyshev synthesis**: Constant sidelobe level, minimal beamwidth (Balanis §6.8.3)
- **Taylor one-parameter distribution**: Near-in constant SLL, far-out monotonic decay (Taylor, 1955)
- **Steering vector**: Phase-based beam steering (Van Trees §2.4)
- **Phase quantization**: Finite-bit phase shifter effects (Mailloux §3.4)
- **Gauss-Jordan inversion**: O(M³) complex matrix inverse for covariance matrices
- **MVDR/Capon beamformer**: Adaptive interference nulling (Capon, 1969)
- **SMI covariance estimation**: Sample matrix inversion with diagonal loading

---

## Canonical Problems (L6)

1. **Linear array beam steering** — 16-element ULA at X-band, steer to 30°
2. **Planar array 2D steering** — 8×8 AESA, azimuth/elevation scanning
3. **Grating lobe analysis** — Detect visible grating lobes for a given scan
4. **Sidelobe suppression** — Compare uniform vs Chebyshev vs Taylor tapers
5. **Monopulse angle estimation** — Sum/difference patterns for tracking radar
6. **Beam squint analysis** — Frequency-dependent beam pointing error

---

## 九校课程映射 (Nine-School Course Alignment)

| School | Course | Topics Covered |
|--------|--------|----------------|
| **MIT** | 6.630 EM Waves | Array factor, pattern multiplication, grating lobes |
| **Stanford** | EE359 Wireless | Beamforming, MIMO arrays, adaptive nulling |
| **Berkeley** | EE117 EM | Antenna arrays, mutual coupling, phased arrays |
| **Michigan** | EECS 411 Microwave | Phased array design, T/R modules, AESA |
| **Georgia Tech** | ECE 6350 EM | Array synthesis, Dolph-Chebyshev, Taylor |
| **TU Munich** | High-Frequency Eng | Beam steering, phase shifters, TTD |
| **ETH** | 227-0455 EM | Array theory, numerical pattern computation |
| **Tsinghua** | 电磁场 | Array antenna fundamentals, phased array radar |

---

## Build and Test

```bash
make          # Build library and test binary
make test     # Run 31 tests (all pass)
make examples # Build 3 end-to-end examples
make clean    # Remove artifacts
```

---

## References

1. Balanis, C.A. (2016) *Antenna Theory: Analysis and Design*, 4th ed. Wiley.
2. Mailloux, R.J. (2005) *Phased Array Antenna Handbook*, 2nd ed. Artech House.
3. Van Trees, H.L. (2002) *Optimum Array Processing*, Part IV. Wiley.
4. Richards, M.A., Scheer, J.A., Holm, W.A. (2010) *Principles of Modern Radar*. SciTech.
5. Skolnik, M.I. (2008) *Radar Handbook*, 3rd ed. McGraw-Hill.
6. Dolph, C.L. (1946) "A Current Distribution for Broadside Arrays...", Proc. IRE.
7. Taylor, T.T. (1955) "Design of Line-Source Antennas...", IRE Trans. AP-3.
8. Capon, J. (1969) "High-Resolution Frequency-Wavenumber Spectrum Analysis", Proc. IEEE.
9. Sherman & Barton (2011) *Monopulse Principles and Techniques*, 2nd ed. Artech House.
10. Brookner, E. (2008) "Phased Arrays Around the World", Microwave Journal.
