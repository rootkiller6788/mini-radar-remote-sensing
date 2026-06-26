# Gap Report — mini-pulse-doppler

## Missing Items by Priority

### High Priority (L1-L6 gaps)
- **None** — All L1-L6 items are covered with implementations.

### Medium Priority (L7 gaps)
| # | Missing Application | Reason | Priority |
|---|--------------------|--------|-----------|
| 1 | SAR (Synthetic Aperture Radar) imaging | Requires 2-D FFT and complex image formation | Medium |
| 2 | Automotive FMCW radar processing | Range-FFT + Doppler-FFT dual processing | Medium |
| 3 | Multi-target tracking (Kalman filter) | Beyond scope of pulse Doppler module | Low |

### Low Priority (L8 gaps)
| # | Missing Advanced Topic | Reason | Priority |
|---|-----------------------|--------|-----------|
| 1 | Full STAP (Space-Time Adaptive Processing) | Requires covariance matrix estimation | Low |
| 2 | MIMO radar beamforming | Requires multi-channel architecture | Low |
| 3 | Track-before-detect (TBD) | Requires particle filters or Hough transform | Low |
| 4 | Compressive sensing radar | Requires sparse recovery algorithms | Low |

### Research (L9)
- MIMO radar: documented only
- Compressive sensing: documented only
- Quantum radar: mentioned only
- AI-based target classification: mentioned only

## Next Steps
1. Add SAR imaging example (range migration algorithm)
2. Add FMCW radar processing (dual-FFT)
3. Implement full STAP with sample matrix inversion
4. Add MIMO beamforming with virtual array concept
