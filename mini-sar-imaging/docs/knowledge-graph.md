# Knowledge Graph -- mini-sar-imaging

## L1: Definitions (COMPLETE)
- SAR parameters: carrier frequency f0, bandwidth B_r, pulse width tau_p, PRF
- Range resolution rho_r = c/(2*B_r)
- Azimuth resolution rho_a = L_a/2 (focused stripmap)
- Chirp/Linear FM waveform: s(t) = exp(j*pi*K_r*t^2)
- Pulse Repetition Frequency (PRF), pulse repetition interval (PRI)
- Range Cell Migration (RCM): range walk + range curvature
- Complex SAR image (SLC): I+jQ, magnitude, phase
- Backscatter coefficient sigma_0 (NRCS)
- Slant range vs ground range
- Doppler centroid f_Dc, Doppler rate f_R
- Synthetic aperture time T_sa, integration angle
- Radar band designations (P, L, S, C, X, Ku, K, Ka, W)
- SAR imaging modes: stripmap, spotlight, ScanSAR, TOPS, ISAR, bistatic
- Interferometric baseline: parallel B_||, perpendicular B_perp
- Interferometric coherence gamma
- Wrapped vs unwrapped phase
- Height ambiguity h_a
- Polarimetric channels: HH, HV, VH, VV
- Coherency matrix T (3x3)

## L2: Core Concepts (COMPLETE)
- Synthetic aperture principle: motion -> virtual array -> fine azimuth resolution
- Pulse compression / matched filtering: max SNR for known signal
- Range-Doppler coupling in chirp signals
- Range Cell Migration Correction (RCMC) principle
- Hyperbolic range equation as SAR phase history
- Speckle: multiplicative noise from coherent imaging
- Multi-looking: incoherent averaging for speckle reduction
- Interferometric SAR (InSAR): phase difference -> topography/deformation
- Phase unwrapping necessity and algorithms
- Antenna beam pattern and azimuth illumination
- Coherent vs incoherent integration
- Earth curvature effects on SAR geometry (effective velocity)

## L3: Mathematical Structures (COMPLETE)
- 2D Fourier transform for frequency-domain SAR processing
- Range-Doppler domain: (tau -> f_tau, eta -> f_eta)
- Wavenumber domain (omega-k): k_r, k_x
- Stationary phase method for SAR signal analysis
- Stolt interpolation: (f_tau, f_eta) -> (k_r, k_x)
- Taylor expansion of hyperbolic range equation
- Chu's relationship: space-Doppler domain mapping
- Hermitian coherency matrix eigenvalue decomposition
- sinc function and its properties for impulse response

## L4: Fundamental Laws (COMPLETE)
- Range resolution: rho_r = c/(2*B_r) -- verified
- Azimuth resolution: rho_a = L_a/2 (focused stripmap) -- verified
- PRF Nyquist constraint: PRF >= 2*v/L_a -- verified
- Range ambiguity constraint: PRF <= c/(2*swath) -- verified
- Doppler centroid: f_Dc = 2*v*sin(theta_sq)/lambda -- verified
- Doppler rate: f_R = -2*v^2*cos^3(theta_sq)/(lambda*R0) -- verified
- InSAR height sensitivity: delta_phi = (4*pi*B_perp*h)/(lambda*R0*sin(theta)) -- verified
- DInSAR displacement: delta_d = lambda*delta_phi_diff/(4*pi) -- verified

## L5: Algorithms/Methods (COMPLETE)
- Range-Doppler Algorithm (RDA): range comp -> azimuth FFT -> RCMC -> azimuth comp -> IFFT
- Chirp Scaling Algorithm (CSA): CS phase -> bulk RCMC -> azimuth comp
- omega-k Algorithm: 2D FFT -> ref multiply -> Stolt -> 2D IFFT
- Backprojection (BP): time-domain coherent summation
- SPECAN: azimuth deramp + FFT for burst-mode
- Phase Gradient Autofocus (PGA): iterative non-parametric
- Map Drift (MD) autofocus: sub-aperture cross-correlation
- Phase unwrapping: Goldstein branch-cut, quality-guided
- Goldstein interferogram filter: adaptive spectral filtering
- Coherence estimation: spatial averaging
- ISTA for compressive sensing SAR

## L6: Canonical Problems (COMPLETE)
- Point target simulation and focusing (example_point_target.c)
- Stripmap SAR image formation (example_stripmap.c)
- RCMC: range cell migration correction (sar_rda_rcmc)
- Impulse response analysis: resolution, PSLR, ISLR
- Multi-looking for speckle reduction (sar_multilook)

## L7: Applications (COMPLETE)
- InSAR DEM generation (example_insar.c)
- DInSAR displacement measurement
- Polarimetric SAR classification (Freeman-Durden, H-Alpha)

## L8: Advanced Topics (PARTIAL)
- Compressive Sensing SAR: ISTA reconstruction -- IMPLEMENTED
- MIMO-SAR: virtual array processing -- IMPLEMENTED
- Bistatic SAR: bistatic range equation + processor -- IMPLEMENTED
- Polarimetric decomposition: Freeman-Durden, H-Alpha -- IMPLEMENTED
- SAR tomography (3D) -- DOCUMENTED ONLY

## L9: Research Frontiers (PARTIAL)
- AI/ML-based SAR imaging and classification -- DOCUMENTED ONLY
- 3D SAR tomography -- DOCUMENTED ONLY
- Quantum SAR concepts -- DOCUMENTED ONLY
- Terahertz SAR -- DOCUMENTED ONLY