# Knowledge Graph — mini-radar-basics

## L1: Definitions (Complete)
- Radar range equation (monostatic, bistatic)
- Radar cross section (RCS), Swerling fluctuation models I-IV
- Pulse width (tau), PRI, PRF, duty cycle
- Range resolution: Delta_R = c*tau/2 or c/(2*B)
- Unambiguous range: R_unamb = c/(2*PRF)
- Unambiguous velocity: v_unamb = lambda*PRF/4
- Doppler shift: f_d = 2*v_r/lambda
- Antenna gain G, effective aperture A_e, beamwidth
- Peak power P_t, average power P_avg
- SNR, noise figure F, noise temperature T_sys
- Radar wavelength lambda, operating band designations
- Waveform types: simple pulse, LFM, Barker, Frank, stepped-freq
- PSL (peak sidelobe level), ISL (integrated sidelobe level)
- CFAR window: CUT, guard cells, reference cells

## L2: Core Concepts (Complete)
- Time-of-flight ranging principle
- Pulse radar vs CW vs FMCW operation
- Power-aperture product
- System noise model: T_sys = T_a + T0*(F-1)
- Radar loss budget
- Matched filter maximizes output SNR
- Pulse compression: long pulse energy + short pulse resolution
- LFM chirp: linear frequency sweep, large TBP
- CFAR principle: adaptive threshold from local noise
- Doppler effect: frequency shift proportional to radial velocity
- Coherent pulse train: phase progression across pulses
- Clutter as distributed (area-extensive) target

## L3: Mathematical Structures (Complete)
- Radar range equation as Friis-like power budget
- Complex baseband / complex envelope representation
- Point target scattering model
- Free-space path loss: L_fs = (4*pi*R/lambda)^2
- Matched filter impulse response: h(t) = s*(T-t)
- Matched filter output = cross-correlation
- Gaussian noise: n ~ N(0, sigma^2)
- Chi-square / gamma for square-law detector
- Marcum Q-function for non-fluctuating targets
- Order statistics for OS-CFAR
- Phase progression: phi[n] = 2*pi*f_d*n*PRI
- Doppler steering vector

## L4: Fundamental Laws (Complete)
- Radar range equation: P_r = P_t*G^2*lambda^2*sigma/((4*pi)^3*R^4*L)
- Maximum detection range (fourth-root dependence)
- Bistatic radar equation
- Unambiguous range/velocity constraints
- Neyman-Pearson lemma (optimal detector for AWGN)
- Albersheim's equation (P_d approximation)
- Shnidman's equation (extended for Swerling targets)
- Nyquist sampling in slow time (PRF >= 2*f_d_max)
- Antenna reciprocity: G = 4*pi*A_e/lambda^2

## L5: Algorithms/Methods (Complete)
- CA-CFAR (Cell-Averaging)
- OS-CFAR (Order-Statistic)
- GO-CFAR, SO-CFAR
- LFM chirp generation: s[n] = exp(j*pi*k*(n/f_s)^2)
- Barker code generation (lengths 2,3,4,5,7,11,13)
- Frank polyphase code generation (N_c x N_c)
- Matched filter: time-domain correlation O(N*M)
- Matched filter: FFT-based fast convolution O(N log N)
- MTI filter design: single/double/triple canceller
- Doppler FFT: 2D range-Doppler map processing
- ULA beamforming: array factor, steering vector
- PSL/ISL computation for sidelobe analysis
- Waveform autocorrelation
- Ambiguity function sampling
- Velocity from phase difference
- Doppler ambiguity resolution (staggered PRF)

## L6: Canonical Problems (Complete)
- Pulse radar range estimation from SNR
- FMCW radar range-velocity joint estimation
- CFAR target detection in homogeneous noise
- CFAR detection with interfering targets
- Multi-target radar return generation
- Clutter + noise simulation
- Range profile computation (matched filter output)
- Range-Doppler map from coherent pulse train
- ROC curve generation (P_d vs P_fa)
- MTI clutter rejection

## L7: Applications (Complete — 10+ items)
- Automotive radar (77 GHz FMCW)
- Weather radar (rain reflectivity)
- Air traffic control (ATC) radar
- Cognitive radar adaptation
- Radar band selection (HF through mm-wave)
- Chaff cloud RCS modeling
- Sea clutter (Nathanson model)
- Ground clutter (constant-gamma model)
- Sensitivity time control (STC)
- Micro-Doppler vibration analysis
- Pulse-Doppler visibility
- Binary integration (M-of-N detection)
- Detection range at specified Pd/Pfa

## L8: Advanced Topics (Partial — 5 items)
- Nonlinear FM (NLFM) waveforms
- MIMO virtual array concept
- MVDR (Capon) adaptive beamforming
- MUSIC DOA estimation
- STAP (Space-Time Adaptive Processing) sampling

## L9: Research Frontiers (Partial)
- Cognitive radar perception-action cycle
- 6G joint communication and sensing (documented)
- Quantum radar concepts (documented)
