# Course Alignment -- mini-sar-imaging

## MIT
- 6.003 Signal Processing: Fourier transforms, convolution (sar_fft2d, sar_pulse_compression)
- 6.630 Electromagnetic Waves: radar equation, RCS, antenna pattern (sar_params_t, sar_antenna_t)

## Stanford
- EE359 Wireless Communications: matched filter, pulse compression (sar_pulse_compression)
- EE368 Radar Remote Sensing: SAR imaging modes, Doppler processing (sar_algorithm.c)

## Berkeley
- EE117 Electromagnetic Fields & Waves: plane wave propagation, RCS (sar_params_init)
- EE123 Digital Signal Processing: FFT, filter design (sar_fft2d, sar_window_*)

## Illinois
- ECE 451 EM Fields: radar cross section, antenna theory (sar_calibrate_sigma0)
- ECE 459 Communications: matched filter, signal detection (sar_matched_filter_coeff)

## Michigan
- EECS 411 Microwave Circuits: radar systems, antenna (sar_antenna_gain)
- EECS 455 Wireless Communications: DSP for radar, SAR processing (sar_rda_process)

## Georgia Tech
- ECE 6350 Applied EM: SAR principles, propagation (sar_range_hyperbolic)
- ECE 6601 Digital Communications: detection theory, matched filter (sar_pulse_compression_fft)

## TU Munich
- Signal Processing: Fourier methods, spectral analysis (sar_spectrum_analyze)
- High-Frequency Engineering: microwave radar, antenna design (sar_antenna_init)
- Communications: matched filtering, optimal detection (sar_matched_filter_coeff)

## ETH Zurich
- 227-0455 EM Waves: radar principles, scattering (sar_calibrate_sigma0)
- 227-0436 Digital Communications: signal processing for radar (sar_csa_process)
- 227-0427 Signal Processing: FFT, convolution, filtering (sar_fft2d)

## Tsinghua
- 信号与系统 (Signals & Systems): Fourier transform, convolution
- 通信原理 (Communication Principles): matched filter, signal detection
- 电磁场 (Electromagnetic Fields): antenna, propagation, RCS
- 雷达信号处理 (Radar Signal Processing): SAR imaging algorithms (full module)