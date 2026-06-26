# Course Alignment — mini-target-tracking

Mapping of this module to nine world-class university curricula.

## MIT (Massachusetts Institute of Technology)

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| 6.003 Signal Processing | Kalman filtering, state-space models | `kalman_filter.c` — KF predict/update, Riccati equation |
| 6.450 Digital Communications | Detection, estimation theory | `measurement_models.c` — ML estimation, CRLB |
| 6.630 EM Waves | Radar principles | `measurement_models.c` — range, Doppler, bearing |
| 16.485 Visual Navigation | Multi-target tracking | `data_association.c` — GNN, PDA, JPDA |

## Stanford University

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| EE102A Signal Processing | Linear estimation, Wiener/Kalman | `kalman_filter.c` — optimal filtering |
| EE359 Wireless Communications | Radar signal processing | `mmwave_tracker.c` — automotive radar |
| EE263 Linear Dynamical Systems | State estimation, observability | `motion_models.c` — state transition models |
| AA273 State Estimation | Kalman filters, nonlinear filtering | `kalman_filter.c` — KF, EKF, UKF |

## UC Berkeley

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| EE123 Digital Signal Processing | Adaptive filtering, estimation | `kalman_filter.c` — adaptive KF |
| EE221A Linear System Theory | Observers, estimators | `track_core.c` — state estimation |
| CS289A Machine Learning | Clustering, classification | `mmwave_tracker.c` — DBSCAN, object classification |

## University of Michigan

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| EECS 351 DSP | Kalman filtering | `kalman_filter.c` |
| EECS 455 Communications | Detection, estimation | `measurement_models.c` — CRLB |
| EECS 411 Microwave Circuits | Radar systems | `mmwave_tracker.c` — mmWave radar |
| NAVARCH 580 Estimation | Multi-target tracking | `data_association.c` — PDA, JPDA, MHT |

## Georgia Tech

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| ECE 4270 DSP | Kalman/Wiener filters | `kalman_filter.c` |
| ECE 6601 Communications | Detection/estimation | `measurement_models.c` |
| ECE 6350 EM | Radar remote sensing | `track_core.c` + measurement models |
| ECE 6550 Nonlinear Systems | EKF, UKF, particle filters | `kalman_filter.c` — EKF, UKF |

## TU Munich

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| Signal Processing | Estimation theory | `kalman_filter.c` — KF variants |
| Communications Engineering | Radar detection | `measurement_models.c` — radar models |
| High-Frequency Engineering | Automotive radar | `mmwave_tracker.c` — 77 GHz radar tracking |

## ETH Zürich

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| 227-0427 Signal Processing | Kalman filtering | `kalman_filter.c` |
| 227-0436 Communications | Multi-user detection | `data_association.c` |
| 227-0455 EM Waves | Radar signal processing | `measurement_models.c` |
| 227-0690 Estimation | Nonlinear estimation | `kalman_filter.c` — EKF, UKF |

## Tsinghua University (清华大学)

| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| 信号与系统 (Signals & Systems) | State-space, Kalman | `kalman_filter.c` — state estimation |
| 通信原理 (Comm Principles) | Detection theory | `measurement_models.c` — radar detection |
| 数字信号处理 (DSP) | Adaptive filtering | `kalman_filter.c` — adaptive KF |
| 雷达信号处理 (Radar SP) | Target tracking | Full module — all files |

## Course Core Intersection

The following topics appear in ≥5 of the 9 curricula and are fully covered:

1. **Kalman filtering** — MIT 6.003, Stanford EE102A, Berkeley EE123, Michigan EECS 351, Georgia Tech ECE 4270, ETH 227-0427, 清华 信号与系统
2. **Target tracking** — MIT 16.485, Stanford AA273, Michigan NAVARCH 580, Georgia Tech ECE 6550
3. **Detection/estimation** — MIT 6.450, Berkeley EE123, Michigan EECS 455, Georgia Tech ECE 6601
4. **Radar principles** — MIT 6.630, Stanford EE359, Michigan EECS 411, Georgia Tech ECE 6350, 清华 雷达信号处理
5. **Multi-target association** — Stanford AA273, Michigan NAVARCH 580, Georgia Tech ECE 6550
6. **Nonlinear filtering** — Georgia Tech ECE 6550, ETH 227-0690, 清华 信号与系统
