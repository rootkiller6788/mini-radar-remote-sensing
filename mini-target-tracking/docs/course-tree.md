# Course Tree — mini-target-tracking

## Prerequisite Dependency Tree

This module depends on the following prerequisite knowledge. Each node lists the prerequisite mini-modules that should be studied first.

```
mini-target-tracking (THIS MODULE)
│
├── mini-radar-basics (雷达基础)
│   └── Radar equation, PRF, range/Doppler ambiguity
│       └── Used in: measurement_models.c (radar sensor model, CRLB)
│
├── mini-signal-basis (信号基础)
│   └── Probability, Gaussian distributions, hypothesis testing
│       └── Used in: kalman_filter.c (Gaussian state-space)
│
├── mini-fourier-analysis (傅里叶分析)
│   └── Fourier transform, frequency domain analysis
│       └── Used in: measurement_models.c (Doppler processing)
│
├── mini-filter-theory (滤波器理论)
│   └── Wiener filter, optimal filtering
│       └── Used in: kalman_filter.c (Kalman optimality)
│
├── mini-system-analysis (系统分析)
│   └── State-space models, linear systems
│       └── Used in: motion_models.c (state transition)
│
├── mini-laplace-z-transform (拉普拉斯/Z变换)
│   └── Discrete-time systems, z-transform
│       └── Used in: motion_models.c (discretization)
│
├── mini-adaptive-filter (自适应滤波)
│   └── LMS, RLS algorithms
│       └── Used in: kalman_filter.c (adaptive KF)
│
├── mini-nonlinear-system (非线性系统)
│   └── Linearization, Jacobian
│       └── Used in: kalman_filter.c (EKF, UKF)
│
├── mini-system-identification (系统辨识)
│   └── Parameter estimation, model selection
│       └── Used in: kalman_filter.c (model likelihood, IMM)
│
└── mini-coherent-detection (相干检测)
    └── Matched filter, detection theory
        └── Used in: measurement_models.c (radar detection)
```

## Forward Dependencies (Modules that depend on THIS)

```
mini-target-tracking (THIS)
│
├── mini-pulse-doppler (脉冲多普勒)
│   └── Uses: measurement_models.c for Doppler processing
│
├── mini-phased-array (相控阵)
│   └── Uses: kalman_filter.c for beam-space tracking
│
├── mini-sar-imaging (SAR成像)
│   └── Uses: motion_models.c for platform motion
│
├── mini-lidar-principle (激光雷达)
│   └── Uses: data_association.c for point cloud association
│
├── mini-infrared-thermal (红外热成像)
│   └── Uses: track_fusion.c for IR+Radar fusion
│
├── mini-hyperspectral (高光谱)
│   └── Uses: measurement_models.c for spectral tracking
│
├── mini-5g-nr-phy (5G NR)
│   └── Uses: data_association.c for beam management
│
└── mini-navigation-positioning (导航定位)
    └── Uses: kalman_filter.c for GPS/INS integration
```

## Internal Dependency Graph

```
track_core.h/c ─────────────── (foundation: types, matrix ops, scoring, gating)
│
├── kalman_filter.h/c ──────── (KF variants, IMM, adaptive)
│   ├── motion_models.h/c ──── (CV, CA, CT, Singer, discretization)
│   ├── measurement_models.h/c (coordinate transforms, CRLB, CMKF)
│   └── track_metrics.h/c ──── (RMSE, NEES, ANEES, OSPA, CLEAR MOT)
│
├── data_association.h/c ───── (GNN, NN, PDA, JPDA, MHT, Auction)
│   └── track_core.h/c (gated cost matrix)
│
├── track_fusion.h/c ───────── (CI, WCF, Info Fusion, tracklet, decorrelation)
│   └── kalman_filter.h/c (matrix operations)
│
└── mmwave_tracker.h/c ─────── (DBSCAN, OGM, classification, ego-motion)
    ├── track_core.h/c
    ├── motion_models.h/c
    └── kalman_filter.h/c
```

## Study Path

### Beginner Path (L1-L3)
1. Start with `track_core.h/c` — understand data structures, matrix operations
2. Study `motion_models.h/c` — learn state transition models
3. Study `measurement_models.h/c` — coordinate systems, radar measurements

### Intermediate Path (L4-L5)
4. Study `kalman_filter.h/c` — KF theory, prediction, update
5. Study `data_association.h/c` — nearest neighbor, Hungarian, PDA
6. Study `track_metrics.h/c` — evaluation metrics

### Advanced Path (L6-L8)
7. Study `kalman_filter.c` — EKF, UKF, SR-KF, adaptive KF, IMM
8. Study `track_fusion.c` — CI, WCF, information fusion
9. Study `mmwave_tracker.c` — DBSCAN, OGM, extended objects
10. Study examples/ — end-to-end scenarios

### Research Path (L9)
11. Explore documented frontiers: JCAS, AI trackers, sub-THz, 4D radar
