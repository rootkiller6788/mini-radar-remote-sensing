# Mini Radar Remote Sensing（雷达与遥感）

一套**从零开始、零外部依赖的 C 语言实现**，涵盖雷达系统、遥感与信号处理理论。每个模块对标 MIT、Stanford 等顶尖高校课程，将教科书中的公式转化为可运行的 C 代码，覆盖从波形生成、天线波束赋形到检测、跟踪与成像的完整雷达信号链路。

## 子模块

| 子模块 | 主题 | 对标课程 |
|--------|------|----------|
| [mini-radar-basics](mini-radar-basics/) | 雷达信号模型、检测理论（CFAR/Marcum/Swerling）、多普勒/MTI 处理、波形设计（LFM/Barker/Costas） | MIT 6.630, Stanford EE359, Skolnik |
| [mini-pulse-doppler](mini-pulse-doppler/) | 雷达距离方程、波形生成、多普勒处理、CFAR 检测（CA/GO/SO/OS）、模糊函数分析 | Stanford EE359, Georgia Tech ECE 6350, Richards |
| [mini-phased-array](mini-phased-array/) | 阵列几何结构、相位/时延波束赋形、自适应波束赋形（SMI/MVDR/LCMV）、方向图综合（Dolph-Chebyshev/Taylor） | Stanford EE252, Van Trees, Balanis |
| [mini-sar-imaging](mini-sar-imaging/) | SAR 几何、聚焦算法（RDA/CSA/ωK）、干涉测量（InSAR/相干性）、压缩感知 SAR、MIMO-SAR | MIT 16.851, Stanford EE355, Georgia Tech ECE 6350 |
| [mini-infrared-thermal](mini-infrared-thermal/) | 红外核心物理与定律、大气传输（MODTRAN）、探测器模型/NUC 标定、辐射测量、目标检测、图像处理 | MIT 22.071, Stanford EE290, ETH 227-0436 |
| [mini-hyperspectral](mini-hyperspectral/) | 光谱学原理、辐射定标与大气校正、降维（PCA/MNF）、分类（SAM/SVM）、光谱解混（LMM） | MIT 12.710, Stanford EE369, Purdue ECE 637 |
| [mini-lidar-principle](mini-lidar-principle/) | LiDAR 检测理论（SNR/Pd/Pfa）、点云几何、全波形处理、扫描模式、ICP 配准、DEM/林业应用 | MIT 16.851, Stanford EE368, 慕尼黑工大高频工程 |
| [mini-target-tracking](mini-target-tracking/) | Kalman 滤波器变体（KF/EKF/UKF/SR-KF）、运动模型（CV/CA/CT/IMM）、数据关联（NN/GNN/PDA/JPDA/MHT）、航迹融合、毫米波跟踪 | Stanford EE363, Georgia Tech ECE 6601, Bar-Shalom |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅依赖 `libc` 和 `libm`
- **独立自包含** — 每个目录拥有自己的 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，附有对标课程的对照说明，引用标准教科书（Skolnik、Balanis、Van Trees、Bar-Shalom）
- **实用演示** — 脉冲多普勒处理链、自适应波束赋形模拟器、SAR 图像生成、多目标跟踪器等

## 构建

每个模块独立运行。进入模块目录并执行：

```bash
cd mini-radar-basics
make all    # 构建所有目标
make test   # 运行测试
```

依赖 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-radar-remote-sensing/
├── mini-radar-basics/          # 雷达基础、检测、多普勒、波形设计
├── mini-pulse-doppler/         # 脉冲多普勒处理、距离方程、模糊函数
├── mini-phased-array/          # 相控阵波束赋形、自适应阵列、方向图综合
├── mini-sar-imaging/           # 合成孔径雷达成像与干涉测量
├── mini-infrared-thermal/      # 红外物理、探测器、辐射测量、图像处理
├── mini-hyperspectral/         # 高光谱光谱学、解混、分类
├── mini-lidar-principle/       # 激光雷达检测、几何、波形、配准
└── mini-target-tracking/       # 卡尔曼滤波、运动模型、数据关联、航迹融合
```

## 许可证

MIT
