# Pangu

English version: [README.en.md](README.en.md)

Pangu 是一个基于 Parthenon 和 Kokkos 的 SRMHD/GRMHD 数值模拟代码。
本文档采用“操作优先”的组织方式，目标是让新用户快速完成构建、运行与分析。

## 作者

- Yuehang Li（原始作者）

## 许可证

本项目采用 MIT License，详见 [LICENSE](LICENSE)。

## 目录

- [1. 项目定位](#1-项目定位)
- [2. 仓库结构](#2-仓库结构)
- [3. 环境依赖](#3-环境依赖)
- [4. 快速开始](#4-快速开始)
- [5. 构建指南](#5-构建指南)
- [6. 运行指南](#6-运行指南)
- [7. 分析与可视化](#7-分析与可视化)
- [8. 算例列表](#8-算例列表)
- [9. 输入文件结构](#9-输入文件结构)
- [10. 常见问题](#10-常见问题)
- [11. 开发说明](#11-开发说明)
- [12. 可复现性清单](#12-可复现性清单)

## 1. 项目定位

Pangu 提供以下能力：

- SRMHD 与 GRMHD 两种求解模式
- CPU 与 CUDA 两种构建目标
- 基于问题目录的初始化机制（problem_generator.cpp）
- 构建、运行、分析三段式脚本流程

标准工作流：

1. 使用 scripts/shell/make.sh 构建
2. 使用 scripts/shell/execute.sh 运行
3. 使用 scripts/shell/analyze.sh 分析

## 2. 仓库结构

关键目录：

- pangu/src：求解器核心与程序入口
- pangu/problem：算例输入与初始化
- scripts/shell：构建/运行/分析脚本
- scripts/python：自定义绘图脚本
- parthenon：框架依赖（仓库内）
- data：运行输出目录（运行时生成）
- pic：图像输出目录（分析时生成）

## 3. 环境依赖

### 3.1 构建依赖

- CMake 3.10+（推荐 3.16+）
- 支持 C++17 的编译器
- Python 3

常见可选依赖：

- MPI（多进程运行）
- CUDA toolkit（GPU 构建）
- HDF5（输出链路常用）

### 3.2 Python 分析依赖

```bash
python3 -m pip install --user -r parthenon/requirements.txt
```

主要包括：h5py、numpy、matplotlib。

## 4. 快速开始

以下命令均在仓库根目录执行。

### 4.1 构建（CPU + SRMHD）

```bash
ENABLE_CUDA=OFF PANGU_GR=OFF PROBLEM=brio_wu_shocktube BUILD_DIR=build ./scripts/shell/make.sh
```

预期结果：

- 可执行文件位于 build/pangu/src/pangu.host
- 生成构建快照 .pangu_build.env

### 4.2 运行

```bash
BUILD_DIR=build ENABLE_CUDA=OFF PROBLEM=brio_wu_shocktube ./scripts/shell/execute.sh -n 1
```

预期结果：

- 输出写入 data/brio_wu_shocktube
- 产生 PHDF 文件

### 4.3 分析

```bash
./scripts/shell/analyze.sh -p brio_wu_shocktube -f density
```

预期结果：

- 图片写入 pic/brio_wu_shocktube

## 5. 构建指南

构建脚本：scripts/shell/make.sh

### 5.1 核心参数

| 变量 | 默认值 | 可选值 | 含义 |
| --- | --- | --- | --- |
| ENABLE_OPENMP | ON | ON/OFF | 是否启用 OpenMP |
| ENABLE_CUDA | ON | ON/OFF | 是否构建 CUDA 目标 |
| BUILD_DIR | build | 路径 | 构建目录 |
| BUILD_TYPE | Release | Debug/Release/... | CMake 构建类型 |
| PROBLEM | brio_wu_shocktube | pangu/problem 下目录名 | 算例源选择 |
| PANGU_GR | OFF | ON/OFF | 选择 GRMHD 或 SRMHD |
| BUILD_JOBS | 4 | 整数 | 并行编译线程数 |
| CMAKE_GENERATOR | 空 | 生成器名 | CMake 生成器覆盖 |
| KOKKOS_ARCH | 空 | 架构标识 | Kokkos 架构开关 |
| CMAKE_EXTRA_ARGS | 空 | 参数串 | 额外 CMake 参数 |
| PROBLEM_PROXY_MODE | ON | ON/OFF | 是否启用问题代理模式 |
| PROBLEM_PROXY_NAME | __active_problem__ | 字符串 | 代理目录名称 |

### 5.2 构建矩阵

| 模式 | ENABLE_CUDA | PANGU_GR | 产物 |
| --- | --- | --- | --- |
| CPU SRMHD | OFF | OFF | pangu.host |
| CPU GRMHD | OFF | ON | pangu.host |
| GPU SRMHD | ON | OFF | pangu.cuda |
| GPU GRMHD | ON | ON | pangu.cuda |

## 6. 运行指南

运行脚本：scripts/shell/execute.sh

### 6.1 命令行参数

| 参数 | 含义 |
| --- | --- |
| -i, --input | 输入文件路径 |
| -b, --build-dir | 构建目录 |
| -p, --problem | 算例名称 |
| -n, --np | MPI 进程数 |

### 6.2 环境变量

| 变量 | 默认值 | 含义 |
| --- | --- | --- |
| BUILD_DIR | build | 构建目录 |
| PROBLEM | brio_wu_shocktube | 用于推断默认输入文件 |
| ENABLE_CUDA | ON | 推断默认可执行文件名 |
| MPI_NP | 1 | 未传 -n 时的进程数 |
| INPUT_FILE | 空 | 为空时使用 pangu/problem/<problem>/inputfile |
| DATA_ROOT | <repo>/data | 运行输出根目录 |

说明：

- 当进程数大于 1 时，系统必须提供 mpirun。
- 运行目录为 DATA_ROOT/<problem_name>。

## 7. 分析与可视化

分析脚本：scripts/shell/analyze.sh

### 7.1 模式

| 模式 | 开关 | 输出 |
| --- | --- | --- |
| contour1d | 默认 | 单张图片 |
| movie2d | --movie2d | 帧图目录 |
| xzplot | --xzplot | x-z 变换帧图目录 |

--movie2d 与 --xzplot 互斥。

### 7.2 常用参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| -p, --problem | 推断或必填 | data 下算例目录 |
| -f, --field | density | 目标字段 |
| -w, --workers | 4 | 并行进程数 |
| --data-root | <repo>/data | 数据根目录 |
| --pic-root | <repo>/pic | 图片根目录 |
| --savename | 自动 | 输出文件或目录名 |
| --colorbar | 字段名 | 色标标题（contour 模式） |

### 7.3 xzplot

自定义脚本：scripts/python/xzplot.py

关键参数：

- --output-directory
- --workers
- --kerr-a
- --kerr-h
- --r0
- --x-max

## 8. 算例列表

当前 pangu/problem 下常用算例：

| 算例 | 典型用途 | 建议模式 |
| --- | --- | --- |
| brio_wu_shocktube | 1D shock tube 基准 | SRMHD |
| brio_wu_shocktube_sr | shock tube 的 SR 变体 | SRMHD |
| kelvin_helmholtz | 剪切不稳定性 | SRMHD |
| orszag_tang_vortex | 2D MHD 涡基准 | SRMHD |
| bondi_flow | 含度规项的吸积流 | GRMHD |
| fm_torus | Kerr 参数下的环流盘 | GRMHD |

每个算例目录应包含：

- pangu/problem/<name>/problem_generator.cpp
- pangu/problem/<name>/inputfile

## 9. 输入文件结构

输入文件使用 Parthenon 分段风格，常见段：

- <parthenon/job>
- <parthenon/mesh>
- <parthenon/meshblock>
- <parthenon/time>
- <core>
- <parthenon/output0>

GR 相关算例通常还包括：

- <metric>
- 算例专属参数段（如 <bondi>、<fm_torus>）

建议从已有 inputfile 拷贝并逐组修改参数。

## 10. 常见问题

| 现象 | 原因 | 处理 |
| --- | --- | --- |
| Problem source not found | PROBLEM 与目录名不一致 | 将 PROBLEM 改为 pangu/problem 下真实目录名 |
| execute.sh 提示找不到可执行文件 | 未完成构建或 BUILD_DIR 不一致 | 重新执行 make.sh 并核对 BUILD_DIR |
| MPI requested but mpirun is not available | 缺少 MPI 运行时 | 安装 MPI，或改为 -n 1 |
| No PHDF files found | 预期目录下无输出数据 | 检查运行目录和输入文件输出配置 |
| --movie2d and --xzplot conflict | 同时开启了两个互斥模式 | 只保留一个模式开关 |

## 11. 开发说明

### 11.1 新增算例

1. 创建 pangu/problem/<new_problem>
2. 添加 problem_generator.cpp
3. 添加 inputfile
4. 使用 PROBLEM=<new_problem> 构建
5. 运行并检查 data/<new_problem> 输出

### 11.2 关键集成点

- 程序入口：pangu/src/main.cc
- 驱动与任务图：pangu/src/task_list
- 算例初始化：pangu/problem/<problem>/problem_generator.cpp
- 构建选择逻辑：pangu/src/CMakeLists.txt

## 12. 可复现性清单

建议在结果记录中保存：

- Git commit hash
- 构建变量快照（.pangu_build.env）
- 使用的 inputfile
- 输出目录快照
- 分析命令与模式

这样可以显著降低复现实验与协作排查成本。