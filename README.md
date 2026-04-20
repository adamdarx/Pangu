# Pangu

Chinese version: [README.md](README.md)

Pangu is an SRMHD/GRMHD simulation code built on top of Parthenon and Kokkos.
This guide is organized in an operations-first way so new users can build, run, and analyze quickly.

## Authors

- Yuehang Li (original author)

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

## Table of Contents

- [1. Project Scope](#1-project-scope)
- [2. Repository Layout](#2-repository-layout)
- [3. Requirements](#3-requirements)
- [4. Quick Start](#4-quick-start)
- [5. Build Guide](#5-build-guide)
- [6. Run Guide](#6-run-guide)
- [7. Analysis and Visualization](#7-analysis-and-visualization)
- [8. Problem Set](#8-problem-set)
- [9. Input File Structure](#9-input-file-structure)
- [10. Troubleshooting](#10-troubleshooting)
- [11. Developer Notes](#11-developer-notes)
- [12. Reproducibility Checklist](#12-reproducibility-checklist)

## 1. Project Scope

Pangu provides:

- SRMHD and GRMHD solver modes
- CPU and CUDA build targets
- Problem-based initialization via problem_generator.cpp
- Unified shell scripts for build, run, and analysis

Standard workflow:

1. Build with scripts/shell/make.sh
2. Run with scripts/shell/execute.sh
3. Analyze with scripts/shell/analyze.sh

## 2. Repository Layout

Key directories:

- pangu/src: solver core and runtime entry
- pangu/problem: problem-specific input files and generators
- scripts/shell: operational scripts (build, run, analyze)
- scripts/python: custom plotting tools
- parthenon: framework dependency (in-repo)
- data: runtime output directory (generated)
- pic: analysis output directory (generated)

## 3. Requirements

### 3.1 Build Requirements

- CMake 3.10+ (3.16+ recommended)
- C++17 compiler
- Python 3

Optional but common:

- MPI (multi-process runs)
- CUDA toolkit (GPU build)
- HDF5 (output workflow)

### 3.2 Python Requirements (Analysis)

```bash
python3 -m pip install --user -r parthenon/requirements.txt
```

Main packages: h5py, numpy, matplotlib.

## 4. Quick Start

Run all commands from the repository root.

### 4.1 Build (CPU + SRMHD)

```bash
ENABLE_CUDA=OFF PANGU_GR=OFF PROBLEM=brio_wu_shocktube BUILD_DIR=build ./scripts/shell/make.sh
```

Expected results:

- Executable at build/pangu/src/pangu.host
- Build snapshot at .pangu_build.env

### 4.2 Run

```bash
BUILD_DIR=build ENABLE_CUDA=OFF PROBLEM=brio_wu_shocktube ./scripts/shell/execute.sh -n 1
```

Expected results:

- Outputs written to data/brio_wu_shocktube
- PHDF files generated

### 4.3 Analyze

```bash
./scripts/shell/analyze.sh -p brio_wu_shocktube -f density
```

Expected results:

- Figures generated under pic/brio_wu_shocktube

## 5. Build Guide

Build script: scripts/shell/make.sh

### 5.1 Core Variables

| Variable | Default | Allowed Values | Meaning |
| --- | --- | --- | --- |
| ENABLE_OPENMP | ON | ON/OFF | Enable OpenMP backend |
| ENABLE_CUDA | ON | ON/OFF | Enable CUDA target |
| BUILD_DIR | build | path | Build directory |
| BUILD_TYPE | Release | Debug/Release/... | CMake build type |
| PROBLEM | brio_wu_shocktube | folder in pangu/problem | Problem source selection |
| PANGU_GR | OFF | ON/OFF | Select GRMHD or SRMHD |
| BUILD_JOBS | 4 | integer | Parallel build jobs |
| CMAKE_GENERATOR | empty | generator name | CMake generator override |
| KOKKOS_ARCH | empty | architecture token | Kokkos architecture selector |
| CMAKE_EXTRA_ARGS | empty | arg string | Extra CMake args |
| PROBLEM_PROXY_MODE | ON | ON/OFF | Enable problem proxy mode |
| PROBLEM_PROXY_NAME | __active_problem__ | string | Proxy folder name |

### 5.2 Build Matrix

| Mode | ENABLE_CUDA | PANGU_GR | Output Binary |
| --- | --- | --- | --- |
| CPU SRMHD | OFF | OFF | pangu.host |
| CPU GRMHD | OFF | ON | pangu.host |
| GPU SRMHD | ON | OFF | pangu.cuda |
| GPU GRMHD | ON | ON | pangu.cuda |

## 6. Run Guide

Run script: scripts/shell/execute.sh

### 6.1 CLI Options

| Option | Meaning |
| --- | --- |
| -i, --input | Input file path |
| -b, --build-dir | Build directory |
| -p, --problem | Problem name |
| -n, --np | MPI process count |

### 6.2 Environment Variables

| Variable | Default | Meaning |
| --- | --- | --- |
| BUILD_DIR | build | Build directory |
| PROBLEM | brio_wu_shocktube | Default input inference |
| ENABLE_CUDA | ON | Default executable inference |
| MPI_NP | 1 | Process count without -n |
| INPUT_FILE | empty | Defaults to pangu/problem/<problem>/inputfile |
| DATA_ROOT | <repo>/data | Output root directory |

Notes:

- If process count is greater than 1, mpirun must be available.
- The script runs in DATA_ROOT/<problem_name>.

## 7. Analysis and Visualization

Analysis script: scripts/shell/analyze.sh

### 7.1 Modes

| Mode | Flag | Output |
| --- | --- | --- |
| contour1d | default | Single image |
| movie2d | --movie2d | Frame directory |
| xzplot | --xzplot | x-z transformed frames |

--movie2d and --xzplot are mutually exclusive.

### 7.2 Common Options

| Option | Default | Meaning |
| --- | --- | --- |
| -p, --problem | inferred or required | Problem under data |
| -f, --field | density | Field to visualize |
| -w, --workers | 4 | Worker count |
| --data-root | <repo>/data | Data root |
| --pic-root | <repo>/pic | Figure root |
| --savename | auto | Output file or directory name |
| --colorbar | field name | Colorbar label (contour mode) |

### 7.3 xzplot

Custom script: scripts/python/xzplot.py

Key arguments:

- --output-directory
- --workers
- --kerr-a
- --kerr-h
- --r0
- --x-max

## 8. Problem Set

Current common directories under pangu/problem:

| Problem | Typical Use | Suggested Mode |
| --- | --- | --- |
| brio_wu_shocktube | 1D shock tube baseline | SRMHD |
| brio_wu_shocktube_sr | SR shock tube variant | SRMHD |
| kelvin_helmholtz | Shear-flow instability | SRMHD |
| orszag_tang_vortex | 2D MHD vortex benchmark | SRMHD |
| bondi_flow | Accretion flow with metric terms | GRMHD |
| fm_torus | Torus with Kerr parameters | GRMHD |

Each problem directory should include:

- pangu/problem/<name>/problem_generator.cpp
- pangu/problem/<name>/inputfile

## 9. Input File Structure

Input files follow Parthenon-style blocks. Common sections:

- <parthenon/job>
- <parthenon/mesh>
- <parthenon/meshblock>
- <parthenon/time>
- <core>
- <parthenon/output0>

GR-focused setups usually also include:

- <metric>
- problem-specific sections (for example <bondi> or <fm_torus>)

Recommendation: start from an existing inputfile and modify one parameter group at a time.

## 10. Troubleshooting

| Symptom | Likely Cause | Fix |
| --- | --- | --- |
| Problem source not found | PROBLEM does not match an existing folder | Set PROBLEM to a valid folder in pangu/problem |
| Executable not found in execute.sh | Build incomplete or wrong BUILD_DIR | Re-run make.sh and verify BUILD_DIR |
| MPI requested but mpirun is not available | Missing MPI runtime | Install MPI or run with -n 1 |
| No PHDF files found | No simulation output in expected directory | Verify run directory and output settings |
| --movie2d and --xzplot conflict | Both mutually-exclusive flags are enabled | Use only one mode flag |

## 11. Developer Notes

### 11.1 Add a New Problem

1. Create pangu/problem/<new_problem>
2. Add problem_generator.cpp
3. Add inputfile
4. Build with PROBLEM=<new_problem>
5. Run and verify data/<new_problem>

### 11.2 Key Integration Points

- Runtime entry: pangu/src/main.cc
- Driver/task graph: pangu/src/task_list
- Problem initialization: pangu/problem/<problem>/problem_generator.cpp
- Build selection: pangu/src/CMakeLists.txt

## 12. Reproducibility Checklist

Before reporting performance or physics results, record:

- Git commit hash
- Build variables (.pangu_build.env)
- Input file used
- Output directory snapshot
- Analysis command and mode

This greatly improves repeatability and collaborative debugging.
