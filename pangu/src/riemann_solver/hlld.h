// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_RIEMANNSOLVER_HLLD_H
#define PANGU_SRC_RIEMANNSOLVER_HLLD_H

#include <memory>
#include <parthenon/package.hpp>

// Computes GRMHD numerical fluxes with the HLLD solver
// across all blocks in the partition.
parthenon::TaskStatus CalculateHLLD(parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_RIEMANNSOLVER_HLLD_H
