// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_RIEMANNSOLVER_HLL_H
#define PANGU_SRC_RIEMANNSOLVER_HLL_H

#include <memory>
#include <parthenon/package.hpp>

// Computes GRMHD numerical fluxes with the HLL solver
// across all blocks in the partition.
parthenon::TaskStatus CalculateHLL(parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_RIEMANNSOLVER_HLL_H
