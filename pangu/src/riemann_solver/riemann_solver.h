// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_RIEMANNSOLVER_RIEMANNSOLVER_H
#define PANGU_SRC_RIEMANNSOLVER_RIEMANNSOLVER_H

#include <memory>
#include <parthenon/package.hpp>

// Dispatches to the configured GRMHD Riemann solver
// across all blocks in the partition.
parthenon::TaskStatus CalculateFluxes(parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_RIEMANNSOLVER_RIEMANNSOLVER_H
