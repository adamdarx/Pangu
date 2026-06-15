// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_RIEMANNSOLVER_SOURCETERM_H
#define PANGU_SRC_RIEMANNSOLVER_SOURCETERM_H

#include <memory>
#include <parthenon/package.hpp>

// Adds geometric source terms to the GRMHD conservative variables
// across all blocks in the partition.
parthenon::TaskStatus AddGeometricSource(parthenon::MeshData<parthenon::Real> *md,
                                         parthenon::Real dt);
#endif  // PANGU_SRC_RIEMANNSOLVER_SOURCETERM_H
