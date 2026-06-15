// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_RIEMANNSOLVER_ELECTRONHEATING_H
#define PANGU_SRC_RIEMANNSOLVER_ELECTRONHEATING_H

#include <memory>
#include <parthenon/package.hpp>

// Applies electron heating using advected and recovered entropy values
// across all blocks in the partition.
parthenon::TaskStatus ApplyElectronHeating(
    parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_RIEMANNSOLVER_ELECTRONHEATING_H