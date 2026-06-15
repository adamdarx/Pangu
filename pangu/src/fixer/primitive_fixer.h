// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_FIXER_PRIMITIVEFIXER_H
#define PANGU_SRC_FIXER_PRIMITIVEFIXER_H

#include <memory>
#include <parthenon/package.hpp>

// Applies GRMHD primitive floors and velocity ceilings using metric data
// across all blocks in the partition.
parthenon::TaskStatus FixPrimitive(parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_FIXER_PRIMITIVEFIXER_H
