// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_RECOVERY_INVERTOR_H
#define PANGU_SRC_RECOVERY_INVERTOR_H
#include <memory>
#include <parthenon/package.hpp>

// Recovers GRMHD primitive variables from conservative variables and metric
// data across all blocks in the partition.
parthenon::TaskStatus Recovery(parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_RECOVERY_INVERTOR_H
