// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#ifndef PANGU_SRC_CONSTRAINEDTRANSPORT_CONSTRAINEDTRANSPORT_H
#define PANGU_SRC_CONSTRAINEDTRANSPORT_CONSTRAINEDTRANSPORT_H

#include <parthenon/package.hpp>

// Updates electric fields and magnetic fluxes with constrained transport
// across all blocks in the partition.
parthenon::TaskStatus ConstraintedTransport(parthenon::MeshData<parthenon::Real> *md);

#endif  // PANGU_SRC_CONSTRAINEDTRANSPORT_CONSTRAINEDTRANSPORT_H
