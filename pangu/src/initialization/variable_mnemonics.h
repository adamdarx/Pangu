// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

// This file in the src/initialization module defines variable_mnemonics.h responsibilities
// for the Pangu runtime. It centers on Index to express core data flow, keep interfaces
// readable, and preserve predictable behavior across task coordination, recovery paths, and
// performance-sensitive execution.

#ifndef PANGU_SRC_INITIALIZATION_VARIABLEMNEMONICS_H
#define PANGU_SRC_INITIALIZATION_VARIABLEMNEMONICS_H

enum Index { RHO, ENY, UX1, UX2, UX3, BX1, BX2, BX3 };

constexpr int NPRIM = 8;

enum Vector3D { X1, X2, X3 };

enum LOC { CENTER, FACEX1, FACEX2, FACEX3 };

#endif  // PANGU_SRC_INITIALIZATION_VARIABLEMNEMONICS_H
