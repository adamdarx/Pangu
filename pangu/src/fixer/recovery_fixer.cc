// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#include "fixer/recovery_fixer.h"

#include <parthenon/package.hpp>
#include <string>
#include <vector>

#include "initialization/variable_mnemonics.h"

parthenon::TaskStatus FixRecovery(parthenon::MeshData<parthenon::Real> *md) {
  using namespace parthenon::package::prelude;
  PARTHENON_INSTRUMENT
  auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
  auto block = IndexRange{0, md->NumBlocks() - 1};

  const auto bound_x1_interior = md->GetBoundsI(IndexDomain::interior);
  const auto bound_x2_interior = md->GetBoundsJ(IndexDomain::interior);
  const auto bound_x3_interior = md->GetBoundsK(IndexDomain::interior);

  PackIndexMap primitiveIndexMap;
  const std::vector<std::string> primitive_tags = {
      "density", "energy", "weighted_velocity", "magnetic_field", "entropy",
      "electron_entropy"
  };
  auto primitive = md->PackVariables(primitive_tags, primitiveIndexMap);
  PackIndexMap flagIndexMap;
  auto flag = md->PackVariables(std::vector<std::string>{"flag"}, flagIndexMap);

  const int offset_x2 = (bound_x2_interior.s != bound_x2_interior.e) ? 1 : 0;
  pmb0->par_for(
      PARTHENON_AUTO_LABEL, block.s, block.e,
      bound_x3_interior.s, bound_x3_interior.e, bound_x2_interior.s, bound_x2_interior.e,
      bound_x1_interior.s, bound_x1_interior.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        if (flag(b, 0, k, j, i) == 0) {
          int pf1 = flag(b, 0, k, j + offset_x2, i - 1);
          int pf2 = flag(b, 0, k, j + offset_x2, i);
          int pf3 = flag(b, 0, k, j + offset_x2, i + 1);
          int pf4 = flag(b, 0, k, j, i + 1);
          int pf5 = flag(b, 0, k, j - offset_x2, i + 1);
          int pf6 = flag(b, 0, k, j - offset_x2, i);
          int pf7 = flag(b, 0, k, j - offset_x2, i - 1);
          int pf8 = flag(b, 0, k, j, i - 1);

          if (pf2 && pf4 && pf6 && pf8) {
            for (int m = 0; m < NPRIM; m++)
              primitive(b, m, k, j, i) =
                  0.25 * (primitive(b, m, k, j + offset_x2, i) +
                          primitive(b, m, k, j - offset_x2, i) +
                          primitive(b, m, k, j, i - 1) +
                          primitive(b, m, k, j, i + 1));
          } else if (pf1 && pf3 && pf5 && pf7) {
            for (int m = 0; m < NPRIM; m++)
              primitive(b, m, k, j, i) =
                  0.25 * (primitive(b, m, k, j + offset_x2, i + 1) +
                          primitive(b, m, k, j + offset_x2, i - 1) +
                          primitive(b, m, k, j - offset_x2, i + 1) +
                          primitive(b, m, k, j - offset_x2, i - 1));
          } else {
            for (int m = 0; m < NPRIM; m++)
              primitive(b, m, k, j, i) =
                  0.125 * (primitive(b, m, k, j + offset_x2, i - 1) +
                           primitive(b, m, k, j + offset_x2, i) +
                           primitive(b, m, k, j + offset_x2, i + 1) +
                           primitive(b, m, k, j, i + 1) +
                           primitive(b, m, k, j - offset_x2, i + 1) +
                           primitive(b, m, k, j - offset_x2, i) +
                           primitive(b, m, k, j - offset_x2, i - 1) +
                           primitive(b, m, k, j, i - 1));
          }
          primitive(b, UX1, k, j, i) = 0.;
          primitive(b, UX2, k, j, i) = 0.;
          primitive(b, UX3, k, j, i) = 0.;
          flag(b, 0, k, j, i) = 1;
        }
      });

  return parthenon::TaskStatus::complete;
}
