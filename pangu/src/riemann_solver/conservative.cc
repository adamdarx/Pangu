// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#include "riemann_solver/conservative.h"

#include <parthenon/package.hpp>
#include <string>
#include <vector>

#include <basic_types.hpp>
#include "initialization/variable_mnemonics.h"
#include "physics/contravariant_flux.h"

parthenon::TaskStatus CalculateConservative(
    parthenon::MeshData<parthenon::Real> *md) {
  using namespace parthenon;
  PARTHENON_INSTRUMENT

  auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
  const auto package_core = pmb0->packages.Get("core");
  const auto kAdiabaticIndex = package_core->Param<Real>("adiabatic_index");

  const auto bound_x1_interior = md->GetBoundsI(IndexDomain::interior);
  const auto bound_x2_interior = md->GetBoundsJ(IndexDomain::interior);
  const auto bound_x3_interior = md->GetBoundsK(IndexDomain::interior);
  auto block = IndexRange{0, md->NumBlocks() - 1};

  PackIndexMap primitiveIndexMap;
  const std::vector<std::string> primitive_tags = {
      "density", "energy", "weighted_velocity", "magnetic_field", "entropy",
      "electron_entropy"
  };
  const auto primitive =
      md->PackVariables(primitive_tags, primitiveIndexMap);

  PackIndexMap conservativeIndexMap;
  const std::vector<std::string> conservative_tags = {"conservative"};
  auto conservative =
      md->PackVariablesAndFluxes(conservative_tags, conservativeIndexMap);

  auto covariant_metric =
      md->PackVariables(std::vector<std::string>{"covariant_metric"});
  auto contravariant_metric =
      md->PackVariables(std::vector<std::string>{"contravariant_metric"});
  auto metric_determinant =
      md->PackVariables(std::vector<std::string>{"metric_determinant"});

  pmb0->par_for(
      PARTHENON_AUTO_LABEL, block.s, block.e,
      bound_x3_interior.s, bound_x3_interior.e,
      bound_x2_interior.s, bound_x2_interior.e,
      bound_x1_interior.s, bound_x1_interior.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        Real gcov[4][4];
        Real gcon[4][4];
        for (int row = 0; row < 4; ++row) {
          for (int col = 0; col < 4; ++col) {
            gcov[row][col] =
                covariant_metric(b, CENTER * 16 + col * 4 + row, k, j, i);
            gcon[row][col] =
                contravariant_metric(b, CENTER * 16 + col * 4 + row, k, j, i);
          }
        }

        Real primitive_c_array[NPRIM];
        for (int index = 0; index < NPRIM; ++index) {
          primitive_c_array[index] = primitive(b, index, k, j, i);
        }

        Real conservative_c_array[NPRIM];
        CalculateContravariantFlux(
            kAdiabaticIndex, primitive_c_array, gcov, gcon,
            metric_determinant(b, CENTER, k, j, i), X0DIR, conservative_c_array);

        for (int index = 0; index < NPRIM; ++index) {
          conservative(b, index, k, j, i) = conservative_c_array[index];
        }
      });

  return TaskStatus::complete;
}
