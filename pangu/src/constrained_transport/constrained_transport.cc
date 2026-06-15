// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

#include "constrained_transport/constrained_transport.h"

#include <parthenon/package.hpp>
#include <string>
#include <vector>

#include <basic_types.hpp>
#include "initialization/variable_mnemonics.h"

parthenon::TaskStatus ConstraintedTransport(
    parthenon::MeshData<parthenon::Real> *md) {
  using namespace parthenon;
  PARTHENON_INSTRUMENT

  auto pmb0 = md->GetBlockData(0)->GetBlockPointer();

  const auto bound_x1_interior = md->GetBoundsI(IndexDomain::interior);
  const auto bound_x2_interior = md->GetBoundsJ(IndexDomain::interior);
  const auto bound_x3_interior = md->GetBoundsK(IndexDomain::interior);
  auto block = IndexRange{0, md->NumBlocks() - 1};

  PackIndexMap conservativeIndexMap;
  const std::vector<std::string> conservative_tags = {"conservative"};
  auto conservative =
      md->PackVariablesAndFluxes(conservative_tags, conservativeIndexMap);
  PackIndexMap electricFieldIndexMap;
  const std::vector<std::string> electric_field_tags = {"electric_field"};
  auto electricField =
      md->PackVariables(electric_field_tags, electricFieldIndexMap);

  const int meshgrid_size_x1 =
      pmb0->cellbounds.ncellsi(IndexDomain::entire);
  const int meshgrid_size_x2 =
      pmb0->cellbounds.ncellsj(IndexDomain::entire);
  const int meshgrid_size_x3 =
      pmb0->cellbounds.ncellsk(IndexDomain::entire);

  const int offset_x1 = (meshgrid_size_x1 > 1) ? 1 : 0;
  const int offset_x2 = (meshgrid_size_x2 > 1) ? 1 : 0;
  const int offset_x3 = (meshgrid_size_x3 > 1) ? 1 : 0;

  const int CalculateElectricFieldX1 = offset_x2 && offset_x3;
  const int CalculateElectricFieldX2 = offset_x3 && offset_x1;
  const int CalculateElectricFieldX3 = offset_x1 && offset_x2;

  pmb0->par_for(
      PARTHENON_AUTO_LABEL, block.s, block.e,
      bound_x3_interior.s, bound_x3_interior.e + offset_x3, bound_x2_interior.s,
      bound_x2_interior.e + offset_x2, bound_x1_interior.s, bound_x1_interior.e + offset_x1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        electricField(b, Vector3D::X1, k, j, i) = 0.0;
        electricField(b, Vector3D::X2, k, j, i) = 0.0;
        electricField(b, Vector3D::X3, k, j, i) = 0.0;

        if (CalculateElectricFieldX1)
          electricField(b, Vector3D::X1, k, j, i) =
              0.25 * ((conservative(b).flux(X2DIR, BX3, k, j, i) +
                       conservative(b).flux(X2DIR, BX3, k - 1, j, i)) -
                      (conservative(b).flux(X3DIR, BX2, k, j, i) +
                       conservative(b).flux(X3DIR, BX2, k, j - 1, i)));
        if (CalculateElectricFieldX2)
          electricField(b, Vector3D::X2, k, j, i) =
              0.25 * ((conservative(b).flux(X3DIR, BX1, k, j, i) +
                       conservative(b).flux(X3DIR, BX1, k, j, i - 1)) -
                      (conservative(b).flux(X1DIR, BX3, k, j, i) +
                       conservative(b).flux(X1DIR, BX3, k - 1, j, i)));
        if (CalculateElectricFieldX3)
          electricField(b, Vector3D::X3, k, j, i) =
              0.25 * ((conservative(b).flux(X1DIR, BX2, k, j, i) +
                       conservative(b).flux(X1DIR, BX2, k, j - 1, i)) -
                      (conservative(b).flux(X2DIR, BX1, k, j, i) +
                       conservative(b).flux(X2DIR, BX1, k, j, i - 1)));
      });

  pmb0->par_for(
      PARTHENON_AUTO_LABEL, block.s, block.e,
      bound_x3_interior.s, bound_x3_interior.e, bound_x2_interior.s, bound_x2_interior.e,
      bound_x1_interior.s, bound_x1_interior.e + offset_x1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        if (offset_x1)
          conservative(b).flux(X1DIR, BX1, k, j, i) = 0;
        if (CalculateElectricFieldX2)
          conservative(b).flux(X1DIR, BX3, k, j, i) =
              -0.5 * (electricField(b, Vector3D::X2, k, j, i) +
                      electricField(b, Vector3D::X2, k + 1, j, i));
        if (CalculateElectricFieldX3)
          conservative(b).flux(X1DIR, BX2, k, j, i) =
              0.5 * (electricField(b, Vector3D::X3, k, j, i) +
                     electricField(b, Vector3D::X3, k, j + 1, i));
      });

  pmb0->par_for(
      PARTHENON_AUTO_LABEL, block.s, block.e,
      bound_x3_interior.s, bound_x3_interior.e, bound_x2_interior.s,
      bound_x2_interior.e + offset_x2, bound_x1_interior.s, bound_x1_interior.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        if (CalculateElectricFieldX1)
          conservative(b).flux(X2DIR, BX3, k, j, i) =
              0.5 * (electricField(b, Vector3D::X1, k, j, i) +
                     electricField(b, Vector3D::X1, k + 1, j, i));
        if (offset_x2)
          conservative(b).flux(X2DIR, BX2, k, j, i) = 0.;
        if (CalculateElectricFieldX3)
          conservative(b).flux(X2DIR, BX1, k, j, i) =
              -0.5 * (electricField(b, Vector3D::X3, k, j, i) +
                      electricField(b, Vector3D::X3, k, j, i + 1));
      });

  pmb0->par_for(
      PARTHENON_AUTO_LABEL, block.s, block.e,
      bound_x3_interior.s, bound_x3_interior.e + offset_x3, bound_x2_interior.s,
      bound_x2_interior.e, bound_x1_interior.s, bound_x1_interior.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        if (CalculateElectricFieldX1)
          conservative(b).flux(X3DIR, BX2, k, j, i) =
              -0.5 * (electricField(b, Vector3D::X1, k, j, i) +
                      electricField(b, Vector3D::X1, k, j + 1, i));
        if (CalculateElectricFieldX2)
          conservative(b).flux(X3DIR, BX1, k, j, i) =
              0.5 * (electricField(b, Vector3D::X2, k, j, i) +
                     electricField(b, Vector3D::X2, k, j, i + 1));
        if (offset_x3)
          conservative(b).flux(X3DIR, BX3, k, j, i) = 0;
      });

  return parthenon::TaskStatus::complete;
}
