// Copyright (c) 2026 Yuehang Li.
// Licensed under the MIT License. See LICENSE in the project root for license information.

// This file in the src/physics module defines two_temperature.h responsibilities for
// the Pangu runtime. It centers on helper functions to express core data flow, keep
// interfaces readable, and preserve predictable behavior across task coordination,
// recovery paths, and performance-sensitive execution.

#ifndef PANGU_SRC_PHYSICS_TWOTEMPERATURE_H
#define PANGU_SRC_PHYSICS_TWOTEMPERATURE_H

#include <parthenon/package.hpp>

#include "basic_types.hpp"

namespace two_temperature {

KOKKOS_INLINE_FUNCTION
parthenon::Real ComputeTotalEntropy(const parthenon::Real gamma,
                                    const parthenon::Real rho,
                                    const parthenon::Real energy) {
  constexpr parthenon::Real kSmall = 1.0e-20;
  return (gamma - 1.0) * energy * Kokkos::pow(Kokkos::max(rho, kSmall), -gamma);
}

KOKKOS_INLINE_FUNCTION
parthenon::Real RecoverAdvectedScalar(const parthenon::Real cons_rho,
                                      const parthenon::Real cons_scalar) {
  constexpr parthenon::Real kSmall = 1.0e-20;
  return cons_scalar / Kokkos::max(cons_rho, kSmall);
}

KOKKOS_INLINE_FUNCTION
parthenon::Real ComputeConstantModelElectronEntropy(
    const parthenon::Real gamma,
    const parthenon::Real gamma_p,
    const parthenon::Real gamma_e,
    const parthenon::Real rho,
    const parthenon::Real bsq,
    const parthenon::Real advected_total_entropy,
    const parthenon::Real recovered_total_entropy,
    const parthenon::Real advected_electron_entropy,
    const parthenon::Real fel_constant,
    const bool limit_kel,
    const parthenon::Real ratio_min,
    const parthenon::Real ratio_max,
    const bool suppress_highb_heat) {
  constexpr parthenon::Real kSmall = 1.0e-20;
  const parthenon::Real rho_safe = Kokkos::max(rho, kSmall);
  const parthenon::Real fel = Kokkos::max(0.0, Kokkos::min(1.0, fel_constant));

  parthenon::Real dissipation =
      (gamma_e - 1.0) / (gamma - 1.0) * Kokkos::pow(rho_safe, gamma - gamma_e) *
      (recovered_total_entropy - advected_total_entropy);
  if (suppress_highb_heat && (bsq / rho_safe > 1.0)) {
    dissipation = 0.0;
  }

  parthenon::Real kel = advected_electron_entropy + fel * dissipation;
  if (!limit_kel) {
    return kel;
  }

  const parthenon::Real kel_max =
      advected_total_entropy * Kokkos::pow(rho_safe, gamma - gamma_e) /
      (ratio_min * (gamma - 1.0) / (gamma_p - 1.0) +
       (gamma - 1.0) / (gamma_e - 1.0));
  const parthenon::Real kel_min =
      advected_total_entropy * Kokkos::pow(rho_safe, gamma - gamma_e) /
      (ratio_max * (gamma - 1.0) / (gamma_p - 1.0) +
       (gamma - 1.0) / (gamma_e - 1.0));
  return Kokkos::max(kel_min, Kokkos::min(kel, kel_max));
}

KOKKOS_INLINE_FUNCTION
parthenon::Real ClampElectronEntropyByRatio(const parthenon::Real total_entropy,
                                            const parthenon::Real electron_entropy,
                                            const parthenon::Real ratio_min,
                                            const parthenon::Real ratio_max) {
  constexpr parthenon::Real kSmall = 1.0e-20;
  const parthenon::Real rmin = Kokkos::max(ratio_min, kSmall);
  const parthenon::Real rmax = Kokkos::max(ratio_max, rmin);

  const parthenon::Real kel_low = total_entropy / (1.0 + rmax);
  const parthenon::Real kel_high = total_entropy / (1.0 + rmin);
  const parthenon::Real kel_bounded = Kokkos::max(0.0, Kokkos::min(electron_entropy, total_entropy));
  return Kokkos::max(kel_low, Kokkos::min(kel_bounded, kel_high));
}

}  // namespace two_temperature

#endif  // PANGU_SRC_PHYSICS_TWOTEMPERATURE_H