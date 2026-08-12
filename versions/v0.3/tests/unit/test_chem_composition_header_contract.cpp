// SPDX-License-Identifier: Apache-2.0

#include "hundun/chem_composition.hpp"
#include "hundun/chem_reports.hpp"

#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_same_v<
              decltype(hundun::chemistry::SpeciesIdentity::element_counts),
              std::vector<std::int32_t>>);
static_assert(std::is_same_v<
              decltype(hundun::chemistry::CompositionIdentity::fingerprint),
              std::uint64_t>);
static_assert(std::is_same_v<
              decltype(hundun::chemistry::ThermochemicalPoint::mass_fractions),
              std::vector<double>>);
static_assert(std::is_same_v<
              decltype(std::declval<const hundun::chemistry::ChemistryIntervalReport &>()
                           .succeeded()),
              bool>);
static_assert(noexcept(
    std::declval<const hundun::chemistry::ChemistryIntervalReport &>()
        .succeeded()));
static_assert(sizeof(hundun::chemistry::ChemistryStatus) ==
              sizeof(std::uint32_t));

int main() { return 0; }
