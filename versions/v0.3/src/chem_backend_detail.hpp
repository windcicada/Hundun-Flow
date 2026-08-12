// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/chem_reports.hpp"

#include <stdexcept>

namespace hundun::chemistry {

class ChemistryBackend {
public:
  virtual ~ChemistryBackend() = default;
  virtual const CompositionIdentity &composition() const noexcept = 0;
  virtual ChemistryIntervalReport
  integrate(const ChemistryIntervalRequest &) = 0;
};

inline bool same_composition_identity(const CompositionIdentity &left,
                                      const CompositionIdentity &right) {
  if (left.fingerprint != right.fingerprint ||
      left.element_names != right.element_names ||
      left.species.size() != right.species.size()) {
    return false;
  }
  for (std::size_t index = 0; index < left.species.size(); ++index) {
    const auto &left_species = left.species[index];
    const auto &right_species = right.species[index];
    if (left_species.name != right_species.name ||
        left_species.molecular_weight_kg_per_kmol !=
            right_species.molecular_weight_kg_per_kmol ||
        left_species.element_counts != right_species.element_counts) {
      return false;
    }
  }
  return true;
}

inline void validate_backend_composition(const CompositionIdentity &expected,
                                         const ChemistryBackend &backend) {
  validate_composition_identity(expected);
  validate_composition_identity(backend.composition());
  if (!same_composition_identity(expected, backend.composition())) {
    throw std::invalid_argument("chemistry backend composition mismatch");
  }
}

} // namespace hundun::chemistry
