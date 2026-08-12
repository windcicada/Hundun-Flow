// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/chem_composition.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::chemistry {
class ThermodynamicsService;
}

namespace hundun::flow {

struct ReactingFieldIds final {
  std::vector<runtime::FieldId> species_density;
  runtime::FieldId total_thermochemical_enthalpy{};
};

ReactingFieldIds
declare_reacting_fields(runtime::FieldRegistry &,
                        const chemistry::CompositionIdentity &);

enum class ReactingLayer : std::uint8_t { history, committed, trial };

struct ReactingLayerValues final {
  std::vector<std::vector<double>> rho_y_kg_per_m3;
  std::vector<double> rho_h_tc_j_per_m3;
};

class ReactingFlowState final {
public:
  static ReactingFlowState
  create(const runtime::FieldRegistry &, runtime::Int3,
         ReactingFieldIds, chemistry::CompositionIdentity);

  ~ReactingFlowState() noexcept;
  ReactingFlowState(ReactingFlowState &&) noexcept;
  ReactingFlowState &operator=(ReactingFlowState &&) noexcept;
  ReactingFlowState(const ReactingFlowState &) = delete;
  ReactingFlowState &operator=(const ReactingFlowState &) = delete;

  void initialize_uniform(const chemistry::ThermodynamicsService &,
                          double p0_pa, double temperature_k,
                          const std::vector<double> &mass_fractions);
  ReactingLayerValues snapshot(ReactingLayer) const;
  void begin_attempt();
  void write_trial_cell(std::size_t,
                        const std::vector<double> &rho_y_kg_per_m3,
                        double rho_h_tc_j_per_m3);
  void rollback_attempt();
  void commit_attempt();
  chemistry::ThermodynamicProperties
  thermodynamics(ReactingLayer, std::size_t,
                 const chemistry::ThermodynamicsService &) const;
  std::uint64_t writer_epoch() const noexcept;
  std::uint64_t cache_epoch() const noexcept;

private:
  struct Impl;
  explicit ReactingFlowState(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
};

} // namespace hundun::flow
