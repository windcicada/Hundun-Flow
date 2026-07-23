// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/momentum_predictor.hpp"
#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/field_storage.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::flow {

namespace detail {
struct FlowStateSolverAccess;
struct AdaptiveTimeControlAccess;
}
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class MaterialDensityTransportTestAccess;
class MaterialDensityPisoTestAccess;
}
#endif

enum class FlowLayer : std::uint8_t { history, committed, trial };

struct FlowFieldIds final {
  runtime::FieldId density{};
  runtime::FieldId velocity{};
  runtime::FieldId mechanical_pressure{};
  runtime::FieldId face_velocity{};
  runtime::FieldId face_mass_flux{};
  std::vector<runtime::FieldId> transported_cell_fields;
};

struct AcceptedStepMetadata final {
  std::uint64_t step{};
  double time_s{};
  double dt_s{};
  double previous_dt_s{};
  MomentumTimeOrder order{MomentumTimeOrder::backward_euler};
};

struct FlowLayerValues final {
  std::vector<double> density;
  std::vector<double> velocity;
  std::vector<double> mechanical_pressure;
  std::vector<double> face_velocity;
  std::vector<double> face_mass_flux;
  std::vector<std::vector<double>> transported_cell_fields;
};

class PisoCoupler;
class FixedStepConstantDensityFlow;
class MaterialDensityTransport;
class MaterialDensityDiagnosticSource;
class FixedStepMaterialDensityFlow;
class MaterialDensityFlowDiagnosticSource;
class IdealGasClosure;
class FixedStepIdealGasFlow;
class IdealGasClosureDiagnosticSource;

class FlowState final {
public:
  static FlowState create(const runtime::FieldRegistry &registry,
                          runtime::FieldLayoutSet layout, FlowFieldIds fields,
                          AcceptedStepMetadata metadata);

  ~FlowState() noexcept;
  FlowState(FlowState &&) noexcept;
  FlowState &operator=(FlowState &&) noexcept;
  FlowState(const FlowState &) = delete;
  FlowState &operator=(const FlowState &) = delete;

  const runtime::FieldStorage &layer(FlowLayer selected) const;
  runtime::FieldStorage &trial_layer();
  const FlowFieldIds &fields() const noexcept;
  AcceptedStepMetadata metadata() const noexcept;

  void seed_accepted_layers(const FlowLayerValues &history,
                            const FlowLayerValues &committed);
  FlowLayerValues snapshot(FlowLayer selected) const;

  void begin_attempt();
  void rollback_attempt();
  void commit_attempt(AcceptedStepMetadata accepted);

private:
  struct Impl;
  explicit FlowState(std::unique_ptr<Impl>) noexcept;
  const runtime::FieldRegistry &solver_registry() const noexcept;
  const runtime::FieldAccessPlan &solver_access_plan() const noexcept;
  runtime::FieldStorage &solver_layer(FlowLayer selected);
  bool attempt_active() const noexcept;
  std::uint64_t attempt_identity() const noexcept;
  std::uint64_t diagnostic_mutation_identity() const noexcept;
  void prepare_commit_attempt(AcceptedStepMetadata accepted);
  void publish_commit_attempt() noexcept;
  std::unique_ptr<Impl> impl_;

  friend class PisoCoupler;
  friend class FixedStepConstantDensityFlow;
  friend class MaterialDensityTransport;
  friend class MaterialDensityDiagnosticSource;
  friend class FixedStepMaterialDensityFlow;
  friend class MaterialDensityFlowDiagnosticSource;
  friend class IdealGasClosure;
  friend class FixedStepIdealGasFlow;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::MaterialDensityTransportTestAccess;
  friend class test::MaterialDensityPisoTestAccess;
#endif
  friend struct detail::FlowStateSolverAccess;
  friend struct detail::AdaptiveTimeControlAccess;
};

} // namespace hundun::flow
