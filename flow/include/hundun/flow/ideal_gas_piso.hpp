// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/ideal_gas_closure.hpp"
#include "hundun/flow/material_density_piso.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string_view>

namespace hundun::flow {

namespace detail {
struct DensityClosureAdapter;
struct DensityClosureDiagnosticAccess;
} // namespace detail

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class IdealGasClosureTestAccess;
}
#endif

class IdealGasStepAttemptReport final {
public:
  IdealGasStepAttemptReport(const IdealGasStepAttemptReport &) = default;
  IdealGasStepAttemptReport(IdealGasStepAttemptReport &&) noexcept = default;
  IdealGasStepAttemptReport &
  operator=(const IdealGasStepAttemptReport &) = default;
  IdealGasStepAttemptReport &
  operator=(IdealGasStepAttemptReport &&) noexcept = default;

  const MaterialDensityStepAttemptReport &flow() const noexcept;
  bool closure_report_available() const noexcept;
  const IdealGasClosureReport &closure_report() const;
  std::uint64_t attempt_identity() const noexcept;

private:
  IdealGasStepAttemptReport();
  MaterialDensityStepAttemptReport flow_;
  std::optional<IdealGasClosureReport> closure_report_;
  std::uint64_t attempt_identity_{};
  std::uint64_t seal_{};
  std::uint64_t compute_seal() const noexcept;
  bool semantic_valid() const noexcept;
  void seal() noexcept;
  bool authenticated() const noexcept;
  friend class FixedStepIdealGasFlow;
  friend class IdealGasClosureDiagnosticSource;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::IdealGasClosureTestAccess;
#endif
};

enum class IdealGasDiagnosticEntity : std::uint8_t { global_scalar, cell };

class IdealGasClosureDiagnosticSource final {
public:
  ~IdealGasClosureDiagnosticSource() noexcept;
  IdealGasClosureDiagnosticSource(IdealGasClosureDiagnosticSource &&) noexcept;
  IdealGasClosureDiagnosticSource &
  operator=(IdealGasClosureDiagnosticSource &&) = delete;
  IdealGasClosureDiagnosticSource(const IdealGasClosureDiagnosticSource &) =
      delete;
  IdealGasClosureDiagnosticSource &
  operator=(const IdealGasClosureDiagnosticSource &) = delete;

  std::size_t fingerprint_field_count() const;
  std::string_view fingerprint_field_id(std::size_t) const;
  std::string_view fingerprint_field_unit(std::size_t) const;
  IdealGasDiagnosticEntity fingerprint_field_entity(std::size_t) const;
  std::size_t fingerprint_field_item_count(std::size_t) const;
  std::uint64_t fingerprint_field_global_id(std::size_t, std::size_t) const;
  double fingerprint_field_value(std::size_t, std::size_t) const;

  std::size_t sample_field_count() const;
  std::string_view sample_field_id(std::size_t) const;
  std::string_view sample_field_unit(std::size_t) const;
  IdealGasDiagnosticEntity sample_field_entity(std::size_t) const;
  std::size_t sample_field_item_count(std::size_t) const;
  std::uint64_t sample_field_global_id(std::size_t, std::size_t) const;
  double sample_field_value(std::size_t, std::size_t) const;

  std::string_view owned_cell_layout_fingerprint() const;
  std::string_view global_cell_layout_fingerprint() const;
  int relative_rank() const;
  std::uint64_t committed_step() const;
  double committed_time_s() const;
  std::size_t owned_cell_count() const;
  double cell_volume_m3(std::size_t) const;
  IdealGasClosureState closure_state() const;
  const IdealGasStepAttemptReport &report() const;

private:
  struct Impl;
  explicit IdealGasClosureDiagnosticSource(std::unique_ptr<Impl>) noexcept;
  void validate() const;
  double committed_cell_value(runtime::FieldId, std::size_t) const;
  std::unique_ptr<Impl> impl_;
  friend class FixedStepIdealGasFlow;
  friend struct detail::DensityClosureDiagnosticAccess;
};

class FixedStepIdealGasFlow final {
public:
  static FixedStepIdealGasFlow
  create(const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
         const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
         const runtime::MpiContext &, execution::ExecutionContext &,
         runtime::HaloExchange &, const linear::LinearSolver &momentum_solver,
         std::array<linear::Preconditioner *, 3> momentum_preconditioners,
         const linear::LinearSolver &pressure_solver,
         linear::Preconditioner &pressure_preconditioner,
         const runtime::FieldRegistry &, FlowFieldIds,
         MaterialDensityTransportSpec, IdealGasClosure &&);

  ~FixedStepIdealGasFlow() noexcept;
  FixedStepIdealGasFlow(FixedStepIdealGasFlow &&) noexcept;
  FixedStepIdealGasFlow &operator=(FixedStepIdealGasFlow &&) = delete;
  FixedStepIdealGasFlow(const FixedStepIdealGasFlow &) = delete;
  FixedStepIdealGasFlow &operator=(const FixedStepIdealGasFlow &) = delete;

  IdealGasStepAttemptReport
  attempt(FlowState &, double mu, const MomentumTimeStencil &,
          const linear::SolveControl &momentum_control,
          const linear::SolveControl &pressure_control) const;

  IdealGasClosureState closure_state() const;
  MaterialDensityFlowDiagnosticSource
  flow_diagnostic_source(const FlowState &,
                         const IdealGasStepAttemptReport &) const;
  IdealGasClosureDiagnosticSource
  closure_diagnostic_source(const FlowState &,
                            const IdealGasStepAttemptReport &) const;

private:
  struct Impl;
  explicit FixedStepIdealGasFlow(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class IdealGasClosureDiagnosticSource;
  friend struct detail::DensityClosureAdapter;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::IdealGasClosureTestAccess;
#endif
};

} // namespace hundun::flow
