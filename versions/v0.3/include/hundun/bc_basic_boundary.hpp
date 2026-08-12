// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_view.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::boundary {

enum class BoundaryKind : std::uint8_t {
  periodic,
  no_slip_wall,
  symmetry,
  velocity_inlet,
  pressure_outlet
};

enum class VelocityRule : std::uint8_t {
  periodic_pair,
  prescribed_zero,
  reflect_normal_copy_tangential,
  prescribed_inlet,
  pure_outflow
};

enum class PressureRule : std::uint8_t {
  periodic_pair,
  zero_normal_gradient,
  prescribed_value
};

enum class TransportRule : std::uint8_t {
  periodic_pair,
  copy_interior,
  zero_normal_diffusive_flux,
  prescribed_value,
  pure_outflow
};

enum class MassFluxRule : std::uint8_t {
  periodic_pair,
  identically_zero,
  prescribed_inlet_state,
  outflow_only
};

struct ResolvedInletState final {
  runtime::Real3 velocity_m_per_s;
  double density_kg_per_m3;
  double enthalpy_J_per_kg;
  std::optional<double> temperature_K;
  std::vector<double> scalar_values;
};

struct VelocityBoundaryValues final {
  runtime::Real3 face;
  runtime::Real3 exterior;
};

struct ScalarBoundaryValues final {
  double face;
  double exterior;
};

class BoundaryDescriptor final {
public:
  BoundaryDescriptor(const BoundaryDescriptor &) = delete;
  BoundaryDescriptor &operator=(const BoundaryDescriptor &) = delete;
  BoundaryDescriptor(BoundaryDescriptor &&) noexcept;
  BoundaryDescriptor &operator=(BoundaryDescriptor &&) noexcept;

  std::uint32_t stable_id() const noexcept;
  std::string_view name() const noexcept;
  BoundaryKind kind() const noexcept;
  VelocityRule velocity_rule() const noexcept;
  PressureRule pressure_rule() const noexcept;
  TransportRule density_rule() const noexcept;
  TransportRule enthalpy_rule() const noexcept;
  TransportRule scalar_rule() const noexcept;
  MassFluxRule mass_flux_rule() const noexcept;
  std::optional<std::uint32_t> paired_patch_id() const noexcept;
  const std::optional<ResolvedInletState> &inlet_state() const noexcept;
  std::optional<double> pressure_value_pa() const noexcept;

private:
  friend class BoundaryRegistry;

  BoundaryDescriptor(std::uint32_t stable_id, std::string name,
                     BoundaryKind kind, VelocityRule velocity_rule,
                     PressureRule pressure_rule, TransportRule density_rule,
                     TransportRule enthalpy_rule, TransportRule scalar_rule,
                     MassFluxRule mass_flux_rule,
                     std::optional<std::uint32_t> paired_patch_id,
                     std::optional<ResolvedInletState> inlet_state,
                     std::optional<double> pressure_value_pa);

  std::uint32_t stable_id_{};
  std::string name_;
  BoundaryKind kind_{BoundaryKind::periodic};
  VelocityRule velocity_rule_{VelocityRule::periodic_pair};
  PressureRule pressure_rule_{PressureRule::periodic_pair};
  TransportRule density_rule_{TransportRule::periodic_pair};
  TransportRule enthalpy_rule_{TransportRule::periodic_pair};
  TransportRule scalar_rule_{TransportRule::periodic_pair};
  MassFluxRule mass_flux_rule_{MassFluxRule::periodic_pair};
  std::optional<std::uint32_t> paired_patch_id_;
  std::optional<ResolvedInletState> inlet_state_;
  std::optional<double> pressure_value_pa_;
};

enum class FinalFluxDecision : std::uint8_t { admissible, outlet_backflow };

struct OutletBackflowEvidence final {
  std::uint32_t patch_id;
  std::uint64_t step;
  double time_s;
  double minimum_outward_mass_flux_kg_per_s;
  mesh::GlobalFaceId global_face_id;
  int lowest_failing_rank;
};

struct FinalFluxAdmissibility final {
  FinalFluxDecision decision;
  std::optional<OutletBackflowEvidence> evidence;
};

class BoundaryRegistry final {
public:
  static BoundaryRegistry create(const config::FlowCaseConfig &resolved,
                                 const mesh::MeshTopology &topology);

  ~BoundaryRegistry() noexcept;
  BoundaryRegistry(BoundaryRegistry &&) noexcept;
  BoundaryRegistry &operator=(BoundaryRegistry &&) = delete;
  BoundaryRegistry(const BoundaryRegistry &) = delete;
  BoundaryRegistry &operator=(const BoundaryRegistry &) = delete;

  const BoundaryDescriptor &patch(std::uint32_t stable_id) const;
  std::size_t scalar_count() const noexcept;
  std::string_view scalar_name(std::size_t scalar) const;
  bool open_domain() const noexcept;
  std::optional<std::uint32_t> velocity_inlet_patch_id() const noexcept;
  std::optional<std::uint32_t> pressure_outlet_patch_id() const noexcept;

  VelocityBoundaryValues
  evaluate_velocity(std::uint32_t patch_id, runtime::Real3 interior_velocity,
                    runtime::Real3 owner_outward_area_vector) const;
  ScalarBoundaryValues evaluate_pressure(std::uint32_t patch_id,
                                         double interior_pi) const;
  ScalarBoundaryValues evaluate_density(std::uint32_t patch_id,
                                        double interior_density) const;
  ScalarBoundaryValues evaluate_enthalpy(std::uint32_t patch_id,
                                         double interior_enthalpy) const;
  ScalarBoundaryValues evaluate_scalar(std::uint32_t patch_id,
                                       std::size_t scalar,
                                       double interior_value) const;

  FinalFluxAdmissibility assess_final_pressure_outlet_flux(
      const mesh::MeshTopology &topology, const runtime::MpiContext &mpi,
      const runtime::FaceFieldView<const double> &final_face_mass_flux,
      std::uint64_t step, double time_s) const;

private:
  struct Impl;
  explicit BoundaryRegistry(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
};

} // namespace hundun::boundary
