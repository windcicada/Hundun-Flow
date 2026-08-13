// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_boundary.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::v04 {

inline constexpr double kUniversalGasConstant = 8314.46261815324;

enum class ThermodynamicsKernel : std::uint8_t { nasa7, constant_cp };

struct ThermalRevisionTuple {
  RevisionToken enthalpy{};
  RevisionToken composition{};
  StorageIdentity enthalpy_storage{};
  StorageIdentity composition_storage{};
  RevisionDomainIdentity enthalpy_revision_domain{};
  RevisionDomainIdentity composition_revision_domain{};
  PlanFingerprint patch_identity{};
};

class ThermalState {
 public:
  ThermalState() noexcept = default;

  double temperature() const noexcept { return temperature_; }
  double cp() const noexcept { return cp_; }
  double gas_constant() const noexcept { return gas_constant_; }
  double gamma() const noexcept { return gamma_; }
  double drho_dp_hY() const noexcept { return drho_dp_hY_; }
  double sound_speed() const noexcept { return sound_speed_; }
  ThermalRevisionTuple revisions() const noexcept { return revisions_; }

 private:
  friend class ThermodynamicsPlan;
  friend class ClosedMassPlan;
  double temperature_{};
  double cp_{};
  double gas_constant_{};
  double gamma_{};
  double drho_dp_hY_{};
  double sound_speed_{};
  PlanFingerprint thermodynamics_{};
  ThermalRevisionTuple revisions_{};
  bool certified_{};
};

struct PressureThermoState {
  double rho{};
  double drho_dp_hY{};
};

struct ThermoState {
  double rho{};
  double temperature{};
  double cp{};
  double gas_constant{};
  double gamma{};
  double drho_dp_hY{};
  double sound_speed{};
  double mach{};
};

class ThermodynamicsPlan {
 public:
  ThermodynamicsPlan() noexcept = default;
  ThermodynamicsPlan(const ThermodynamicsPlan&) = delete;
  ThermodynamicsPlan& operator=(const ThermodynamicsPlan&) = delete;
  ThermodynamicsPlan(ThermodynamicsPlan&&) noexcept = default;
  ThermodynamicsPlan& operator=(ThermodynamicsPlan&&) noexcept = default;

  static Status compile(const ThermophysicalSpec& spec,
                        Span<const TransportedScalarSpec> scalar_catalog,
                        ThermodynamicsPlan& out) noexcept;

  Status evaluate(double p_abs, double h,
                  Span<const double> independent_mass_fractions,
                  Real3 velocity, ThermoState& out,
                  double temperature_hint =
                      std::numeric_limits<double>::quiet_NaN()) const noexcept;
  Status evaluate_thermal(
      double h, Span<const double> independent_mass_fractions,
      ThermalRevisionTuple revisions, ThermalState& out,
      double temperature_hint =
          std::numeric_limits<double>::quiet_NaN()) const noexcept;
  Status complete_state(double p_abs, const ThermalState& thermal,
                        ThermalRevisionTuple revisions, Real3 velocity,
                        ThermoState& out) const noexcept;
  Status evaluate_pressure(double p_abs, const ThermalState& thermal,
                           ThermalRevisionTuple revisions,
                           PressureThermoState& out) const noexcept;
  Status evaluate_from_reference_pressure(
      double p_ref, double pi, double h,
      Span<const double> independent_mass_fractions, Real3 velocity,
      ThermoState& out,
      double temperature_hint =
          std::numeric_limits<double>::quiet_NaN()) const noexcept;
  Status mixture_enthalpy(double temperature,
                          Span<const double> independent_mass_fractions,
                          double& enthalpy, double& cp,
                          double& gas_constant) const noexcept;

  std::size_t species_count() const noexcept {
    return inverse_molecular_weight_.size();
  }
  std::size_t independent_species_count() const noexcept {
    return independent_to_species_.size();
  }
  std::size_t dependent_species_index() const noexcept {
    return dependent_species_;
  }
  ThermodynamicsKernel kernel() const noexcept { return kernel_; }
  double minimum_temperature() const noexcept {
    return minimum_temperature_;
  }
  double maximum_temperature() const noexcept {
    return maximum_temperature_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class TransportPlan;
  friend class ClosedMassPlan;
  Status evaluate_thermal_impl(
      double h, Span<const double> independent_mass_fractions,
      ThermalRevisionTuple revisions, bool require_certificate,
      ThermalState& out, double temperature_hint) const noexcept;
  Status evaluate_pressure_impl(double p_abs, const ThermalState& thermal,
                                ThermalRevisionTuple revisions,
                                bool require_certificate,
                                PressureThermoState& out) const noexcept;
  Status complete_state_impl(double p_abs, const ThermalState& thermal,
                             ThermalRevisionTuple revisions,
                             bool require_certificate, Real3 velocity,
                             ThermoState& out) const noexcept;
  Status composition(Span<const double> independent_mass_fractions,
                     double& dependent) const noexcept;
  Status mixture_properties(double temperature,
                            Span<const double> independent_mass_fractions,
                            double dependent, double& enthalpy, double& cp,
                            double& gas_constant) const noexcept;

  std::vector<double> inverse_molecular_weight_;
  std::array<std::vector<double>, 7U> nasa_low_;
  std::array<std::vector<double>, 7U> nasa_high_;
  std::vector<double> temperature_switch_;
  std::vector<std::uint16_t> independent_to_species_;
  std::size_t dependent_species_{};
  double minimum_temperature_{};
  double maximum_temperature_{};
  double relative_tolerance_{};
  std::uint32_t maximum_iterations_{};
  double constant_enthalpy_offset_{};
  PlanFingerprint source_fingerprint_{};
  ThermodynamicsKernel kernel_{ThermodynamicsKernel::nasa7};
  PlanFingerprint fingerprint_{};
};

enum class TransportKernel : std::uint8_t { constant, sutherland_wilke };

struct MolecularTransportState {
  double viscosity{};
  double conductivity{};
};

class TransportPlan {
 public:
  TransportPlan() noexcept = default;
  TransportPlan(const TransportPlan&) = delete;
  TransportPlan& operator=(const TransportPlan&) = delete;
  TransportPlan(TransportPlan&&) noexcept = default;
  TransportPlan& operator=(TransportPlan&&) noexcept = default;

  static Status compile(const ThermophysicalSpec& spec,
                        const ThermodynamicsPlan& thermodynamics,
                        TransportPlan& out) noexcept;
  Status evaluate(double temperature,
                  Span<const double> independent_mass_fractions,
                  MolecularTransportState& out) const noexcept;

  TransportKernel kernel() const noexcept { return kernel_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  std::vector<double> molecular_weight_;
  std::vector<double> viscosity_reference_;
  std::vector<double> reference_temperature_;
  std::vector<double> sutherland_temperature_;
  std::vector<double> prandtl_;
  std::vector<double> conductivity_;
  std::array<std::vector<double>, 7U> nasa_low_;
  std::array<std::vector<double>, 7U> nasa_high_;
  std::vector<double> temperature_switch_;
  std::vector<double> molecular_weight_ratio_quarter_;
  std::vector<double> wilke_denominator_reciprocal_;
  std::vector<double> constant_wilke_phi_;
  std::vector<std::uint16_t> independent_to_species_;
  std::size_t dependent_species_{};
  double minimum_temperature_{};
  double maximum_temperature_{};
  TransportKernel kernel_{TransportKernel::constant};
  PlanFingerprint fingerprint_{};
};

class ThermophysicalCompiler {
 public:
  static Status load_and_compile(
      MPI_Comm communicator, const ValidatedModel& model,
      ThermophysicalSpec& spec,
      ThermodynamicsPlan& thermodynamics,
      TransportPlan& transport,
      struct ThermophysicalCompileDiagnostics* diagnostics = nullptr) noexcept;
};

struct ThermophysicalCompileDiagnostics {
  int lowest_failing_rank{-1};
};

struct VelocityGradient {
  std::array<double, 9U> value{};
};

struct DerivedRevisionTuple {
  RevisionToken velocity{};
  RevisionToken geometry{};
  RevisionToken boundary{};
  RevisionToken turbulence{};
  StorageIdentity velocity_storage{};
  RevisionDomainIdentity velocity_revision_domain{};
  PlanFingerprint patch_identity{};
  PlanFingerprint geometry_plan_identity{};
  PlanFingerprint boundary_plan_identity{};
  PlanFingerprint turbulence_plan_identity{};
};

DerivedRevisionTuple make_derived_revision_tuple(
    ConstFieldView velocity, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, RevisionToken boundary,
    PlanFingerprint boundary_plan_identity, RevisionToken turbulence,
    PlanFingerprint turbulence_plan_identity) noexcept;

class DerivedFieldPlan {
 public:
  DerivedFieldPlan() noexcept = default;
  DerivedFieldPlan(const DerivedFieldPlan&) = delete;
  DerivedFieldPlan& operator=(const DerivedFieldPlan&) = delete;
  DerivedFieldPlan(DerivedFieldPlan&&) noexcept = default;
  DerivedFieldPlan& operator=(DerivedFieldPlan&&) noexcept = default;

  static Status compile(FieldId velocity_field,
                        Span<const FieldId> declared_fields,
                        RevisionSlotId gradient_cache_slot,
                        RevisionSourceId geometry_source,
                        RevisionSourceId boundary_source,
                        RevisionSourceId turbulence_source,
                        std::size_t cell_capacity,
                        DerivedFieldPlan& out) noexcept;
  Status prepare_velocity_gradient(
      ConstFieldView velocity, const CartesianGeometryPlan& geometry,
      const MeshPatch& patch,
      DerivedRevisionTuple revisions, AttemptTransaction& transaction,
      Span<const VelocityGradient>& out) noexcept;
  Status finalize_velocity_gradient(
      const AttemptTransaction& transaction) noexcept;
  Status velocity_gradient(DerivedRevisionTuple revisions,
                           Span<const VelocityGradient>& out) const noexcept;

  std::uint64_t gradient_compute_count() const noexcept {
    return gradient_compute_count_;
  }
  RevisionToken cache_revision() const noexcept { return cache_revision_; }

 private:
  std::vector<VelocityGradient> active_gradient_;
  std::vector<VelocityGradient> pending_gradient_;
  DerivedRevisionTuple active_revisions_{};
  DerivedRevisionTuple pending_revisions_{};
  FieldId velocity_field_{};
  RevisionSlotId cache_slot_{};
  RevisionSourceId geometry_source_{};
  RevisionSourceId boundary_source_{};
  RevisionSourceId turbulence_source_{};
  RevisionToken cache_revision_{};
  RevisionToken pending_cache_revision_{};
  RevisionToken next_cache_revision_{1U};
  std::uint64_t pending_attempt_identity_{};
  std::uint64_t gradient_compute_count_{};
  bool active_valid_{};
  bool pending_valid_{};
};

class ContributionRegistry;

class EffectiveViscosityAuthority {
 public:
  Status claim(FieldId output, StageId stage,
               ContributionRegistry& registry) noexcept;
  bool claimed() const noexcept { return claimed_; }
  FieldId output() const noexcept { return output_; }
  StageId stage() const noexcept { return stage_; }

 private:
  FieldId output_{};
  StageId stage_{};
  bool claimed_{};
};

struct UnitDimension {
  std::array<std::int8_t, 7U> si_exponents{};
  friend bool operator==(UnitDimension left,
                         UnitDimension right) noexcept {
    return left.si_exponents == right.si_exponents;
  }
};

enum class ContributionCapability : std::uint8_t {
  inert_source,
  chemistry,
  reacting
};

struct ContributionSpec {
  FieldId conserved_quantity{};
  UnitDimension units{};
  StageId stage{};
  Span<const FieldId> reads{};
  FieldId explicit_source{};
  FieldId implicit_diagonal{};
  bool supplies_implicit_diagonal{};
  ContributionCapability capability{ContributionCapability::inert_source};
};

struct CompiledContribution {
  FieldId conserved_quantity{};
  UnitDimension units{};
  StageId stage{};
  std::uint32_t read_begin{};
  std::uint16_t read_count{};
  FieldId explicit_source{};
  FieldId implicit_diagonal{};
  bool supplies_implicit_diagonal{};
  std::uint32_t registration_ordinal{};
};

class ContributionRegistry {
 public:
  Status configure(Span<const FieldId> declared_fields) noexcept;
  Status register_contribution(const ContributionSpec& spec) noexcept;
  Status freeze() noexcept;
  Span<const CompiledContribution> contributions() const noexcept {
    return {contributions_.data(), contributions_.size()};
  }
  Span<const FieldId> reads() const noexcept {
    return {reads_.data(), reads_.size()};
  }
  bool frozen() const noexcept { return frozen_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class EffectiveViscosityAuthority;
  std::vector<CompiledContribution> contributions_;
  std::vector<FieldId> reads_;
  std::vector<FieldId> declared_fields_;
  std::uint32_t next_ordinal_{};
  PlanFingerprint fingerprint_{};
  bool frozen_{};
  bool effective_viscosity_claimed_{};
  FieldId effective_viscosity_output_{};
  StageId effective_viscosity_stage_{};
};

struct ClosedMassCellView {
  Span<const double> pressure_perturbation;
  Span<const double> enthalpy;
  // Species-major SoA: independent species s starts at s*cell_count.
  Span<const double> independent_mass_fractions;
  Span<const double> volume;
  Span<const std::uint8_t> active;
};

struct ClosedMassResult {
  double pressure_reference{};
  double mass{};
  double residual{};
  std::uint32_t iterations{};
  int lowest_failing_rank{-1};
};

class ClosedMassPlan {
 public:
  static Status compile(PressureReferenceKind authority,
                        const ThermophysicalSpec& spec,
                        ClosedMassPlan& out) noexcept;
  Status solve(MPI_Comm communicator,
               const ThermodynamicsPlan& thermodynamics,
               const ClosedMassCellView& cells, double target_mass,
               double current_pressure_reference,
               ClosedMassResult& out) const noexcept;
  PressureReferenceKind authority() const noexcept { return authority_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  PressureReferenceKind authority_{PressureReferenceKind::boundary_absolute};
  double relative_tolerance_{};
  double maximum_relative_step_{};
  std::uint32_t maximum_iterations_{};
  PlanFingerprint fingerprint_{};
};

}  // namespace hundun::v04
