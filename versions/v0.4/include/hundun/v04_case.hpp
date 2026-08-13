// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_status.hpp"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace hundun::v04 {

enum class GeometryKind : std::uint8_t { uniform, tensor_stretched };
enum class TurbulenceKind : std::uint8_t {
  none,
  wale,
  vreman_wall_function
};
enum class TimeControlKind : std::uint8_t {
  fixed,
  adaptive_flow,
  adaptive_acoustic
};
enum class PressureReferenceKind : std::uint8_t {
  boundary_absolute,
  closed_mass
};
enum class CartesianFace : std::uint8_t {
  x_min,
  x_max,
  y_min,
  y_max,
  z_min,
  z_max
};
enum class BoundaryKind : std::uint8_t {
  none,
  velocity_inlet,
  mass_flow_inlet,
  static_state_inlet,
  total_state_inlet,
  pressure_outlet,
  nscbc_inlet,
  nscbc_outlet,
  no_slip_wall,
  moving_wall,
  slip,
  symmetry,
  periodic,
  adiabatic_wall,
  isothermal_wall,
  heat_flux_wall
};
enum class ScalarBoundaryKind : std::uint8_t {
  dirichlet,
  normal_flux,
  zero_gradient,
  convective
};
enum class TransportedScalarRole : std::uint8_t { species, passive_scalar };
enum class ConvectionScheme : std::uint8_t {
  central2,
  limited_central2,
  tvd2
};
enum class DiffusionScheme : std::uint8_t { central2 };
enum class TimeScheme : std::uint8_t { backward_euler, variable_bdf2 };
enum class TransportLaw : std::uint8_t { constant, sutherland };

struct CaseSpec {
  std::filesystem::path root;
};

struct FocusRegionSpec {
  Real3 lower{};
  Real3 upper{};
  Real3 target_spacing{};
};

struct MeshLimits {
  std::uint64_t max_global_cells{};
  std::uint64_t max_memory_bytes_per_rank{};
};

struct CartesianMeshSpec {
  GeometryKind kind{GeometryKind::uniform};
  Real3 lower{};
  Real3 upper{};
  bool has_exact_cells{};
  Int3 exact_cells{};
  bool has_base_spacing{};
  Real3 base_spacing{};
  Real3 minimum_spacing{};
  double max_growth_ratio{1.0};
  std::vector<FocusRegionSpec> focus_regions;
  MeshLimits limits{};
};

struct ScalarBoundarySpec {
  std::string stable_name;
  ScalarBoundaryKind kind{ScalarBoundaryKind::zero_gradient};
  double value{};
  // Used only when an outlet explicitly permits inflow. The ordinary outlet
  // closure remains independent so the hot resolver can select by face-cell
  // flow direction without reparsing case data.
  ScalarBoundaryKind backflow_kind{ScalarBoundaryKind::zero_gradient};
  double backflow_value{};
};

struct TransportedScalarSpec {
  std::string stable_name;
  TransportedScalarRole role{TransportedScalarRole::passive_scalar};
};

struct BoundaryFaceSpec {
  BoundaryKind flow_kind{BoundaryKind::symmetry};
  BoundaryKind thermal_kind{BoundaryKind::none};
  Real3 velocity{};
  Real3 direction{};
  Real3 backflow_velocity{};
  double mass_flow_rate{};
  double pressure{};
  double temperature{};
  double total_pressure{};
  double total_temperature{};
  double backflow_temperature{};
  double heat_flux{};
  double relaxation{};
  double mach_limit{0.95};
  bool allow_backflow{};
  std::vector<ScalarBoundarySpec> scalars;
};

struct SchemeSpec {
  ConvectionScheme momentum{ConvectionScheme::limited_central2};
  ConvectionScheme enthalpy{ConvectionScheme::limited_central2};
  ConvectionScheme species{ConvectionScheme::tvd2};
  ConvectionScheme passive_scalar{ConvectionScheme::tvd2};
  DiffusionScheme diffusion{DiffusionScheme::central2};
  double limiter{1.0};
};

struct TimeControlSpec {
  TimeControlKind control{TimeControlKind::adaptive_flow};
  TimeScheme scheme{TimeScheme::variable_bdf2};
  double initial_dt{1.0e-4};
  double minimum_dt{1.0e-10};
  double maximum_dt{1.0};
  double convective_cfl{0.8};
  double viscous_cfl{0.5};
  double thermal_cfl{0.5};
  double species_cfl{0.5};
  double acoustic_cfl{0.8};
  double maximum_growth{1.25};
  double retry_factor{0.5};
  std::uint32_t maximum_retries{8U};
  double minimum_bdf_ratio{0.2};
  double maximum_bdf_ratio{5.0};
};

struct SpeciesThermophysicalSpec {
  std::string stable_name;
  double molecular_weight{};
  double temperature_switch{};
  std::array<double, 7U> nasa7_low{};
  std::array<double, 7U> nasa7_high{};
  TransportLaw transport_law{TransportLaw::constant};
  double viscosity_reference{};
  double transport_reference_temperature{};
  double sutherland_temperature{};
  double prandtl{};
  double conductivity{};
};

struct ThermophysicalSpec {
  std::filesystem::path data_file;
  double minimum_temperature{};
  double maximum_temperature{};
  double temperature_relative_tolerance{};
  std::uint32_t maximum_temperature_iterations{};
  double closed_mass_relative_tolerance{};
  std::uint32_t maximum_closed_mass_iterations{};
  double maximum_closed_mass_relative_step{};
  std::vector<SpeciesThermophysicalSpec> species;
};

struct ValidatedModel {
  CartesianMeshSpec mesh;
  TurbulenceKind turbulence{TurbulenceKind::vreman_wall_function};
  PressureReferenceKind pressure_reference{
      PressureReferenceKind::boundary_absolute};
  std::array<BoundaryFaceSpec, 6U> boundaries;
  SchemeSpec schemes;
  TimeControlSpec time;
  ThermophysicalSpec thermophysics;
  std::vector<TransportedScalarSpec> transported_scalars;
  std::vector<std::filesystem::path> data_files;
  std::optional<std::filesystem::path> stl_file;
  PlanFingerprint fingerprint{};
};

class CaseCompiler {
 public:
  static Status load_and_compile(MPI_Comm communicator,
                                 const std::filesystem::path& case_root,
                                 ValidatedModel& out);
};

}  // namespace hundun::v04
