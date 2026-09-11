// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_status.hpp"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace hundun::v04 {

enum class LinearAlgorithm : std::uint8_t { pcg, fgmres, bicgstab };
enum class CouplingKind : std::uint8_t { piso, simple };
enum class MgCorrectionScaling : std::uint8_t {
  residual_minimizing,
  unit_linear
};

enum class GeometryKind : std::uint8_t {
  uniform,
  tensor_stretched,
  coast_runtime_axes_v1
};
enum class TurbulenceKind : std::uint8_t {
  none,
  wale,
  vreman_wall_function,
  vreman
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
  heat_flux_wall,
  // Zero normal gradients and signed outlet-flux mass closure (COAST -2).
  zero_gradient_mass_outlet
};
inline bool is_candidate_transport_outlet(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::pressure_outlet ||
         kind == BoundaryKind::zero_gradient_mass_outlet;
}
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
enum class TimeScheme : std::uint8_t {
  backward_euler,
  variable_bdf2,
  // CN midpoint momentum; BE mass, species and total energy.
  cn_be
};
enum class TransportLaw : std::uint8_t {
  constant,
  sutherland,
  coast_native_air,
  coast_perry
};
enum class ImmersedFluidSide : std::uint8_t { outside, inside };
enum class IbmReconstructionPolicy : std::uint8_t {
  strict_quadratic,
  adaptive_order
};

enum class ReactionMode : std::uint8_t {
  none,
  finite_rate_mean,
  pasr_algebraic_v1,
  esf_tpdf
};

enum class TcrMode : std::uint8_t { off, shadow, experimental, validated };
struct TcrSpec {
  TcrMode mode{TcrMode::off};
  std::vector<std::string> reactants;
  std::vector<double> progress_weights;
  int initialization_sign{};
  double weak_rate_threshold{1e-12};
};
struct EsfSpec {
  std::uint32_t fields{2};
  std::uint64_t seed{};
  // Field-major independent-species offsets about the gas initial mean.
  // Empty initializes every stochastic field to that mean.
  std::vector<double> initial_species_offsets;
  TcrSpec tcr;
};

struct ReactionSpec {
  ReactionMode mode{ReactionMode::none};
  std::string mechanism_sha256;
  std::string phase;
  enum class Representation : std::uint8_t {
    external_provider, analytic_isomer, direct_cantera
  };
  Representation representation{Representation::external_provider};
  std::filesystem::path mechanism_file;
  double analytic_rate_s{2.0};
  double analytic_cp_j_per_kg_k{1000.0};
  double relative_tolerance{1e-8};
  double absolute_tolerance{1e-14};
  std::uint32_t maximum_internal_steps{2000};
  double mixing_c_z{1.0};
  double turbulent_schmidt{0.7};
  std::optional<EsfSpec> esf;
};

struct SprayInjectionSpec {
  std::uint64_t id{};
  Real3 origin_m{};
  Real3 axis{1.0, 0.0, 0.0};
  double cone_half_angle_rad{}; // zero is a point injector
  double speed_m_per_s{};
  double mass_flow_rate_kg_per_s{};
  double represented_mass_per_parcel_kg{};
  double droplet_diameter_m{};
  double temperature_k{};
};

struct SpraySpec {
  std::filesystem::path liquid_file;
  std::uint64_t liquid_fingerprint{}; // FNV-1a64 of exact liquid asset bytes
  std::uint64_t seed{};
  std::uint32_t maximum_local_parcels{1024};
  std::uint32_t maximum_local_segments{8192};
  double maximum_substep_s{1e-4};
  double minimum_substep_s{1e-12};
  double relative_tolerance{1e-6};
  bool tab_breakup{};
  std::vector<SprayInjectionSpec> injectors;
};

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
  std::filesystem::path axes_file;
  std::array<std::vector<double>, 3U> coast_runtime_faces;
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

struct ImmersedBoundarySpec {
  std::filesystem::path stl_file;
  ImmersedFluidSide fluid_side{ImmersedFluidSide::outside};
  IbmReconstructionPolicy reconstruction_policy{
      IbmReconstructionPolicy::strict_quadratic};
  // Optional frozen binary 0=solid, 1=fluid Cartesian cell authority.
  // This selects explicit Cartesian face wall geometry instead of STL scans.
  std::optional<std::filesystem::path> marker_file;
  PlanFingerprint marker_fingerprint{};
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
  // Molecular and turbulent scalar diffusivities are formed from the
  // corresponding dynamic viscosity divided by these dimensionless
  // Schmidt numbers.  JSON input requires both values explicitly; the
  // positive defaults keep programmatic typed fixtures well formed.
  double molecular_schmidt{1.0};
  double turbulent_schmidt{1.0};
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

struct PatchInletSpec {
  std::int32_t label{};
  // x_min, x_max, y_min, y_max, z_min, z_max; for an immersed patch
  // this is the direction from its fluid owner toward the solid donor.
  std::uint8_t face{};
  bool immersed{};
  BoundaryFaceSpec boundary;
};

struct PatchInletsSpec {
  // Dense little-endian int32 labels, global owned cells in x-fast order.
  // Positive labels select fluid owners; immersed sources require the
  // corresponding negative label in the immediately adjacent solid cell.
  std::filesystem::path labels_file;
  PlanFingerprint labels_fingerprint{};
  std::vector<PatchInletSpec> patches;
};

struct SchemeSpec {
  ConvectionScheme momentum{ConvectionScheme::central2};
  ConvectionScheme enthalpy{ConvectionScheme::limited_central2};
  ConvectionScheme species{ConvectionScheme::tvd2};
  ConvectionScheme passive_scalar{ConvectionScheme::tvd2};
  DiffusionScheme diffusion{DiffusionScheme::central2};
  double limiter{1.0};
};

struct PressureLinearSolverSpec {
  double absolute_tolerance{1.0e-13};
  double relative_tolerance{1.0e-13};
  std::uint32_t maximum_iterations{400U};
  std::uint32_t true_residual_interval{4U};
  std::uint32_t krylov_restart{12U};
  LinearAlgorithm algorithm{LinearAlgorithm::fgmres};
  MgCorrectionScaling mg_correction_scaling{
      MgCorrectionScaling::residual_minimizing};
};

struct SolverToleranceSpec {
  double eos{1.0e-10};
  double continuity{1.0e-10};
  double closed_mass{1.0e-10};
  double gauge{1.0e-10};
};

// COAST ewt/ewt_pdf: freeze max accepted conservative variables per attempt,
// divide by this physical reference time, and test equation residuals.
struct ColdStoppingSpec {
  double reference_time{};
  double momentum{1.0e-4};
  double enthalpy{1.0e-4};
  double species{1.0e-4};

  bool valid() const noexcept {
    return std::isfinite(reference_time) && reference_time > 0.0 &&
           std::isfinite(momentum) && momentum > 0.0 && momentum < 1.0 &&
           std::isfinite(enthalpy) && enthalpy > 0.0 && enthalpy < 1.0 &&
           std::isfinite(species) && species > 0.0 && species < 1.0;
  }
};

struct SolverSpec {
  CouplingKind coupling{CouplingKind::piso};
  PressureLinearSolverSpec pressure;
  SolverToleranceSpec terminal;
  std::optional<ColdStoppingSpec> cold_stopping;
};

struct TimeControlSpec {
  TimeControlKind control{TimeControlKind::adaptive_flow};
  TimeScheme scheme{TimeScheme::cn_be};
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
  double critical_temperature{};  // K, COAST/Perry corresponding states.
  double critical_pressure{};     // atm.
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
  SolverSpec solver;
  SchemeSpec schemes;
  TimeControlSpec time;
  ThermophysicalSpec thermophysics;
  std::vector<TransportedScalarSpec> transported_scalars;
  std::vector<std::filesystem::path> data_files;
  std::optional<ImmersedBoundarySpec> immersed_boundary;
  ReactionSpec reaction;
  std::optional<SpraySpec> spray;
  std::optional<PatchInletsSpec> patch_inlets;
  PlanFingerprint fingerprint{};
  // Same physical case with variable_bdf2, for explicit history rebuild only.
  // Zero for methods which do not migrate from that time scheme.
  PlanFingerprint legacy_time_fingerprint{};
};

class CaseCompiler {
 public:
  static Status load_and_compile(MPI_Comm communicator,
                                 const std::filesystem::path& case_root,
                                 ValidatedModel& out);
  // Re-read both cold cases for Sutherland -> Perry transport, or a Perry
  // zero-gradient outlet -> static-pressure outlet with backflow transition.
  // The selected change preserves every other control and referenced byte.
  static Status validate_transport_change(MPI_Comm communicator,
      const std::filesystem::path& source_root, const ValidatedModel& source,
      const std::filesystem::path& target_root, const ValidatedModel& target);
};

}  // namespace hundun::v04
