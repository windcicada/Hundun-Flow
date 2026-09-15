// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_combustion.hpp"
#include "hundun/v04_portable.hpp"
#include <array>
#include <vector>

namespace hundun::v04::esf::detail {
// HUNDUN Philox v1 address: no cell/rank/species/retry cursor.
struct CounterAddress {
  std::uint64_t seed{}, accepted_step{};
  std::uint32_t stochastic_stage{}, field_pair{}, spatial_direction{},
      purpose{1};
};
std::array<std::uint32_t, 4>
    philox_words(std::array<std::uint32_t, 4>,
                 std::array<std::uint32_t, 2>) noexcept;
std::array<std::uint32_t, 4> philox(const CounterAddress &,
                                    std::uint32_t lane) noexcept;
struct WienerReport {
  portable::Status status{portable::Status::invalid_input};
  std::array<std::array<double, 3>, 4> increments{};
};
WienerReport balanced_wiener(std::size_t fields, double dt,
                             const CounterAddress &) noexcept;
// Full species; fluxes are kg/(m2 s), correction J_i-Y_i sum(J).
portable::Status correct_species_flux(const double *y, const double *raw,
                                      std::size_t species,
                                      double *corrected) noexcept;
class Workspace;
struct View {
  portable::Revision revision{};
  std::uint64_t composition_fingerprint{};
  std::size_t fields{}, species{};
  // Field-major tuples [Y_0,...Y_Ns-1,h_tc(J/kg)].
  const double *values{};
  // Borrowed until next advance (even failed) or owner destruction.
  const Workspace *owner{};
  std::uint64_t generation{};
};
// Equal-weight Favre samples define physical Y/h; the independently
// transported auxiliary tuple defines pressure-coupling thermodynamics.
// Auxiliary species can be signed; their sum remains one. They are kept
// separate from the physical mass-fraction simplex.
// Densities are EOS results for these same tuples, pressure and revision.
// The reducer preserves every input and publishes caller-owned moments only
// after validating the complete request. It performs no state correction.
struct DualStateRequest {
  View stochastic{}, auxiliary{}; // auxiliary.fields == 1
  portable::Revision expected_revision{};
  const double* stochastic_densities_kg_per_m3{};
  double auxiliary_density_kg_per_m3{};
};
struct DualStateOutput {
  double* physical_mean{}; // Ns full mass fractions followed by h [J/kg]
  double* physical_variance{};
  std::size_t capacity{};
};
struct DualStateReport {
  portable::Status status{portable::Status::invalid_input};
  portable::Revision revision{};
  std::uint64_t composition_fingerprint{};
  std::uint64_t model_identity{0x4553464455410001ULL};
  double statistical_density_kg_per_m3{}; // 1 / mean(1 / rho_f)
  double pressure_density_kg_per_m3{}; // auxiliary EOS
  double max_species_mean_gap{};
  double enthalpy_mean_gap_j_per_kg{}; // physical mean minus auxiliary
};
DualStateReport dual_state_moments(const DualStateRequest&,
                                  DualStateOutput) noexcept;
// COAST evaluates auxiliary EOS with positive species weights, retaining
// the raw field0 tuple. A normalized PH query uses Y_s=max(Y0_s,0)/M and
// h=h0/M, M=sum(max(Y0_s,0)); its density divided by M is the auxiliary
// pressure density. This is an algebraic coordinate change for the same
// extensive mixture EOS. Physical stochastic fields retain their own EOS.
struct AuxiliaryCoordinatesReport {
  portable::Status status{portable::Status::invalid_input};
  double positive_weight_sum{};
  double added_positive_weight{};
  double query_enthalpy_j_per_kg{};
  portable::Revision revision{};
};
AuxiliaryCoordinatesReport auxiliary_eos_coordinates(
    const View& raw, portable::Revision expected_revision,
    double* query_mass_fractions, std::size_t capacity) noexcept;
struct AuxiliaryPressureState {
  portable::Revision revision{};
  std::uint64_t composition_fingerprint{};
  double density_kg_per_m3{};
  double density_pressure_derivative{}; // fixed raw h0 and Y0
  double density_enthalpy_derivative{}; // fixed p and raw Y0; h0 in J/kg
};
// Derivatives of the positive-weight ideal-mixture EOS in the raw field0
// coordinates. Mapping h0 to a physical-mean energy correction is a separate
// scheduling decision; this query does not assume their increments coincide.
portable::Status auxiliary_pressure_state(const AuxiliaryCoordinatesReport&,
    const portable::GasSample& normalized, portable::Revision expected_revision,
    std::uint64_t expected_composition, AuxiliaryPressureState&) noexcept;
struct Request {
  View accepted{};
  portable::Revision expected_revision{};
  const combustion::ChemistryIdentity *identity{};
  double dt_s{}, mixing_time_s{}, tcr_control{1}, turbulent_diffusivity_m2_s{};
  // Already assembled deterministic rates; molecular diffusion only here.
  const double *deterministic_rates{};
  // Field-major/component-major xyz gradients; never includes diffusivity.
  const double *gradients{};
  CounterAddress random{};
};
// IEM contribution to a finite-volume scalar equation before multiplication
// by cell volume: beta [kg/(m3 s)] and beta*mean [kg/(m3 s) * scalar].
// Density and time storage remain owned by the common equation assembler.
struct IemSource {
  double explicit_source_density{}, implicit_sink_density{};
};
portable::Status iem_source(double volume_m3, double molecular_viscosity,
    double turbulent_viscosity, double cd, double control, double mean,
    IemSource&) noexcept;
// Frozen stochastic RHS per volume. All species and enthalpy in one field
// share one attenuation factor. Bounds are global for that field/component;
// Wiener increments are the already drawn, step-owned sqrt(s) values.
struct StochasticSourceRequest {
  double density{}, molecular_viscosity{}, turbulent_viscosity{};
  double molecular_schmidt{.7}, turbulent_schmidt{.5}, dt{};
  std::array<double,3> wiener{};
  const double *values{}, *gradients{}, *lower{}, *upper{};
  std::size_t components{};
  // Canonicalize unresolved increments and available bound distances before
  // common limiting. Accepted values and the caller repair budget stay intact.
  bool canonicalize_roundoff{};
};
struct StochasticSourceReport {
  portable::Status status{portable::Status::invalid_input};
  double attenuation{1};
};
// One attenuation for the full tuple. Optional canonicalization maps raw
// increments within 32 FP64 eps at max(1, |lower|, |upper|) to zero before
// bounding the complete tuple. Resolvable increments retain the same factor.
StochasticSourceReport stochastic_source(const StochasticSourceRequest&,
    double* source_density) noexcept;
struct Report {
  portable::Status status{portable::Status::invalid_input};
  std::size_t failure_field{}, failure_component{};
  View candidate{};
  double relaxation_factor{1}, max_mean_residual{}, max_element_residual{};
  const double *means{};
  const double *variances{};
  std::uint64_t model_identity{0x4553465048490001ULL};
  std::uint32_t chemistry_call_count{};
  std::array<double, 4> final_densities_kg_per_m3{};
  double ensemble_heat_release_j_per_m3{};
  // Equal-weight sum of validated chemical intervals, Ns entries.
  // Borrowed under candidate generation; null on a failed reaction.
  const double* mean_integrated_species_density_delta_kg_per_m3{};
  // Favre fields carry equal-weight mass fractions. The flow equation
  // multiplies this endpoint increment by its shared carrier density.
  const double* mean_integrated_mass_fraction_delta{};
};
enum class ReactionIntervals : std::uint8_t { full, two_halves };
struct ReactionRequest {
  View accepted{};
  portable::Revision expected_revision{};
  const combustion::ChemistryIdentity *chemistry_identity{};
  const portable::GasIdentity *gas_identity{};
  const double *pressures_pa{};
  const double *initial_densities_kg_per_m3{};
  double start_time_s{}, duration_s{};
  ReactionIntervals intervals{ReactionIntervals::full};
};
class Workspace {
public:
  explicit Workspace(std::size_t species_capacity);
  std::uint64_t owned_bytes() const noexcept {
    std::uint64_t bytes = sizeof(*this);
    for (const auto *v : {&candidate_, &transported_, &means_, &variances_,
                          &next_y_, &species_delta_, &mean_species_delta_, &mean_fraction_delta_})
      bytes += v->capacity() * sizeof(double);
    return bytes;
  }
  Report advance(const Request &) noexcept;
  // Translate to a prescribed physical mean, then apply the largest common
  // contraction that keeps every species in its mass-fraction simplex.
  // The same contraction acts on enthalpy fluctuations and all fields.
  Report recenter(const View&, const double* target_mean) noexcept;
  // Persistent field chemistry after transport/mixing. The default invokes
  // one full interval per field; two contiguous halves remain an explicit
  // integration comparison, with each provider increment validated.
  Report react(const ReactionRequest &,
               portable::GasAdvanceProvider &) noexcept;
  bool valid(const View &v) const noexcept {
    return v.owner == this && v.generation == generation_ &&
           v.values != nullptr;
  }

private:
  std::uint64_t generation_{};
  std::size_t capacity_{};
  std::vector<double> candidate_, transported_, means_, variances_;
  std::vector<double> next_y_, species_delta_, mean_species_delta_, mean_fraction_delta_;
};
// Shared IEM/TCR relaxation: exp[-cbrt(control) dt/(2 tau)].
portable::Status iem_factor(double dt, double tau, double control,
                            double &factor) noexcept;
} // namespace hundun::v04::esf::detail
