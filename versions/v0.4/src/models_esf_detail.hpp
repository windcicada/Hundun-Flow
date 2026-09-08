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
  // Equal-weight sum of both validated chemical half-intervals, Ns entries.
  // Borrowed under candidate generation; null on a failed reaction.
  const double* mean_integrated_species_density_delta_kg_per_m3{};
};
struct ReactionRequest {
  View accepted{};
  portable::Revision expected_revision{};
  const combustion::ChemistryIdentity *chemistry_identity{};
  const portable::GasIdentity *gas_identity{};
  const double *pressures_pa{};
  const double *initial_densities_kg_per_m3{};
  double start_time_s{}, duration_s{};
};
class Workspace {
public:
  explicit Workspace(std::size_t species_capacity);
  Report advance(const Request &) noexcept;
  // Persistent field chemistry, exactly two half-intervals per field.
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
  std::vector<double> next_y_, species_delta_, mean_species_delta_;
};
// Shared IEM/TCR relaxation: exp[-cbrt(control) dt/(2 tau)].
portable::Status iem_factor(double dt, double tau, double control,
                            double &factor) noexcept;
} // namespace hundun::v04::esf::detail
