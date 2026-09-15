// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_portable.hpp"
#include "models_spray_transfer_detail.hpp"
#include <filesystem>
#include <string>
#include <vector>
namespace hundun::v04::spray::detail {
struct LiquidAsset {
  LiquidPropertyPack pack;
  double reference_temperature_k{}, reference_liquid_enthalpy_j_per_kg{};
  double vapor_molecular_weight_kg_per_kmol{};
  std::size_t vapor_species_index{};
  std::string vapor_species_name, source;
  portable::GasIdentity gas_identity;
  // FNV1a64 of exact file bytes, a versioned content identity, not SHA256.
  std::uint64_t content_fingerprint{};
};
struct LiquidAssetReport {
  portable::Status status{portable::Status::invalid_input};
  LiquidAsset asset;
  bool available{};
};
// Cold path. Strict SI format, bounded file length; exceptions become status.
LiquidAssetReport load_liquid_asset(const std::filesystem::path &,
                                    std::uint64_t expected_fnv1a64,
                                    const portable::GasIdentity &) noexcept;
struct LiquidEnthalpyReport {
  portable::Status status{portable::Status::invalid_input};
  double liquid_enthalpy_j_per_kg{}, cp_j_per_kg_k{};
  bool available{};
};
LiquidEnthalpyReport evaluate_liquid_enthalpy(const LiquidAsset &,
                                              double temperature_k) noexcept;
// Subcritical kerosene surrogate phase relations used by THICK_EX. The
// sensible increment integrates liquid cp; its absolute reference belongs
// to the material asset and the gas-liquid exchange ledger.
// Callers admit temperatures through their material/mechanism range first.
struct KerosenePhaseReport {
  portable::Status status{portable::Status::invalid_input};
  double saturation_pressure_pa{}, boiling_temperature_k{};
  double liquid_density_kg_per_m3{};
  double liquid_cp_j_per_kg_k{}, latent_heat_j_per_kg{};
  double sensible_enthalpy_increment_j_per_kg{};
  bool available{};
};
KerosenePhaseReport evaluate_kerosene_phase(double temperature_k,
    double pressure_pa, double reference_temperature_k) noexcept;
struct FilmTransportReport {
  portable::Status status{portable::Status::unavailable};
  portable::Revision revision{};
  double dynamic_viscosity_pa_s{};
};
// Queries the admitted molecular transport law at the sampled far PT/Y.
class FilmTransportProvider {
 public:
  virtual ~FilmTransportProvider() = default;
  virtual FilmTransportReport query(const portable::GasQuery&) const noexcept = 0;
};
struct KeroseneFilmInput {
  portable::GasQuery far_gas;
  portable::Revision expected_revision{}, transport_revision{};
  double surface_temperature_k{}, far_dynamic_viscosity_pa_s{};
  const FilmTransportProvider* transport{};
};
struct KeroseneFilmReport {
  portable::Status status{portable::Status::invalid_input};
  portable::Revision revision{};
  KerosenePhaseReport liquid;
  double far_temperature_k{}, film_temperature_k{};
  double gas_cp_j_per_kg_k{}, gas_dynamic_viscosity_pa_s{};
  double vapor_cp_j_per_kg_k{}, vapor_prandtl_number{};
  double vapor_absolute_enthalpy_j_per_kg{}, surface_vapor_mass_fraction{};
  double far_density_kg_per_m3{}, far_dynamic_viscosity_pa_s{};
  bool available{};
};
// THICK_EX evaluates cp at film T with far composition and Pr from vapor.
// The admitted flow material supplies viscosity with the same revision.
class KeroseneFilmWorkspace {
 public:
  explicit KeroseneFilmWorkspace(std::size_t species);
  KeroseneFilmReport query(const LiquidAsset&, portable::GasQueryProvider&,
                          const KeroseneFilmInput&) noexcept;
  std::size_t owned_payload_bytes() const noexcept;
 private:
  std::vector<double> pure_y_, d_, h_, w_;
};
struct FilmQueryInput {
  portable::GasQuery far_gas;
  portable::Revision expected_revision;
  Vector3 velocity_m_per_s{};
  double surface_temperature_k{};
};
struct FilmQueryReport {
  portable::Status status{portable::Status::invalid_input};
  OneThirdFilmSample film;
  double liquid_enthalpy_j_per_kg{}, vapor_enthalpy_j_per_kg{};
  double latent_heat_consistency_residual_j_per_kg{};
  double far_gas_density_kg_per_m3{};
  bool available{};
};
class FilmQueryWorkspace {
public:
  explicit FilmQueryWorkspace(std::size_t species);
  FilmQueryReport query(const LiquidAsset &, portable::GasQueryProvider &,
                        const FilmQueryInput &) noexcept;
  // Full Ns film composition; unavailable after every failed query.
  const double *film_mass_fractions() const noexcept;

private:
  bool available_{};
  std::vector<double> film_y_, surface_y_, pure_y_, d_, h_, w_;
};
} // namespace hundun::v04::spray::detail
