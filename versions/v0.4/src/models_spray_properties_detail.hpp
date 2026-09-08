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
