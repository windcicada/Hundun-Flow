// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_collective_status.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::flow::detail {

struct ReactingLayerState final {
  std::vector<double> rho_y_kg_per_m3;
  std::vector<double> rho_h_tc_j_per_m3;
  double p0_pa{};
};

bool operator==(const ReactingLayerState &left,
                const ReactingLayerState &right) noexcept;
bool operator!=(const ReactingLayerState &left,
                const ReactingLayerState &right) noexcept;

class ReactingAttemptState final {
public:
  ReactingAttemptState(std::size_t cell_count, std::size_t species_count,
                       ReactingLayerState history,
                       ReactingLayerState committed);

  std::size_t cell_count() const noexcept;
  std::size_t species_count() const noexcept;
  const ReactingLayerState &history() const noexcept;
  const ReactingLayerState &committed() const noexcept;
  bool attempt_active() const noexcept;
  void begin_attempt();
  void rollback_attempt();

private:
  friend class ReactingSourceTransaction;
  void publish_attempt();
  std::size_t cell_count_{};
  std::size_t species_count_{};
  ReactingLayerState history_;
  ReactingLayerState committed_;
  ReactingLayerState trial_;
  bool attempt_active_{};
};

enum class ReactingSourceKind : std::uint32_t {
  chemistry = 0U,
  transport = 1U,
  boundary = 2U
};

struct ReactingSourceIdentity final {
  std::string id;
  ReactingSourceKind kind{ReactingSourceKind::transport};
};

enum class ReactingSourceQuantity : std::uint32_t {
  species_density = 0U,
  enthalpy_density = 1U
};

struct ReactingSourceRecord final {
  ReactingSourceIdentity source;
  ReactingSourceQuantity quantity{ReactingSourceQuantity::species_density};
  std::size_t cell{};
  std::size_t component{};
  double delta{};
  std::string_view units;
};

class ReactingSourceTransaction final {
public:
  explicit ReactingSourceTransaction(ReactingAttemptState &state);
  ReactingSourceTransaction(const ReactingSourceTransaction &) = delete;
  ReactingSourceTransaction &
  operator=(const ReactingSourceTransaction &) = delete;

  void add_species(const ReactingSourceIdentity &source, std::size_t cell,
                   std::size_t species, double delta_kg_per_m3);
  void add_enthalpy(const ReactingSourceIdentity &source, std::size_t cell,
                    double delta_j_per_m3);
  const std::vector<double> &species_delta_kg_per_m3() const noexcept;
  const std::vector<double> &enthalpy_delta_j_per_m3() const noexcept;
  const std::vector<ReactingSourceRecord> &records() const noexcept;
  bool commit(const runtime::CollectiveStatus &status);
  void rollback();

private:
  void require_open() const;
  void rollback_noexcept() noexcept;
  ReactingAttemptState *state_{};
  std::vector<double> species_delta_;
  std::vector<double> enthalpy_delta_;
  std::vector<double> chemistry_mass_delta_;
  std::vector<ReactingSourceRecord> records_;
  bool closed_{};
};

} // namespace hundun::flow::detail
