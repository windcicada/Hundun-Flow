// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_types.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::test::stage3 {

struct ScientificMetric final {
  bool available{};
  double value{};
};

struct ScientificIdentity final {
  bool available{};
  std::uint64_t value{};
};

struct ScientificRow final {
  std::string row_id;
  runtime::Int3 cells{};
  int ranks{};
  runtime::Int3 process_grid{};
  std::uint64_t steps{};
  double dt_s{};
  double final_time_s{};
  ScientificMetric velocity_l2;
  ScientificMetric pressure_l2;
  ScientificMetric nu_t_l2;
  ScientificMetric continuity;
  ScientificMetric conservation;
  ScientificMetric closure;
  ScientificMetric force;
  ScientificIdentity wale_identity;
  std::string status;
};

bool validate_scientific_row(const ScientificRow &) noexcept;
std::string serialize_scientific_row(const ScientificRow &);
ScientificRow parse_scientific_row(std::string_view);

enum class CellRestriction : std::uint8_t {
  cell_average,
  point_sample,
};

struct CellAverageSnapshot final {
  runtime::Int3 cells{};
  std::vector<double> velocity;
  std::vector<double> pressure;
  std::vector<double> nu_t;
};

CellAverageSnapshot restrict_cell_averages(const CellAverageSnapshot &,
                                           runtime::Int3 coarse_cells,
                                           CellRestriction);

struct TgvConvergenceInput final {
  std::array<double, 3> velocity_l2{};
  std::array<double, 3> pressure_l2{};
  std::array<double, 3> nu_t_l2{};
  std::array<double, 3> final_time_s{};
  std::array<std::uint64_t, 3> density_revision{};
  std::array<std::uint64_t, 3> wale_density_revision{};
  std::array<std::uint64_t, 3> wale_evaluation_count{};
  CellRestriction restriction{CellRestriction::cell_average};
};

struct TgvConvergenceReport final {
  std::array<double, 2> velocity_order{};
  std::array<double, 2> pressure_order{};
  std::array<double, 2> nu_t_order{};
  bool accepted{};
};

TgvConvergenceReport assess_tgv_convergence(const TgvConvergenceInput &);

} // namespace hundun::test::stage3
