// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/stage3_scientific_row.hpp"

#include "hundun/rt_error.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace hundun::test::stage3 {
namespace {

constexpr std::string_view kPrefix = "STAGE3_SCIENTIFIC_ROW ";

bool positive(runtime::Int3 value) noexcept {
  return value.x > 0 && value.y > 0 && value.z > 0;
}

bool safe_row_id(std::string_view value) noexcept {
  return !value.empty() &&
         std::all_of(value.begin(), value.end(), [](char character) {
           return (character >= 'a' && character <= 'z') ||
                  (character >= '0' && character <= '9') || character == '-' ||
                  character == '_' || character == '.';
         });
}

bool metric_valid(const ScientificMetric &metric) noexcept {
  if (!metric.available)
    return metric.value == 0.0 && !std::signbit(metric.value);
  return std::isfinite(metric.value) && metric.value >= 0.0;
}

std::string dimensions(runtime::Int3 value) {
  return std::to_string(value.x) + 'x' + std::to_string(value.y) + 'x' +
         std::to_string(value.z);
}

runtime::Int3 parse_dimensions(const std::string &value) {
  const auto first = value.find('x');
  const auto second = first == std::string::npos ? std::string::npos
                                                 : value.find('x', first + 1U);
  if (first == std::string::npos || second == std::string::npos ||
      value.find('x', second + 1U) != std::string::npos)
    throw runtime::Error("Stage 3 scientific row dimensions are invalid");
  try {
    std::size_t used_x{}, used_y{}, used_z{};
    const int x = std::stoi(value.substr(0U, first), &used_x);
    const int y =
        std::stoi(value.substr(first + 1U, second - first - 1U), &used_y);
    const int z = std::stoi(value.substr(second + 1U), &used_z);
    if (used_x != first || used_y != second - first - 1U ||
        used_z != value.size() - second - 1U)
      throw runtime::Error("Stage 3 scientific row dimensions are invalid");
    return {x, y, z};
  } catch (const runtime::Error &) {
    throw;
  } catch (...) {
    throw runtime::Error("Stage 3 scientific row dimensions are invalid");
  }
}

std::uint64_t parse_u64(const std::string &value, const char *message) {
  try {
    std::size_t used{};
    const auto parsed = std::stoull(value, &used);
    if (used != value.size())
      throw runtime::Error(message);
    return parsed;
  } catch (const runtime::Error &) {
    throw;
  } catch (...) {
    throw runtime::Error(message);
  }
}

int parse_int(const std::string &value, const char *message) {
  const auto parsed = parse_u64(value, message);
  if (parsed > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw runtime::Error(message);
  return static_cast<int>(parsed);
}

double parse_double(const std::string &value, const char *message) {
  try {
    std::size_t used{};
    const double parsed = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(parsed))
      throw runtime::Error(message);
    return parsed;
  } catch (const runtime::Error &) {
    throw;
  } catch (...) {
    throw runtime::Error(message);
  }
}

bool parse_available(const std::map<std::string, std::string> &values,
                     const std::string &name) {
  const auto found = values.find(name + "_available");
  if (found == values.end() || (found->second != "0" && found->second != "1"))
    throw runtime::Error("Stage 3 scientific availability bit is invalid");
  return found->second == "1";
}

ScientificMetric parse_metric(const std::map<std::string, std::string> &values,
                              const std::string &name) {
  const bool available = parse_available(values, name);
  const auto found = values.find(name);
  if (available != (found != values.end()))
    throw runtime::Error("Stage 3 scientific metric availability is invalid");
  return {available, available
                         ? parse_double(found->second,
                                        "Stage 3 scientific metric is invalid")
                         : 0.0};
}

std::size_t cell_count(runtime::Int3 cells) {
  if (!positive(cells))
    throw runtime::Error("Stage 3 scientific snapshot extent is invalid");
  const auto x = static_cast<std::size_t>(cells.x);
  const auto y = static_cast<std::size_t>(cells.y);
  const auto z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y ||
      x * y > std::numeric_limits<std::size_t>::max() / z)
    throw runtime::Error("Stage 3 scientific snapshot extent overflows");
  return x * y * z;
}

std::size_t offset(runtime::Int3 cells, int i, int j, int k) noexcept {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(cells.y) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(cells.x) +
         static_cast<std::size_t>(i);
}

std::array<double, 2> observed_orders(const std::array<double, 3> &errors,
                                      bool &valid) noexcept {
  constexpr double minimum_error =
      1024.0 * std::numeric_limits<double>::epsilon();
  for (double value : errors)
    valid = valid && std::isfinite(value) && value > minimum_error;
  valid = valid && errors[0] > errors[1] && errors[1] > errors[2];
  if (!valid)
    return {};
  const std::array<double, 2> result{
      std::log(errors[0] / errors[1]) / std::log(2.0),
      std::log(errors[1] / errors[2]) / std::log(2.0)};
  valid = std::all_of(result.begin(), result.end(), [](double value) {
    return std::isfinite(value) && value >= 1.8;
  });
  return result;
}

} // namespace

bool validate_scientific_row(const ScientificRow &row) noexcept {
  if (!safe_row_id(row.row_id) || !positive(row.cells) || row.ranks <= 0 ||
      !positive(row.process_grid) || row.steps == 0U || !(row.dt_s > 0.0) ||
      !std::isfinite(row.dt_s) || !(row.final_time_s > 0.0) ||
      !std::isfinite(row.final_time_s) ||
      (row.status != "pass" && row.status != "fail"))
    return false;
  const std::int64_t grid = static_cast<std::int64_t>(row.process_grid.x) *
                            row.process_grid.y * row.process_grid.z;
  if (grid != row.ranks)
    return false;
  const double expected_time = static_cast<double>(row.steps) * row.dt_s;
  const double time_scale =
      std::max({std::abs(expected_time), std::abs(row.final_time_s), 1.0});
  if (std::abs(expected_time - row.final_time_s) >
      16.0 * std::numeric_limits<double>::epsilon() * time_scale)
    return false;
  if (!metric_valid(row.velocity_l2) || !metric_valid(row.pressure_l2) ||
      !metric_valid(row.nu_t_l2) || !metric_valid(row.continuity) ||
      !metric_valid(row.conservation) || !metric_valid(row.closure) ||
      !metric_valid(row.force))
    return false;
  return row.wale_identity.available ? row.wale_identity.value != 0U
                                     : row.wale_identity.value == 0U;
}

std::string serialize_scientific_row(const ScientificRow &row) {
  if (!validate_scientific_row(row))
    throw runtime::Error("Stage 3 scientific row is invalid");
  std::ostringstream output;
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << kPrefix << "row_id=" << row.row_id
         << " cells=" << dimensions(row.cells) << " ranks=" << row.ranks
         << " process_grid=" << dimensions(row.process_grid)
         << " steps=" << row.steps << " dt=" << row.dt_s
         << " final_time=" << row.final_time_s;
  const auto append_metric = [&](const char *name,
                                 const ScientificMetric &metric) {
    output << ' ' << name << "_available=" << (metric.available ? 1 : 0);
    if (metric.available)
      output << ' ' << name << '=' << metric.value;
  };
  append_metric("velocity_l2", row.velocity_l2);
  append_metric("pressure_l2", row.pressure_l2);
  append_metric("nu_t_l2", row.nu_t_l2);
  append_metric("continuity", row.continuity);
  append_metric("conservation", row.conservation);
  append_metric("closure", row.closure);
  append_metric("force", row.force);
  output << " wale_identity_available="
         << (row.wale_identity.available ? 1 : 0);
  if (row.wale_identity.available)
    output << " wale_identity=" << row.wale_identity.value;
  output << " status=" << row.status;
  return output.str();
}

ScientificRow parse_scientific_row(std::string_view text) {
  if (text.substr(0U, kPrefix.size()) != kPrefix)
    throw runtime::Error("Stage 3 scientific row prefix is invalid");
  std::istringstream input(std::string{text.substr(kPrefix.size())});
  std::map<std::string, std::string> values;
  std::string token;
  while (input >> token) {
    const auto split = token.find('=');
    if (split == std::string::npos || split == 0U ||
        split + 1U == token.size() ||
        !values.emplace(token.substr(0U, split), token.substr(split + 1U))
             .second)
      throw runtime::Error("Stage 3 scientific row token is invalid");
  }
  const auto required = [&](const char *name) -> const std::string & {
    const auto found = values.find(name);
    if (found == values.end())
      throw runtime::Error("Stage 3 scientific row field is missing");
    return found->second;
  };
  ScientificRow row;
  row.row_id = required("row_id");
  row.cells = parse_dimensions(required("cells"));
  row.ranks = parse_int(required("ranks"), "Stage 3 ranks are invalid");
  row.process_grid = parse_dimensions(required("process_grid"));
  row.steps = parse_u64(required("steps"), "Stage 3 steps are invalid");
  row.dt_s = parse_double(required("dt"), "Stage 3 dt is invalid");
  row.final_time_s =
      parse_double(required("final_time"), "Stage 3 final time is invalid");
  row.velocity_l2 = parse_metric(values, "velocity_l2");
  row.pressure_l2 = parse_metric(values, "pressure_l2");
  row.nu_t_l2 = parse_metric(values, "nu_t_l2");
  row.continuity = parse_metric(values, "continuity");
  row.conservation = parse_metric(values, "conservation");
  row.closure = parse_metric(values, "closure");
  row.force = parse_metric(values, "force");
  row.wale_identity.available = parse_available(values, "wale_identity");
  const auto wale = values.find("wale_identity");
  if (row.wale_identity.available != (wale != values.end()))
    throw runtime::Error("Stage 3 WALE identity availability is invalid");
  row.wale_identity.value =
      row.wale_identity.available
          ? parse_u64(wale->second, "Stage 3 WALE identity is invalid")
          : 0U;
  row.status = required("status");
  const std::size_t expected_fields =
      8U + 7U + static_cast<std::size_t>(row.velocity_l2.available) +
      static_cast<std::size_t>(row.pressure_l2.available) +
      static_cast<std::size_t>(row.nu_t_l2.available) +
      static_cast<std::size_t>(row.continuity.available) +
      static_cast<std::size_t>(row.conservation.available) +
      static_cast<std::size_t>(row.closure.available) +
      static_cast<std::size_t>(row.force.available) + 1U +
      static_cast<std::size_t>(row.wale_identity.available);
  if (values.size() != expected_fields || !validate_scientific_row(row))
    throw runtime::Error("Stage 3 scientific row is invalid");
  return row;
}

CellAverageSnapshot restrict_cell_averages(const CellAverageSnapshot &fine,
                                           runtime::Int3 coarse,
                                           CellRestriction restriction) {
  const auto fine_count = cell_count(fine.cells);
  const auto coarse_count = cell_count(coarse);
  const bool has_nu_t = fine.nu_t.size() == fine_count;
  if (fine.velocity.size() != fine_count * 3U ||
      fine.pressure.size() != fine_count || (!fine.nu_t.empty() && !has_nu_t) ||
      fine.cells.x % coarse.x != 0 || fine.cells.y % coarse.y != 0 ||
      fine.cells.z % coarse.z != 0)
    throw runtime::Error("Stage 3 cell-average restriction layout is invalid");
  const runtime::Int3 ratio{fine.cells.x / coarse.x, fine.cells.y / coarse.y,
                            fine.cells.z / coarse.z};
  if (!positive(ratio))
    throw runtime::Error("Stage 3 cell-average restriction ratio is invalid");
  CellAverageSnapshot result;
  result.cells = coarse;
  result.velocity.assign(coarse_count * 3U, 0.0);
  result.pressure.assign(coarse_count, 0.0);
  if (has_nu_t)
    result.nu_t.assign(coarse_count, 0.0);
  const double divisor = restriction == CellRestriction::cell_average
                             ? static_cast<double>(ratio.x * ratio.y * ratio.z)
                             : 1.0;
  for (int k = 0; k < coarse.z; ++k) {
    for (int j = 0; j < coarse.y; ++j) {
      for (int i = 0; i < coarse.x; ++i) {
        const auto coarse_offset = offset(coarse, i, j, k);
        const int dk_end =
            restriction == CellRestriction::cell_average ? ratio.z : 1;
        const int dj_end =
            restriction == CellRestriction::cell_average ? ratio.y : 1;
        const int di_end =
            restriction == CellRestriction::cell_average ? ratio.x : 1;
        for (int dk = 0; dk < dk_end; ++dk) {
          for (int dj = 0; dj < dj_end; ++dj) {
            for (int di = 0; di < di_end; ++di) {
              const auto fine_offset =
                  offset(fine.cells, i * ratio.x + di, j * ratio.y + dj,
                         k * ratio.z + dk);
              result.pressure[coarse_offset] += fine.pressure[fine_offset];
              if (has_nu_t)
                result.nu_t[coarse_offset] += fine.nu_t[fine_offset];
              for (std::size_t component = 0U; component < 3U; ++component)
                result.velocity[coarse_offset * 3U + component] +=
                    fine.velocity[fine_offset * 3U + component];
            }
          }
        }
        result.pressure[coarse_offset] /= divisor;
        if (has_nu_t)
          result.nu_t[coarse_offset] /= divisor;
        for (std::size_t component = 0U; component < 3U; ++component)
          result.velocity[coarse_offset * 3U + component] /= divisor;
      }
    }
  }
  return result;
}

TgvConvergenceReport assess_tgv_convergence(const TgvConvergenceInput &input) {
  TgvConvergenceReport result;
  bool valid = input.restriction == CellRestriction::cell_average;
  const double time_scale = std::max({std::abs(input.final_time_s[0]),
                                      std::abs(input.final_time_s[1]),
                                      std::abs(input.final_time_s[2]), 1.0});
  for (double time : input.final_time_s)
    valid = valid && std::isfinite(time) && time > 0.0 &&
            std::abs(time - input.final_time_s[0]) <=
                16.0 * std::numeric_limits<double>::epsilon() * time_scale;
  for (std::size_t level = 0U; level < 3U; ++level) {
    valid =
        valid && input.density_revision[level] != 0U &&
        input.density_revision[level] == input.wale_density_revision[level] &&
        input.wale_evaluation_count[level] == 1U;
  }
  result.velocity_order = observed_orders(input.velocity_l2, valid);
  result.pressure_order = observed_orders(input.pressure_l2, valid);
  result.nu_t_order = observed_orders(input.nu_t_l2, valid);
  result.accepted = valid;
  return result;
}

} // namespace hundun::test::stage3
