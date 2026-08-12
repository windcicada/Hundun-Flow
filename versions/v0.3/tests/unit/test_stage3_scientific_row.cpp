// SPDX-License-Identifier: Apache-2.0

#include "tests/support/stage3_scientific_row.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/rt_error.hpp"

#include <cmath>
#include <limits>
#include <string>

namespace {

using namespace hundun;

test::stage3::ScientificRow row_fixture() {
  return {"wale-tgv-n12-r1",
          {12, 12, 12},
          1,
          {1, 1, 1},
          1U,
          1.0e-4,
          1.0e-4,
          {true, 1.0e-3},
          {true, 2.0e-3},
          {true, 3.0e-3},
          {true, 4.0e-11},
          {true, 5.0e-12},
          {false, 0.0},
          {false, 0.0},
          {true, 42U},
          "pass"};
}

void run() {
  using namespace test::stage3;
  const auto row = row_fixture();
  HUNDUN_CHECK(validate_scientific_row(row));
  const auto encoded = serialize_scientific_row(row);
  HUNDUN_CHECK(encoded.rfind("STAGE3_SCIENTIFIC_ROW ", 0U) == 0U);
  HUNDUN_CHECK(serialize_scientific_row(parse_scientific_row(encoded)) ==
               encoded);
  HUNDUN_CHECK(encoded.find("closure_available=0") != std::string::npos);
  HUNDUN_CHECK(encoded.find(" closure=") == std::string::npos);

  auto final_time = row;
  final_time.final_time_s = 2.0e-4;
  HUNDUN_CHECK(!validate_scientific_row(final_time));
  auto wrong_grid = row;
  wrong_grid.process_grid = {2, 1, 1};
  HUNDUN_CHECK(!validate_scientific_row(wrong_grid));
  auto wrong_ranks = row;
  wrong_ranks.ranks = 2;
  HUNDUN_CHECK(!validate_scientific_row(wrong_ranks));
  auto zero_identity = row;
  zero_identity.wale_identity.value = 0U;
  HUNDUN_CHECK(!validate_scientific_row(zero_identity));

  auto missing_availability = encoded;
  const auto marker = missing_availability.find("pressure_l2_available=1 ");
  HUNDUN_CHECK(marker != std::string::npos);
  missing_availability.erase(marker,
                             std::string("pressure_l2_available=1 ").size());
  bool missing_rejected = false;
  try {
    static_cast<void>(parse_scientific_row(missing_availability));
  } catch (const runtime::Error &) {
    missing_rejected = true;
  }
  HUNDUN_CHECK(missing_rejected);

  CellAverageSnapshot fine;
  fine.cells = {4, 4, 4};
  fine.velocity.resize(4U * 4U * 4U * 3U);
  fine.pressure.resize(4U * 4U * 4U);
  fine.nu_t.resize(4U * 4U * 4U);
  for (std::size_t cell = 0U; cell < fine.pressure.size(); ++cell) {
    fine.pressure[cell] = static_cast<double>(cell + 1U);
    fine.nu_t[cell] = 0.5 * static_cast<double>(cell + 1U);
    for (std::size_t component = 0U; component < 3U; ++component)
      fine.velocity[cell * 3U + component] =
          static_cast<double>((cell + 1U) * (component + 1U));
  }
  const auto averaged =
      restrict_cell_averages(fine, {2, 2, 2}, CellRestriction::cell_average);
  const auto sampled =
      restrict_cell_averages(fine, {2, 2, 2}, CellRestriction::point_sample);
  HUNDUN_CHECK(averaged.cells.x == 2 && averaged.pressure.size() == 8U);
  HUNDUN_CHECK(averaged.pressure != sampled.pressure);
  auto without_nu_t = fine;
  without_nu_t.nu_t.clear();
  HUNDUN_CHECK(restrict_cell_averages(without_nu_t, {2, 2, 2},
                                      CellRestriction::cell_average)
                   .nu_t.empty());
  auto partial_nu_t = fine;
  partial_nu_t.nu_t.pop_back();
  bool partial_nu_t_rejected = false;
  try {
    static_cast<void>(restrict_cell_averages(partial_nu_t, {2, 2, 2},
                                             CellRestriction::cell_average));
  } catch (const runtime::Error &) {
    partial_nu_t_rejected = true;
  }
  HUNDUN_CHECK(partial_nu_t_rejected);

  TgvConvergenceInput convergence;
  convergence.velocity_l2 = {4.0e-3, 1.0e-3, 2.5e-4};
  convergence.pressure_l2 = {8.0e-3, 2.0e-3, 5.0e-4};
  convergence.nu_t_l2 = {2.0e-3, 5.0e-4, 1.25e-4};
  convergence.final_time_s = {1.0e-4, 1.0e-4, 1.0e-4};
  convergence.density_revision = {1U, 2U, 3U};
  convergence.wale_density_revision = {1U, 2U, 3U};
  convergence.wale_evaluation_count = {1U, 1U, 1U};
  const auto accepted = assess_tgv_convergence(convergence);
  HUNDUN_CHECK(accepted.accepted);
  for (const auto order : accepted.velocity_order)
    HUNDUN_CHECK(order >= 1.8);

  const auto rejected = [](TgvConvergenceInput changed) {
    HUNDUN_CHECK(!assess_tgv_convergence(changed).accepted);
  };
  auto mismatch = convergence;
  mismatch.final_time_s[2] = 2.0e-4;
  rejected(mismatch);
  auto point_sample = convergence;
  point_sample.restriction = CellRestriction::point_sample;
  rejected(point_sample);
  auto hidden_segment = convergence;
  hidden_segment.velocity_l2 = {4.0e-3, 3.9e-3, 2.5e-4};
  rejected(hidden_segment);
  auto zero_error = convergence;
  zero_error.pressure_l2[2] = 0.0;
  rejected(zero_error);
  auto epsilon_error = convergence;
  epsilon_error.nu_t_l2[2] = std::numeric_limits<double>::epsilon();
  rejected(epsilon_error);
  auto stale_density = convergence;
  stale_density.wale_density_revision[1] = 1U;
  rejected(stale_density);
  auto second_evaluation = convergence;
  second_evaluation.wale_evaluation_count[2] = 2U;
  rejected(second_evaluation);
}

} // namespace

int main() { return hundun::test::run(run); }
