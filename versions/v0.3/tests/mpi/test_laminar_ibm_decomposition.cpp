// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/stage3_decomposition_equality.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_test_contracts.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

namespace {

using namespace hundun;

int parse_positive(const char *text, const char *name) {
  if (text == nullptr || *text == '\0')
    throw runtime::Error(std::string("missing ") + name);
  char *end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (*end != '\0' || value <= 0 || value > std::numeric_limits<int>::max())
    throw runtime::Error(std::string("invalid ") + name);
  return static_cast<int>(value);
}

runtime::Int3 default_process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported Task 11 decomposition rank count");
}

test::stage3::ManufacturedCase make_case(int cells,
                                         runtime::Int3 process_grid) {
  test::stage3::ManufacturedCase result;
  result.body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  result.cells = cells;
  result.mapping = test::stage3::ManufacturedMapping::uniform;
  result.fluid_side = config::ImmersedFluidSide::outside;
  result.process_grid = process_grid;
  result.collect_force = true;
  return result;
}

void run(int argc, char **argv, const runtime::MpiContext &world) {
  HUNDUN_CHECK(
      test::stage3::decomposition_equality_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(argc == 2 || argc == 5);
  const std::string selector = argv[1];
  HUNDUN_CHECK(selector == "fast" || selector == "acceptance");
  const int cells = selector == "acceptance" ? 24 : 16;
  runtime::Int3 process_grid = default_process_grid(world.size());
  if (argc == 5)
    process_grid = {parse_positive(argv[2], "process-grid x"),
                    parse_positive(argv[3], "process-grid y"),
                    parse_positive(argv[4], "process-grid z")};
  HUNDUN_CHECK(process_grid.x * process_grid.y * process_grid.z ==
               world.size());

  const auto distributed = test::stage3::run_manufactured_case(
      world, make_case(cells, process_grid));
  std::optional<test::stage3::ManufacturedRunResult> reference;
  int reference_ok = 1;
  if (world.rank() == 0) {
    try {
      const auto self = runtime::MpiContext::duplicate(MPI_COMM_SELF);
      reference.emplace(test::stage3::run_manufactured_case(
          self, make_case(cells, {1, 1, 1})));
    } catch (...) {
      reference_ok = 0;
    }
  }
  HUNDUN_CHECK(MPI_Bcast(&reference_ok, 1, MPI_INT, 0, world.comm()) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(reference_ok == 1);
  int comparison_ok = 1;
  if (world.rank() == 0) {
    try {
      HUNDUN_CHECK(reference.has_value());
      HUNDUN_CHECK(reference->global_active_cell_ids ==
                   distributed.global_active_cell_ids);
      HUNDUN_CHECK(reference->active_cell_count ==
                   distributed.active_cell_count);
      HUNDUN_CHECK(reference->immersed_link_count ==
                   distributed.immersed_link_count);
      HUNDUN_CHECK(reference->wall_point_count == distributed.wall_point_count);
      HUNDUN_CHECK(reference->classification_fingerprint ==
                   distributed.classification_fingerprint);
      HUNDUN_CHECK(reference->surface_coverage_fingerprint ==
                   distributed.surface_coverage_fingerprint);
      HUNDUN_CHECK(reference->ghost_plan_fingerprint ==
                   distributed.ghost_plan_fingerprint);
      HUNDUN_CHECK(reference->wall_plan_fingerprint ==
                   distributed.wall_plan_fingerprint);
      HUNDUN_CHECK(reference->operator_structure_fingerprint ==
                   distributed.operator_structure_fingerprint);
      HUNDUN_CHECK(reference->pressure_error_extremum.global_cell_id ==
                   distributed.pressure_error_extremum.global_cell_id);
      HUNDUN_CHECK(reference->pressure_error_extremum.logical_cell.x ==
                   distributed.pressure_error_extremum.logical_cell.x);
      HUNDUN_CHECK(reference->pressure_error_extremum.logical_cell.y ==
                   distributed.pressure_error_extremum.logical_cell.y);
      HUNDUN_CHECK(reference->pressure_error_extremum.logical_cell.z ==
                   distributed.pressure_error_extremum.logical_cell.z);
      test::stage3::require_decomposition_field(
          "pressure error extremum",
          {reference->pressure_error_extremum.signed_error_pa,
           reference->pressure_error_extremum.absolute_error_pa,
           reference->pressure_error_extremum.wall_distance_m},
          {distributed.pressure_error_extremum.signed_error_pa,
           distributed.pressure_error_extremum.absolute_error_pa,
           distributed.pressure_error_extremum.wall_distance_m});
      test::stage3::require_decomposition_field(
          "pressure", reference->global_pressure, distributed.global_pressure);
      test::stage3::require_decomposition_field(
          "velocity", reference->global_velocity, distributed.global_velocity);
      test::stage3::require_decomposition_force("operator force",
                                                reference->operator_force,
                                                distributed.operator_force);
      test::stage3::require_decomposition_force("budget reaction",
                                                reference->budget_reaction,
                                                distributed.budget_reaction);
      test::stage3::require_decomposition_force("surface traction",
                                                reference->surface_traction,
                                                distributed.surface_traction);
      test::stage3::require_decomposition_force("operator/surface consistency",
                                                reference->consistency,
                                                distributed.consistency);
    } catch (const std::exception &error) {
      std::cerr << "Task 11 decomposition comparison failed: " << error.what()
                << '\n';
      comparison_ok = 0;
    }
  }
  HUNDUN_CHECK(MPI_Bcast(&comparison_ok, 1, MPI_INT, 0, world.comm()) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(comparison_ok == 1);

  auto failure_case = make_case(cells, process_grid);
  failure_case.final_force_failure_rank = world.size() - 1;
  const auto failure = test::stage3::run_manufactured_case(world, failure_case);
  HUNDUN_CHECK(!failure.committed);
  HUNDUN_CHECK(failure.rollback_bitwise_equal);
  HUNDUN_CHECK(failure.lowest_failing_rank == world.size() - 1);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  const auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run(argc, argv, mpi); });
}
