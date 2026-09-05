// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>
#include <unistd.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"

namespace {
using namespace hundun::v04;
bool same_fields(Span<const RestartFieldView> actual,
                 const std::vector<RestartImageField>& expected) {
  if (actual.size != expected.size()) return false;
  for (std::size_t f = 0; f < actual.size; ++f) {
    const auto& view = actual.data[f].values;
    const auto& field = expected[f];
    if (actual.data[f].role != field.role || view.field != field.field ||
        view.components != field.components)
      return false;
    std::size_t i = 0;
    for (int z = 0; z < view.interior.z; ++z)
      for (int y = 0; y < view.interior.y; ++y)
        for (int x = 0; x < view.interior.x; ++x)
          for (std::uint8_t c = 0; c < view.components; ++c)
            if (view.unchecked({x, y, z}, c) != field.values[i++]) return false;
    if (i != field.values.size()) return false;
  }
  return true;
}
bool same_flux(ConstFaceFluxView actual,
               const std::array<std::vector<double>, 3>& expected) {
  const std::array<ConstFaceFieldView, 3> views{actual.x, actual.y, actual.z};
  for (std::size_t axis = 0; axis < views.size(); ++axis) {
    const auto& v = views[axis];
    std::size_t i = 0;
    for (int z = 0; z < v.extents.z; ++z)
      for (int y = 0; y < v.extents.y; ++y)
        for (int x = 0; x < v.extents.x; ++x)
          if (v.unchecked({x, y, z}) != expected[axis][i++]) return false;
    if (i != expected[axis].size()) return false;
  }
  return true;
}
}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0, id = static_cast<int>(getpid());
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Bcast(&id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  const auto root = std::filesystem::temp_directory_path() /
                    ("hundun-restart-storage-test-" + std::to_string(id));
  bool passed = true;
  const auto agree = [&] {
    int local = passed ? 1 : 0, global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    passed = global != 0;
  };
  {
    auto model = test::product_model({8, 8, 8});
    const auto create = [&](ProductDriver& driver) {
      CompiledCasePlan plan;
      Status s = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
      if (s) s = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
      return s;
    };
    ProductDriver seed, restored;
    Status s = create(seed);
    RestartExpected strict, compatible;
    if (s) s = seed.restart_expected(strict);
    if (s)
      s = seed.restart_expected(
          compatible, RestartStorageCompatibility::mg_bundle_ghost_v1);
    // Portable synthetic payload exercises the public compatibility contract;
    // an independently generated 69d8eee checkpoint is replayed separately.
    passed = s && strict.compatible_storage_plan == 0 &&
             compatible.compatible_storage_plan != 0 &&
             compatible.compatible_storage_schema != 0 &&
             compatible.compatible_storage_plan != strict.plan &&
             compatible.compatible_storage_schema != strict.schema;
    if (!passed && rank == 0)
      std::cerr << "missing explicit legacy MG storage identity\n";
    agree();
    DriverInitialState initial;
    initial.velocity = {0.11, -0.03, 0.02};
    if (passed) s = seed.initialize(initial);
    DriverStepReport step;
    for (int n = 0; n < 2 && passed && s; ++n)
      s = seed.advance({1, 1, 1, 1, 1}, step);
    RestartSnapshot snapshot;
    if (passed && s) s = seed.committed_restart_snapshot(snapshot);
    if (passed && s) {
      snapshot.plan = compatible.compatible_storage_plan;
      snapshot.schema = compatible.compatible_storage_schema;
      s = RestartWriter::write(MPI_COMM_WORLD, root, snapshot);
    }
    passed &= static_cast<bool>(s);
    RestartImage image;
    if (passed) {
      image.step = 991;
      const auto rejected =
          RestartReader::load(MPI_COMM_WORLD, root, strict, image);
      passed &= !rejected && rejected.detail == 10307U && image.step == 991;
      s = RestartReader::load(MPI_COMM_WORLD, root, compatible, image);
      passed &=
          s && image.storage_layout_migrated && image.plan == snapshot.plan &&
          image.schema == snapshot.schema && !image.backward_euler_recovery;
    }
    if (passed) s = create(restored);
    if (passed && s) {
      passed &= !restored.initialize_restart(image);
      auto wrong = image;
      wrong.plan ^= 1U;
      passed &= !restored.initialize_restart(
          wrong, RestartStorageCompatibility::mg_bundle_ghost_v1);
      s = restored.initialize_restart(
          image, RestartStorageCompatibility::mg_bundle_ghost_v1);
      passed &= static_cast<bool>(s);
    }
    RestartSnapshot migrated;
    if (passed) {
      s = restored.committed_restart_snapshot(migrated);
      passed &=
          s && migrated.plan == strict.plan &&
          migrated.schema == strict.schema && migrated.time == image.time &&
          migrated.dt == image.dt && migrated.step == image.step &&
          migrated.controller_state == image.controller_state &&
          migrated.pressure_reference == image.pressure_reference &&
          migrated.previous_pressure_reference ==
              image.previous_pressure_reference &&
          migrated.closed_mass_target == image.closed_mass_target &&
          same_fields(migrated.fields, image.fields) &&
          same_fields(migrated.previous_fields, image.previous_fields) &&
          same_fields(migrated.accepted_rate_fields,
                      image.accepted_rate_fields) &&
          same_fields(migrated.previous_rate_fields,
                      image.previous_rate_fields) &&
          same_flux(migrated.final_mass_flux, image.final_mass_flux) &&
          same_flux(migrated.previous_mass_flux, image.previous_mass_flux);
    }
    agree();
    if (passed) {
      const auto next_root = root / "current-layout";
      s = RestartWriter::write(MPI_COMM_WORLD, next_root, migrated);
      RestartImage next;
      if (s) s = RestartReader::load(MPI_COMM_WORLD, next_root, strict, next);
      passed &=
          s && !next.storage_layout_migrated && next.plan == strict.plan &&
          next.schema == strict.schema &&
          next.controller_state == image.controller_state &&
          next.fields[2].values == image.fields[2].values &&
          next.previous_fields[2].values == image.previous_fields[2].values;
    }
    agree();
    if (passed) {
      s = restored.advance({1, 1, 1, 1, 1}, step);
      passed &= s && step.accepted && step.accepted_step == image.step + 1 &&
                !step.temporal_method_fallback;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::error_code error;
    std::filesystem::remove_all(root, error);
  }
  int local = passed ? 1 : 0, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "restart storage compatibility " << (global ? "PASS" : "FAIL")
              << '\n';
  MPI_Finalize();
  return global ? 0 : 1;
}
