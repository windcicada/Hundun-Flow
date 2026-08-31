// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_io.hpp"

#include "../support/piso_fixture.hpp"
#include "io_restart_detail.hpp"

#include <mpi.h>
#include <unistd.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

constexpr Int3 kGlobal{7, 5, 3};
constexpr PlanFingerprint kPlan = UINT64_C(0x19010001);
constexpr PlanFingerprint kSchema = UINT64_C(0x19010002);
constexpr PlanFingerprint kGeometry = UINT64_C(0x19010003);

double cell_value(std::size_t field, std::uint8_t component,
                  Int3 global) noexcept {
  return 1000.0 * static_cast<double>(field + 1U) +
         100.0 * component + global.x + 0.01 * global.y + 0.0001 * global.z;
}

double face_value(std::size_t axis, Int3 global) noexcept {
  return 10000.0 * static_cast<double>(axis + 1U) + global.x +
         0.01 * global.y + 0.0001 * global.z;
}

struct Fixture {
  MeshPatch patch{};
  std::array<test::OwnedField, 4U> fields;
  std::array<RestartFieldView, 4U> restart_fields{};
  FaceFluxStorage provisional;
  FaceFluxView provisional_view{};
  FaceFluxStorage final_storage;
  StateLayers layers;
  AttemptTransaction transaction;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
  ConstFaceFluxView committed_flux{};
  ConstFaceFluxView previous_flux{};

  bool initialize(MPI_Comm communicator) {
    CartesianMeshSpec mesh;
    mesh.kind = GeometryKind::uniform;
    mesh.lower = {0.0, 0.0, 0.0};
    mesh.upper = {1.0, 1.0, 1.0};
    mesh.has_exact_cells = true;
    mesh.exact_cells = kGlobal;
    mesh.minimum_spacing = {1.0 / kGlobal.x, 1.0 / kGlobal.y,
                            1.0 / kGlobal.z};
    mesh.max_growth_ratio = 1.0;
    mesh.limits.max_global_cells = 105U;
    mesh.limits.max_memory_bytes_per_rank = UINT64_C(134217728);
    CartesianGeometryPlan geometry;
    if (!CartesianGeometryCompiler::compile(communicator, mesh, {}, geometry,
                                             patch))
      return false;
    const std::array<std::uint8_t, 4U> components{{3U, 1U, 1U, 1U}};
    const std::array<RestartFieldRole, 4U> roles{{
        RestartFieldRole::velocity,
        RestartFieldRole::pressure_perturbation,
        RestartFieldRole::enthalpy,
        RestartFieldRole::independent_species}};
    for (std::size_t field = 0U; field < fields.size(); ++field) {
      fields[field] = test::make_field(
          static_cast<FieldId>(field), patch.cells, components[field], 0U,
          static_cast<RevisionToken>(100U + field),
          static_cast<StorageIdentity>(200U + field));
      for (std::int32_t z = 0; z < patch.cells.z; ++z)
        for (std::int32_t y = 0; y < patch.cells.y; ++y)
          for (std::int32_t x = 0; x < patch.cells.x; ++x)
            for (std::uint8_t component = 0U;
                 component < components[field]; ++component) {
              const Int3 global{patch.begin.x + x, patch.begin.y + y,
                                patch.begin.z + z};
              fields[field].view.unchecked({x, y, z}, component) =
                  cell_value(field, component, global);
            }
      restart_fields[field] = {roles[field], as_const(fields[field].view)};
    }
    if (!FaceFluxStorage::allocate_workspace(patch.cells, 1U, provisional) ||
        !provisional.workspace_view(0U, 301U, provisional_view))
      return false;
    FaceFieldView faces[3]{provisional_view.x, provisional_view.y,
                           provisional_view.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis)
      for (std::int32_t z = 0; z < faces[axis].extents.z; ++z)
        for (std::int32_t y = 0; y < faces[axis].extents.y; ++y)
          for (std::int32_t x = 0; x < faces[axis].extents.x; ++x)
            faces[axis].unchecked({x, y, z}) = face_value(
                axis, {patch.begin.x + x, patch.begin.y + y,
                       patch.begin.z + z});

    FieldRegistry registry;
    FieldSchema schema;
    FieldId dummy = 0U;
    if (!registry.declare_field("restart_dummy", 1U, 0U, dummy) ||
        !registry.freeze(schema))
      return false;
    const std::array requests{
        ArenaFieldRequest{dummy, patch.cells, {0U}, FieldLifetime::state_layer}};
    ArenaLayout layout;
    if (!ArenaLayout::compile(schema, {requests.data(), requests.size()}, layout) ||
        !StateLayers::allocate(layout, layers) ||
        !AttemptTransaction::create(layers.field_count(), 1U, 1U,
                                    transaction) ||
        !authority.claim(50U, 0U, transaction, writer) ||
        !FaceFluxStorage::allocate_final(patch.cells, final_storage) ||
        !writer.initialize_committed(final_storage,
                                     as_const(provisional_view)) ||
        !writer.committed(final_storage, committed_flux))
      return false;
    return true;
  }

  bool promote_flux(MPI_Comm communicator) {
    if (!transaction.begin(layers) || !transaction.revise_trial(0U))
      return false;
    PendingFaceFluxView pending;
    const std::array dependencies{RevisionDependency{
        AttemptTransaction::field_revision_source(0U),
        transaction.trial_revision(0U)}};
    if (!writer.begin_pending(transaction, final_storage, pending) ||
        !detail::overwrite_pending_face_flux_for_test(pending, 777.0) ||
        !writer.publish_pending(
            {dependencies.data(), dependencies.size()}, pending) ||
        !transaction.collective_finish(communicator, Status{}) ||
        !writer.committed(final_storage, committed_flux) ||
        !writer.committed_previous(final_storage, previous_flux))
      return false;
    return committed_flux.revision == 2U && previous_flux.revision == 1U;
  }

  RestartSnapshot snapshot(std::uint64_t step) const noexcept {
    return {kGlobal,
            patch,
            kPlan,
            kSchema,
            kGeometry,
            0.01 * static_cast<double>(step),
            0.01,
            101325.0,
            step,
            UINT64_C(0x55aa),
            {restart_fields.data(), restart_fields.size()},
            committed_flux};
  }
};

struct ExactFixture {
  Fixture base;
  std::array<test::OwnedField, 4U> previous_fields;
  std::array<test::OwnedField, 2U> accepted_rates;
  std::array<test::OwnedField, 2U> previous_rates;
  std::array<RestartFieldView, 4U> previous_restart_fields{};
  std::array<RestartFieldView, 2U> accepted_restart_rates{};
  std::array<RestartFieldView, 2U> previous_restart_rates{};

  bool initialize(MPI_Comm communicator) {
    if (!base.initialize(communicator) || !base.promote_flux(communicator))
      return false;
    const std::array<std::uint8_t, 4U> components{{3U, 1U, 1U, 1U}};
    const std::array<RestartFieldRole, 4U> roles{{
        RestartFieldRole::velocity,
        RestartFieldRole::pressure_perturbation,
        RestartFieldRole::enthalpy,
        RestartFieldRole::independent_species}};
    for (std::size_t field = 0U; field < previous_fields.size(); ++field) {
      previous_fields[field] = test::make_field(
          static_cast<FieldId>(field), base.patch.cells, components[field],
          0U, static_cast<RevisionToken>(500U + field),
          static_cast<StorageIdentity>(600U + field));
      for (std::int32_t z = 0; z < base.patch.cells.z; ++z)
        for (std::int32_t y = 0; y < base.patch.cells.y; ++y)
          for (std::int32_t x = 0; x < base.patch.cells.x; ++x)
            for (std::uint8_t component = 0U;
                 component < components[field]; ++component) {
              const Int3 global{base.patch.begin.x + x,
                                base.patch.begin.y + y,
                                base.patch.begin.z + z};
              previous_fields[field].view.unchecked({x, y, z}, component) =
                  cell_value(field, component, global) + 50000.0;
            }
      previous_restart_fields[field] = {
          roles[field], as_const(previous_fields[field].view)};
    }
    const std::array<RestartFieldRole, 2U> rate_roles{{
        RestartFieldRole::enthalpy_nonadvective_rate,
        RestartFieldRole::scalar_nonadvective_rate}};
    for (std::size_t field = 0U; field < accepted_rates.size(); ++field) {
      const FieldId id = static_cast<FieldId>(10U + field);
      accepted_rates[field] = test::make_field(
          id, base.patch.cells, 1U, 0U,
          static_cast<RevisionToken>(700U + field),
          static_cast<StorageIdentity>(800U + field));
      previous_rates[field] = test::make_field(
          id, base.patch.cells, 1U, 0U,
          static_cast<RevisionToken>(900U + field),
          static_cast<StorageIdentity>(1000U + field));
      for (std::int32_t z = 0; z < base.patch.cells.z; ++z)
        for (std::int32_t y = 0; y < base.patch.cells.y; ++y)
          for (std::int32_t x = 0; x < base.patch.cells.x; ++x) {
            const Int3 global{base.patch.begin.x + x,
                              base.patch.begin.y + y,
                              base.patch.begin.z + z};
            accepted_rates[field].view.unchecked({x, y, z}, 0U) =
                70000.0 + 1000.0 * field + global.x + 0.01 * global.y;
            previous_rates[field].view.unchecked({x, y, z}, 0U) =
                80000.0 + 1000.0 * field + global.x + 0.01 * global.y;
          }
      accepted_restart_rates[field] = {
          rate_roles[field], as_const(accepted_rates[field].view)};
      previous_restart_rates[field] = {
          rate_roles[field], as_const(previous_rates[field].view)};
    }
    return true;
  }

  RestartSnapshot snapshot() const noexcept {
    RestartSnapshot value = base.snapshot(50U);
    value.controller_state = 51U;
    value.previous_fields = {previous_restart_fields.data(),
                             previous_restart_fields.size()};
    value.accepted_rate_fields = {accepted_restart_rates.data(),
                                  accepted_restart_rates.size()};
    value.previous_rate_fields = {previous_restart_rates.data(),
                                  previous_restart_rates.size()};
    value.previous_mass_flux = base.previous_flux;
    value.previous_pressure_reference = 101300.0;
    value.closed_mass_target = 3.5;
    return value;
  }
};

bool verify(const RestartImage& image, const MeshPatch& patch) {
  if (image.fields.size() != 4U || !image.backward_euler_recovery ||
      image.step == 0U || image.patch.begin.x != patch.begin.x ||
      image.patch.cells.x != patch.cells.x)
    return false;
  for (std::size_t index = 0U; index < kRuntimeSha256HexCharacters;
       ++index) {
    const char value = image.source_manifest_sha256[index];
    if (!((value >= '0' && value <= '9') ||
          (value >= 'a' && value <= 'f')))
      return false;
  }
  if (image.source_manifest_sha256[kRuntimeSha256HexCharacters] != '\0')
    return false;
  const auto index = [](Int3 local, Int3 cells) {
    return (static_cast<std::size_t>(local.z) * cells.y + local.y) * cells.x +
           local.x;
  };
  for (std::size_t field = 0U; field < image.fields.size(); ++field)
    for (std::int32_t z = 0; z < patch.cells.z; ++z)
      for (std::int32_t y = 0; y < patch.cells.y; ++y)
        for (std::int32_t x = 0; x < patch.cells.x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{patch.begin.x + x, patch.begin.y + y,
                            patch.begin.z + z};
          const std::size_t cell = index(local, patch.cells);
          for (std::uint8_t component = 0U;
               component < image.fields[field].components; ++component) {
            if (image.fields[field]
                    .values[cell * image.fields[field].components + component] !=
                cell_value(field, component, global))
              return false;
          }
        }
  const std::array<Int3, 3U> extents{{
      {patch.cells.x + 1, patch.cells.y, patch.cells.z},
      {patch.cells.x, patch.cells.y + 1, patch.cells.z},
      {patch.cells.x, patch.cells.y, patch.cells.z + 1}}};
  for (std::size_t axis = 0U; axis < 3U; ++axis)
    for (std::int32_t z = 0; z < extents[axis].z; ++z)
      for (std::int32_t y = 0; y < extents[axis].y; ++y)
        for (std::int32_t x = 0; x < extents[axis].x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{patch.begin.x + x, patch.begin.y + y,
                            patch.begin.z + z};
          if (image.final_mass_flux[axis][index(local, extents[axis])] !=
              face_value(axis, global))
            return false;
        }
  return true;
}

bool exact_transition(MPI_Comm communicator, const fs::path& directory) {
  ExactFixture fixture;
  bool passed = fixture.initialize(communicator);
  Status status;
  if (passed)
    status = RestartWriter::write(communicator, directory,
                                  fixture.snapshot(), {1U});
  passed = passed && static_cast<bool>(status);
  const std::array<RestartExpectedField, 4U> expected_fields{{
      {RestartFieldRole::velocity, 0U, 3U},
      {RestartFieldRole::pressure_perturbation, 1U, 1U},
      {RestartFieldRole::enthalpy, 2U, 1U},
      {RestartFieldRole::independent_species, 3U, 1U}}};
  const std::array<RestartExpectedField, 2U> expected_rates{{
      {RestartFieldRole::enthalpy_nonadvective_rate, 10U, 1U},
      {RestartFieldRole::scalar_nonadvective_rate, 11U, 1U}}};
  const RestartExpected expected{
      kGlobal,
      fixture.base.patch,
      kPlan,
      kSchema,
      kGeometry,
      {expected_fields.data(), expected_fields.size()},
      {expected_rates.data(), expected_rates.size()}};
  RestartImage image;
  if (passed)
    status = RestartReader::load(communicator, directory, expected, image);
  passed = passed && static_cast<bool>(status) &&
           !image.backward_euler_recovery &&
           image.controller_state == 51U &&
           image.previous_pressure_reference == 101300.0 &&
           image.closed_mass_target == 3.5 &&
           image.final_mass_flux_revision == 2U &&
           image.previous_mass_flux_revision == 1U &&
           image.previous_fields.size() == expected_fields.size() &&
           image.accepted_rate_fields.size() == expected_rates.size() &&
           image.previous_rate_fields.size() == expected_rates.size();
  const auto dense_index = [](Int3 local, Int3 cells) {
    return (static_cast<std::size_t>(local.z) * cells.y + local.y) * cells.x +
           local.x;
  };
  for (std::size_t field = 0U;
       field < image.previous_fields.size() && passed; ++field)
    for (std::int32_t z = 0; z < fixture.base.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < fixture.base.patch.cells.y; ++y)
        for (std::int32_t x = 0; x < fixture.base.patch.cells.x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{fixture.base.patch.begin.x + x,
                            fixture.base.patch.begin.y + y,
                            fixture.base.patch.begin.z + z};
          const std::size_t cell =
              dense_index(local, fixture.base.patch.cells);
          for (std::uint8_t component = 0U;
               component < image.previous_fields[field].components;
               ++component)
            passed = passed &&
                     image.previous_fields[field]
                             .values[cell *
                                         image.previous_fields[field]
                                             .components +
                                     component] ==
                         cell_value(field, component, global) + 50000.0;
        }
  for (std::size_t field = 0U;
       field < image.accepted_rate_fields.size() && passed; ++field)
    for (std::int32_t z = 0; z < fixture.base.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < fixture.base.patch.cells.y; ++y)
        for (std::int32_t x = 0; x < fixture.base.patch.cells.x; ++x) {
          const Int3 global{fixture.base.patch.begin.x + x,
                            fixture.base.patch.begin.y + y,
                            fixture.base.patch.begin.z + z};
          const std::size_t cell = dense_index({x, y, z},
                                               fixture.base.patch.cells);
          passed = image.accepted_rate_fields[field].values[cell] ==
                       70000.0 + 1000.0 * field + global.x +
                           0.01 * global.y &&
                   image.previous_rate_fields[field].values[cell] ==
                       80000.0 + 1000.0 * field + global.x +
                           0.01 * global.y;
        }
  const std::array<Int3, 3U> extents{{
      {fixture.base.patch.cells.x + 1, fixture.base.patch.cells.y,
       fixture.base.patch.cells.z},
      {fixture.base.patch.cells.x, fixture.base.patch.cells.y + 1,
       fixture.base.patch.cells.z},
      {fixture.base.patch.cells.x, fixture.base.patch.cells.y,
       fixture.base.patch.cells.z + 1}}};
  for (std::size_t axis = 0U; axis < extents.size() && passed; ++axis)
    for (std::int32_t z = 0; z < extents[axis].z; ++z)
      for (std::int32_t y = 0; y < extents[axis].y; ++y)
        for (std::int32_t x = 0; x < extents[axis].x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{fixture.base.patch.begin.x + x,
                            fixture.base.patch.begin.y + y,
                            fixture.base.patch.begin.z + z};
          const std::size_t index = dense_index(local, extents[axis]);
          passed = passed &&
                   image.final_mass_flux[axis][index] == 777.0 &&
                   image.previous_mass_flux[axis][index] ==
                       face_value(axis, global);
        }
  const int local = passed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
  return global != 0;
}

bool transition(MPI_Comm world, int writer_size, int reader_size,
                const fs::path& directory, std::uint64_t step) {
  int world_rank = 0;
  MPI_Comm_rank(world, &world_rank);
  MPI_Comm writer = MPI_COMM_NULL;
  MPI_Comm_split(world, world_rank < writer_size ? 0 : MPI_UNDEFINED,
                 world_rank, &writer);
  bool local = true;
  if (writer != MPI_COMM_NULL) {
    Fixture fixture;
    local = fixture.initialize(writer);
    Status write_status;
    if (local)
      write_status = RestartWriter::write(writer, directory,
                                          fixture.snapshot(step), {1U});
    if (!write_status) {
      std::cerr << "world rank " << world_rank << " transition "
                << writer_size << "->" << reader_size << " write status="
                << static_cast<unsigned>(write_status.code)
                << " detail=" << write_status.detail << '\n';
    }
    local = local && static_cast<bool>(write_status);
    MPI_Comm_free(&writer);
  }
  MPI_Barrier(world);
  MPI_Comm reader = MPI_COMM_NULL;
  MPI_Comm_split(world, world_rank < reader_size ? 0 : MPI_UNDEFINED,
                 world_rank, &reader);
  if (reader != MPI_COMM_NULL) {
    Fixture fixture;
    local = local && fixture.initialize(reader);
    const std::array<RestartExpectedField, 4U> expected_fields{{
        {RestartFieldRole::velocity, 0U, 3U},
        {RestartFieldRole::pressure_perturbation, 1U, 1U},
        {RestartFieldRole::enthalpy, 2U, 1U},
        {RestartFieldRole::independent_species, 3U, 1U}}};
    RestartExpected expected{kGlobal,
                             fixture.patch,
                             kPlan,
                             kSchema,
                             kGeometry,
                             {expected_fields.data(), expected_fields.size()}};
    RestartImage image;
    Status read_status;
    if (local)
      read_status = RestartReader::load(reader, directory, expected, image);
    if (!read_status) {
      std::cerr << "world rank " << world_rank << " transition "
                << writer_size << "->" << reader_size << " read status="
                << static_cast<unsigned>(read_status.code)
                << " detail=" << read_status.detail << '\n';
    }
    local = local && static_cast<bool>(read_status) &&
            verify(image, fixture.patch);
    MPI_Comm_free(&reader);
  }
  MPI_Barrier(world);
  const int value = local ? 1 : 0;
  int passed = 0;
  MPI_Allreduce(&value, &passed, 1, MPI_INT, MPI_MIN, world);
  return passed != 0;
}

bool failure_boundaries(MPI_Comm communicator, const fs::path& directory) {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  Fixture fixture;
  bool passed = fixture.initialize(communicator);
  Status status;
  if (passed)
    status = RestartWriter::write(communicator, directory,
                                  fixture.snapshot(20U), {1U});
  passed = passed && static_cast<bool>(status);
  const std::array points{
      detail::RestartFailurePoint::after_directory,
      detail::RestartFailurePoint::after_rank_file,
      detail::RestartFailurePoint::after_manifest,
      detail::RestartFailurePoint::after_generation_rename,
      detail::RestartFailurePoint::after_current_switch};
  const std::array<RestartExpectedField, 4U> expected_fields{{
      {RestartFieldRole::velocity, 0U, 3U},
      {RestartFieldRole::pressure_perturbation, 1U, 1U},
      {RestartFieldRole::enthalpy, 2U, 1U},
      {RestartFieldRole::independent_species, 3U, 1U}}};
  const RestartExpected expected{kGlobal,
                                 fixture.patch,
                                 kPlan,
                                 kSchema,
                                 kGeometry,
                                 {expected_fields.data(),
                                  expected_fields.size()}};
  for (std::size_t index = 0U; index < points.size() && passed; ++index) {
    detail::set_restart_failure_for_test(points[index], 0);
    status = RestartWriter::write(communicator, directory,
                                  fixture.snapshot(21U + index), {1U});
    detail::clear_restart_failure_for_test();
    passed &= status.code == StatusCode::io_failure;
    RestartImage image;
    status = RestartReader::load(communicator, directory, expected, image);
    passed &= static_cast<bool>(status) && verify(image, fixture.patch);
  }
  status = RestartWriter::write(communicator, directory,
                                fixture.snapshot(40U), {1U});
  passed &= static_cast<bool>(status);
  if (rank == 0) {
    std::error_code error;
    std::size_t generations = 0U;
    std::size_t pending = 0U;
    for (fs::directory_iterator iterator(directory, error), end;
         !error && iterator != end; iterator.increment(error)) {
      const std::string name = iterator->path().filename().string();
      if (iterator->is_directory(error) &&
          name.rfind("generation-", 0U) == 0U) {
        ++generations;
        if (name.find("-pending") != std::string::npos) ++pending;
      }
    }
    passed &= !error && generations == 1U && pending == 0U;
  }
  RestartExpected mismatch = expected;
  mismatch.plan = kPlan + 1U;
  RestartImage unchanged;
  unchanged.plan = UINT64_C(0xdeadbeef);
  status = RestartReader::load(communicator, directory, mismatch, unchanged);
  passed &= status.code == StatusCode::invalid_plan &&
            unchanged.plan == UINT64_C(0xdeadbeef);
  if (rank == 0) {
    std::ifstream pointer(directory / "current", std::ios::binary);
    std::string generation;
    std::getline(pointer, generation);
    std::fstream rank_file(directory / generation / "rank-00000000.bin",
                           std::ios::in | std::ios::out | std::ios::binary);
    rank_file.seekg(20);
    char byte = 0;
    rank_file.read(&byte, 1);
    byte ^= 0x5a;
    rank_file.seekp(20);
    rank_file.write(&byte, 1);
    rank_file.flush();
  }
  MPI_Barrier(communicator);
  unchanged.plan = UINT64_C(0xcafebabe);
  status = RestartReader::load(communicator, directory, expected, unchanged);
  passed &= status.code == StatusCode::io_failure &&
            unchanged.plan == UINT64_C(0xcafebabe);
  const int value = passed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, communicator);
  return global != 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 4) {
    if (rank == 0) std::cerr << "io_restart_mpi_test requires four ranks\n";
    MPI_Finalize();
    return 2;
  }
  std::string root;
  if (rank == 0)
    root = (fs::temp_directory_path() /
            ("hundun-v04-restart-" + std::to_string(::getpid())))
               .string();
  std::uint64_t length = root.size();
  MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (rank != 0) root.resize(length);
  MPI_Bcast(root.data(), static_cast<int>(length), MPI_CHAR, 0, MPI_COMM_WORLD);
  const fs::path base(root);
  if (rank == 0) {
    std::error_code error;
    fs::remove_all(base, error);
    fs::create_directories(base);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  bool passed = true;
  passed &= transition(MPI_COMM_WORLD, 1, 2, base / "one-to-two", 1U);
  passed &= transition(MPI_COMM_WORLD, 2, 4, base / "two-to-four", 2U);
  passed &= transition(MPI_COMM_WORLD, 4, 1, base / "four-to-one", 3U);
  passed &= exact_transition(MPI_COMM_WORLD, base / "exact-four-to-four");
  passed &= failure_boundaries(MPI_COMM_WORLD, base / "failure-boundaries");
  if (rank == 0) {
    std::error_code error;
    fs::remove_all(base, error);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  const int value = passed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0 && global == 0) std::cerr << "restart rank-change failure\n";
  MPI_Finalize();
  return global != 0 ? 0 : 1;
}
