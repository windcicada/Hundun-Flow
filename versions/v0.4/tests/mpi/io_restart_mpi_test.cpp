// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_io.hpp"

#include "../support/piso_fixture.hpp"
#include "core_spray_history_detail.hpp"
#include "core_tcr_history_detail.hpp"
#include "io_restart_detail.hpp"

#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/syscall.h>
#include <cerrno>
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <new>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Avoid sys/stat.h's optimized inline fstat definition: this test interposes
// fstat/__fxstat themselves below. Use the ordinary libc FIFO entry point.
extern "C" int mkfifo(const char*, mode_t) noexcept;

namespace restart_allocation_probe {
bool armed = false;
unsigned oversized = 0U;
}
void* operator new(std::size_t bytes) {
  if (restart_allocation_probe::armed && bytes > 16U * 1024U * 1024U) {
    ++restart_allocation_probe::oversized;
    throw std::bad_alloc();  // Never exhaust memory to test an untrusted length.
  }
  if (void* value = std::malloc(bytes == 0U ? 1U : bytes)) return value;
  throw std::bad_alloc();
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* value) noexcept { std::free(value); }
void operator delete[](void* value) noexcept { std::free(value); }
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }

namespace restart_syscall_probe {
bool enabled = false, fired = false, close_fired = false, current_renamed = false;
int mode = 0, descriptor = -1;
std::string prefix;
}
extern "C" int open(const char* path, int flags, ...) {
  mode_t permissions = 0;
  if ((flags & O_CREAT) != 0) {
    va_list args;
    va_start(args, flags);
    permissions = va_arg(args, int);
    va_end(args);
  }
  using namespace restart_syscall_probe;
  const bool watched = enabled && ((flags & O_DIRECTORY) == 0 ||
                                   (mode == 10 && current_renamed)) &&
      std::strncmp(path, prefix.c_str(), prefix.size()) == 0;
  if (watched && !fired && (mode == 0 || mode == 4)) {
    fired = true;
    errno = mode == 0 ? EINTR : EACCES;
    return -1;
  }
  const int result = static_cast<int>(
      ::syscall(SYS_openat, AT_FDCWD, path, flags, permissions));
  if (watched) descriptor = result;
  return result;
}
extern "C" ssize_t write(int fd, const void* data, std::size_t count) {
  using namespace restart_syscall_probe;
  if (enabled && fd == descriptor && !fired &&
      (mode == 1 || mode == 3 || mode == 5 || mode == 8 || mode == 9)) {
    fired = true;
    if (mode == 3) return ::syscall(SYS_write, fd, data, std::max(std::size_t{1}, count / 2));
    if (mode == 9) return 0;
    errno = mode == 1 ? EINTR : ENOSPC;
    return -1;
  }
  return ::syscall(SYS_write, fd, data, count);
}
extern "C" int fsync(int fd) {
  using namespace restart_syscall_probe;
  if (enabled && fd == descriptor && !fired &&
      (mode == 2 || mode == 6 || (mode == 10 && current_renamed))) {
    fired = true;
    errno = mode == 2 ? EINTR : EIO;
    return -1;
  }
  return static_cast<int>(::syscall(SYS_fsync, fd));
}
extern "C" int close(int fd) {
  using namespace restart_syscall_probe;
  const bool fail = enabled && fd == descriptor && !close_fired &&
                    (mode == 7 || mode == 8);
  const int result = static_cast<int>(::syscall(SYS_close, fd));
  if (fd == descriptor) descriptor = -1;
  if (fail) {
    fired = close_fired = true;
    errno = EIO;
    return -1;
  }
  return result;
}
extern "C" ssize_t read(int fd, void* data, std::size_t count) {
  using namespace restart_syscall_probe;
  if (enabled && fd == descriptor && !fired && (mode == 11 || mode == 12)) {
    fired = true;
    errno = mode == 11 ? EINTR : EIO;
    return -1;
  }
  return ::syscall(SYS_read, fd, data, count);
}
static int observed_fstat(int fd, struct stat* out) {
  using namespace restart_syscall_probe;
  if (enabled && fd == descriptor && !fired && (mode == 13 || mode == 14)) {
    fired = true;
    errno = mode == 13 ? EINTR : EIO;
    return -1;
  }
  const int result = static_cast<int>(::syscall(SYS_fstat, fd, out));
  if (result == 0 && enabled && fd == descriptor && !fired && mode == 15) {
    fired = true;
    out->st_size += static_cast<off_t>(UINT64_C(1) << 40U);
  }
  return result;
}
extern "C" int fstat(int fd, struct stat* out) {
  return observed_fstat(fd, out);
}
extern "C" int __fxstat(int, int fd, struct stat* out) {
  return observed_fstat(fd, out);
}
extern "C" int rename(const char* from, const char* to) {
  const int result = static_cast<int>(::syscall(SYS_renameat, AT_FDCWD, from, AT_FDCWD, to));
  using namespace restart_syscall_probe;
  if (result == 0 && enabled && std::strncmp(to, prefix.c_str(), prefix.size()) == 0 &&
      std::strlen(to) >= 8U && std::strcmp(to + std::strlen(to) - 8U, "/current") == 0)
    current_renamed = true;
  return result;
}

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

bool read_budget_boundaries(MPI_Comm communicator, const fs::path& directory,
                            const RestartExpected& expected, RestartImage& image) {
  int rank = 0, ranks = 0;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &ranks);
  RestartReadReport measured;
  bool passed = static_cast<bool>(RestartReader::load(communicator, directory, expected, image, &measured));
  passed &= measured.retained_image_bytes == measured.new_image_bytes &&
            measured.peak_bulk_bytes > 2U * measured.new_image_bytes;
  RestartReadLimits limits;
  limits.maximum_bulk_bytes = measured.peak_bulk_bytes;
  passed &= static_cast<bool>(RestartReader::load(communicator, directory, expected, image, nullptr, limits));
  const auto saved = image;
  const auto same_payload = [&] {
    bool equal =
        image.step == saved.step && image.time == saved.time &&
        image.dt == saved.dt &&
        image.closed_mass_target == saved.closed_mass_target &&
        image.plan == saved.plan && image.schema == saved.schema &&
        image.geometry == saved.geometry &&
        image.controller_state == saved.controller_state &&
        image.cell_record_identity == saved.cell_record_identity &&
        image.cell_record_bytes == saved.cell_record_bytes &&
        image.cell_record_lengths == saved.cell_record_lengths &&
        image.cell_records == saved.cell_records &&
        image.source_manifest_sha256 == saved.source_manifest_sha256 &&
        image.method_history_signature == saved.method_history_signature &&
        image.final_mass_flux == saved.final_mass_flux &&
        image.previous_mass_flux == saved.previous_mass_flux;
    const std::array<const std::vector<RestartImageField>*, 4U> a{{&image.fields, &image.previous_fields,
        &image.accepted_rate_fields, &image.previous_rate_fields}};
    const std::array<const std::vector<RestartImageField>*, 4U> b{{&saved.fields, &saved.previous_fields,
        &saved.accepted_rate_fields, &saved.previous_rate_fields}};
    for (std::size_t layer = 0; layer < a.size(); ++layer) {
      equal &= a[layer]->size() == b[layer]->size();
      for (std::size_t f = 0; equal && f < a[layer]->size(); ++f)
        equal &= (*a[layer])[f].role == (*b[layer])[f].role &&
            (*a[layer])[f].field == (*b[layer])[f].field &&
            (*a[layer])[f].components == (*b[layer])[f].components &&
            (*a[layer])[f].values == (*b[layer])[f].values;
    }
    return equal;
  };
  for (int mode = 0; mode < 4; ++mode) {
    limits = {};
    if (rank == ranks - 1) {
      if (mode == 0) limits.maximum_bulk_bytes = measured.peak_bulk_bytes - 1U;
      if (mode == 1) limits.maximum_manifest_bytes = 1U;
      if (mode == 2) limits.maximum_source_ranks = static_cast<std::uint32_t>(ranks - 1);
      if (mode == 3) limits.maximum_bulk_bytes = 0U;
    }
    const auto* before = image.fields[0U].values.data();
    const Status denied = RestartReader::load(communicator, directory, expected, image, nullptr, limits);
    const std::uint64_t wire = (std::uint64_t(denied.code) << 32U) | denied.detail;
    auto low = wire, high = wire;
    MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_UINT64_T, MPI_MIN, communicator);
    MPI_Allreduce(MPI_IN_PLACE, &high, 1, MPI_UINT64_T, MPI_MAX, communicator);
    passed &= !denied && low == high && before == image.fields[0U].values.data() && same_payload();
    if (mode == 0) passed &= denied.code == StatusCode::allocation_failure && denied.detail == 10310U;
  }
  if (rank == 0) std::cout << "read_budget version=" << image.source_format_version
      << " new=" << measured.new_image_bytes << " retained=" << measured.retained_image_bytes
      << " peak=" << measured.peak_bulk_bytes << " passed=" << passed << '\n';
  return passed;
}

std::vector<std::uint8_t> model_records(const MeshPatch &patch) {
  std::vector<std::uint8_t> records(std::size_t(patch.cells.x) * patch.cells.y *
                                    patch.cells.z * 8U);
  std::size_t cell = 0U;
  for (int z = 0; z < patch.cells.z; ++z)
    for (int y = 0; y < patch.cells.y; ++y)
      for (int x = 0; x < patch.cells.x; ++x, ++cell) {
        const auto global =
            (std::uint64_t(z + patch.begin.z) * kGlobal.y + y + patch.begin.y) *
                kGlobal.x +
            x + patch.begin.x;
        const auto counter = UINT64_MAX - global;
        for (unsigned byte = 0; byte < 8; ++byte)
          records[8U * cell + byte] = std::uint8_t(counter >> (8U * byte));
      }
  return records;
}

bool exact_transition(MPI_Comm communicator, const fs::path &directory,
                      PlanFingerprint signature = 0U,
                      bool cell_records = false) {
  ExactFixture fixture;
  bool passed = fixture.initialize(communicator);
  Status status;
  auto snapshot = fixture.snapshot();
  snapshot.method_history_signature = signature;
  std::vector<std::uint8_t> records;
  if (cell_records) {
    records = model_records(fixture.base.patch);
    snapshot.cell_records = {
        UINT64_C(0x5443525245433031), 8U, {records.data(), records.size()}};
  }
  if (passed)
    status = RestartWriter::write(communicator, directory,
                                  snapshot, {1U});
  passed = passed && static_cast<bool>(status);
  const std::array<RestartExpectedField, 4U> expected_fields{{
      {RestartFieldRole::velocity, 0U, 3U},
      {RestartFieldRole::pressure_perturbation, 1U, 1U},
      {RestartFieldRole::enthalpy, 2U, 1U},
      {RestartFieldRole::independent_species, 3U, 1U}}};
  const std::array<RestartExpectedField, 2U> expected_rates{{
      {RestartFieldRole::enthalpy_nonadvective_rate, 10U, 1U},
      {RestartFieldRole::scalar_nonadvective_rate, 11U, 1U}}};
  RestartExpected expected{kGlobal,
                           fixture.base.patch,
                           kPlan,
                           kSchema,
                           kGeometry,
                           {expected_fields.data(), expected_fields.size()},
                           {expected_rates.data(), expected_rates.size()}};
  if (cell_records) {
    expected.cell_record_identity = snapshot.cell_records.identity;
    expected.cell_record_bytes = 8U;
  }
  RestartImage image;
  if (passed)
    status = RestartReader::load(communicator, directory, expected, image);
  passed = passed && static_cast<bool>(status) &&
           image.source_format_version ==
               (cell_records ? 4U : (signature == 0U ? 2U : 3U)) &&
           image.method_history_signature == signature &&
           image.history_compatibility(signature) ==
               (signature == 0U ? RestartHistoryCompatibility::unknown
                                : RestartHistoryCompatibility::compatible) &&
           !image.backward_euler_recovery && image.controller_state == 51U &&
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
  if (cell_records)
    passed &= image.cell_record_identity == snapshot.cell_records.identity &&
              image.cell_record_bytes == 8U && image.cell_records == records;
  passed &= read_budget_boundaries(communicator, directory, expected, image);
  const int local = passed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
  return global != 0;
}

std::vector<std::uint8_t> variable_records(const MeshPatch &patch,
                                           std::vector<std::uint32_t> &lengths,
                                           bool empty = false) {
  auto fixed = model_records(patch);
  std::vector<std::uint8_t> values;
  lengths.clear();
  std::size_t cell = 0;
  for (int z = 0; z < patch.cells.z; ++z)
    for (int y = 0; y < patch.cells.y; ++y)
      for (int x = 0; x < patch.cells.x; ++x, ++cell) {
        const auto global =
            (std::uint64_t(z + patch.begin.z) * kGlobal.y + y + patch.begin.y) *
                kGlobal.x +
            x + patch.begin.x;
        const unsigned records = empty ? 0 : global % 5;
        lengths.push_back(8U * records);
        for (unsigned r = 0; r < records; ++r)
          values.insert(values.end(), fixed.begin() + 8 * cell,
                        fixed.begin() + 8 * cell + 8);
      }
  return values;
}

struct SprayHistoryFixture {
  detail::ProductSprayHistory history;
  detail::ProductTcrHistory tcr;
  spray::detail::DeterministicInjector injector;
  std::vector<detail::ProductSprayHistory::Parcel> parcels;
  bool owns_injector{};
  bool initialize(MeshPatch patch, bool advance) {
    const auto count =
        std::size_t(patch.cells.x) * patch.cells.y * patch.cells.z;
    tcr.configure(991, count, -1);
    owns_injector =
        patch.begin.x == 0 && patch.begin.y == 0 && patch.begin.z == 0;
    spray::detail::DeterministicInjector *pointer = &injector;
    if (owns_injector) {
      spray::detail::InjectorSpec spec;
      spec.seed = 93;
      spec.injector_id = 77;
      spec.axis = {1, 0, 0};
      spec.mass_flow_rate_kg_per_s = .5;
      spec.represented_mass_per_parcel_kg = 1;
      spec.droplet_mass_kg = 1e-9;
      spec.droplet_diameter_m = 1e-4;
      spec.temperature_k = 300;
      spec.liquid_material_fingerprint = 123;
      spec.owner_global_cell = 0;
      if (!injector.reserve(2) ||
          !injector.configure(spec, {.25, UINT64_C(9007199254740993)}))
        return false;
    }
    if (!history.configure(
            990, patch, kGlobal, 2 * count, 123,
            {owns_injector ? &pointer : nullptr, owns_injector ? 1U : 0U},
            tcr.snapshot(), 1U << 24))
      return false;
    if (!advance)
      return true;
    for (int z = 0; z < patch.cells.z; ++z)
      for (int y = 0; y < patch.cells.y; ++y)
        for (int x = 0; x < patch.cells.x; ++x) {
          const std::uint64_t global =
              (x + patch.begin.x) +
              kGlobal.x *
                  ((y + patch.begin.y) + kGlobal.y * (z + patch.begin.z));
          for (unsigned j = 0; j < global % 3; ++j) {
            detail::ProductSprayHistory::Parcel v;
            v.parcel.id = {UINT64_C(9007199254741001) + global, j + 1};
            v.parcel.position_m = {(x + patch.begin.x + .5) / kGlobal.x,
                                   (y + patch.begin.y + .5) / kGlobal.y,
                                   (z + patch.begin.z + .5) / kGlobal.z};
            v.parcel.velocity_m_per_s = {.1, -.2, .3};
            v.parcel.droplet_mass_kg = 1e-9;
            v.parcel.droplet_diameter_m = 1e-4;
            v.parcel.multiplicity = 2;
            v.parcel.temperature_k = 300 + .1 * global;
            v.parcel.liquid_material_fingerprint = 123;
            v.parcel.owner_global_cell = global;
            v.parcel.age_s = .01;
            v.tab_deformation = -.25;
            v.tab_deformation_rate_per_s = global + .5;
            v.breakup_ordinal = UINT64_MAX - global - j;
            parcels.push_back(v);
          }
        }
    for (std::size_t i = 0; i < count; ++i) {
      tcr::detail::TrialRequest request;
      request.expected_revision = tcr.accepted(i).revision;
      request.mode = tcr::detail::Mode::experimental;
      request.initialization_sign = -1;
      request.mapping = {tcr::detail::Status::success,
                         {.25, 1},
                         tcr::detail::kReactantMoleFractionMappingIdentity};
      if (!tcr.stage(i, tcr::detail::prepare(tcr.accepted(i), request), 0))
        return false;
    }
    tcr.seal();
    if (owns_injector && !injector.begin_trial(0, 1).succeeded())
      return false;
    if (!history.stage_next({parcels.data(), parcels.size()}, 0,
                            tcr.prepared_snapshot()) ||
        !history.preflight_commit())
      return false;
    history.commit();
    tcr.commit();
    return true;
  }
};

bool spray_history_repartition(MPI_Comm world, int writers, int readers,
                               const fs::path &directory) {
  int rank = 0;
  MPI_Comm_rank(world, &rank);
  MPI_Comm writer = MPI_COMM_NULL, reader = MPI_COMM_NULL;
  MPI_Comm_split(world, rank < writers ? 0 : MPI_UNDEFINED, rank, &writer);
  bool passed = true;
  if (writer != MPI_COMM_NULL) {
    ExactFixture fields;
    SprayHistoryFixture spray;
    passed =
        fields.initialize(writer) && spray.initialize(fields.base.patch, true);
    int all_ready = passed ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &all_ready, 1, MPI_INT, MPI_MIN, writer);
    passed = all_ready != 0;
    if (passed) {
      auto snapshot = fields.snapshot();
      snapshot.step = 1;
      snapshot.time = .01;
      snapshot.method_history_signature = 993;
      snapshot.cell_records = spray.history.snapshot();
      passed = bool(RestartWriter::write(writer, directory, snapshot));
    }
    MPI_Comm_free(&writer);
  }
  MPI_Barrier(world);
  MPI_Comm_split(world, rank < readers ? 0 : MPI_UNDEFINED, rank, &reader);
  if (reader != MPI_COMM_NULL) {
    ExactFixture fixture;
    SprayHistoryFixture restored, reference;
    passed &= fixture.initialize(reader) &&
              restored.initialize(fixture.base.patch, false) &&
              reference.initialize(fixture.base.patch, true);
    const std::array<RestartExpectedField, 4> fields{
        {{RestartFieldRole::velocity, 0, 3},
         {RestartFieldRole::pressure_perturbation, 1, 1},
         {RestartFieldRole::enthalpy, 2, 1},
         {RestartFieldRole::independent_species, 3, 1}}};
    const std::array<RestartExpectedField, 2> rates{
        {{RestartFieldRole::enthalpy_nonadvective_rate, 10, 1},
         {RestartFieldRole::scalar_nonadvective_rate, 11, 1}}};
    RestartExpected expected{kGlobal,
                             fixture.base.patch,
                             kPlan,
                             kSchema,
                             kGeometry,
                             {fields.data(), fields.size()},
                             {rates.data(), rates.size()}};
    expected.cell_record_identity = 990;
    expected.cell_record_bytes = 0;
    RestartImage image;
    const auto read = RestartReader::load(reader, directory, expected, image);
    passed &= read && image.history_compatibility(993) ==
                          RestartHistoryCompatibility::compatible;
    if (read) {
      const auto prepared = restored.history.stage_restore(image);
      const auto tcr = restored.tcr.stage_restore_records(
          restored.history.prepared_tcr_records(), image.step);
      passed &= prepared && tcr && restored.history.preflight_commit();
      if (prepared && tcr && restored.history.preflight_commit()) {
        restored.history.commit();
        restored.tcr.commit();
        const auto a = restored.history.snapshot(),
                   b = reference.history.snapshot();
        passed &=
            a.values.size == b.values.size &&
            a.variable_cell_bytes.size == b.variable_cell_bytes.size &&
            std::equal(a.values.data, a.values.data + a.values.size,
                       b.values.data) &&
            std::equal(a.variable_cell_bytes.data,
                       a.variable_cell_bytes.data + a.variable_cell_bytes.size,
                       b.variable_cell_bytes.data);
        if (restored.owns_injector)
          passed &= restored.injector.committed_state() ==
                    reference.injector.committed_state();
      }
    }
    MPI_Comm_free(&reader);
  }
  MPI_Barrier(world);
  int okay = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, world);
  if (!okay && rank == 0)
    std::cerr << "FAIL: native spray/TCR variable-cell repartition " << writers
              << "->" << readers << '\n';
  return okay != 0;
}

bool record_repartition(MPI_Comm world, int writer_size, int reader_size,
                        const fs::path &directory, bool variable = false,
                        bool empty = false) {
  int rank = 0;
  MPI_Comm_rank(world, &rank);
  MPI_Comm writer = MPI_COMM_NULL, reader = MPI_COMM_NULL;
  MPI_Comm_split(world, rank < writer_size ? 0 : MPI_UNDEFINED, rank, &writer);
  bool passed = true;
  constexpr PlanFingerprint identity = UINT64_C(0x5443525245433031);
  if (writer != MPI_COMM_NULL) {
    ExactFixture fixture;
    passed = fixture.initialize(writer);
    auto records = model_records(fixture.base.patch);
    std::vector<std::uint32_t> lengths;
    if (variable)
      records = variable_records(fixture.base.patch, lengths, empty);
    auto snapshot = fixture.snapshot();
    snapshot.method_history_signature = UINT64_C(0x391003);
    snapshot.cell_records = {
        identity, variable ? 0U : 8U, {records.data(), records.size()}};
    if (variable)
      snapshot.cell_records.variable_cell_bytes = {lengths.data(),
                                                   lengths.size()};
    if (variable) {
      int writer_rank = 0;
      MPI_Comm_rank(writer, &writer_rank);
      if (writer_rank == writer_size - 1)
        ++lengths[0];
      passed &= !RestartWriter::write(writer, directory, snapshot);
      if (writer_rank == writer_size - 1)
        --lengths[0];
    }
    if (passed)
      passed =
          static_cast<bool>(RestartWriter::write(writer, directory, snapshot));
    MPI_Comm_free(&writer);
  }
  MPI_Barrier(world);
  MPI_Comm_split(world, rank < reader_size ? 0 : MPI_UNDEFINED, rank, &reader);
  if (reader != MPI_COMM_NULL) {
    ExactFixture fixture;
    passed &= fixture.initialize(reader);
    const std::array<RestartExpectedField, 4U> fields{
        {{RestartFieldRole::velocity, 0U, 3U},
         {RestartFieldRole::pressure_perturbation, 1U, 1U},
         {RestartFieldRole::enthalpy, 2U, 1U},
         {RestartFieldRole::independent_species, 3U, 1U}}};
    const std::array<RestartExpectedField, 2U> rates{
        {{RestartFieldRole::enthalpy_nonadvective_rate, 10U, 1U},
         {RestartFieldRole::scalar_nonadvective_rate, 11U, 1U}}};
    RestartExpected expected{kGlobal,
                             fixture.base.patch,
                             kPlan,
                             kSchema,
                             kGeometry,
                             {fields.data(), fields.size()},
                             {rates.data(), rates.size()}};
    expected.cell_record_identity = identity;
    expected.cell_record_bytes = variable ? 0U : 8U;
    std::vector<std::uint32_t> lengths;
    RestartImage image;
    const auto status = RestartReader::load(reader, directory, expected, image);
    passed &=
        status &&
        image.cell_records ==
            (variable ? variable_records(fixture.base.patch, lengths, empty)
                      : model_records(fixture.base.patch)) &&
        image.cell_record_identity == identity &&
        image.cell_record_bytes == (variable ? 0U : 8U) &&
        image.history_compatibility(UINT64_C(0x391003)) ==
            RestartHistoryCompatibility::compatible;
    if (variable)
      passed &= image.cell_record_lengths == lengths &&
                image.source_format_version == 5U;
    if (variable && writer_size == reader_size)
      passed &= read_budget_boundaries(reader, directory, expected, image);
    const auto saved = image.cell_records;
    for (int mismatch = 0; mismatch < 2; ++mismatch) {
      auto wrong = expected;
      if (rank == reader_size - 1) {
        if (mismatch == 0)
          wrong.cell_record_identity ^= 1U;
        else
          ++wrong.cell_record_bytes;
      }
      const auto denied = RestartReader::load(reader, directory, wrong, image);
      passed &= !denied && image.cell_records == saved;
    }
    MPI_Comm_free(&reader);
  }
  MPI_Barrier(world);
  int okay = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, world);
  return okay != 0;
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
    RestartReadReport observation;
    Status read_status;
    if (local)
      read_status = RestartReader::load(reader, directory, expected, image, &observation);
    if (!read_status) {
      std::cerr << "world rank " << world_rank << " transition "
                << writer_size << "->" << reader_size << " read status="
                << static_cast<unsigned>(read_status.code)
                << " detail=" << read_status.detail << '\n';
    }
    local = local && static_cast<bool>(read_status) &&
            verify(image, fixture.patch);
    unsigned integrity_count = 0U, skipped = 0U;
    const unsigned local_skipped = writer_size - observation.restoration_blocks;
    MPI_Allreduce(&observation.integrity_blocks, &integrity_count, 1,
                  MPI_UNSIGNED, MPI_SUM, reader);
    MPI_Allreduce(&local_skipped, &skipped, 1, MPI_UNSIGNED, MPI_SUM, reader);
    local &= integrity_count == static_cast<unsigned>(writer_size);
    if (writer_size == 4 && reader_size == 4) {
      if (skipped == 0U && world_rank == 0)
        std::cerr << "restore prefilter decoded every source on every rank\n";
      local &= skipped > 0U;
    }
    MPI_Comm_free(&reader);
  }
  MPI_Barrier(world);
  const int value = local ? 1 : 0;
  int passed = 0;
  MPI_Allreduce(&value, &passed, 1, MPI_INT, MPI_MIN, world);
  return passed != 0;
}

bool restart_syscalls(MPI_Comm communicator, const fs::path& directory) {
  int rank = 0, ranks = 0;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &ranks);
  Fixture fixture;
  bool passed = fixture.initialize(communicator);
  const std::array<RestartExpectedField, 4U> fields{{
      {RestartFieldRole::velocity, 0U, 3U},
      {RestartFieldRole::pressure_perturbation, 1U, 1U},
      {RestartFieldRole::enthalpy, 2U, 1U},
      {RestartFieldRole::independent_species, 3U, 1U}}};
  const RestartExpected expected{kGlobal, fixture.patch, kPlan, kSchema, kGeometry,
                                 {fields.data(), fields.size()}};
  for (int target : {0, ranks - 1}) {
    for (int mode = 0; mode < 11; ++mode) {
      if (mode == 10 && target != 0) continue;
      const auto root = directory / (std::to_string(target) + "-" + std::to_string(mode));
      passed &= static_cast<bool>(RestartWriter::write(communicator, root, fixture.snapshot(29U)));
      restart_syscall_probe::prefix = root.string();
      restart_syscall_probe::mode = mode;
      restart_syscall_probe::fired = false;
      restart_syscall_probe::close_fired = restart_syscall_probe::current_renamed = false;
      restart_syscall_probe::enabled = rank == target;
      RestartWriteReport report;
      const auto written = RestartWriter::write(communicator, root, fixture.snapshot(30U), {1U, &report});
      restart_syscall_probe::enabled = false;
      int okay = rank != target || restart_syscall_probe::fired;
      if (mode < 4) {
        okay &= written && !report.failure.valid &&
                 report.publication == RestartPublicationState::durable;
      } else {
        const auto operation = mode == 4 ? IoFailureOperation::open
            : (mode == 5 || mode == 8 || mode == 9) ? IoFailureOperation::write
            : mode == 7 ? IoFailureOperation::close : IoFailureOperation::sync;
        const int error = mode == 4 ? EACCES : (mode == 5 || mode == 8) ? ENOSPC : EIO;
        okay &= written.code == StatusCode::io_failure && report.failure.valid &&
                 report.failure.operation == operation && report.failure.system_error == error &&
                 report.failure.rank == target &&
                 report.publication == (mode == 10 ? RestartPublicationState::visible_not_durable
                                                    : RestartPublicationState::not_switched);
      }
      RestartImage restored;
      const auto read = RestartReader::load(communicator, root, expected, restored);
      okay &= read && restored.step == (mode < 4 || mode == 10 ? 30U : 29U) &&
               verify(restored, fixture.patch);
      MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, communicator);
      if (rank == 0)
        std::cout << "restart syscall target=" << target << " mode=" << mode
                  << " passed=" << okay << '\n';
      passed &= okay != 0;
    }
    for (int mode : {7, 11, 12, 13, 14}) {
      const auto root = directory / ("reader-" + std::to_string(target) + "-" + std::to_string(mode));
      passed &= static_cast<bool>(RestartWriter::write(communicator, root, fixture.snapshot(29U)));
      restart_syscall_probe::prefix = root.string();
      restart_syscall_probe::mode = mode;
      restart_syscall_probe::fired = restart_syscall_probe::close_fired = false;
      restart_syscall_probe::enabled = rank == target;
      RestartImage image;
      image.step = 999U;
      RestartReadReport report;
      const auto loaded = RestartReader::load(communicator, root, expected, image, &report);
      restart_syscall_probe::enabled = false;
      int okay = rank != target || restart_syscall_probe::fired;
      if (mode == 11 || mode == 13)
        okay &= loaded && !report.failure.valid && image.step == 29U && verify(image, fixture.patch);
      else
        okay &= loaded.code == StatusCode::io_failure && image.step == 999U &&
                 report.failure.valid && report.failure.rank == target &&
                 report.failure.system_error == EIO &&
                 report.failure.operation == (mode == 7 ? IoFailureOperation::close
                     : mode == 12 ? IoFailureOperation::read : IoFailureOperation::stat);
      MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, communicator);
      if (rank == 0) std::cout << "restart reader syscall target=" << target
                              << " mode=" << mode << " passed=" << okay << '\n';
      passed &= okay != 0;
    }
  }
  const auto budget_root = directory / "budget";
  RestartWriteReport budget;
  passed &= static_cast<bool>(RestartWriter::write(communicator, budget_root, fixture.snapshot(40U), {2U, &budget}));
  const auto exact_capacity = budget.peak_bulk_staging_bytes;
  passed &= exact_capacity >= budget.rank_payload_bytes && exact_capacity != 0U;
  const auto denied = RestartWriter::write(communicator, budget_root, fixture.snapshot(41U),
      {2U, &budget, exact_capacity - 1U});
  passed &= denied.code == StatusCode::allocation_failure &&
            budget.publication == RestartPublicationState::not_switched;
  RestartImage image;
  passed &= static_cast<bool>(RestartReader::load(communicator, budget_root, expected, image)) && image.step == 40U;
  passed &= static_cast<bool>(RestartWriter::write(communicator, budget_root, fixture.snapshot(41U),
      {2U, &budget, exact_capacity}));
  restart_syscall_probe::prefix = budget_root.string();
  restart_syscall_probe::mode = 15;
  restart_syscall_probe::fired = false;
  restart_syscall_probe::enabled = rank == 0;
  restart_allocation_probe::oversized = 0U;
  restart_allocation_probe::armed = true;
  const auto oversized_verify = RestartWriter::write(communicator, budget_root, fixture.snapshot(42U),
      {2U, &budget, exact_capacity});
  restart_allocation_probe::armed = restart_syscall_probe::enabled = false;
  passed &= oversized_verify.code == StatusCode::io_failure &&
      budget.publication == RestartPublicationState::not_switched &&
      restart_allocation_probe::oversized == 0U && (rank != 0 || restart_syscall_probe::fired);
  passed &= static_cast<bool>(RestartReader::load(communicator, budget_root, expected, image)) && image.step == 41U;
  if (rank == 0) std::cout << "writer oversized verification preflight=" << passed << '\n';
  return passed;
}

bool retention_order(MPI_Comm communicator, const fs::path& directory) {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  Fixture fixture;
  bool passed = fixture.initialize(communicator);
  for (std::uint32_t keep = 1U; keep <= 3U; ++keep) {
    const auto root = directory / std::to_string(keep);
    if (rank == 0) {
      fs::create_directories(root / "unrelated-pending");
      fs::create_directories(root / "generation-invalid-pending");
      fs::create_directory_symlink(root / "unrelated-pending", root / "generation-999-999-pending");
    }
    MPI_Barrier(communicator);
    for (std::uint64_t step : {8U, 9U, 10U, 11U, 98U, 99U, 100U, 101U})
      passed &= static_cast<bool>(RestartWriter::write(
          communicator, root, fixture.snapshot(step), {keep}));
    if (rank == 0) {
      std::vector<std::uint64_t> retained;
      for (const auto& entry : fs::directory_iterator(root)) {
        const auto name = entry.path().filename().string();
        if (entry.is_directory() && !entry.is_symlink() &&
            name.rfind("generation-", 0U) == 0U &&
            name.size() > 11U && name[11U] >= '0' && name[11U] <= '9')
          retained.push_back(std::stoull(name.substr(11U)));
      }
      std::sort(retained.begin(), retained.end());
      const std::vector<std::uint64_t> expected =
          keep == 1U ? std::vector<std::uint64_t>{101U}
          : keep == 2U ? std::vector<std::uint64_t>{100U, 101U}
                       : std::vector<std::uint64_t>{99U, 100U, 101U};
      if (retained != expected) {
        std::cerr << "retention keep=" << keep << " wrong steps:";
        for (auto step : retained) std::cerr << ' ' << step;
        std::cerr << '\n';
        passed = false;
      }
      passed &= fs::exists(root / "unrelated-pending") &&
                fs::exists(root / "generation-invalid-pending") &&
                fs::is_symlink(root / "generation-999-999-pending");
    }
    for (std::uint64_t step = 102U - keep; step <= 101U; ++step) {
      const auto inspect = directory / ("read-" + std::to_string(keep) + "-" + std::to_string(step));
      if (rank == 0) {
        fs::create_directory(inspect);
        for (const auto& entry : fs::directory_iterator(root)) {
          const auto name = entry.path().filename().string();
          if (name.rfind("generation-" + std::to_string(step) + "-", 0U) == 0U) {
            fs::create_directory_symlink(entry.path(), inspect / name);
            std::ofstream current(inspect / "current");
            current << name << '\n';
          }
        }
      }
      MPI_Barrier(communicator);
      const std::array<RestartExpectedField, 4U> catalog{{
          {RestartFieldRole::velocity,0U,3U}, {RestartFieldRole::pressure_perturbation,1U,1U},
          {RestartFieldRole::enthalpy,2U,1U}, {RestartFieldRole::independent_species,3U,1U}}};
      RestartExpected expected{kGlobal, fixture.patch, kPlan, kSchema, kGeometry,
                               {catalog.data(), catalog.size()}};
      RestartImage image;
      passed &= static_cast<bool>(RestartReader::load(communicator, inspect, expected, image)) &&
                image.step == step && verify(image, fixture.patch);
    }
  }
  return passed;
}

bool read_preallocation_bounds(MPI_Comm communicator, const fs::path& directory) {
  int rank = 0, ranks = 0;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &ranks);
  Fixture fixture;
  if (!fixture.initialize(communicator)) return false;
  const std::array<RestartExpectedField, 4U> fields{{
      {RestartFieldRole::velocity, 0U, 3U}, {RestartFieldRole::pressure_perturbation, 1U, 1U},
      {RestartFieldRole::enthalpy, 2U, 1U}, {RestartFieldRole::independent_species, 3U, 1U}}};
  const RestartExpected expected{kGlobal, fixture.patch, kPlan, kSchema, kGeometry,
                                  {fields.data(), fields.size()}};
  if (!RestartWriter::write(communicator, directory, fixture.snapshot(60U))) return false;
  RestartImage image;
  if (!RestartReader::load(communicator, directory, expected, image)) return false;
  bool passed = read_budget_boundaries(communicator, directory, expected, image);
  const double* original = image.fields[0U].values.data();
  std::string generation;
  { std::ifstream current(directory / "current"); std::getline(current, generation); }
  const auto file = directory / generation /
      ("rank-0000000" + std::to_string(ranks - 1) + ".bin");
  if (rank == 0) {
    const int fd = ::open(file.c_str(), O_WRONLY);
    if (fd < 0 || ::ftruncate(fd, static_cast<off_t>(UINT64_C(1) << 40U)) != 0)
      MPI_Abort(communicator, 2);
    ::close(fd);
  }
  MPI_Barrier(communicator);
  restart_allocation_probe::oversized = 0U;
  restart_allocation_probe::armed = true;
  RestartReadReport report;
  const Status status = RestartReader::load(communicator, directory, expected, image, &report);
  restart_allocation_probe::armed = false;
  passed &= status.code == StatusCode::io_failure &&
      restart_allocation_probe::oversized == 0U &&
      image.step == 60U && image.fields[0U].values.data() == original &&
      verify(image, fixture.patch) && fs::file_size(file) == (UINT64_C(1) << 40U);
  std::cerr << "read_preallocation rank=" << rank << " status=" << unsigned(status.code)
            << '/' << status.detail << " oversized_allocations="
            << restart_allocation_probe::oversized << " passed=" << passed << '\n';
  return passed;
}

bool read_malformed_metadata(MPI_Comm communicator, const fs::path& directory) {
  int rank = 0, ranks = 0;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &ranks);
  Fixture fixture;
  if (!fixture.initialize(communicator)) return false;
  const std::array<RestartExpectedField, 4U> fields{{
      {RestartFieldRole::velocity, 0U, 3U}, {RestartFieldRole::pressure_perturbation, 1U, 1U},
      {RestartFieldRole::enthalpy, 2U, 1U}, {RestartFieldRole::independent_species, 3U, 1U}}};
  const RestartExpected expected{kGlobal, fixture.patch, kPlan, kSchema, kGeometry,
                                  {fields.data(), fields.size()}};
  const auto original = directory / "original";
  if (!RestartWriter::write(communicator, original, fixture.snapshot(61U))) return false;
  RestartImage image;
  if (!RestartReader::load(communicator, original, expected, image)) return false;
  bool passed = true;
  std::string generation;
  { std::ifstream current(original / "current"); std::getline(current, generation); }
  const auto rank_file = "rank-0000000" + std::to_string(ranks - 1) + ".bin";
  for (int mode = 0; mode < 7; ++mode) {
    const auto root = directory / std::to_string(mode);
    if (rank == 0) {
      fs::copy(original, root, fs::copy_options::recursive);
      const auto manifest_path = root / generation / "manifest.bin";
      if (mode <= 1) {
        const auto file = mode == 0 ? root / "current" : manifest_path;
        const int fd = ::open(file.c_str(), O_WRONLY);
        if (fd < 0 || ::ftruncate(fd, static_cast<off_t>(UINT64_C(1) << 40U)) != 0)
          MPI_Abort(communicator, 2);
        ::close(fd);
      } else if (mode <= 3) {
        std::ifstream input(manifest_path, std::ios::binary);
        std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(input)), {});
        input.close();
        if (mode == 2) {
          // Legal integrity hash but an incorrect last-rank expected length.
          const std::size_t offset = bytes.size() - 8U - 40U + 24U;
          std::uint64_t size = 0U;
          for (unsigned n = 0; n < 8; ++n) size |= std::uint64_t(bytes[offset + n]) << (8U * n);
          ++size;
          for (unsigned n = 0; n < 8; ++n) bytes[offset + n] = (size >> (8U * n)) & 255U;
        } else {
          for (unsigned n = 12U; n < 16U; ++n) bytes[n] = 255U;
        }
        std::uint64_t hash = UINT64_C(1469598103934665603);
        for (std::size_t n = 0; n < bytes.size() - 8U; ++n) {
          hash ^= bytes[n]; hash *= UINT64_C(1099511628211);
        }
        if (hash == 0U) hash = 1U;
        for (unsigned n = 0; n < 8; ++n) bytes[bytes.size() - 8U + n] = (hash >> (8U * n)) & 255U;
        std::ofstream output(manifest_path, std::ios::binary | std::ios::trunc);
        output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
      } else {
        const auto file = mode == 6 ? root / "current" : root / generation / rank_file;
        fs::remove(file);
        if (mode == 4) {
          if (::mkfifo(file.c_str(), 0600) != 0) MPI_Abort(communicator, 2);
        } else if (mode == 5) {
          std::ofstream empty(file);
        } else fs::create_directory(file);
      }
    }
    MPI_Barrier(communicator);
    const auto* before = image.fields[0U].values.data();
    restart_allocation_probe::oversized = 0U;
    restart_allocation_probe::armed = true;
    const Status status = RestartReader::load(communicator, root, expected, image);
    restart_allocation_probe::armed = false;
    passed &= status.code == StatusCode::io_failure && restart_allocation_probe::oversized == 0U &&
        image.step == 61U && before == image.fields[0U].values.data() && verify(image, fixture.patch);
    if (rank == 0) std::cout << "read_malformed mode=" << mode << " status=" << unsigned(status.code)
        << '/' << status.detail << " passed=" << passed << '\n';
  }
  // Only independent copies were malformed; the original checkpoint is still readable.
  const Status intact = RestartReader::load(communicator, original, expected, image);
  return passed && intact && image.step == 61U && verify(image, fixture.patch);
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
    RestartWriteReport publication;
    status = RestartWriter::write(communicator, directory,
                                  fixture.snapshot(21U + index), {1U, &publication});
    detail::clear_restart_failure_for_test();
    passed &= status.code == StatusCode::io_failure;
    const auto expected_publication =
        points[index] == detail::RestartFailurePoint::after_current_switch
            ? RestartPublicationState::durable
            : RestartPublicationState::not_switched;
    passed &= publication.publication == expected_publication;
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
  passed &= read_preallocation_bounds(MPI_COMM_WORLD, base / "read-bounds");
  passed &= read_malformed_metadata(MPI_COMM_WORLD, base / "read-malformed");
  if (argc > 1 && std::strcmp(argv[1], "--read-bounds-only") == 0) {
    const int local = passed ? 1 : 0;
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) fs::remove_all(base);
    MPI_Finalize();
    return global ? 0 : 1;
  }
  passed &= spray_history_repartition(MPI_COMM_WORLD, 1, 4,
                                      base / "spray-one-to-four");
  passed &= spray_history_repartition(MPI_COMM_WORLD, 4, 1,
                                      base / "spray-four-to-one");
  passed &=
      record_repartition(MPI_COMM_WORLD, 1, 4, base / "records-one-to-four");
  passed &=
      record_repartition(MPI_COMM_WORLD, 4, 1, base / "records-four-to-one");
  passed &= record_repartition(MPI_COMM_WORLD, 1, 4,
                               base / "variable-one-to-four", true);
  passed &= record_repartition(MPI_COMM_WORLD, 4, 1,
                               base / "variable-four-to-one", true);
  passed &= record_repartition(MPI_COMM_WORLD, 4, 4,
                               base / "variable-four-to-four", true);
  passed &= record_repartition(MPI_COMM_WORLD, 4, 4, base / "variable-empty",
                               true, true);
  passed &= transition(MPI_COMM_WORLD, 1, 2, base / "one-to-two", 1U);
  passed &= transition(MPI_COMM_WORLD, 2, 4, base / "two-to-four", 2U);
  passed &= transition(MPI_COMM_WORLD, 4, 1, base / "four-to-one", 3U);
  passed &= transition(MPI_COMM_WORLD, 4, 4, base / "four-to-four", 4U);
  passed &= exact_transition(MPI_COMM_WORLD, base / "exact-four-to-four");
  passed &= exact_transition(MPI_COMM_WORLD, base / "signed-four-to-four", UINT64_C(0x391002));
  passed &= exact_transition(MPI_COMM_WORLD, base / "records-four-to-four",
                             UINT64_C(0x391003), true);
  passed &= failure_boundaries(MPI_COMM_WORLD, base / "failure-boundaries");
  passed &= retention_order(MPI_COMM_WORLD, base / "retention-order");
  passed &= restart_syscalls(MPI_COMM_WORLD, base / "syscalls");
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
