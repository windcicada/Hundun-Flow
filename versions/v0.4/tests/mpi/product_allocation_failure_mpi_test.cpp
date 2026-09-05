// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>
#include <sys/resource.h>
#include <unistd.h>

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>

#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"

namespace {
long fail_after = -1;
bool injected = false;
bool counting = false;
std::uint64_t allocations = 0U;
bool observing = false;
struct LiveAllocation {
  void* pointer{};
  std::size_t bytes{};
  bool aligned{};
};
std::array<LiveAllocation, 32768U> live{};
std::size_t live_count = 0U, live_bytes = 0U, peak_bytes = 0U;
std::array<MPI_Comm, 128U> owned_comms{};
std::array<MPI_Request, 4096U> owned_requests{};
std::size_t comm_count = 0U, request_count = 0U;
// Independent observation of the blocking collective primitives used by the
// exercised product paths. This is not a replacement for production counters.
std::array<std::uint64_t, 8U> collective_calls{};
void record(void* pointer, std::size_t bytes, bool aligned = false) {
  if (!observing) return;
  const auto hash =
      (reinterpret_cast<std::uintptr_t>(pointer) >> 4U) % live.size();
  for (std::size_t i = 0U; i < live.size(); ++i) {
    auto& slot = live[(hash + i) % live.size()];
    if (slot.pointer == nullptr) {
      slot = {pointer, bytes, aligned};
      ++live_count;
      live_bytes += bytes;
      if (live_bytes > peak_bytes) peak_bytes = live_bytes;
      return;
    }
  }
  std::abort();  // Observer overflow must never masquerade as a passing test.
}
void release(void* pointer) noexcept {
  if (pointer != nullptr && live_count != 0U) {
    const auto hash =
        (reinterpret_cast<std::uintptr_t>(pointer) >> 4U) % live.size();
    // Deletions leave holes; do not stop the search at an empty slot.
    for (std::size_t i = 0U; i < live.size(); ++i) {
      auto& slot = live[(hash + i) % live.size()];
      if (slot.pointer != pointer) continue;
      --live_count;
      live_bytes -= slot.bytes;
      slot = {};
      break;
    }
  }
  std::free(pointer);
}
template <class Handle, std::size_t N>
void acquired(Handle handle, std::array<Handle, N>& values,
              std::size_t& count) {
  if (!observing) return;
  if (count == N) std::abort();
  values[count++] = handle;
}
template <class Handle, std::size_t N>
void released(Handle handle, std::array<Handle, N>& values,
              std::size_t& count) {
  for (std::size_t i = 0U; i < count; ++i)
    if (values[i] == handle) {
      values[i] = values[--count];
      break;
    }
}
void checkpoint() {
  if (counting) ++allocations;
  if (fail_after >= 0 && fail_after-- == 0) {
    injected = true;
    throw std::bad_alloc{};
  }
}
}  // namespace
void* operator new(std::size_t size) {
  checkpoint();
  if (void* pointer = std::malloc(size ? size : 1U)) {
    record(pointer, size);
    return pointer;
  }
  throw std::bad_alloc{};
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return ::operator new(size);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return ::operator new(size);
  } catch (...) {
    return nullptr;
  }
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  checkpoint();
  void* pointer = nullptr;
  if (::posix_memalign(&pointer, static_cast<std::size_t>(alignment),
                       size ? size : 1U) == 0) {
    record(pointer, size, true);
    return pointer;
  }
  throw std::bad_alloc{};
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return ::operator new(size, alignment);
}
void operator delete(void* pointer) noexcept { release(pointer); }
void operator delete[](void* pointer) noexcept { release(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { release(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  release(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  release(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  release(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  release(pointer);
}
void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept {
  release(pointer);
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  release(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  release(pointer);
}

extern "C" int MPI_Comm_dup(MPI_Comm input, MPI_Comm* output) {
  if (observing) ++collective_calls[7U];
  const int result = PMPI_Comm_dup(input, output);
  if (result == MPI_SUCCESS) acquired(*output, owned_comms, comm_count);
  return result;
}
extern "C" int MPI_Comm_split_type(MPI_Comm input, int type, int key,
                                   MPI_Info info, MPI_Comm* output) {
  const int result = PMPI_Comm_split_type(input, type, key, info, output);
  if (result == MPI_SUCCESS && *output != MPI_COMM_NULL)
    acquired(*output, owned_comms, comm_count);
  return result;
}
extern "C" int MPI_Comm_free(MPI_Comm* comm) {
  const MPI_Comm previous = *comm;
  const int result = PMPI_Comm_free(comm);
  if (result == MPI_SUCCESS) released(previous, owned_comms, comm_count);
  return result;
}
extern "C" int MPI_Send_init(const void* buffer, int count,
                             MPI_Datatype datatype, int peer, int tag,
                             MPI_Comm comm, MPI_Request* request) {
  const int result =
      PMPI_Send_init(buffer, count, datatype, peer, tag, comm, request);
  if (result == MPI_SUCCESS) acquired(*request, owned_requests, request_count);
  return result;
}
extern "C" int MPI_Recv_init(void* buffer, int count, MPI_Datatype datatype,
                             int peer, int tag, MPI_Comm comm,
                             MPI_Request* request) {
  const int result =
      PMPI_Recv_init(buffer, count, datatype, peer, tag, comm, request);
  if (result == MPI_SUCCESS) acquired(*request, owned_requests, request_count);
  return result;
}
extern "C" int MPI_Request_free(MPI_Request* request) {
  const MPI_Request previous = *request;
  const int result = PMPI_Request_free(request);
  if (result == MPI_SUCCESS) released(previous, owned_requests, request_count);
  return result;
}

extern "C" int MPI_Allreduce(const void* send, void* receive, int count,
                             MPI_Datatype type, MPI_Op op, MPI_Comm comm) {
  if (observing) ++collective_calls[0U];
  return PMPI_Allreduce(send, receive, count, type, op, comm);
}
extern "C" int MPI_Bcast(void* buffer, int count, MPI_Datatype type, int root,
                         MPI_Comm comm) {
  if (observing) ++collective_calls[1U];
  return PMPI_Bcast(buffer, count, type, root, comm);
}
extern "C" int MPI_Allgather(const void* send, int send_count,
                             MPI_Datatype send_type, void* receive,
                             int receive_count, MPI_Datatype receive_type,
                             MPI_Comm comm) {
  if (observing) ++collective_calls[2U];
  return PMPI_Allgather(send, send_count, send_type, receive, receive_count,
                        receive_type, comm);
}
extern "C" int MPI_Allgatherv(const void* send, int count,
                              MPI_Datatype send_type, void* receive,
                              const int* counts, const int* offsets,
                              MPI_Datatype receive_type, MPI_Comm comm) {
  if (observing) ++collective_calls[3U];
  return PMPI_Allgatherv(send, count, send_type, receive, counts, offsets,
                         receive_type, comm);
}
extern "C" int MPI_Alltoall(const void* send, int send_count,
                            MPI_Datatype send_type, void* receive,
                            int receive_count, MPI_Datatype receive_type,
                            MPI_Comm comm) {
  if (observing) ++collective_calls[4U];
  return PMPI_Alltoall(send, send_count, send_type, receive, receive_count,
                       receive_type, comm);
}
extern "C" int MPI_Alltoallv(const void* send, const int* send_counts,
                             const int* send_offsets, MPI_Datatype send_type,
                             void* receive, const int* receive_counts,
                             const int* receive_offsets,
                             MPI_Datatype receive_type, MPI_Comm comm) {
  if (observing) ++collective_calls[5U];
  return PMPI_Alltoallv(send, send_counts, send_offsets, send_type, receive,
                        receive_counts, receive_offsets, receive_type, comm);
}
extern "C" int MPI_Barrier(MPI_Comm comm) {
  if (observing) ++collective_calls[6U];
  return PMPI_Barrier(comm);
}

namespace {
hundun::v04::Status prepare_restart_image(
    const hundun::v04::ValidatedModel& model,
    const std::filesystem::path& case_root, hundun::v04::RestartImage& image,
    std::filesystem::path* retained_root = nullptr) {
  using namespace hundun::v04;
  int rank = 0, id = static_cast<int>(getpid());
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Bcast(&id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  const auto root = std::filesystem::temp_directory_path() /
                    ("hundun-restart-allocation-seed-" + std::to_string(id));
  CompiledCasePlan plan;
  ProductDriver driver;
  RestartExpected expected;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, case_root, plan);
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  if (status) status = driver.restart_expected(expected);
  if (status) status = driver.initialize({});
  DriverStepReport step;
  if (status) status = driver.advance({1, 1, 1, 1, 1}, step);
  RestartSnapshot snapshot;
  if (status) status = driver.committed_restart_snapshot(snapshot);
  if (status) status = RestartWriter::write(MPI_COMM_WORLD, root, snapshot);
  if (status) status = RestartReader::load(MPI_COMM_WORLD, root, expected, image);
  MPI_Barrier(MPI_COMM_WORLD);
  if (retained_root != nullptr) *retained_root = root;
  if (rank == 0 && retained_root == nullptr) {
    std::error_code error;
    std::filesystem::remove_all(root, error);
  }
  return status;
}

// Requested, simultaneously live C++ storage: no views, allocator overhead,
// libc/MPI malloc, stacks, input model or observer table are counted twice.
bool memory_profile(bool immersed, int rank) {
  using namespace hundun::v04;
  ValidatedModel model =
      test::product_model(immersed ? Int3{16, 16, 16} : Int3{8, 8, 8});
  std::filesystem::path case_root;
  if (immersed) {
    model.mesh.lower = {-2.0, -2.0, -2.0};
    model.mesh.upper = {2.0, 2.0, 2.0};
    model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
    model.immersed_boundary =
        ImmersedBoundarySpec{"cylinder_ascii.stl", ImmersedFluidSide::outside};
    case_root = HUNDUN_V04_TEST_DATA_ROOT;
  }
  int root_pid = rank == 0 ? ::getpid() : 0;
  MPI_Bcast(&root_pid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  const auto root = std::filesystem::temp_directory_path() /
                    ("hundun-memory-profile-" + std::to_string(root_pid));
  const auto visit_path = root / "Visit";
  const auto restart_path = root / "Restart";
  const auto row = [&](const char* phase, Status status) {
    std::size_t aligned_bytes = 0U;
    for (const auto& value : live)
      if (value.pointer && value.aligned) aligned_bytes += value.bytes;
    long rss_pages = 0, virtual_pages = 0;
    if (FILE* file = std::fopen("/proc/self/statm", "r")) {
      if (std::fscanf(file, "%ld %ld", &virtual_pages, &rss_pages) != 2)
        rss_pages = 0;
      std::fclose(file);
    }
    struct rusage usage {};
    (void)::getrusage(RUSAGE_SELF, &usage);
    std::printf(
        "MEMORY {\"geometry\":\"%s\",\"rank\":%d,\"phase\":\"%s\","
        "\"status_code\":%u,\"status_detail\":%u,"
        "\"cpp_live_bytes\":%zu,\"cpp_phase_peak_bytes\":%zu,"
        "\"aligned_live_bytes\":%zu,\"cpp_live_objects\":%zu,"
        "\"sampled_rss_bytes\":%llu,\"historical_maxrss_bytes\":%llu}\n",
        immersed ? "immersed" : "cartesian", rank, phase,
        static_cast<unsigned>(status.code), status.detail, live_bytes,
        peak_bytes, aligned_bytes, live_count,
        static_cast<unsigned long long>(rss_pages) *
            static_cast<unsigned long long>(::sysconf(_SC_PAGESIZE)),
        static_cast<unsigned long long>(usage.ru_maxrss) * 1024ULL);
    return aligned_bytes;
  };
  const auto stage = [&](const char* name, auto&& operation) {
    peak_bytes = live_bytes;
    collective_calls = {};
    observing = true;
    const Status status = operation();
    observing = false;
    row(name, status);
    return status;
  };
  row("before_compile", {});
  bool passed = true;
  {
    CompiledCasePlan plan;
    ProductDriver driver;
    IoServicePlan services;
    RestartExpected expected;
    Status status = stage("compile", [&] {
      return ProductCompiler::compile(MPI_COMM_WORLD, model, case_root, plan);
    });
    if (status) {
      services = *plan.io_services();
      std::size_t arena_bytes = 0U;
      for (const auto& value : live)
        if (value.pointer && value.aligned) arena_bytes += value.bytes;
      const Int3 n = plan.summary().local_cells;
      const auto align8 = [](std::size_t value) {
        return (value + 7U) & ~std::size_t{7U};
      };
      // FaceFluxStorage owns a separate three-replica final/history bundle
      // and six workspace replicas. They are aligned, but NOT in the arena.
      const std::size_t face_replica_doubles = align8(n.x + 1U) * n.y * n.z +
                                               align8(n.x) * (n.y + 1U) * n.z +
                                               align8(n.x) * n.y * (n.z + 1U);
      const std::size_t face_bytes = 9U * face_replica_doubles * sizeof(double);
      const std::size_t arena_theory =
          plan.summary().arena_doubles * sizeof(double);
      passed &= arena_bytes == arena_theory + face_bytes;
      std::printf(
          "MEMORY_CAPACITY {\"geometry\":\"%s\",\"rank\":%d,"
          "\"arena_bytes\":%zu,\"face_flux_bytes\":%zu,"
          "\"observed_aligned_bytes\":%zu,\"matched\":%s}\n",
          immersed ? "immersed" : "cartesian", rank, arena_theory, face_bytes,
          arena_bytes, passed ? "true" : "false");
      status = stage("create", [&] {
        Status created =
            ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
        return created ? driver.restart_expected(expected) : created;
      });
    }
    if (status)
      status = stage("initialize", [&] { return driver.initialize({}); });
    DriverStepReport report;
    for (unsigned step_index = 0; step_index < 2U && status; ++step_index) {
      status = stage(step_index == 0U ? "advance_cold" : "advance_warm", [&] {
        return driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, report);
      });
      if (status) {
        std::uint64_t total = 0U;
        for (auto calls : collective_calls) total += calls;
        const auto& resources = report.resources;
        const std::uint64_t legacy = resources.reduction_collectives +
                                     resources.mg_blocking_collectives +
                                     resources.predictor_blocking_collectives;
        const std::uint64_t reported =
            legacy + resources.structured_control_collectives +
            resources.ibm_control_collectives;
        passed &= total >= reported &&
                  resources.structured_control_collectives > 0U &&
                  (!immersed || resources.ibm_control_collectives > 0U);
        std::printf(
            "COLLECTIVES {\"geometry\":\"%s\",\"rank\":%d,"
            "\"step\":%u,\"legacy_solver_partial\":%llu,\"halo_control\":%llu,"
            "\"donor_control\":%llu,"
            "\"pmpi_advance_total\":%llu,\"reported_solver_partial\":%llu,"
            "\"unattributed\":%llu,\"allreduce\":%llu,\"bcast\":%llu,"
            "\"allgather\":%llu,\"allgatherv\":%llu,\"comm_dup\":%llu}\n",
            immersed ? "immersed" : "cartesian", rank, step_index + 1U,
            static_cast<unsigned long long>(legacy),
            static_cast<unsigned long long>(
                resources.structured_control_collectives),
            static_cast<unsigned long long>(resources.ibm_control_collectives),
            static_cast<unsigned long long>(total),
            static_cast<unsigned long long>(reported),
            static_cast<unsigned long long>(total - reported),
            static_cast<unsigned long long>(collective_calls[0]),
            static_cast<unsigned long long>(collective_calls[1]),
            static_cast<unsigned long long>(collective_calls[2]),
            static_cast<unsigned long long>(collective_calls[3]),
            static_cast<unsigned long long>(collective_calls[7]));
      }
    }
    if (status)
      status = stage("visit", [&] {
        CommittedOutputSnapshot snapshot;
        Status local = driver.committed_output_snapshot(snapshot);
        return local ? VisitWriter::write(MPI_COMM_WORLD, visit_path, services,
                                          snapshot)
                     : local;
      });
    if (status)
      status = stage("restart_write", [&] {
        RestartSnapshot snapshot;
        Status local = driver.committed_restart_snapshot(snapshot);
        return local ? RestartWriter::write(MPI_COMM_WORLD, restart_path,
                                            snapshot)
                     : local;
      });
    if (status)
      status = stage("restart_read", [&] {
        RestartImage image;
        return RestartReader::load(MPI_COMM_WORLD, restart_path, expected,
                                   image);
      });
    passed &= status && report.accepted;
  }
  peak_bytes = live_bytes;
  row("released", {});
  passed &= live_count == 0U && live_bytes == 0U && comm_count == 0U &&
            request_count == 0U;
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::error_code error;
    std::filesystem::remove_all(root, error);
    passed &= !error;
  }
  return passed;
}
}  // namespace

int main(int argc, char** argv) {
  using namespace hundun::v04;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (argc > 1 && std::strcmp(argv[1], "--memory-profile") == 0) {
    int passed =
        memory_profile(argc > 2 && std::strcmp(argv[2], "immersed") == 0, rank)
            ? 1
            : 0;
    MPI_Allreduce(MPI_IN_PLACE, &passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  const long first = argc > 1 ? std::strtol(argv[1], nullptr, 10) : 0;
  const long requested_last =
      argc > 2 ? std::strtol(argv[2], nullptr, 10) : first;
  const bool create = argc > 3 && std::strcmp(argv[3], "create") == 0;
  const bool initialize = argc > 3 && std::strcmp(argv[3], "initialize") == 0;
  const bool restore = argc > 3 && std::strcmp(argv[3], "initialize_restart") == 0;
  const bool read_restart =
      argc > 3 && std::strcmp(argv[3], "read_restart") == 0;
  const bool immersed = argc > 4 && std::strcmp(argv[4], "immersed") == 0;
  const char* operation = read_restart ? "read_restart"
                          : restore    ? "initialize_restart"
                          : initialize ? "initialize"
                          : create     ? "create"
                                       : "compile";
  bool passed = true;
  const int target_stride = ranks > 1 ? ranks - 1 : 1;
  for (int target = 0; target < ranks; target += target_stride) {
    if (!passed) break;
    ValidatedModel model =
        test::product_model(immersed ? Int3{16, 16, 16} : Int3{8, 8, 8});
    std::filesystem::path case_root;
    if (immersed) {
      model.mesh.lower = {-2.0, -2.0, -2.0};
      model.mesh.upper = {2.0, 2.0, 2.0};
      model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
      model.immersed_boundary = ImmersedBoundarySpec{
          "cylinder_ascii.stl", ImmersedFluidSide::outside};
      case_root = HUNDUN_V04_TEST_DATA_ROOT;
    }
    DriverInitialState initial;
    initial.velocity = {0.1, 0.0, 0.0};
    RestartImage image;
    std::filesystem::path restart_root;
    if ((restore || read_restart) &&
        !prepare_restart_image(model, case_root, image,
                               read_restart ? &restart_root : nullptr))
      MPI_Abort(MPI_COMM_WORLD, 2);
    std::uint64_t count = 0U;
    if (requested_last < 0) {
      CompiledCasePlan baseline;
      ProductDriver baseline_driver;
      if ((create || initialize || restore || read_restart) &&
          !ProductCompiler::compile(MPI_COMM_WORLD, model, case_root, baseline))
        MPI_Abort(MPI_COMM_WORLD, 2);
      if ((initialize || restore || read_restart) &&
          !ProductDriver::create(MPI_COMM_WORLD, std::move(baseline),
                                 baseline_driver))
        MPI_Abort(MPI_COMM_WORLD, 2);
      RestartExpected expected;
      RestartImage loaded;
      if (read_restart && !baseline_driver.restart_expected(expected))
        MPI_Abort(MPI_COMM_WORLD, 2);
      allocations = 0U;
      counting = rank == target;
      const Status compiled =
          read_restart ? RestartReader::load(MPI_COMM_WORLD, restart_root,
                                             expected, loaded)
          : restore    ? baseline_driver.initialize_restart(image)
          : initialize ? baseline_driver.initialize(initial)
          : create ? ProductDriver::create(MPI_COMM_WORLD, std::move(baseline),
                                           baseline_driver)
                   : ProductCompiler::compile(MPI_COMM_WORLD, model, case_root,
                                              baseline);
      counting = false;
      count = allocations;
      MPI_Bcast(&count, 1, MPI_UINT64_T, target, MPI_COMM_WORLD);
      if (!compiled || (!initialize && count == 0U) || count > 1000000U) {
        if (rank == 0)
          std::cerr << "baseline failed status="
                    << static_cast<int>(compiled.code) << '/' << compiled.detail
                    << '\n';
        MPI_Abort(MPI_COMM_WORLD, 2);
        return 2;
      }
      if (rank == 0)
        std::cerr << operation << " allocation sites=" << count << '\n';
    }
    const long last =
        requested_last < 0 ? static_cast<long>(count) - 1 : requested_last;
    for (long index = first; index <= last && passed; ++index) {
      CompiledCasePlan plan;
      ProductDriver driver;
      if ((create || initialize || restore || read_restart) &&
          !ProductCompiler::compile(MPI_COMM_WORLD, model, case_root, plan))
        MPI_Abort(MPI_COMM_WORLD, 2);
      if ((initialize || restore || read_restart) &&
          !ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver))
        MPI_Abort(MPI_COMM_WORLD, 2);
      const PlanFingerprint original_plan = plan.fingerprint();
      RestartExpected expected;
      RestartImage loaded;
      loaded.step = 999U;
      if (read_restart && !driver.restart_expected(expected))
        MPI_Abort(MPI_COMM_WORLD, 2);
      if (rank == 0)
        std::cerr << operation << " allocation target=" << target
                  << " index=" << index << '\n';
      injected = false;
      observing = true;
      peak_bytes = 0U;
      fail_after = rank == target ? index : -1;
      const Status status =
          read_restart ? RestartReader::load(MPI_COMM_WORLD, restart_root,
                                             expected, loaded)
          : restore    ? driver.initialize_restart(image)
          : initialize ? driver.initialize(initial)
          : create
              ? ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver)
              : ProductCompiler::compile(MPI_COMM_WORLD, model, case_root,
                                         plan);
      fail_after = -1;
      observing = false;
      const int local_failure_rank = injected ? rank : ranks;
      int lowest_failure_rank = ranks;
      MPI_Allreduce(&local_failure_rank, &lowest_failure_rank, 1, MPI_INT,
                    MPI_MIN, MPI_COMM_WORLD);
      const std::uint64_t wire =
          (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
      std::uint64_t low = 0U, high = 0U;
      MPI_Allreduce(&wire, &low, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&wire, &high, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
      int okay = lowest_failure_rank == target && low == high &&
                         status.code == StatusCode::allocation_failure &&
                         plan.fingerprint() == original_plan &&
                         plan.summary().sealed == create &&
                         !driver.initialized()
                     ? 1
                     : 0;
      if (restore) {
        // Every recoverable allocation failure leaves this driver reusable;
        // retry the real image with injection disabled, without re-creating it.
        const Status retried = driver.initialize_restart(image);
        RestartSnapshot restored;
        const Status captured = retried ? driver.committed_restart_snapshot(restored)
                                        : retried;
        okay &= captured && restored.step == image.step &&
                restored.time == image.time && restored.dt == image.dt &&
                restored.controller_state == image.controller_state;
      }
      if (read_restart) {
        okay &= loaded.step == 999U && loaded.fields.empty();
        const Status retried =
            RestartReader::load(MPI_COMM_WORLD, restart_root, expected, loaded);
        okay &= retried && loaded.step == image.step &&
                loaded.time == image.time &&
                loaded.fields.size() == image.fields.size() &&
                loaded.final_mass_flux == image.final_mass_flux &&
                loaded.previous_mass_flux == image.previous_mass_flux;
        loaded = RestartImage{};
      }
      driver = ProductDriver{};
      plan = CompiledCasePlan{};
      okay &= live_count == 0U && live_bytes == 0U && comm_count == 0U &&
              request_count == 0U;
      MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if (rank == 0)
        std::cout << operation << " failure consistent=" << okay
                  << " lowest_failure_rank=" << lowest_failure_rank
                  << " status=" << static_cast<int>(status.code) << '/'
                  << status.detail << " live_cpp=" << live_count
                  << " live_comms=" << comm_count
                  << " live_requests=" << request_count << '\n';
      passed = okay != 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (read_restart && rank == 0) {
      std::error_code error;
      std::filesystem::remove_all(restart_root, error);
    }
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
