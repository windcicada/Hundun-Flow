// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>
#include <unistd.h>

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <new>
#include <string>

#include "hundun/v04_app.hpp"

namespace {
long fail_after = -1;
std::uint64_t allocations = 0U;
bool counting = false;
bool injected = false;
}  // namespace

void* operator new(std::size_t size) {
  if (counting) ++allocations;
  if (fail_after >= 0 && fail_after-- == 0) {
    injected = true;
    throw std::bad_alloc{};
  }
  if (void* memory = std::malloc(size == 0U ? 1U : size)) return memory;
  throw std::bad_alloc{};
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }

int main(int argc, char** argv) {
  using namespace hundun::v04;
  namespace fs = std::filesystem;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  std::array<char, 256U> root_text{};
  if (rank == 0) {
    const fs::path root =
        fs::temp_directory_path() /
        ("hundun-case-allocation-" + std::to_string(::getpid()));
    const fs::path data = HUNDUN_V04_TEST_DATA;
    fs::create_directories(root);
    std::ifstream input(data / "case_minimal_valid.json");
    std::string json{std::istreambuf_iterator<char>(input), {}};
    const std::string empty = "\"scalars\": []";
    const std::string scalar =
        R"json("scalars": [{"stable_name":"tracer","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0}])json";
    for (std::size_t position = 0U;
         (position = json.find(empty, position)) != std::string::npos;
         position += scalar.size())
      json.replace(position, empty.size(), scalar);
    const std::string catalog = "\"transported_scalars\": []";
    json.replace(
        json.find(catalog), catalog.size(),
        R"json("transported_scalars": [{"stable_name":"tracer","role":"passive_scalar","molecular_schmidt":1,"turbulent_schmidt":1}])json");
    std::ofstream(root / "case.json") << json;
    fs::copy_file(data / "thermophysics.d", root / "thermophysics.d",
                  fs::copy_options::overwrite_existing);
    const std::string path = root.string();
    path.copy(root_text.data(), root_text.size() - 1U);
  }
  MPI_Bcast(root_text.data(), static_cast<int>(root_text.size()), MPI_CHAR, 0,
            MPI_COMM_WORLD);
  const fs::path root{root_text.data()};
  bool passed = true;
  for (int target = 0; target < ranks && passed; ++target) {
    ValidatedModel baseline;
    allocations = 0U;
    counting = rank == target;
    const Status initial =
        CaseCompiler::load_and_compile(MPI_COMM_WORLD, root, baseline);
    counting = false;
    std::uint64_t count = allocations;
    MPI_Bcast(&count, 1, MPI_UINT64_T, target, MPI_COMM_WORLD);
    if (!initial || count == 0U || count > 4096U) {
      passed = false;
      break;
    }
    std::uint64_t tested = 0U;
    for (std::uint64_t index = 0U; index < count; ++index) {
      ValidatedModel output;
      output.fingerprint = 123U;
      injected = false;
      fail_after = rank == target ? static_cast<long>(index) : -1;
      const Status status =
          CaseCompiler::load_and_compile(MPI_COMM_WORLD, root, output);
      fail_after = -1;
      int fired = rank == target && injected ? 1 : 0;
      MPI_Bcast(&fired, 1, MPI_INT, target, MPI_COMM_WORLD);
      const std::uint64_t packed =
          (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
      std::uint64_t low = 0U, high = 0U;
      MPI_Allreduce(&packed, &low, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&packed, &high, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
      const bool local =
          low == high && (!fired || (!status && output.fingerprint == 123U));
      int ok = local ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if (!ok) {
        if (rank == 0)
          std::cerr << "FAIL case allocation target=" << target
                    << " index=" << index << " status=" << packed << '\n';
        passed = false;
        break;
      }
      tested += fired != 0;
    }
    if (rank == 0)
      std::cout << "case allocation target=" << target << " tested=" << tested
                << " sites=" << count << '\n';
  }
  if (rank == 0 && passed) {
    const fs::path run = root.string() + "-run";
    const fs::path source = HUNDUN_V04_TEST_DATA;
    allocations = 0U;
    counting = true;
    const Status initial =
        ApplicationService::validate_run_directories(root, run, source);
    counting = false;
    const std::uint64_t count = allocations;
    passed &= initial && count > 0U;
    for (std::uint64_t index = 0U; index < count && passed; ++index) {
      injected = false;
      fail_after = static_cast<long>(index);
      const Status status =
          ApplicationService::validate_run_directories(root, run, source);
      fail_after = -1;
      passed &= !injected || !status;
    }
    std::cout << "application path allocation sites=" << count << '\n';
  }
  int success = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  passed = success != 0;
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) fs::remove_all(root);
  MPI_Finalize();
  return passed ? 0 : 1;
}
