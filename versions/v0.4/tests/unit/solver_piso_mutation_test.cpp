// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/piso_fixture.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  if (void* pointer = std::malloc(size == 0U ? 1U : size)) {
    return pointer;
  }
  throw std::bad_alloc{};
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* pointer = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&pointer, alignment, requested) == 0 &&
      pointer != nullptr) {
    return pointer;
  }
  throw std::bad_alloc{};
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool test_hot_numeric_refresh_and_atomic_failure() {
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(8), "PISO hot fixture compiles cold");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField accepted = make_field(50U, cells, 1U, 0U, 6101U, 7101U);
  OwnedField previous = make_field(51U, cells, 1U, 0U, 6102U, 7102U);
  OwnedField drho_dp = make_field(52U, cells, 1U, 0U, 6103U, 7103U);
  OwnedField diagonal = make_field(53U, cells, 1U, 0U, 6104U, 7104U);
  OwnedField rhs = make_field(54U, cells, 1U, 0U, 6105U, 7105U);
  constexpr FieldId correction_field = 90U;
  OwnedField correction = make_field(correction_field, cells, 1U, 1U, 6106U,
                                     7106U);
  OwnedField applied = make_field(91U, cells, 1U, 0U, 6107U, 7107U);
  fill(accepted, 0.9);
  fill(previous, 0.7);
  fill(drho_dp, 0.02);
  fill(correction, 1.0);

  const BdfCoefficients bdf{10.0, -15.0, 5.0, 2U};
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 6201U);
  PisoIntermediateCertificate intermediate;
  passed &= expect(static_cast<bool>(fixture.coupler.refresh(
                       intermediate_input, intermediate)),
                   "PISO hot intermediate warms");
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(fixture.density.view);
  pressure_input.density_accepted = as_const(accepted.view);
  pressure_input.density_previous = as_const(previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = bdf;
  pressure_input.time = intermediate_input.momentum.time;
  pressure_input.geometry = intermediate_input.momentum.geometry;
  pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
  const PressureCorrectionSystemView system{diagonal.view, rhs.view};
  PressureCorrectionCertificate pressure;
  passed &= expect(static_cast<bool>(fixture.coupler.assemble_pressure_system(
                       pressure_input, system, pressure)),
                   "PISO pressure numeric warms");

  const std::array<HaloFieldSpec, 1U> operator_fields{{
      {correction_field, 1U, 1U}}};
  HaloEngine operator_halo;
  passed &= expect(static_cast<bool>(operator_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {operator_fields.data(), operator_fields.size()},
                       fixture.boundary.halo_topology())),
                   "pressure operator halo reserves cold");
  PressureLinearOperator pressure_operator;
  passed &= expect(static_cast<bool>(fixture.coupler.bind_pressure_operator(
                       {MPI_COMM_SELF, &operator_halo, 6202U,
                        correction_field},
                       system, pressure_operator)) &&
                       static_cast<bool>(pressure_operator.refresh(
                           {pressure,
                            {6203U, 6204U, 6205U, 6206U, 6207U},
                            6208U})) &&
                       static_cast<bool>(pressure_operator.apply(
                           correction.view, applied.view)),
                   "pressure operator hot path warms");
  if (!passed) {
    return false;
  }

  const std::uintptr_t workspace_address =
      fixture.coupler.workspace_storage_address();
  const std::uintptr_t coefficient_address =
      pressure_operator.coefficient_storage_address();
  Status hot_status;
  {
    allocation_observer::Guard guard;
    for (std::uint64_t repetition = 0U; repetition < 32U; ++repetition) {
      ++fixture.density.view.revision;
      ++intermediate_input.momentum.time;
      intermediate_input.predictor.time = intermediate_input.momentum.time;
      intermediate_input.pressure_reference.time =
          intermediate_input.momentum.time;
      intermediate_input.density = fixture.density.view;
      intermediate_input.predictor.predicted_density =
          fixture.density.view.revision;
      Status status = fixture.coupler.refresh(intermediate_input, intermediate);
      if (status) {
        pressure_input.intermediate = intermediate;
        pressure_input.pressure_reference =
            intermediate_input.pressure_reference;
        pressure_input.density_trial = as_const(fixture.density.view);
        pressure_input.time = intermediate_input.momentum.time;
        status = fixture.coupler.assemble_pressure_system(
            pressure_input, system, pressure);
      }
      if (status) {
        status = pressure_operator.refresh(
            {pressure,
             {6301U + repetition, 6302U + repetition, 6303U + repetition,
              6304U + repetition, 6305U + repetition},
             6306U + repetition});
      }
      if (status) {
        ++correction.view.revision;
        status = pressure_operator.apply(correction.view, applied.view);
      }
      if (!status && hot_status) {
        hot_status = status;
      }
    }
  }
  passed &= expect(static_cast<bool>(hot_status) &&
                       allocation_observer::count.load(
                           std::memory_order_relaxed) == 0U,
                   "PISO numeric refresh/operator hot path performs zero allocations");
  passed &= expect(fixture.coupler.workspace_storage_address() ==
                           workspace_address &&
                       pressure_operator.coefficient_storage_address() ==
                           coefficient_address,
                   "hot refresh preserves workspace and coefficient addresses");

  const PressureCorrectionCertificate marker = pressure;
  PressureCorrectionInput foreign_thermophysical = pressure_input;
  ++foreign_thermophysical.intermediate
        .thermophysical_boundary_rank_local_binding;
  fill(diagonal, -67.0);
  fill(rhs, -69.0);
  PressureCorrectionCertificate foreign_rejected = marker;
  const Status foreign_failure = fixture.coupler.assemble_pressure_system(
      foreign_thermophysical, system, foreign_rejected);
  passed &= expect(
      foreign_failure.code == StatusCode::invalid_plan &&
          diagonal.view.unchecked({1, 1, 1}, 0U) == -67.0 &&
          rhs.view.unchecked({1, 1, 1}, 0U) == -69.0 &&
          foreign_rejected.state == marker.state,
      "pressure assembly rejects a foreign thermophysical local binding "
      "atomically");
  const double saved = drho_dp.view.unchecked({1, 1, 1}, 0U);
  drho_dp.view.unchecked({1, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  fill(diagonal, -71.0);
  fill(rhs, -73.0);
  PressureCorrectionCertificate rejected = marker;
  const Status failure = fixture.coupler.assemble_pressure_system(
      pressure_input, system, rejected);
  passed &= expect(failure.code == StatusCode::numerical_failure &&
                       diagonal.view.unchecked({1, 1, 1}, 0U) == -71.0 &&
                       rhs.view.unchecked({1, 1, 1}, 0U) == -73.0 &&
                       rejected.state == marker.state,
                   "failed numeric refresh rolls back system and certificate");
  drho_dp.view.unchecked({1, 1, 1}, 0U) = saved;
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_hot_numeric_refresh_and_atomic_failure();
  MPI_Finalize();
  return passed ? 0 : 1;
}
