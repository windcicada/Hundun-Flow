// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh_uniform_structured.hpp"

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "hundun/rt_types.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

using hundun::mesh::UniformStructuredMesh;
using hundun::runtime::Box3;
using hundun::runtime::CollectiveStatus;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Error;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;
using hundun::runtime::collective_status;

constexpr Int3 kGlobalCells{4, 2, 2};
constexpr Real3 kOrigin{-1.0, 0.0, 2.0};
constexpr Real3 kLength{2.0, 1.0, 4.0};
constexpr std::array<bool, 3> kNonperiodic{false, false, false};
constexpr double kTolerance = 1.0e-12;

static_assert(noexcept(std::declval<const UniformStructuredMesh&>().spacing_m()));
static_assert(
    noexcept(std::declval<const UniformStructuredMesh&>().cell_volume_m3()));
static_assert(
    noexcept(std::declval<const UniformStructuredMesh&>().local_extent()));
static_assert(noexcept(
    std::declval<const UniformStructuredMesh&>().owned_global_box()));

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(Box3 lhs, Box3 rhs) {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

bool finite_near(double actual, double expected, double tolerance) noexcept {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      !std::isfinite(tolerance) || tolerance < 0.0) {
    return false;
  }
  const double difference = std::abs(actual - expected);
  return std::isfinite(difference) && difference <= tolerance;
}

class LocalCheckState final {
 public:
  void require(bool condition, const char* expression, const char* file,
               int line) noexcept {
    if (!condition) {
      record(file, line, expression);
    }
  }

  void near(double actual, double expected, double tolerance,
            const char* file, int line) noexcept {
    if (finite_near(actual, expected, tolerance)) {
      return;
    }
    std::array<char, 224> detail{};
    const int count = std::snprintf(
        detail.data(), detail.size(),
        "numerical check failed: actual %.17g, expected %.17g, tolerance %.17g",
        actual, expected, tolerance);
    record(file, line,
           count < 0 ? "numerical check failed" : detail.data());
  }

  void unexpected(const std::exception& error, const char* file,
                  int line) noexcept {
    std::array<char, 320> detail{};
    const int count = std::snprintf(detail.data(), detail.size(),
                                    "unexpected exception: %s", error.what());
    record(file, line,
           count < 0 ? "unexpected exception" : detail.data());
  }

  void unexpected(const char* file, int line) noexcept {
    record(file, line, "unexpected non-standard exception");
  }

  bool passed() const noexcept { return passed_; }

  std::string_view message() const noexcept {
    return passed_ ? std::string_view{} : std::string_view{message_.data()};
  }

 private:
  void record(const char* file, int line, const char* detail) noexcept {
    if (!passed_) {
      return;
    }
    const int count = std::snprintf(message_.data(), message_.size(),
                                    "%s:%d: %s", file, line, detail);
    if (count < 0) {
      static_cast<void>(std::snprintf(message_.data(), message_.size(),
                                      "mesh test failure"));
    }
    passed_ = false;
  }

  bool passed_{true};
  std::array<char, 512> message_{};
};

#define MESH_LOCAL_CHECK(state, expression)                              \
  (state).require(static_cast<bool>(expression), #expression, __FILE__, \
                  __LINE__)

#define MESH_LOCAL_NEAR(state, actual, expected, tolerance) \
  (state).near((actual), (expected), (tolerance), __FILE__, __LINE__)

void check_near(LocalCheckState& state, Real3 actual,
                Real3 expected) noexcept {
  MESH_LOCAL_NEAR(state, actual.x, expected.x, kTolerance);
  MESH_LOCAL_NEAR(state, actual.y, expected.y, kTolerance);
  MESH_LOCAL_NEAR(state, actual.z, expected.z, kTolerance);
}

void collective_checkpoint(const MpiContext& mpi,
                           const LocalCheckState& state) {
  const CollectiveStatus status =
      collective_status(mpi, state.passed(), state.message());
  if (!status.ok) {
    throw std::runtime_error("rank " + std::to_string(status.failing_rank) +
                             ": " + status.message);
  }
}

Int3 with_component(Int3 value, int axis, int component) {
  if (axis == 0) {
    value.x = component;
  } else if (axis == 1) {
    value.y = component;
  } else {
    value.z = component;
  }
  return value;
}

Real3 with_component(Real3 value, int axis, double component) {
  if (axis == 0) {
    value.x = component;
  } else if (axis == 1) {
    value.y = component;
  } else {
    value.z = component;
  }
  return value;
}

template <class Function>
void expect_error(LocalCheckState& state, Function&& function,
                  const char* file, int line) noexcept {
  bool threw = false;
  try {
    function();
  } catch (const Error& error) {
    threw = true;
    state.require(error.what() != nullptr && error.what()[0] != '\0',
                  "runtime error message is nonempty", file, line);
  } catch (const std::exception& error) {
    state.unexpected(error, file, line);
    return;
  } catch (...) {
    state.unexpected(file, line);
    return;
  }
  state.require(threw, "expected hundun::runtime::Error", file, line);
}

StructuredDecomposition make_decomposition(const MpiContext& mpi) {
  return StructuredDecomposition::create(
      mpi, kGlobalCells, kNonperiodic,
      DecompositionOptions{Int3{2, 1, 1}});
}

void check_volume_case(LocalCheckState& state,
                       const StructuredDecomposition& decomposition,
                       Real3 spacing, double expected_volume) noexcept {
  const Real3 length{spacing.x * static_cast<double>(kGlobalCells.x),
                     spacing.y * static_cast<double>(kGlobalCells.y),
                     spacing.z * static_cast<double>(kGlobalCells.z)};
  try {
    const UniformStructuredMesh mesh(kGlobalCells, Real3{}, length,
                                     decomposition);
    check_near(state, mesh.spacing_m(), spacing);
    MESH_LOCAL_NEAR(state, mesh.cell_volume_m3(), expected_volume, 0.0);
  } catch (const std::exception& error) {
    state.unexpected(error, __FILE__, __LINE__);
  } catch (...) {
    state.unexpected(__FILE__, __LINE__);
  }
}

enum class VolumeCaseOutcome {
  value,
  runtime_error,
  unexpected_exception,
  unexpected_nonstandard_exception,
};

struct VolumeCaseObservation final {
  VolumeCaseOutcome outcome;
  double volume;
  bool exact_input;
  bool nonempty_error;
};

VolumeCaseObservation observe_volume_case(
    const StructuredDecomposition& decomposition, Real3 spacing) noexcept {
  const Real3 length{spacing.x * static_cast<double>(kGlobalCells.x),
                     spacing.y * static_cast<double>(kGlobalCells.y),
                     spacing.z * static_cast<double>(kGlobalCells.z)};
  const bool exact_input =
      std::isfinite(length.x) && std::isfinite(length.y) &&
      std::isfinite(length.z) &&
      length.x / static_cast<double>(kGlobalCells.x) == spacing.x &&
      length.y / static_cast<double>(kGlobalCells.y) == spacing.y &&
      length.z / static_cast<double>(kGlobalCells.z) == spacing.z;
  try {
    const UniformStructuredMesh mesh(kGlobalCells, Real3{}, length,
                                     decomposition);
    return VolumeCaseObservation{VolumeCaseOutcome::value,
                                 mesh.cell_volume_m3(), exact_input, true};
  } catch (const Error& error) {
    return VolumeCaseObservation{
        VolumeCaseOutcome::runtime_error, 0.0, exact_input,
        error.what() != nullptr && error.what()[0] != '\0'};
  } catch (const std::exception&) {
    return VolumeCaseObservation{VolumeCaseOutcome::unexpected_exception,
                                 0.0, exact_input, false};
  } catch (...) {
    return VolumeCaseObservation{
        VolumeCaseOutcome::unexpected_nonstandard_exception, 0.0,
        exact_input, false};
  }
}

const char* volume_case_outcome_name(VolumeCaseOutcome outcome) noexcept {
  switch (outcome) {
    case VolumeCaseOutcome::value:
      return "value";
    case VolumeCaseOutcome::runtime_error:
      return "runtime-error";
    case VolumeCaseOutcome::unexpected_exception:
      return "unexpected-exception";
    case VolumeCaseOutcome::unexpected_nonstandard_exception:
      return "unexpected-nonstandard-exception";
  }
  return "unknown";
}

bool accepted_as(const VolumeCaseObservation& observation,
                 double expected) noexcept {
  return observation.exact_input &&
         observation.outcome == VolumeCaseOutcome::value &&
         observation.volume == expected;
}

bool rejected(const VolumeCaseObservation& observation) noexcept {
  return observation.exact_input &&
         observation.outcome == VolumeCaseOutcome::runtime_error &&
         observation.nonempty_error;
}

void test_numerical_oracle(const MpiContext& mpi) {
  LocalCheckState state;
  const double infinity = std::numeric_limits<double>::infinity();
  const double quiet_nan = std::numeric_limits<double>::quiet_NaN();
  MESH_LOCAL_CHECK(state, finite_near(1.0, 1.0, 0.0));
  MESH_LOCAL_CHECK(state, !finite_near(quiet_nan, 0.0, kTolerance));
  MESH_LOCAL_CHECK(state, !finite_near(infinity, 0.0, kTolerance));
  MESH_LOCAL_CHECK(state, !finite_near(-infinity, 0.0, kTolerance));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, quiet_nan, kTolerance));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, infinity, kTolerance));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, -infinity, kTolerance));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, 0.0, quiet_nan));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, 0.0, infinity));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, 0.0, -infinity));
  MESH_LOCAL_CHECK(state, !finite_near(0.0, 0.0, -kTolerance));
  collective_checkpoint(mpi, state);
}

void test_mixed_scale_volumes(const MpiContext& mpi) {
  LocalCheckState state;
  {
    auto decomposition = make_decomposition(mpi);
    const double large = std::ldexp(1.0, 700);
    const double small = std::ldexp(1.0, -700);
    check_volume_case(state, decomposition, Real3{large, large, small},
                      large);
    check_volume_case(state, decomposition, Real3{small, small, large},
                      small);
    check_volume_case(
        state, decomposition,
        Real3{small, small, std::ldexp(1.0, 330)},
        std::ldexp(1.0, -1070));
  }
  collective_checkpoint(mpi, state);
}

void test_range_boundary_volumes(const MpiContext& mpi) {
  LocalCheckState state;
  {
    auto decomposition = make_decomposition(mpi);
    const VolumeCaseObservation upper_finite = observe_volume_case(
        decomposition,
        Real3{0x1.128fedbab2a43p+341, 0x1.e05f53e74d73dp+341,
              0x1.fcd0fe4fe3b83p+340});
    const VolumeCaseObservation upper_overflow = observe_volume_case(
        decomposition,
        Real3{0x1.8ee2bfb0af1e3p+341, 0x1.78cea81633acfp+341,
              0x1.be7d8ae880d00p+340});

    const double lower_base = std::ldexp(1.0, -500);
    const double lower_last = std::ldexp(1.0, -75);
    const VolumeCaseObservation lower_tie = observe_volume_case(
        decomposition, Real3{lower_base, lower_base, lower_last});
    const VolumeCaseObservation lower_above = observe_volume_case(
        decomposition,
        Real3{lower_base, lower_base,
              std::nextafter(lower_last,
                             std::numeric_limits<double>::infinity())});
    const VolumeCaseObservation lower_below = observe_volume_case(
        decomposition,
        Real3{lower_base, lower_base, std::nextafter(lower_last, 0.0)});

    std::array<char, 448> detail{};
    const int count = std::snprintf(
        detail.data(), detail.size(),
        "range observations: upper-finite=%s(%a), upper-overflow=%s(%a); "
        "lower-tie=%s(%a), lower-above=%s(%a), lower-below=%s(%a)",
        volume_case_outcome_name(upper_finite.outcome), upper_finite.volume,
        volume_case_outcome_name(upper_overflow.outcome), upper_overflow.volume,
        volume_case_outcome_name(lower_tie.outcome), lower_tie.volume,
        volume_case_outcome_name(lower_above.outcome), lower_above.volume,
        volume_case_outcome_name(lower_below.outcome), lower_below.volume);
    const char* const diagnostic =
        count < 0 ? "range-boundary observation failed" : detail.data();

    state.require(accepted_as(upper_finite,
                              std::numeric_limits<double>::max()),
                  diagnostic, __FILE__, __LINE__);
    state.require(rejected(upper_overflow), diagnostic, __FILE__, __LINE__);
    state.require(rejected(lower_tie), diagnostic, __FILE__, __LINE__);
    state.require(
        accepted_as(lower_above,
                    std::numeric_limits<double>::denorm_min()),
        diagnostic, __FILE__, __LINE__);
    state.require(rejected(lower_below), diagnostic, __FILE__, __LINE__);
  }
  collective_checkpoint(mpi, state);
}

void test_geometry_and_ownership(const MpiContext& mpi) {
  LocalCheckState state;
  {
    auto decomposition = make_decomposition(mpi);
    const Int3 expected_extent = decomposition.local_extent();
    const Box3 expected_box = decomposition.owned_box();
    std::optional<UniformStructuredMesh> mesh;
    try {
      mesh.emplace(kGlobalCells, kOrigin, kLength, decomposition);
    } catch (const std::exception& error) {
      state.unexpected(error, __FILE__, __LINE__);
    } catch (...) {
      state.unexpected(__FILE__, __LINE__);
    }

    if (mesh.has_value()) {
      check_near(state, mesh->spacing_m(), Real3{0.5, 0.5, 2.0});
      MESH_LOCAL_NEAR(state, mesh->cell_volume_m3(), 0.5, kTolerance);
      MESH_LOCAL_CHECK(state, same(mesh->local_extent(), expected_extent));
      MESH_LOCAL_CHECK(state,
                       same(mesh->owned_global_box(), expected_box));
    }

    int local_first = 0;
    int local_last = 0;
    for (int k = 0; k < expected_extent.z; ++k) {
      for (int j = 0; j < expected_extent.y; ++j) {
        for (int i = 0; i < expected_extent.x; ++i) {
          const Int3 local{i, j, k};
          const Int3 global{expected_box.begin.x + i,
                            expected_box.begin.y + j,
                            expected_box.begin.z + k};
          const Real3 expected_center{
              kOrigin.x + (static_cast<double>(global.x) + 0.5) * 0.5,
              kOrigin.y + (static_cast<double>(global.y) + 0.5) * 0.5,
              kOrigin.z + (static_cast<double>(global.z) + 0.5) * 2.0};
          std::optional<Real3> actual;
          if (mesh.has_value()) {
            try {
              actual = mesh->cell_center(local);
            } catch (const std::exception& error) {
              state.unexpected(error, __FILE__, __LINE__);
            } catch (...) {
              state.unexpected(__FILE__, __LINE__);
            }
          }
          if (actual.has_value()) {
            check_near(state, *actual, expected_center);
          }

          if (same(global, Int3{0, 0, 0})) {
            if (actual.has_value()) {
              check_near(state, *actual, Real3{-0.75, 0.25, 3.0});
            }
            local_first = 1;
          }
          if (same(global, Int3{3, 1, 1})) {
            if (actual.has_value()) {
              check_near(state, *actual, Real3{0.75, 0.75, 5.0});
            }
            local_last = 1;
          }
        }
      }
    }

    const int local_nonzero_begin =
        expected_box.begin.x != 0 || expected_box.begin.y != 0 ||
                expected_box.begin.z != 0
            ? 1
            : 0;
    int global_first = 0;
    int global_last = 0;
    int ranks_with_nonzero_begin = 0;
    const int first_result = MPI_Allreduce(&local_first, &global_first, 1,
                                           MPI_INT, MPI_SUM, mpi.comm());
    const int last_result = MPI_Allreduce(&local_last, &global_last, 1,
                                          MPI_INT, MPI_SUM, mpi.comm());
    const int nonzero_result =
        MPI_Allreduce(&local_nonzero_begin, &ranks_with_nonzero_begin, 1,
                      MPI_INT, MPI_SUM, mpi.comm());
    MESH_LOCAL_CHECK(state, first_result == MPI_SUCCESS);
    MESH_LOCAL_CHECK(state, last_result == MPI_SUCCESS);
    MESH_LOCAL_CHECK(state, nonzero_result == MPI_SUCCESS);
    MESH_LOCAL_CHECK(state, global_first == 1);
    MESH_LOCAL_CHECK(state, global_last == 1);
    MESH_LOCAL_CHECK(state, ranks_with_nonzero_begin >= 1);

    if (mesh.has_value()) {
      for (int axis = 0; axis < 3; ++axis) {
        expect_error(
            state,
            [&] {
              static_cast<void>(mesh->cell_center(
                  with_component(Int3{0, 0, 0}, axis, -1)));
            },
            __FILE__, __LINE__);
        expect_error(
            state,
            [&] {
              static_cast<void>(mesh->cell_center(with_component(
                  Int3{0, 0, 0}, axis,
                  axis == 0 ? expected_extent.x
                            : (axis == 1 ? expected_extent.y
                                         : expected_extent.z))));
            },
            __FILE__, __LINE__);
      }
    }
  }
  collective_checkpoint(mpi, state);
}

void test_constructor_validation(const MpiContext& mpi) {
  LocalCheckState state;
  {
    auto decomposition = make_decomposition(mpi);

    for (int axis = 0; axis < 3; ++axis) {
      for (int invalid : {0, -1}) {
        expect_error(
            state,
            [&] {
              static_cast<void>(UniformStructuredMesh(
                  with_component(kGlobalCells, axis, invalid), kOrigin,
                  kLength, decomposition));
            },
            __FILE__, __LINE__);
      }
    }

    expect_error(
        state,
        [&] {
          static_cast<void>(UniformStructuredMesh(
              Int3{5, 2, 2}, kOrigin, kLength, decomposition));
        },
        __FILE__, __LINE__);

    for (int axis = 0; axis < 3; ++axis) {
      for (double invalid : {0.0, -1.0}) {
        expect_error(
            state,
            [&] {
              static_cast<void>(UniformStructuredMesh(
                  kGlobalCells, kOrigin,
                  with_component(kLength, axis, invalid), decomposition));
            },
            __FILE__, __LINE__);
      }
    }

    const double infinity = std::numeric_limits<double>::infinity();
    const double quiet_nan = std::numeric_limits<double>::quiet_NaN();
    for (int axis = 0; axis < 3; ++axis) {
      for (double invalid : {infinity, quiet_nan}) {
        expect_error(
            state,
            [&] {
              static_cast<void>(UniformStructuredMesh(
                  kGlobalCells, with_component(kOrigin, axis, invalid),
                  kLength, decomposition));
            },
            __FILE__, __LINE__);
        expect_error(
            state,
            [&] {
              static_cast<void>(UniformStructuredMesh(
                  kGlobalCells, kOrigin,
                  with_component(kLength, axis, invalid), decomposition));
            },
            __FILE__, __LINE__);
      }
    }

    expect_error(
        state,
        [&] {
          static_cast<void>(UniformStructuredMesh(
              kGlobalCells, kOrigin,
              Real3{std::numeric_limits<double>::denorm_min(), 1.0, 1.0},
              decomposition));
        },
        __FILE__, __LINE__);
    expect_error(
        state,
        [&] {
          static_cast<void>(UniformStructuredMesh(
              kGlobalCells, kOrigin, Real3{1.0e-200, 1.0e-200, 1.0e-200},
              decomposition));
        },
        __FILE__, __LINE__);
    expect_error(
        state,
        [&] {
          static_cast<void>(UniformStructuredMesh(
              kGlobalCells, kOrigin, Real3{1.0e200, 1.0e200, 1.0e200},
              decomposition));
        },
        __FILE__, __LINE__);
    expect_error(
        state,
        [&] {
          const double maximum = std::numeric_limits<double>::max();
          static_cast<void>(UniformStructuredMesh(
              kGlobalCells, Real3{maximum, 0.0, 0.0},
              Real3{maximum, 1.0e-100, 1.0e-100}, decomposition));
        },
        __FILE__, __LINE__);
  }
  collective_checkpoint(mpi, state);
}

void test_value_state_lifetime(const MpiContext& mpi) {
  LocalCheckState state;
  Int3 expected_extent{};
  Box3 expected_box{};
  std::optional<UniformStructuredMesh> mesh;
  {
    auto decomposition = make_decomposition(mpi);
    expected_extent = decomposition.local_extent();
    expected_box = decomposition.owned_box();
    try {
      mesh.emplace(kGlobalCells, kOrigin, kLength, decomposition);
    } catch (const std::exception& error) {
      state.unexpected(error, __FILE__, __LINE__);
    } catch (...) {
      state.unexpected(__FILE__, __LINE__);
    }
  }

  if (mesh.has_value()) {
    check_near(state, mesh->spacing_m(), Real3{0.5, 0.5, 2.0});
    MESH_LOCAL_NEAR(state, mesh->cell_volume_m3(), 0.5, kTolerance);
    MESH_LOCAL_CHECK(state, same(mesh->local_extent(), expected_extent));
    MESH_LOCAL_CHECK(state, same(mesh->owned_global_box(), expected_box));

    std::optional<Real3> actual;
    try {
      actual = mesh->cell_center(Int3{0, 0, 0});
    } catch (const std::exception& error) {
      state.unexpected(error, __FILE__, __LINE__);
    } catch (...) {
      state.unexpected(__FILE__, __LINE__);
    }
    if (actual.has_value()) {
      check_near(
          state, *actual,
          Real3{kOrigin.x +
                    (static_cast<double>(expected_box.begin.x) + 0.5) * 0.5,
                kOrigin.y +
                    (static_cast<double>(expected_box.begin.y) + 0.5) * 0.5,
                kOrigin.z +
                    (static_cast<double>(expected_box.begin.z) + 0.5) * 2.0});
    }
  }

  collective_checkpoint(mpi, state);
}

void run_tests(const MpiContext& mpi) {
  LocalCheckState state;
  MESH_LOCAL_CHECK(state, mpi.size() == 2);
  collective_checkpoint(mpi, state);
  test_numerical_oracle(mpi);
  test_geometry_and_ownership(mpi);
  test_mixed_scale_volumes(mpi);
  test_range_boundary_volumes(mpi);
  test_constructor_validation(mpi);
  test_value_state_lifetime(mpi);
}

#undef MESH_LOCAL_NEAR
#undef MESH_LOCAL_CHECK

}  // namespace

int main(int argc, char** argv) {
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    result = hundun::test::run([&] { run_tests(mpi); });
  }
  return result;
}
