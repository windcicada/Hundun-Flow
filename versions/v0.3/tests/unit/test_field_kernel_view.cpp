// SPDX-License-Identifier: Apache-2.0

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "hundun/rt_error.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_kernel_field_view.hpp"
#include "tests/support/test_main.hpp"

namespace {

using hundun::runtime::AccessMode;
using hundun::runtime::ActorId;
using hundun::runtime::Error;
using hundun::runtime::FaceFieldView;
using hundun::runtime::FieldAccessPlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldLayoutSet;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FieldView;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::KernelCellView;
using hundun::runtime::KernelFaceView;
using hundun::runtime::OutputPolicy;
using hundun::runtime::PhaseId;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::with_kernel_cell_view;
using hundun::runtime::with_kernel_face_view;

constexpr PhaseId kPhase = 12U;
constexpr ActorId kActor = 34U;
constexpr const char *kDeadOwner = "field view owner is no longer alive";
constexpr const char *kStaleGeneration = "field view generation is stale";

template <class T>
constexpr ScalarType scalar_type();

template <>
constexpr ScalarType scalar_type<double>() {
  return ScalarType::float64;
}

template <>
constexpr ScalarType scalar_type<std::int32_t>() {
  return ScalarType::int32;
}

template <>
constexpr ScalarType scalar_type<std::uint8_t>() {
  return ScalarType::uint8;
}

FieldDescriptor descriptor(std::string name, FunctionSpace space,
                           ScalarType type, std::uint32_t components,
                           int ghost_width) {
  return FieldDescriptor{std::move(name),
                         "1",
                         "stage2_kernel_test",
                         space,
                         type,
                         components,
                         ghost_width,
                         true,
                         RestartPolicy::persistent,
                         OutputPolicy::selected};
}

template <class Function>
void expect_error_message(Function &&function, const char *expected) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()) == expected);
  }
  HUNDUN_CHECK(threw);
}

template <class T>
T sample_value(std::size_t linear) {
  if constexpr (std::is_same_v<T, double>) {
    return static_cast<double>(linear) + 0.25;
  } else if constexpr (std::is_same_v<T, std::int32_t>) {
    return static_cast<std::int32_t>(linear) - 97;
  } else {
    return static_cast<std::uint8_t>((linear * 17U + 3U) % 251U);
  }
}

template <class T>
struct ReturningCellCallback {
  int operator()(KernelCellView<T>) const { return 1; }
};

template <class T>
struct ReturningFaceCallback {
  int operator()(KernelFaceView<T>) const { return 1; }
};

template <class T>
constexpr void assert_kernel_representation_contracts() {
  using MutableCell = KernelCellView<T>;
  using ConstCell = KernelCellView<const T>;
  using MutableFace = KernelFaceView<T>;
  using ConstFace = KernelFaceView<const T>;

  static_assert(std::is_trivially_copyable_v<MutableCell>);
  static_assert(std::is_trivially_copyable_v<ConstCell>);
  static_assert(std::is_trivially_copyable_v<MutableFace>);
  static_assert(std::is_trivially_copyable_v<ConstFace>);
  static_assert(std::is_trivially_destructible_v<MutableCell>);
  static_assert(std::is_trivially_destructible_v<ConstCell>);
  static_assert(std::is_trivially_destructible_v<MutableFace>);
  static_assert(std::is_trivially_destructible_v<ConstFace>);
  static_assert(std::is_standard_layout_v<MutableCell>);
  static_assert(std::is_standard_layout_v<ConstCell>);
  static_assert(std::is_standard_layout_v<MutableFace>);
  static_assert(std::is_standard_layout_v<ConstFace>);
  static_assert(!std::is_default_constructible_v<MutableCell>);
  static_assert(!std::is_default_constructible_v<ConstCell>);
  static_assert(!std::is_default_constructible_v<MutableFace>);
  static_assert(!std::is_default_constructible_v<ConstFace>);

  static_assert(!std::is_constructible_v<
                MutableCell, std::byte *, Int3, int, std::uint32_t,
                std::size_t, std::size_t, std::size_t>);
  static_assert(!std::is_constructible_v<
                ConstCell, const std::byte *, Int3, int, std::uint32_t,
                std::size_t, std::size_t, std::size_t>);
  static_assert(!std::is_constructible_v<MutableCell, FieldView<T>>);
  static_assert(!std::is_constructible_v<ConstCell, FieldView<const T>>);
  static_assert(!std::is_constructible_v<
                MutableFace, std::byte *, std::size_t, std::uint32_t>);
  static_assert(!std::is_constructible_v<
                ConstFace, const std::byte *, std::size_t, std::uint32_t>);
  static_assert(!std::is_constructible_v<MutableFace, FaceFieldView<T>>);
  static_assert(!std::is_constructible_v<ConstFace,
                                         FaceFieldView<const T>>);

  static_assert(std::is_void_v<decltype(with_kernel_cell_view(
      std::declval<const FieldView<T> &>(), ReturningCellCallback<T>{}))>);
  static_assert(std::is_void_v<decltype(with_kernel_face_view(
      std::declval<const FaceFieldView<T> &>(), ReturningFaceCallback<T>{}))>);
}

void test_static_representation_and_public_construction_contracts() {
  assert_kernel_representation_contracts<double>();
  assert_kernel_representation_contracts<std::int32_t>();
  assert_kernel_representation_contracts<std::uint8_t>();
}

template <class T>
void test_cell_numerical_equivalence() {
  constexpr Int3 extent{3, 2, 2};
  constexpr int ghost = 2;
  constexpr std::uint32_t components = 3U;
  constexpr std::size_t padded_x = 7U;
  constexpr std::size_t padded_y = 6U;

  FieldRegistry registry;
  const FieldId field = registry.declare_field(descriptor(
      "cell", FunctionSpace::cell_average, scalar_type<T>(), components,
      ghost));
  registry.freeze();

  FieldAccessPlan plan(registry);
  plan.declare_access(kPhase, kActor, field, AccessMode::read_write);
  plan.freeze();
  FieldStorage storage(registry, FieldLayoutSet{extent, 0U});
  auto checked = storage.acquire_write<T>(plan, kPhase, kActor, field);

  int mutable_callbacks = 0;
  with_kernel_cell_view(checked, [&](auto kernel) {
    using Actual = decltype(kernel);
    static_assert(std::is_same_v<Actual, KernelCellView<T>>);
    ++mutable_callbacks;
    const auto kernel_extent = kernel.interior_extent();
    HUNDUN_CHECK(kernel_extent.x == extent.x);
    HUNDUN_CHECK(kernel_extent.y == extent.y);
    HUNDUN_CHECK(kernel_extent.z == extent.z);
    HUNDUN_CHECK(kernel.ghost_width() == ghost);
    HUNDUN_CHECK(kernel.components() == components);
    std::size_t linear = 0U;
    for (int k = -ghost; k < extent.z + ghost; ++k) {
      for (int j = -ghost; j < extent.y + ghost; ++j) {
        for (int i = -ghost; i < extent.x + ghost; ++i) {
          for (int component = 0;
               component < static_cast<int>(components);
               ++component, ++linear) {
            kernel(i, j, k, component) = sample_value<T>(linear);
            HUNDUN_CHECK(&kernel(i, j, k, component) ==
                         &checked(i, j, k, component));
          }
        }
      }
    }

    const auto base = reinterpret_cast<std::uintptr_t>(
        &kernel(-ghost, -ghost, -ghost, 0));
    HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(
                     &kernel(-ghost, -ghost, -ghost, 1)) ==
                 base + sizeof(T));
    HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(
                     &kernel(-ghost + 1, -ghost, -ghost, 0)) ==
                 base + components * sizeof(T));
    HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(
                     &kernel(-ghost, -ghost + 1, -ghost, 0)) ==
                 base + padded_x * components * sizeof(T));
    HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(
                     &kernel(-ghost, -ghost, -ghost + 1, 0)) ==
                 base + padded_x * padded_y * components * sizeof(T));
  });
  HUNDUN_CHECK(mutable_callbacks == 1);

  const FieldStorage &const_storage = storage;
  auto reader = const_storage.acquire_read<T>(plan, kPhase, kActor, field);
  int const_callbacks = 0;
  with_kernel_cell_view(reader, [&](auto kernel) {
    using Actual = decltype(kernel);
    static_assert(std::is_same_v<Actual, KernelCellView<const T>>);
    static_assert(std::is_same_v<decltype(kernel(0, 0, 0, 0)), const T &>);
    ++const_callbacks;
    std::size_t linear = 0U;
    for (int k = -ghost; k < extent.z + ghost; ++k) {
      for (int j = -ghost; j < extent.y + ghost; ++j) {
        for (int i = -ghost; i < extent.x + ghost; ++i) {
          for (int component = 0;
               component < static_cast<int>(components);
               ++component, ++linear) {
            HUNDUN_CHECK(kernel(i, j, k, component) ==
                         sample_value<T>(linear));
            HUNDUN_CHECK(&kernel(i, j, k, component) ==
                         &reader(i, j, k, component));
          }
        }
      }
    }
  });
  HUNDUN_CHECK(const_callbacks == 1);
}

template <class T>
void test_face_numerical_equivalence() {
  constexpr std::size_t face_count = 5U;
  constexpr std::uint32_t components = 3U;

  FieldRegistry registry;
  const FieldId field = registry.declare_field(descriptor(
      "face", FunctionSpace::face_value, scalar_type<T>(), components, 0));
  registry.freeze();

  FieldAccessPlan plan(registry);
  plan.declare_access(kPhase, kActor, field, AccessMode::read_write);
  plan.freeze();
  FieldStorage storage(registry,
                       FieldLayoutSet{Int3{1, 1, 1}, face_count});
  auto checked =
      storage.acquire_face_write<T>(plan, kPhase, kActor, field);

  int mutable_callbacks = 0;
  with_kernel_face_view(checked, [&](auto kernel) {
    using Actual = decltype(kernel);
    static_assert(std::is_same_v<Actual, KernelFaceView<T>>);
    ++mutable_callbacks;
    HUNDUN_CHECK(kernel.face_count() == face_count);
    HUNDUN_CHECK(kernel.components() == components);
    std::size_t linear = 0U;
    for (std::size_t face = 0U; face < face_count; ++face) {
      for (int component = 0; component < static_cast<int>(components);
           ++component, ++linear) {
        kernel(face, component) = sample_value<T>(linear);
        HUNDUN_CHECK(&kernel(face, component) == &checked(face, component));
      }
    }
    const auto base =
        reinterpret_cast<std::uintptr_t>(&kernel(std::size_t{0}, 0));
    HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(
                     &kernel(std::size_t{0}, 1)) ==
                 base + sizeof(T));
    HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(
                     &kernel(std::size_t{1}, 0)) ==
                 base + components * sizeof(T));
  });
  HUNDUN_CHECK(mutable_callbacks == 1);

  const FieldStorage &const_storage = storage;
  auto reader =
      const_storage.acquire_face_read<T>(plan, kPhase, kActor, field);
  int const_callbacks = 0;
  with_kernel_face_view(reader, [&](auto kernel) {
    using Actual = decltype(kernel);
    static_assert(std::is_same_v<Actual, KernelFaceView<const T>>);
    static_assert(std::is_same_v<decltype(kernel(std::size_t{0}, 0)),
                                 const T &>);
    ++const_callbacks;
    std::size_t linear = 0U;
    for (std::size_t face = 0U; face < face_count; ++face) {
      for (int component = 0; component < static_cast<int>(components);
           ++component, ++linear) {
        HUNDUN_CHECK(kernel(face, component) == sample_value<T>(linear));
        HUNDUN_CHECK(&kernel(face, component) == &reader(face, component));
      }
    }
  });
  HUNDUN_CHECK(const_callbacks == 1);
}

class CallbackFailure final : public std::runtime_error {
 public:
  CallbackFailure() : std::runtime_error("kernel callback failure") {}
};

void test_callback_count_return_and_exception_behavior() {
  FieldRegistry registry;
  const auto cell = registry.declare_field(descriptor(
      "cell", FunctionSpace::cell_average, ScalarType::float64, 1U, 0));
  const auto face = registry.declare_field(descriptor(
      "face", FunctionSpace::face_value, ScalarType::float64, 1U, 0));
  registry.freeze();

  FieldAccessPlan plan(registry);
  plan.declare_access(kPhase, kActor, cell, AccessMode::read_write);
  plan.declare_access(kPhase, kActor, face, AccessMode::read_write);
  plan.freeze();
  FieldStorage storage(registry, FieldLayoutSet{Int3{1, 1, 1}, 1U});
  auto cell_view =
      storage.acquire_write<double>(plan, kPhase, kActor, cell);
  auto face_view =
      storage.acquire_face_write<double>(plan, kPhase, kActor, face);

  int cell_calls = 0;
  with_kernel_cell_view(cell_view, [&](KernelCellView<double>) {
    ++cell_calls;
    return 73;
  });
  HUNDUN_CHECK(cell_calls == 1);

  int face_calls = 0;
  with_kernel_face_view(face_view, [&](KernelFaceView<double>) {
    ++face_calls;
    return 79;
  });
  HUNDUN_CHECK(face_calls == 1);

  bool cell_threw = false;
  try {
    with_kernel_cell_view(cell_view, [](KernelCellView<double>) -> void {
      throw CallbackFailure{};
    });
  } catch (const CallbackFailure &error) {
    cell_threw = std::string(error.what()) == "kernel callback failure";
  }
  HUNDUN_CHECK(cell_threw);

  bool face_threw = false;
  try {
    with_kernel_face_view(face_view, [](KernelFaceView<double>) -> void {
      throw CallbackFailure{};
    });
  } catch (const CallbackFailure &error) {
    face_threw = std::string(error.what()) == "kernel callback failure";
  }
  HUNDUN_CHECK(face_threw);
}

void test_generation_transition_rejects_before_callback() {
  FieldRegistry registry;
  const auto cell = registry.declare_field(descriptor(
      "cell", FunctionSpace::cell_average, ScalarType::float64, 1U, 0));
  const auto face = registry.declare_field(descriptor(
      "face", FunctionSpace::face_value, ScalarType::float64, 1U, 0));
  registry.freeze();
  FieldAccessPlan plan(registry);
  plan.declare_access(kPhase, kActor, cell, AccessMode::read_write);
  plan.declare_access(kPhase, kActor, face, AccessMode::read_write);
  plan.freeze();

  FieldStorage storage(registry, FieldLayoutSet{Int3{1, 1, 1}, 2U});
  auto old_cell =
      storage.acquire_write<double>(plan, kPhase, kActor, cell);
  auto old_face =
      storage.acquire_face_write<double>(plan, kPhase, kActor, face);
  old_cell(0, 0, 0, 0) = 83.0;
  old_face(1U, 0) = 89.0;
  storage.begin_rebuild();

  int stale_cell_calls = 0;
  expect_error_message(
      [&] {
        with_kernel_cell_view(old_cell, [&](auto) { ++stale_cell_calls; });
      },
      kStaleGeneration);
  HUNDUN_CHECK(stale_cell_calls == 0);

  int stale_face_calls = 0;
  expect_error_message(
      [&] {
        with_kernel_face_view(old_face, [&](auto) { ++stale_face_calls; });
      },
      kStaleGeneration);
  HUNDUN_CHECK(stale_face_calls == 0);

  auto fresh_cell =
      storage.acquire_write<double>(plan, kPhase, kActor, cell);
  auto fresh_face =
      storage.acquire_face_write<double>(plan, kPhase, kActor, face);
  with_kernel_cell_view(fresh_cell, [](auto kernel) {
    HUNDUN_CHECK_NEAR(kernel(0, 0, 0, 0), 83.0, 0.0);
  });
  with_kernel_face_view(fresh_face, [](auto kernel) {
    HUNDUN_CHECK_NEAR(kernel(1U, 0), 89.0, 0.0);
  });
}

void test_owner_destruction_rejects_before_callback() {
  FieldRegistry registry;
  const auto cell = registry.declare_field(descriptor(
      "cell", FunctionSpace::cell_average, ScalarType::float64, 1U, 0));
  const auto face = registry.declare_field(descriptor(
      "face", FunctionSpace::face_value, ScalarType::float64, 1U, 0));
  registry.freeze();
  FieldAccessPlan plan(registry);
  plan.declare_access(kPhase, kActor, cell, AccessMode::read_write);
  plan.declare_access(kPhase, kActor, face, AccessMode::read_write);
  plan.freeze();

  std::optional<FieldView<double>> dead_cell;
  std::optional<FaceFieldView<double>> dead_face;
  {
    FieldStorage storage(registry, FieldLayoutSet{Int3{1, 1, 1}, 1U});
    dead_cell.emplace(
        storage.acquire_write<double>(plan, kPhase, kActor, cell));
    dead_face.emplace(
        storage.acquire_face_write<double>(plan, kPhase, kActor, face));
  }

  int cell_calls = 0;
  expect_error_message(
      [&] {
        with_kernel_cell_view(*dead_cell, [&](auto) { ++cell_calls; });
      },
      kDeadOwner);
  HUNDUN_CHECK(cell_calls == 0);

  int face_calls = 0;
  expect_error_message(
      [&] {
        with_kernel_face_view(*dead_face, [&](auto) { ++face_calls; });
      },
      kDeadOwner);
  HUNDUN_CHECK(face_calls == 0);
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_static_representation_and_public_construction_contracts();
    test_cell_numerical_equivalence<double>();
    test_cell_numerical_equivalence<std::int32_t>();
    test_cell_numerical_equivalence<std::uint8_t>();
    test_face_numerical_equivalence<double>();
    test_face_numerical_equivalence<std::int32_t>();
    test_face_numerical_equivalence<std::uint8_t>();
    test_callback_count_return_and_exception_behavior();
    test_generation_transition_rejects_before_callback();
    test_owner_destruction_rejects_before_callback();
  });
}
