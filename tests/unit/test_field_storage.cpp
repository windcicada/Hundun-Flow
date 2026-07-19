// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "runtime/src/field_epoch_test_access.hpp"
#include "tests/support/test_main.hpp"

namespace {

using hundun::runtime::Error;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FieldView;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::OutputPolicy;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::detail::FieldEpochTestAccess;

constexpr const char *kDeadOwnerMessage =
    "field view owner is no longer alive";
constexpr const char *kStaleGenerationMessage =
    "field view generation is stale";
constexpr const char *kGenerationWrapMessage =
    "field storage generation would wrap";

FieldDescriptor descriptor(std::string name = "passive_scalar",
                           ScalarType scalar_type = ScalarType::float64,
                           std::uint32_t components = 1, int ghost_width = 2,
                           FunctionSpace space = FunctionSpace::cell_average) {
  return FieldDescriptor{std::move(name),
                         "1",
                         "passive_scalar_solver",
                         space,
                         scalar_type,
                         components,
                         ghost_width,
                         true,
                         RestartPolicy::persistent,
                         OutputPolicy::selected};
}

template <class Function>
void expect_error(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()).empty() == false);
  }
  HUNDUN_CHECK(threw);
}

template <class Function>
void expect_error_message(Function &&function, const char *expected_message) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()) == expected_message);
  }
  HUNDUN_CHECK(threw);
}

FieldRegistry make_epoch_registry() {
  FieldRegistry registry;
  registry.declare_field(
      descriptor("epoch_value", ScalarType::float64, 1, 1));
  registry.freeze();
  return registry;
}

void test_registry_validation_and_stable_ids() {
  FieldRegistry registry;
  HUNDUN_CHECK(registry.frozen() == false);
  HUNDUN_CHECK(registry.size() == 0U);

  auto empty_name = descriptor("");
  expect_error([&] { registry.declare_field(std::move(empty_name)); });
  auto missing_unit = descriptor();
  missing_unit.unit.clear();
  expect_error([&] { registry.declare_field(std::move(missing_unit)); });
  auto missing_owner = descriptor();
  missing_owner.owner.clear();
  expect_error([&] { registry.declare_field(std::move(missing_owner)); });
  expect_error([&] {
    registry.declare_field(descriptor("zero", ScalarType::float64, 0));
  });
  expect_error([&] {
    registry.declare_field(
        descriptor("negative_ghost", ScalarType::float64, 1, -1));
  });

  for (const auto space :
       {static_cast<FunctionSpace>(-1), static_cast<FunctionSpace>(6)}) {
    auto invalid = descriptor("invalid_space");
    invalid.space = space;
    expect_error([&] { registry.declare_field(std::move(invalid)); });
  }
  for (const auto scalar :
       {static_cast<ScalarType>(-1), static_cast<ScalarType>(3)}) {
    auto invalid = descriptor("invalid_scalar");
    invalid.scalar_type = scalar;
    expect_error([&] { registry.declare_field(std::move(invalid)); });
  }
  for (const auto policy :
       {static_cast<RestartPolicy>(-1), static_cast<RestartPolicy>(2)}) {
    auto invalid = descriptor("invalid_restart");
    invalid.restart = policy;
    expect_error([&] { registry.declare_field(std::move(invalid)); });
  }
  for (const auto policy :
       {static_cast<OutputPolicy>(-1), static_cast<OutputPolicy>(3)}) {
    auto invalid = descriptor("invalid_output");
    invalid.output = policy;
    expect_error([&] { registry.declare_field(std::move(invalid)); });
  }
  HUNDUN_CHECK(registry.size() == 0U);

  const FieldId scalar_id = registry.declare_field(descriptor());
  const FieldId marker_id =
      registry.declare_field(descriptor("marker", ScalarType::uint8, 1, 0));
  HUNDUN_CHECK(scalar_id == 0U);
  HUNDUN_CHECK(marker_id == 1U);
  HUNDUN_CHECK(registry.size() == 2U);
  HUNDUN_CHECK(registry.field_id("passive_scalar") == scalar_id);
  HUNDUN_CHECK(registry.field_id("marker") == marker_id);
  HUNDUN_CHECK(registry.descriptor(scalar_id).name == "passive_scalar");
  HUNDUN_CHECK(registry.descriptor(marker_id).scalar_type == ScalarType::uint8);

  expect_error([&] { registry.declare_field(descriptor()); });
  HUNDUN_CHECK(registry.size() == 2U);
  expect_error([&] { static_cast<void>(registry.field_id("unknown")); });
  expect_error([&] {
    static_cast<void>(registry.descriptor(std::numeric_limits<FieldId>::max()));
  });

  registry.freeze();
  registry.freeze();
  HUNDUN_CHECK(registry.frozen());
  expect_error(
      [&] { registry.declare_field(descriptor("declared_too_late")); });
  HUNDUN_CHECK(registry.size() == 2U);
}

void test_storage_preconditions_and_function_spaces() {
  FieldRegistry unfrozen;
  unfrozen.declare_field(descriptor());
  expect_error([&] { FieldStorage storage(unfrozen, Int3{1, 1, 1}); });

  FieldRegistry registry;
  registry.declare_field(descriptor());
  registry.freeze();
  for (const auto extent : {Int3{0, 1, 1}, Int3{-1, 1, 1}, Int3{1, 0, 1},
                            Int3{1, -1, 1}, Int3{1, 1, 0}, Int3{1, 1, -1}}) {
    expect_error([&] { FieldStorage storage(registry, extent); });
  }

  for (const auto space :
       {FunctionSpace::face_value, FunctionSpace::vertex_value,
        FunctionSpace::element_dof, FunctionSpace::quadrature_point,
        FunctionSpace::particle}) {
    FieldRegistry unsupported;
    unsupported.declare_field(
        descriptor("unsupported", ScalarType::float64, 1, 0, space));
    unsupported.freeze();
    expect_error([&] { FieldStorage storage(unsupported, Int3{1, 1, 1}); });
  }
}

void test_typed_views_layout_bounds_and_constness() {
  FieldRegistry registry;
  const FieldId q_id =
      registry.declare_field(descriptor("q", ScalarType::float64, 2, 1));
  const FieldId index_id =
      registry.declare_field(descriptor("index", ScalarType::int32, 1, 2));
  const FieldId mask_id =
      registry.declare_field(descriptor("mask", ScalarType::uint8, 1, 0));
  registry.freeze();

  FieldStorage storage(registry, Int3{2, 2, 2});
  HUNDUN_CHECK(storage.interior_extent().x == 2);
  HUNDUN_CHECK(storage.interior_extent().y == 2);
  HUNDUN_CHECK(storage.interior_extent().z == 2);

  auto q = storage.view<double>(q_id);
  static_assert(std::is_same_v<decltype(q(-1, -1, -1, 0)), double &>,
                "mutable storage must return mutable references");
  HUNDUN_CHECK(q.interior_extent().x == 2);
  HUNDUN_CHECK(q.interior_extent().y == 2);
  HUNDUN_CHECK(q.interior_extent().z == 2);
  HUNDUN_CHECK(q.ghost_width() == 1);
  HUNDUN_CHECK(q.components() == 2U);

  q(-1, -1, -1, 1) = 11.0;
  q(0, -1, -1, 0) = 12.0;
  q(-1, 0, -1, 0) = 13.0;
  q(-1, -1, 0, 0) = 14.0;
  q(1, 1, 1, 1) = 16.0;
  q(2, 2, 2, 1) = 18.0;

  struct CornerSample {
    Int3 coordinate;
    double value;
  };
  constexpr std::array<CornerSample, 8> ghost_corners{{
      {{-1, -1, -1}, 101.0},
      {{2, -1, -1}, 102.0},
      {{-1, 2, -1}, 103.0},
      {{2, 2, -1}, 104.0},
      {{-1, -1, 2}, 105.0},
      {{2, -1, 2}, 106.0},
      {{-1, 2, 2}, 107.0},
      {{2, 2, 2}, 108.0},
  }};
  constexpr std::array<CornerSample, 8> interior_corners{{
      {{0, 0, 0}, 201.0},
      {{1, 0, 0}, 202.0},
      {{0, 1, 0}, 203.0},
      {{1, 1, 0}, 204.0},
      {{0, 0, 1}, 205.0},
      {{1, 0, 1}, 206.0},
      {{0, 1, 1}, 207.0},
      {{1, 1, 1}, 208.0},
  }};
  for (const auto &sample : ghost_corners) {
    q(sample.coordinate.x, sample.coordinate.y, sample.coordinate.z, 0) =
        sample.value;
  }
  for (const auto &sample : interior_corners) {
    q(sample.coordinate.x, sample.coordinate.y, sample.coordinate.z, 0) =
        sample.value;
  }

  HUNDUN_CHECK_NEAR(q(-1, -1, -1, 1), 11.0, 0.0);
  HUNDUN_CHECK_NEAR(q(0, -1, -1, 0), 12.0, 0.0);
  HUNDUN_CHECK_NEAR(q(-1, 0, -1, 0), 13.0, 0.0);
  HUNDUN_CHECK_NEAR(q(-1, -1, 0, 0), 14.0, 0.0);
  HUNDUN_CHECK_NEAR(q(1, 1, 1, 1), 16.0, 0.0);
  HUNDUN_CHECK_NEAR(q(2, 2, 2, 1), 18.0, 0.0);
  for (const auto &sample : ghost_corners) {
    HUNDUN_CHECK_NEAR(
        q(sample.coordinate.x, sample.coordinate.y, sample.coordinate.z, 0),
        sample.value, 0.0);
  }
  for (const auto &sample : interior_corners) {
    HUNDUN_CHECK_NEAR(
        q(sample.coordinate.x, sample.coordinate.y, sample.coordinate.z, 0),
        sample.value, 0.0);
  }

  const auto base = reinterpret_cast<std::uintptr_t>(&q(-1, -1, -1, 0));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&q(-1, -1, -1, 1)) ==
               base + sizeof(double));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&q(0, -1, -1, 0)) ==
               base + 2U * sizeof(double));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&q(-1, 0, -1, 0)) ==
               base + 8U * sizeof(double));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&q(-1, -1, 0, 0)) ==
               base + 32U * sizeof(double));
  HUNDUN_CHECK(base % alignof(double) == 0U);

  auto index = storage.view<std::int32_t>(index_id);
  HUNDUN_CHECK(index.ghost_width() == 2);
  index(-2, -2, -2, 0) = -123;
  index(3, 3, 3, 0) = 456;
  HUNDUN_CHECK(index(-2, -2, -2, 0) == -123);
  HUNDUN_CHECK(index(3, 3, 3, 0) == 456);
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&index(-2, -2, -2, 0)) %
                   alignof(std::int32_t) ==
               0U);

  auto mask = storage.view<std::uint8_t>(mask_id);
  mask(0, 0, 0, 0) = static_cast<std::uint8_t>(7);
  mask(1, 1, 1, 0) = static_cast<std::uint8_t>(251);
  HUNDUN_CHECK(mask(0, 0, 0, 0) == static_cast<std::uint8_t>(7));
  HUNDUN_CHECK(mask(1, 1, 1, 0) == static_cast<std::uint8_t>(251));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&mask(0, 0, 0, 0)) %
                   alignof(std::uint8_t) ==
               0U);

  const FieldStorage &const_storage = storage;
  auto const_q = const_storage.view<double>(q_id);
  static_assert(
      std::is_same_v<decltype(const_q(-1, -1, -1, 0)), const double &>,
      "const storage must return read-only references");
  HUNDUN_CHECK_NEAR(const_q(2, 2, 2, 1), 18.0, 0.0);
  HUNDUN_CHECK_NEAR(const_q(2, -1, 2, 0), 106.0, 0.0);

  expect_error([&] { static_cast<void>(storage.view<std::int32_t>(q_id)); });
  expect_error([&] { static_cast<void>(storage.view<double>(index_id)); });
  expect_error([&] { static_cast<void>(storage.view<double>(mask_id)); });
  expect_error([&] {
    static_cast<void>(
        storage.view<double>(std::numeric_limits<FieldId>::max()));
  });

  for (const auto coordinate : {Int3{-2, 0, 0}, Int3{3, 0, 0}, Int3{0, -2, 0},
                                Int3{0, 3, 0}, Int3{0, 0, -2}, Int3{0, 0, 3}}) {
    expect_error([&] {
      static_cast<void>(q(coordinate.x, coordinate.y, coordinate.z, 0));
    });
  }
  expect_error([&] { static_cast<void>(q(0, 0, 0, -1)); });
  expect_error([&] { static_cast<void>(q(0, 0, 0, 2)); });
}

void test_checked_size_overflow() {
  FieldRegistry element_overflow;
  element_overflow.declare_field(
      descriptor("element_overflow", ScalarType::uint8,
                 std::numeric_limits<std::uint32_t>::max(),
                 std::numeric_limits<int>::max()));
  element_overflow.freeze();
  expect_error([&] {
    FieldStorage storage(
        element_overflow,
        Int3{std::numeric_limits<int>::max(), std::numeric_limits<int>::max(),
             std::numeric_limits<int>::max()});
  });

  FieldRegistry byte_overflow;
  byte_overflow.declare_field(
      descriptor("byte_overflow", ScalarType::float64, 1, 0));
  byte_overflow.freeze();
  expect_error([&] {
    FieldStorage storage(byte_overflow,
                         Int3{std::numeric_limits<int>::max(),
                              std::numeric_limits<int>::max(), 1});
  });
}

void test_unindexable_component_count() {
  constexpr auto max_indexable_components =
      static_cast<std::uintmax_t>(std::numeric_limits<int>::max()) + 1U;
  constexpr auto first_unindexable_component_count =
      max_indexable_components + 1U;

  if constexpr (first_unindexable_component_count <=
                static_cast<std::uintmax_t>(
                    std::numeric_limits<std::uint32_t>::max())) {
    FieldRegistry registry;
    registry.declare_field(descriptor(
        "unindexable_components", ScalarType::uint8,
        static_cast<std::uint32_t>(first_unindexable_component_count), 0));
    registry.freeze();
    expect_error([&] { FieldStorage storage(registry, Int3{1, 1, 1}); });
  }
}

void test_unindexable_upper_ghost_coordinate() {
  FieldRegistry registry;
  registry.declare_field(
      descriptor("unindexable_upper_ghost", ScalarType::uint8, 1, 2));
  registry.freeze();

  constexpr std::array<Int3, 3> first_unindexable_extents{{
      {std::numeric_limits<int>::max(), 1, 1},
      {1, std::numeric_limits<int>::max(), 1},
      {1, 1, std::numeric_limits<int>::max()},
  }};
  for (const auto extent : first_unindexable_extents) {
    bool rejected_by_coordinate_domain_guard = false;
    try {
      FieldStorage storage(registry, extent);
    } catch (const Error &error) {
      rejected_by_coordinate_domain_guard =
          std::string(error.what()) ==
          "field storage upper ghost coordinate exceeds the view index range";
    }
    HUNDUN_CHECK(rejected_by_coordinate_domain_guard);
  }
}

void test_view_rejects_access_after_owner_destruction() {
  auto registry = make_epoch_registry();
  std::optional<FieldView<double>> surviving_view;
  {
    FieldStorage storage(registry, Int3{2, 1, 1});
    auto view = storage.view<double>(0U);
    view(1, 0, 0, 0) = 19.0;
    surviving_view.emplace(view);
  }

  HUNDUN_CHECK(surviving_view->interior_extent().x == 2);
  HUNDUN_CHECK(surviving_view->ghost_width() == 1);
  HUNDUN_CHECK(surviving_view->components() == 1U);
  expect_error_message(
      [&] { static_cast<void>((*surviving_view)(1, 0, 0, 0)); },
      kDeadOwnerMessage);
}

void test_lifecycle_transitions_invalidate_old_views() {
  auto registry = make_epoch_registry();
  FieldStorage storage(registry, Int3{2, 1, 1});
  auto current = storage.view<double>(0U);
  current(0, 0, 0, 0) = 23.0;

  using Transition = void (FieldStorage::*)();
  constexpr std::array<Transition, 3> transitions{{
      &FieldStorage::begin_rebuild,
      &FieldStorage::begin_repartition,
      &FieldStorage::begin_restart_v2_read_transaction,
  }};

  for (const auto transition : transitions) {
    const auto previous = current;
    (storage.*transition)();
    expect_error_message(
        [&] { static_cast<void>(previous(0, 0, 0, 0)); },
        kStaleGenerationMessage);

    current = storage.view<double>(0U);
    HUNDUN_CHECK_NEAR(current(0, 0, 0, 0), 23.0, 0.0);
    current(0, 0, 0, 0) = 23.0;
  }
}

void test_move_construction_preserves_source_views() {
  static_assert(std::is_nothrow_move_constructible_v<FieldStorage>);

  auto registry = make_epoch_registry();
  FieldStorage source(registry, Int3{2, 1, 1});
  auto source_view = source.view<double>(0U);
  source_view(1, 0, 0, 0) = 29.0;

  FieldStorage destination(std::move(source));
  HUNDUN_CHECK_NEAR(source_view(1, 0, 0, 0), 29.0, 0.0);
  auto destination_view = destination.view<double>(0U);
  HUNDUN_CHECK_NEAR(destination_view(1, 0, 0, 0), 29.0, 0.0);
}

void test_move_assignment_invalidates_target_and_preserves_source_views() {
  static_assert(std::is_nothrow_move_assignable_v<FieldStorage>);

  auto registry = make_epoch_registry();
  FieldStorage target(registry, Int3{2, 1, 1});
  FieldStorage source(registry, Int3{2, 1, 1});
  auto target_old_view = target.view<double>(0U);
  auto source_old_view = source.view<double>(0U);
  target_old_view(0, 0, 0, 0) = 31.0;
  source_old_view(0, 0, 0, 0) = 37.0;

  target = std::move(source);

  expect_error_message(
      [&] { static_cast<void>(target_old_view(0, 0, 0, 0)); },
      kDeadOwnerMessage);
  HUNDUN_CHECK_NEAR(source_old_view(0, 0, 0, 0), 37.0, 0.0);
  auto target_new_view = target.view<double>(0U);
  HUNDUN_CHECK_NEAR(target_new_view(0, 0, 0, 0), 37.0, 0.0);
}

void test_generation_is_monotonic_and_wrap_is_rejected() {
  auto registry = make_epoch_registry();
  FieldStorage storage(registry, Int3{2, 1, 1});

  const auto initial = FieldEpochTestAccess::generation(storage);
  storage.begin_rebuild();
  HUNDUN_CHECK(FieldEpochTestAccess::generation(storage) == initial + 1U);
  storage.begin_repartition();
  HUNDUN_CHECK(FieldEpochTestAccess::generation(storage) == initial + 2U);
  storage.begin_restart_v2_read_transaction();
  HUNDUN_CHECK(FieldEpochTestAccess::generation(storage) == initial + 3U);

  constexpr auto maximum = std::numeric_limits<std::uint64_t>::max();
  FieldEpochTestAccess::force_generation(storage, maximum);
  auto active_view = storage.view<double>(0U);
  active_view(0, 0, 0, 0) = 41.0;

  expect_error_message([&] { storage.begin_rebuild(); },
                       kGenerationWrapMessage);
  HUNDUN_CHECK(FieldEpochTestAccess::generation(storage) == maximum);
  HUNDUN_CHECK_NEAR(active_view(0, 0, 0, 0), 41.0, 0.0);
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_registry_validation_and_stable_ids();
    test_storage_preconditions_and_function_spaces();
    test_typed_views_layout_bounds_and_constness();
    test_checked_size_overflow();
    test_unindexable_component_count();
    test_unindexable_upper_ghost_coordinate();
    test_view_rejects_access_after_owner_destruction();
    test_lifecycle_transitions_invalidate_old_views();
    test_move_construction_preserves_source_views();
    test_move_assignment_invalidates_target_and_preserves_source_views();
    test_generation_is_monotonic_and_wrap_is_rejected();
  });
}
