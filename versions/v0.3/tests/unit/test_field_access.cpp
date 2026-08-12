// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

#include "hundun/rt_error.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
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
using hundun::runtime::OutputPolicy;
using hundun::runtime::PhaseId;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;

constexpr const char *kPlanRequiresFrozenRegistry =
    "field access plan requires a frozen registry";
constexpr const char *kPlanFrozen = "field access plan is frozen";
constexpr const char *kUnknownField =
    "field access declaration references an unknown field";
constexpr const char *kInvalidMode =
    "field access declaration mode is unrecognized";
constexpr const char *kDuplicateDeclaration =
    "field access declaration is duplicated";
constexpr const char *kWriterConflict =
    "field access plan already has a writer for this phase and field";
constexpr const char *kUnfrozenAcquisition =
    "field access acquisition requires a frozen plan";
constexpr const char *kDomainMismatch =
    "field access plan field domain does not match storage";
constexpr const char *kUndeclaredAccess = "field access was not declared";
constexpr const char *kReadNotPermitted =
    "field access declaration does not permit reading";
constexpr const char *kWriteNotPermitted =
    "field access declaration does not permit writing";
constexpr const char *kScalarMismatch =
    "field view scalar type does not match its descriptor";
constexpr const char *kCellSpaceMismatch =
    "field acquisition requires a cell_average field";
constexpr const char *kFaceSpaceMismatch =
    "field acquisition requires a face_value field";
constexpr const char *kFaceCountRequired =
    "face field storage requires a positive face count";
constexpr const char *kFaceGhostUnsupported =
    "face_value fields require zero ghost width";
constexpr const char *kLegacyCellOnly =
    "passive-scalar field storage supports cell_average fields only";
constexpr const char *kUnsupportedSpace =
    "field storage does not support this function space";
constexpr const char *kDeadOwner = "field view owner is no longer alive";
constexpr const char *kStaleGeneration = "field view generation is stale";

FieldDescriptor descriptor(
    std::string name, FunctionSpace space = FunctionSpace::cell_average,
    ScalarType scalar_type = ScalarType::float64,
    std::uint32_t components = 1U, int ghost_width = 0) {
  return FieldDescriptor{std::move(name),
                         "1",
                         "stage2_field_test",
                         space,
                         scalar_type,
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

template <class Function>
void expect_error(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &) {
    threw = true;
  }
  HUNDUN_CHECK(threw);
}

FieldRegistry one_cell_registry() {
  FieldRegistry registry;
  registry.declare_field(descriptor("cell"));
  registry.freeze();
  return registry;
}

FieldRegistry one_face_registry(int ghost_width = 0,
                                std::uint32_t components = 2U,
                                ScalarType scalar_type = ScalarType::float64) {
  FieldRegistry registry;
  registry.declare_field(descriptor("face", FunctionSpace::face_value,
                                    scalar_type, components, ghost_width));
  registry.freeze();
  return registry;
}

void test_plan_validation_and_freeze() {
  FieldRegistry unfrozen;
  const FieldId field = unfrozen.declare_field(descriptor("q"));
  expect_error_message([&] { FieldAccessPlan plan(unfrozen); },
                       kPlanRequiresFrozenRegistry);

  unfrozen.freeze();
  FieldAccessPlan plan(unfrozen);
  HUNDUN_CHECK(!plan.frozen());

  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{0}, ActorId{1},
                            std::numeric_limits<FieldId>::max(),
                            AccessMode::read);
      },
      kUnknownField);
  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{0}, ActorId{1}, field,
                            static_cast<AccessMode>(-1));
      },
      kInvalidMode);
  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{0}, ActorId{1}, field,
                            static_cast<AccessMode>(3));
      },
      kInvalidMode);

  plan.declare_access(PhaseId{0}, ActorId{1}, field, AccessMode::read);
  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{0}, ActorId{1}, field,
                            AccessMode::read_write);
      },
      kDuplicateDeclaration);

  plan.freeze();
  plan.freeze();
  HUNDUN_CHECK(plan.frozen());
  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{1}, ActorId{2}, field,
                            AccessMode::read);
      },
      kPlanFrozen);
}

void test_single_writer_and_reader_rules() {
  auto registry = one_cell_registry();
  const FieldId field = 0U;
  FieldAccessPlan plan(registry);

  plan.declare_access(PhaseId{2}, ActorId{1}, field, AccessMode::read);
  plan.declare_access(PhaseId{2}, ActorId{2}, field, AccessMode::read);
  plan.declare_access(PhaseId{2}, ActorId{3}, field, AccessMode::write);
  plan.declare_access(PhaseId{3}, ActorId{4}, field, AccessMode::read_write);

  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{2}, ActorId{5}, field,
                            AccessMode::read_write);
      },
      kWriterConflict);
  expect_error_message(
      [&] {
        plan.declare_access(PhaseId{3}, ActorId{5}, field,
                            AccessMode::write);
      },
      kWriterConflict);

  plan.declare_access(PhaseId{2}, ActorId{6}, field, AccessMode::read);
  plan.declare_access(PhaseId{3}, ActorId{6}, field, AccessMode::read);
  plan.freeze();
}

void test_cell_acquisition_matrix_and_domain_checks() {
  FieldRegistry registry;
  const FieldId read_id = registry.declare_field(
      descriptor("read", FunctionSpace::cell_average, ScalarType::float64));
  const FieldId write_id = registry.declare_field(
      descriptor("write", FunctionSpace::cell_average, ScalarType::int32));
  const FieldId read_write_id = registry.declare_field(descriptor(
      "read_write", FunctionSpace::cell_average, ScalarType::float64));
  registry.freeze();

  constexpr PhaseId phase = 7U;
  constexpr ActorId reader = 11U;
  constexpr ActorId writer = 12U;
  constexpr ActorId both = 13U;
  FieldAccessPlan plan(registry);
  plan.declare_access(phase, reader, read_id, AccessMode::read);
  plan.declare_access(phase, writer, write_id, AccessMode::write);
  plan.declare_access(phase, both, read_write_id, AccessMode::read_write);

  FieldStorage storage(registry, Int3{2, 1, 1});
  expect_error_message(
      [&] {
        static_cast<void>(
            storage.acquire_read<double>(plan, phase, reader, read_id));
      },
      kUnfrozenAcquisition);
  plan.freeze();

  auto read_view = storage.acquire_read<double>(plan, phase, reader, read_id);
  auto write_view =
      storage.acquire_write<std::int32_t>(plan, phase, writer, write_id);
  auto read_write_read =
      storage.acquire_read<double>(plan, phase, both, read_write_id);
  auto read_write_write =
      storage.acquire_write<double>(plan, phase, both, read_write_id);
  static_assert(
      std::is_same_v<decltype(read_view(0, 0, 0, 0)), const double &>);
  static_assert(std::is_same_v<decltype(write_view(0, 0, 0, 0)),
                               std::int32_t &>);
  static_assert(std::is_same_v<decltype(read_write_read(0, 0, 0, 0)),
                               const double &>);
  static_assert(std::is_same_v<decltype(read_write_write(0, 0, 0, 0)),
                               double &>);

  write_view(1, 0, 0, 0) = 17;
  read_write_write(0, 0, 0, 0) = 19.0;
  HUNDUN_CHECK(write_view(1, 0, 0, 0) == 17);
  HUNDUN_CHECK_NEAR(read_write_read(0, 0, 0, 0), 19.0, 0.0);

  const FieldStorage &const_storage = storage;
  auto const_read =
      const_storage.acquire_read<double>(plan, phase, reader, read_id);
  static_assert(
      std::is_same_v<decltype(const_read(0, 0, 0, 0)), const double &>);

  expect_error_message(
      [&] {
        static_cast<void>(
            storage.acquire_write<double>(plan, phase, reader, read_id));
      },
      kWriteNotPermitted);
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_read<std::int32_t>(
            plan, phase, writer, write_id));
      },
      kReadNotPermitted);
  expect_error_message(
      [&] {
        static_cast<void>(
            storage.acquire_read<double>(plan, phase, ActorId{99}, read_id));
      },
      kUndeclaredAccess);
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_read<double>(
            plan, PhaseId{99}, reader, read_id));
      },
      kUndeclaredAccess);
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_read<double>(plan, phase, reader,
                                                       read_write_id));
      },
      kUndeclaredAccess);
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_read<std::int32_t>(
            plan, phase, reader, read_id));
      },
      kScalarMismatch);

  FieldRegistry shorter_registry;
  shorter_registry.declare_field(descriptor("only"));
  shorter_registry.freeze();
  FieldStorage shorter_storage(shorter_registry, Int3{2, 1, 1});
  expect_error_message(
      [&] {
        static_cast<void>(shorter_storage.acquire_read<double>(
            plan, phase, reader, read_id));
      },
      kDomainMismatch);
}

void test_plan_and_registry_lifetimes_are_not_borrowed() {
  std::optional<FieldStorage> storage;
  std::optional<FieldAccessPlan> surviving_plan;
  {
    auto registry = one_cell_registry();
    storage.emplace(registry, Int3{1, 1, 1});
    surviving_plan.emplace(registry);
    surviving_plan->declare_access(PhaseId{1}, ActorId{2}, FieldId{0},
                                   AccessMode::read_write);
    surviving_plan->freeze();
  }

  auto writer = storage->acquire_write<double>(
      *surviving_plan, PhaseId{1}, ActorId{2}, FieldId{0});
  writer(0, 0, 0, 0) = 23.0;

  std::optional<FieldView<const double>> surviving_view;
  {
    FieldAccessPlan temporary_plan = *surviving_plan;
    surviving_view.emplace(storage->acquire_read<double>(
        temporary_plan, PhaseId{1}, ActorId{2}, FieldId{0}));
  }
  HUNDUN_CHECK_NEAR((*surviving_view)(0, 0, 0, 0), 23.0, 0.0);
}

void test_legacy_constructor_stays_cell_only() {
  auto cell_registry = one_cell_registry();
  FieldStorage legacy(cell_registry, Int3{2, 2, 1});
  auto cell = legacy.view<double>(0U);
  cell(1, 1, 0, 0) = 29.0;
  HUNDUN_CHECK_NEAR(cell(1, 1, 0, 0), 29.0, 0.0);

  auto face_registry = one_face_registry();
  expect_error_message(
      [&] { FieldStorage rejected(face_registry, Int3{2, 2, 1}); },
      kLegacyCellOnly);
}

void test_mixed_cell_face_layout_and_acquisition() {
  FieldRegistry registry;
  const FieldId cell_id = registry.declare_field(descriptor(
      "cell", FunctionSpace::cell_average, ScalarType::float64, 2U, 1));
  const FieldId face_id = registry.declare_field(descriptor(
      "face", FunctionSpace::face_value, ScalarType::float64, 3U, 0));
  const FieldId marker_id = registry.declare_field(descriptor(
      "marker", FunctionSpace::face_value, ScalarType::uint8, 1U, 0));
  registry.freeze();

  const FieldLayoutSet layout{Int3{2, 1, 1}, 4U};
  FieldStorage storage(registry, layout);
  const auto actual_layout = storage.layout_set();
  HUNDUN_CHECK(actual_layout.cell_interior_extent.x == 2);
  HUNDUN_CHECK(actual_layout.cell_interior_extent.y == 1);
  HUNDUN_CHECK(actual_layout.cell_interior_extent.z == 1);
  HUNDUN_CHECK(actual_layout.face_count == 4U);
  HUNDUN_CHECK(storage.interior_extent().x == 2);

  constexpr PhaseId phase = 4U;
  constexpr ActorId actor = 8U;
  FieldAccessPlan plan(registry);
  plan.declare_access(phase, actor, cell_id, AccessMode::read_write);
  plan.declare_access(phase, actor, face_id, AccessMode::read_write);
  plan.declare_access(phase, actor, marker_id, AccessMode::read_write);
  plan.freeze();

  auto cell = storage.acquire_write<double>(plan, phase, actor, cell_id);
  auto face = storage.acquire_face_write<double>(plan, phase, actor, face_id);
  auto marker =
      storage.acquire_face_write<std::uint8_t>(plan, phase, actor, marker_id);
  static_assert(
      std::is_same_v<decltype(face(std::size_t{0}, 0)), double &>);
  HUNDUN_CHECK(face.face_count() == 4U);
  HUNDUN_CHECK(face.components() == 3U);

  cell(-1, -1, -1, 1) = 31.0;
  face(0U, 0) = 37.0;
  face(0U, 1) = 41.0;
  face(1U, 0) = 43.0;
  face(3U, 2) = 47.0;
  marker(2U, 0) = static_cast<std::uint8_t>(251);
  HUNDUN_CHECK_NEAR(cell(-1, -1, -1, 1), 31.0, 0.0);
  HUNDUN_CHECK_NEAR(face(3U, 2), 47.0, 0.0);
  HUNDUN_CHECK(marker(2U, 0) == static_cast<std::uint8_t>(251));

  const auto base = reinterpret_cast<std::uintptr_t>(&face(0U, 0));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&face(0U, 1)) ==
               base + sizeof(double));
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(&face(1U, 0)) ==
               base + 3U * sizeof(double));
  HUNDUN_CHECK(base % alignof(double) == 0U);

  const FieldStorage &const_storage = storage;
  auto const_face = const_storage.acquire_face_read<double>(
      plan, phase, actor, face_id);
  static_assert(std::is_same_v<decltype(const_face(std::size_t{0}, 0)),
                               const double &>);
  HUNDUN_CHECK_NEAR(const_face(1U, 0), 43.0, 0.0);

  expect_error_message(
      [&] { static_cast<void>(face(4U, 0)); },
      "face field view face index is out of bounds");
  expect_error_message(
      [&] { static_cast<void>(face(0U, -1)); },
      "face field view component index is out of bounds");
  expect_error_message(
      [&] { static_cast<void>(face(0U, 3)); },
      "face field view component index is out of bounds");
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_read<double>(plan, phase, actor,
                                                       face_id));
      },
      kCellSpaceMismatch);
  expect_error_message(
      [&] { static_cast<void>(storage.view<double>(face_id)); },
      kCellSpaceMismatch);
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_face_read<double>(
            plan, phase, actor, cell_id));
      },
      kFaceSpaceMismatch);
  expect_error_message(
      [&] {
        static_cast<void>(storage.acquire_face_read<std::int32_t>(
            plan, phase, actor, face_id));
      },
      kScalarMismatch);
}

void test_face_layout_preconditions_and_overflow() {
  auto face_registry = one_face_registry();
  expect_error_message(
      [&] {
        FieldStorage storage(face_registry,
                             FieldLayoutSet{Int3{1, 1, 1}, 0U});
      },
      kFaceCountRequired);

  auto ghosted_face_registry = one_face_registry(1);
  expect_error_message(
      [&] {
        FieldStorage storage(ghosted_face_registry,
                             FieldLayoutSet{Int3{1, 1, 1}, 2U});
      },
      kFaceGhostUnsupported);

  for (const auto extent : {Int3{0, 1, 1}, Int3{-1, 1, 1}, Int3{1, 0, 1},
                            Int3{1, -1, 1}, Int3{1, 1, 0},
                            Int3{1, 1, -1}}) {
    expect_error([&] {
      FieldStorage storage(face_registry, FieldLayoutSet{extent, 2U});
    });
  }

  for (const auto space :
       {FunctionSpace::vertex_value, FunctionSpace::element_dof,
        FunctionSpace::quadrature_point, FunctionSpace::particle}) {
    FieldRegistry unsupported;
    unsupported.declare_field(descriptor("unsupported", space));
    unsupported.freeze();
    expect_error_message(
        [&] {
          FieldStorage storage(unsupported,
                               FieldLayoutSet{Int3{1, 1, 1}, 1U});
        },
        kUnsupportedSpace);
  }

  auto component_overflow = one_face_registry(
      0, std::numeric_limits<std::uint32_t>::max(), ScalarType::uint8);
  expect_error([&] {
    FieldStorage storage(component_overflow,
                         FieldLayoutSet{Int3{1, 1, 1},
                                        std::numeric_limits<std::size_t>::max()});
  });

  auto element_overflow = one_face_registry(0, 2U, ScalarType::uint8);
  expect_error([&] {
    FieldStorage storage(element_overflow,
                         FieldLayoutSet{Int3{1, 1, 1},
                                        std::numeric_limits<std::size_t>::max()});
  });

  auto byte_overflow = one_face_registry(0, 1U, ScalarType::float64);
  const auto first_overflowing_double_count =
      std::numeric_limits<std::size_t>::max() / sizeof(double) + 1U;
  expect_error([&] {
    FieldStorage storage(
        byte_overflow,
        FieldLayoutSet{Int3{1, 1, 1}, first_overflowing_double_count});
  });

  auto capacity_overflow = one_face_registry(0, 1U, ScalarType::uint8);
  expect_error([&] {
    FieldStorage storage(
        capacity_overflow,
        FieldLayoutSet{Int3{1, 1, 1},
                       std::numeric_limits<std::size_t>::max()});
  });
}

FieldAccessPlan make_face_plan(const FieldRegistry &registry) {
  FieldAccessPlan plan(registry);
  plan.declare_access(PhaseId{5}, ActorId{6}, FieldId{0},
                      AccessMode::read_write);
  plan.freeze();
  return plan;
}

void test_face_view_epoch_transitions_and_destruction() {
  auto registry = one_face_registry();
  auto plan = make_face_plan(registry);
  std::optional<FaceFieldView<double>> surviving;
  {
    FieldStorage storage(registry, FieldLayoutSet{Int3{1, 1, 1}, 3U});
    auto face = storage.acquire_face_write<double>(
        plan, PhaseId{5}, ActorId{6}, FieldId{0});
    face(2U, 1) = 53.0;
    surviving.emplace(face);
  }
  HUNDUN_CHECK(surviving->face_count() == 3U);
  HUNDUN_CHECK(surviving->components() == 2U);
  expect_error_message([&] { static_cast<void>((*surviving)(2U, 1)); },
                       kDeadOwner);

  FieldStorage storage(registry, FieldLayoutSet{Int3{1, 1, 1}, 3U});
  auto current = storage.acquire_face_write<double>(
      plan, PhaseId{5}, ActorId{6}, FieldId{0});
  current(1U, 0) = 59.0;

  using Transition = void (FieldStorage::*)();
  constexpr std::array<Transition, 3> transitions{{
      &FieldStorage::begin_rebuild,
      &FieldStorage::begin_repartition,
      &FieldStorage::begin_restart_v2_read_transaction,
  }};
  for (const auto transition : transitions) {
    const auto previous = current;
    (storage.*transition)();
    HUNDUN_CHECK(previous.face_count() == 3U);
    HUNDUN_CHECK(previous.components() == 2U);
    expect_error_message([&] { static_cast<void>(previous(1U, 0)); },
                         kStaleGeneration);
    current = storage.acquire_face_write<double>(
        plan, PhaseId{5}, ActorId{6}, FieldId{0});
    HUNDUN_CHECK_NEAR(current(1U, 0), 59.0, 0.0);
  }
}

void test_face_view_move_replacement_behavior() {
  static_assert(std::is_nothrow_move_constructible_v<FieldStorage>);
  static_assert(std::is_nothrow_move_assignable_v<FieldStorage>);

  auto registry = one_face_registry();
  auto plan = make_face_plan(registry);
  const FieldLayoutSet layout{Int3{1, 1, 1}, 2U};

  FieldStorage source(registry, layout);
  auto source_view = source.acquire_face_write<double>(
      plan, PhaseId{5}, ActorId{6}, FieldId{0});
  source_view(1U, 1) = 61.0;
  FieldStorage constructed(std::move(source));
  HUNDUN_CHECK_NEAR(source_view(1U, 1), 61.0, 0.0);
  auto constructed_view = constructed.acquire_face_write<double>(
      plan, PhaseId{5}, ActorId{6}, FieldId{0});
  HUNDUN_CHECK_NEAR(constructed_view(1U, 1), 61.0, 0.0);

  FieldStorage target(registry, layout);
  FieldStorage replacement(registry, layout);
  auto target_old = target.acquire_face_write<double>(
      plan, PhaseId{5}, ActorId{6}, FieldId{0});
  auto replacement_old = replacement.acquire_face_write<double>(
      plan, PhaseId{5}, ActorId{6}, FieldId{0});
  target_old(0U, 0) = 67.0;
  replacement_old(0U, 0) = 71.0;
  target = std::move(replacement);

  expect_error_message([&] { static_cast<void>(target_old(0U, 0)); },
                       kDeadOwner);
  HUNDUN_CHECK_NEAR(replacement_old(0U, 0), 71.0, 0.0);
  auto target_new = target.acquire_face_write<double>(
      plan, PhaseId{5}, ActorId{6}, FieldId{0});
  HUNDUN_CHECK_NEAR(target_new(0U, 0), 71.0, 0.0);
  HUNDUN_CHECK(target.layout_set().face_count == 2U);
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_plan_validation_and_freeze();
    test_single_writer_and_reader_rules();
    test_cell_acquisition_matrix_and_domain_checks();
    test_plan_and_registry_lifetimes_are_not_borrowed();
    test_legacy_constructor_stays_cell_only();
    test_mixed_cell_face_layout_and_acquisition();
    test_face_layout_preconditions_and_overflow();
    test_face_view_epoch_transitions_and_destruction();
    test_face_view_move_replacement_behavior();
  });
}
