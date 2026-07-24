// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/flow_state.hpp"
#include "adaptive_time_control_detail.hpp"
#include "checkpoint_v2_detail.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "material_density_transport_test_access.hpp"
#include "material_density_piso_test_access.hpp"
#endif

#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace hundun::flow {
namespace {

constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;

std::size_t cell_count(runtime::Int3 extent) {
  const auto x = static_cast<std::size_t>(extent.x);
  const auto y = static_cast<std::size_t>(extent.y);
  const auto z = static_cast<std::size_t>(extent.z);
  if (y != 0U && x > std::numeric_limits<std::size_t>::max() / y) {
    throw runtime::Error("flow-state cell count overflows");
  }
  const std::size_t xy = x * y;
  if (z != 0U && xy > std::numeric_limits<std::size_t>::max() / z) {
    throw runtime::Error("flow-state cell count overflows");
  }
  return xy * z;
}

void require_cell_field(const runtime::FieldRegistry &registry,
                        runtime::FieldId id, std::uint32_t components,
                        const char *name) {
  const auto &descriptor = registry.descriptor(id);
  if (descriptor.space != runtime::FunctionSpace::cell_average ||
      descriptor.scalar_type != runtime::ScalarType::float64 ||
      descriptor.components != components || descriptor.ghost_width < 2) {
    throw runtime::Error(std::string("flow-state ") + name +
                         " descriptor is invalid");
  }
}

void require_face_field(const runtime::FieldRegistry &registry,
                        runtime::FieldId id, std::uint32_t components,
                        const char *name) {
  const auto &descriptor = registry.descriptor(id);
  if (descriptor.space != runtime::FunctionSpace::face_value ||
      descriptor.scalar_type != runtime::ScalarType::float64 ||
      descriptor.components != components || descriptor.ghost_width != 0) {
    throw runtime::Error(std::string("flow-state ") + name +
                         " descriptor is invalid");
  }
}

std::vector<runtime::FieldId> ordered_fields(const FlowFieldIds &fields) {
  std::vector<runtime::FieldId> result{
      fields.density, fields.velocity, fields.mechanical_pressure,
      fields.face_velocity, fields.face_mass_flux};
  result.insert(result.end(), fields.transported_cell_fields.begin(),
                fields.transported_cell_fields.end());
  return result;
}

void validate_fields(const runtime::FieldRegistry &registry,
                     const FlowFieldIds &fields) {
  if (!registry.frozen()) {
    throw runtime::Error("FlowState requires a frozen field registry");
  }
  require_cell_field(registry, fields.density, 1U, "density");
  require_cell_field(registry, fields.velocity, 3U, "velocity");
  require_cell_field(registry, fields.mechanical_pressure, 1U,
                     "mechanical-pressure");
  require_face_field(registry, fields.face_velocity, 3U, "face-velocity");
  finite_volume::require_face_mass_flux_field(registry, fields.face_mass_flux);
  for (const auto field : fields.transported_cell_fields) {
    require_cell_field(registry, field, 1U, "transported-field");
  }
  const auto ids = ordered_fields(fields);
  const std::set<runtime::FieldId> unique(ids.begin(), ids.end());
  if (unique.size() != ids.size()) {
    throw runtime::Error("flow-state field identifiers must be unique");
  }
}

void validate_metadata(AcceptedStepMetadata metadata) {
  if (!std::isfinite(metadata.time_s) || metadata.time_s < 0.0 ||
      !(metadata.dt_s > 0.0) || !std::isfinite(metadata.dt_s)) {
    throw runtime::Error("flow-state accepted metadata is invalid");
  }
  static_cast<void>(make_momentum_time_stencil(metadata.order, metadata.dt_s,
                                               metadata.previous_dt_s));
}

void require_size(const std::vector<double> &values, std::size_t expected,
                  const char *name) {
  if (values.size() != expected) {
    throw runtime::Error(std::string("flow-state ") + name +
                         " value count is invalid");
  }
  if (!std::all_of(values.begin(), values.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw runtime::Error(std::string("flow-state ") + name +
                         " contains a non-finite value");
  }
}

template <class Function>
void for_each_cell(runtime::Int3 extent, Function &&function) {
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        function(i, j, k);
      }
    }
  }
}

} // namespace

FlowState::Impl::Impl(const runtime::FieldRegistry &supplied_registry,
                      runtime::FieldLayoutSet supplied_layout,
                      FlowFieldIds supplied_fields,
                      AcceptedStepMetadata supplied_metadata)
    : registry(&supplied_registry), layout(supplied_layout),
      fields(std::move(supplied_fields)), metadata(supplied_metadata),
      access(supplied_registry), history(supplied_registry, supplied_layout),
      committed(supplied_registry, supplied_layout),
      trial(supplied_registry, supplied_layout),
      rollback_snapshot(supplied_registry, supplied_layout) {
  for (const auto field : ordered_fields(fields)) {
    access.declare_access(kStatePhase, kStateActor, field,
                          runtime::AccessMode::read_write);
  }
  access.freeze();
}

static_assert(std::is_nothrow_move_constructible_v<runtime::FieldStorage>);
static_assert(std::is_nothrow_move_assignable_v<runtime::FieldStorage>);
static_assert(std::is_nothrow_copy_assignable_v<AcceptedStepMetadata>);

namespace {

template <class StateImplementation>
void write_cell_values(StateImplementation &impl,
                       runtime::FieldStorage &storage, runtime::FieldId field,
                       const std::vector<double> &values,
                       std::uint32_t components) {
  const std::size_t count = cell_count(impl.layout.cell_interior_extent);
  if (components != 0U &&
      count > std::numeric_limits<std::size_t>::max() / components) {
    throw runtime::Error("flow-state cell value count overflows");
  }
  require_size(values, count * components, "cell");
  auto view = storage.acquire_write<double>(impl.access, kStatePhase,
                                            kStateActor, field);
  std::size_t offset = 0U;
  for_each_cell(impl.layout.cell_interior_extent, [&](int i, int j, int k) {
    for (std::uint32_t component = 0; component < components; ++component) {
      view(i, j, k, static_cast<int>(component)) = values[offset++];
    }
  });
}

template <class StateImplementation>
std::vector<double> read_cell_values(const StateImplementation &impl,
                                     const runtime::FieldStorage &storage,
                                     runtime::FieldId field,
                                     std::uint32_t components) {
  std::vector<double> result;
  result.reserve(cell_count(impl.layout.cell_interior_extent) * components);
  const auto view = storage.acquire_read<double>(impl.access, kStatePhase,
                                                 kStateActor, field);
  for_each_cell(impl.layout.cell_interior_extent, [&](int i, int j, int k) {
    for (std::uint32_t component = 0; component < components; ++component) {
      result.push_back(view(i, j, k, static_cast<int>(component)));
    }
  });
  return result;
}

template <class StateImplementation>
void write_face_values(StateImplementation &impl,
                       runtime::FieldStorage &storage, runtime::FieldId field,
                       const std::vector<double> &values,
                       std::uint32_t components) {
  if (impl.layout.face_count >
      std::numeric_limits<std::size_t>::max() / components) {
    throw runtime::Error("flow-state face value count overflows");
  }
  require_size(values, impl.layout.face_count * components, "face");
  auto view = storage.acquire_face_write<double>(impl.access, kStatePhase,
                                                 kStateActor, field);
  for (std::size_t face = 0; face < impl.layout.face_count; ++face) {
    for (std::uint32_t component = 0; component < components; ++component) {
      view(face, static_cast<int>(component)) =
          values[face * components + component];
    }
  }
}

template <class StateImplementation>
std::vector<double> read_face_values(const StateImplementation &impl,
                                     const runtime::FieldStorage &storage,
                                     runtime::FieldId field,
                                     std::uint32_t components) {
  std::vector<double> result(impl.layout.face_count * components);
  const auto view = storage.acquire_face_read<double>(impl.access, kStatePhase,
                                                      kStateActor, field);
  for (std::size_t face = 0; face < impl.layout.face_count; ++face) {
    for (std::uint32_t component = 0; component < components; ++component) {
      result[face * components + component] =
          view(face, static_cast<int>(component));
    }
  }
  return result;
}

template <class StateImplementation>
void validate_layer_values(const StateImplementation &impl,
                           const FlowLayerValues &values) {
  const std::size_t cells = cell_count(impl.layout.cell_interior_extent);
  require_size(values.density, cells, "density");
  require_size(values.velocity, cells * 3U, "velocity");
  require_size(values.mechanical_pressure, cells, "mechanical-pressure");
  require_size(values.face_velocity, impl.layout.face_count * 3U,
               "face-velocity");
  require_size(values.face_mass_flux, impl.layout.face_count, "face-mass-flux");
  if (values.transported_cell_fields.size() !=
      impl.fields.transported_cell_fields.size()) {
    throw runtime::Error("flow-state transported-field count is invalid");
  }
  for (const auto &transported : values.transported_cell_fields) {
    require_size(transported, cells, "transported-field");
  }
}

template <class StateImplementation>
void write_layer(StateImplementation &impl, runtime::FieldStorage &storage,
                 const FlowLayerValues &values) {
  write_cell_values(impl, storage, impl.fields.density, values.density, 1U);
  write_cell_values(impl, storage, impl.fields.velocity, values.velocity, 3U);
  write_cell_values(impl, storage, impl.fields.mechanical_pressure,
                    values.mechanical_pressure, 1U);
  write_face_values(impl, storage, impl.fields.face_velocity,
                    values.face_velocity, 3U);
  write_face_values(impl, storage, impl.fields.face_mass_flux,
                    values.face_mass_flux, 1U);
  for (std::size_t index = 0;
       index < impl.fields.transported_cell_fields.size(); ++index) {
    write_cell_values(impl, storage, impl.fields.transported_cell_fields[index],
                      values.transported_cell_fields[index], 1U);
  }
}

template <class StateImplementation>
FlowLayerValues read_layer(const StateImplementation &impl,
                           const runtime::FieldStorage &storage) {
  FlowLayerValues values;
  values.density = read_cell_values(impl, storage, impl.fields.density, 1U);
  values.velocity = read_cell_values(impl, storage, impl.fields.velocity, 3U);
  values.mechanical_pressure =
      read_cell_values(impl, storage, impl.fields.mechanical_pressure, 1U);
  values.face_velocity =
      read_face_values(impl, storage, impl.fields.face_velocity, 3U);
  values.face_mass_flux =
      read_face_values(impl, storage, impl.fields.face_mass_flux, 1U);
  for (const auto field : impl.fields.transported_cell_fields) {
    values.transported_cell_fields.push_back(
        read_cell_values(impl, storage, field, 1U));
  }
  return values;
}

template <class StateImplementation>
void copy_layer(StateImplementation &impl, const runtime::FieldStorage &source,
                runtime::FieldStorage &destination) {
  const auto copy_cell = [&](runtime::FieldId field, std::uint32_t components) {
    const auto input = source.template acquire_read<double>(
        impl.access, kStatePhase, kStateActor, field);
    auto output = destination.template acquire_write<double>(
        impl.access, kStatePhase, kStateActor, field);
    for_each_cell(impl.layout.cell_interior_extent, [&](int i, int j, int k) {
      for (std::uint32_t component = 0; component < components; ++component) {
        output(i, j, k, static_cast<int>(component)) =
            input(i, j, k, static_cast<int>(component));
      }
    });
  };
  const auto copy_face = [&](runtime::FieldId field, std::uint32_t components) {
    const auto input = source.template acquire_face_read<double>(
        impl.access, kStatePhase, kStateActor, field);
    auto output = destination.template acquire_face_write<double>(
        impl.access, kStatePhase, kStateActor, field);
    for (std::size_t face = 0; face < impl.layout.face_count; ++face) {
      for (std::uint32_t component = 0; component < components; ++component) {
        output(face, static_cast<int>(component)) =
            input(face, static_cast<int>(component));
      }
    }
  };
  copy_cell(impl.fields.density, 1U);
  copy_cell(impl.fields.velocity, 3U);
  copy_cell(impl.fields.mechanical_pressure, 1U);
  copy_face(impl.fields.face_velocity, 3U);
  copy_face(impl.fields.face_mass_flux, 1U);
  for (const runtime::FieldId field : impl.fields.transported_cell_fields) {
    copy_cell(field, 1U);
  }
}

} // namespace

FlowState FlowState::create(const runtime::FieldRegistry &registry,
                            runtime::FieldLayoutSet layout, FlowFieldIds fields,
                            AcceptedStepMetadata metadata) {
  validate_fields(registry, fields);
  validate_metadata(metadata);
  if (layout.face_count == 0U) {
    throw runtime::Error("FlowState requires a positive local face count");
  }
  static_cast<void>(cell_count(layout.cell_interior_extent));
  return FlowState(
      std::make_unique<Impl>(registry, layout, std::move(fields), metadata));
}

FlowState::FlowState(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
FlowState::~FlowState() noexcept = default;
FlowState::FlowState(FlowState &&) noexcept = default;
FlowState &FlowState::operator=(FlowState &&) noexcept = default;

const runtime::FieldRegistry &FlowState::solver_registry() const noexcept {
  return *impl_->registry;
}

const runtime::FieldAccessPlan &FlowState::solver_access_plan() const noexcept {
  return impl_->access;
}

runtime::FieldStorage &FlowState::solver_layer(FlowLayer selected) {
  switch (selected) {
  case FlowLayer::history:
    return impl_->history;
  case FlowLayer::committed:
    return impl_->committed;
  case FlowLayer::trial:
    return impl_->trial;
  }
  throw runtime::Error("invalid flow-state layer");
}

bool FlowState::attempt_active() const noexcept {
  return impl_->attempt_active;
}

std::uint64_t FlowState::attempt_identity() const noexcept {
  return impl_->attempt_identity;
}

std::uint64_t FlowState::diagnostic_mutation_identity() const noexcept {
  return impl_ == nullptr ? 0U : impl_->diagnostic_mutation_identity;
}

const runtime::FieldStorage &FlowState::layer(FlowLayer selected) const {
  switch (selected) {
  case FlowLayer::history:
    return impl_->history;
  case FlowLayer::committed:
    return impl_->committed;
  case FlowLayer::trial:
    return impl_->trial;
  }
  throw runtime::Error("invalid flow-state layer");
}

runtime::FieldStorage &FlowState::trial_layer() { return impl_->trial; }

const FlowFieldIds &FlowState::fields() const noexcept { return impl_->fields; }

AcceptedStepMetadata FlowState::metadata() const noexcept {
  return impl_->metadata;
}

const FlowState *
detail::AdaptiveTimeControlAccess::state_identity(const FlowState &state) noexcept {
  return &state;
}

bool detail::AdaptiveTimeControlAccess::state_live(
    const FlowState &state) noexcept {
  return static_cast<bool>(state.impl_);
}

bool detail::AdaptiveTimeControlAccess::state_layout_matches(
    const FlowState &state, const mesh::MeshTopology &topology) noexcept {
  if (!state.impl_)
    return false;
  const auto box = topology.owned_global_box();
  const runtime::Int3 expected{box.end.x - box.begin.x,
                              box.end.y - box.begin.y,
                              box.end.z - box.begin.z};
  const auto actual = state.impl_->layout.cell_interior_extent;
  return actual.x == expected.x && actual.y == expected.y &&
         actual.z == expected.z &&
         state.impl_->layout.face_count == topology.local_face_count();
}

std::uint64_t detail::AdaptiveTimeControlAccess::diagnostic_identity(
    const FlowState &state) noexcept {
  return state.diagnostic_mutation_identity();
}

std::vector<double>
detail::AdaptiveTimeControlAccess::committed_density(const FlowState &state) {
  const auto extent = state.impl_->layout.cell_interior_extent;
  std::vector<double> result;
  result.reserve(cell_count(extent));
  auto view = state.impl_->committed.acquire_read<double>(
      state.impl_->access, kStatePhase, kStateActor, state.impl_->fields.density);
  for_each_cell(extent, [&](int i, int j, int k) {
    result.push_back(view(i, j, k, 0));
  });
  return result;
}

std::vector<double>
detail::AdaptiveTimeControlAccess::committed_face_mass_flux(
    const FlowState &state) {
  auto view = state.impl_->committed.acquire_face_read<double>(
      state.impl_->access, kStatePhase, kStateActor,
      state.impl_->fields.face_mass_flux);
  std::vector<double> result(view.face_count());
  for (std::size_t face = 0; face < view.face_count(); ++face)
    result[face] = view(face, 0);
  return result;
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void test::MaterialDensityTransportTestAccess::force_attempt_identity_wrap(
    FlowState &state) noexcept {
  state.impl_->attempt_identity = std::numeric_limits<std::uint64_t>::max();
}
void test::MaterialDensityTransportTestAccess::set_accepted_transport_value(
    FlowState &state, bool history, std::size_t field, std::size_t local_cell,
    double value) {
  if (field >= state.impl_->fields.transported_cell_fields.size() ||
      local_cell >= cell_count(state.impl_->layout.cell_interior_extent))
    throw runtime::Error("material test transport index is invalid");
  const auto extent = state.impl_->layout.cell_interior_extent;
  const int i =
      static_cast<int>(local_cell % static_cast<std::size_t>(extent.x));
  const std::size_t yz = local_cell / static_cast<std::size_t>(extent.x);
  const int j = static_cast<int>(yz % static_cast<std::size_t>(extent.y));
  const int k = static_cast<int>(yz / static_cast<std::size_t>(extent.y));
  auto &storage = history ? state.impl_->history : state.impl_->committed;
  auto view = storage.acquire_write<double>(
      state.impl_->access, kStatePhase, kStateActor,
      state.impl_->fields.transported_cell_fields[field]);
  view(i, j, k, 0) = value;
}
void test::MaterialDensityPisoTestAccess::force_state_diagnostic_identity(
    FlowState &state, std::uint64_t value) noexcept {
  state.impl_->diagnostic_mutation_identity = value;
}
std::uint64_t
test::MaterialDensityPisoTestAccess::state_diagnostic_identity(
    const FlowState &state) noexcept {
  return state.diagnostic_mutation_identity();
}

bool test::MaterialDensityPisoTestAccess::state_attempt_active(
    const FlowState &state) noexcept {
  return state.attempt_active();
}
std::uint64_t test::MaterialDensityPisoTestAccess::state_attempt_identity(
    const FlowState &state) noexcept {
  return state.attempt_identity();
}
std::uint64_t test::MaterialDensityPisoTestAccess::state_allocation_identity(
    const FlowState &state, FlowLayer selected) {
  const runtime::FieldStorage *storage{};
  switch (selected) {
  case FlowLayer::history:
    storage = &state.impl_->history;
    break;
  case FlowLayer::committed:
    storage = &state.impl_->committed;
    break;
  case FlowLayer::trial:
    storage = &state.impl_->trial;
    break;
  }
  if (storage == nullptr)
    throw runtime::Error("flow-state test layer is invalid");
  const auto view = storage->acquire_read<double>(
      state.impl_->access, kStatePhase, kStateActor,
      state.impl_->fields.density);
  return static_cast<std::uint64_t>(
      reinterpret_cast<std::uintptr_t>(&view(0, 0, 0, 0)));
}
void test::MaterialDensityPisoTestAccess::set_accepted_face_mass_flux(
    FlowState &state, std::size_t face, double value) {
  auto view = state.impl_->committed.acquire_face_write<double>(
      state.impl_->access, kStatePhase, kStateActor,
      state.impl_->fields.face_mass_flux);
  if (face >= view.face_count())
    throw runtime::Error("flow-state test face is out of range");
  view(face, 0) = value;
}
void test::MaterialDensityPisoTestAccess::force_state_metadata(
    FlowState &state, AcceptedStepMetadata metadata) noexcept {
  state.impl_->metadata = metadata;
}
#endif

void FlowState::seed_accepted_layers(const FlowLayerValues &history,
                                     const FlowLayerValues &committed) {
  if (impl_->attempt_active) {
    throw runtime::Error("cannot seed FlowState during an active attempt");
  }
  validate_layer_values(*impl_, history);
  validate_layer_values(*impl_, committed);
  if (impl_->diagnostic_mutation_identity ==
      std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error("FlowState diagnostic mutation identity would wrap");
  }
  ++impl_->diagnostic_mutation_identity;
  impl_->history.begin_rebuild();
  impl_->committed.begin_rebuild();
  impl_->trial.begin_rebuild();
  write_layer(*impl_, impl_->history, history);
  write_layer(*impl_, impl_->committed, committed);
  write_layer(*impl_, impl_->trial, committed);
}

FlowLayerValues FlowState::snapshot(FlowLayer selected) const {
  return read_layer(*impl_, layer(selected));
}

bool detail::FlowStateCheckpointAccess::live(
    const FlowState &state) noexcept {
  return static_cast<bool>(state.impl_);
}
const runtime::FieldRegistry &
detail::FlowStateCheckpointAccess::registry(const FlowState &state) {
  if (!state.impl_)
    throw runtime::Error("Checkpoint v2 FlowState has been moved from");
  return *state.impl_->registry;
}
runtime::FieldLayoutSet
detail::FlowStateCheckpointAccess::layout(const FlowState &state) {
  if (!state.impl_)
    throw runtime::Error("Checkpoint v2 FlowState has been moved from");
  return state.impl_->layout;
}
bool detail::FlowStateCheckpointAccess::attempt_active(
    const FlowState &state) noexcept {
  return state.impl_ && state.impl_->attempt_active;
}
bool detail::FlowStateCheckpointAccess::diagnostic_identity_can_advance(
    const FlowState &state) noexcept {
  return state.impl_ &&
         state.impl_->diagnostic_mutation_identity !=
             std::numeric_limits<std::uint64_t>::max();
}
bool detail::FlowStateCheckpointAccess::read_transaction_ready(
    FlowState &state) noexcept {
  if (!diagnostic_identity_can_advance(state))
    return false;
  std::array<runtime::FieldStorage *, 4> storages{
      &state.impl_->history, &state.impl_->committed, &state.impl_->trial,
      &state.impl_->rollback_snapshot};
  return runtime::FieldStorage::restart_v2_read_transactions_ready(
      storages.data(), storages.size());
}
void detail::FlowStateCheckpointAccess::enter_read_transaction(
    FlowState &state) {
  if (!diagnostic_identity_can_advance(state))
    throw runtime::Error("Checkpoint v2 diagnostic identity would wrap");
  std::array<runtime::FieldStorage *, 4> storages{
      &state.impl_->history, &state.impl_->committed, &state.impl_->trial,
      &state.impl_->rollback_snapshot};
  runtime::FieldStorage::begin_restart_v2_read_transactions(storages.data(),
                                                            storages.size());
  ++state.impl_->diagnostic_mutation_identity;
}
std::uint64_t detail::FlowStateCheckpointAccess::diagnostic_identity(
    const FlowState &state) noexcept {
  return state.impl_ ? state.impl_->diagnostic_mutation_identity : 0U;
}
FlowState detail::FlowStateCheckpointAccess::prepare_replacement(
    const FlowState &state, const FlowLayerValues &history,
    const FlowLayerValues &committed, AcceptedStepMetadata metadata) {
  if (!state.impl_)
    throw runtime::Error("Checkpoint v2 FlowState has been moved from");
  FlowState replacement = FlowState::create(
      *state.impl_->registry, state.impl_->layout, state.impl_->fields,
      metadata);
  replacement.seed_accepted_layers(history, committed);
  return replacement;
}
void detail::FlowStateCheckpointAccess::publish_replacement(
    FlowState &state, FlowState &&replacement) noexcept {
  const auto identity = state.impl_->diagnostic_mutation_identity;
  replacement.impl_->diagnostic_mutation_identity = identity;
  state.impl_ = std::move(replacement.impl_);
}

void FlowState::begin_attempt() {
  if (impl_->attempt_active) {
    throw runtime::Error("FlowState attempt is already active");
  }
  if (impl_->attempt_identity == std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error("FlowState attempt identity would wrap");
  }
  if (impl_->diagnostic_mutation_identity ==
      std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error("FlowState diagnostic mutation identity would wrap");
  }
  ++impl_->diagnostic_mutation_identity;
  ++impl_->attempt_identity;

  // Preserve the exact inactive trial bytes in a fourth preallocated buffer.
  // Both control blocks are advanced before rotation, so neither a view of
  // the former inactive trial nor a view of the reusable buffer can become
  // live again when storage roles are restored.
  impl_->trial.begin_rebuild();
  impl_->rollback_snapshot.begin_rebuild();
  runtime::FieldStorage reusable(std::move(impl_->rollback_snapshot));
  impl_->rollback_snapshot = std::move(impl_->trial);
  impl_->trial = std::move(reusable);
  impl_->trial.begin_rebuild();
  impl_->rollback_snapshot_valid = true;
  try {
    copy_layer(*impl_, impl_->committed, impl_->trial);
  } catch (...) {
    impl_->trial.begin_rebuild();
    runtime::FieldStorage failed_buffer(std::move(impl_->trial));
    impl_->trial = std::move(impl_->rollback_snapshot);
    impl_->rollback_snapshot = std::move(failed_buffer);
    impl_->rollback_snapshot_valid = false;
    throw;
  }
  impl_->attempt_active = true;
}

void FlowState::rollback_attempt() {
  if (!impl_->attempt_active) {
    throw runtime::Error("FlowState has no active attempt to roll back");
  }
  if (!impl_->rollback_snapshot_valid) {
    throw runtime::Error("FlowState rollback snapshot is unavailable");
  }
  // Invalidate views of the active candidate even when commit preparation
  // already invalidated them, then restore the exact pre-attempt storage.
  impl_->trial.begin_rebuild();
  runtime::FieldStorage rejected_candidate(std::move(impl_->trial));
  impl_->trial = std::move(impl_->rollback_snapshot);
  impl_->rollback_snapshot = std::move(rejected_candidate);
  impl_->rollback_snapshot_valid = false;
  impl_->commit_prepared = false;
  impl_->attempt_active = false;
}

void FlowState::prepare_commit_attempt(AcceptedStepMetadata accepted) {
  if (!impl_->attempt_active) {
    throw runtime::Error("FlowState has no active attempt to commit");
  }
  if (impl_->commit_prepared) {
    throw runtime::Error("FlowState commit is already prepared");
  }
  validate_metadata(accepted);
  if (accepted.step != impl_->metadata.step + 1U ||
      accepted.time_s != impl_->metadata.time_s + accepted.dt_s) {
    throw runtime::Error("FlowState commit metadata does not advance once");
  }
  impl_->prepared_metadata = accepted;
  impl_->commit_prepared = true;
  // These are the only potentially throwing operations in commit
  // preparation.  They run before the caller's collective decision and
  // invalidate every checked view before the no-throw storage rotation.
  impl_->trial.begin_rebuild();
  impl_->history.begin_rebuild();
  impl_->committed.begin_rebuild();
}

void FlowState::publish_commit_attempt() noexcept {
  // FieldStorage move construction/assignment is statically required to be
  // noexcept above.  Rotating the already allocated layers publishes trial
  // as committed and the old committed layer as history without allocation.
  runtime::FieldStorage recycled_history(std::move(impl_->history));
  impl_->history = std::move(impl_->committed);
  impl_->committed = std::move(impl_->trial);
  impl_->trial = std::move(recycled_history);
  impl_->metadata = impl_->prepared_metadata;
  impl_->rollback_snapshot_valid = false;
  impl_->commit_prepared = false;
  impl_->attempt_active = false;
}

void FlowState::commit_attempt(AcceptedStepMetadata accepted) {
  prepare_commit_attempt(accepted);
  publish_commit_attempt();
}

} // namespace hundun::flow
