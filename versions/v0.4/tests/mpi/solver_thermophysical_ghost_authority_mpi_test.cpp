// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_execution.hpp"
#include "hundun/v04_flow.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr FieldId kVelocity = 0U;
constexpr FieldId kPressure = 1U;
constexpr FieldId kEnthalpy = 2U;
constexpr FieldId kSpecies = 3U;
constexpr std::uint8_t kReach = 2U;
constexpr double kGhostSentinel = -777777.0;
constexpr double kInletEnthalpy = 700.0;
constexpr double kInletSpecies = 0.25;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local, MPI_Comm communicator = MPI_COMM_WORLD) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         output != 0;
}

std::uint64_t packed_status(Status status) noexcept {
  return (static_cast<std::uint64_t>(status.code) << 32U) |
         static_cast<std::uint64_t>(status.detail);
}

bool same_status(Status status, MPI_Comm communicator) {
  const std::uint64_t local = packed_status(status);
  std::uint64_t minimum = local;
  std::uint64_t maximum = local;
  return MPI_Allreduce(&local, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(&local, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {2.0, 8.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {2, 8, 2};
  mesh.minimum_spacing = {1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 32U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 20U;
  return mesh;
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec face;
  face.flow_kind = BoundaryKind::no_slip_wall;
  face.thermal_kind = BoundaryKind::adiabatic_wall;
  face.temperature = 300.0;
  face.mach_limit = 0.95;
  return face;
}

ValidatedModel authority_model() {
  ValidatedModel model;
  model.fingerprint = UINT64_C(0x544845524d4f4155);
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  model.mesh = mesh_spec();
  model.transported_scalars.push_back(
      {"Y", TransportedScalarRole::species, 1.0, 1.0});
  for (BoundaryFaceSpec& face : model.boundaries) {
    face = wall();
  }

  // Rank decomposition for 2 ranks is 1x2x1 for this anisotropic mesh. The
  // inlet therefore meets a true MPI y interface in the local x-y corner.
  BoundaryFaceSpec& x_min = model.boundaries[0U];
  x_min.flow_kind = BoundaryKind::velocity_inlet;
  x_min.thermal_kind = BoundaryKind::none;
  x_min.velocity = {1.0, 0.0, 0.0};
  x_min.temperature = 300.0;
  x_min.scalars.clear();
  x_min.scalars.push_back(
      {"Y", ScalarBoundaryKind::dirichlet, kInletSpecies,
       ScalarBoundaryKind::zero_gradient, 0.0});

  BoundaryFaceSpec& x_max = model.boundaries[1U];
  x_max.flow_kind = BoundaryKind::pressure_outlet;
  x_max.thermal_kind = BoundaryKind::none;
  x_max.pressure = 101325.0;

  // y remains physical so only the internal y interface is supplied by MPI.
  // z is periodic to keep the fixture focused on the x-inlet/y-MPI seam.
  model.boundaries[4U].flow_kind = BoundaryKind::periodic;
  model.boundaries[4U].thermal_kind = BoundaryKind::none;
  model.boundaries[5U].flow_kind = BoundaryKind::periodic;
  model.boundaries[5U].thermal_kind = BoundaryKind::none;
  for (std::size_t index = 0U; index < model.boundaries.size(); ++index) {
    if (model.boundaries[index].flow_kind != BoundaryKind::periodic &&
        index != 0U) {
      model.boundaries[index].scalars.push_back(
          {"Y", ScalarBoundaryKind::zero_gradient, 0.0,
           ScalarBoundaryKind::zero_gradient, 0.0});
    }
  }
  return model;
}

struct GhostBundle {
  ThermophysicalGhostHistory enthalpy{};
  std::array<ThermophysicalGhostHistory, 1U> species{};
};

ThermophysicalGhostAuthority authority_for(ConstFieldView view,
                                            const CartesianGeometryPlan& geometry,
                                            const BoundaryPlan& boundary,
                                            const HaloEngine& halo) noexcept {
  return {halo.instance_identity(),
          view.field,
          view.revision,
          view.storage_identity,
          view.revision_domain,
          geometry.topology_revision(),
          boundary.revision(),
          kReach};
}

bool authority_matches(ThermophysicalGhostAuthority authority,
                       ConstFieldView view,
                       const CartesianGeometryPlan& geometry,
                       const BoundaryPlan& boundary) noexcept {
  return authority.valid() && authority.field == view.field &&
         authority.state == view.revision &&
         authority.storage == view.storage_identity &&
         authority.revision_domain == view.revision_domain &&
         authority.geometry == geometry.topology_revision() &&
         authority.boundary == boundary.revision() && authority.reach >= kReach;
}

void clear_padded(FieldView view, double value) {
  for (std::int32_t z = -view.ghosts.z;
       z < view.interior.z + view.ghosts.z; ++z) {
    for (std::int32_t y = -view.ghosts.y;
         y < view.interior.y + view.ghosts.y; ++y) {
      for (std::int32_t x = -view.ghosts.x;
           x < view.interior.x + view.ghosts.x; ++x) {
        view.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
}

void fill_interior(FieldView view, double value) {
  for (std::int32_t z = 0; z < view.interior.z; ++z) {
    for (std::int32_t y = 0; y < view.interior.y; ++y) {
      for (std::int32_t x = 0; x < view.interior.x; ++x) {
        view.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
}

bool make_plans(MPI_Comm communicator, CartesianGeometryPlan& geometry,
                MeshPatch& patch, FieldRegistry& registry,
                BoundaryPlan& boundary, SchemePlan& schemes,
                TimeSchemePlan& time, FieldSchema& schema, FieldId& enthalpy,
                FieldId& species, int rank) {
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          communicator, mesh_spec(), GeometryBudget{0U, 1U}, geometry,
          patch)),
      rank, "anisotropic authority geometry compiles");
  if (passed) {
    passed &= expect(patch.process_grid.x == 1 && patch.process_grid.y == 2 &&
                         patch.process_grid.z == 1,
                     rank, "2-rank decomposition is along y");
  }
  ValidatedModel model = authority_model();
  if (passed) {
    const Status boundary_status = BoundaryCompiler::compile(
        communicator, model, geometry, patch, registry, boundary, schemes,
        time);
    if (!boundary_status) {
      std::cerr << "rank " << rank << " boundary compile status="
                << static_cast<unsigned>(boundary_status.code) << "/"
                << boundary_status.detail << '\n';
    }
    passed &= expect(static_cast<bool>(boundary_status),
                     rank, "x-inlet/y-MPI boundary plan compiles");
  }
  if (passed) {
    passed &= expect(static_cast<bool>(registry.freeze(schema)), rank,
                     "authority field schema freezes");
  }
  if (passed) {
    for (const FieldDescriptor& field : schema) {
      if (field.stable_name == "h") enthalpy = field.id;
      if (field.stable_name == "Y") species = field.id;
    }
    passed &= expect(enthalpy == kEnthalpy && species == kSpecies, rank,
                     "boundary registry preserves h/Y field identities");
  }
  return all_true(passed, communicator);
}

bool make_layers(const FieldSchema& schema, MeshPatch patch, StateLayers& layers,
                 int rank) {
  const std::array<ArenaFieldRequest, 4U> requests{{
      {kVelocity, patch.cells, {0U}, FieldLifetime::state_layer},
      {kPressure, patch.cells, {0U}, FieldLifetime::state_layer},
      {kEnthalpy, patch.cells, {0U}, FieldLifetime::state_layer},
      {kSpecies, patch.cells, {0U}, FieldLifetime::state_layer},
  }};
  ArenaLayout layout;
  const bool passed = expect(
      static_cast<bool>(ArenaLayout::compile(
          schema, {requests.data(), requests.size()}, layout)) &&
          static_cast<bool>(StateLayers::allocate(layout, layers)),
      rank, "state-layer arena allocates h/Y with two ghosts");
  return all_true(passed);
}

bool exchange_state(HaloEngine& halo, StageId stage,
                    std::array<FieldView, 2U>& fields, int rank) {
  HaloTicket ticket;
  const Status begun = halo.begin(
      stage, {fields.data(), fields.size()}, ticket);
  bool passed = expect(static_cast<bool>(begun) && ticket.active(), rank,
                       "state ghost exchange begins");
  const Status finished =
      begun ? halo.finish(ticket, {fields.data(), fields.size()}) : begun;
  passed &= expect(static_cast<bool>(finished) && !ticket.active(), rank,
                   "state ghost exchange finishes");
  passed &= expect(halo.ghost_revision(kEnthalpy) == fields[0U].revision &&
                       halo.ghost_revision(kSpecies) == fields[1U].revision,
                   rank, "halo publishes exact state revisions");
  return all_true(passed);
}

bool apply_physical_boundaries(const BoundaryPlan& boundary,
                               FieldView enthalpy, FieldView species,
                               double inlet_enthalpy, int rank) {
  std::vector<double> resolved(boundary.resolved_scalar_count(),
                               inlet_enthalpy);
  const BoundaryResolvedValues values{
      {resolved.data(), resolved.size()}, {}, {}};
  bool passed = expect(
      static_cast<bool>(apply_boundary_ghosts(
          BoundaryStage::enthalpy, boundary, {&enthalpy, 1U}, values)),
      rank, "x-inlet enthalpy mirror applies after MPI exchange");
  passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                       BoundaryStage::scalar, boundary, {&species, 1U}, {})),
                   rank, "x-inlet scalar mirror applies after MPI exchange");
  return all_true(passed);
}

bool verify_y_interface(ConstFieldView view, int rank, double rank_zero_label,
                        double rank_one_label, std::string_view description) {
  const double expected = rank == 0 ? rank_one_label : rank_zero_label;
  const std::int32_t y = rank == 0 ? view.interior.y : 0;
  bool passed = true;
  for (std::int32_t z = 0; z < view.interior.z; ++z) {
    for (std::int32_t x = 0; x < view.interior.x; ++x) {
      for (std::int32_t layer = 1; layer <= kReach; ++layer) {
        const std::int32_t ghost_y = rank == 0 ? y - 1 + layer : -layer;
        passed &= expect(view.unchecked({x, ghost_y, z}, 0U) == expected,
                         rank, description);
      }
    }
  }
  return passed;
}

bool verify_x_inlet(ConstFieldView view, int rank, double interior_label,
                    double face_value, std::string_view description) {
  bool passed = true;
  for (std::int32_t z = 0; z < view.interior.z; ++z) {
    for (std::int32_t y = 0; y < view.interior.y; ++y) {
      for (std::int32_t layer = 1; layer <= kReach; ++layer) {
        const double expected =
            2.0 * face_value - (interior_label + static_cast<double>(rank));
        passed &= expect(view.unchecked({-layer, y, z}, 0U) == expected,
                         rank, description);
      }
    }
  }
  return passed;
}

bool verify_identity(ThermophysicalGhostAuthority authority,
                     ThermophysicalGhostAuthority expected, int rank,
                     std::string_view description) {
  return expect(authority.exchange_plan == expected.exchange_plan &&
                    authority.field == expected.field &&
                    authority.state == expected.state &&
                    authority.storage == expected.storage &&
                    authority.revision_domain == expected.revision_domain &&
                    authority.geometry == expected.geometry &&
                    authority.boundary == expected.boundary &&
                    authority.reach == expected.reach,
                rank, description);
}

bool run_authority_case(int rank, int size) {
  if (size != 2) {
    return expect(false, rank, "authority test requires exactly two ranks");
  }

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  FieldRegistry registry;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  FieldSchema schema;
  FieldId enthalpy = 0U;
  FieldId species = 0U;
  bool passed = make_plans(MPI_COMM_WORLD, geometry, patch, registry, boundary,
                           schemes, time, schema, enthalpy, species, rank);
  StateLayers layers;
  if (passed) passed = make_layers(schema, patch, layers, rank);

  AttemptTransaction transaction;
  if (passed) {
    passed &= expect(static_cast<bool>(AttemptTransaction::create(
                             layers.field_count(), 1U,
                             layers.field_count(), transaction)),
                     rank, "authority transaction allocates");
  }

  const std::array<HaloFieldSpec, 2U> halo_specs{{
      {enthalpy, kReach, 1U},
      {species, kReach, 1U},
  }};
  HaloEngine halo;
  if (passed) {
    passed &= expect(static_cast<bool>(halo.reserve(
                             MPI_COMM_WORLD, patch,
                             {halo_specs.data(), halo_specs.size()},
                             boundary.halo_topology())),
                     rank, "production HaloEngine reserves h/Y exchange");
  }
  passed = all_true(passed);
  if (!passed) return false;

  FieldView previous_h;
  FieldView previous_y;
  FieldView accepted_h;
  FieldView accepted_y;
  passed &= expect(static_cast<bool>(layers.view(
                           StateRole::accepted_n_minus_one, enthalpy,
                           previous_h)) &&
                       static_cast<bool>(layers.view(
                           StateRole::accepted_n_minus_one, species,
                           previous_y)) &&
                       static_cast<bool>(layers.view(StateRole::accepted_n,
                                                     enthalpy, accepted_h)) &&
                       static_cast<bool>(layers.view(StateRole::accepted_n,
                                                     species, accepted_y)),
                   rank, "accepted/previous state views bind");
  if (!all_true(passed)) return false;

  clear_padded(previous_h, kGhostSentinel);
  clear_padded(previous_y, kGhostSentinel);
  clear_padded(accepted_h, kGhostSentinel);
  clear_padded(accepted_y, kGhostSentinel);
  fill_interior(previous_h, 2000.0 + rank);
  fill_interior(previous_y, 4000.0 + rank);
  fill_interior(accepted_h, 1000.0 + rank);
  fill_interior(accepted_y, 3000.0 + rank);

  std::array<FieldView, 2U> previous_fields{previous_h, previous_y};
  std::array<FieldView, 2U> accepted_fields{accepted_h, accepted_y};
  passed &= exchange_state(halo, 10U, previous_fields, rank);
  passed &= apply_physical_boundaries(boundary, previous_h, previous_y, 800.0,
                                      rank);
  passed &= exchange_state(halo, 11U, accepted_fields, rank);
  passed &= apply_physical_boundaries(boundary, accepted_h, accepted_y,
                                      kInletEnthalpy, rank);

  GhostBundle ghosts;
  ghosts.enthalpy.previous = authority_for(
      as_const(previous_h), geometry, boundary, halo);
  ghosts.enthalpy.accepted =
      authority_for(as_const(accepted_h), geometry, boundary, halo);
  ghosts.species[0U].previous = authority_for(
      as_const(previous_y), geometry, boundary, halo);
  ghosts.species[0U].accepted =
      authority_for(as_const(accepted_y), geometry, boundary, halo);
  passed &= expect(ghosts.enthalpy.previous.exchange_plan ==
                       ghosts.enthalpy.accepted.exchange_plan &&
                       ghosts.species[0U].previous.exchange_plan ==
                           ghosts.enthalpy.accepted.exchange_plan,
                   rank, "accepted/previous authorities share one halo plan");
  passed &= verify_y_interface(as_const(previous_h), rank, 2000.0, 2001.0,
                              "previous h has both MPI y ghost layers");
  passed &= verify_y_interface(as_const(accepted_h), rank, 1000.0, 1001.0,
                              "accepted h has both MPI y ghost layers");
  passed &= verify_y_interface(as_const(previous_y), rank, 4000.0, 4001.0,
                              "previous Y has both MPI y ghost layers");
  passed &= verify_y_interface(as_const(accepted_y), rank, 3000.0, 3001.0,
                              "accepted Y has both MPI y ghost layers");
  passed &= verify_x_inlet(as_const(previous_h), rank, 2000.0, 800.0,
                           "previous h has two-layer x-inlet mirror");
  passed &= verify_x_inlet(as_const(accepted_h), rank, 1000.0,
                           kInletEnthalpy,
                           "accepted h has two-layer x-inlet mirror");
  passed &= verify_x_inlet(as_const(previous_y), rank, 4000.0,
                           kInletSpecies,
                           "previous Y has two-layer x-inlet mirror");
  passed &= verify_x_inlet(as_const(accepted_y), rank, 3000.0,
                           kInletSpecies,
                           "accepted Y has two-layer x-inlet mirror");
  passed = all_true(passed);
  if (!passed) return false;

  const std::array<std::size_t, 3U> old_handles{
      layers.handle(StateRole::accepted_n),
      layers.handle(StateRole::accepted_n_minus_one),
      layers.handle(StateRole::trial)};

  // A rank-one stale previous authority must become one collective failure;
  // no role handle or published authority is allowed to rotate.
  passed &= expect(static_cast<bool>(transaction.begin(layers)), rank,
                   "stale-authority attempt begins");
  for (FieldId field = 0U; field < layers.field_count() && passed; ++field) {
    passed &= expect(static_cast<bool>(transaction.revise_trial(field)), rank,
                     "stale-authority attempt revises every state field");
  }
  ThermophysicalGhostHistory stale_enthalpy = ghosts.enthalpy;
  if (rank == 1) ++stale_enthalpy.previous.state;
  const bool authority_ok =
      authority_matches(stale_enthalpy.accepted, as_const(accepted_h),
                        geometry, boundary) &&
      authority_matches(stale_enthalpy.previous, as_const(previous_h),
                        geometry, boundary) &&
      authority_matches(ghosts.species[0U].accepted, as_const(accepted_y),
                        geometry, boundary) &&
      authority_matches(ghosts.species[0U].previous, as_const(previous_y),
                        geometry, boundary);
  const Status stale_status = transaction.collective_finish(
      MPI_COMM_WORLD,
      authority_ok ? Status{} : Status{StatusCode::invalid_plan, 9201U});
  passed &= expect(stale_status.code == StatusCode::invalid_plan &&
                       transaction.lowest_failing_rank() == 1 &&
                       !transaction.committed(),
                   rank, "stale authority is rejected collectively");
  passed &= expect(same_status(stale_status, MPI_COMM_WORLD), rank,
                   "stale authority status is identical on both ranks");
  passed &= expect(layers.handle(StateRole::accepted_n) == old_handles[0U] &&
                       layers.handle(StateRole::accepted_n_minus_one) ==
                           old_handles[1U] &&
                       layers.handle(StateRole::trial) == old_handles[2U],
                   rank, "stale authority rollback preserves role handles");
  passed = all_true(passed);
  if (!passed) return false;

  // A second rank-one fault uses a different exchange-plan token. This is the
  // cross-history mismatch rejected by the production predictor contract;
  // it must also be collective and leave all role handles untouched.
  passed &= expect(static_cast<bool>(transaction.begin(layers)), rank,
                   "mismatched-authority attempt begins");
  for (FieldId field = 0U; field < layers.field_count() && passed; ++field) {
    passed &= expect(static_cast<bool>(transaction.revise_trial(field)), rank,
                     "mismatched-authority attempt revises every state field");
  }
  ThermophysicalGhostHistory mismatched_enthalpy = ghosts.enthalpy;
  if (rank == 1) ++mismatched_enthalpy.previous.exchange_plan;
  const bool mismatch_authority_ok =
      authority_matches(mismatched_enthalpy.accepted, as_const(accepted_h),
                         geometry, boundary) &&
      authority_matches(mismatched_enthalpy.previous, as_const(previous_h),
                         geometry, boundary) &&
      mismatched_enthalpy.previous.exchange_plan ==
          mismatched_enthalpy.accepted.exchange_plan &&
      authority_matches(ghosts.species[0U].accepted, as_const(accepted_y),
                         geometry, boundary) &&
      authority_matches(ghosts.species[0U].previous, as_const(previous_y),
                         geometry, boundary) &&
      ghosts.species[0U].accepted.exchange_plan ==
          mismatched_enthalpy.accepted.exchange_plan &&
      ghosts.species[0U].previous.exchange_plan ==
          mismatched_enthalpy.accepted.exchange_plan;
  const Status mismatch_status = transaction.collective_finish(
      MPI_COMM_WORLD,
      mismatch_authority_ok ? Status{}
                            : Status{StatusCode::invalid_plan, 9202U});
  passed &= expect(mismatch_status.code == StatusCode::invalid_plan &&
                       transaction.lowest_failing_rank() == 1 &&
                       !transaction.committed(),
                   rank, "mismatched ghost authority is rejected collectively");
  passed &= expect(same_status(mismatch_status, MPI_COMM_WORLD), rank,
                   "mismatched authority status is identical on both ranks");
  passed &= expect(layers.handle(StateRole::accepted_n) == old_handles[0U] &&
                       layers.handle(StateRole::accepted_n_minus_one) ==
                           old_handles[1U] &&
                       layers.handle(StateRole::trial) == old_handles[2U],
                   rank, "mismatched authority rollback preserves role handles");
  passed = all_true(passed);
  if (!passed) return false;

  // A fresh attempt exchanges the new trial state, then commits. The explicit
  // authority co-rotation mirrors the product driver's published seam: old
  // accepted becomes previous and trial becomes accepted only after commit.
  passed &= expect(static_cast<bool>(transaction.begin(layers)), rank,
                   "co-rotation attempt begins after rollback");
  for (FieldId field = 0U; field < layers.field_count() && passed; ++field) {
    passed &= expect(static_cast<bool>(transaction.revise_trial(field)), rank,
                     "co-rotation attempt revises every state field");
  }
  FieldView trial_h;
  FieldView trial_y;
  passed &= expect(static_cast<bool>(layers.view(StateRole::trial, enthalpy,
                                                  trial_h)) &&
                       static_cast<bool>(layers.view(StateRole::trial, species,
                                                     trial_y)),
                   rank, "trial state views bind");
  clear_padded(trial_h, kGhostSentinel);
  clear_padded(trial_y, kGhostSentinel);
  fill_interior(trial_h, 5000.0 + rank);
  fill_interior(trial_y, 6000.0 + rank);
  std::array<FieldView, 2U> trial_fields{trial_h, trial_y};
  passed &= exchange_state(halo, 12U, trial_fields, rank);
  passed &= apply_physical_boundaries(boundary, trial_h, trial_y,
                                      kInletEnthalpy, rank);
  const ThermophysicalGhostAuthority trial_h_authority = authority_for(
      as_const(trial_h), geometry, boundary, halo);
  const ThermophysicalGhostAuthority trial_y_authority = authority_for(
      as_const(trial_y), geometry, boundary, halo);
  const Status committed = transaction.collective_finish(MPI_COMM_WORLD,
                                                          Status{});
  passed &= expect(static_cast<bool>(committed) && transaction.committed(),
                   rank, "valid trial commits collectively");
  passed &= expect(same_status(committed, MPI_COMM_WORLD), rank,
                   "successful commit status is identical on both ranks");
  if (committed) {
    ghosts.enthalpy.previous = ghosts.enthalpy.accepted;
    ghosts.enthalpy.accepted = trial_h_authority;
    ghosts.species[0U].previous = ghosts.species[0U].accepted;
    ghosts.species[0U].accepted = trial_y_authority;
  }
  FieldView committed_h;
  FieldView committed_previous_h;
  FieldView committed_y;
  FieldView committed_previous_y;
  passed &= expect(static_cast<bool>(layers.view(StateRole::accepted_n,
                                                 enthalpy, committed_h)) &&
                       static_cast<bool>(layers.view(
                           StateRole::accepted_n_minus_one, enthalpy,
                           committed_previous_h)) &&
                       static_cast<bool>(layers.view(StateRole::accepted_n,
                                                     species, committed_y)) &&
                       static_cast<bool>(layers.view(
                           StateRole::accepted_n_minus_one, species,
                           committed_previous_y)),
                   rank, "committed state views rebind after rotation");
  passed &= verify_identity(ghosts.enthalpy.accepted, trial_h_authority, rank,
                            "accepted authority co-rotates to committed trial");
  passed &= verify_identity(ghosts.enthalpy.previous,
                            authority_for(as_const(committed_previous_h),
                                          geometry, boundary, halo),
                            rank,
                            "previous authority co-rotates to old accepted");
  passed &= verify_identity(ghosts.species[0U].accepted, trial_y_authority,
                            rank, "species accepted authority co-rotates");
  passed &= verify_identity(ghosts.species[0U].previous,
                            authority_for(as_const(committed_previous_y),
                                          geometry, boundary, halo),
                            rank, "species previous authority co-rotates");
  passed &= verify_y_interface(as_const(committed_h), rank, 5000.0, 5001.0,
                              "committed accepted h retains trial MPI ghosts");
  passed &= verify_y_interface(as_const(committed_previous_h), rank, 1000.0,
                              1001.0,
                              "committed previous h retains old MPI ghosts");
  passed &= verify_y_interface(as_const(committed_y), rank, 6000.0, 6001.0,
                              "committed accepted Y retains trial MPI ghosts");
  passed &= verify_y_interface(as_const(committed_previous_y), rank, 3000.0,
                              3001.0,
                              "committed previous Y retains old MPI ghosts");
  passed &= verify_x_inlet(as_const(committed_h), rank, 5000.0,
                           kInletEnthalpy,
                           "committed accepted h retains inlet mirror");
  passed &= verify_x_inlet(as_const(committed_previous_h), rank, 1000.0,
                           kInletEnthalpy,
                           "committed previous h retains inlet mirror");
  passed &= verify_x_inlet(as_const(committed_y), rank, 6000.0,
                           kInletSpecies,
                           "committed accepted Y retains inlet mirror");
  passed &= verify_x_inlet(as_const(committed_previous_y), rank, 3000.0,
                           kInletSpecies,
                           "committed previous Y retains inlet mirror");
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const bool passed = run_authority_case(rank, size);
  if (rank == 0 && passed) {
    std::cout << "v0.4 thermophysical ghost authority MPI test passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
