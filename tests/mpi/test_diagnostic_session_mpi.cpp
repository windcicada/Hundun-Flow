// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/diagnostics/diagnostic_session.hpp"
#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"
#include "hundun/diagnostics/stage2_module_diagnostics.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/matrix_free_poisson.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/linear/ghosted_vector_halo.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

class CaptureSink final : public hundun::diagnostics::DiagnosticSink {
 public:
  void submit(const hundun::diagnostics::DiagnosticRecord& record) override {
    records.push_back(record);
  }
  std::vector<hundun::diagnostics::DiagnosticRecord> records;
};

class RankFailingSink final : public hundun::diagnostics::DiagnosticSink {
 public:
  RankFailingSink(int rank, int failing_rank)
      : rank_(rank), failing_rank_(failing_rank) {}
  void submit(const hundun::diagnostics::DiagnosticRecord&) override {
    if (rank_ == failing_rank_)
      throw std::runtime_error("injected rank-local sink failure");
  }

 private:
  int rank_{};
  int failing_rank_{};
};

template <class Source>
void check_fingerprint_ids(const Source &source,
                           std::initializer_list<std::string_view> expected) {
  const auto actual =
      hundun::diagnostics::diagnostic_fingerprint_field_ids(source);
  HUNDUN_CHECK(std::is_sorted(actual.begin(), actual.end()));
  HUNDUN_CHECK(std::adjacent_find(actual.begin(), actual.end()) ==
               actual.end());
  HUNDUN_CHECK(std::vector<std::string_view>(expected) == actual);
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back({"alpha", 0.0});
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t i = 0; i < names.size(); ++i) {
    config.boundaries[i].patch = names[i];
    config.boundaries[i].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::diagnostics::DiagnosticBatch batch_for_rank(int rank) {
  using namespace hundun::diagnostics;
  const DiagnosticDescriptor descriptor{
      kDiagnosticRecordSchemaV1, DiagnosticModuleKind::flow_driver,
      "hundun.application.flow_driver", "flow-driver",
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary)};
  const DiagnosticRequest request{
      DiagnosticLevel::summary, DiagnosticScope::local,
      {rank, 2U, 0.25, "accepted-step"}, {}, 0U};
  DiagnosticRecord record;
  record.module_kind = descriptor.module_kind;
  record.module_id = std::string(descriptor.module_id);
  record.instance_id = std::string(descriptor.instance_id);
  record.rank = rank;
  record.step = 2U;
  record.time_s = describe_fp64(0.25);
  record.phase = "accepted-step";
  record.metrics.push_back(
      {"attempt_count", DiagnosticMetricKind::state_summary, "count",
       describe_fp64(1.0)});
  DiagnosticFingerprintAccumulator fingerprint;
  fingerprint.add("attempt_count", 0U, 0U, describe_fp64(1.0));
  record.state_fingerprint = fingerprint.finish();
  DiagnosticBatch result;
  DiagnosticBatchSink sink(descriptor, request, result);
  sink.submit(record);
  return result;
}

std::filesystem::path rank_file(const std::filesystem::path& directory,
                                int rank) {
  std::ostringstream name;
  name << "diagnostics.v1.rank-" << std::setw(6) << std::setfill('0') << rank
       << ".step-00000000000000000002.jsonl";
  return directory / name.str();
}

std::filesystem::path mesh_file(const std::filesystem::path& directory,
                                int rank) {
  std::ostringstream name;
  name << "meshdiag.v2.rank-" << std::setw(6) << std::setfill('0') << rank
       << ".bin";
  return directory / name.str();
}

std::vector<char> read_bytes(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  HUNDUN_CHECK(static_cast<bool>(stream));
  return {std::istreambuf_iterator<char>(stream),
          std::istreambuf_iterator<char>()};
}

struct IncidentAreaEvidence {
  std::uint64_t remote_owner_face_count{};
  std::optional<std::uint64_t> mutation_cell;
  double old_owned_face_area_sum{};
  double complete_incident_area_sum{};
};

IncidentAreaEvidence check_incident_face_areas(
    const hundun::mesh::MeshTopology& topology,
    const hundun::mesh::MeshGeometry& geometry,
    const hundun::diagnostics::MeshDiagnosticV2& decoded) {
  std::vector<hundun::mesh::LocalFaceId> incident_faces;
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    if (topology.cell_ownership(owner) ==
            hundun::mesh::EntityOwnership::owned ||
        (neighbour &&
         topology.cell_ownership(*neighbour) ==
             hundun::mesh::EntityOwnership::owned)) {
      incident_faces.push_back(face);
    }
  }
  std::sort(incident_faces.begin(), incident_faces.end(),
            [&](const auto left, const auto right) {
              return topology.global_face_id(left) <
                     topology.global_face_id(right);
            });
  HUNDUN_CHECK(incident_faces.size() == decoded.faces.size());

  std::map<std::uint64_t, double> expected_area_sums;
  std::map<std::uint64_t, double> decoded_area_sums;
  for (hundun::mesh::LocalCellId cell = 0;
       cell < topology.owned_cell_count(); ++cell) {
    const auto global_id = topology.global_cell_id(cell);
    expected_area_sums.emplace(global_id, 0.0);
    decoded_area_sums.emplace(global_id, 0.0);
  }

  IncidentAreaEvidence evidence;
  for (std::size_t index = 0; index < incident_faces.size(); ++index) {
    const auto local_face = incident_faces[index];
    const auto owner = topology.owner(local_face);
    const auto neighbour = topology.neighbour(local_face);
    const auto owner_global = topology.global_cell_id(owner);
    const auto& actual = decoded.faces[index];
    HUNDUN_CHECK(actual.global_id == topology.global_face_id(local_face));
    HUNDUN_CHECK(actual.owner_global_cell == owner_global);
    HUNDUN_CHECK(
        actual.neighbour_global_cell ==
        (neighbour
             ? std::optional<std::uint64_t>{
                   topology.global_cell_id(*neighbour)}
             : std::nullopt));
    HUNDUN_CHECK(actual.area_m2 == geometry.face_area_m2(local_face));

    if (const auto found = expected_area_sums.find(owner_global);
        found != expected_area_sums.end())
      found->second += geometry.face_area_m2(local_face);
    if (neighbour) {
      if (const auto found =
              expected_area_sums.find(topology.global_cell_id(*neighbour));
          found != expected_area_sums.end())
        found->second += geometry.face_area_m2(local_face);
    }
    if (const auto found =
            decoded_area_sums.find(actual.owner_global_cell);
        found != decoded_area_sums.end())
      found->second += actual.area_m2;
    if (actual.neighbour_global_cell) {
      if (const auto found =
              decoded_area_sums.find(*actual.neighbour_global_cell);
          found != decoded_area_sums.end())
        found->second += actual.area_m2;
    }

    if (neighbour &&
        topology.cell_ownership(owner) ==
            hundun::mesh::EntityOwnership::ghost &&
        topology.cell_ownership(*neighbour) ==
            hundun::mesh::EntityOwnership::owned) {
      ++evidence.remote_owner_face_count;
      if (!evidence.mutation_cell)
        evidence.mutation_cell = topology.global_cell_id(*neighbour);
    }
  }
  HUNDUN_CHECK(expected_area_sums == decoded_area_sums);

  if (evidence.mutation_cell) {
    for (const auto local_face : incident_faces) {
      const auto owner = topology.owner(local_face);
      const auto neighbour = topology.neighbour(local_face);
      const bool incident =
          topology.global_cell_id(owner) == *evidence.mutation_cell ||
          (neighbour &&
           topology.global_cell_id(*neighbour) ==
               *evidence.mutation_cell);
      if (!incident)
        continue;
      evidence.complete_incident_area_sum +=
          geometry.face_area_m2(local_face);
      if (topology.cell_ownership(owner) ==
          hundun::mesh::EntityOwnership::owned) {
        evidence.old_owned_face_area_sum +=
            geometry.face_area_m2(local_face);
      }
    }
    HUNDUN_CHECK(evidence.complete_incident_area_sum >
                 evidence.old_owned_face_area_sum);
  }
  return evidence;
}

}  // namespace

int main(int argc, char** argv) {
  return hundun::test::run([&] {
    hundun::runtime::MpiEnvironment environment(argc, argv);
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    check_fingerprint_ids(mpi,
                          {"collective_calls", "logical_payload_bytes", "rank",
                           "reduced_scalars", "size", "thread_level"});
    const auto counters_before = mpi.fp64_reduction_counters();
    const hundun::diagnostics::DiagnosticRequest local_request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local,
        {mpi.rank(), 2U, 0.25, "accepted-step"},
        {},
        0U};
    CaptureSink local_sink;
    hundun::diagnostics::collect_diagnostics(
        mpi, local_request, local_sink);
    HUNDUN_CHECK(local_sink.records.size() == 1U);
    const auto counters_after_local = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(counters_after_local.collective_calls ==
                 counters_before.collective_calls);
    HUNDUN_CHECK(counters_after_local.reduced_scalars ==
                 counters_before.reduced_scalars);
    HUNDUN_CHECK(counters_after_local.logical_payload_bytes ==
                 counters_before.logical_payload_bytes);

    auto collective_request = local_request;
    collective_request.scope =
        hundun::diagnostics::DiagnosticScope::collective;
    CaptureSink collective_sink;
    hundun::diagnostics::collect_diagnostics(
        mpi, mpi, collective_request, collective_sink);
    HUNDUN_CHECK(collective_sink.records.size() == 1U);
    auto root_fingerprint =
        collective_sink.records.front().state_fingerprint.hex;
    std::uint64_t root_fingerprint_size =
        mpi.rank() == 0
            ? static_cast<std::uint64_t>(root_fingerprint.size())
            : 0U;
    HUNDUN_CHECK(MPI_Bcast(&root_fingerprint_size, 1, MPI_UINT64_T, 0,
                           mpi.comm()) == MPI_SUCCESS);
    root_fingerprint.resize(
        static_cast<std::size_t>(root_fingerprint_size));
    HUNDUN_CHECK(MPI_Bcast(root_fingerprint.data(),
                           static_cast<int>(root_fingerprint.size()),
                           MPI_CHAR, 0, mpi.comm()) == MPI_SUCCESS);
    HUNDUN_CHECK(
        collective_sink.records.front().state_fingerprint.hex ==
        root_fingerprint);
    const auto counters_after_collective = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(counters_after_collective.collective_calls ==
                 counters_before.collective_calls);
    HUNDUN_CHECK(counters_after_collective.reduced_scalars ==
                 counters_before.reduced_scalars);
    HUNDUN_CHECK(counters_after_collective.logical_payload_bytes ==
                 counters_before.logical_payload_bytes);
    if (mpi.size() > 1) {
      const auto expect_collective_rejection =
          [&](hundun::diagnostics::DiagnosticRequest mismatched,
              int expected_lowest) {
            CaptureSink sink;
            bool rejected = false;
            int lowest = -1;
            try {
              hundun::diagnostics::collect_diagnostics(
                  mpi, mpi, mismatched, sink);
            } catch (
                const hundun::diagnostics::DiagnosticCollectionError& error) {
              rejected = true;
              lowest = error.lowest_failing_rank();
            }
            HUNDUN_CHECK(rejected);
            HUNDUN_CHECK(lowest == expected_lowest);
          };
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.rank = 1;
        expect_collective_rejection(std::move(mismatched), 0);
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 1)
          mismatched.level = hundun::diagnostics::DiagnosticLevel::counters;
        expect_collective_rejection(std::move(mismatched), 1);
      }
      {
        auto mismatched = collective_request;
        if (mpi.size() == 4 && mpi.rank() >= 2)
          mismatched.frame.step += 1U;
        if (mpi.size() == 4)
          expect_collective_rejection(std::move(mismatched), 2);
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.time_s = 0.5;
        expect_collective_rejection(std::move(mismatched), 1);
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.phase = "different-phase";
        expect_collective_rejection(std::move(mismatched), 1);
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 1)
          mismatched.selected_fields = {"rank"};
        expect_collective_rejection(std::move(mismatched), 1);
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 1)
          mismatched.sample_budget = 1U;
        expect_collective_rejection(std::move(mismatched), 1);
      }

      RankFailingSink failing_sink(mpi.rank(), 0);
      bool rejected = false;
      int lowest = -1;
      try {
        hundun::diagnostics::collect_diagnostics(
            mpi, mpi, collective_request, failing_sink);
      } catch (
          const hundun::diagnostics::DiagnosticCollectionError& error) {
        rejected =
            error.classification() ==
            hundun::diagnostics::DiagnosticFailureClass::sink_failure;
        lowest = error.lowest_failing_rank();
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(lowest == 0);
    }

    const auto root = std::filesystem::temp_directory_path() /
                      ("hundun-task24-session-mpi-" +
                       std::to_string(mpi.size()));
    if (mpi.rank() == 0)
      std::filesystem::remove_all(root);
    mpi.barrier();

    {
      auto decomposition =
          hundun::runtime::StructuredDecomposition::create(
              mpi, {8, 6, 4}, {true, true, true});
      const auto exchange_plan =
          hundun::runtime::ExchangePlan::create(
              decomposition, decomposition.local_extent(), 2);
      hundun::mesh::MeshTopology topology(decomposition);
      hundun::mesh::MeshGeometry geometry(
          topology,
          hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
      HUNDUN_CHECK(geometry.local_face_count() == topology.local_face_count());
      if (mpi.size() > 1) {
        std::uint64_t local_ghost_faces = topology.ghost_face_count();
        std::uint64_t global_ghost_faces = 0U;
        HUNDUN_CHECK(MPI_Allreduce(&local_ghost_faces, &global_ghost_faces, 1,
                                   MPI_UINT64_T, MPI_SUM,
                                   mpi.comm()) == MPI_SUCCESS);
        HUNDUN_CHECK(global_ghost_faces > 0U);
        if (local_ghost_faces > 0U)
          HUNDUN_CHECK(geometry.local_face_count() >
                       topology.owned_face_count());
      }
      check_fingerprint_ids(decomposition, {"decomposition", "periodicity"});
      check_fingerprint_ids(exchange_plan, {"ghost_width", "regions"});
      check_fingerprint_ids(topology, {"cells.global_id",
                                       "cells.ownership",
                                       "cells.x",
                                       "cells.y",
                                       "cells.z",
                                       "faces.axis",
                                       "faces.global_id",
                                       "faces.neighbour",
                                       "faces.neighbour.present",
                                       "faces.owner",
                                       "faces.ownership",
                                       "faces.pair",
                                       "faces.pair.present",
                                       "faces.patch",
                                       "faces.patch.present",
                                       "faces.x",
                                       "faces.y",
                                       "faces.z",
                                       "mesh_metadata",
                                       "patches.face",
                                       "patches.id",
                                       "patches.name",
                                       "patches.name.length",
                                       "patches.paired",
                                       "patches.paired.present",
                                       "patches.pairing"});
      check_fingerprint_ids(geometry, {"cells", "faces", "faces.local_count",
                                       "mapping.extent", "mapping.kind",
                                       "mapping.metric"});
      hundun::execution::CpuReferenceContext execution;
      hundun::linear::GhostedVector vector(
          execution,
          hundun::linear::VectorLayout::from_topology(topology));
      auto vector_halo = hundun::linear::GhostedVectorHalo::create(
          decomposition, topology, execution);
      check_fingerprint_ids(vector, {"allocation_identity", "backend_identity",
                                     "epoch", "ghost_count", "layout.global_id",
                                     "layout.global_id.count", "local_count",
                                     "owned_count", "space"});
      check_fingerprint_ids(vector_halo,
                            {"ghost_count", "owned_count", "path",
                             "receive_value_count", "send_value_count"});
      auto boundaries =
          hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);
      check_fingerprint_ids(
          boundaries,
          {"inlet_patch.present", "open_domain", "outlet_patch.present",
           "patch.density_rule", "patch.enthalpy_rule", "patch.id",
           "patch.inlet.present", "patch.kind", "patch.mass_flux_rule",
           "patch.name", "patch.name.length", "patch.paired",
           "patch.paired.present", "patch.pressure.present",
           "patch.pressure_rule", "patch.scalar_rule", "patch.velocity_rule",
           "scalar", "scalar.length"});
      hundun::execution::Buffer gamma(execution, topology.local_face_count() *
                                                     sizeof(double));
      auto gamma_values = gamma.view(0U, topology.local_face_count());
      for (std::size_t face = 0; face < gamma_values.size(); ++face)
        gamma_values[face] = 1.0;
      const auto &const_gamma = gamma;
      auto poisson = hundun::finite_volume::MatrixFreePoissonOperator::create(
          decomposition, topology, geometry, execution,
          const_gamma.view(0U, topology.local_face_count()),
          {hundun::finite_volume::PressureConstraintMode::constant_nullspace,
           std::nullopt});
      check_fingerprint_ids(
          poisson, {"constraint_mode", "context.backend_identity",
                    "context.backend_name", "context.backend_name.length",
                    "context.space", "diagonal_available", "domain.global_id",
                    "domain.owned", "pressure_reference_patch",
                    "pressure_reference_patch.present", "range.global_id",
                    "range.owned", "revision", "solver_family"});
      const auto adapter_counters_before = mpi.fp64_reduction_counters();
      for (const auto* source_name :
           {"hundun.runtime.structured_decomposition",
            "hundun.runtime.halo",
            "hundun.mesh.topology", "hundun.mesh.geometry",
            "hundun.linear.ghosted_vector",
            "hundun.linear.ghosted_vector_halo",
            "hundun.finite_volume.poisson"}) {
        CaptureSink sink;
        if (std::string_view(source_name) ==
            "hundun.runtime.structured_decomposition")
          hundun::diagnostics::collect_diagnostics(
              decomposition, local_request, sink);
        else if (std::string_view(source_name) == "hundun.runtime.halo")
          hundun::diagnostics::collect_diagnostics(
              exchange_plan, local_request, sink);
        else if (std::string_view(source_name) == "hundun.mesh.topology")
          hundun::diagnostics::collect_diagnostics(topology, local_request,
                                                   sink);
        else if (std::string_view(source_name) == "hundun.mesh.geometry")
          hundun::diagnostics::collect_diagnostics(geometry, local_request,
                                                   sink);
        else if (std::string_view(source_name) ==
                 "hundun.linear.ghosted_vector")
          hundun::diagnostics::collect_diagnostics(vector, local_request,
                                                   sink);
        else if (std::string_view(source_name) ==
                 "hundun.linear.ghosted_vector_halo")
          hundun::diagnostics::collect_diagnostics(vector_halo, local_request,
                                                   sink);
        else
          hundun::diagnostics::collect_diagnostics(poisson, local_request,
                                                   sink);
        HUNDUN_CHECK(sink.records.size() == 1U);
        HUNDUN_CHECK(sink.records.front().module_id == source_name);
      }
      const auto adapter_counters_after = mpi.fp64_reduction_counters();
      HUNDUN_CHECK(adapter_counters_after.collective_calls ==
                   adapter_counters_before.collective_calls);
      HUNDUN_CHECK(adapter_counters_after.reduced_scalars ==
                   adapter_counters_before.reduced_scalars);
      HUNDUN_CHECK(adapter_counters_after.logical_payload_bytes ==
                   adapter_counters_before.logical_payload_bytes);

      const auto uniform_record =
          hundun::diagnostics::make_mesh_diagnostic_v2(
              mpi.rank(), mpi.size(), topology, geometry);
      const auto uniform_bytes =
          hundun::diagnostics::encode_mesh_diagnostic_v2(uniform_record);
      hundun::diagnostics::write_mesh_diagnostic_v2(
          mpi, root / "mesh-uniform", uniform_record);
      const auto uniform_read =
          hundun::diagnostics::read_mesh_diagnostic_v2_file(
              mesh_file(root / "mesh-uniform", mpi.rank()));
      HUNDUN_CHECK(
          hundun::diagnostics::encode_mesh_diagnostic_v2(uniform_read) ==
          uniform_bytes);
      const auto uniform_incident_evidence =
          check_incident_face_areas(topology, geometry, uniform_read);
      std::uint64_t global_remote_owner_faces = 0U;
      HUNDUN_CHECK(
          MPI_Allreduce(
              &uniform_incident_evidence.remote_owner_face_count,
              &global_remote_owner_faces, 1, MPI_UINT64_T, MPI_SUM,
              mpi.comm()) == MPI_SUCCESS);
      if (mpi.size() > 1)
        HUNDUN_CHECK(global_remote_owner_faces > 0U);

      std::uint64_t local_mutation_count = 0U;
      if (uniform_incident_evidence.mutation_cell) {
        auto incident_area_mutation = uniform_read;
        const double old_limit =
            256.0 * std::numeric_limits<double>::epsilon() *
            uniform_incident_evidence.old_owned_face_area_sum;
        const double complete_limit =
            256.0 * std::numeric_limits<double>::epsilon() *
            uniform_incident_evidence.complete_incident_area_sum;
        HUNDUN_CHECK(complete_limit > old_limit);
        const double mutation =
            old_limit + 0.5 * (complete_limit - old_limit);
        bool found_mutation_cell = false;
        double maximum_closure = 0.0;
        for (auto& cell : incident_area_mutation.cells) {
          if (cell.global_id ==
              *uniform_incident_evidence.mutation_cell) {
            cell.closure_m2 = {mutation, 0.0, 0.0};
            found_mutation_cell = true;
          }
          maximum_closure =
              std::max(maximum_closure,
                       std::sqrt(cell.closure_m2.x * cell.closure_m2.x +
                                 cell.closure_m2.y * cell.closure_m2.y +
                                 cell.closure_m2.z * cell.closure_m2.z));
        }
        HUNDUN_CHECK(found_mutation_cell);
        incident_area_mutation.maximum_cell_closure_norm =
            maximum_closure;
        static_cast<void>(
            hundun::diagnostics::encode_mesh_diagnostic_v2(
                incident_area_mutation));
        local_mutation_count = 1U;
      }
      std::uint64_t global_mutation_count = 0U;
      HUNDUN_CHECK(
          MPI_Allreduce(&local_mutation_count, &global_mutation_count, 1,
                        MPI_UINT64_T, MPI_SUM,
                        mpi.comm()) == MPI_SUCCESS);
      if (mpi.size() > 1)
        HUNDUN_CHECK(global_mutation_count > 0U);

      auto wrong_in_range_owner = uniform_record;
      bool mutated_owner = false;
      for (auto& face : wrong_in_range_owner.faces) {
        if (face.neighbour_global_cell.has_value() &&
            face.owner_global_cell != *face.neighbour_global_cell) {
          face.owner_global_cell = *face.neighbour_global_cell;
          mutated_owner = true;
          break;
        }
      }
      HUNDUN_CHECK(mutated_owner);
      bool rejected_owner = false;
      try {
        static_cast<void>(
            hundun::diagnostics::encode_mesh_diagnostic_v2(
                wrong_in_range_owner));
      } catch (const hundun::runtime::Error&) {
        rejected_owner = true;
      }
      HUNDUN_CHECK(rejected_owner);

      hundun::mesh::MeshGeometry warped_geometry(
          topology, hundun::mesh::AnalyticWarpedBoxMapping(
                        {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0},
                        {0.02, -0.015, 0.01}));
      auto warped_request = local_request;
      warped_request.level =
          hundun::diagnostics::DiagnosticLevel::invariants;
      CaptureSink warped_sink;
      hundun::diagnostics::collect_diagnostics(
          warped_geometry, warped_request, warped_sink);
      HUNDUN_CHECK(warped_sink.records.size() == 1U);
      HUNDUN_CHECK(warped_sink.records.front().module_id ==
                   "hundun.mesh.geometry");
      HUNDUN_CHECK(!warped_sink.records.front().invariants.empty());
      const auto warped_record =
          hundun::diagnostics::make_mesh_diagnostic_v2(
              mpi.rank(), mpi.size(), topology, warped_geometry);
      const auto warped_bytes =
          hundun::diagnostics::encode_mesh_diagnostic_v2(warped_record);
      hundun::diagnostics::write_mesh_diagnostic_v2(
          mpi, root / "mesh-warped", warped_record);
      const auto warped_read =
          hundun::diagnostics::read_mesh_diagnostic_v2_file(
              mesh_file(root / "mesh-warped", mpi.rank()));
      HUNDUN_CHECK(
          hundun::diagnostics::encode_mesh_diagnostic_v2(warped_read) ==
          warped_bytes);
      static_cast<void>(
          check_incident_face_areas(topology, warped_geometry, warped_read));
    }

    const auto directory = root / "records";
    hundun::diagnostics::DiagnosticSession session(directory, 2, mpi.rank());
    const auto batch = batch_for_rank(mpi.rank());
    session.publish(mpi, 2U, batch);
    const auto final_path = rank_file(directory, mpi.rank());
    HUNDUN_CHECK(std::filesystem::is_regular_file(final_path));
    const auto first_bytes = read_bytes(final_path);
    session.publish(mpi, 2U, batch);
    HUNDUN_CHECK(std::filesystem::is_regular_file(final_path));
    HUNDUN_CHECK(read_bytes(final_path) == first_bytes);
    HUNDUN_CHECK(
        !std::filesystem::exists(final_path.string() + ".tmp"));

    const auto blocked = root / "blocked";
    if (mpi.rank() == 0) {
      std::ofstream stream(blocked);
      stream << "not a directory";
    }
    mpi.barrier();
    bool rejected = false;
    try {
      hundun::diagnostics::DiagnosticSession failing(
          blocked / "child", 2, mpi.rank());
      failing.publish(mpi, 2U, batch);
    } catch (const hundun::runtime::Error&) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(
        !std::filesystem::exists(rank_file(blocked / "child", mpi.rank())));

    if (mpi.size() > 1) {
      const auto equivalent_path =
          mpi.rank() == 0 ? root / "normalized" / "." : root / "normalized";
      hundun::diagnostics::DiagnosticSession normalized(
          equivalent_path, 2, mpi.rank());
      normalized.publish(mpi, 2U, batch);
      HUNDUN_CHECK(std::filesystem::is_regular_file(
          rank_file(root / "normalized", mpi.rank())));

      rejected = false;
      try {
        hundun::diagnostics::DiagnosticSession disagreeing(
            root / "rank",
            2, mpi.rank() == 0 ? mpi.rank() + 1 : mpi.rank());
        disagreeing.publish(mpi, 2U, batch);
      } catch (const hundun::runtime::Error&) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(!std::filesystem::exists(root / "rank"));

      rejected = false;
      try {
        hundun::diagnostics::DiagnosticSession disagreeing(
            root / (mpi.rank() == 0 ? "left" : "right"), 2, mpi.rank());
        disagreeing.publish(mpi, 2U, batch);
      } catch (const hundun::runtime::Error&) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);

      rejected = false;
      try {
        hundun::diagnostics::DiagnosticSession disagreeing(
            root / "schedule", mpi.rank() == 0 ? 2 : 3, mpi.rank());
        disagreeing.publish(mpi, 2U, batch);
      } catch (const hundun::runtime::Error&) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(!std::filesystem::exists(root / "schedule"));

      rejected = false;
      try {
        hundun::diagnostics::DiagnosticSession disagreeing(
            root / "step", 2, mpi.rank());
        disagreeing.publish(mpi, mpi.rank() == 0 ? 2U : 4U, batch);
      } catch (const hundun::runtime::Error&) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(!std::filesystem::exists(root / "step"));
    }

    const auto rename_directory = root / "rename";
    const auto rename_final = rank_file(rename_directory, mpi.rank());
    std::filesystem::create_directories(rename_final);
    rejected = false;
    try {
      hundun::diagnostics::DiagnosticSession rename_failing(
          rename_directory, 2, mpi.rank());
      rename_failing.publish(mpi, 2U, batch);
    } catch (const hundun::runtime::Error&) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(std::filesystem::is_directory(rename_final));
    HUNDUN_CHECK(
        !std::filesystem::exists(rename_final.string() + ".tmp"));

    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(root);
    mpi.barrier();
  });
}
