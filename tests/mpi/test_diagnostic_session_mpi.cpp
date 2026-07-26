// SPDX-License-Identifier: Apache-2.0

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

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
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

}  // namespace

int main(int argc, char** argv) {
  return hundun::test::run([&] {
    hundun::runtime::MpiEnvironment environment(argc, argv);
    auto mpi =
        hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
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
          [&](hundun::diagnostics::DiagnosticRequest mismatched) {
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
            HUNDUN_CHECK(lowest == 0);
          };
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.rank = 1;
        expect_collective_rejection(std::move(mismatched));
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.level =
              hundun::diagnostics::DiagnosticLevel::counters;
        expect_collective_rejection(std::move(mismatched));
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.step += 1U;
        expect_collective_rejection(std::move(mismatched));
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.time_s = 0.5;
        expect_collective_rejection(std::move(mismatched));
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.frame.phase = "different-phase";
        expect_collective_rejection(std::move(mismatched));
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.selected_fields = {"rank"};
        expect_collective_rejection(std::move(mismatched));
      }
      {
        auto mismatched = collective_request;
        if (mpi.rank() == 0)
          mismatched.sample_budget = 1U;
        expect_collective_rejection(std::move(mismatched));
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
          topology, hundun::mesh::UniformBoxMapping(
                        {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
      hundun::execution::CpuReferenceContext execution;
      hundun::linear::GhostedVector vector(
          execution,
          hundun::linear::VectorLayout::from_topology(topology));
      auto vector_halo = hundun::linear::GhostedVectorHalo::create(
          decomposition, topology, execution);
      hundun::execution::Buffer gamma(
          execution, topology.local_face_count() * sizeof(double));
      auto gamma_values = gamma.view(0U, topology.local_face_count());
      for (std::size_t face = 0; face < gamma_values.size(); ++face)
        gamma_values[face] = 1.0;
      const auto& const_gamma = gamma;
      auto poisson =
          hundun::finite_volume::MatrixFreePoissonOperator::create(
              decomposition, topology, geometry, execution,
              const_gamma.view(0U, topology.local_face_count()),
              {hundun::finite_volume::PressureConstraintMode::
                   constant_nullspace,
               std::nullopt});
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
