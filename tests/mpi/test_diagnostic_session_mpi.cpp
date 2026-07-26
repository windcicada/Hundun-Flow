// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/diagnostic_session.hpp"
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
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

class CaptureSink final : public hundun::diagnostics::DiagnosticSink {
 public:
  void submit(const hundun::diagnostics::DiagnosticRecord& record) override {
    records.push_back(record);
  }
  std::vector<hundun::diagnostics::DiagnosticRecord> records;
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
    const auto counters_after_collective = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(counters_after_collective.collective_calls ==
                 counters_before.collective_calls);
    HUNDUN_CHECK(counters_after_collective.reduced_scalars ==
                 counters_before.reduced_scalars);
    HUNDUN_CHECK(counters_after_collective.logical_payload_bytes ==
                 counters_before.logical_payload_bytes);

    {
      auto decomposition =
          hundun::runtime::StructuredDecomposition::create(
              mpi, {8, 6, 4}, {true, true, true});
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
            "hundun.mesh.topology", "hundun.mesh.geometry",
            "hundun.linear.ghosted_vector",
            "hundun.linear.ghosted_vector_halo",
            "hundun.finite_volume.poisson"}) {
        CaptureSink sink;
        if (std::string_view(source_name) ==
            "hundun.runtime.structured_decomposition")
          hundun::diagnostics::collect_diagnostics(
              decomposition, local_request, sink);
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
    }
    const auto root = std::filesystem::temp_directory_path() /
                      ("hundun-task24-session-mpi-" +
                       std::to_string(mpi.size()));
    if (mpi.rank() == 0)
      std::filesystem::remove_all(root);
    mpi.barrier();

    const auto directory = root / "records";
    hundun::diagnostics::DiagnosticSession session(directory, 2, mpi.rank());
    const auto batch = batch_for_rank(mpi.rank());
    session.publish(mpi, 2U, batch);
    HUNDUN_CHECK(std::filesystem::is_regular_file(
        rank_file(directory, mpi.rank())));
    session.publish(mpi, 2U, batch);
    HUNDUN_CHECK(std::filesystem::is_regular_file(
        rank_file(directory, mpi.rank())));

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

    if (mpi.size() > 1) {
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

    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(root);
    mpi.barrier();
  });
}
