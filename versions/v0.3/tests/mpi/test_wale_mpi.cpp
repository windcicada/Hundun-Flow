// SPDX-License-Identifier: Apache-2.0

#include "hundun/les_wale.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {

using namespace hundun;

runtime::Int3 process_grid(int ranks) {
  return ranks == 1 ? runtime::Int3{1, 1, 1} : runtime::Int3{2, 1, 1};
}

runtime::FieldDescriptor field(const char *name, std::uint32_t components) {
  return {name, "1", "task12-wale", runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64, components, 1, false,
          runtime::RestartPolicy::transient, runtime::OutputPolicy::never};
}

struct Index final { int i{}, j{}, k{}; };

class DeviceContext final : public execution::ExecutionContext {
public:
  std::string_view backend_name() const noexcept override { return "device"; }
  execution::BackendIdentity backend_identity() const noexcept override {
    return 99U;
  }
  execution::ExecutionSpace space() const noexcept override {
    return execution::ExecutionSpace::device;
  }
  bool ordered() const noexcept override { return true; }
  bool supports(execution::ExecutionCapability) const noexcept override {
    return true;
  }
};

Index map(runtime::Int3 global, runtime::Box3 box, runtime::Int3 extent) {
  const runtime::Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                            box.end.z - box.begin.z};
  const auto axis = [](int x, int begin, int end, int n, int local_n) {
    if (x >= begin && x < end) return x - begin;
    if (x == begin - 1 || (begin == 0 && x == n - 1)) return -1;
    if (x == end || (end == n && x == 0)) return local_n;
    throw runtime::Error("test WALE mapping failed");
  };
  return {axis(global.x, box.begin.x, box.end.x, extent.x, local.x),
          axis(global.y, box.begin.y, box.end.y, extent.y, local.y),
          axis(global.z, box.begin.z, box.end.z, extent.z, local.z)};
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {8, 6, 4}, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {2.0, 1.5, 1.0}));
  std::vector<mesh::GlobalCellId> active;
  active.reserve(topology.local_cell_count());
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    active.push_back(topology.global_cell_id(cell));
  for (mesh::LocalCellId cell = topology.owned_cell_count();
       cell < topology.local_cell_count(); ++cell)
    active.push_back(topology.global_cell_id(cell));

  runtime::FieldRegistry registry;
  const auto gradient_id = registry.declare_field(field("gradient", 9U));
  const auto density_id = registry.declare_field(field("density", 1U));
  registry.freeze();
  const auto box = topology.owned_global_box();
  const runtime::Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                            box.end.z - box.begin.z};
  runtime::FieldStorage storage(registry, local);
  auto gradient = storage.view<double>(gradient_id);
  auto density = storage.view<double>(density_id);
  const std::array<double, 9> tensor{0.0, 1.0, 0.0, 0.25, 0.0,
                                     0.0, 0.0, 0.0, 0.0};
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell) {
    const auto index = map(topology.global_cell(cell), box,
                           topology.global_extent());
    for (std::size_t c = 0U; c < tensor.size(); ++c)
      gradient(index.i, index.j, index.k, static_cast<int>(c)) = tensor[c];
    density(index.i, index.j, index.k, 0) = 2.0;
  }

  execution::CpuReferenceContext execution;
  bool invalid_control_rejected = false;
  try {
    static_cast<void>(les::WaleModel::create(
        {0.0, 0.9, 0.7}, topology, geometry, topology.owned_cell_count(),
        active, execution));
  } catch (const runtime::Error &) {
    invalid_control_rejected = true;
  }
  HUNDUN_CHECK(invalid_control_rejected);
  DeviceContext device;
  bool device_rejected = false;
  try {
    static_cast<void>(les::WaleModel::create(
        {0.5, 0.9, 0.7}, topology, geometry, topology.owned_cell_count(),
        active, device));
  } catch (const runtime::Error &) {
    device_rejected = true;
  }
  HUNDUN_CHECK(device_rejected);
  auto model = les::WaleModel::create(
      {0.5, 0.9, 0.7}, topology, geometry, topology.owned_cell_count(), active,
      execution);
  const auto make_input = [&](std::uint64_t density_fingerprint) {
    return les::WaleAttemptInput{
        3U, 0.01, les::WaleTimeOrder::bdf2, 11U, 12U, 13U,
        density_fingerprint,
        static_cast<const runtime::FieldStorage &>(storage)
            .view<double>(gradient_id),
        static_cast<const runtime::FieldStorage &>(storage)
            .view<double>(density_id)};
  };
  auto first = model.evaluate(make_input(14U));
  HUNDUN_CHECK(first.owned_active_count() == topology.owned_cell_count());
  HUNDUN_CHECK(first.local_active_count() == topology.local_cell_count());
  HUNDUN_CHECK(first.identity().value != 0U);
  const auto nu = first.nu_t_m2_per_s();
  const auto mu = first.mu_sgs_pa_s();
  for (std::size_t row = 0U; row < active.size(); ++row) {
    HUNDUN_CHECK(std::isfinite(nu[row]) && nu[row] > 0.0);
    HUNDUN_CHECK(mu[row] == 2.0 * nu[row]);
  }
  auto identity_mutation = model.evaluate(make_input(15U));
  HUNDUN_CHECK(identity_mutation.identity() != first.identity());

  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell) {
    const auto index = map(topology.global_cell(cell), box,
                           topology.global_extent());
    density(index.i, index.j, index.k, 0) = 3.0;
  }
  auto fresh_density = model.evaluate(make_input(16U));
  const auto fresh_mu = fresh_density.mu_sgs_pa_s();
  for (std::size_t row = 0U; row < active.size(); ++row)
    HUNDUN_CHECK(fresh_mu[row] == 3.0 * nu[row]);

  const auto stale = make_input(17U);
  storage.begin_rebuild();
  bool stale_rejected = false;
  try { static_cast<void>(model.evaluate(stale)); }
  catch (const runtime::Error &) { stale_rejected = true; }
  HUNDUN_CHECK(stale_rejected);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run(run);
}
