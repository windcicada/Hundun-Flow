// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#include "../../src/core_patch_inlets_detail.hpp"
#include "../support/ibm_force_fixture.hpp"
#include "../support/product_fixture.hpp"

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace hundun::v04;
using namespace hundun::v04::test;
namespace {
bool check(bool value, const char* message) {
  if (!value) std::cerr << "FAIL: " << message << '\n';
  return value;
}
bool check(Status value, const char* message) {
  if (!value) std::cerr << "status=" << static_cast<unsigned>(value.code) << '/' << value.detail << '\n';
  return check(static_cast<bool>(value), message);
}
PlanFingerprint write_labels(const std::filesystem::path& path,
                            const std::vector<std::int32_t>& labels) {
  std::ofstream file(path, std::ios::binary | std::ios::trunc);
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto label : labels)
    for (unsigned byte = 0U; byte < 4U; ++byte) {
      const auto value = static_cast<unsigned char>(static_cast<std::uint32_t>(label) >> (8U * byte));
      file.put(static_cast<char>(value));
      hash ^= value;
      hash *= UINT64_C(1099511628211);
    }
  return hash == 0U ? 1U : hash;
}
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool ok = true;
  char directory[] = "/tmp/hundun-patch-test.XXXXXX";
  const char* created = ::mkdtemp(directory);
  ok &= check(created != nullptr, "scratch directory");
  if (!ok) { MPI_Finalize(); return 1; }
  const std::filesystem::path root{created};
  constexpr int n = 16;
  IbmForceFixture fixture;
  ok &= check(fixture.initialize(MPI_COMM_SELF, n), "geometry fixture");
  ValidatedModel model = product_model({n,n,n});
  model.mesh = force_mesh(n);
  BoundaryFaceSpec air;
  air.flow_kind = BoundaryKind::mass_flow_inlet;
  air.mass_flow_rate = 0.01;
  air.pressure = 100000;
  air.temperature = 295;
  air.direction = {1,0,0};
  air.relaxation = 1;
  model.boundaries[0] = air;
  BoundaryFaceSpec fuel = air;
  fuel.mass_flow_rate = 0.002;
  fuel.direction = {0,1,0};
  model.patch_inlets.emplace();
  model.patch_inlets->labels_file = "labels.d";
  model.patch_inlets->patches = {{7,0,false,air},{11,2,true,fuel}};
  std::vector<std::int32_t> labels(static_cast<std::size_t>(n)*n*n, 0);
  labels[0] = 7;
  labels[static_cast<std::size_t>(n)] = 7;
  const auto links = fixture.topology.links();
  const ImmersedLink* source = nullptr;
  for (std::size_t i = 0U; i < links.size; ++i)
    if (links.data[i].direction == ImmersedFaceDirection::y_negative) {
      source = &links.data[i]; break;
    }
  ok &= check(source != nullptr, "fuel face found");
  ThermodynamicsPlan thermo;
  ok &= check(ThermodynamicsPlan::compile(model.thermophysics, {}, thermo), "thermodynamics");
  if (ok) {
    labels[detail::patch_flat({n,n,n}, source->fluid_global_index)] = 11;
    labels[detail::patch_flat({n,n,n}, source->solid_global_index)] = -11;
    model.patch_inlets->labels_fingerprint = write_labels(root / "labels.d", labels);
    detail::CompiledPatchInlets plan;
    ok &= check(detail::compile_patch_inlets(model, root, fixture.geometry, fixture.patch,
        &fixture.topology, thermo, plan), "real fluid owners and solid backing bind");
    if (!plan.immersed_states.empty()) {
      const auto& state = plan.immersed_states.front();
      ok &= check(plan.immersed_states.size() == 1U && state.global_link == source->global_link &&
          std::abs(state.face_mass_flux - 0.002) < 1e-18 && state.velocity.y > 0.0 &&
          state.enthalpy > 0.0, "fuel state and oriented exact mass source");
      ok &= check(plan.face_patch[0][0] == 0 && plan.face_patch[0][1] == 0 &&
          plan.face_patch[0][2] == -1, "air support does not extend to unlabelled cells");
    } else ok &= check(false, "fuel source not dropped");
    const auto identity = model.patch_inlets->labels_fingerprint;
    labels[0] = 0;
    write_labels(root / "labels.d", labels);
    detail::CompiledPatchInlets rejected;
    ok &= check(detail::compile_patch_inlets(model, root, fixture.geometry, fixture.patch,
        &fixture.topology, thermo, rejected).detail == 15903U,
        "file changed after parse cannot be frozen");
    labels[0] = 7;
    labels[detail::patch_flat({n,n,n}, source->solid_global_index)] = 0;
    model.patch_inlets->labels_fingerprint = write_labels(root / "labels.d", labels);
    ok &= check(!detail::compile_patch_inlets(model, root, fixture.geometry, fixture.patch,
        &fixture.topology, thermo, rejected), "fuel without negative backing label rejects");
    labels[detail::patch_flat({n,n,n}, source->solid_global_index)] = -11;
    model.patch_inlets->labels_fingerprint = write_labels(root / "labels.d", labels);
    ok &= check(model.patch_inlets->labels_fingerprint == identity, "restored binary identity");
    model.patch_inlets->patches[0].boundary.temperature = 300;
    ok &= check(!detail::compile_patch_inlets(model, root, fixture.geometry, fixture.patch,
        &fixture.topology, thermo, rejected), "unsupported heterogeneous outer thermal target rejects");
  }
  std::filesystem::remove(root / "labels.d");
  std::filesystem::remove(root);
  MPI_Finalize();
  return ok ? 0 : 1;
}
