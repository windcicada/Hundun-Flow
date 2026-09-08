// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/candidate_boundary_fixture.hpp"
#include "../../src/solver_scalar_mass_remap_detail.hpp"

#include <mpi.h>
#include <array>
#include <iostream>
#include <limits>

using namespace hundun::v04;
using namespace hundun::v04::test;

bool run() {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.multispecies = true;
  if (!fixture.initialize(MPI_COMM_WORLD, spec)) return false;
  std::array<FieldId, 2> fields{};
  std::array<FieldView, 2> scalars{};
  const std::array<TransportedScalarRole, 2> roles{
      TransportedScalarRole::species, TransportedScalarRole::species};
  for (std::size_t s=0; s<2; ++s) {
    fields[s] = fixture.independent_species[s].view.field;
    scalars[s] = fixture.independent_species[s].view;
    fill(fixture.independent_species[s], 0.0);
  }
  const auto cells = fixture.patch.cells;
  for (int z=0; z<cells.z; ++z)
    for (int y=0; y<cells.y; ++y)
      for (int x=0; x<cells.x; ++x)
        fixture.density.view.unchecked({x,y,z},0U) =
            1.0 / fixture.cell_volume({x,y,z});
  detail::ScalarMassRemap remap;
  Status status = remap.allocate(fixture.patch, {fields.data(),2U},
                                 {roles.data(),2U}, 2U);
  if (status) status = remap.bind(MPI_COMM_WORLD, fixture.boundary);
  FaceFluxStorage storage;
  FaceFluxView base, final;
  if (status) status = FaceFluxStorage::allocate_workspace(cells,2U,storage);
  if (status) status = storage.workspace_view(0U,81U,base);
  if (status) status = storage.workspace_view(1U,82U,final);
  if (!status) return false;
  for (auto flux : {base,final})
    for (auto face : {flux.x,flux.y,flux.z})
      for (int z=0; z<face.extents.z; ++z)
        for (int y=0; y<face.extents.y; ++y)
          for (int x=0; x<face.extents.x; ++x) face.unchecked({x,y,z})=0.0;
  status = remap.capture(fixture.kernels, as_const(fixture.density.view),
                         {scalars.data(),2U}, as_const(base), 1.0);
  // Fixed physical donor is 0.2/0.3. The extended-precision incoming
  // quantity is nonzero, but its correctly rounded binary64 q is zero.
  // A relative residual of one must not demand an unrepresentable update.
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y)
        final.x.unchecked({0,y,z}) = std::numeric_limits<double>::denorm_min();
  const BoundaryResolvedValues values{
      {fixture.boundary_scalar_values.data(),fixture.boundary_scalar_values.size()},
      {fixture.boundary_vector_values.data(),fixture.boundary_vector_values.size()},
      {fixture.boundary_normal_gradient_values.data(),fixture.boundary_normal_gradient_values.size()}};
  detail::ScalarMassRemap::Report report;
  if (status) status = remap.solve(fixture.kernels, as_const(fixture.density.view),
      as_const(fixture.velocity.view), as_const(final), {scalars.data(),2U}, {},
      values, fixture.reductions, report);
  if (!status) std::cerr << "FAIL: representable trace remap "
      << unsigned(status.code) << '/' << status.detail << " residual="
      << report.residual << " iterations=" << report.iterations << '\n';
  bool passed = status && report.residual == 1.0 &&
      report.convergence_residual == 0.0 && report.iterations == 0U;
  remap.copy_solution({scalars.data(),2U}, TransportedScalarRole::species);
  for (auto scalar : scalars)
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y)
        for (int x=0; x<cells.x; ++x)
          passed &= scalar.unchecked({x,y,z},0U) == 0.0;
  // Ordinary representable transport still requires a real update; the
  // quantization gate must not mistake an initially zero field for closure.
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y) final.x.unchecked({0,y,z}) = 1.0e-3;
  status = remap.solve(fixture.kernels, as_const(fixture.density.view),
      as_const(fixture.velocity.view), as_const(final), {scalars.data(),2U}, {},
      values, fixture.reductions, report);
  passed &= status && report.iterations > 0U &&
      report.residual <= detail::ScalarMassRemap::tolerance;
  remap.copy_solution({scalars.data(),2U}, TransportedScalarRole::species);
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y) {
        passed &= std::abs(scalars[0].unchecked({0,y,z},0U)-2.0e-4) < 1e-18;
        passed &= std::abs(scalars[1].unchecked({0,y,z},0U)-3.0e-4) < 1e-18;
      }
  if (!passed) std::cerr << "FAIL: quantization and ordinary update contracts\n";
  return passed;
}

int main(int argc, char** argv) {
  if (MPI_Init(&argc,&argv) != MPI_SUCCESS) return 2;
  int local=run() ? 1 : 0, global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
