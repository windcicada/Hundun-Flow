// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include "../../src/solver_conservative_energy_detail.hpp"
#include "hundun/v04_boundary.hpp"
#include "hundun/v04_execution.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <mpi.h>
using namespace hundun::v04;
CartesianMeshSpec mesh_spec(std::int32_t n, bool stretched) {
  CartesianMeshSpec mesh;
  mesh.kind =
      stretched ? GeometryKind::tensor_stretched : GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  // Refine one smooth-family monitor: the permitted adjacent ratio tends to
  // one with h, so the 16/32/64 sequence samples a smooth tensor mapping
  // instead of preserving a fixed grid-spacing jump.
  mesh.max_growth_ratio = stretched ? 1.0 + 0.8 / static_cast<double>(n) : 1.0;
  if (stretched) {
    const double inverse_n = 1.0 / static_cast<double>(n);
    mesh.has_base_spacing = true;
    mesh.base_spacing = {1.08 * inverse_n, 1.08 * inverse_n, 1.08 * inverse_n};
    mesh.minimum_spacing = {0.92 * inverse_n, 0.92 * inverse_n,
                            0.92 * inverse_n};
    mesh.focus_regions.push_back(
        {{0.35, 0.35, 0.35},
         {0.65, 0.65, 0.65},
         {0.96 * inverse_n, 0.96 * inverse_n, 0.96 * inverse_n}});
  } else {
    mesh.minimum_spacing = {0.25 / n, 0.25 / n, 0.25 / n};
  }
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(n) *
                                 static_cast<std::uint64_t>(n) *
                                 static_cast<std::uint64_t>(n);
  mesh.limits.max_memory_bytes_per_rank = 1U << 30U;
  return mesh;
}

ValidatedModel periodic_model(
    const CartesianMeshSpec &mesh,
    ConvectionScheme convection_scheme = ConvectionScheme::central2) {
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x94c17e32U;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec &face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
    face.mach_limit = 0.95;
  }
  model.schemes.momentum = convection_scheme;
  model.schemes.enthalpy = convection_scheme;
  model.schemes.species = convection_scheme;
  model.schemes.passive_scalar = convection_scheme;
  model.schemes.diffusion = DiffusionScheme::central2;
  model.schemes.limiter = 1.0;
  return model;
}

struct KernelFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  CartesianKernelPlan kernels;
};

bool make_fixture(
    std::int32_t n, bool stretched, KernelFixture &fixture,
    ConvectionScheme convection_scheme = ConvectionScheme::central2) {
  const CartesianMeshSpec mesh = mesh_spec(n, stretched);
  FieldRegistry registry;
  TimeSchemePlan time;
  return CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh,
                                            GeometryBudget{}, fixture.geometry,
                                            fixture.patch) &&
         BoundaryCompiler::compile(MPI_COMM_SELF,
                                   periodic_model(mesh, convection_scheme),
                                   fixture.geometry, fixture.patch, registry,
                                   fixture.boundary, fixture.schemes, time) &&
         CartesianKernelPlan::compile(fixture.schemes, fixture.geometry,
                                      fixture.patch, fixture.boundary,
                                      fixture.kernels);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int failures = 0, cases = 0;
  for (bool stretched : {false, true})
    for (auto scheme :
         {ConvectionScheme::central2, ConvectionScheme::limited_central2,
          ConvectionScheme::tvd2}) {
      KernelFixture f;
      if (!make_fixture(17, stretched, f, scheme)) {
        ++failures;
        continue;
      }
      for (auto axis : {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z})
        for (double mass : {-1.0, 1.0})
          for (auto q : {std::array<double, 4>{1, 2, 4, 7},
                         std::array<double, 4>{7, 4, 2, 1},
                         std::array<double, 4>{1, 3, 2, 4},
                         std::array<double, 4>{-14, -13.5, -9, 15}}) {
            const std::array<double, 4> d{0.13, -0.29, 0.41, -0.17};
            const Int3 face{8, 8, 8};
            double analytical = detail::sampled_convection_direction(
                f.kernels, scheme, q, d, axis, face, mass);
            auto hi = q, lo = q;
            constexpr double eps = 1e-5;
            for (int j = 0; j < 4; ++j) {
              hi[j] += eps * d[j];
              lo[j] -= eps * d[j];
            }
            double fd = (detail::sampled_convection_face(f.kernels, scheme, hi,
                                                         axis, face, mass) -
                         detail::sampled_convection_face(f.kernels, scheme, lo,
                                                         axis, face, mass)) /
                        (2 * eps);
            double linear = 0;
            for (int j = 0; j < 4; ++j) {
              std::array<double, 4> b{};
              b[j] = 1;
              linear += d[j] * detail::sampled_convection_direction(
                                   f.kernels, scheme, q, b, axis, face, mass);
            }
            ++cases;
            if (!std::isfinite(analytical) ||
                std::abs(analytical - fd) > 1e-7 ||
                std::abs(analytical - linear) > 1e-12) {
              ++failures;
              printf("FAIL scheme=%d stretched=%d mass=%g analytic=%g fd=%g "
                     "linear=%g\n",
                     int(scheme), stretched, mass, analytical, fd, linear);
            }
          }
    }
  printf("sampled derivative cases=%d failures=%d\n", cases, failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
