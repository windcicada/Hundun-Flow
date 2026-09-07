// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/ibm_force_fixture.hpp"
#include "../support/product_fixture.hpp"

#include <iomanip>
#include <iostream>
#include <numeric>
#include <string_view>

namespace {
using namespace hundun::v04;
using namespace hundun::v04::test;
constexpr double pi = 3.14159265358979323846;
constexpr double tau = 0.02;
constexpr double radius = 0.65;

bool check(Status s, const char* phase) {
  if (!s) std::cerr << phase << " status=" << static_cast<int>(s.code)
                    << '/' << s.detail << '\n';
  return static_cast<bool>(s);
}

std::vector<TriangleInput> surface_for(int shape, int n, int cylinder_facets = 0) {
  if (shape != 2) {
    const auto cube = force_cube();
    std::vector<TriangleInput> result(cube.begin(), cube.end());
    if (shape == 1) {
      const auto rotate = [](Real3 p) {
        return Real3{0.8 * p.x - 0.6 * p.y,
                     0.6 * p.x + 0.8 * p.y, p.z};
      };
      for (auto& t : result) t = {rotate(t.a), rotate(t.b), rotate(t.c)};
    }
    return result;
  }
  std::vector<TriangleInput> result;
  const int facets = cylinder_facets == 0 ? 4 * n : cylinder_facets;
  for (int i = 0; i < facets; ++i) {
    const double a = 2.0 * pi * i / facets;
    // Weld the seam by reusing the first angular coordinate exactly; sin(2pi)
    // is not bitwise zero and would create an open surface in the fixture.
    const double b = 2.0 * pi * ((i + 1) % facets) / facets;
    const Real3 p{radius * std::cos(a), radius * std::sin(a), -0.45};
    const Real3 q{radius * std::cos(b), radius * std::sin(b), -0.45};
    const Real3 r{p.x, p.y, 0.45}, s{q.x, q.y, 0.45};
    result.push_back({p, q, s});
    result.push_back({p, s, r});
    result.push_back({{0.0, 0.0, -0.45}, q, p});
    result.push_back({{0.0, 0.0, 0.45}, r, s});
  }
  return result;
}

// Exact Cartesian cell averages of q and Laplacian(q), not point samples.
// The analytic cube/cylinder walls have grad(q).normal=0. For the cylinder,
// the default polygonal STL error decreases with h. The separate geometry
// experiment freezes its facets instead; neither is a curved cut-volume scheme.
std::array<double, 2U> exact(int shape, Real3 p, double h) {
  if (shape != 2) {
    const double c = shape == 0 ? 1.0 : 0.8;
    const double s = shape == 0 ? 0.0 : 0.6;
    const auto sinc = [](double v) { return v == 0.0 ? 1.0 : std::sin(v) / v; };
    const double q = std::cos(2.0 * pi * (c * p.x + s * p.y)) *
                     sinc(pi * h * c) * sinc(pi * h * s);
    return {q, -4.0 * pi * pi * q};
  }
  const double x2 = p.x * p.x + h * h / 12.0;
  const double y2 = p.y * p.y + h * h / 12.0;
  const auto fourth = [h](double x) {
    return x * x * x * x + 0.5 * x * x * h * h + std::pow(h, 4) / 80.0;
  };
  return {fourth(p.x) + fourth(p.y) + 2.0 * x2 * y2 -
              2.0 * radius * radius * (x2 + y2) + std::pow(radius, 4),
          16.0 * (x2 + y2) - 8.0 * radius * radius};
}

bool run(int shape, int n, double& solution_error, int cylinder_facets = 0) {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  EBTopology topology;
  BoundaryStencilPlan boundary;
  ImmersedPlanLimits limits;
  limits.stencil.policy = IbmReconstructionPolicy::adaptive_order;
  const auto triangles = surface_for(shape, n, cylinder_facets);
  if (!check(CartesianGeometryCompiler::compile(MPI_COMM_SELF, force_mesh(n), {},
      geometry, patch), "mesh") ||
      !check(StlScanCompiler::compile_triangles(geometry, patch,
          {triangles.data(), triangles.size()}, CartesianAxis::y, kForceScanBudget, scan), "scan") ||
      !check(ImmersedSurfaceCompiler::compile(scan, surface), "surface") ||
      !check(EBTopologyCompiler::compile(MPI_COMM_SELF, geometry, patch, scan, surface,
          ImmersedFluidSide::inside, limits, topology), "topology") ||
      !check(BoundaryStencilCompiler::compile(MPI_COMM_SELF, geometry, patch, surface,
          topology, limits, boundary), "boundary")) return false;
  auto model = product_model({n, n, n});
  model.mesh = force_mesh(n);
  FieldRegistry registry;
  BoundaryPlan physical;
  SchemePlan schemes;
  TimeSchemePlan time;
  CartesianKernelPlan kernels;
  IbmEquationInterfacePlan ibm;
  if (!check(BoundaryCompiler::compile(MPI_COMM_SELF, model, geometry, patch,
      registry, physical, schemes, time), "physical") ||
      !check(CartesianKernelPlan::compile(schemes, geometry, patch, physical, kernels), "kernels") ||
      !check(IbmEquationInterfacePlan::compile(kernels, topology, boundary,
          topology.interface_metric(), ibm), "interface")) return false;
  const Int3 cells = patch.cells;
  const double h = 3.0 / n;
  const auto region = topology.region();
  auto q = make_force_field(1U, cells, 1U, 1U, 1U, 1U);
  auto d = make_force_field(2U, cells, 1U, 1U, 1U, 2U);
  auto rate = make_force_field(3U, cells, 1U, 0U, 1U, 3U);
  std::fill(d.storage.begin(), d.storage.end(), 1.0);
  std::vector<Int3> active;
  std::vector<double> reference, laplacian;
  for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
    if (region.data[x + n * (y + n * z)] == 0U) continue;
    active.push_back({x, y, z});
    const auto value = exact(shape, {-1.5 + (x + 0.5) * h,
        -1.5 + (y + 0.5) * h, -1.5 + (z + 0.5) * h}, h);
    reference.push_back(value[0U]);
    laplacian.push_back(value[1U]);
  }
  if (active.empty()) return false;
  const auto apply = [&](const std::vector<double>& values, std::vector<double>& out) {
    for (std::size_t i = 0; i < active.size(); ++i) q.view.unchecked(active[i], 0U) = values[i];
    const std::array<ConstFieldView, 1U> reads{as_const(q.view)};
    const std::array<FieldView, 1U> writes{rate.view};
    if (!check(cartesian_diffusion(kernels, as_const(d.view),
        {{reads.data(), reads.size()}, {writes.data(), writes.size()},
         {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr}), "diffusion") ||
        !check(ibm.correct_impermeable_scalar_diffusion(as_const(q.view),
            as_const(d.view), rate.view), "impermeable")) return false;
    for (std::size_t i = 0; i < active.size(); ++i)
      out[i] = values[i] - tau * rate.view.unchecked(active[i], 0U);
    return true;
  };
  const std::size_t size = active.size();
  std::vector<double> out(size), ones(size, 1.0), rhs(size), solution(size), residual(size), direction(size);
  if (!apply(ones, out)) return false;
  double constant_error = 0.0;
  for (double v : out) constant_error = std::max(constant_error, std::abs(v - 1.0));
  if (!apply(reference, out)) return false;
  double flux_sum = 0.0, flux_l1 = 0.0, row_error = 0.0, truncation = 0.0;
  const std::array<Int3, 6U> offsets{{{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}}};
  for (std::size_t i = 0; i < size; ++i) {
    const Int3 p = active[i];
    const double value = rate.view.unchecked(p, 0U);
    double paired = 0.0;
    for (const auto step : offsets) {
      const Int3 neighbor{p.x + step.x, p.y + step.y, p.z + step.z};
      if (region.data[neighbor.x + n * (neighbor.y + n * neighbor.z)] != 0U)
        paired += (q.view.unchecked(neighbor, 0U) - reference[i]) / (h * h);
    }
    row_error = std::max(row_error, std::abs(value - paired));
    flux_sum += value * h * h * h;
    flux_l1 += std::abs(value) * h * h * h;
    truncation += (value - laplacian[i]) * (value - laplacian[i]);
    rhs[i] = reference[i] - tau * laplacian[i];
  }
  // Test-only CG for the symmetric uniform-grid Helmholtz operator. This is
  // not a change of the production pressure/energy Krylov method.
  residual = rhs;
  direction = rhs;
  const auto dot = [](const auto& a, const auto& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  };
  double rr = dot(residual, residual);
  const double initial_rr = rr;
  int iterations = 0;
  while (rr > 1.0e-24 * initial_rr && iterations < 1000) {
    if (!apply(direction, out)) return false;
    const double denominator = dot(direction, out);
    if (!(denominator > 0.0)) return false;
    const double alpha = rr / denominator;
    for (std::size_t i = 0; i < size; ++i) {
      solution[i] += alpha * direction[i];
      residual[i] -= alpha * out[i];
    }
    const double next = dot(residual, residual);
    for (std::size_t i = 0; i < size; ++i) direction[i] = residual[i] + (next / rr) * direction[i];
    rr = next;
    ++iterations;
  }
  if (!apply(solution, out)) return false;
  double error = 0.0, true_residual = 0.0;
  for (std::size_t i = 0; i < size; ++i) {
    error += std::pow(solution[i] - reference[i], 2);
    true_residual += std::pow(out[i] - rhs[i], 2);
  }
  solution_error = std::sqrt(error / size);
  if (shape == 2) {
    const int facets = cylinder_facets == 0 ? 4 * n : cylinder_facets;
    const double analytic_volume = pi * radius * radius * 0.9;
    const double stl_volume = 0.5 * facets * radius * radius *
                              std::sin(2.0 * pi / facets) * 0.9;
    const double binary_volume = size * h * h * h;
    // Analytic q obeys Neumann data on the circle, not exactly on the polygon.
    // This face-midpoint mismatch is reported separately from the discrete
    // zero-flux row check; it must not be mistaken for solver non-conservation.
    const double midpoint_radius = radius * std::cos(pi / facets);
    const double midpoint_normal_gradient =
        4.0 * (midpoint_radius * midpoint_radius - radius * radius) * midpoint_radius;
    std::cout << std::setprecision(17) << "IBM_SCALAR_GEOMETRY n=" << n
              << " facets=" << facets << " binary_volume=" << binary_volume
              << " stl_volume=" << stl_volume << " analytic_volume=" << analytic_volume
              << " binary_minus_stl=" << binary_volume - stl_volume
              << " stl_minus_analytic=" << stl_volume - analytic_volume
              << " analytic_midface_normal_gradient=" << midpoint_normal_gradient << '\n';
  }
  std::cout << std::setprecision(17) << "IBM_SCALAR_MMS shape=" << shape << " n=" << n
            << " active=" << size << " l2=" << solution_error
            << " truncation_l2=" << std::sqrt(truncation / size)
            << " paired_row_error=" << row_error << " net_flux=" << flux_sum
            << " constant_error=" << constant_error << " iterations=" << iterations
            << " true_relative=" << std::sqrt(true_residual / initial_rr) << '\n';
  return std::isfinite(solution_error) && row_error < 1.0e-9 && constant_error < 1.0e-12 &&
      std::abs(flux_sum) < 1.0e-12 * std::max(1.0, flux_l1) && true_residual < 1.0e-20 * initial_rr;
}
}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  bool passed = true;
  if (argc == 2 && std::string_view(argv[1]) == "--geometry-separation") {
    // Independent experiment, not a replacement for the existing MMS gates.
    // First freeze the STL while refining h, then freeze h and vary the STL.
    double previous = 0.0;
    for (int n : {24, 48, 96}) {
      double error = 0.0;
      if (!run(2, n, error, 1024)) { passed = false; break; }
      if (previous > 0.0)
        std::cout << "IBM_SCALAR_FIXED_STL_ORDER facets=1024 n=" << n
                  << " order=" << std::log(previous / error) / std::log(2.0) << '\n';
      previous = error;
    }
    for (int facets : {48, 96, 192}) {
      double error = 0.0;
      if (!run(2, 48, error, facets)) { passed = false; break; }
    }
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc != 1) { MPI_Finalize(); return 2; }
  for (int shape = 0; shape < 3; ++shape) {
    double previous = 0.0;
    for (int n : {24, 48, 96}) {
      double error = 0.0;
      if (!run(shape, n, error)) { passed = false; break; }
      if (previous > 0.0) {
        const double order = std::log(previous / error) / std::log(2.0);
        std::cout << "IBM_SCALAR_ORDER shape=" << shape << " n=" << n << " order=" << order << '\n';
        // Only the grid-aligned flat-wall second-order property is a gate.
        // Oblique/curved error is measured, never labeled second order just
        // because cut-face inventory is conservative.
        if (shape == 0 && order < 1.8) passed = false;
      }
      previous = error;
    }
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
