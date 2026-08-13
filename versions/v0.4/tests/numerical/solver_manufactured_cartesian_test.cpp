// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_execution.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr double kPi = 3.141592653589793238462643383279502884;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

CartesianMeshSpec mesh_spec(std::int32_t n, bool stretched) {
  CartesianMeshSpec mesh;
  mesh.kind = stretched ? GeometryKind::tensor_stretched
                        : GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  // Refine one smooth-family monitor: the permitted adjacent ratio tends to
  // one with h, so the 16/32/64 sequence samples a smooth tensor mapping
  // instead of preserving a fixed grid-spacing jump.
  mesh.max_growth_ratio =
      stretched ? 1.0 + 0.8 / static_cast<double>(n) : 1.0;
  if (stretched) {
    const double inverse_n = 1.0 / static_cast<double>(n);
    mesh.has_base_spacing = true;
    mesh.base_spacing = {1.08 * inverse_n, 1.08 * inverse_n,
                         1.08 * inverse_n};
    mesh.minimum_spacing = {0.92 * inverse_n, 0.92 * inverse_n,
                            0.92 * inverse_n};
    mesh.focus_regions.push_back(
        {{0.35, 0.35, 0.35}, {0.65, 0.65, 0.65},
         {0.96 * inverse_n, 0.96 * inverse_n, 0.96 * inverse_n}});
  } else {
    mesh.minimum_spacing = {0.25 / n, 0.25 / n, 0.25 / n};
  }
  mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(n) * static_cast<std::uint64_t>(n) *
      static_cast<std::uint64_t>(n);
  mesh.limits.max_memory_bytes_per_rank = 1U << 30U;
  return mesh;
}

ValidatedModel periodic_model(
    const CartesianMeshSpec& mesh,
    ConvectionScheme convection_scheme = ConvectionScheme::central2) {
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x94c17e32U;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& face : model.boundaries) {
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
    std::int32_t n, bool stretched, KernelFixture& fixture,
    ConvectionScheme convection_scheme = ConvectionScheme::central2) {
  const CartesianMeshSpec mesh = mesh_spec(n, stretched);
  FieldRegistry registry;
  TimeSchemePlan time;
  return CartesianGeometryCompiler::compile(
             MPI_COMM_SELF, mesh, GeometryBudget{}, fixture.geometry,
             fixture.patch) &&
         BoundaryCompiler::compile(
             MPI_COMM_SELF, periodic_model(mesh, convection_scheme), fixture.geometry,
             fixture.patch, registry, fixture.boundary, fixture.schemes,
             time) &&
         CartesianKernelPlan::compile(
             fixture.schemes, fixture.geometry, fixture.patch,
             fixture.boundary, fixture.kernels);
}

struct OwnedField {
  std::vector<double> allocation;
  FieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision) {
  OwnedField field;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field.allocation.assign(nx * ny * nz * components, 0.0);
  field.view.base = field.allocation.data() + ghosts + ghosts * nx +
                    ghosts * nx * ny;
  field.view.interior = cells;
  field.view.ghosts = {ghosts, ghosts, ghosts};
  field.view.components = components;
  field.view.stride_y = nx;
  field.view.stride_z = nx * ny;
  field.view.component_stride = nx * ny * nz;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = static_cast<StorageIdentity>(id) + 1001U;
  field.view.revision_domain = 7001U;
  return field;
}

double extended_centre(const AxisMetrics& axis, std::int32_t index) {
  const Span<const double> centres = axis.centres();
  const Span<const double> widths = axis.widths();
  if (index < 0) {
    return centres.data[0U] + static_cast<double>(index) * widths.data[0U];
  }
  if (static_cast<std::size_t>(index) >= centres.size) {
    const auto last = static_cast<std::int32_t>(centres.size - 1U);
    return centres.data[centres.size - 1U] +
           static_cast<double>(index - last) * widths.data[widths.size - 1U];
  }
  return centres.data[static_cast<std::size_t>(index)];
}

double q_value(double x, double y, double z) {
  return 1.0 + 0.20 * std::sin(2.0 * kPi * x) +
         0.15 * std::cos(2.0 * kPi * y) +
         0.10 * std::sin(2.0 * kPi * z);
}

Real3 q_gradient(double x, double y, double z) {
  return {0.40 * kPi * std::cos(2.0 * kPi * x),
          -0.30 * kPi * std::sin(2.0 * kPi * y),
          0.20 * kPi * std::cos(2.0 * kPi * z)};
}

double q_laplacian(double x, double y, double z) {
  const double wave = 2.0 * kPi;
  return -wave * wave *
         (0.20 * std::sin(wave * x) + 0.15 * std::cos(wave * y) +
          0.10 * std::sin(wave * z));
}

double gamma_value(double x, double y, double z) {
  return 1.5 + 0.25 * x + 0.15 * y + 0.10 * z;
}

double variable_diffusion_oracle(double x, double y, double z) {
  const Real3 gradient = q_gradient(x, y, z);
  return gamma_value(x, y, z) * q_laplacian(x, y, z) +
         0.25 * gradient.x + 0.15 * gradient.y + 0.10 * gradient.z;
}

double monotone_q(double x, double y, double z) {
  return std::exp(0.32 * x + 0.23 * y + 0.17 * z);
}

Real3 monotone_q_gradient(double x, double y, double z) {
  const double q = monotone_q(x, y, z);
  return {0.32 * q, 0.23 * q, 0.17 * q};
}

Real3 divergent_velocity(double x, double y, double z) {
  return {std::sin(2.0 * kPi * x), 0.5 * std::cos(2.0 * kPi * y),
          0.25 * std::sin(2.0 * kPi * z)};
}

double velocity_divergence(double x, double y, double z) {
  return 2.0 * kPi * std::cos(2.0 * kPi * x) -
         kPi * std::sin(2.0 * kPi * y) +
         0.5 * kPi * std::cos(2.0 * kPi * z);
}

template <class Callback>
void fill_cells(FieldView field, const CartesianGeometryPlan& geometry,
                const MeshPatch& patch, std::uint8_t ghosts,
                Callback&& callback) {
  for (std::int32_t z = -ghosts; z < patch.cells.z + ghosts; ++z) {
    const double zc = extended_centre(geometry.z(), patch.begin.z + z);
    for (std::int32_t y = -ghosts; y < patch.cells.y + ghosts; ++y) {
      const double yc = extended_centre(geometry.y(), patch.begin.y + y);
      for (std::int32_t x = -ghosts; x < patch.cells.x + ghosts; ++x) {
        const double xc = extended_centre(geometry.x(), patch.begin.x + x);
        callback(field, Int3{x, y, z}, xc, yc, zc);
      }
    }
  }
}

double volume(const CartesianGeometryPlan& geometry, const MeshPatch& patch,
              Int3 local) {
  return geometry.x().widths().data[patch.begin.x + local.x] *
         geometry.y().widths().data[patch.begin.y + local.y] *
         geometry.z().widths().data[patch.begin.z + local.z];
}

template <class Oracle>
double l2_error(ConstFieldView field, std::uint8_t component,
                const CartesianGeometryPlan& geometry, const MeshPatch& patch,
                Oracle&& oracle) {
  long double error = 0.0L;
  long double measure = 0.0L;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    const double zc = geometry.z().centres().data[patch.begin.z + z];
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      const double yc = geometry.y().centres().data[patch.begin.y + y];
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const double xc = geometry.x().centres().data[patch.begin.x + x];
        const Int3 cell{x, y, z};
        const long double weight = volume(geometry, patch, cell);
        const long double difference =
            field.unchecked(cell, component) - oracle(xc, yc, zc);
        error += difference * difference * weight;
        measure += weight;
      }
    }
  }
  return std::sqrt(static_cast<double>(error / measure));
}

template <class Oracle>
double interior_l2_error(ConstFieldView field, std::uint8_t component,
                         const CartesianGeometryPlan& geometry,
                         const MeshPatch& patch, std::int32_t trim,
                         Oracle&& oracle) {
  long double error = 0.0L;
  long double measure = 0.0L;
  const Int3 end{patch.cells.x - trim, patch.cells.y - trim,
                 patch.cells.z - trim};
  for (std::int32_t z = trim; z < end.z; ++z) {
    const double zc = geometry.z().centres().data[patch.begin.z + z];
    for (std::int32_t y = trim; y < end.y; ++y) {
      const double yc = geometry.y().centres().data[patch.begin.y + y];
      for (std::int32_t x = trim; x < end.x; ++x) {
        const double xc = geometry.x().centres().data[patch.begin.x + x];
        const Int3 cell{x, y, z};
        const long double weight = volume(geometry, patch, cell);
        const long double difference =
            field.unchecked(cell, component) - oracle(xc, yc, zc);
        error += difference * difference * weight;
        measure += weight;
      }
    }
  }
  return std::sqrt(static_cast<double>(error / measure));
}

void fill_constant_velocity_flux(FaceFluxView flux,
                                 const CartesianGeometryPlan& geometry,
                                 Real3 velocity) {
  for (std::int32_t z = 0; z < flux.x.extents.z; ++z) {
    const double dz = geometry.z().widths().data[z];
    for (std::int32_t y = 0; y < flux.x.extents.y; ++y) {
      const double area = geometry.y().widths().data[y] * dz;
      for (std::int32_t x = 0; x < flux.x.extents.x; ++x) {
        flux.x.unchecked({x, y, z}) = velocity.x * area;
      }
    }
  }
  for (std::int32_t z = 0; z < flux.y.extents.z; ++z) {
    const double dz = geometry.z().widths().data[z];
    for (std::int32_t y = 0; y < flux.y.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.y.extents.x; ++x) {
        flux.y.unchecked({x, y, z}) =
            velocity.y * geometry.x().widths().data[x] * dz;
      }
    }
  }
  for (std::int32_t z = 0; z < flux.z.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.z.extents.y; ++y) {
      const double dy = geometry.y().widths().data[y];
      for (std::int32_t x = 0; x < flux.z.extents.x; ++x) {
        flux.z.unchecked({x, y, z}) =
            velocity.z * geometry.x().widths().data[x] * dy;
      }
    }
  }
}

struct Errors {
  double gradient{};
  double divergence{};
  double convection{};
  double diffusion{};
};

bool run_manufactured(std::int32_t n, bool stretched, Errors& errors) {
  KernelFixture fixture;
  if (!expect(make_fixture(n, stretched, fixture),
              "manufactured kernel fixture compiles")) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  constexpr std::uint8_t ghosts = 2U;
  OwnedField q = make_field(20U, cells, 1U, ghosts, 11U);
  OwnedField gradient = make_field(21U, cells, 3U, 0U, 12U);
  OwnedField rho = make_field(22U, cells, 1U, ghosts, 13U);
  OwnedField velocity = make_field(23U, cells, 3U, ghosts, 14U);
  OwnedField divergence = make_field(24U, cells, 1U, 0U, 15U);
  OwnedField convection = make_field(25U, cells, 1U, 0U, 16U);
  OwnedField gamma = make_field(26U, cells, 1U, ghosts, 17U);
  OwnedField diffusion = make_field(27U, cells, 1U, 0U, 18U);

  fill_cells(q.view, fixture.geometry, fixture.patch, ghosts,
             [](FieldView field, Int3 cell, double x, double y, double z) {
               field.unchecked(cell, 0U) = q_value(x, y, z);
             });
  fill_cells(rho.view, fixture.geometry, fixture.patch, ghosts,
             [](FieldView field, Int3 cell, double, double, double) {
               field.unchecked(cell, 0U) = 1.0;
             });
  fill_cells(gamma.view, fixture.geometry, fixture.patch, ghosts,
             [](FieldView field, Int3 cell, double x, double y, double z) {
               field.unchecked(cell, 0U) = gamma_value(x, y, z);
             });
  fill_cells(velocity.view, fixture.geometry, fixture.patch, ghosts,
             [](FieldView field, Int3 cell, double x, double y, double z) {
               const Real3 value = divergent_velocity(x, y, z);
               field.unchecked(cell, 0U) = value.x;
               field.unchecked(cell, 1U) = value.y;
               field.unchecked(cell, 2U) = value.z;
             });

  const KernelBox all{{0, 0, 0}, cells};
  const std::array<ConstFieldView, 1U> q_read{as_const(q.view)};
  const std::array<FieldView, 1U> gradient_write{gradient.view};
  const KernelInvocation gradient_call{
      {q_read.data(), q_read.size()}, {gradient_write.data(), 1U}, all,
      0U, 0U, 1U, 0U, nullptr};
  bool passed = expect(static_cast<bool>(cartesian_gradient(
                           fixture.kernels, gradient_call)),
                       "manufactured gradient evaluates");

  FaceFluxStorage face_storage;
  FaceFluxView flux;
  passed &= expect(FaceFluxStorage::allocate_workspace(cells, 1U,
                                                       face_storage) &&
                       face_storage.workspace_view(0U, 51U, flux),
                   "manufactured face arena allocates");
  const std::array<ConstFieldView, 2U> mass_reads{as_const(rho.view),
                                                 as_const(velocity.view)};
  const KernelInvocation mass_call{{mass_reads.data(), mass_reads.size()}, {},
                                   all, 0U, 0U, 1U, 0U, nullptr};
  passed &= expect(static_cast<bool>(reconstruct_mass_flux(
                       fixture.kernels, mass_call, flux)),
                   "manufactured mass flux reconstructs");
  const std::array<FieldView, 1U> divergence_write{divergence.view};
  const KernelInvocation divergence_call{
      {}, {divergence_write.data(), 1U}, all, 0U, 0U, 1U, flux.revision,
      nullptr};
  passed &= expect(static_cast<bool>(cartesian_provisional_face_divergence(
                       fixture.kernels, as_const(flux), divergence_call)),
                   "manufactured divergence evaluates");

  constexpr Real3 transport_velocity{0.7, -0.4, 0.25};
  fill_constant_velocity_flux(flux, fixture.geometry, transport_velocity);
  const std::array<FieldView, 1U> convection_write{convection.view};
  const KernelInvocation convection_call{
      {q_read.data(), 1U}, {convection_write.data(), 1U}, all, 0U, 0U, 1U,
      flux.revision, nullptr};
  passed &= expect(static_cast<bool>(cartesian_provisional_convection(
                       fixture.kernels, ConvectionScheme::central2,
                       as_const(flux), convection_call)),
                   "manufactured convection evaluates");

  const std::array<FieldView, 1U> diffusion_write{diffusion.view};
  const KernelInvocation diffusion_call{
      {q_read.data(), 1U}, {diffusion_write.data(), 1U}, all, 0U, 0U, 1U,
      0U, nullptr};
  passed &= expect(static_cast<bool>(cartesian_diffusion(
                       fixture.kernels, as_const(gamma.view), diffusion_call)),
                   "manufactured diffusion evaluates");
  if (!passed) {
    return false;
  }

  const double ex = l2_error(
      as_const(gradient.view), 0U, fixture.geometry, fixture.patch,
      [](double x, double y, double z) { return q_gradient(x, y, z).x; });
  const double ey = l2_error(
      as_const(gradient.view), 1U, fixture.geometry, fixture.patch,
      [](double x, double y, double z) { return q_gradient(x, y, z).y; });
  const double ez = l2_error(
      as_const(gradient.view), 2U, fixture.geometry, fixture.patch,
      [](double x, double y, double z) { return q_gradient(x, y, z).z; });
  errors.gradient = std::sqrt(ex * ex + ey * ey + ez * ez);
  errors.divergence = l2_error(
      as_const(divergence.view), 0U, fixture.geometry, fixture.patch,
      [](double x, double y, double z) {
        return velocity_divergence(x, y, z);
      });
  errors.convection = l2_error(
      as_const(convection.view), 0U, fixture.geometry, fixture.patch,
      [](double x, double y, double z) {
        const Real3 gradient_value = q_gradient(x, y, z);
        return transport_velocity.x * gradient_value.x +
               transport_velocity.y * gradient_value.y +
               transport_velocity.z * gradient_value.z;
      });
  errors.diffusion = l2_error(
      as_const(diffusion.view), 0U, fixture.geometry, fixture.patch,
      [](double x, double y, double z) {
        return variable_diffusion_oracle(x, y, z);
      });
  return std::isfinite(errors.gradient) && std::isfinite(errors.divergence) &&
         std::isfinite(errors.convection) && std::isfinite(errors.diffusion);
}

double order(double coarse, double fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

bool check_orders(bool stretched) {
  std::array<Errors, 3U> values{};
  const std::array<std::int32_t, 3U> resolutions{16, 32, 64};
  bool passed = true;
  for (std::size_t level = 0U; level < resolutions.size(); ++level) {
    passed &= run_manufactured(resolutions[level], stretched, values[level]);
  }
  const auto gate = [&](double e16, double e32, double e64,
                        std::string_view description) {
    const bool above_roundoff = e16 > 1.0e-12 && e32 > 1.0e-13 &&
                                e64 > 1.0e-14;
    const double first_order = order(e16, e32);
    const double second_order = order(e32, e64);
    const bool accepted = above_roundoff && first_order >= 1.8 &&
                          second_order >= 1.8;
    if (!accepted) {
      std::cerr << description << ": errors=" << e16 << ',' << e32 << ','
                << e64 << " orders=" << first_order << ',' << second_order
                << '\n';
    }
    return expect(accepted, description);
  };
  passed &= gate(values[0U].gradient, values[1U].gradient,
                 values[2U].gradient, "gradient retains second order");
  passed &= gate(values[0U].divergence, values[1U].divergence,
                 values[2U].divergence, "divergence retains second order");
  passed &= gate(values[0U].convection, values[1U].convection,
                 values[2U].convection, "convection retains second order");
  passed &= gate(values[0U].diffusion, values[1U].diffusion,
                 values[2U].diffusion, "diffusion retains second order");
  return passed;
}

bool run_limited_convection(std::int32_t n, bool stretched,
                            ConvectionScheme scheme, double& error) {
  KernelFixture fixture;
  if (!expect(make_fixture(n, stretched, fixture, scheme),
              "limited convection fixture compiles")) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  constexpr std::uint8_t ghosts = 2U;
  OwnedField transported = make_field(40U, cells, 1U, ghosts, 41U);
  OwnedField convection = make_field(41U, cells, 1U, 0U, 42U);
  fill_cells(
      transported.view, fixture.geometry, fixture.patch, ghosts,
      [](FieldView field, Int3 cell, double x, double y, double z) {
        field.unchecked(cell, 0U) = monotone_q(x, y, z);
      });

  FaceFluxStorage face_storage;
  FaceFluxView flux;
  bool passed = expect(
      FaceFluxStorage::allocate_workspace(cells, 1U, face_storage) &&
          face_storage.workspace_view(0U, 71U, flux),
      "limited convection face arena allocates");
  constexpr Real3 transport_velocity{0.7, -0.4, 0.25};
  fill_constant_velocity_flux(flux, fixture.geometry, transport_velocity);
  const std::array<ConstFieldView, 1U> reads{as_const(transported.view)};
  const std::array<FieldView, 1U> writes{convection.view};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {writes.data(), writes.size()},
      {{0, 0, 0}, cells}, 0U, 0U, 1U, flux.revision, nullptr};
  passed &= expect(static_cast<bool>(cartesian_provisional_convection(
                       fixture.kernels, scheme, as_const(flux), invocation)),
                   "limited convection evaluates");
  if (!passed) {
    return false;
  }
  error = interior_l2_error(
      as_const(convection.view), 0U, fixture.geometry, fixture.patch, 2,
      [](double x, double y, double z) {
        constexpr Real3 velocity{0.7, -0.4, 0.25};
        const Real3 gradient = monotone_q_gradient(x, y, z);
        return velocity.x * gradient.x + velocity.y * gradient.y +
               velocity.z * gradient.z;
      });
  return std::isfinite(error);
}

bool check_limited_convection_orders(bool stretched,
                                     ConvectionScheme scheme,
                                     std::string_view description) {
  const std::array<std::int32_t, 3U> resolutions{16, 32, 64};
  std::array<double, 3U> errors{};
  bool passed = true;
  for (std::size_t level = 0U; level < resolutions.size(); ++level) {
    passed &= run_limited_convection(resolutions[level], stretched, scheme,
                                     errors[level]);
  }
  const double order_16_32 = order(errors[0U], errors[1U]);
  const double order_32_64 = order(errors[1U], errors[2U]);
  const bool accepted = errors[0U] > 1.0e-12 && errors[1U] > 1.0e-13 &&
                        errors[2U] > 1.0e-14 && order_16_32 >= 1.8 &&
                        order_32_64 >= 1.8;
  if (!accepted) {
    std::cerr << description << ": errors=" << errors[0U] << ','
              << errors[1U] << ',' << errors[2U] << " orders="
              << order_16_32 << ',' << order_32_64 << '\n';
  }
  passed &= expect(accepted, description);
  return passed;
}

bool same_counters(const KernelCounters& left,
                   const KernelCounters& right) {
  return left.invocations == right.invocations && left.cells == right.cells &&
         left.faces == right.faces &&
         left.logical_bytes_read == right.logical_bytes_read &&
         left.logical_bytes_written == right.logical_bytes_written;
}

bool test_noncentral_rejects_one_ghost(ConvectionScheme scheme,
                                       std::string_view description) {
  KernelFixture fixture;
  if (!expect(make_fixture(16, false, fixture, scheme),
              "one-ghost mutation fixture compiles")) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField transported = make_field(50U, cells, 1U, 1U, 51U);
  OwnedField output = make_field(51U, cells, 1U, 0U, 52U);
  fill_cells(
      transported.view, fixture.geometry, fixture.patch, 1U,
      [](FieldView field, Int3 cell, double x, double y, double z) {
        field.unchecked(cell, 0U) = monotone_q(x, y, z);
      });
  std::fill(output.allocation.begin(), output.allocation.end(), -8125.25);
  const std::vector<double> output_before = output.allocation;

  FaceFluxStorage face_storage;
  FaceFluxView flux;
  bool passed = expect(
      FaceFluxStorage::allocate_workspace(cells, 1U, face_storage) &&
          face_storage.workspace_view(0U, 81U, flux),
      "one-ghost mutation face arena allocates");
  fill_constant_velocity_flux(flux, fixture.geometry, {0.7, -0.4, 0.25});
  KernelCounters counters{3U, 5U, 7U, 11U, 13U};
  const KernelCounters counters_before = counters;
  const std::array<ConstFieldView, 1U> reads{as_const(transported.view)};
  const std::array<FieldView, 1U> writes{output.view};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {writes.data(), writes.size()},
      {{0, 0, 0}, cells}, 0U, 0U, 1U, flux.revision, &counters};
  const Status result = cartesian_provisional_convection(fixture.kernels, scheme,
                                               as_const(flux), invocation);
  passed &= expect(result.code == StatusCode::invalid_plan, description);
  passed &= expect(output.allocation == output_before,
                   "one-ghost rejection leaves output unchanged");
  passed &= expect(same_counters(counters, counters_before),
                   "one-ghost rejection leaves counters unchanged");
  return passed;
}

bool test_stretched_affine_and_harmonic_mutation() {
  KernelFixture fixture;
  bool passed = expect(make_fixture(32, true, fixture),
                       "stretched mutation fixture compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  constexpr std::uint8_t ghosts = 2U;
  OwnedField affine = make_field(60U, cells, 1U, ghosts, 61U);
  OwnedField gradient = make_field(61U, cells, 3U, 0U, 62U);
  fill_cells(affine.view, fixture.geometry, fixture.patch, ghosts,
             [](FieldView field, Int3 cell, double x, double y, double z) {
               field.unchecked(cell, 0U) = 1.0 + 2.0 * x - 3.0 * y + 4.0 * z;
             });
  const std::array<ConstFieldView, 1U> reads{as_const(affine.view)};
  const std::array<FieldView, 1U> writes{gradient.view};
  const KernelInvocation gradient_call{{reads.data(), 1U}, {writes.data(), 1U},
                                       {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U,
                                       nullptr};
  passed &= expect(static_cast<bool>(cartesian_gradient(fixture.kernels,
                                                       gradient_call)),
                   "stretched affine gradient evaluates");
  double maximum_affine_error = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        maximum_affine_error =
            std::max(maximum_affine_error,
                     std::abs(gradient.view.unchecked(cell, 0U) - 2.0));
        maximum_affine_error =
            std::max(maximum_affine_error,
                     std::abs(gradient.view.unchecked(cell, 1U) + 3.0));
        maximum_affine_error =
            std::max(maximum_affine_error,
                     std::abs(gradient.view.unchecked(cell, 2U) - 4.0));
      }
    }
  }
  passed &= expect(maximum_affine_error <= 2.0e-12,
                   "stretched three-point gradient is affine exact");

  OwnedField q = make_field(62U, cells, 1U, ghosts, 63U);
  OwnedField gamma = make_field(63U, cells, 1U, ghosts, 64U);
  OwnedField output = make_field(64U, cells, 1U, 0U, 65U);
  const std::int32_t interface_cell = cells.x / 2 - 1;
  fill_cells(q.view, fixture.geometry, fixture.patch, ghosts,
             [](FieldView field, Int3 cell, double x, double, double) {
               field.unchecked(cell, 0U) = x * x;
             });
  fill_cells(gamma.view, fixture.geometry, fixture.patch, ghosts,
             [interface_cell](FieldView field, Int3 cell, double, double,
                              double) {
               field.unchecked(cell, 0U) =
                   cell.x <= interface_cell ? 1.0 : 9.0;
             });
  const std::array<ConstFieldView, 1U> q_reads{as_const(q.view)};
  const std::array<FieldView, 1U> output_writes{output.view};
  const KernelInvocation diffusion_call{
      {q_reads.data(), 1U}, {output_writes.data(), 1U}, {{0, 0, 0}, cells},
      0U, 0U, 1U, 0U, nullptr};
  passed &= expect(static_cast<bool>(cartesian_diffusion(
                       fixture.kernels, as_const(gamma.view), diffusion_call)),
                   "positive discontinuous diffusion remains finite");

  const std::int32_t i = interface_cell;
  const double x_l = fixture.geometry.x().centres().data[i];
  const double x_r = fixture.geometry.x().centres().data[i + 1];
  const double x_m = fixture.geometry.x().centres().data[i - 1];
  const double face_r = fixture.geometry.x().faces().data[i + 1];
  const double face_l = fixture.geometry.x().faces().data[i];
  const double right_harmonic =
      ((x_r * x_r) - (x_l * x_l)) /
      ((face_r - x_l) / 1.0 + (x_r - face_r) / 9.0);
  const double left_rate =
      ((x_l * x_l) - (x_m * x_m)) /
      ((face_l - x_m) / 1.0 + (x_l - face_l) / 1.0);
  const double expected = (right_harmonic - left_rate) /
                          fixture.geometry.x().widths().data[i];
  const double arithmetic_gamma =
      ((x_r - face_r) * 1.0 + (face_r - x_l) * 9.0) / (x_r - x_l);
  const double arithmetic_right =
      arithmetic_gamma * ((x_r * x_r) - (x_l * x_l)) / (x_r - x_l);
  const double arithmetic = (arithmetic_right - left_rate) /
                            fixture.geometry.x().widths().data[i];
  const double actual = output.view.unchecked({i, cells.y / 2, cells.z / 2},
                                               0U);
  const double tolerance = 2.0e-12 * std::max(1.0, std::abs(expected));
  passed &= expect(std::isfinite(actual) &&
                       std::abs(actual - expected) <= tolerance &&
                       std::abs(actual - arithmetic) > 100.0 * tolerance,
                   "diffusion uses distance-weighted harmonic face conductance");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = true;
  passed &= check_orders(false);
  passed &= check_orders(true);
  passed &= check_limited_convection_orders(
      false, ConvectionScheme::limited_central2,
      "uniform limited-central convection retains second order");
  passed &= check_limited_convection_orders(
      true, ConvectionScheme::limited_central2,
      "stretched limited-central convection retains second order");
  passed &= check_limited_convection_orders(
      false, ConvectionScheme::tvd2,
      "uniform TVD convection retains second order");
  passed &= check_limited_convection_orders(
      true, ConvectionScheme::tvd2,
      "stretched TVD convection retains second order");
  passed &= test_noncentral_rejects_one_ghost(
      ConvectionScheme::limited_central2,
      "limited-central rejects a one-ghost transported field");
  passed &= test_noncentral_rejects_one_ghost(
      ConvectionScheme::tvd2, "TVD rejects a one-ghost transported field");
  passed &= test_stretched_affine_and_harmonic_mutation();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 manufactured Cartesian kernel tests passed\n";
  return 0;
}
