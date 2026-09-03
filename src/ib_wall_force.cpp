// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_wall_force.hpp"

#include "ib_quadratic_reconstruction_detail.hpp"
#include "ib_wall_force_detail.hpp"

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::immersed {
namespace {

using runtime::Real3;

thread_local const std::vector<detail::WallPressureNormalGradient>
    *scoped_wall_pressure_normal_gradients = nullptr;

#ifdef HUNDUN_IMMERSED_ENABLE_TEST_ACCESS
thread_local std::vector<detail::WallForcePointSnapshot>
    *scoped_wall_force_point_snapshots = nullptr;

class WallForceTraceScope final {
public:
  explicit WallForceTraceScope(
      std::vector<detail::WallForcePointSnapshot> &snapshots) {
    if (scoped_wall_force_point_snapshots != nullptr)
      throw runtime::Error("wall force trace scope is nested");
    scoped_wall_force_point_snapshots = &snapshots;
  }
  ~WallForceTraceScope() noexcept { scoped_wall_force_point_snapshots = nullptr; }
};
#endif

class WallPressureAuthorityScope final {
public:
  explicit WallPressureAuthorityScope(
      const std::vector<detail::WallPressureNormalGradient> &gradients) {
    if (scoped_wall_pressure_normal_gradients != nullptr)
      throw runtime::Error("wall force pressure-authority scope is nested");
    scoped_wall_pressure_normal_gradients = &gradients;
  }
  ~WallPressureAuthorityScope() noexcept {
    scoped_wall_pressure_normal_gradients = nullptr;
  }
};

bool same(runtime::Int3 lhs, runtime::Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 add(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 multiply(double factor, Real3 value) noexcept {
  return {factor * value.x, factor * value.y, factor * value.z};
}

double dot(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Real3 cross(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

void require_finite(Real3 value, const char *message) {
  if (!finite(value)) {
    throw runtime::Error(message);
  }
}

void add_checked(Real3 &target, Real3 increment, const char *message) {
  const auto candidate = add(target, increment);
  require_finite(candidate, message);
  target = candidate;
}

void check_mpi(int code, const char *operation) {
  if (code == MPI_SUCCESS) {
    return;
  }
  char text[MPI_MAX_ERROR_STRING]{};
  int length = 0;
  const int text_code = MPI_Error_string(code, text, &length);
  if (text_code != MPI_SUCCESS || length < 0 || length > MPI_MAX_ERROR_STRING) {
    throw runtime::Error(std::string(operation) + ": MPI operation failed");
  }
  throw runtime::Error(std::string(operation) + ": " +
                       std::string(text, static_cast<std::size_t>(length)));
}

void reduce_real3(Real3 &value, MPI_Comm communicator, const char *operation) {
  std::array<double, 3> buffer{value.x, value.y, value.z};
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, buffer.data(),
                          static_cast<int>(buffer.size()), MPI_DOUBLE, MPI_SUM,
                          communicator),
            operation);
  value = {buffer[0], buffer[1], buffer[2]};
}

void validate_layouts(const runtime::FieldView<const double> &pressure,
                      const runtime::FieldView<const double> &velocity,
                      const runtime::FieldView<const double> &gradient,
                      const runtime::FieldView<const double> &viscosity) {
  if (pressure.components() != 1U || velocity.components() != 3U ||
      gradient.components() != 9U || viscosity.components() != 1U) {
    throw runtime::Error("wall force field component layout is invalid");
  }
  const auto extent = pressure.interior_extent();
  const int ghost_width = pressure.ghost_width();
  if (!same(velocity.interior_extent(), extent) ||
      !same(gradient.interior_extent(), extent) ||
      !same(viscosity.interior_extent(), extent) ||
      velocity.ghost_width() != ghost_width ||
      gradient.ghost_width() != ghost_width ||
      viscosity.ghost_width() != ghost_width) {
    throw runtime::Error("wall force field spatial layouts do not match");
  }
}

struct LocalSample final {
  ForceComponents force;
  MomentComponents moment;
  Real3 area_closure{};
  std::uint64_t count{};
};

Real3 viscous_traction(const std::array<double, 9> &gradient, double viscosity,
                       Real3 normal) {
  const double divergence = gradient[0] + gradient[4] + gradient[8];
  const std::array<double, 3> n{normal.x, normal.y, normal.z};
  std::array<double, 3> traction{};
  for (std::size_t row = 0U; row < 3U; ++row) {
    for (std::size_t column = 0U; column < 3U; ++column) {
      const double strain = gradient[row * 3U + column] +
                            gradient[column * 3U + row] -
                            (row == column ? (2.0 / 3.0) * divergence : 0.0);
      const double candidate = traction[row] + viscosity * strain * n[column];
      if (!std::isfinite(candidate)) {
        throw runtime::Error("wall force viscous traction is non-finite");
      }
      traction[row] = candidate;
    }
  }
  return {traction[0], traction[1], traction[2]};
}

LocalSample
integrate_local(const WallQuadraturePlan &plan, int rank,
                const std::vector<std::size_t> &local_order,
                const std::vector<detail::WallPressureNormalGradient>
                    *wall_pressure_normal_gradients,
                const runtime::FieldView<const double> &pressure,
                const runtime::FieldView<const double> &velocity,
                const runtime::FieldView<const double> &gradient_field,
                const runtime::FieldView<const double> &viscosity_field) {
  validate_layouts(pressure, velocity, gradient_field, viscosity_field);
  LocalSample result;
  const auto &points = plan.local_points();
  if (local_order.size() != points.size()) {
    throw runtime::Error("wall force canonical ordering is invalid");
  }
  for (const auto index : local_order) {
    if (index >= points.size() || points[index].owner_rank != rank) {
      throw runtime::Error("wall force quadrature ownership is invalid");
    }
    const auto &point = points[index];
    const auto pressure_link =
        detail::boundary_authority_link(point.reconstruction);
    const auto pressure_value_reconstruction =
        wall_pressure_normal_gradients == nullptr
            ? point.reconstruction
            : detail::boundary_authority_reconstruction(point.reconstruction);
    const auto &velocity_reconstruction = point.reconstruction;
    if (!finite(point.position_m) || !finite(point.solid_to_fluid_normal) ||
        !std::isfinite(point.weight_m2) || point.weight_m2 <= 0.0) {
      throw runtime::Error("wall force quadrature geometry is invalid");
    }
    const double normal_norm2 =
        dot(point.solid_to_fluid_normal, point.solid_to_fluid_normal);
    const double normal_tolerance =
        1024.0 * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(normal_norm2) ||
        std::abs(normal_norm2 - 1.0) > normal_tolerance) {
      throw runtime::Error("wall force quadrature normal is not unit length");
    }

    double pressure_normal_gradient = 0.0;
    if (wall_pressure_normal_gradients == nullptr) {
      const auto independent_pressure_gradient =
          point.reconstruction.gradient(point.position_m, pressure, 0U);
      pressure_normal_gradient =
          dot(independent_pressure_gradient, point.solid_to_fluid_normal);
    } else {
      const auto found =
          std::lower_bound(wall_pressure_normal_gradients->begin(),
                           wall_pressure_normal_gradients->end(), pressure_link,
                           [](const auto &candidate, ImmersedLinkId link) {
                             return candidate.link < link;
                           });
      if (found == wall_pressure_normal_gradients->end() ||
          found->link != pressure_link || !std::isfinite(found->value))
        throw runtime::Error(
            "wall force pressure-gradient authority is unavailable");
      pressure_normal_gradient = found->value;
    }
    // Both paths are quadratic-exact. The supplied datum is constrained along
    // its associated link normal; the resulting scalar pressure is still
    // integrated against the real surface normal below.
    const double pi = detail::value_with_origin_normal_gradient(
        pressure_value_reconstruction, point.position_m, pressure, 0U,
        pressure_normal_gradient);
    std::array<double, 3> reconstructed_velocity{};
    for (std::size_t component = 0U; component < 3U; ++component) {
      reconstructed_velocity[component] =
          velocity_reconstruction.value(point.position_m, velocity, component);
      if (!std::isfinite(reconstructed_velocity[component])) {
        throw runtime::Error("wall force reconstructed velocity is non-finite");
      }
    }
    std::array<double, 9> gradient{};
    for (std::size_t component = 0U; component < gradient.size(); ++component) {
      const double provided_gradient = velocity_reconstruction.value(
          point.position_m, gradient_field, component);
      if (!std::isfinite(provided_gradient)) {
        throw runtime::Error(
            "wall force reconstructed velocity gradient is non-finite");
      }
    }
    for (std::size_t component = 0U; component < 3U; ++component) {
      const auto direct_gradient = detail::gradient_with_origin_constraint(
          velocity_reconstruction, point.position_m, velocity, component, 0.0);
      const double normal_derivative =
          dot(direct_gradient, point.solid_to_fluid_normal);
      gradient[component * 3U] =
          normal_derivative * point.solid_to_fluid_normal.x;
      gradient[component * 3U + 1U] =
          normal_derivative * point.solid_to_fluid_normal.y;
      gradient[component * 3U + 2U] =
          normal_derivative * point.solid_to_fluid_normal.z;
    }
    const double viscosity =
        point.reconstruction.value(point.position_m, viscosity_field, 0U);
    if (!std::isfinite(pi)) {
      throw runtime::Error(
          "wall force reconstructed mechanical pressure is non-finite");
    }
    if (!std::isfinite(viscosity) || viscosity < 0.0) {
      throw runtime::Error(
          "wall force reconstructed effective viscosity is invalid");
    }

    const auto pressure_force =
        multiply(-pi * point.weight_m2, point.solid_to_fluid_normal);
    const auto viscous_force = multiply(
        point.weight_m2,
        viscous_traction(gradient, viscosity, point.solid_to_fluid_normal));
    const auto total_force = add(pressure_force, viscous_force);
    require_finite(pressure_force,
                   "wall force pressure contribution is non-finite");
    require_finite(viscous_force,
                   "wall force viscous contribution is non-finite");
    require_finite(total_force, "wall force total contribution is non-finite");

#ifdef HUNDUN_IMMERSED_ENABLE_TEST_ACCESS
    if (scoped_wall_force_point_snapshots != nullptr)
      scoped_wall_force_point_snapshots->push_back(
          {pressure_link, point.triangle, point.point_index, point.position_m,
           point.solid_to_fluid_normal, point.weight_m2, pressure_force,
           viscous_force, total_force});
#endif

    add_checked(result.force.pressure_N, pressure_force,
                "wall force pressure accumulation is non-finite");
    add_checked(result.force.viscous_N, viscous_force,
                "wall force viscous accumulation is non-finite");
    add_checked(result.force.total_N, total_force,
                "wall force total accumulation is non-finite");
    add_checked(result.moment.pressure_N_m,
                cross(point.position_m, pressure_force),
                "wall force pressure moment is non-finite");
    add_checked(result.moment.viscous_N_m,
                cross(point.position_m, viscous_force),
                "wall force viscous moment is non-finite");
    add_checked(result.moment.total_N_m, cross(point.position_m, total_force),
                "wall force total moment is non-finite");
    add_checked(result.area_closure,
                multiply(point.weight_m2, point.solid_to_fluid_normal),
                "wall force area closure is non-finite");
    if (result.count == std::numeric_limits<std::uint64_t>::max()) {
      throw runtime::Error("wall force quadrature count would overflow");
    }
    ++result.count;
  }
  return result;
}

bool valid_reduced(const WallForceSample &sample) noexcept {
  return finite(sample.surface_traction.pressure_N) &&
         finite(sample.surface_traction.viscous_N) &&
         finite(sample.surface_traction.total_N) &&
         finite(sample.moment_about_global_origin.pressure_N_m) &&
         finite(sample.moment_about_global_origin.viscous_N_m) &&
         finite(sample.moment_about_global_origin.total_N_m) &&
         finite(sample.area_vector_closure_m2);
}

std::vector<detail::BoundaryAuthorityOwner>
collective_authority_catalog(const WallQuadraturePlan &plan,
                             const runtime::MpiContext &mpi) {
  int local_root = plan.local_points().empty() ? mpi.size() : mpi.rank();
  int root = mpi.size();
  check_mpi(MPI_Allreduce(&local_root, &root, 1, MPI_INT, MPI_MIN, mpi.comm()),
            "MPI_Allreduce wall pressure-authority catalog root");
  if (root < 0 || root >= mpi.size())
    throw runtime::Error("wall force pressure-authority catalog is empty");
  int count = 0;
  const std::vector<detail::BoundaryAuthorityOwner> *root_catalog = nullptr;
  if (mpi.rank() == root) {
    root_catalog = &detail::boundary_authority_catalog(
        plan.local_points().front().reconstruction);
    count = static_cast<int>(root_catalog->size());
  }
  check_mpi(MPI_Bcast(&count, 1, MPI_INT, root, mpi.comm()),
            "MPI_Bcast wall pressure-authority catalog count");
  if (count <= 0)
    throw runtime::Error("wall force pressure-authority catalog is invalid");
  std::vector<std::uint64_t> links(static_cast<std::size_t>(count));
  std::vector<int> owners(static_cast<std::size_t>(count));
  if (mpi.rank() == root) {
    for (int index = 0; index < count; ++index) {
      links[static_cast<std::size_t>(index)] =
          (*root_catalog)[static_cast<std::size_t>(index)].link;
      owners[static_cast<std::size_t>(index)] =
          (*root_catalog)[static_cast<std::size_t>(index)].owner_rank;
    }
  }
  check_mpi(MPI_Bcast(links.data(), count, MPI_UINT64_T, root, mpi.comm()),
            "MPI_Bcast wall pressure-authority catalog links");
  check_mpi(MPI_Bcast(owners.data(), count, MPI_INT, root, mpi.comm()),
            "MPI_Bcast wall pressure-authority catalog owners");
  std::vector<detail::BoundaryAuthorityOwner> result(
      static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index)
    result[static_cast<std::size_t>(index)] = {
        links[static_cast<std::size_t>(index)],
        owners[static_cast<std::size_t>(index)]};
  return result;
}

} // namespace

WallForceIntegrator
WallForceIntegrator::create(const WallQuadraturePlan &plan,
                            const runtime::MpiContext &mpi) {
  bool local_ok = true;
  std::string local_message;
  std::vector<std::size_t> local_order;
  try {
    for (const auto &point : plan.local_points()) {
      if (point.owner_rank != mpi.rank()) {
        throw runtime::Error("wall force plan ownership is invalid");
      }
      if (!finite(point.position_m) || !finite(point.solid_to_fluid_normal) ||
          !std::isfinite(point.weight_m2) || point.weight_m2 <= 0.0) {
        throw runtime::Error("wall force plan geometry is invalid");
      }
      const auto link = detail::boundary_authority_link(point.reconstruction);
      const int owner =
          detail::boundary_authority_owner_rank(point.reconstruction);
      const auto &catalog =
          detail::boundary_authority_catalog(point.reconstruction);
      if (catalog.empty() ||
          catalog.size() >
              static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
          owner < 0 || owner >= mpi.size() ||
          !std::is_sorted(catalog.begin(), catalog.end(),
                          [](const auto &left, const auto &right) {
                            return left.link < right.link;
                          }))
        throw runtime::Error(
            "wall force pressure-authority catalog is invalid");
      const auto found =
          std::lower_bound(catalog.begin(), catalog.end(), link,
                           [](const auto &candidate, ImmersedLinkId target) {
                             return candidate.link < target;
                           });
      if (found == catalog.end() || found->link != link ||
          found->owner_rank != owner)
        throw runtime::Error(
            "wall force pressure-authority association is invalid");
    }
    local_order.resize(plan.local_points().size());
    std::iota(local_order.begin(), local_order.end(), std::size_t{0});
    const auto &points = plan.local_points();
    std::sort(local_order.begin(), local_order.end(),
              [&points](std::size_t lhs, std::size_t rhs) {
                const auto &left = points[lhs];
                const auto &right = points[rhs];
                const auto left_key = std::tie(
                    left.position_m.x, left.position_m.y, left.position_m.z,
                    left.solid_to_fluid_normal.x, left.solid_to_fluid_normal.y,
                    left.solid_to_fluid_normal.z, left.weight_m2,
                    left.reconstruction.quality().pivot_fingerprint,
                    left.point_index);
                const auto right_key =
                    std::tie(right.position_m.x, right.position_m.y,
                             right.position_m.z, right.solid_to_fluid_normal.x,
                             right.solid_to_fluid_normal.y,
                             right.solid_to_fluid_normal.z, right.weight_m2,
                             right.reconstruction.quality().pivot_fingerprint,
                             right.point_index);
                return left_key < right_key;
              });
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "wall force plan validation failed";
  }
  const auto status = runtime::collective_status(mpi, local_ok, local_message);
  if (!status.ok) {
    throw runtime::Error(status.message + " (lowest failing rank " +
                         std::to_string(status.failing_rank) + ")");
  }
  return WallForceIntegrator(plan, mpi, std::move(local_order));
}

WallForceSample WallForceIntegrator::integrate(
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FieldView<const double> &mu_eff_by_cell) const {
  const auto failed = [](int rank) {
    WallForceSample result;
    result.lowest_failing_rank = rank;
    return result;
  };
  bool local_layout_ok = true;
  std::string local_layout_message;
  try {
    validate_layouts(mechanical_pressure, velocity, velocity_gradient,
                     mu_eff_by_cell);
    const auto required_reach = plan_->maximum_halo_reach();
    if (required_reach == 0U ||
        mechanical_pressure.ghost_width() <
            static_cast<int>(required_reach) ||
        velocity.ghost_width() < static_cast<int>(required_reach) ||
        velocity_gradient.ghost_width() < static_cast<int>(required_reach) ||
        mu_eff_by_cell.ghost_width() < static_cast<int>(required_reach))
      throw runtime::Error("wall force field Halo width is insufficient");
  } catch (const std::exception &error) {
    local_layout_ok = false;
    local_layout_message = error.what();
  } catch (...) {
    local_layout_ok = false;
    local_layout_message = "wall force field layout validation failed";
  }
  const auto layout_status = runtime::collective_status(
      *mpi_, local_layout_ok, local_layout_message);
  if (!layout_status.ok)
    return failed(layout_status.failing_rank);

  const auto *wall_pressure_normal_gradients =
      scoped_wall_pressure_normal_gradients;
  std::vector<detail::WallPressureNormalGradient> authoritative;
  if (wall_pressure_normal_gradients != nullptr) {
    struct GradientWire final {
      std::uint64_t link{};
      double value{};
      int source_rank{};
    };
    static_assert(std::is_trivially_copyable_v<GradientWire>);
    const bool local_input_ok =
        wall_pressure_normal_gradients->size() <=
            static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                sizeof(GradientWire) &&
        std::all_of(
            wall_pressure_normal_gradients->begin(),
            wall_pressure_normal_gradients->end(),
            [](const auto &value) { return std::isfinite(value.value); });
    const auto input_status = runtime::collective_status(
        *mpi_, local_input_ok, "wall force pressure-gradient input is invalid");
    if (!input_status.ok)
      return failed(input_status.failing_rank);
    std::vector<GradientWire> local;
    local.reserve(wall_pressure_normal_gradients->size());
    for (const auto &gradient : *wall_pressure_normal_gradients)
      local.push_back({gradient.link, gradient.value, mpi_->rank()});
    const int local_bytes =
        static_cast<int>(local.size() * sizeof(GradientWire));
    std::vector<int> counts(static_cast<std::size_t>(mpi_->size()));
    check_mpi(MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                            mpi_->comm()),
              "MPI_Allgather wall pressure-gradient counts");
    std::vector<int> offsets(counts.size(), 0);
    std::size_t total_bytes = 0U;
    bool payload_ok = true;
    for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
      if (counts[rank] < 0 ||
          total_bytes >
              static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                  static_cast<std::size_t>(counts[rank])) {
        payload_ok = false;
        break;
      }
      offsets[rank] = static_cast<int>(total_bytes);
      total_bytes += static_cast<std::size_t>(counts[rank]);
    }
    if (!payload_ok || total_bytes % sizeof(GradientWire) != 0U)
      return failed(0);
    std::vector<GradientWire> gathered(total_bytes / sizeof(GradientWire));
    check_mpi(MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE,
                             gathered.data(), counts.data(), offsets.data(),
                             MPI_BYTE, mpi_->comm()),
              "MPI_Allgatherv wall pressure gradients");
    std::sort(gathered.begin(), gathered.end(),
              [](const auto &left, const auto &right) {
                return std::tie(left.link, left.source_rank) <
                       std::tie(right.link, right.source_rank);
              });
    int duplicate_rank = mpi_->size();
    for (std::size_t index = 1U; index < gathered.size(); ++index)
      if (gathered[index - 1U].link == gathered[index].link)
        duplicate_rank =
            std::min({duplicate_rank, gathered[index - 1U].source_rank,
                      gathered[index].source_rank});
    if (duplicate_rank < mpi_->size())
      return failed(duplicate_rank);

    const auto expected = collective_authority_catalog(*plan_, *mpi_);
    int catalog_failure_rank = mpi_->size();
    std::size_t gathered_index = 0U;
    std::size_t expected_index = 0U;
    while (gathered_index < gathered.size() ||
           expected_index < expected.size()) {
      if (expected_index == expected.size() ||
          (gathered_index < gathered.size() &&
           gathered[gathered_index].link < expected[expected_index].link)) {
        catalog_failure_rank = std::min(catalog_failure_rank,
                                        gathered[gathered_index].source_rank);
        ++gathered_index;
        continue;
      }
      if (gathered_index == gathered.size() ||
          expected[expected_index].link < gathered[gathered_index].link) {
        catalog_failure_rank =
            std::min(catalog_failure_rank, expected[expected_index].owner_rank);
        ++expected_index;
        continue;
      }
      if (gathered[gathered_index].source_rank !=
          expected[expected_index].owner_rank)
        catalog_failure_rank = std::min({catalog_failure_rank,
                                         gathered[gathered_index].source_rank,
                                         expected[expected_index].owner_rank});
      ++gathered_index;
      ++expected_index;
    }
    if (catalog_failure_rank < mpi_->size())
      return failed(catalog_failure_rank);

    authoritative.reserve(gathered.size());
    for (const auto &gradient : gathered)
      authoritative.push_back({gradient.link, gradient.value});
    wall_pressure_normal_gradients = &authoritative;
  }
  LocalSample local;
  bool local_ok = true;
  std::string local_message;
  try {
    local = integrate_local(*plan_, mpi_->rank(), local_order_,
                            wall_pressure_normal_gradients, mechanical_pressure,
                            velocity, velocity_gradient, mu_eff_by_cell);
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "wall force local integration failed";
  }
  const auto preparation =
      runtime::collective_status(*mpi_, local_ok, local_message);
  if (!preparation.ok) {
    WallForceSample failure;
    failure.lowest_failing_rank = preparation.failing_rank;
    return failure;
  }
  const auto communicator_size = static_cast<std::uint64_t>(mpi_->size());
  const bool count_safe =
      communicator_size > 0U &&
      local.count <=
          std::numeric_limits<std::uint64_t>::max() / communicator_size;
  const auto count_status = runtime::collective_status(
      *mpi_, count_safe, "wall force quadrature count sum would overflow");
  if (!count_status.ok) {
    WallForceSample failure;
    failure.lowest_failing_rank = count_status.failing_rank;
    return failure;
  }

  WallForceSample result;
  result.surface_traction = local.force;
  result.moment_about_global_origin = local.moment;
  result.area_vector_closure_m2 = local.area_closure;
  result.quadrature_point_count = local.count;
  const auto communicator = mpi_->comm();

  reduce_real3(result.surface_traction.pressure_N, communicator,
               "MPI_Allreduce wall pressure force");
  reduce_real3(result.surface_traction.viscous_N, communicator,
               "MPI_Allreduce wall viscous force");
  reduce_real3(result.surface_traction.total_N, communicator,
               "MPI_Allreduce wall total force");
  reduce_real3(result.moment_about_global_origin.pressure_N_m, communicator,
               "MPI_Allreduce wall pressure moment");
  reduce_real3(result.moment_about_global_origin.viscous_N_m, communicator,
               "MPI_Allreduce wall viscous moment");
  reduce_real3(result.moment_about_global_origin.total_N_m, communicator,
               "MPI_Allreduce wall total moment");
  reduce_real3(result.area_vector_closure_m2, communicator,
               "MPI_Allreduce wall area closure");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &result.quadrature_point_count, 1,
                          MPI_UINT64_T, MPI_SUM, communicator),
            "MPI_Allreduce wall quadrature count");

  const auto reduced_status = runtime::collective_status(
      *mpi_, valid_reduced(result), "wall force reduced result is non-finite");
  if (!reduced_status.ok) {
    WallForceSample failure;
    failure.lowest_failing_rank = reduced_status.failing_rank;
    return failure;
  }
  return result;
}

namespace detail {

#ifdef HUNDUN_IMMERSED_ENABLE_TEST_ACCESS
WallForceTrace trace_wall_force_for_test(
    const WallForceIntegrator &integrator,
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FieldView<const double> &mu_eff_by_cell) {
  WallForceTrace result;
  WallForceTraceScope scope(result.local_points);
  result.reduced = integrator.integrate(mechanical_pressure, velocity,
                                        velocity_gradient, mu_eff_by_cell);
  return result;
}
#endif

WallForceSample integrate_with_wall_pressure_authority(
    const WallForceIntegrator &integrator,
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FieldView<const double> &mu_eff_by_cell,
    const std::vector<WallPressureNormalGradient> &gradients) {
  WallPressureAuthorityScope scope(gradients);
  return integrator.integrate(mechanical_pressure, velocity, velocity_gradient,
                              mu_eff_by_cell);
}

} // namespace detail

} // namespace hundun::immersed
