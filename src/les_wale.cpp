// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/les_wale.hpp"

#include "hundun/rt_error.hpp"
#include "hundun/rt_kernel_field_view.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <utility>

namespace hundun::les {
namespace {

constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
    hash ^= (value >> (8U * byte)) & UINT64_C(0xff);
    hash *= kFnvPrime;
  }
  return hash;
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double wale_kinematic_viscosity(const std::array<double, 9> &gradient,
                                double coefficient,
                                double filter_width_m) {
  double scale = 0.0;
  for (const double value : gradient) {
    if (!std::isfinite(value))
      throw runtime::Error("WALE velocity gradient is non-finite");
    scale = std::max(scale, std::abs(value));
  }
  if (scale == 0.0)
    return 0.0;

  std::array<double, 9> g{};
  for (std::size_t index = 0U; index < g.size(); ++index)
    g[index] = gradient[index] / scale;
  const auto at = [](const auto &tensor, std::size_t i, std::size_t j) {
    return tensor[3U * i + j];
  };
  std::array<double, 9> strain{};
  std::array<double, 9> square{};
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      strain[3U * i + j] = 0.5 * (at(g, i, j) + at(g, j, i));
      for (std::size_t k = 0U; k < 3U; ++k)
        square[3U * i + j] += at(g, i, k) * at(g, k, j);
    }
  }
  const double trace = at(square, 0U, 0U) + at(square, 1U, 1U) +
                       at(square, 2U, 2U);
  std::array<double, 9> deviatoric_square{};
  double strain_invariant = 0.0;
  double square_invariant = 0.0;
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      const std::size_t index = 3U * i + j;
      deviatoric_square[index] =
          0.5 * (at(square, i, j) + at(square, j, i)) -
          (i == j ? trace / 3.0 : 0.0);
      strain_invariant += strain[index] * strain[index];
      square_invariant +=
          deviatoric_square[index] * deviatoric_square[index];
    }
  }
  if (strain_invariant == 0.0 && square_invariant == 0.0)
    return 0.0;
  const double numerator =
      square_invariant * std::sqrt(square_invariant);
  const double denominator =
      strain_invariant * strain_invariant * std::sqrt(strain_invariant) +
      square_invariant * std::sqrt(std::sqrt(square_invariant));
  if (denominator == 0.0)
    return 0.0;
  const double length = coefficient * filter_width_m;
  const double result =
      length * length * scale * numerator / denominator;
  if (!(result >= 0.0) || !std::isfinite(result))
    throw runtime::Error("WALE kinematic viscosity is non-finite");
  return result == 0.0 ? 0.0 : result;
}

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

StructuredIndex map_cell(runtime::Int3 global, runtime::Box3 owned,
                         runtime::Int3 global_extent) {
  const runtime::Int3 local{owned.end.x - owned.begin.x,
                            owned.end.y - owned.begin.y,
                            owned.end.z - owned.begin.z};
  const auto axis = [](int coordinate, int begin, int end, int global_n,
                       int local_n) {
    if (coordinate >= begin && coordinate < end)
      return coordinate - begin;
    if (coordinate == begin - 1 || (begin == 0 && coordinate == global_n - 1))
      return -1;
    if (coordinate == end || (end == global_n && coordinate == 0))
      return local_n;
    throw runtime::Error("WALE active cell has no structured mapping");
  };
  return {axis(global.x, owned.begin.x, owned.end.x, global_extent.x, local.x),
          axis(global.y, owned.begin.y, owned.end.y, global_extent.y, local.y),
          axis(global.z, owned.begin.z, owned.end.z, global_extent.z, local.z)};
}

std::size_t coefficient_bytes(std::size_t count) {
  if (count > std::numeric_limits<std::size_t>::max() /
                  (2U * sizeof(double)))
    throw runtime::Error("WALE coefficient storage size overflows");
  return 2U * count * sizeof(double);
}

} // namespace

struct WaleAttemptCoefficients::Impl final {
  Impl(execution::ExecutionContext &execution, std::size_t count,
       std::size_t owned)
      : coefficients(execution, coefficient_bytes(count)), local_count(count),
        owned_count(owned) {}

  execution::Buffer coefficients;
  std::size_t local_count{};
  std::size_t owned_count{};
  WaleSummary summary;
};

struct WaleModel::Impl final {
  WaleControl control;
  execution::ExecutionContext *execution{};
  runtime::Int3 local_extent{};
  std::size_t owned_count{};
  std::vector<StructuredIndex> indices;
  std::vector<double> filter_width_m;
  std::uint64_t model_fingerprint{};
  mutable diagnostics::Stage3PerformanceCounters performance;
};

WaleAttemptCoefficients::WaleAttemptCoefficients(
    std::unique_ptr<Impl> implementation) noexcept
    : impl_(std::move(implementation)) {}
WaleAttemptCoefficients::~WaleAttemptCoefficients() noexcept = default;
WaleAttemptCoefficients::WaleAttemptCoefficients(
    WaleAttemptCoefficients &&) noexcept = default;

WaleCoefficientIdentity WaleAttemptCoefficients::identity() const noexcept {
  return impl_ ? impl_->summary.identity : WaleCoefficientIdentity{};
}
const WaleSummary &WaleAttemptCoefficients::summary() const noexcept {
  static const WaleSummary empty{};
  return impl_ ? impl_->summary : empty;
}
std::size_t WaleAttemptCoefficients::owned_active_count() const noexcept {
  return impl_ ? impl_->owned_count : 0U;
}
std::size_t WaleAttemptCoefficients::local_active_count() const noexcept {
  return impl_ ? impl_->local_count : 0U;
}
execution::VectorView<const double>
WaleAttemptCoefficients::nu_t_m2_per_s() const {
  if (!impl_)
    throw runtime::Error("WALE coefficients have been moved from");
  return static_cast<const execution::Buffer &>(impl_->coefficients)
      .view(0U, impl_->local_count);
}
execution::VectorView<const double>
WaleAttemptCoefficients::mu_sgs_pa_s() const {
  if (!impl_)
    throw runtime::Error("WALE coefficients have been moved from");
  return static_cast<const execution::Buffer &>(impl_->coefficients)
      .view(impl_->local_count, impl_->local_count);
}

WaleModel::WaleModel(std::unique_ptr<Impl> implementation) noexcept
    : impl_(std::move(implementation)) {}
WaleModel::~WaleModel() noexcept = default;
WaleModel::WaleModel(WaleModel &&) noexcept = default;

WaleModel WaleModel::create(
    WaleControl control, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry, std::size_t owned_active_count,
    const std::vector<mesh::GlobalCellId> &ordered_active,
    execution::ExecutionContext &execution) {
  if (!std::isfinite(control.coefficient) || control.coefficient < 1.0e-6 ||
      control.coefficient > 1.0 ||
      !std::isfinite(control.turbulent_prandtl) ||
      control.turbulent_prandtl < 0.1 || control.turbulent_prandtl > 10.0 ||
      !std::isfinite(control.turbulent_schmidt) ||
      control.turbulent_schmidt < 0.1 || control.turbulent_schmidt > 10.0)
    throw runtime::Error("WALE control is invalid");
  if (execution.backend_identity() == 0U ||
      execution.space() != execution::ExecutionSpace::host ||
      !execution.supports(execution::ExecutionCapability::buffer_allocation) ||
      !execution.supports(execution::ExecutionCapability::host_access))
    throw runtime::Error("WALE requires a host-access execution backend");
  geometry.require_compatible(topology);
  if (owned_active_count > ordered_active.size())
    throw runtime::Error("WALE owned active count is invalid");
  std::set<mesh::GlobalCellId> unique;
  auto impl = std::make_unique<Impl>();
  impl->control = control;
  impl->execution = &execution;
  const auto box = topology.owned_global_box();
  impl->local_extent = {box.end.x - box.begin.x, box.end.y - box.begin.y,
                        box.end.z - box.begin.z};
  impl->owned_count = owned_active_count;
  impl->indices.reserve(ordered_active.size());
  impl->filter_width_m.reserve(ordered_active.size());
  std::uint64_t fingerprint = kFnvOffset;
  fingerprint = mix(fingerprint, bits(control.coefficient));
  fingerprint = mix(fingerprint, bits(control.turbulent_prandtl));
  fingerprint = mix(fingerprint, bits(control.turbulent_schmidt));
  for (std::size_t row = 0U; row < ordered_active.size(); ++row) {
    if (!unique.insert(ordered_active[row]).second)
      throw runtime::Error("WALE active global cell IDs are not unique");
    const auto local = topology.find_local_cell(ordered_active[row]);
    if (!local.has_value())
      throw runtime::Error("WALE active global cell is not local");
    const bool owned = topology.cell_ownership(*local) ==
                       mesh::EntityOwnership::owned;
    if (owned != (row < owned_active_count))
      throw runtime::Error("WALE active cells are not owned-first");
    const double volume = geometry.cell_volume_m3(*local);
    if (!(volume > 0.0) || !std::isfinite(volume))
      throw runtime::Error("WALE active cell volume is invalid");
    impl->indices.push_back(map_cell(topology.global_cell(*local), box,
                                     topology.global_extent()));
    impl->filter_width_m.push_back(std::cbrt(volume));
    fingerprint = mix(fingerprint, ordered_active[row]);
    fingerprint = mix(fingerprint, bits(volume));
  }
  impl->model_fingerprint = fingerprint == 0U ? 1U : fingerprint;
  return WaleModel(std::move(impl));
}

WaleAttemptCoefficients WaleModel::evaluate(
    const WaleAttemptInput &input) const {
  if (!impl_)
    throw runtime::Error("WALE model has been moved from");
  if (!(input.attempted_dt_s > 0.0) ||
      !std::isfinite(input.attempted_dt_s) ||
      (input.order != WaleTimeOrder::backward_euler &&
       input.order != WaleTimeOrder::bdf2) ||
      input.lagged_velocity_gradient.components() != 9U ||
      input.rho_attempt.components() != 1U ||
      input.lagged_velocity_gradient.interior_extent().x !=
          impl_->local_extent.x ||
      input.lagged_velocity_gradient.interior_extent().y !=
          impl_->local_extent.y ||
      input.lagged_velocity_gradient.interior_extent().z !=
          impl_->local_extent.z ||
      input.rho_attempt.interior_extent().x != impl_->local_extent.x ||
      input.rho_attempt.interior_extent().y != impl_->local_extent.y ||
      input.rho_attempt.interior_extent().z != impl_->local_extent.z)
    throw runtime::Error("WALE attempt input is invalid");
  const auto covered = [](StructuredIndex index, runtime::Int3 extent,
                          int ghost_width) {
    return index.i >= -ghost_width && index.i < extent.x + ghost_width &&
           index.j >= -ghost_width && index.j < extent.y + ghost_width &&
           index.k >= -ghost_width && index.k < extent.z + ghost_width;
  };
  for (const auto index : impl_->indices)
    if (!covered(index, input.lagged_velocity_gradient.interior_extent(),
                 input.lagged_velocity_gradient.ghost_width()) ||
        !covered(index, input.rho_attempt.interior_extent(),
                 input.rho_attempt.ghost_width()))
      throw runtime::Error("WALE active cell exceeds the input ghost extent");

  const auto add_work = [](std::uint64_t &counter, std::uint64_t value) {
    if (value > std::numeric_limits<std::uint64_t>::max() - counter)
      throw runtime::Error("WALE performance counter would overflow");
    counter += value;
  };
  add_work(impl_->performance.step_wale_evaluations, 1U);

  std::unique_ptr<WaleAttemptCoefficients::Impl> output;
  runtime::with_kernel_cell_view(
      input.lagged_velocity_gradient, [&](auto gradient) {
        runtime::with_kernel_cell_view(input.rho_attempt, [&](auto density) {
          output = std::make_unique<WaleAttemptCoefficients::Impl>(
              *impl_->execution, impl_->indices.size(), impl_->owned_count);
          auto nu_view = output->coefficients.view(0U, impl_->indices.size());
          auto mu_view = output->coefficients.view(impl_->indices.size(),
                                                   impl_->indices.size());
          double *const nu = nu_view.data();
          double *const mu = mu_view.data();
          double minimum = std::numeric_limits<double>::infinity();
          double maximum = 0.0;
          double l2 = 0.0;
          std::uint64_t zeros{};
          for (std::size_t row = 0U; row < impl_->indices.size(); ++row) {
            const auto index = impl_->indices[row];
            std::array<double, 9> tensor{};
            for (std::size_t component = 0U; component < tensor.size();
                 ++component)
              tensor[component] = gradient(index.i, index.j, index.k,
                                            static_cast<int>(component));
            if (row < impl_->owned_count)
              add_work(impl_->performance.step_wale_gradient_cells, 1U);
            const double rho = density(index.i, index.j, index.k, 0);
            if (!(rho > 0.0) || !std::isfinite(rho))
              throw runtime::Error("WALE attempt density is invalid");
            nu[row] = wale_kinematic_viscosity(
                tensor, impl_->control.coefficient,
                impl_->filter_width_m[row]);
            mu[row] = rho * nu[row];
            if (!std::isfinite(mu[row]))
              throw runtime::Error("WALE dynamic viscosity is non-finite");
            if (row < impl_->owned_count) {
              minimum = std::min(minimum, nu[row]);
              maximum = std::max(maximum, nu[row]);
              l2 = std::hypot(l2, nu[row]);
              if (nu[row] == 0.0 && !std::signbit(nu[row]))
                ++zeros;
            }
          }
          std::uint64_t identity = impl_->model_fingerprint;
          identity = mix(identity, input.step);
          identity = mix(identity, bits(input.attempted_dt_s));
          identity = mix(identity, static_cast<std::uint64_t>(input.order));
          identity = mix(identity, input.committed_state_fingerprint);
          identity = mix(identity, input.history_state_fingerprint);
          identity = mix(identity, input.lagged_gradient_fingerprint);
          identity = mix(identity, input.density_fingerprint);
          if (identity == 0U)
            identity = 1U;
          output->summary = {
              {identity},
              impl_->owned_count == 0U ? 0.0 : minimum,
              impl_->owned_count == 0U ? 0.0 : maximum,
              l2,
              zeros,
              static_cast<std::uint64_t>(impl_->owned_count)};
        });
      });
  return WaleAttemptCoefficients(std::move(output));
}

WaleControl WaleModel::control() const noexcept {
  return impl_ ? impl_->control : WaleControl{};
}

diagnostics::Stage3PerformanceCounters
WaleModel::performance_counters() const noexcept {
  return impl_ ? impl_->performance
               : diagnostics::Stage3PerformanceCounters{};
}

#ifdef HUNDUN_LES_ENABLE_TEST_ACCESS
namespace test {
double wale_kinematic_viscosity_for_test(
    const std::array<double, 9> &gradient, double coefficient,
    double filter_width_m) {
  if (!std::isfinite(coefficient) || coefficient < 1.0e-6 ||
      coefficient > 1.0 || !(filter_width_m > 0.0) ||
      !std::isfinite(filter_width_m))
    throw runtime::Error("WALE test input is invalid");
  return wale_kinematic_viscosity(gradient, coefficient, filter_width_m);
}
} // namespace test
#endif

} // namespace hundun::les
