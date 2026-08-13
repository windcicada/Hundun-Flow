// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include "hundun/v04_case.hpp"
#include "hundun/v04_mesh.hpp"

#include "solver_cartesian_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTransportKernel = 921U;
constexpr std::uint32_t kTransportNumerical = 922U;
constexpr std::int32_t kXBlock = 16;

template <std::size_t Axis>
inline Int3 offset_axis(Int3 value, std::int32_t amount) noexcept {
  if constexpr (Axis == 0U) {
    value.x += amount;
  } else if constexpr (Axis == 1U) {
    value.y += amount;
  } else {
    value.z += amount;
  }
  return value;
}

template <std::size_t Axis>
inline std::int32_t axis_index(Int3 value) noexcept {
  if constexpr (Axis == 0U) {
    return value.x;
  } else if constexpr (Axis == 1U) {
    return value.y;
  } else {
    return value.z;
  }
}

inline double minmod(double left, double right) noexcept {
  if (left * right <= 0.0) {
    return 0.0;
  }
  return std::copysign(std::min(std::abs(left), std::abs(right)), left);
}

template <bool Uniform, std::size_t Axis>
inline double limited_slope_values(const CartesianKernelPlan& plan,
                                   std::int32_t normal, double q_minus,
                                   double q_centre, double q_plus) noexcept {
  const double centre_coordinate =
      detail::metric_centre<Uniform>(plan, Axis, normal);
  const double left_distance =
      centre_coordinate -
      detail::metric_centre<Uniform>(plan, Axis, normal - 1);
  const double right_distance =
      detail::metric_centre<Uniform>(plan, Axis, normal + 1) -
      centre_coordinate;
  const double delta_left = (q_centre - q_minus) / left_distance;
  const double delta_right = (q_plus - q_centre) / right_distance;
  const double centred =
      (right_distance * delta_left + left_distance * delta_right) /
      (left_distance + right_distance);
  const double monotone = minmod(2.0 * delta_left, 2.0 * delta_right);
  return plan.limiter() * minmod(centred, monotone);
}

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
inline double reconstructed_face(const CartesianKernelPlan& plan,
                                 ConstFieldView field,
                                 std::uint8_t component_id, Int3 face,
                                 double mass_rate) noexcept {
  const std::int32_t normal = axis_index<Axis>(face);
  const Int3 left = offset_axis<Axis>(face, -1);
  const Int3 right = face;
  const double q_left = field.unchecked(left, component_id);
  const double q_right = field.unchecked(right, component_id);
  if constexpr (Scheme == ConvectionScheme::central2) {
    return detail::metric_interpolate_face<Uniform>(
        plan, Axis, normal, q_left, q_right);
  }

  const double face_coordinate =
      detail::metric_face<Uniform>(plan, Axis, normal);
  const double q_left_left =
      field.unchecked(offset_axis<Axis>(left, -1), component_id);
  const double q_right_right =
      field.unchecked(offset_axis<Axis>(right, 1), component_id);
  const double left_reconstructed =
      q_left + limited_slope_values<Uniform, Axis>(
                   plan, normal - 1, q_left_left, q_left, q_right) *
                   (face_coordinate -
                    detail::metric_centre<Uniform>(plan, Axis, normal - 1));
  const double right_reconstructed =
      q_right + limited_slope_values<Uniform, Axis>(
                    plan, normal, q_left, q_right, q_right_right) *
          (face_coordinate -
           detail::metric_centre<Uniform>(plan, Axis, normal));
  if constexpr (Scheme == ConvectionScheme::tvd2) {
    return mass_rate >= 0.0 ? left_reconstructed : right_reconstructed;
  }
  return 0.5 * (left_reconstructed + right_reconstructed);
}

template <std::size_t Axis>
inline double face_rate(ConstFaceFluxView flux, Int3 face) noexcept {
  if constexpr (Axis == 0U) {
    return flux.x.unchecked(face);
  } else if constexpr (Axis == 1U) {
    return flux.y.unchecked(face);
  } else {
    return flux.z.unchecked(face);
  }
}

template <bool Uniform, std::size_t Axis>
inline double diffusion_transmissibility(const CartesianKernelPlan& plan,
                                         ConstFieldView diffusivity,
                                         Int3 face) noexcept {
  const std::int32_t normal = axis_index<Axis>(face);
  const Int3 left = offset_axis<Axis>(face, -1);
  const double gamma_left = diffusivity.unchecked(left, 0U);
  const double gamma_right = diffusivity.unchecked(face, 0U);
  const double face_coordinate =
      detail::metric_face<Uniform>(plan, Axis, normal);
  const double distance_left =
      face_coordinate -
      detail::metric_centre<Uniform>(plan, Axis, normal - 1);
  const double distance_right =
      detail::metric_centre<Uniform>(plan, Axis, normal) - face_coordinate;
  if (!std::isfinite(gamma_left) || !std::isfinite(gamma_right) ||
      gamma_left <= 0.0 || gamma_right <= 0.0 || distance_left <= 0.0 ||
      distance_right <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return detail::metric_face_area<Uniform>(plan, Axis, face) /
         (distance_left / gamma_left + distance_right / gamma_right);
}

bool valid_transport_invocation(const CartesianKernelPlan& plan,
                                const KernelInvocation& invocation,
                                std::uint8_t required_ghost_width,
                                bool face_flux_required) noexcept {
  return plan.fingerprint() != 0U && invocation.reads.data != nullptr &&
         invocation.reads.size == 1U && invocation.writes.data != nullptr &&
         invocation.writes.size == 1U &&
         detail::valid_kernel_box(invocation.box, plan.cells()) &&
         detail::valid_cell_view(invocation.reads.data[0U], plan.cells(),
                                 invocation.read_component_begin,
                                 invocation.component_count,
                                 required_ghost_width) &&
         detail::valid_cell_view(invocation.writes.data[0U], plan.cells(),
                                 invocation.write_component_begin,
                                 invocation.component_count) &&
         (!face_flux_required ||
          invocation.required_face_flux_revision != 0U);
}

template <bool Uniform, ConvectionScheme Scheme>
Status convection_kernel(const CartesianKernelPlan& plan,
                         ConstFaceFluxView flux,
                         const KernelInvocation& invocation,
                         ConstFieldView transported,
                         FieldView output) noexcept {
  const std::uint64_t components = invocation.component_count;
  const std::uint64_t scalar_reads =
      Scheme == ConvectionScheme::central2 ? 12U : 24U;
  KernelCounters prepared{};
  if (invocation.counters != nullptr) {
    prepared = *invocation.counters;
    const Status counted = detail::add_kernel_counters(
        &prepared, invocation.box, 6U,
        6U + scalar_reads * components, components);
    if (!counted) {
      return counted;
    }
  }
  const Int3 end{invocation.box.begin.x + invocation.box.cells.x,
                 invocation.box.begin.y + invocation.box.cells.y,
                 invocation.box.begin.z + invocation.box.cells.z};
  for (std::int32_t z = invocation.box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = invocation.box.begin.y; y < end.y; ++y) {
      alignas(64) std::array<double, kXBlock> mxm;
      alignas(64) std::array<double, kXBlock> mxp;
      alignas(64) std::array<double, kXBlock> mym;
      alignas(64) std::array<double, kXBlock> myp;
      alignas(64) std::array<double, kXBlock> mzm;
      alignas(64) std::array<double, kXBlock> mzp;
      alignas(64) std::array<double, kXBlock> inverse_volume;
      for (std::int32_t x_begin = invocation.box.begin.x; x_begin < end.x;) {
        const std::int32_t lanes =
            std::min(kXBlock, static_cast<std::int32_t>(end.x - x_begin));
        for (std::int32_t lane = 0; lane < lanes; ++lane) {
          const std::int32_t x = x_begin + lane;
          const Int3 cell{x, y, z};
          mxm[static_cast<std::size_t>(lane)] = face_rate<0U>(flux, cell);
          mxp[static_cast<std::size_t>(lane)] =
              face_rate<0U>(flux, {x + 1, y, z});
          mym[static_cast<std::size_t>(lane)] = face_rate<1U>(flux, cell);
          myp[static_cast<std::size_t>(lane)] =
              face_rate<1U>(flux, {x, y + 1, z});
          mzm[static_cast<std::size_t>(lane)] = face_rate<2U>(flux, cell);
          mzp[static_cast<std::size_t>(lane)] =
              face_rate<2U>(flux, {x, y, z + 1});
          inverse_volume[static_cast<std::size_t>(lane)] =
              detail::metric_inverse_volume<Uniform>(plan, cell);
        }
        for (std::uint8_t component_id = 0U;
             component_id < invocation.component_count; ++component_id) {
          const std::uint8_t read_component = static_cast<std::uint8_t>(
              invocation.read_component_begin + component_id);
          const std::uint8_t write_component = static_cast<std::uint8_t>(
              invocation.write_component_begin + component_id);
          for (std::int32_t lane = 0; lane < lanes; ++lane) {
            const std::size_t selected = static_cast<std::size_t>(lane);
            const std::int32_t x = x_begin + lane;
            const Int3 cell{x, y, z};
            const Int3 xp{x + 1, y, z};
            const Int3 yp{x, y + 1, z};
            const Int3 zp{x, y, z + 1};
            const double raw =
                mxp[selected] * reconstructed_face<Uniform, 0U, Scheme>(
                                    plan, transported, read_component, xp,
                                    mxp[selected]) -
                mxm[selected] * reconstructed_face<Uniform, 0U, Scheme>(
                                    plan, transported, read_component, cell,
                                    mxm[selected]) +
                myp[selected] * reconstructed_face<Uniform, 1U, Scheme>(
                                    plan, transported, read_component, yp,
                                    myp[selected]) -
                mym[selected] * reconstructed_face<Uniform, 1U, Scheme>(
                                    plan, transported, read_component, cell,
                                    mym[selected]) +
                mzp[selected] * reconstructed_face<Uniform, 2U, Scheme>(
                                    plan, transported, read_component, zp,
                                    mzp[selected]) -
                mzm[selected] * reconstructed_face<Uniform, 2U, Scheme>(
                                    plan, transported, read_component, cell,
                                    mzm[selected]);
            const double value = raw * inverse_volume[selected];
            if (!std::isfinite(value)) {
              return {StatusCode::numerical_failure, kTransportNumerical};
            }
            output.unchecked(cell, write_component) = value;
          }
        }
        x_begin += lanes;
      }
    }
  }
  if (invocation.counters != nullptr) {
    *invocation.counters = prepared;
  }
  return {};
}

template <bool Uniform>
Status diffusion_kernel(const CartesianKernelPlan& plan,
                        ConstFieldView diffusivity,
                        const KernelInvocation& invocation,
                        ConstFieldView transported,
                        FieldView output) noexcept {
  const std::uint64_t components = invocation.component_count;
  KernelCounters prepared{};
  if (invocation.counters != nullptr) {
    prepared = *invocation.counters;
    const Status counted = detail::add_kernel_counters(
        &prepared, invocation.box, 6U, 12U + 12U * components, components);
    if (!counted) {
      return counted;
    }
  }
  const Int3 end{invocation.box.begin.x + invocation.box.cells.x,
                 invocation.box.begin.y + invocation.box.cells.y,
                 invocation.box.begin.z + invocation.box.cells.z};
  for (std::int32_t z = invocation.box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = invocation.box.begin.y; y < end.y; ++y) {
      alignas(64) std::array<double, kXBlock> txm;
      alignas(64) std::array<double, kXBlock> txp;
      alignas(64) std::array<double, kXBlock> tym;
      alignas(64) std::array<double, kXBlock> typ;
      alignas(64) std::array<double, kXBlock> tzm;
      alignas(64) std::array<double, kXBlock> tzp;
      alignas(64) std::array<double, kXBlock> inverse_volume;
      for (std::int32_t x_begin = invocation.box.begin.x; x_begin < end.x;) {
        const std::int32_t lanes =
            std::min(kXBlock, static_cast<std::int32_t>(end.x - x_begin));
        for (std::int32_t lane = 0; lane < lanes; ++lane) {
          const std::size_t selected = static_cast<std::size_t>(lane);
          const std::int32_t x = x_begin + lane;
          const Int3 cell{x, y, z};
          txm[selected] = diffusion_transmissibility<Uniform, 0U>(
              plan, diffusivity, cell);
          txp[selected] = diffusion_transmissibility<Uniform, 0U>(
              plan, diffusivity, {x + 1, y, z});
          tym[selected] = diffusion_transmissibility<Uniform, 1U>(
              plan, diffusivity, cell);
          typ[selected] = diffusion_transmissibility<Uniform, 1U>(
              plan, diffusivity, {x, y + 1, z});
          tzm[selected] = diffusion_transmissibility<Uniform, 2U>(
              plan, diffusivity, cell);
          tzp[selected] = diffusion_transmissibility<Uniform, 2U>(
              plan, diffusivity, {x, y, z + 1});
          inverse_volume[selected] =
              detail::metric_inverse_volume<Uniform>(plan, cell);
          if (!std::isfinite(txm[selected]) ||
              !std::isfinite(txp[selected]) ||
              !std::isfinite(tym[selected]) ||
              !std::isfinite(typ[selected]) ||
              !std::isfinite(tzm[selected]) ||
              !std::isfinite(tzp[selected])) {
            return {StatusCode::numerical_failure, kTransportNumerical};
          }
        }
        for (std::uint8_t component_id = 0U;
             component_id < invocation.component_count; ++component_id) {
          const std::uint8_t read_component = static_cast<std::uint8_t>(
              invocation.read_component_begin + component_id);
          const std::uint8_t write_component = static_cast<std::uint8_t>(
              invocation.write_component_begin + component_id);
          for (std::int32_t lane = 0; lane < lanes; ++lane) {
            const std::size_t selected = static_cast<std::size_t>(lane);
            const std::int32_t x = x_begin + lane;
            const Int3 cell{x, y, z};
            const double q = transported.unchecked(cell, read_component);
            const double raw =
                txp[selected] *
                    (transported.unchecked({x + 1, y, z}, read_component) -
                     q) -
                txm[selected] *
                    (q - transported.unchecked({x - 1, y, z},
                                               read_component)) +
                typ[selected] *
                    (transported.unchecked({x, y + 1, z}, read_component) -
                     q) -
                tym[selected] *
                    (q - transported.unchecked({x, y - 1, z},
                                               read_component)) +
                tzp[selected] *
                    (transported.unchecked({x, y, z + 1}, read_component) -
                     q) -
                tzm[selected] *
                    (q - transported.unchecked({x, y, z - 1},
                                               read_component));
            const double value = raw * inverse_volume[selected];
            if (!std::isfinite(value)) {
              return {StatusCode::numerical_failure, kTransportNumerical};
            }
            output.unchecked(cell, write_component) = value;
          }
        }
        x_begin += lanes;
      }
    }
  }
  if (invocation.counters != nullptr) {
    *invocation.counters = prepared;
  }
  return {};
}

template <ConvectionScheme Scheme>
Status dispatch_convection_geometry(const CartesianKernelPlan& plan,
                                    ConstFaceFluxView flux,
                                    const KernelInvocation& invocation,
                                    ConstFieldView transported,
                                    FieldView output) noexcept {
  return plan.geometry_kind() == GeometryKind::uniform
             ? convection_kernel<true, Scheme>(plan, flux, invocation,
                                                transported, output)
             : convection_kernel<false, Scheme>(plan, flux, invocation,
                                                 transported, output);
}

}  // namespace

Status cartesian_provisional_convection(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView flux, const KernelInvocation& invocation) noexcept {
  const std::uint8_t required_ghost_width =
      scheme == ConvectionScheme::central2 ? 1U : 2U;
  if (!valid_transport_invocation(plan, invocation, required_ghost_width,
                                  true) ||
      static_cast<std::uint8_t>(scheme) >
          static_cast<std::uint8_t>(ConvectionScheme::tvd2) ||
      !detail::valid_flux_view(flux, plan.cells(),
                               invocation.required_face_flux_revision) ||
      detail::component_ranges_overlap(
          invocation.reads.data[0U], invocation.read_component_begin,
          invocation.component_count, invocation.writes.data[0U],
          invocation.write_component_begin, invocation.component_count)) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  const ConstFieldView transported = invocation.reads.data[0U];
  const FieldView output = invocation.writes.data[0U];
  switch (scheme) {
    case ConvectionScheme::central2:
      return dispatch_convection_geometry<ConvectionScheme::central2>(
          plan, flux, invocation, transported, output);
    case ConvectionScheme::limited_central2:
      return dispatch_convection_geometry<
          ConvectionScheme::limited_central2>(plan, flux, invocation,
                                               transported, output);
    case ConvectionScheme::tvd2:
      return dispatch_convection_geometry<ConvectionScheme::tvd2>(
          plan, flux, invocation, transported, output);
  }
  return {StatusCode::invalid_plan, kTransportKernel};
}

Status cartesian_convection(const CartesianKernelPlan& plan,
                            ConvectionScheme scheme, ConstFaceFluxView flux,
                            const KernelInvocation& invocation) noexcept {
  if (!flux.certificate.matches(flux)) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  return cartesian_provisional_convection(plan, scheme, flux, invocation);
}

Status cartesian_diffusion(const CartesianKernelPlan& plan,
                           ConstFieldView diffusivity,
                           const KernelInvocation& invocation) noexcept {
  if (!valid_transport_invocation(plan, invocation, 1U, false) ||
      invocation.required_face_flux_revision != 0U ||
      !detail::valid_cell_view(diffusivity, plan.cells(), 0U, 1U, 1U) ||
      detail::component_ranges_overlap(
          invocation.reads.data[0U], invocation.read_component_begin,
          invocation.component_count, invocation.writes.data[0U],
          invocation.write_component_begin, invocation.component_count) ||
      detail::component_ranges_overlap(
          diffusivity, 0U, 1U, invocation.writes.data[0U],
          invocation.write_component_begin, invocation.component_count)) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  const ConstFieldView transported = invocation.reads.data[0U];
  const FieldView output = invocation.writes.data[0U];
  return plan.geometry_kind() == GeometryKind::uniform
             ? diffusion_kernel<true>(plan, diffusivity, invocation,
                                      transported, output)
             : diffusion_kernel<false>(plan, diffusivity, invocation,
                                       transported, output);
}

}  // namespace hundun::v04
