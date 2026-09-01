// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include "hundun/v04_case.hpp"
#include "hundun/v04_mesh.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTransportKernel = 921U;
constexpr std::uint32_t kTransportNumerical = 922U;
constexpr std::uint32_t kTransportNondifferentiable = 923U;
constexpr std::int32_t kXBlock = 16;
constexpr std::uint64_t kFrozenConvectionSchema =
    UINT64_C(0x7630346672636f6e);
constexpr std::uint64_t kFrozenConvectionDerivativeSchema =
    UINT64_C(0x7630346672646572);
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t nonzero_hash(std::uint64_t value) noexcept {
  return value == 0U ? 1U : value;
}

bool same_bits(double left, double right) noexcept {
  std::uint64_t left_bits = 0U;
  std::uint64_t right_bits = 0U;
  std::memcpy(&left_bits, &left, sizeof(left_bits));
  std::memcpy(&right_bits, &right, sizeof(right_bits));
  return left_bits == right_bits;
}

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
  if (left == 0.0 || right == 0.0 ||
      std::signbit(left) != std::signbit(right)) {
    return 0.0;
  }
  return std::copysign(std::min(std::abs(left), std::abs(right)), left);
}

enum class LimitedSlopeBranch : std::uint8_t {
  zero = 1U,
  centred = 2U,
  left_delta = 3U,
  right_delta = 4U,
  nondifferentiable = 5U
};

struct LimitedSlopeEvaluation {
  double value{};
  LimitedSlopeBranch branch{LimitedSlopeBranch::nondifferentiable};
};

bool branch_tie(double left, double right) noexcept {
  const double scale = std::max(std::abs(left), std::abs(right));
  return left == right ||
         (scale > 0.0 &&
          std::abs(left - right) <=
              64.0 * std::numeric_limits<double>::epsilon() * scale);
}

template <bool Uniform, std::size_t Axis>
inline LimitedSlopeEvaluation evaluate_limited_slope_values(
    const CartesianKernelPlan& plan, std::int32_t normal, double q_minus,
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
  LimitedSlopeEvaluation result{
      plan.limiter() * minmod(centred, monotone),
      LimitedSlopeBranch::nondifferentiable};
  if (!std::isfinite(result.value) || !std::isfinite(delta_left) ||
      !std::isfinite(delta_right) || !std::isfinite(centred) ||
      !std::isfinite(monotone)) {
    return result;
  }
  if (plan.limiter() == 0.0) {
    result.branch = LimitedSlopeBranch::zero;
    return result;
  }
  if (delta_left == 0.0 || delta_right == 0.0) {
    return result;
  }
  if (std::signbit(delta_left) != std::signbit(delta_right)) {
    result.branch = LimitedSlopeBranch::zero;
    return result;
  }
  const double centred_magnitude = std::abs(centred);
  const double monotone_magnitude = std::abs(monotone);
  if (branch_tie(centred_magnitude, monotone_magnitude)) {
    return result;
  }
  if (centred_magnitude < monotone_magnitude) {
    result.branch = LimitedSlopeBranch::centred;
    return result;
  }
  const double left_magnitude = std::abs(delta_left);
  const double right_magnitude = std::abs(delta_right);
  if (branch_tie(left_magnitude, right_magnitude)) {
    return result;
  }
  result.branch = left_magnitude < right_magnitude
                      ? LimitedSlopeBranch::left_delta
                      : LimitedSlopeBranch::right_delta;
  return result;
}

template <bool Uniform, std::size_t Axis>
inline double limited_slope_values(const CartesianKernelPlan& plan,
                                   std::int32_t normal, double q_minus,
                                   double q_centre, double q_plus) noexcept {
  return evaluate_limited_slope_values<Uniform, Axis>(
             plan, normal, q_minus, q_centre, q_plus)
      .value;
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

struct FaceBranchSelection {
  LimitedSlopeEvaluation left{};
  LimitedSlopeEvaluation right{};
  std::uint64_t code{};
  bool differentiable{};
};

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
inline FaceBranchSelection select_face_branches(
    const CartesianKernelPlan& plan, ConstFieldView target,
    std::uint8_t component, Int3 face, double mass_rate) noexcept {
  if constexpr (Scheme == ConvectionScheme::central2) {
    return {{}, {}, UINT64_C(0x101), true};
  }
  const std::int32_t normal = axis_index<Axis>(face);
  const Int3 left = offset_axis<Axis>(face, -1);
  const Int3 right = face;
  const double q_left = target.unchecked(left, component);
  const double q_right = target.unchecked(right, component);
  const LimitedSlopeEvaluation left_slope =
      evaluate_limited_slope_values<Uniform, Axis>(
          plan, normal - 1,
          target.unchecked(offset_axis<Axis>(left, -1), component), q_left,
          q_right);
  const LimitedSlopeEvaluation right_slope =
      evaluate_limited_slope_values<Uniform, Axis>(
          plan, normal, q_left, q_right,
          target.unchecked(offset_axis<Axis>(right, 1), component));
  if constexpr (Scheme == ConvectionScheme::tvd2) {
    const bool left_donor = mass_rate >= 0.0;
    const LimitedSlopeBranch selected =
        left_donor ? left_slope.branch : right_slope.branch;
    const std::uint64_t code =
        UINT64_C(0x300) |
        (left_donor ? UINT64_C(0x80) : UINT64_C(0x40)) |
        static_cast<std::uint8_t>(selected);
    return {left_slope, right_slope, code,
            selected != LimitedSlopeBranch::nondifferentiable};
  }
  const std::uint64_t code =
      UINT64_C(0x200) |
      static_cast<std::uint8_t>(left_slope.branch) |
      (static_cast<std::uint64_t>(
           static_cast<std::uint8_t>(right_slope.branch))
       << 4U);
  return {left_slope, right_slope, code,
          left_slope.branch != LimitedSlopeBranch::nondifferentiable &&
              right_slope.branch !=
                  LimitedSlopeBranch::nondifferentiable};
}

template <bool Uniform, std::size_t Axis>
inline double directional_limited_slope_values(
    const CartesianKernelPlan& plan, std::int32_t normal, double q_minus,
    double q_centre, double q_plus, LimitedSlopeBranch branch) noexcept {
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
  switch (branch) {
    case LimitedSlopeBranch::zero:
      return 0.0;
    case LimitedSlopeBranch::centred:
      return plan.limiter() *
             (right_distance * delta_left + left_distance * delta_right) /
             (left_distance + right_distance);
    case LimitedSlopeBranch::left_delta:
      return plan.limiter() * 2.0 * delta_left;
    case LimitedSlopeBranch::right_delta:
      return plan.limiter() * 2.0 * delta_right;
    case LimitedSlopeBranch::nondifferentiable:
      break;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

struct DirectionalFaceEvaluation {
  double value{};
  std::uint64_t branch_code{};
  bool differentiable{};
  bool generalized{};
};

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
inline DirectionalFaceEvaluation reconstructed_face_direction(
    const CartesianKernelPlan& plan, ConstFieldView target,
    std::uint8_t target_component, ConstFieldView variation,
    std::uint8_t variation_component, Int3 face,
    double mass_rate, FrozenConvectionLinearizationPolicy policy) noexcept {
  const FaceBranchSelection selection = select_face_branches<
      Uniform, Axis, Scheme>(plan, target, target_component, face, mass_rate);
  const std::int32_t normal = axis_index<Axis>(face);
  const Int3 left = offset_axis<Axis>(face, -1);
  const Int3 right = face;
  const double v_left = variation.unchecked(left, variation_component);
  const double v_right = variation.unchecked(right, variation_component);
  if constexpr (Scheme == ConvectionScheme::central2) {
    return {detail::metric_interpolate_face<Uniform>(
                plan, Axis, normal, v_left, v_right),
            selection.code, true, false};
  }
  if (!selection.differentiable) {
    if (policy != FrozenConvectionLinearizationPolicy::
                      semismooth_generalized_zero_slope) {
      return {std::numeric_limits<double>::quiet_NaN(), selection.code, false,
              false};
    }
    if constexpr (Scheme == ConvectionScheme::tvd2) {
      return {mass_rate >= 0.0 ? v_left : v_right, selection.code, true,
              true};
    }
    return {0.5 * (v_left + v_right), selection.code, true, true};
  }
  const double face_coordinate =
      detail::metric_face<Uniform>(plan, Axis, normal);
  const double left_slope = directional_limited_slope_values<Uniform, Axis>(
      plan, normal - 1,
      variation.unchecked(offset_axis<Axis>(left, -1), variation_component),
      v_left, v_right, selection.left.branch);
  const double right_slope = directional_limited_slope_values<Uniform, Axis>(
      plan, normal, v_left, v_right,
      variation.unchecked(offset_axis<Axis>(right, 1), variation_component),
      selection.right.branch);
  const double left_reconstructed =
      v_left + left_slope *
                   (face_coordinate - detail::metric_centre<Uniform>(
                                          plan, Axis, normal - 1));
  const double right_reconstructed =
      v_right + right_slope *
                    (face_coordinate - detail::metric_centre<Uniform>(
                                           plan, Axis, normal));
  if constexpr (Scheme == ConvectionScheme::tvd2) {
    return {mass_rate >= 0.0 ? left_reconstructed : right_reconstructed,
            selection.code, true, false};
  }
  return {0.5 * (left_reconstructed + right_reconstructed), selection.code,
          true, false};
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

Int3 expected_face_extents(Int3 cells, CartesianAxis axis) noexcept {
  if (axis == CartesianAxis::x) {
    ++cells.x;
  } else if (axis == CartesianAxis::y) {
    ++cells.y;
  } else {
    ++cells.z;
  }
  return cells;
}

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template <class T>
bool valid_frozen_output(BasicFaceFieldView<T> view, CartesianAxis axis,
                         Int3 cells) noexcept {
  detail::FieldStorageInterval interval;
  return detail::face_storage_interval(view, interval) && view.axis == axis &&
         same_shape(view.extents, expected_face_extents(cells, axis)) &&
         view.storage_identity != 0U && view.revision_domain != 0U;
}

template <class T>
std::uint64_t mix_face_view(std::uint64_t hash,
                            BasicFaceFieldView<T> view) noexcept {
  hash = hash_mix(hash, reinterpret_cast<std::uintptr_t>(view.base));
  hash = hash_mix(hash, static_cast<std::uint64_t>(view.extents.x));
  hash = hash_mix(hash, static_cast<std::uint64_t>(view.extents.y));
  hash = hash_mix(hash, static_cast<std::uint64_t>(view.extents.z));
  hash = hash_mix(hash, view.stride_y);
  hash = hash_mix(hash, view.stride_z);
  hash = hash_mix(hash, static_cast<std::uint8_t>(view.axis));
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

std::uint64_t mix_field_view(std::uint64_t hash,
                             ConstFieldView view) noexcept {
  hash = hash_mix(hash, reinterpret_cast<std::uintptr_t>(view.base));
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.revision);
  hash = hash_mix(hash, view.storage_identity);
  hash = hash_mix(hash, view.revision_domain);
  hash = hash_mix(hash, view.stride_y);
  hash = hash_mix(hash, view.stride_z);
  hash = hash_mix(hash, view.component_stride);
  return hash_mix(hash, view.replica);
}

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
bool preflight_frozen_axis(const CartesianKernelPlan& plan,
                           ConstFaceFluxView flux,
                           ConstFieldView transported,
                           std::uint8_t component,
                           FaceFieldView output) noexcept {
  for (std::int32_t z = 0; z < output.extents.z; ++z) {
    for (std::int32_t y = 0; y < output.extents.y; ++y) {
      for (std::int32_t x = 0; x < output.extents.x; ++x) {
        const Int3 face{x, y, z};
        const double rate = face_rate<Axis>(flux, face);
        const double value = reconstructed_face<Uniform, Axis, Scheme>(
            plan, transported, component, face, rate);
        if (!std::isfinite(rate) || !std::isfinite(value)) return false;
      }
    }
  }
  return true;
}

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
void commit_frozen_axis(const CartesianKernelPlan& plan,
                        ConstFaceFluxView flux,
                        ConstFieldView transported,
                        std::uint8_t component,
                        FaceFieldView output) noexcept {
  for (std::int32_t z = 0; z < output.extents.z; ++z) {
    for (std::int32_t y = 0; y < output.extents.y; ++y) {
      for (std::int32_t x = 0; x < output.extents.x; ++x) {
        const Int3 face{x, y, z};
        output.unchecked(face) = reconstructed_face<Uniform, Axis, Scheme>(
            plan, transported, component, face,
            face_rate<Axis>(flux, face));
      }
    }
  }
}

template <bool Uniform, ConvectionScheme Scheme>
bool preflight_frozen_faces(const CartesianKernelPlan& plan,
                            ConstFaceFluxView flux,
                            ConstFieldView transported,
                            std::uint8_t component,
                            FrozenConvectionFaceOutput output) noexcept {
  return preflight_frozen_axis<Uniform, 0U, Scheme>(
             plan, flux, transported, component, output.x) &&
         preflight_frozen_axis<Uniform, 1U, Scheme>(
             plan, flux, transported, component, output.y) &&
         preflight_frozen_axis<Uniform, 2U, Scheme>(
             plan, flux, transported, component, output.z);
}

template <bool Uniform, ConvectionScheme Scheme>
void commit_frozen_faces(const CartesianKernelPlan& plan,
                         ConstFaceFluxView flux,
                         ConstFieldView transported,
                         std::uint8_t component,
                         FrozenConvectionFaceOutput output) noexcept {
  commit_frozen_axis<Uniform, 0U, Scheme>(plan, flux, transported, component,
                                          output.x);
  commit_frozen_axis<Uniform, 1U, Scheme>(plan, flux, transported, component,
                                          output.y);
  commit_frozen_axis<Uniform, 2U, Scheme>(plan, flux, transported, component,
                                          output.z);
}

template <ConvectionScheme Scheme>
bool dispatch_frozen_preflight(const CartesianKernelPlan& plan,
                               ConstFaceFluxView flux,
                               ConstFieldView transported,
                               std::uint8_t component,
                               FrozenConvectionFaceOutput output) noexcept {
  return plan.geometry_kind() == GeometryKind::uniform
             ? preflight_frozen_faces<true, Scheme>(
                   plan, flux, transported, component, output)
             : preflight_frozen_faces<false, Scheme>(
                   plan, flux, transported, component, output);
}

template <ConvectionScheme Scheme>
void dispatch_frozen_commit(const CartesianKernelPlan& plan,
                            ConstFaceFluxView flux,
                            ConstFieldView transported,
                            std::uint8_t component,
                            FrozenConvectionFaceOutput output) noexcept {
  if (plan.geometry_kind() == GeometryKind::uniform) {
    commit_frozen_faces<true, Scheme>(plan, flux, transported, component,
                                      output);
  } else {
    commit_frozen_faces<false, Scheme>(plan, flux, transported, component,
                                       output);
  }
}

std::uint64_t frozen_reconstruction_identity(
    ConvectionScheme scheme, FrozenConvectionContext context) noexcept {
  std::uint64_t reconstruction = hash_mix(kFnvOffset, kFrozenConvectionSchema);
  reconstruction = hash_mix(reconstruction, context.collective_semantics);
  reconstruction = hash_mix(reconstruction, static_cast<std::uint8_t>(scheme));
  return nonzero_hash(reconstruction);
}

std::uint64_t frozen_revision_identity(
    std::uint64_t reconstruction, FrozenConvectionContext context,
    ConstFaceFluxView target_flux, ConstFieldView transported,
    std::uint8_t component) noexcept {
  std::uint64_t revision = hash_mix(kFnvOffset, reconstruction);
  revision = hash_mix(revision, context.closure);
  revision = hash_mix(revision, target_flux.revision);
  revision = hash_mix(revision, transported.field);
  revision = hash_mix(revision, transported.revision);
  revision = hash_mix(revision, component);
  return nonzero_hash(revision);
}

template <class OutputFace>
std::uint64_t frozen_local_binding_identity(
    std::uint64_t revision, std::uint64_t reconstruction,
    const CartesianKernelPlan& plan, ConstFaceFluxView target_flux,
    ConstFieldView transported,
    const std::array<OutputFace, 3U>& output_faces) noexcept {
  const ConstFaceFieldView flux_faces[]{target_flux.x, target_flux.y,
                                       target_flux.z};
  std::uint64_t local = hash_mix(kFnvOffset, revision);
  local = hash_mix(local, reconstruction);
  local = hash_mix(local, plan.fingerprint());
  local = mix_field_view(local, transported);
  for (ConstFaceFieldView face : flux_faces)
    local = mix_face_view(local, face);
  for (OutputFace face : output_faces) local = mix_face_view(local, face);
  return nonzero_hash(local);
}

bool valid_linearization_policy(
    FrozenConvectionLinearizationPolicy policy) noexcept {
  return policy ==
             FrozenConvectionLinearizationPolicy::classical_active_branch ||
         policy == FrozenConvectionLinearizationPolicy::
                       semismooth_generalized_zero_slope;
}

std::uint64_t directional_reconstruction_identity(
    std::uint64_t frozen_reconstruction,
    FrozenConvectionLinearizationPolicy policy) noexcept {
  std::uint64_t identity =
      hash_mix(kFnvOffset, kFrozenConvectionDerivativeSchema);
  identity = hash_mix(identity, frozen_reconstruction);
  identity = hash_mix(identity, static_cast<std::uint8_t>(policy));
  return nonzero_hash(identity);
}

enum class DirectionalPreflight : std::uint8_t {
  success,
  nonfinite,
  nondifferentiable,
  stale_numeric
};

struct DirectionalPreflightState {
  DirectionalPreflight result{DirectionalPreflight::success};
  std::uint64_t branch_hash{};
  std::uint64_t generalized_face_count{};
};

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
void preflight_direction_axis(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    ConstFieldView variation, std::uint8_t variation_component,
    FrozenConvectionLinearizationPolicy policy, ConstFaceFieldView frozen,
    FaceFieldView output, DirectionalPreflightState& state) noexcept {
  state.branch_hash = hash_mix(state.branch_hash, Axis);
  for (std::int32_t z = 0; z < output.extents.z; ++z) {
    for (std::int32_t y = 0; y < output.extents.y; ++y) {
      for (std::int32_t x = 0; x < output.extents.x; ++x) {
        const Int3 face{x, y, z};
        const double rate = face_rate<Axis>(flux, face);
        const double target_value = reconstructed_face<Uniform, Axis, Scheme>(
            plan, target, target_component, face, rate);
        const DirectionalFaceEvaluation directional =
            reconstructed_face_direction<Uniform, Axis, Scheme>(
                plan, target, target_component, variation,
                variation_component, face, rate, policy);
        if (!std::isfinite(rate) || !std::isfinite(target_value) ||
            !std::isfinite(frozen.unchecked(face))) {
          state.result = DirectionalPreflight::nonfinite;
          return;
        }
        if (!same_bits(target_value, frozen.unchecked(face))) {
          state.result = DirectionalPreflight::stale_numeric;
          return;
        }
        if (!directional.differentiable) {
          state.result = DirectionalPreflight::nondifferentiable;
          return;
        }
        if (!std::isfinite(directional.value)) {
          state.result = DirectionalPreflight::nonfinite;
          return;
        }
        state.branch_hash = hash_mix(state.branch_hash,
                                     directional.branch_code);
        state.branch_hash = hash_mix(state.branch_hash,
                                     directional.generalized ? 1U : 0U);
        if (directional.generalized) {
          if (state.generalized_face_count ==
              std::numeric_limits<std::uint64_t>::max()) {
            state.result = DirectionalPreflight::nonfinite;
            return;
          }
          ++state.generalized_face_count;
        }
      }
    }
  }
}

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
void commit_direction_axis(const CartesianKernelPlan& plan,
                           ConstFaceFluxView flux, ConstFieldView target,
                           std::uint8_t target_component,
                           ConstFieldView variation,
                           std::uint8_t variation_component,
                           FrozenConvectionLinearizationPolicy policy,
                           FaceFieldView output) noexcept {
  for (std::int32_t z = 0; z < output.extents.z; ++z) {
    for (std::int32_t y = 0; y < output.extents.y; ++y) {
      for (std::int32_t x = 0; x < output.extents.x; ++x) {
        const Int3 face{x, y, z};
        output.unchecked(face) =
            reconstructed_face_direction<Uniform, Axis, Scheme>(
                plan, target, target_component, variation,
                variation_component, face, face_rate<Axis>(flux, face), policy)
                .value;
      }
    }
  }
}

template <bool Uniform, ConvectionScheme Scheme>
void preflight_direction_faces(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    ConstFieldView variation, std::uint8_t variation_component,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen,
    FrozenConvectionFaceOutput output,
    DirectionalPreflightState& state) noexcept {
  preflight_direction_axis<Uniform, 0U, Scheme>(
      plan, flux, target, target_component, variation, variation_component,
      policy, frozen.x, output.x, state);
  if (state.result == DirectionalPreflight::success)
    preflight_direction_axis<Uniform, 1U, Scheme>(
        plan, flux, target, target_component, variation, variation_component,
        policy, frozen.y, output.y, state);
  if (state.result == DirectionalPreflight::success)
    preflight_direction_axis<Uniform, 2U, Scheme>(
        plan, flux, target, target_component, variation, variation_component,
        policy, frozen.z, output.z, state);
}

template <bool Uniform, ConvectionScheme Scheme>
void commit_direction_faces(const CartesianKernelPlan& plan,
                            ConstFaceFluxView flux, ConstFieldView target,
                            std::uint8_t target_component,
                            ConstFieldView variation,
                            std::uint8_t variation_component,
                            FrozenConvectionLinearizationPolicy policy,
                            FrozenConvectionFaceOutput output) noexcept {
  commit_direction_axis<Uniform, 0U, Scheme>(
      plan, flux, target, target_component, variation, variation_component,
      policy, output.x);
  commit_direction_axis<Uniform, 1U, Scheme>(
      plan, flux, target, target_component, variation, variation_component,
      policy, output.y);
  commit_direction_axis<Uniform, 2U, Scheme>(
      plan, flux, target, target_component, variation, variation_component,
      policy, output.z);
}

template <ConvectionScheme Scheme>
void dispatch_direction_preflight(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    ConstFieldView variation, std::uint8_t variation_component,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen, FrozenConvectionFaceOutput output,
    DirectionalPreflightState& state) noexcept {
  if (plan.geometry_kind() == GeometryKind::uniform) {
    preflight_direction_faces<true, Scheme>(
        plan, flux, target, target_component, variation, variation_component,
        policy, frozen, output, state);
  } else {
    preflight_direction_faces<false, Scheme>(
        plan, flux, target, target_component, variation, variation_component,
        policy, frozen, output, state);
  }
}

template <ConvectionScheme Scheme>
void dispatch_direction_commit(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    ConstFieldView variation, std::uint8_t variation_component,
    FrozenConvectionLinearizationPolicy policy,
    FrozenConvectionFaceOutput output) noexcept {
  if (plan.geometry_kind() == GeometryKind::uniform) {
    commit_direction_faces<true, Scheme>(
        plan, flux, target, target_component, variation, variation_component,
        policy, output);
  } else {
    commit_direction_faces<false, Scheme>(
        plan, flux, target, target_component, variation, variation_component,
        policy, output);
  }
}

constexpr std::uint16_t kCompiledBranchGeneralized = UINT16_C(0x8000);

std::size_t face_value_count(Int3 extents) noexcept {
  return static_cast<std::size_t>(extents.x) *
         static_cast<std::size_t>(extents.y) *
         static_cast<std::size_t>(extents.z);
}

std::size_t compiled_branch_count(Int3 cells) noexcept {
  return face_value_count(expected_face_extents(cells, CartesianAxis::x)) +
         face_value_count(expected_face_extents(cells, CartesianAxis::y)) +
         face_value_count(expected_face_extents(cells, CartesianAxis::z));
}

template <std::size_t Axis>
std::size_t compiled_branch_offset(Int3 cells, Int3 face) noexcept {
  const CartesianAxis axis = static_cast<CartesianAxis>(Axis);
  const Int3 extents = expected_face_extents(cells, axis);
  std::size_t base = 0U;
  if constexpr (Axis >= 1U)
    base += face_value_count(
        expected_face_extents(cells, CartesianAxis::x));
  if constexpr (Axis >= 2U)
    base += face_value_count(
        expected_face_extents(cells, CartesianAxis::y));
  return base + static_cast<std::size_t>(face.x) +
         static_cast<std::size_t>(extents.x) *
             (static_cast<std::size_t>(face.y) +
              static_cast<std::size_t>(extents.y) *
                  static_cast<std::size_t>(face.z));
}

template <bool Uniform, std::size_t Axis, bool Commit>
void compile_limited_branch_axis(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionLinearizationPolicy policy, ConstFaceFieldView frozen,
    std::uint16_t* output_values,
    DirectionalPreflightState& state) noexcept {
  state.branch_hash = hash_mix(state.branch_hash, Axis);
  for (std::int32_t z = 0; z < frozen.extents.z; ++z) {
    for (std::int32_t y = 0; y < frozen.extents.y; ++y) {
      for (std::int32_t x = 0; x < frozen.extents.x; ++x) {
        const Int3 face{x, y, z};
        const double rate = face_rate<Axis>(flux, face);
        const double target_value =
            reconstructed_face<Uniform, Axis,
                               ConvectionScheme::limited_central2>(
                plan, target, target_component, face, rate);
        const FaceBranchSelection selection =
            select_face_branches<Uniform, Axis,
                                 ConvectionScheme::limited_central2>(
                plan, target, target_component, face, rate);
        const bool generalized = !selection.differentiable;
        if (!std::isfinite(rate) || !std::isfinite(target_value) ||
            !std::isfinite(frozen.unchecked(face))) {
          state.result = DirectionalPreflight::nonfinite;
          return;
        }
        if (!same_bits(target_value, frozen.unchecked(face))) {
          state.result = DirectionalPreflight::stale_numeric;
          return;
        }
        if (generalized &&
            policy != FrozenConvectionLinearizationPolicy::
                          semismooth_generalized_zero_slope) {
          state.result = DirectionalPreflight::nondifferentiable;
          return;
        }
        if (selection.code >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::uint16_t>::max() &
                ~kCompiledBranchGeneralized)) {
          state.result = DirectionalPreflight::nonfinite;
          return;
        }
        state.branch_hash = hash_mix(state.branch_hash, selection.code);
        state.branch_hash = hash_mix(state.branch_hash, generalized ? 1U : 0U);
        if (generalized) {
          if (state.generalized_face_count ==
              std::numeric_limits<std::uint64_t>::max()) {
            state.result = DirectionalPreflight::nonfinite;
            return;
          }
          ++state.generalized_face_count;
        }
        if constexpr (Commit) {
          output_values[compiled_branch_offset<Axis>(plan.cells(), face)] =
              static_cast<std::uint16_t>(selection.code) |
              (generalized ? kCompiledBranchGeneralized : UINT16_C(0));
        }
      }
    }
  }
}

template <bool Uniform, bool Commit>
void compile_limited_branch_faces(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen,
    std::uint16_t* output_values,
    DirectionalPreflightState& state) noexcept {
  compile_limited_branch_axis<Uniform, 0U, Commit>(
      plan, flux, target, target_component, policy, frozen.x, output_values,
      state);
  if (state.result == DirectionalPreflight::success)
    compile_limited_branch_axis<Uniform, 1U, Commit>(
        plan, flux, target, target_component, policy, frozen.y, output_values,
        state);
  if (state.result == DirectionalPreflight::success)
    compile_limited_branch_axis<Uniform, 2U, Commit>(
        plan, flux, target, target_component, policy, frozen.z, output_values,
        state);
}

template <bool Uniform, std::size_t Axis>
bool apply_limited_branch_axis(
    const CartesianKernelPlan& plan,
    const FrozenConvectionBranchPlan& branches, ConstFieldView variation,
    std::uint8_t variation_component, FaceFieldView output) noexcept {
  for (std::int32_t z = 0; z < output.extents.z; ++z) {
    for (std::int32_t y = 0; y < output.extents.y; ++y) {
      for (std::int32_t x = 0; x < output.extents.x; ++x) {
        const Int3 face{x, y, z};
        const std::uint16_t encoded =
            branches.values.data[compiled_branch_offset<Axis>(plan.cells(),
                                                               face)];
        const std::uint16_t branch_code =
            encoded & static_cast<std::uint16_t>(~kCompiledBranchGeneralized);
        const bool generalized =
            (encoded & kCompiledBranchGeneralized) != 0U;
        const LimitedSlopeBranch left_branch =
            static_cast<LimitedSlopeBranch>(branch_code & UINT16_C(0x0f));
        const LimitedSlopeBranch right_branch = static_cast<LimitedSlopeBranch>(
            (branch_code >> 4U) & UINT16_C(0x0f));
        if ((branch_code & UINT16_C(0xff00)) != UINT16_C(0x0200) ||
            static_cast<std::uint8_t>(left_branch) <
                static_cast<std::uint8_t>(LimitedSlopeBranch::zero) ||
            static_cast<std::uint8_t>(left_branch) >
                static_cast<std::uint8_t>(
                    LimitedSlopeBranch::nondifferentiable) ||
            static_cast<std::uint8_t>(right_branch) <
                static_cast<std::uint8_t>(LimitedSlopeBranch::zero) ||
            static_cast<std::uint8_t>(right_branch) >
                static_cast<std::uint8_t>(
                    LimitedSlopeBranch::nondifferentiable)) {
          return false;
        }
        const Int3 left = offset_axis<Axis>(face, -1);
        const Int3 right = face;
        const double v_left = variation.unchecked(left, variation_component);
        const double v_right = variation.unchecked(right, variation_component);
        double value = 0.0;
        if (generalized) {
          if (branches.policy != FrozenConvectionLinearizationPolicy::
                                     semismooth_generalized_zero_slope ||
              left_branch != LimitedSlopeBranch::nondifferentiable &&
                  right_branch != LimitedSlopeBranch::nondifferentiable) {
            return false;
          }
          value = 0.5 * (v_left + v_right);
        } else {
          if (left_branch == LimitedSlopeBranch::nondifferentiable ||
              right_branch == LimitedSlopeBranch::nondifferentiable) {
            return false;
          }
          const std::int32_t normal = axis_index<Axis>(face);
          const double face_coordinate =
              detail::metric_face<Uniform>(plan, Axis, normal);
          const double left_slope =
              directional_limited_slope_values<Uniform, Axis>(
                  plan, normal - 1,
                  variation.unchecked(offset_axis<Axis>(left, -1),
                                      variation_component),
                  v_left, v_right, left_branch);
          const double right_slope =
              directional_limited_slope_values<Uniform, Axis>(
                  plan, normal, v_left, v_right,
                  variation.unchecked(offset_axis<Axis>(right, 1),
                                      variation_component),
                  right_branch);
          const double left_reconstructed =
              v_left + left_slope *
                           (face_coordinate -
                            detail::metric_centre<Uniform>(plan, Axis,
                                                           normal - 1));
          const double right_reconstructed =
              v_right + right_slope *
                            (face_coordinate -
                             detail::metric_centre<Uniform>(plan, Axis,
                                                            normal));
          value = 0.5 * (left_reconstructed + right_reconstructed);
        }
        if (!std::isfinite(value)) return false;
        output.unchecked(face) = value;
      }
    }
  }
  return true;
}

template <bool Uniform, std::size_t Axis>
bool verify_limited_branch_axis(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionBranchPlan& branches,
    ConstFaceFieldView frozen) noexcept {
  for (std::int32_t z = 0; z < frozen.extents.z; ++z) {
    for (std::int32_t y = 0; y < frozen.extents.y; ++y) {
      for (std::int32_t x = 0; x < frozen.extents.x; ++x) {
        const Int3 face{x, y, z};
        const double rate = face_rate<Axis>(flux, face);
        const FaceBranchSelection selection =
            select_face_branches<Uniform, Axis,
                                 ConvectionScheme::limited_central2>(
                plan, target, target_component, face, rate);
        const bool generalized = !selection.differentiable;
        if (generalized &&
            policy != FrozenConvectionLinearizationPolicy::
                          semismooth_generalized_zero_slope) {
          return false;
        }
        const std::uint16_t expected =
            static_cast<std::uint16_t>(selection.code) |
            (generalized ? kCompiledBranchGeneralized : UINT16_C(0));
        if (branches.values.data[compiled_branch_offset<Axis>(plan.cells(),
                                                              face)] !=
            expected) {
          return false;
        }
      }
    }
  }
  return true;
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

struct PointDiagnosticAccumulator {
  double convection_integral{};
  double mass_integral{};
  double maximum_violation{};
  double envelope_scale{1.0};
  double selected_face{};
  double selected_minimum{};
  double selected_maximum{};
  bool sampled{};
};

template <bool Uniform, std::size_t Axis, ConvectionScheme Scheme>
Status accumulate_diagnostic_face(const CartesianKernelPlan& plan,
                                  ConstFaceFluxView flux,
                                  ConstFieldView transported,
                                  std::uint8_t component, Int3 face,
                                  double sign,
                                  PointDiagnosticAccumulator& out) noexcept {
  const double rate = face_rate<Axis>(flux, face);
  const double value = reconstructed_face<Uniform, Axis, Scheme>(
      plan, transported, component, face, rate);
  const Int3 left = offset_axis<Axis>(face, -1);
  const double left_value = transported.unchecked(left, component);
  const double right_value = transported.unchecked(face, component);
  const double donor_minimum = std::min(left_value, right_value);
  const double donor_maximum = std::max(left_value, right_value);
  const double violation =
      std::max({0.0, donor_minimum - value, value - donor_maximum});
  if (!std::isfinite(rate) || !std::isfinite(value) ||
      !std::isfinite(donor_minimum) || !std::isfinite(donor_maximum) ||
      !std::isfinite(violation)) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }
  out.convection_integral += sign * rate * value;
  out.mass_integral += sign * rate;
  out.envelope_scale =
      std::max({out.envelope_scale, std::abs(value),
                std::abs(donor_minimum), std::abs(donor_maximum)});
  if (!out.sampled || violation > out.maximum_violation) {
    out.maximum_violation = violation;
    out.selected_face = value;
    out.selected_minimum = donor_minimum;
    out.selected_maximum = donor_maximum;
    out.sampled = true;
  }
  return {};
}

template <bool Uniform, ConvectionScheme Scheme>
Status diagnose_convection_geometry(const CartesianKernelPlan& plan,
                                    ConstFaceFluxView flux,
                                    ConstFieldView transported,
                                    std::uint8_t component, Int3 cell,
                                    ConvectionPointDiagnostic& out) noexcept {
  PointDiagnosticAccumulator accumulated;
  Status status = accumulate_diagnostic_face<Uniform, 0U, Scheme>(
      plan, flux, transported, component, cell, -1.0, accumulated);
  if (status)
    status = accumulate_diagnostic_face<Uniform, 0U, Scheme>(
        plan, flux, transported, component,
        {cell.x + 1, cell.y, cell.z}, 1.0, accumulated);
  if (status)
    status = accumulate_diagnostic_face<Uniform, 1U, Scheme>(
        plan, flux, transported, component, cell, -1.0, accumulated);
  if (status)
    status = accumulate_diagnostic_face<Uniform, 1U, Scheme>(
        plan, flux, transported, component,
        {cell.x, cell.y + 1, cell.z}, 1.0, accumulated);
  if (status)
    status = accumulate_diagnostic_face<Uniform, 2U, Scheme>(
        plan, flux, transported, component, cell, -1.0, accumulated);
  if (status)
    status = accumulate_diagnostic_face<Uniform, 2U, Scheme>(
        plan, flux, transported, component,
        {cell.x, cell.y, cell.z + 1}, 1.0, accumulated);
  if (!status) return status;

  const double inverse_volume =
      detail::metric_inverse_volume<Uniform>(plan, cell);
  ConvectionPointDiagnostic candidate;
  candidate.divergence = accumulated.convection_integral * inverse_volume;
  candidate.mass_divergence = accumulated.mass_integral * inverse_volume;
  candidate.maximum_face_envelope_violation =
      accumulated.maximum_violation;
  candidate.selected_face_value = accumulated.selected_face;
  candidate.selected_donor_minimum = accumulated.selected_minimum;
  candidate.selected_donor_maximum = accumulated.selected_maximum;
  candidate.face_envelope_checked = Scheme == ConvectionScheme::tvd2;
  candidate.face_envelope_valid =
      !candidate.face_envelope_checked ||
      candidate.maximum_face_envelope_violation <=
          128.0 * std::numeric_limits<double>::epsilon() *
              accumulated.envelope_scale;
  if (!std::isfinite(candidate.divergence) ||
      !std::isfinite(candidate.mass_divergence)) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }
  out = candidate;
  return {};
}

}  // namespace

Status freeze_cartesian_target_convection_faces(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView target_flux, ConstFieldView transported,
    std::uint8_t component, FrozenConvectionContext context,
    FrozenConvectionFaceOutput output,
    FrozenConvectionFaceField& frozen) noexcept {
  const Int3 cells = plan.cells();
  const std::uint8_t required_ghost_width =
      scheme == ConvectionScheme::central2 ? 1U : 2U;
  const bool valid_scheme =
      static_cast<std::uint8_t>(scheme) <=
      static_cast<std::uint8_t>(ConvectionScheme::tvd2);
  const ConstFaceFieldView flux_faces[]{target_flux.x, target_flux.y,
                                       target_flux.z};
  const FaceFieldView output_faces[]{output.x, output.y, output.z};
  bool aliases = false;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    aliases = aliases ||
              detail::cell_face_views_overlap(transported,
                                               output_faces[axis]);
    for (std::size_t flux_axis = 0U; flux_axis < 3U; ++flux_axis)
      aliases = aliases || detail::face_views_overlap(
                               output_faces[axis], flux_faces[flux_axis]);
    for (std::size_t prior = 0U; prior < axis; ++prior)
      aliases = aliases || detail::face_views_overlap(
                               output_faces[axis], output_faces[prior]);
  }
  if (plan.fingerprint() == 0U || !valid_scheme ||
      context.collective_semantics == 0U || context.closure == 0U ||
      !detail::valid_cell_view(transported, cells, component, 1U,
                               required_ghost_width) ||
      !detail::valid_flux_view(target_flux, cells, target_flux.revision) ||
      target_flux.certificate.valid() ||
      !valid_frozen_output(output.x, CartesianAxis::x, cells) ||
      !valid_frozen_output(output.y, CartesianAxis::y, cells) ||
      !valid_frozen_output(output.z, CartesianAxis::z, cells) ||
      output.x.storage_identity != output.y.storage_identity ||
      output.x.storage_identity != output.z.storage_identity ||
      output.x.revision_domain != output.y.revision_domain ||
      output.x.revision_domain != output.z.revision_domain || aliases) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }

  bool finite = false;
  switch (scheme) {
    case ConvectionScheme::central2:
      finite = dispatch_frozen_preflight<ConvectionScheme::central2>(
          plan, target_flux, transported, component, output);
      break;
    case ConvectionScheme::limited_central2:
      finite =
          dispatch_frozen_preflight<ConvectionScheme::limited_central2>(
              plan, target_flux, transported, component, output);
      break;
    case ConvectionScheme::tvd2:
      finite = dispatch_frozen_preflight<ConvectionScheme::tvd2>(
          plan, target_flux, transported, component, output);
      break;
  }
  if (!finite) return {StatusCode::numerical_failure, kTransportNumerical};

  switch (scheme) {
    case ConvectionScheme::central2:
      dispatch_frozen_commit<ConvectionScheme::central2>(
          plan, target_flux, transported, component, output);
      break;
    case ConvectionScheme::limited_central2:
      dispatch_frozen_commit<ConvectionScheme::limited_central2>(
          plan, target_flux, transported, component, output);
      break;
    case ConvectionScheme::tvd2:
      dispatch_frozen_commit<ConvectionScheme::tvd2>(
          plan, target_flux, transported, component, output);
      break;
  }

  const std::uint64_t reconstruction =
      frozen_reconstruction_identity(scheme, context);
  const std::uint64_t revision = frozen_revision_identity(
      reconstruction, context, target_flux, transported, component);
  const std::array<FaceFieldView, 3U> local_output_faces{
      output.x, output.y, output.z};
  const std::uint64_t local = frozen_local_binding_identity(
      revision, reconstruction, plan, target_flux, transported,
      local_output_faces);

  FrozenConvectionFaceField candidate{
      as_const(output.x), as_const(output.y), as_const(output.z), revision,
      reconstruction, local, true};
  frozen = candidate;
  return {};
}

Status differentiate_frozen_cartesian_target_convection_faces(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView target_flux, ConstFieldView target,
    std::uint8_t target_component, FrozenConvectionContext context,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen, ConstFieldView variation,
    std::uint8_t variation_component, FrozenConvectionFaceOutput output,
    FrozenConvectionFaceDirectionalDerivative& derivative) noexcept {
  const Int3 cells = plan.cells();
  const bool valid_scheme =
      static_cast<std::uint8_t>(scheme) <=
      static_cast<std::uint8_t>(ConvectionScheme::tvd2);
  const std::uint8_t required_ghost_width =
      scheme == ConvectionScheme::central2 ? 1U : 2U;
  const ConstFaceFieldView flux_faces[]{target_flux.x, target_flux.y,
                                       target_flux.z};
  const ConstFaceFieldView frozen_faces[]{frozen.x, frozen.y, frozen.z};
  const FaceFieldView output_faces[]{output.x, output.y, output.z};
  bool aliases = false;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    aliases = aliases ||
              detail::cell_face_views_overlap(target, output_faces[axis]) ||
              detail::cell_face_views_overlap(variation,
                                               output_faces[axis]);
    for (std::size_t input_axis = 0U; input_axis < 3U; ++input_axis) {
      aliases = aliases || detail::face_views_overlap(
                               output_faces[axis], flux_faces[input_axis]) ||
                detail::face_views_overlap(output_faces[axis],
                                           frozen_faces[input_axis]);
    }
    for (std::size_t prior = 0U; prior < axis; ++prior)
      aliases = aliases || detail::face_views_overlap(
                               output_faces[axis], output_faces[prior]);
  }
  if (plan.fingerprint() == 0U || !valid_scheme ||
      !valid_linearization_policy(policy) ||
      context.collective_semantics == 0U || context.closure == 0U ||
      !detail::valid_cell_view(target, cells, target_component, 1U,
                               required_ghost_width) ||
      !detail::valid_cell_view(variation, cells, variation_component, 1U,
                               required_ghost_width) ||
      !detail::valid_flux_view(target_flux, cells, target_flux.revision) ||
      target_flux.certificate.valid() || !frozen.valid() ||
      !valid_frozen_output(frozen.x, CartesianAxis::x, cells) ||
      !valid_frozen_output(frozen.y, CartesianAxis::y, cells) ||
      !valid_frozen_output(frozen.z, CartesianAxis::z, cells) ||
      frozen.x.storage_identity != frozen.y.storage_identity ||
      frozen.x.storage_identity != frozen.z.storage_identity ||
      frozen.x.revision_domain != frozen.y.revision_domain ||
      frozen.x.revision_domain != frozen.z.revision_domain ||
      !valid_frozen_output(output.x, CartesianAxis::x, cells) ||
      !valid_frozen_output(output.y, CartesianAxis::y, cells) ||
      !valid_frozen_output(output.z, CartesianAxis::z, cells) ||
      output.x.storage_identity != output.y.storage_identity ||
      output.x.storage_identity != output.z.storage_identity ||
      output.x.revision_domain != output.y.revision_domain ||
      output.x.revision_domain != output.z.revision_domain || aliases) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }

  const std::uint64_t frozen_reconstruction =
      frozen_reconstruction_identity(scheme, context);
  const std::uint64_t frozen_revision = frozen_revision_identity(
      frozen_reconstruction, context, target_flux, target, target_component);
  const std::array<ConstFaceFieldView, 3U> frozen_output_faces{
      frozen.x, frozen.y, frozen.z};
  const std::uint64_t frozen_local = frozen_local_binding_identity(
      frozen_revision, frozen_reconstruction, plan, target_flux, target,
      frozen_output_faces);
  if (frozen.revision != frozen_revision ||
      frozen.reconstruction != frozen_reconstruction ||
      frozen.local_binding != frozen_local ||
      !frozen.exact_target_reconstruction) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }

  const std::uint64_t reconstruction =
      directional_reconstruction_identity(frozen_reconstruction, policy);
  DirectionalPreflightState preflight;
  preflight.branch_hash = hash_mix(kFnvOffset, reconstruction);
  preflight.branch_hash = hash_mix(preflight.branch_hash, frozen_revision);
  switch (scheme) {
    case ConvectionScheme::central2:
      dispatch_direction_preflight<ConvectionScheme::central2>(
          plan, target_flux, target, target_component, variation,
          variation_component, policy, frozen, output, preflight);
      break;
    case ConvectionScheme::limited_central2:
      dispatch_direction_preflight<
          ConvectionScheme::limited_central2>(
          plan, target_flux, target, target_component, variation,
          variation_component, policy, frozen, output, preflight);
      break;
    case ConvectionScheme::tvd2:
      dispatch_direction_preflight<ConvectionScheme::tvd2>(
          plan, target_flux, target, target_component, variation,
          variation_component, policy, frozen, output, preflight);
      break;
  }
  if (preflight.result == DirectionalPreflight::stale_numeric) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  if (preflight.result == DirectionalPreflight::nondifferentiable) {
    return {StatusCode::numerical_failure, kTransportNondifferentiable};
  }
  if (preflight.result != DirectionalPreflight::success) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }

  switch (scheme) {
    case ConvectionScheme::central2:
      dispatch_direction_commit<ConvectionScheme::central2>(
          plan, target_flux, target, target_component, variation,
          variation_component, policy, output);
      break;
    case ConvectionScheme::limited_central2:
      dispatch_direction_commit<ConvectionScheme::limited_central2>(
          plan, target_flux, target, target_component, variation,
          variation_component, policy, output);
      break;
    case ConvectionScheme::tvd2:
      dispatch_direction_commit<ConvectionScheme::tvd2>(
          plan, target_flux, target, target_component, variation,
          variation_component, policy, output);
      break;
  }

  const std::uint64_t branch_authority =
      nonzero_hash(preflight.branch_hash);
  std::uint64_t derivative_revision =
      hash_mix(kFnvOffset, kFrozenConvectionDerivativeSchema);
  derivative_revision = hash_mix(derivative_revision, reconstruction);
  derivative_revision = hash_mix(derivative_revision, frozen_revision);
  derivative_revision = hash_mix(derivative_revision, variation.field);
  derivative_revision = hash_mix(derivative_revision, variation.revision);
  derivative_revision = hash_mix(derivative_revision, variation_component);
  derivative_revision = nonzero_hash(derivative_revision);
  std::uint64_t local = hash_mix(kFnvOffset, derivative_revision);
  local = hash_mix(local, reconstruction);
  local = hash_mix(local, branch_authority);
  local = hash_mix(local, static_cast<std::uint8_t>(policy));
  local = hash_mix(local, preflight.generalized_face_count);
  local = hash_mix(local, plan.fingerprint());
  local = mix_field_view(local, target);
  local = mix_field_view(local, variation);
  for (ConstFaceFieldView face : flux_faces)
    local = mix_face_view(local, face);
  for (ConstFaceFieldView face : frozen_faces)
    local = mix_face_view(local, face);
  for (FaceFieldView face : output_faces)
    local = mix_face_view(local, face);
  local = nonzero_hash(local);
  const bool classical_everywhere =
      preflight.generalized_face_count == 0U;
  derivative = {as_const(output.x), as_const(output.y), as_const(output.z),
                derivative_revision, reconstruction, branch_authority, local,
                policy, preflight.generalized_face_count,
                classical_everywhere};
  return {};
}

Status compile_frozen_limited_central2_branches(
    const CartesianKernelPlan& plan, ConstFaceFluxView target_flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionContext context,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen,
    FrozenConvectionBranchOutput output,
    FrozenConvectionBranchPlan& branches) noexcept {
  branches = {};
  const Int3 cells = plan.cells();
  const ConstFaceFieldView flux_faces[]{target_flux.x, target_flux.y,
                                       target_flux.z};
  const ConstFaceFieldView frozen_faces[]{frozen.x, frozen.y, frozen.z};
  if (plan.fingerprint() == 0U || !valid_linearization_policy(policy) ||
      context.collective_semantics == 0U || context.closure == 0U ||
      !detail::valid_cell_view(target, cells, target_component, 1U, 2U) ||
      !detail::valid_flux_view(target_flux, cells, target_flux.revision) ||
      target_flux.certificate.valid() || !frozen.valid() ||
      !valid_frozen_output(frozen.x, CartesianAxis::x, cells) ||
      !valid_frozen_output(frozen.y, CartesianAxis::y, cells) ||
      !valid_frozen_output(frozen.z, CartesianAxis::z, cells) ||
      frozen.x.storage_identity != frozen.y.storage_identity ||
      frozen.x.storage_identity != frozen.z.storage_identity ||
      frozen.x.revision_domain != frozen.y.revision_domain ||
      frozen.x.revision_domain != frozen.z.revision_domain ||
      output.values.data == nullptr ||
      output.values.size != compiled_branch_count(cells)) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }

  const std::uint64_t frozen_reconstruction = frozen_reconstruction_identity(
      ConvectionScheme::limited_central2, context);
  const std::uint64_t frozen_revision = frozen_revision_identity(
      frozen_reconstruction, context, target_flux, target, target_component);
  const std::array<ConstFaceFieldView, 3U> frozen_output_faces{
      frozen.x, frozen.y, frozen.z};
  const std::uint64_t frozen_local = frozen_local_binding_identity(
      frozen_revision, frozen_reconstruction, plan, target_flux, target,
      frozen_output_faces);
  if (frozen.revision != frozen_revision ||
      frozen.reconstruction != frozen_reconstruction ||
      frozen.local_binding != frozen_local ||
      !frozen.exact_target_reconstruction) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }

  const std::uint64_t reconstruction =
      directional_reconstruction_identity(frozen_reconstruction, policy);
  DirectionalPreflightState preflight;
  preflight.branch_hash = hash_mix(kFnvOffset, reconstruction);
  preflight.branch_hash = hash_mix(preflight.branch_hash, frozen_revision);
  if (plan.geometry_kind() == GeometryKind::uniform) {
    compile_limited_branch_faces<true, false>(
        plan, target_flux, target, target_component, policy, frozen, nullptr,
        preflight);
  } else {
    compile_limited_branch_faces<false, false>(
        plan, target_flux, target, target_component, policy, frozen, nullptr,
        preflight);
  }
  if (preflight.result == DirectionalPreflight::stale_numeric)
    return {StatusCode::invalid_plan, kTransportKernel};
  if (preflight.result == DirectionalPreflight::nondifferentiable)
    return {StatusCode::numerical_failure, kTransportNondifferentiable};
  if (preflight.result != DirectionalPreflight::success)
    return {StatusCode::numerical_failure, kTransportNumerical};

  DirectionalPreflightState commit;
  commit.branch_hash = hash_mix(kFnvOffset, reconstruction);
  commit.branch_hash = hash_mix(commit.branch_hash, frozen_revision);
  if (plan.geometry_kind() == GeometryKind::uniform) {
    compile_limited_branch_faces<true, true>(
        plan, target_flux, target, target_component, policy, frozen,
        output.values.data, commit);
  } else {
    compile_limited_branch_faces<false, true>(
        plan, target_flux, target, target_component, policy, frozen,
        output.values.data, commit);
  }
  if (commit.result != DirectionalPreflight::success ||
      commit.branch_hash != preflight.branch_hash ||
      commit.generalized_face_count != preflight.generalized_face_count) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }

  const std::uint64_t branch_authority = nonzero_hash(preflight.branch_hash);
  std::uint64_t revision =
      hash_mix(kFnvOffset, kFrozenConvectionDerivativeSchema);
  revision = hash_mix(revision, reconstruction);
  revision = hash_mix(revision, frozen_revision);
  revision = hash_mix(revision, target.field);
  revision = hash_mix(revision, target.revision);
  revision = nonzero_hash(revision);
  std::uint64_t local = hash_mix(kFnvOffset, revision);
  local = hash_mix(local, reconstruction);
  local = hash_mix(local, branch_authority);
  local = hash_mix(local, static_cast<std::uint8_t>(policy));
  local = hash_mix(local, preflight.generalized_face_count);
  local = hash_mix(local, plan.fingerprint());
  local = mix_field_view(local, target);
  for (ConstFaceFieldView face : flux_faces) local = mix_face_view(local, face);
  for (ConstFaceFieldView face : frozen_faces)
    local = mix_face_view(local, face);
  local = hash_mix(local,
                   reinterpret_cast<std::uintptr_t>(output.values.data));
  local = hash_mix(local, output.values.size);
  local = nonzero_hash(local);

  FrozenConvectionBranchPlan candidate{
      {output.values.data, output.values.size},
      cells,
      plan.fingerprint(),
      revision,
      reconstruction,
      branch_authority,
      local,
      policy,
      preflight.generalized_face_count,
      preflight.generalized_face_count == 0U};
  if (!candidate.valid())
    return {StatusCode::invalid_plan, kTransportKernel};
  branches = candidate;
  return {};
}

Status validate_frozen_limited_central2_branches(
    const CartesianKernelPlan& plan, ConstFaceFluxView target_flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionContext context,
    const FrozenConvectionFaceField& frozen,
    const FrozenConvectionBranchPlan& branches) noexcept {
  const Int3 cells = plan.cells();
  if (!branches.valid() || !same_shape(branches.cells, cells) ||
      branches.kernels != plan.fingerprint() ||
      branches.values.size != compiled_branch_count(cells) ||
      !detail::valid_cell_view(target, cells, target_component, 1U, 2U) ||
      !detail::valid_flux_view(target_flux, cells, target_flux.revision) ||
      target_flux.certificate.valid() || !frozen.valid()) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  const std::uint64_t frozen_reconstruction = frozen_reconstruction_identity(
      ConvectionScheme::limited_central2, context);
  const std::uint64_t frozen_revision = frozen_revision_identity(
      frozen_reconstruction, context, target_flux, target, target_component);
  const std::array<ConstFaceFieldView, 3U> frozen_output_faces{
      frozen.x, frozen.y, frozen.z};
  const std::uint64_t frozen_local = frozen_local_binding_identity(
      frozen_revision, frozen_reconstruction, plan, target_flux, target,
      frozen_output_faces);
  if (frozen.revision != frozen_revision ||
      frozen.reconstruction != frozen_reconstruction ||
      frozen.local_binding != frozen_local ||
      !frozen.exact_target_reconstruction) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  const std::uint64_t reconstruction = directional_reconstruction_identity(
      frozen_reconstruction, branches.policy);
  DirectionalPreflightState preflight;
  preflight.branch_hash = hash_mix(kFnvOffset, reconstruction);
  preflight.branch_hash = hash_mix(preflight.branch_hash, frozen_revision);
  if (plan.geometry_kind() == GeometryKind::uniform) {
    compile_limited_branch_faces<true, false>(
        plan, target_flux, target, target_component, branches.policy, frozen,
        nullptr, preflight);
  } else {
    compile_limited_branch_faces<false, false>(
        plan, target_flux, target, target_component, branches.policy, frozen,
        nullptr, preflight);
  }
  if (preflight.result == DirectionalPreflight::stale_numeric)
    return {StatusCode::invalid_plan, kTransportKernel};
  if (preflight.result == DirectionalPreflight::nondifferentiable)
    return {StatusCode::numerical_failure, kTransportNondifferentiable};
  if (preflight.result != DirectionalPreflight::success)
    return {StatusCode::numerical_failure, kTransportNumerical};
  const std::uint64_t branch_authority = nonzero_hash(preflight.branch_hash);
  std::uint64_t revision =
      hash_mix(kFnvOffset, kFrozenConvectionDerivativeSchema);
  revision = hash_mix(revision, reconstruction);
  revision = hash_mix(revision, frozen_revision);
  revision = hash_mix(revision, target.field);
  revision = hash_mix(revision, target.revision);
  revision = nonzero_hash(revision);
  const ConstFaceFieldView flux_faces[]{target_flux.x, target_flux.y,
                                       target_flux.z};
  const ConstFaceFieldView frozen_faces[]{frozen.x, frozen.y, frozen.z};
  std::uint64_t local = hash_mix(kFnvOffset, revision);
  local = hash_mix(local, reconstruction);
  local = hash_mix(local, branch_authority);
  local = hash_mix(local, static_cast<std::uint8_t>(branches.policy));
  local = hash_mix(local, preflight.generalized_face_count);
  local = hash_mix(local, plan.fingerprint());
  local = mix_field_view(local, target);
  for (ConstFaceFieldView face : flux_faces) local = mix_face_view(local, face);
  for (ConstFaceFieldView face : frozen_faces)
    local = mix_face_view(local, face);
  local = hash_mix(local,
                   reinterpret_cast<std::uintptr_t>(branches.values.data));
  local = hash_mix(local, branches.values.size);
  local = nonzero_hash(local);
  const bool identity = revision == branches.revision &&
                        reconstruction == branches.reconstruction &&
                        branch_authority == branches.branch_authority &&
                        local == branches.local_binding &&
                        preflight.generalized_face_count ==
                            branches.generalized_face_count &&
                        (preflight.generalized_face_count == 0U) ==
                            branches.classical_everywhere;
  if (!identity)
    return {StatusCode::invalid_plan, kTransportKernel};
  const bool verified = plan.geometry_kind() == GeometryKind::uniform
                            ? verify_limited_branch_axis<true, 0U>(
                                  plan, target_flux, target, target_component,
                                  branches.policy, branches, frozen.x) &&
                                  verify_limited_branch_axis<true, 1U>(
                                      plan, target_flux, target,
                                      target_component, branches.policy,
                                      branches, frozen.y) &&
                                  verify_limited_branch_axis<true, 2U>(
                                      plan, target_flux, target,
                                      target_component, branches.policy,
                                      branches, frozen.z)
                            : verify_limited_branch_axis<false, 0U>(
                                  plan, target_flux, target, target_component,
                                  branches.policy, branches, frozen.x) &&
                                  verify_limited_branch_axis<false, 1U>(
                                      plan, target_flux, target,
                                      target_component, branches.policy,
                                      branches, frozen.y) &&
                                  verify_limited_branch_axis<false, 2U>(
                                      plan, target_flux, target,
                                      target_component, branches.policy,
                                      branches, frozen.z);
  return verified ? Status{}
                  : Status{StatusCode::invalid_plan, kTransportKernel};
}

Status apply_frozen_limited_central2_branches(
    const CartesianKernelPlan& plan,
    const FrozenConvectionBranchPlan& branches, ConstFieldView variation,
    std::uint8_t variation_component,
    FrozenConvectionFaceOutput output) noexcept {
  const Int3 cells = plan.cells();
  const FaceFieldView outputs[]{output.x, output.y, output.z};
  bool aliases = false;
  for (FaceFieldView face : outputs)
    aliases = aliases || detail::cell_face_views_overlap(variation, face);
  for (std::size_t left = 0U; left < 3U; ++left)
    for (std::size_t right = left + 1U; right < 3U; ++right)
      aliases = aliases || detail::face_views_overlap(outputs[left],
                                                      outputs[right]);
  if (!branches.valid() || !same_shape(branches.cells, cells) ||
      branches.kernels != plan.fingerprint() ||
      branches.values.size != compiled_branch_count(cells) ||
      !detail::valid_cell_view(variation, cells, variation_component, 1U, 2U) ||
      !valid_frozen_output(output.x, CartesianAxis::x, cells) ||
      !valid_frozen_output(output.y, CartesianAxis::y, cells) ||
      !valid_frozen_output(output.z, CartesianAxis::z, cells) || aliases) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  const bool finite =
      plan.geometry_kind() == GeometryKind::uniform
          ? apply_limited_branch_axis<true, 0U>(
                plan, branches, variation, variation_component, output.x) &&
                apply_limited_branch_axis<true, 1U>(
                    plan, branches, variation, variation_component, output.y) &&
                apply_limited_branch_axis<true, 2U>(
                    plan, branches, variation, variation_component, output.z)
          : apply_limited_branch_axis<false, 0U>(
                plan, branches, variation, variation_component, output.x) &&
                apply_limited_branch_axis<false, 1U>(
                    plan, branches, variation, variation_component, output.y) &&
                apply_limited_branch_axis<false, 2U>(
                    plan, branches, variation, variation_component, output.z);
  return finite ? Status{}
                : Status{StatusCode::numerical_failure, kTransportNumerical};
}

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

Status cartesian_target_convection(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView flux, const KernelInvocation& invocation) noexcept {
  if (flux.certificate.valid()) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  return cartesian_provisional_convection(plan, scheme, flux, invocation);
}

Status diagnose_cartesian_convection_point(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView flux, ConstFieldView transported,
    std::uint8_t component, Int3 cell,
    ConvectionPointDiagnostic& diagnostic) noexcept {
  const std::uint8_t required_ghost_width =
      scheme == ConvectionScheme::central2 ? 1U : 2U;
  const Int3 cells = plan.cells();
  if (plan.fingerprint() == 0U ||
      static_cast<std::uint8_t>(scheme) >
          static_cast<std::uint8_t>(ConvectionScheme::tvd2) ||
      component >= transported.components || cell.x < 0 || cell.y < 0 ||
      cell.z < 0 || cell.x >= cells.x || cell.y >= cells.y ||
      cell.z >= cells.z ||
      !detail::valid_cell_view(transported, cells, component, 1U,
                               required_ghost_width) ||
      !detail::valid_flux_view(flux, cells, flux.revision) ||
      !flux.certificate.matches(flux)) {
    return {StatusCode::invalid_plan, kTransportKernel};
  }
  const bool uniform = plan.geometry_kind() == GeometryKind::uniform;
  switch (scheme) {
    case ConvectionScheme::central2:
      return uniform
                 ? diagnose_convection_geometry<
                       true, ConvectionScheme::central2>(
                       plan, flux, transported, component, cell, diagnostic)
                 : diagnose_convection_geometry<
                       false, ConvectionScheme::central2>(
                       plan, flux, transported, component, cell, diagnostic);
    case ConvectionScheme::limited_central2:
      return uniform
                 ? diagnose_convection_geometry<
                       true, ConvectionScheme::limited_central2>(
                       plan, flux, transported, component, cell, diagnostic)
                 : diagnose_convection_geometry<
                       false, ConvectionScheme::limited_central2>(
                       plan, flux, transported, component, cell, diagnostic);
    case ConvectionScheme::tvd2:
      return uniform
                 ? diagnose_convection_geometry<true,
                                                  ConvectionScheme::tvd2>(
                       plan, flux, transported, component, cell, diagnostic)
                 : diagnose_convection_geometry<false,
                                                  ConvectionScheme::tvd2>(
                       plan, flux, transported, component, cell, diagnostic);
  }
  return {StatusCode::invalid_plan, kTransportKernel};
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
