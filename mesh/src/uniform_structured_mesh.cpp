// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/uniform_structured_mesh.hpp"

#include "geometry_arithmetic.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::mesh {
namespace {

bool same(runtime::Int3 lhs, runtime::Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool positive(runtime::Int3 value) noexcept {
  return value.x > 0 && value.y > 0 && value.z > 0;
}

bool positive(runtime::Real3 value) noexcept {
  return value.x > 0.0 && value.y > 0.0 && value.z > 0.0;
}

static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::radix == 2);
static_assert(std::numeric_limits<double>::digits == 53);
static_assert(std::numeric_limits<double>::min_exponent == -1021);
static_assert(std::numeric_limits<double>::max_exponent == 1024);
static_assert(std::numeric_limits<double>::has_denorm ==
              std::denorm_present);
static_assert(std::numeric_limits<double>::has_infinity);
static_assert(std::numeric_limits<std::uint32_t>::digits == 32);
static_assert(std::numeric_limits<std::uint64_t>::digits == 64);

constexpr int kSignificandBits = std::numeric_limits<double>::digits;
constexpr int kMinimumNormalExponent =
    std::numeric_limits<double>::min_exponent - 1;
constexpr int kMaximumFiniteExponent =
    std::numeric_limits<double>::max_exponent - 1;
constexpr int kSubnormalQuantumExponent =
    kMinimumNormalExponent - (kSignificandBits - 1);
constexpr int kLimbBits = std::numeric_limits<std::uint32_t>::digits;

using SignificandLimbs = std::array<std::uint32_t, 2>;
using ProductLimbs = std::array<std::uint32_t, 5>;

struct BinaryFactor final {
  std::uint64_t significand;
  int exponent;
};

BinaryFactor decompose(double value) noexcept {
  int exponent = 0;
  const double fraction = std::frexp(value, &exponent);
  const auto significand = static_cast<std::uint64_t>(
      std::ldexp(fraction, kSignificandBits));
  return BinaryFactor{significand, exponent - kSignificandBits};
}

SignificandLimbs split(std::uint64_t value) noexcept {
  return SignificandLimbs{static_cast<std::uint32_t>(value),
                          static_cast<std::uint32_t>(value >> kLimbBits)};
}

template <std::size_t LeftSize, std::size_t RightSize>
std::array<std::uint32_t, LeftSize + RightSize> multiply_limbs(
    const std::array<std::uint32_t, LeftSize>& lhs,
    const std::array<std::uint32_t, RightSize>& rhs) noexcept {
  std::array<std::uint32_t, LeftSize + RightSize> result{};
  for (std::size_t left = 0; left < LeftSize; ++left) {
    std::uint64_t carry = 0;
    for (std::size_t right = 0; right < RightSize; ++right) {
      const std::uint64_t accumulator =
          static_cast<std::uint64_t>(lhs[left]) * rhs[right] +
          result[left + right] + carry;
      result[left + right] = static_cast<std::uint32_t>(accumulator);
      carry = accumulator >> kLimbBits;
    }
    result[left + RightSize] = static_cast<std::uint32_t>(carry);
  }
  return result;
}

ProductLimbs multiply_significands(std::uint64_t x, std::uint64_t y,
                                   std::uint64_t z) noexcept {
  const auto xy = multiply_limbs(split(x), split(y));
  const auto xyz = multiply_limbs(xy, split(z));
  // Three 53-bit integers occupy at most 159 bits, so xyz[5] is zero.
  return ProductLimbs{xyz[0], xyz[1], xyz[2], xyz[3], xyz[4]};
}

bool product_bit(const ProductLimbs& product, int position) noexcept {
  if (position < 0 ||
      position >= static_cast<int>(product.size()) * kLimbBits) {
    return false;
  }
  return ((product[static_cast<std::size_t>(position / kLimbBits)] >>
           (position % kLimbBits)) &
          std::uint32_t{1}) != 0;
}

int product_bit_length(const ProductLimbs& product) noexcept {
  for (std::size_t index = product.size(); index > 0; --index) {
    std::uint32_t word = product[index - 1];
    if (word == 0) {
      continue;
    }
    int word_bits = 0;
    while (word != 0) {
      ++word_bits;
      word >>= 1;
    }
    return static_cast<int>((index - 1) * kLimbBits) + word_bits;
  }
  return 0;
}

bool any_product_bit_below(const ProductLimbs& product,
                           int exclusive_end) noexcept {
  if (exclusive_end <= 0) {
    return false;
  }
  const int complete_limbs = exclusive_end / kLimbBits;
  for (int index = 0;
       index < complete_limbs &&
       index < static_cast<int>(product.size());
       ++index) {
    if (product[static_cast<std::size_t>(index)] != 0) {
      return true;
    }
  }
  const int remaining_bits = exclusive_end % kLimbBits;
  if (remaining_bits == 0 ||
      complete_limbs >= static_cast<int>(product.size())) {
    return false;
  }
  const std::uint32_t mask =
      (std::uint32_t{1} << remaining_bits) - std::uint32_t{1};
  return (product[static_cast<std::size_t>(complete_limbs)] & mask) != 0;
}

std::uint64_t round_right_shift(const ProductLimbs& product,
                                int shift) noexcept {
  std::uint64_t quotient = 0;
  for (int position = product_bit_length(product) - 1; position >= shift;
       --position) {
    quotient = (quotient << 1) |
               static_cast<std::uint64_t>(product_bit(product, position));
  }
  const bool guard = product_bit(product, shift - 1);
  const bool sticky = any_product_bit_below(product, shift - 1);
  if (guard && (sticky || (quotient & std::uint64_t{1}) != 0)) {
    ++quotient;
  }
  return quotient;
}

double range_safe_product_impl(runtime::Real3 factors) noexcept {
  const BinaryFactor x = decompose(factors.x);
  const BinaryFactor y = decompose(factors.y);
  const BinaryFactor z = decompose(factors.z);
  const ProductLimbs product =
      multiply_significands(x.significand, y.significand, z.significand);
  const int product_exponent = x.exponent + y.exponent + z.exponent;
  const int leading_exponent =
      product_exponent + product_bit_length(product) - 1;
  if (leading_exponent > kMaximumFiniteExponent) {
    return std::numeric_limits<double>::infinity();
  }

  const int quantum_exponent =
      leading_exponent >= kMinimumNormalExponent
          ? leading_exponent - (kSignificandBits - 1)
          : kSubnormalQuantumExponent;
  const std::uint64_t rounded =
      round_right_shift(product, quantum_exponent - product_exponent);
  if (rounded == 0) {
    return 0.0;
  }
  if (leading_exponent == kMaximumFiniteExponent &&
      rounded == (std::uint64_t{1} << kSignificandBits)) {
    return std::numeric_limits<double>::infinity();
  }
  return std::ldexp(static_cast<double>(rounded), quantum_exponent);
}

}  // namespace

namespace detail {

double range_safe_product(runtime::Real3 factors) noexcept {
  return range_safe_product_impl(factors);
}

}  // namespace detail

UniformStructuredMesh::UniformStructuredMesh(
    runtime::Int3 global_cells, runtime::Real3 origin_m,
    runtime::Real3 length_m,
    const runtime::StructuredDecomposition& decomposition) {
  if (!positive(global_cells)) {
    throw runtime::Error("mesh global cell counts must be positive");
  }
  if (!same(global_cells, decomposition.global_extent())) {
    throw runtime::Error(
        "mesh global cell counts must match the decomposition extent");
  }
  if (!finite(origin_m)) {
    throw runtime::Error("mesh origin components must be finite");
  }
  if (!finite(length_m) || !positive(length_m)) {
    throw runtime::Error("mesh lengths must be finite and positive");
  }

  const runtime::Real3 spacing{
      length_m.x / static_cast<double>(global_cells.x),
      length_m.y / static_cast<double>(global_cells.y),
      length_m.z / static_cast<double>(global_cells.z)};
  if (!finite(spacing) || !positive(spacing)) {
    throw runtime::Error("mesh spacing must be finite and positive");
  }

  const double cell_volume = detail::range_safe_product(spacing);
  if (!std::isfinite(cell_volume) || cell_volume <= 0.0) {
    throw runtime::Error("mesh cell volume must be finite and positive");
  }

  const runtime::Real3 upper_endpoint{origin_m.x + length_m.x,
                                      origin_m.y + length_m.y,
                                      origin_m.z + length_m.z};
  if (!finite(upper_endpoint)) {
    throw runtime::Error("mesh upper endpoints must be finite");
  }

  origin_m_ = origin_m;
  spacing_m_ = spacing;
  cell_volume_m3_ = cell_volume;
  local_extent_ = decomposition.local_extent();
  owned_global_box_ = decomposition.owned_box();
}

runtime::Real3 UniformStructuredMesh::spacing_m() const noexcept {
  return spacing_m_;
}

runtime::Real3 UniformStructuredMesh::cell_center(
    runtime::Int3 local_cell) const {
  if (local_cell.x < 0 || local_cell.x >= local_extent_.x ||
      local_cell.y < 0 || local_cell.y >= local_extent_.y ||
      local_cell.z < 0 || local_cell.z >= local_extent_.z) {
    throw runtime::Error("local cell is outside the mesh owned box");
  }

  const runtime::Int3 global_cell{owned_global_box_.begin.x + local_cell.x,
                                  owned_global_box_.begin.y + local_cell.y,
                                  owned_global_box_.begin.z + local_cell.z};
  return runtime::Real3{
      origin_m_.x + (static_cast<double>(global_cell.x) + 0.5) * spacing_m_.x,
      origin_m_.y + (static_cast<double>(global_cell.y) + 0.5) * spacing_m_.y,
      origin_m_.z + (static_cast<double>(global_cell.z) + 0.5) * spacing_m_.z};
}

double UniformStructuredMesh::cell_volume_m3() const noexcept {
  return cell_volume_m3_;
}

runtime::Int3 UniformStructuredMesh::local_extent() const noexcept {
  return local_extent_;
}

runtime::Box3 UniformStructuredMesh::owned_global_box() const noexcept {
  return owned_global_box_;
}

}  // namespace hundun::mesh
