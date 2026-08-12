// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_exchange_plan.hpp"

#include "hundun/rt_error.hpp"
#include "rt_halo_detail.hpp"

#include <limits>
#include <utility>
#include <vector>

namespace hundun::runtime {
namespace {

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

int checked_upper_endpoint(int extent, int ghost_width) {
  if (ghost_width > std::numeric_limits<int>::max() - extent) {
    throw Error("halo plan upper ghost endpoint exceeds the int range");
  }
  return extent + ghost_width;
}

struct AxisBoxes {
  int send_begin{};
  int send_end{};
  int receive_begin{};
  int receive_end{};
};

AxisBoxes make_axis_boxes(int offset, int extent, int ghost_width,
                          int upper_endpoint) {
  if (offset == -1) {
    return AxisBoxes{0, ghost_width, -ghost_width, 0};
  }
  if (offset == 0) {
    return AxisBoxes{0, extent, 0, extent};
  }
  return AxisBoxes{extent - ghost_width, extent, extent, upper_endpoint};
}

}  // namespace

ExchangePlan ExchangePlan::create(
    const StructuredDecomposition& decomposition, Int3 local_extent,
    int ghost_width) {
  if (local_extent.x <= 0 || local_extent.y <= 0 || local_extent.z <= 0) {
    throw Error("halo plan local extents must be positive");
  }
  if (!same(local_extent, decomposition.local_extent())) {
    throw Error("halo plan local extent does not match the decomposition");
  }
  if (ghost_width < 0) {
    throw Error("halo plan ghost width must not be negative");
  }
  if (local_extent.x < ghost_width || local_extent.y < ghost_width ||
      local_extent.z < ghost_width) {
    throw Error("halo plan ghost width exceeds a local dimension");
  }

  const Int3 upper{
      checked_upper_endpoint(local_extent.x, ghost_width),
      checked_upper_endpoint(local_extent.y, ghost_width),
      checked_upper_endpoint(local_extent.z, ghost_width)};

  std::vector<ExchangeRegion> regions;
  regions.reserve(26U);
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        const Int3 offset{x, y, z};
        if (x == 0 && y == 0 && z == 0) {
          continue;
        }
        static_cast<void>(detail::halo_offset_code(offset));
        const AxisBoxes x_boxes =
            make_axis_boxes(x, local_extent.x, ghost_width, upper.x);
        const AxisBoxes y_boxes =
            make_axis_boxes(y, local_extent.y, ghost_width, upper.y);
        const AxisBoxes z_boxes =
            make_axis_boxes(z, local_extent.z, ghost_width, upper.z);
        regions.push_back(ExchangeRegion{
            offset,
            decomposition.neighbor_rank(offset),
            Box3{Int3{x_boxes.send_begin, y_boxes.send_begin,
                      z_boxes.send_begin},
                 Int3{x_boxes.send_end, y_boxes.send_end, z_boxes.send_end}},
            Box3{Int3{x_boxes.receive_begin, y_boxes.receive_begin,
                      z_boxes.receive_begin},
                 Int3{x_boxes.receive_end, y_boxes.receive_end,
                      z_boxes.receive_end}}});
      }
    }
  }
  return ExchangePlan(local_extent, ghost_width, std::move(regions));
}

ExchangePlan::ExchangePlan(Int3 local_extent, int ghost_width,
                           std::vector<ExchangeRegion> regions) noexcept
    : local_extent_(local_extent),
      ghost_width_(ghost_width),
      regions_(std::move(regions)) {}

const std::vector<ExchangeRegion>& ExchangePlan::regions() const noexcept {
  return regions_;
}

int ExchangePlan::ghost_width() const noexcept { return ghost_width_; }

}  // namespace hundun::runtime
