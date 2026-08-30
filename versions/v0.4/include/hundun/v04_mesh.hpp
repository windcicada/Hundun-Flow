// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_case.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

namespace hundun::v04 {

enum class CartesianAxis : std::uint8_t { x, y, z };
enum class RegionFlag : std::uint8_t { solid = 0U, fluid = 1U };

struct GeometryBudget {
  std::uint64_t fixed_bytes_per_rank{};
  std::uint64_t bytes_per_owned_cell_upper_bound{};
};

class AxisMetrics {
 public:
  AxisMetrics() noexcept = default;
  AxisMetrics(const AxisMetrics&) = delete;
  AxisMetrics& operator=(const AxisMetrics&) = delete;
  AxisMetrics(AxisMetrics&&) noexcept = default;
  AxisMetrics& operator=(AxisMetrics&&) noexcept = default;

  Span<const double> faces() const noexcept;
  Span<const double> centres() const noexcept;
  Span<const double> widths() const noexcept;
  Span<const double> inverse_widths() const noexcept;
  Span<const double> inverse_centre_distances() const noexcept;
  bool uniform() const noexcept { return uniform_; }
  double uniform_width() const noexcept { return uniform_width_; }
  double uniform_inverse_width() const noexcept {
    return uniform_inverse_width_;
  }

 private:
  friend class CartesianGeometryCompiler;
  std::vector<double> faces_;
  std::vector<double> centres_;
  std::vector<double> widths_;
  std::vector<double> inverse_widths_;
  std::vector<double> inverse_centre_distances_;
  bool uniform_{};
  double uniform_width_{};
  double uniform_inverse_width_{};
};

class CartesianGeometryPlan {
 public:
  CartesianGeometryPlan() noexcept = default;
  CartesianGeometryPlan(const CartesianGeometryPlan&) = delete;
  CartesianGeometryPlan& operator=(const CartesianGeometryPlan&) = delete;
  CartesianGeometryPlan(CartesianGeometryPlan&&) noexcept = default;
  CartesianGeometryPlan& operator=(CartesianGeometryPlan&&) noexcept = default;

  GeometryKind kind() const noexcept { return kind_; }
  Int3 global_cells() const noexcept { return global_cells_; }
  Real3 lower() const noexcept { return lower_; }
  Real3 upper() const noexcept { return upper_; }
  RevisionToken topology_revision() const noexcept {
    return topology_revision_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  const AxisMetrics& axis(CartesianAxis selected) const noexcept;
  const AxisMetrics& x() const noexcept { return axes_[0]; }
  const AxisMetrics& y() const noexcept { return axes_[1]; }
  const AxisMetrics& z() const noexcept { return axes_[2]; }

 private:
  friend class CartesianGeometryCompiler;
  GeometryKind kind_{GeometryKind::uniform};
  Int3 global_cells_{};
  Real3 lower_{};
  Real3 upper_{};
  AxisMetrics axes_[3];
  RevisionToken topology_revision_{};
  PlanFingerprint fingerprint_{};
};

struct MeshPatch {
  Int3 begin{};
  Int3 cells{};
  Int3 process_grid{};
  Int3 process_coord{};
};

struct CpuTile {
  Int3 begin{};
  Int3 cells{};
};

class CartesianGeometryCompiler {
 public:
  static Status compile(MPI_Comm communicator, const CartesianMeshSpec& mesh,
                        GeometryBudget budget, CartesianGeometryPlan& geometry,
                        MeshPatch& patch) noexcept;
};

Status build_cpu_tiles(const MeshPatch& patch, Int3 target_cells,
                       std::vector<CpuTile>& out) noexcept;

struct TriangleInput {
  Real3 a{};
  Real3 b{};
  Real3 c{};
};

class TriangleSoA {
 public:
  TriangleSoA() noexcept = default;
  TriangleSoA(const TriangleSoA&) = delete;
  TriangleSoA& operator=(const TriangleSoA&) = delete;
  TriangleSoA(TriangleSoA&&) noexcept = default;
  TriangleSoA& operator=(TriangleSoA&&) noexcept = default;

  std::size_t size() const noexcept { return ax_.size(); }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Span<const double> ax() const noexcept { return {ax_.data(), ax_.size()}; }
  Span<const double> ay() const noexcept { return {ay_.data(), ay_.size()}; }
  Span<const double> az() const noexcept { return {az_.data(), az_.size()}; }
  Span<const double> bx() const noexcept { return {bx_.data(), bx_.size()}; }
  Span<const double> by() const noexcept { return {by_.data(), by_.size()}; }
  Span<const double> bz() const noexcept { return {bz_.data(), bz_.size()}; }
  Span<const double> cx() const noexcept { return {cx_.data(), cx_.size()}; }
  Span<const double> cy() const noexcept { return {cy_.data(), cy_.size()}; }
  Span<const double> cz() const noexcept { return {cz_.data(), cz_.size()}; }
  Span<const double> e1x() const noexcept { return {e1x_.data(), e1x_.size()}; }
  Span<const double> e1y() const noexcept { return {e1y_.data(), e1y_.size()}; }
  Span<const double> e1z() const noexcept { return {e1z_.data(), e1z_.size()}; }
  Span<const double> e2x() const noexcept { return {e2x_.data(), e2x_.size()}; }
  Span<const double> e2y() const noexcept { return {e2y_.data(), e2y_.size()}; }
  Span<const double> e2z() const noexcept { return {e2z_.data(), e2z_.size()}; }
  Span<const double> nx() const noexcept { return {nx_.data(), nx_.size()}; }
  Span<const double> ny() const noexcept { return {ny_.data(), ny_.size()}; }
  Span<const double> nz() const noexcept { return {nz_.data(), nz_.size()}; }
  Span<const double> min_x() const noexcept {
    return {min_x_.data(), min_x_.size()};
  }
  Span<const double> min_y() const noexcept {
    return {min_y_.data(), min_y_.size()};
  }
  Span<const double> min_z() const noexcept {
    return {min_z_.data(), min_z_.size()};
  }
  Span<const double> max_x() const noexcept {
    return {max_x_.data(), max_x_.size()};
  }
  Span<const double> max_y() const noexcept {
    return {max_y_.data(), max_y_.size()};
  }
  Span<const double> max_z() const noexcept {
    return {max_z_.data(), max_z_.size()};
  }

 private:
  friend class StlScanCompiler;
  std::vector<double> ax_, ay_, az_;
  // Preserve the exact parsed STL vertices for cold topology construction.
  // Reconstructing b/c from a + edge can change a valid binary32 coordinate
  // after promotion to binary64 and destroy exact shared-edge identity.
  std::vector<double> bx_, by_, bz_;
  std::vector<double> cx_, cy_, cz_;
  std::vector<double> e1x_, e1y_, e1z_;
  std::vector<double> e2x_, e2y_, e2z_;
  std::vector<double> nx_, ny_, nz_;
  std::vector<double> min_x_, min_y_, min_z_;
  std::vector<double> max_x_, max_y_, max_z_;
  PlanFingerprint fingerprint_{};
};

class StlScanPlan {
 public:
  StlScanPlan() noexcept = default;
  StlScanPlan(const StlScanPlan&) = delete;
  StlScanPlan& operator=(const StlScanPlan&) = delete;
  StlScanPlan(StlScanPlan&&) noexcept = default;
  StlScanPlan& operator=(StlScanPlan&&) noexcept = default;

  CartesianAxis scan_axis() const noexcept { return scan_axis_; }
  std::size_t triangle_count() const noexcept { return triangles_.size(); }
  std::size_t scan_line_count() const noexcept {
    return line_offsets_.empty() ? 0U : line_offsets_.size() - 1U;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  const TriangleSoA& triangles() const noexcept { return triangles_; }

  Status classify(const CartesianGeometryPlan& geometry,
                  const MeshPatch& patch,
                  Span<std::uint8_t> region) const noexcept;

 private:
  friend class StlScanCompiler;
  TriangleSoA triangles_;
  CartesianAxis scan_axis_{CartesianAxis::y};
  Int3 patch_begin_{};
  Int3 patch_cells_{};
  std::vector<std::size_t> line_offsets_;
  std::vector<double> span_bounds_;
  PlanFingerprint geometry_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

struct StlScanBudget {
  std::uint64_t max_persistent_bytes_per_rank{};
  std::uint64_t max_peak_bytes_per_rank{};
  std::uint64_t max_bin_references{};
  std::uint64_t max_events_per_line{};
  std::uint32_t worker_threads{1U};
};

class StlScanCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const std::filesystem::path& case_root,
                        const std::optional<std::filesystem::path>& stl_file,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch, CartesianAxis scan_axis,
                        StlScanBudget budget,
                        StlScanPlan& out) noexcept;

  static Status compile_triangles(const CartesianGeometryPlan& geometry,
                                  const MeshPatch& patch,
                                  Span<const TriangleInput> triangles,
                                  CartesianAxis scan_axis,
                                  StlScanBudget budget,
                                  StlScanPlan& out) noexcept;
};

}  // namespace hundun::v04
