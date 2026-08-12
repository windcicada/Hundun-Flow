// SPDX-License-Identifier: Apache-2.0
#pragma once

#if !defined(HUNDUN_IMMERSED_ENABLE_TEST_ACCESS)
#error "immersed_test_access.hpp is available only to approved test targets"
#endif

#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"

#include "ib_surface_bvh_detail.hpp"

#include <cstddef>
#include <filesystem>
#include <utility>
#include <vector>

namespace hundun::immersed::detail {

bool coplanar_contact_is_forbidden_for_test(
    const SurfaceTriangle &target, runtime::Real3 first, runtime::Real3 second,
    bool first_is_canonical_shared_vertex,
    bool second_is_canonical_shared_vertex) noexcept;

} // namespace hundun::immersed::detail

namespace hundun::immersed::test {

class ImmersedTestAccess final {
public:
  static ImmersedSurface
  load_collective_with_chunk_limit(const std::filesystem::path &path,
                                   double length_scale_to_m,
                                   const runtime::MpiContext &mpi, int root,
                                   std::size_t maximum_chunk_bytes) {
    return ImmersedSurface::load_collective_impl(path, length_scale_to_m, mpi,
                                                 root, maximum_chunk_bytes);
  }

  static std::vector<TriangleId> bounded_candidates(const SurfaceQuery &query,
                                                    runtime::Real3 minimum_m,
                                                    runtime::Real3 maximum_m) {
    return detail::bounded_candidates(*query.storage_, minimum_m, maximum_m);
  }

  static std::vector<TriangleId> bvh_order(const SurfaceQuery &query) {
    return query.storage_->triangle_order;
  }

  static SurfaceQuery
  create_query_with_initial_order(const ImmersedSurface &surface,
                                  std::vector<TriangleId> initial_order) {
    return SurfaceQuery(detail::build_surface_query_with_initial_order(
        surface.storage_, std::move(initial_order)));
  }

  static bool resolve_parity(std::array<bool, 3> votes) {
    return detail::require_consistent_parity(votes);
  }

  static bool coplanar_contact_is_forbidden(
      const SurfaceTriangle &target, runtime::Real3 first,
      runtime::Real3 second, bool first_is_canonical_shared_vertex,
      bool second_is_canonical_shared_vertex) noexcept {
    return detail::coplanar_contact_is_forbidden_for_test(
        target, first, second, first_is_canonical_shared_vertex,
        second_is_canonical_shared_vertex);
  }
};

} // namespace hundun::immersed::test
