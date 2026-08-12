// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_quadratic_reconstruction.hpp"

#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_view.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::QuadraticReconstruction;
using hundun::immersed::ReconstructionQuality;
using hundun::immersed::WeightedDonor;
using hundun::mesh::GlobalCellId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::FieldView;
using hundun::runtime::Int3;
using hundun::runtime::Real3;

using Create = QuadraticReconstruction (*)(Real3, Real3, Real3, Real3, double,
                                           Int3, const std::vector<Int3> &,
                                           const MeshTopology &,
                                           const MeshGeometry &);
using Value = double (QuadraticReconstruction::*)(
    Real3, const FieldView<const double> &, std::size_t) const;
using Gradient = Real3 (QuadraticReconstruction::*)(
    Real3, const FieldView<const double> &, std::size_t) const;
using ConstrainedValue = double (QuadraticReconstruction::*)(
    Real3, const FieldView<const double> &, std::size_t, double) const;
using ConstrainedGradient = Real3 (QuadraticReconstruction::*)(
    Real3, const FieldView<const double> &, std::size_t, double) const;
template <class T, class = void>
struct HasOriginConstraint : std::false_type {};
template <class T>
struct HasOriginConstraint<
    T, std::void_t<decltype(&T::value_with_origin_constraint),
                   decltype(&T::gradient_with_origin_constraint),
                   decltype(&T::value_with_origin_normal_gradient),
                   decltype(&T::gradient_with_origin_normal_gradient)>>
    : std::true_type {};

static_assert(std::is_final_v<WeightedDonor>);
static_assert(
    std::is_same_v<decltype(WeightedDonor::global_cell), GlobalCellId>);
static_assert(std::is_same_v<decltype(WeightedDonor::weight), double>);
static_assert(std::is_final_v<ReconstructionQuality>);
static_assert(
    std::is_same_v<decltype(ReconstructionQuality::rank), std::uint32_t>);
static_assert(std::is_same_v<
              decltype(ReconstructionQuality::condition_estimate), double>);
static_assert(
    std::is_same_v<decltype(ReconstructionQuality::halo_reach), std::uint32_t>);
static_assert(std::is_same_v<decltype(ReconstructionQuality::pivot_fingerprint),
                             std::uint64_t>);
static_assert(
    std::is_same_v<Create, decltype(&QuadraticReconstruction::create)>);
static_assert(std::is_same_v<Value, decltype(&QuadraticReconstruction::value)>);
static_assert(
    std::is_same_v<Gradient, decltype(&QuadraticReconstruction::gradient)>);
static_assert(!HasOriginConstraint<QuadraticReconstruction>::value);
static_assert(
    std::is_same_v<
        decltype(std::declval<const QuadraticReconstruction &>().quality()),
        const ReconstructionQuality &>);
static_assert(
    noexcept(std::declval<const QuadraticReconstruction &>().quality()));
static_assert(std::is_copy_constructible_v<QuadraticReconstruction>);
static_assert(std::is_copy_assignable_v<QuadraticReconstruction>);
static_assert(std::is_nothrow_move_constructible_v<QuadraticReconstruction>);
static_assert(std::is_nothrow_move_assignable_v<QuadraticReconstruction>);
static_assert(!std::is_default_constructible_v<QuadraticReconstruction>);

} // namespace

int main() { return 0; }
