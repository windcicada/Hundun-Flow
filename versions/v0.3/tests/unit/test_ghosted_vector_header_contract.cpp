// SPDX-License-Identifier: Apache-2.0

#include "hundun/lin_ghosted_vector.hpp"
#include "hundun/lin_ghosted_vector_halo.hpp"

#include <type_traits>

static_assert(std::is_default_constructible_v<hundun::linear::VectorLayout>);
static_assert(std::is_copy_constructible_v<hundun::linear::VectorLayout>);
static_assert(!std::is_copy_constructible_v<hundun::linear::GhostedVector>);
static_assert(std::is_nothrow_move_constructible_v<
              hundun::linear::GhostedVector>);
static_assert(!std::is_copy_constructible_v<
              hundun::linear::GhostedVectorHalo>);
static_assert(std::is_nothrow_move_constructible_v<
              hundun::linear::GhostedVectorHalo>);
static_assert(!std::is_move_assignable_v<hundun::linear::GhostedVectorHalo>);

int main() { return 0; }
