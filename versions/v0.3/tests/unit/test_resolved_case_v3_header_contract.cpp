// SPDX-License-Identifier: Apache-2.0

#include "src/app_resolved_case_v3_broadcast_detail.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/cfg_resolved_case_v3_loader.hpp"

#include <cstdint>
#include <filesystem>
#include <mpi.h>
#include <type_traits>
#include <utility>
#include <variant>

static_assert(std::variant_size_v<hundun::config::ResolvedCase> == 2);
static_assert(std::variant_size_v<hundun::config::ResolvedCaseV3> == 3);
static_assert(std::is_same_v<
              std::underlying_type_t<hundun::config::ImmersedBoundaryModel>,
              std::uint8_t>);
static_assert(
    std::is_same_v<std::underlying_type_t<hundun::config::ImmersedFluidSide>,
                   std::uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<hundun::config::LesModel>,
                             std::uint8_t>);
static_assert(
    std::is_same_v<decltype(hundun::config::load_resolved_case_v3(
                       std::declval<const std::filesystem::path &>())),
                   hundun::config::ResolvedCaseV3>);
static_assert(
    std::is_same_v<decltype(hundun::config::to_resolved_json_v3(
                       std::declval<const hundun::config::ResolvedCaseV3 &>())),
                   std::string>);
static_assert(
    std::is_same_v<decltype(hundun::config::broadcast_resolved_case_v3(
                       std::declval<MPI_Comm>(), 0,
                       std::declval<const hundun::config::ResolvedCaseV3 *>())),
                   hundun::config::ResolvedCaseV3>);

int main() { return 0; }
