// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/flow_checkpoint_v3.hpp"

#include "tests/support/test_main.hpp"

#include <cstdint>
#include <type_traits>

static_assert(std::is_enum_v<hundun::flow::CheckpointV3Operation>);
static_assert(std::is_enum_v<hundun::flow::CheckpointV3Disposition>);
static_assert(std::is_enum_v<hundun::flow::CheckpointV3FailureReason>);
static_assert(std::is_enum_v<hundun::flow::CheckpointV3Phase>);
static_assert(std::is_enum_v<hundun::flow::CheckpointV3Presence>);
static_assert(hundun::flow::CheckpointV3Phase::rank_payload !=
              hundun::flow::CheckpointV3Phase::manifest);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::CheckpointV3Presence::constant_static_ibm) ==
              1U);
static_assert(
    static_cast<std::uint8_t>(
        hundun::flow::CheckpointV3Presence::constant_body_fitted_wale) == 2U);
static_assert(
    static_cast<std::uint8_t>(
        hundun::flow::CheckpointV3Presence::constant_static_ibm_wale) == 3U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::CheckpointV3Presence::material_static_ibm) ==
              4U);
static_assert(
    static_cast<std::uint8_t>(
        hundun::flow::CheckpointV3Presence::material_body_fitted_wale) == 5U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::CheckpointV3Presence::material_static_ibm_wale) ==
              6U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::CheckpointV3Presence::ideal_gas_static_ibm) ==
              7U);
static_assert(
    static_cast<std::uint8_t>(
        hundun::flow::CheckpointV3Presence::ideal_gas_body_fitted_wale) == 8U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::CheckpointV3Presence::ideal_gas_static_ibm_wale) ==
              9U);

static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3ControlState::last_retry_count),
              std::uint32_t>);
static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3ControlState::proposed_next_dt_s),
              double>);

static_assert(std::is_same_v<
              decltype(&hundun::flow::CheckpointV3Report::operation),
              hundun::flow::CheckpointV3Operation (
                  hundun::flow::CheckpointV3Report::*)() const noexcept>);
static_assert(std::is_same_v<
              decltype(&hundun::flow::CheckpointV3Report::disposition),
              hundun::flow::CheckpointV3Disposition (
                  hundun::flow::CheckpointV3Report::*)() const noexcept>);
static_assert(std::is_same_v<
              decltype(&hundun::flow::CheckpointV3Report::presence),
              hundun::flow::CheckpointV3Presence (
                  hundun::flow::CheckpointV3Report::*)() const noexcept>);

static_assert(std::is_same_v<
              decltype(&hundun::flow::CheckpointV3ReadResult::restored),
              bool (hundun::flow::CheckpointV3ReadResult::*)() const noexcept>);
static_assert(std::is_same_v<
              decltype(&hundun::flow::CheckpointV3ReadResult::control_state),
              const hundun::flow::CheckpointV3ControlState &(
                  hundun::flow::CheckpointV3ReadResult::*)() const>);
static_assert(std::is_same_v<
              decltype(&hundun::flow::checkpoint_v3_performance_counters),
              hundun::diagnostics::Stage3PerformanceCounters (*)() noexcept>);

static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3WriteModules::wale),
              const hundun::les::WaleModel *>);
static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3WriteModules::flow),
              const hundun::flow::FixedStepImmersedFlow *>);
static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3WriteModules::ideal_gas),
              const hundun::flow::IdealGasClosure *>);
static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3ReadModules::wale),
              const hundun::les::WaleModel *>);
static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3ReadModules::flow),
              hundun::flow::FixedStepImmersedFlow *>);
static_assert(std::is_same_v<
              decltype(hundun::flow::CheckpointV3ReadModules::ideal_gas),
              hundun::flow::IdealGasClosure *>);

using WriteModulesOverload = hundun::flow::CheckpointV3Report (*)(
    const hundun::runtime::MpiContext &,
    const hundun::runtime::StructuredDecomposition &,
    const hundun::mesh::MeshTopology &, const hundun::mesh::MeshGeometry &,
    const hundun::boundary::BoundaryRegistry &,
    const hundun::config::ImmersedFlowCaseConfig &,
    const hundun::flow::CheckpointV3WriteModules &,
    const hundun::flow::FlowState &, hundun::flow::CheckpointV3ControlState,
    const std::filesystem::path &);
using ReadModulesOverload = hundun::flow::CheckpointV3ReadResult (*)(
    const hundun::runtime::MpiContext &,
    const hundun::runtime::StructuredDecomposition &,
    const hundun::mesh::MeshTopology &, const hundun::mesh::MeshGeometry &,
    const hundun::boundary::BoundaryRegistry &,
    const hundun::config::ImmersedFlowCaseConfig &,
    const hundun::flow::CheckpointV3ReadModules &, hundun::flow::FlowState &,
    const std::filesystem::path &);

static_assert(std::is_same_v<
              decltype(static_cast<WriteModulesOverload>(
                  &hundun::flow::write_checkpoint_v3)),
              WriteModulesOverload>);
static_assert(std::is_same_v<
              decltype(static_cast<ReadModulesOverload>(
                  &hundun::flow::read_checkpoint_v3)),
              ReadModulesOverload>);

int main() { return hundun::test::run([] {}); }
