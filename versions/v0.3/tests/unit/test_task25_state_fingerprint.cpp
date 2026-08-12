// SPDX-License-Identifier: Apache-2.0

#include "src/app_performance_fingerprint_detail.hpp"
#include "src/app_performance_crc64_detail.hpp"
#include "tests/support/test_main.hpp"

#include <optional>

int main() {
  return hundun::test::run([] {
    HUNDUN_CHECK(
        hundun::application::detail::performance_crc64_ecma("123456789") ==
        UINT64_C(0x6c40df5f0b497347));
    hundun::flow::FlowLayerValues original;
    original.density = {1.0, 2.0};
    original.velocity = {3.0, 4.0, 5.0};
    original.mechanical_pressure = {6.0};
    original.face_velocity = {7.0, 8.0, 9.0};
    original.face_mass_flux = {10.0};
    original.transported_cell_fields = {{11.0, 12.0}, {13.0}};

    const auto expected =
        hundun::application::detail::performance_flow_layer_encoding(original);
    auto copy = original;
    HUNDUN_CHECK(
        hundun::application::detail::performance_flow_layer_encoding(copy) ==
        expected);

    copy.velocity[1] = 4.5;
    HUNDUN_CHECK(
        hundun::application::detail::performance_flow_layer_encoding(copy) !=
        expected);

    copy = original;
    copy.transported_cell_fields[1][0] = 13.5;
    HUNDUN_CHECK(
        hundun::application::detail::performance_flow_layer_encoding(copy) !=
        expected);

    copy = original;
    copy.transported_cell_fields.push_back({});
    HUNDUN_CHECK(
        hundun::application::detail::performance_flow_layer_encoding(copy) !=
        expected);

    copy = original;
    copy.transported_cell_fields[0].pop_back();
    HUNDUN_CHECK(
        hundun::application::detail::performance_flow_layer_encoding(copy) !=
        expected);

    const std::optional<hundun::flow::IdealGasClosureState> no_closure;
    const auto none_encoding =
        hundun::application::detail::performance_closure_state_encoding(
            no_closure);
    hundun::flow::IdealGasClosureState closure;
    closure.mode = hundun::flow::IdealGasPressureMode::closed_dynamic;
    closure.thermodynamic_pressure_pa = 101325.0;
    closure.target_mass_kg = 2.5;
    closure.revision = 7U;
    const auto closure_encoding =
        hundun::application::detail::performance_closure_state_encoding(
            closure);
    HUNDUN_CHECK(closure_encoding != none_encoding);
    const auto closure_digest =
        hundun::application::detail::tagged_performance_crc64(
            "hundun-performance-state-fp-v1", closure_encoding);

    auto mutated_closure = closure;
    mutated_closure.thermodynamic_pressure_pa = 101326.0;
    HUNDUN_CHECK(
        hundun::application::detail::tagged_performance_crc64(
            "hundun-performance-state-fp-v1",
            hundun::application::detail::performance_closure_state_encoding(
                mutated_closure)) != closure_digest);
    mutated_closure = closure;
    ++mutated_closure.revision;
    HUNDUN_CHECK(
        hundun::application::detail::tagged_performance_crc64(
            "hundun-performance-state-fp-v1",
            hundun::application::detail::performance_closure_state_encoding(
                mutated_closure)) != closure_digest);
  });
}
