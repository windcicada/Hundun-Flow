// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_case.hpp"

#include <string_view>

namespace hundun::v04::detail {

// Validates only the typed thermophysical model. The source path belongs to
// CaseCompiler's file-identity contract and is intentionally not required here.
bool valid_thermophysical_spec(const ThermophysicalSpec& spec) noexcept;

// Replaces only the high-interval NASA enthalpy integration constant with the
// unique value implied by the accepted low interval at T_switch.  The
// operation is cold, deterministic, and idempotent; all other coefficients
// remain input authorities.
Status canonicalize_thermophysical_spec(
    ThermophysicalSpec& spec) noexcept;

// Hashes the complete typed payload, including its direct-root source path.
// Invalid models or empty source paths return zero.
PlanFingerprint thermophysical_spec_fingerprint(
    const ThermophysicalSpec& spec) noexcept;

Status parse_thermophysical_text(std::string_view text,
                                 ThermophysicalSpec& out) noexcept;

}  // namespace hundun::v04::detail
