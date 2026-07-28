#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#
# HUNDUN_STAGE2_LIST_ONLY=1 validates and prints the exact registered
# inventory without executing numerical tests. It exists only for the
# Task 26 contract test; the normal no-argument path always runs the matrix.

set -euo pipefail

fail() {
  printf 'stage2 acceptance: %s\n' "$*" >&2
  exit 2
}

if test "$#" -ne 0; then
  fail "no command-line arguments are accepted"
fi

source_root=$(
  CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P
)
requested_build=${HUNDUN_STAGE2_BUILD_DIR:-"${source_root}/build/debug"}
if ! test -d "${requested_build}"; then
  fail "build directory not found: ${requested_build}"
fi
build_root=$(CDPATH= cd -- "${requested_build}" && pwd -P)

cache="${build_root}/CMakeCache.txt"
if ! test -f "${cache}"; then
  fail "CMake cache not found: ${build_root}/CMakeCache.txt"
fi

cache_source=$(sed -n 's/^HundunFlow_SOURCE_DIR:STATIC=//p' "${cache}")
if test "${cache_source}" != "${source_root}"; then
  fail "build source mismatch: expected ${source_root}, found ${cache_source:-<missing>}"
fi
cache_type=$(sed -n 's/^CMAKE_BUILD_TYPE:STRING=//p' "${cache}")
if test "${cache_type}" != "Debug"; then
  fail "build type mismatch: expected Debug, found ${cache_type:-<missing>}"
fi
cache_tests=$(sed -n 's/^HUNDUN_BUILD_TESTS:BOOL=//p' "${cache}")
if test "${cache_tests}" != "ON"; then
  fail "build test registry is disabled"
fi

build_identity="${build_root}/generated/applications/hundun/stage2_performance_build.hpp"
if ! test -f "${build_identity}"; then
  fail "build identity not found: ${build_identity}"
fi
build_commit=$(
  awk '/performance_source_commit[ \t]*=/{getline; gsub(/[";]/, ""); gsub(/^[ \t]+|[ \t]+$/, ""); print; exit}' \
    "${build_identity}"
)
if ! source_commit=$(git -C "${source_root}" rev-parse HEAD 2>/dev/null); then
  fail "source HEAD is unavailable"
fi
if test "${build_commit}" != "${source_commit}"; then
  fail "build HEAD mismatch: expected ${source_commit}, found ${build_commit:-<missing>}"
fi
if ! test -f "${source_root}/docs/numerics/stage2-capability-ledger.md"; then
  fail "capability ledger not found"
fi

ctest_command=${CTEST_COMMAND:-ctest}
case "${ctest_command}" in
  */*)
    test -x "${ctest_command}" ||
      fail "CTest command is not executable: ${ctest_command}"
    ;;
  *)
    command -v "${ctest_command}" >/dev/null 2>&1 ||
      fail "CTest command is unavailable: ${ctest_command}"
    ;;
esac

required_tests=(
  test_stage2_task2_dispatch
  test_mesh_geometry_1_rank
  test_mesh_geometry_2_rank
  test_mesh_geometry_4_rank
  test_conjugate_gradient_1_rank
  test_conjugate_gradient_2_rank
  test_conjugate_gradient_4_rank
  test_bicgstab_1_rank
  test_bicgstab_2_rank
  test_bicgstab_4_rank
  test_matrix_free_poisson_1_rank
  test_matrix_free_poisson_2_rank
  test_matrix_free_poisson_4_rank
  test_basic_boundary_1_rank
  test_basic_boundary_2_rank
  test_basic_boundary_4_rank
  test_cell_centered_fvm_1_rank
  test_cell_centered_fvm_2_rank
  test_cell_centered_fvm_4_rank
  test_fixed_step_piso_1_rank
  test_fixed_step_piso_2_rank
  test_fixed_step_piso_4_rank
  test_taylor_green_piso_1_rank
  test_taylor_green_piso_2_rank
  test_taylor_green_piso_4_rank
  test_material_density_wave_1_rank
  test_material_density_wave_2_rank
  test_material_density_wave_4_rank
  test_variable_density_vortex_full_1_rank
  test_variable_density_vortex_full_2_rank
  test_variable_density_vortex_full_4_rank
  test_ideal_gas_closed_heating_full_1_rank
  test_ideal_gas_closed_heating_full_2_rank
  test_ideal_gas_closed_heating_full_4_rank
  test_ideal_gas_open_plug_full_1_rank
  test_ideal_gas_open_plug_full_2_rank
  test_ideal_gas_open_plug_full_4_rank
  test_adaptive_time_control_1_rank_acceptance
  test_adaptive_time_control_2_rank_acceptance
  test_adaptive_time_control_4_rank_acceptance
  test_checkpoint_v2_acceptance_1_rank
  test_checkpoint_v2_acceptance_2_rank
  test_checkpoint_v2_acceptance_4_rank
  test_checkpoint_v2_diagnostics_acceptance_1_rank
  test_checkpoint_v2_diagnostics_acceptance_2_rank
  test_checkpoint_v2_diagnostics_acceptance_4_rank
  test_task24_flow_models_1rank
  test_task24_flow_models_2rank
  test_task24_flow_models_4rank
  test_task24_restart_1rank
  test_task24_restart_2rank
  test_task24_restart_4rank
  test_diagnostic_session_mpi_1rank
  test_diagnostic_session_mpi_2rank
  test_diagnostic_session_mpi_4rank
  test_stage2_module_diagnostics
  test_task25_performance_1rank
  test_task25_performance_2rank
  test_task25_performance_4rank
)
required_count=${#required_tests[@]}
if test "${required_count}" -ne 59; then
  fail "internal inventory cardinality mismatch"
fi

selection_regex='^('
separator=
for test_name in "${required_tests[@]}"; do
  selection_regex+="${separator}${test_name}"
  separator='|'
done
selection_regex+=')$'

if ! inventory_output=$(
  "${ctest_command}" --test-dir "${build_root}" -N -R "${selection_regex}" 2>&1
); then
  fail "CTest inventory query failed"
fi
case "${inventory_output}" in
  *"Could not find executable"*)
    fail "selected test executable is missing"
    ;;
esac
actual_tests=$(
  printf '%s\n' "${inventory_output}" |
    sed -n \
      's/^[[:space:]]*Test[[:space:]]*#[0-9][0-9]*:[[:space:]]*//p'
)
actual_count=$(
  printf '%s\n' "${actual_tests}" |
    awk 'NF { count += 1 } END { print count + 0 }'
)
if test "${actual_count}" -ne "${required_count}"; then
  fail "registered inventory cardinality mismatch: expected ${required_count}, found ${actual_count}"
fi
expected_sorted=$(
  printf '%s\n' "${required_tests[@]}" | LC_ALL=C sort
)
actual_sorted=$(
  printf '%s\n' "${actual_tests}" | LC_ALL=C sort
)
if test "${actual_sorted}" != "${expected_sorted}"; then
  fail "registered inventory name mismatch"
fi

if test "${HUNDUN_STAGE2_LIST_ONLY:-0}" = "1"; then
  printf 'STAGE2_ACCEPTANCE_INVENTORY cardinality=%s\n' "${required_count}"
  printf '%s\n' "${required_tests[@]}"
  exit 0
fi
if test "${HUNDUN_STAGE2_LIST_ONLY:-0}" != "0"; then
  fail "HUNDUN_STAGE2_LIST_ONLY must be 0 or 1"
fi

printf 'STAGE2_ACCEPTANCE_RUN cardinality=%s build=%s\n' \
  "${required_count}" "${build_root}"
"${ctest_command}" --test-dir "${build_root}" --output-on-failure \
  -j1 -R "${selection_regex}"
