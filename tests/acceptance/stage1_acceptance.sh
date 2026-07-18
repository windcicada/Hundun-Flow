#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
project_root=$(CDPATH= cd -- "${script_dir}/../.." && pwd -P)
cmake_command=${CMAKE_COMMAND:-cmake}
ctest_command=${CTEST_COMMAND:-ctest}
mpiexec_command=${MPIEXEC_COMMAND:-mpiexec}
ldd_command=${LDD_COMMAND:-ldd}
work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-stage1-acceptance.XXXXXX")
build_dir="${work_root}/build"
case_dir="${work_root}/case"
failure_case_dir="${work_root}/failure-case"
uneven_valid_case_dir="${work_root}/uneven-valid-partition"
uneven_automatic_case_dir="${work_root}/uneven-automatic-partition"
minimum_valid_case_dir="${work_root}/minimum-valid-partition"

cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

export OMPI_MCA_rmaps_base_oversubscribe=1

check_normal_log() {
  local path=$1
  local -a lines=()
  if ! mapfile -t lines <"${path}"; then
    echo "unable to read normal numerical log" >&2
    return 1
  fi
  if test "${#lines[@]}" -ne 5; then
    echo "normal numerical log must contain exactly five lines" >&2
    return 1
  fi

  local step_10='^STEP 10 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$'
  local step_20='^STEP 20 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$'
  local finished='^FINISHED step=20 time_s=[0-9.eE+-]+$'
  if test "${lines[0]}" != "HUNDUN-FLOW 0.0.0-stage1"; then
    echo "normal numerical log line 1 is not the Stage 1 banner" >&2
    return 1
  fi
  if test "${lines[1]}" != \
      "CASE name=periodic_passive_scalar ranks=2 cells=64x8x8"; then
    echo "normal numerical log line 2 is not the canonical CASE record" >&2
    return 1
  fi
  if [[ ! ${lines[2]} =~ ${step_10} ]]; then
    echo "normal numerical log line 3 is not the STEP 10 record" >&2
    return 1
  fi
  if [[ ! ${lines[3]} =~ ${step_20} ]]; then
    echo "normal numerical log line 4 is not the STEP 20 record" >&2
    return 1
  fi
  if [[ ! ${lines[4]} =~ ${finished} ]]; then
    echo "normal numerical log line 5 is not the FINISHED record" >&2
    return 1
  fi
}

require_bounded_diagnostic() {
  local path=$1
  local expected=$2
  local matches
  matches=$(grep -Fxc -- "${expected}" "${path}" || true)
  if test "${matches}" -ne 1; then
    echo "output failure did not emit the required diagnostic" >&2
    return 1
  fi
  local bytes
  bytes=$(wc -c <"${path}")
  if test "${bytes}" -gt 256; then
    echo "output failure diagnostic is not bounded" >&2
    return 1
  fi
}

require_single_newline_record() {
  local path=$1
  if ! test -s "${path}"; then
    echo "resolved output record is empty" >&2
    return 1
  fi
  local newline_count
  newline_count=$(LC_ALL=C tr -cd '\n' <"${path}" | wc -c)
  if test "${newline_count}" -ne 1; then
    echo "resolved output must contain exactly one newline" >&2
    return 1
  fi
  local final_byte
  final_byte=$(LC_ALL=C tail -c 1 -- "${path}" | od -An -tu1 | \
    tr -d '[:space:]')
  if test "${final_byte}" != "10"; then
    echo "resolved output must end with a newline" >&2
    return 1
  fi
}

require_no_numerical_artifacts() {
  local directory=$1
  local output_name=$2
  local restart_name=$3
  local artifact
  if test -e "${directory}/${output_name}" || \
      test -e "${directory}/${restart_name}"; then
    echo "zero-step partition case created a configured output directory" >&2
    return 1
  fi
  if ! artifact=$(find "${directory}" -type f \
      \( -name '*.vtk' -o -name 'restart.rank*.bin' -o -name '*.tmp' \) \
      -print -quit); then
    echo "unable to inspect zero-step partition artifacts" >&2
    return 1
  fi
  if test -n "${artifact}"; then
    echo "zero-step partition case created a numerical artifact" >&2
    return 1
  fi
}

check_invalid_local_width_case() {
  local label=$1
  local directory=$2
  local output_name=$3
  local restart_name=$4
  local validation_output="${work_root}/${label}.validation.txt"
  local normal_output="${work_root}/${label}.stdout"
  local normal_error="${work_root}/${label}.stderr"

  set +e
  timeout --kill-after=2 15 "${mpiexec_command}" -n 2 \
    "${build_dir}/hundun" "${directory}/case.json" --validate \
    >"${validation_output}"
  local validation_status=$?
  set -e
  if test "${validation_status}" -ne 0; then
    echo "${label} validation did not complete successfully" >&2
    return 1
  fi
  cmp --silent "${work_root}/validation.expected" "${validation_output}"

  : >"${normal_output}"
  : >"${normal_error}"
  set +e
  timeout --kill-after=2 5 "${mpiexec_command}" -n 2 sh -c \
    'normal_output=$1
normal_error=$2
shift 2
exec "$@" >>"${normal_output}" 2>>"${normal_error}"' \
    sh "${normal_output}" "${normal_error}" "${build_dir}/hundun" \
    "${directory}/case.json" >>"${normal_error}" 2>&1
  local normal_status=$?
  set -e
  case "${normal_status}" in
  1) ;;
  124 | 137)
    echo "${label} normal execution reached the time limit with status ${normal_status}" >&2
    return 1
    ;;
  *)
    echo "${label} normal execution returned unknown status ${normal_status}" >&2
    return 1
    ;;
  esac
  if test -s "${normal_output}"; then
    echo "${label} normal failure emitted stdout" >&2
    return 1
  fi
  local diagnostic_count
  diagnostic_count=$(grep -Fxc \
    'halo plan ghost width exceeds a local dimension' \
    "${normal_error}" || true)
  if test "${diagnostic_count}" -ne 1; then
    echo "${label} did not emit exactly one complete project diagnostic" >&2
    return 1
  fi
  local diagnostic_bytes
  diagnostic_bytes=$(wc -c <"${normal_error}")
  if test "${diagnostic_bytes}" -gt 4096; then
    echo "${label} stderr exceeded the fixed diagnostic bound" >&2
    return 1
  fi
  require_no_numerical_artifacts "${directory}" "${output_name}" \
    "${restart_name}"
}

check_minimum_valid_partition() {
  local normal_output="${work_root}/minimum-valid-partition.stdout"
  local normal_error="${work_root}/minimum-valid-partition.stderr"
  local expected_output="${work_root}/minimum-valid-partition.expected"

  set +e
  timeout --kill-after=2 15 "${mpiexec_command}" -n 2 \
    "${build_dir}/hundun" "${minimum_valid_case_dir}/case.json" \
    >"${normal_output}" 2>"${normal_error}"
  local normal_status=$?
  set -e
  if test "${normal_status}" -ne 0; then
    echo "minimum valid partition returned status ${normal_status}" >&2
    return 1
  fi
  printf '%s\n' \
    'HUNDUN-FLOW 0.0.0-stage1' \
    'CASE name=minimum_valid_partition ranks=2 cells=4x2x2' \
    'FINISHED step=0 time_s=0' >"${expected_output}"
  cmp --silent "${expected_output}" "${normal_output}"
  if test -s "${normal_error}"; then
    echo "minimum valid partition emitted stderr" >&2
    return 1
  fi
  require_no_numerical_artifacts "${minimum_valid_case_dir}" \
    'output.minimum' 'Restart.minimum'
}

artifact_inspection_missing="${work_root}/missing-artifact-inspection"
artifact_inspection_error="${work_root}/artifact-inspection.stderr"
set +e
require_no_numerical_artifacts "${artifact_inspection_missing}" \
  'output.artifact-inspection' 'Restart.artifact-inspection' \
  2>"${artifact_inspection_error}"
artifact_inspection_status=$?
set -e
if test "${artifact_inspection_status}" -eq 0; then
  echo "artifact inspection traversal error was incorrectly accepted" >&2
  exit 1
fi

if ! artifact_inspection_bytes=$(wc -c <"${artifact_inspection_error}"); then
  echo "unable to measure artifact inspection diagnostic" >&2
  exit 1
fi
if [[ ! ${artifact_inspection_bytes} =~ ^[0-9]+$ ]]; then
  echo "artifact inspection diagnostic size is invalid" >&2
  exit 1
fi
if test "${artifact_inspection_bytes}" -gt 1024; then
  echo "artifact inspection diagnostic is not bounded" >&2
  exit 1
fi

set +e
artifact_inspection_matches=$(grep -Fxc -- \
  'unable to inspect zero-step partition artifacts' \
  "${artifact_inspection_error}")
artifact_inspection_match_status=$?
set -e
case "${artifact_inspection_match_status}" in
0 | 1) ;;
*)
  echo "unable to count artifact inspection diagnostic" >&2
  exit 1
  ;;
esac
if test "${artifact_inspection_matches}" -ne 1; then
  echo "artifact inspection did not emit exactly one project diagnostic" >&2
  exit 1
fi

mkdir -p -- "${case_dir}" "${failure_case_dir}" \
  "${uneven_valid_case_dir}" "${uneven_automatic_case_dir}" \
  "${minimum_valid_case_dir}" \
  "${work_root}/launch-a" "${work_root}/launch-b"
cp -- "${project_root}/cases/passive_scalar/case.json" "${case_dir}/case.json"
cp -- "${project_root}/cases/passive_scalar/case_restart.json" \
  "${case_dir}/case_restart.json"
cp -- "${project_root}/cases/passive_scalar/case.json" \
  "${failure_case_dir}/case.json"

cat >"${uneven_valid_case_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "uneven_valid_partition"},
  "resources": {
    "expected_ranks": 2,
    "process_grid": [2, 1, 1]
  },
  "mesh": {
    "cells": [3, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 0.125, 0.125],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 0.003125, "steps": 0},
  "transport": {
    "velocity_m_per_s": [1.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart.uneven-valid"
  },
  "output": {
    "directory": "output.uneven-valid",
    "write_interval": 10,
    "restart_interval": 10
  }
}
EOF

cat >"${uneven_automatic_case_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "uneven_automatic_partition"},
  "resources": {"expected_ranks": 2},
  "mesh": {
    "cells": [3, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 0.125, 0.125],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 0.003125, "steps": 0},
  "transport": {
    "velocity_m_per_s": [1.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart.uneven-automatic"
  },
  "output": {
    "directory": "output.uneven-automatic",
    "write_interval": 10,
    "restart_interval": 10
  }
}
EOF

cat >"${minimum_valid_case_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "minimum_valid_partition"},
  "resources": {
    "expected_ranks": 2,
    "process_grid": [2, 1, 1]
  },
  "mesh": {
    "cells": [4, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 0.125, 0.125],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 0.003125, "steps": 0},
  "transport": {
    "velocity_m_per_s": [1.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart.minimum"
  },
  "output": {
    "directory": "output.minimum",
    "write_interval": 10,
    "restart_interval": 10
  }
}
EOF

"${cmake_command}" -S "${project_root}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=ON
"${cmake_command}" --build "${build_dir}" -j 2

"${build_dir}/hundun" --version >"${work_root}/version.txt"
printf '%s\n' "HUNDUN-FLOW 0.0.0-stage1" \
  >"${work_root}/version.expected"
cmp --silent "${work_root}/version.expected" "${work_root}/version.txt"
set +e
"${build_dir}/hundun" --version >/dev/full \
  2>"${work_root}/version-output-failure.log"
version_output_status=$?
set -e
if test "${version_output_status}" -eq 0; then
  echo "version output unexpectedly succeeded after a write failure" >&2
  exit 1
fi
require_bounded_diagnostic "${work_root}/version-output-failure.log" \
  "unable to write version output"

"${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case.json" --validate >"${work_root}/validation.txt"
printf '%s\n' "VALID" >"${work_root}/validation.expected"
cmp --silent "${work_root}/validation.expected" \
  "${work_root}/validation.txt"
set +e
timeout 15 "${mpiexec_command}" -n 2 sh -c \
  'exec "$@" >/dev/full' sh "${build_dir}/hundun" \
  "${case_dir}/case.json" --validate \
  2>"${work_root}/validation-output-failure.log"
validation_output_status=$?
set -e
if test "${validation_output_status}" -eq 0; then
  echo "validation output unexpectedly succeeded after a write failure" >&2
  exit 1
fi
if test "${validation_output_status}" -eq 124; then
  echo "validation output failure did not converge before timeout" >&2
  exit 1
fi
require_bounded_diagnostic "${work_root}/validation-output-failure.log" \
  "unable to write validation output"

(
  cd -- "${work_root}/launch-a"
  "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
    "${case_dir}/case.json" --print-resolved
) >"${work_root}/resolved-a.json"
(
  cd -- "${work_root}/launch-b"
  "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
    "${case_dir}/case.json" --print-resolved
) >"${work_root}/resolved-b.json"
require_single_newline_record "${work_root}/resolved-a.json"
require_single_newline_record "${work_root}/resolved-b.json"
cmp --silent "${work_root}/resolved-a.json" "${work_root}/resolved-b.json"
grep -F '"process_grid":[2,1,1]' "${work_root}/resolved-a.json" >/dev/null
grep -F '"restart":{"read":false,"write_directory":"Restart"}' \
  "${work_root}/resolved-a.json" >/dev/null

set +e
timeout 15 "${mpiexec_command}" -n 2 sh -c \
  'exec "$@" >/dev/full' sh "${build_dir}/hundun" \
  "${case_dir}/case.json" --print-resolved \
  2>"${work_root}/resolved-output-failure.log"
resolved_output_status=$?
set -e
if test "${resolved_output_status}" -eq 0; then
  echo "resolved output unexpectedly succeeded after a write failure" >&2
  exit 1
fi
if test "${resolved_output_status}" -eq 124; then
  echo "resolved output failure did not converge before timeout" >&2
  exit 1
fi
grep -F 'unable to write resolved case configuration' \
  "${work_root}/resolved-output-failure.log" >/dev/null

set +e
timeout 30 "${mpiexec_command}" -n 2 sh -c \
  'exec "$@" >/dev/full' sh "${build_dir}/hundun" \
  "${failure_case_dir}/case.json" \
  2>"${work_root}/normal-output-failure.log"
normal_output_status=$?
set -e
if test "${normal_output_status}" -eq 0; then
  echo "normal output unexpectedly succeeded after a write failure" >&2
  exit 1
fi
if test "${normal_output_status}" -eq 124; then
  echo "normal output failure did not converge before timeout" >&2
  exit 1
fi
require_bounded_diagnostic "${work_root}/normal-output-failure.log" \
  "unable to write Stage 1 output"

check_invalid_local_width_case 'uneven-valid-partition' \
  "${uneven_valid_case_dir}" 'output.uneven-valid' \
  'Restart.uneven-valid'
check_invalid_local_width_case 'uneven-automatic-partition' \
  "${uneven_automatic_case_dir}" 'output.uneven-automatic' \
  'Restart.uneven-automatic'
check_minimum_valid_partition

mkdir -p -- "${case_dir}/output"
"${mpiexec_command}" -n 2 "${build_dir}/hundun" "${case_dir}/case.json" \
  >"${case_dir}/output/run.log" 2>&1
check_normal_log "${case_dir}/output/run.log"

check_checkpoint() {
  local directory=$1
  test -s "${directory}/manifest.v1.bin"
  test -s "${directory}/COMPLETED"
  test "$(find "${directory}" -maxdepth 1 -type f \
    -name 'restart.rank*.bin' | wc -l)" -eq 2
  test -s "${directory}/restart.rank000000.bin"
  test -s "${directory}/restart.rank000001.bin"
  local records
  records=$(od -An -tu4 -j 75 -N 4 "${directory}/manifest.v1.bin" | \
    tr -d '[:space:]')
  test "${records}" = "2"
}

check_checkpoint "${case_dir}/Restart/step00000010"
check_checkpoint "${case_dir}/Restart/step00000020"
test "$(find "${case_dir}/output" -maxdepth 1 -type f \
  -name 'scalar.step*.rank*.vtk' | wc -l)" -eq 4

mv -- "${case_dir}/Restart/step00000010/COMPLETED" \
  "${case_dir}/Restart/step00000010/COMPLETED.saved"
if "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case_restart.json" >"${work_root}/missing-marker.log" 2>&1; then
  echo "restart unexpectedly accepted a missing COMPLETED marker" >&2
  exit 1
fi
mv -- "${case_dir}/Restart/step00000010/COMPLETED.saved" \
  "${case_dir}/Restart/step00000010/COMPLETED"

cp -p -- "${case_dir}/Restart/step00000010/restart.rank000001.bin" \
  "${case_dir}/Restart/step00000010/restart.rank000001.bin.saved"
: >"${case_dir}/Restart/step00000010/restart.rank000001.bin"
if "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case_restart.json" >"${work_root}/truncated-rank.log" 2>&1; then
  echo "restart unexpectedly accepted a truncated rank file" >&2
  exit 1
fi
mv -- "${case_dir}/Restart/step00000010/restart.rank000001.bin.saved" \
  "${case_dir}/Restart/step00000010/restart.rank000001.bin"

mkdir -p -- "${case_dir}/output.resumed"
"${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case_restart.json" >"${case_dir}/output.resumed/run.log" 2>&1
grep -E '^FINISHED step=20 time_s=[0-9.eE+-]+$' \
  "${case_dir}/output.resumed/run.log" >/dev/null
check_checkpoint "${case_dir}/Restart.resumed/step00000020"

for rank in 000000 000001; do
  "${build_dir}/compare_scalar_vtk" \
    "${case_dir}/output/scalar.step00000020.rank${rank}.vtk" \
    "${case_dir}/output.resumed/scalar.step00000020.rank${rank}.vtk" \
    1e-13
done
if "${build_dir}/compare_scalar_vtk" \
  "${case_dir}/output/scalar.step00000010.rank000000.vtk" \
  "${case_dir}/output/scalar.step00000020.rank000000.vtk" 1e-13; then
  echo "VTK comparator unexpectedly accepted distinct generated fields" >&2
  exit 1
fi

"${ctest_command}" --test-dir "${build_dir}" -j 1 --output-on-failure
"${ctest_command}" --test-dir "${build_dir}" -L provenance -j 1 \
  --output-on-failure

linkage_output="${work_root}/hundun.ldd"
if ! "${ldd_command}" "${build_dir}/hundun" >"${linkage_output}"; then
  echo "unable to inspect hundun dynamic linkage" >&2
  exit 1
fi
set +e
grep -E -i 'python|not found|coast|boffin' "${linkage_output}" >/dev/null
linkage_scan_status=$?
set -e
case "${linkage_scan_status}" in
0)
  echo "hundun has forbidden or missing dynamic linkage" >&2
  exit 1
  ;;
1) ;;
*)
  echo "unable to scan hundun dynamic linkage" >&2
  exit 1
  ;;
esac
if find "${case_dir}" -type f -name '*.tmp' -print -quit | grep -q .; then
  echo "generated case tree contains a residual temporary file" >&2
  exit 1
fi
