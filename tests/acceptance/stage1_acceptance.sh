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

mkdir -p -- "${case_dir}" "${failure_case_dir}" \
  "${work_root}/launch-a" "${work_root}/launch-b"
cp -- "${project_root}/cases/passive_scalar/case.json" "${case_dir}/case.json"
cp -- "${project_root}/cases/passive_scalar/case_restart.json" \
  "${case_dir}/case_restart.json"
cp -- "${project_root}/cases/passive_scalar/case.json" \
  "${failure_case_dir}/case.json"

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
