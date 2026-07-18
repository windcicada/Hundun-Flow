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
mpmd_root_a_dir="${work_root}/mpmd-root-a"
mpmd_root_b_dir="${work_root}/mpmd-root-b"
time_overflow_case_dir="${work_root}/time-overflow"
time_finite_case_dir="${work_root}/time-finite"
state_nonfinite_case_dir="${work_root}/state-nonfinite"

cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

export OMPI_MCA_rmaps_base_oversubscribe=1

check_normal_log() {
  local stdout_path=$1
  local stderr_path=$2
  if ! test -f "${stderr_path}"; then
    echo "normal numerical stderr capture is missing" >&2
    return 1
  fi
  if test -s "${stderr_path}"; then
    echo "normal numerical run emitted stderr" >&2
    return 1
  fi

  local newline_count
  if ! newline_count=$(LC_ALL=C tr -cd '\n' <"${stdout_path}" | wc -c); then
    echo "unable to count normal numerical stdout LF bytes" >&2
    return 1
  fi
  if [[ ! ${newline_count} =~ ^[0-9]+$ ]]; then
    echo "normal numerical stdout LF count is invalid" >&2
    return 1
  fi
  if test "${newline_count}" -ne 5; then
    echo "normal numerical stdout must contain exactly five LF bytes" >&2
    return 1
  fi

  local final_byte
  if ! final_byte=$(LC_ALL=C tail -c 1 -- "${stdout_path}" | od -An -tu1); then
    echo "unable to inspect normal numerical stdout final byte" >&2
    return 1
  fi
  final_byte=${final_byte//[[:space:]]/}
  if test "${final_byte}" != "10"; then
    echo "normal numerical stdout must end with LF" >&2
    return 1
  fi

  local -a lines=()
  if ! mapfile -t lines <"${stdout_path}"; then
    echo "unable to read normal numerical stdout" >&2
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

captured_generated_files=()

capture_generated_files() {
  local find_command=$1
  shift
  local captured
  captured_generated_files=()
  if ! captured=$("${find_command}" "$@"); then
    return 1
  fi
  if test -z "${captured}"; then
    return 0
  fi
  if ! mapfile -t captured_generated_files <<<"${captured}"; then
    captured_generated_files=()
    return 1
  fi
}

require_exact_canonical_vtk_set() {
  local directory=$1
  local enumerator=$2
  local -a expected_basenames=(
    'scalar.step00000010.rank000000.vtk'
    'scalar.step00000010.rank000001.vtk'
    'scalar.step00000020.rank000000.vtk'
    'scalar.step00000020.rank000001.vtk')

  if ! capture_generated_files "${enumerator}" "${directory}" -maxdepth 1 \
      -type f -name 'scalar.step*.rank*.vtk' -print; then
    echo "unable to enumerate canonical VTK files" >&2
    return 1
  fi
  if test "${#captured_generated_files[@]}" -ne \
      "${#expected_basenames[@]}"; then
    echo "canonical output does not contain exactly four VTK files" >&2
    return 1
  fi

  local expected
  for expected in "${expected_basenames[@]}"; do
    if ! test -s "${directory}/${expected}"; then
      echo "canonical output is missing nonempty VTK file ${expected}" >&2
      return 1
    fi
    local matches=0
    local generated
    for generated in "${captured_generated_files[@]}"; do
      if test "${generated##*/}" = "${expected}"; then
        matches=$((matches + 1))
      fi
    done
    if test "${matches}" -ne 1; then
      echo "canonical VTK enumeration does not contain ${expected} once" >&2
      return 1
    fi
  done

  local generated
  for generated in "${captured_generated_files[@]}"; do
    local basename=${generated##*/}
    local allowed=0
    for expected in "${expected_basenames[@]}"; do
      if test "${basename}" = "${expected}"; then
        allowed=1
        break
      fi
    done
    if test "${allowed}" -ne 1; then
      echo "canonical output contains stray VTK file ${basename}" >&2
      return 1
    fi
  done
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
  if ! newline_count=$(LC_ALL=C tr -cd '\n' <"${path}" | wc -c); then
    echo "unable to count resolved output LF bytes" >&2
    return 1
  fi
  if [[ ! ${newline_count} =~ ^[0-9]+$ ]]; then
    echo "resolved output LF count is invalid" >&2
    return 1
  fi
  if test "${newline_count}" -ne 1; then
    echo "resolved output must contain exactly one newline" >&2
    return 1
  fi
  local final_byte
  if ! final_byte=$(LC_ALL=C tail -c 1 -- "${path}" | od -An -tu1); then
    echo "unable to inspect resolved output final byte" >&2
    return 1
  fi
  final_byte=${final_byte//[[:space:]]/}
  if test "${final_byte}" != "10"; then
    echo "resolved output must end with a newline" >&2
    return 1
  fi
}

require_no_numerical_artifacts() {
  local directory=$1
  local output_name=$2
  local restart_name=$3
  if test -e "${directory}/${output_name}" || \
      test -e "${directory}/${restart_name}"; then
    echo "zero-step partition case created a configured output directory" >&2
    return 1
  fi
  if ! capture_generated_files find "${directory}" -type f \
      \( -name '*.vtk' -o -name 'restart.rank*.bin' -o -name '*.tmp' \) \
      -print; then
    echo "unable to inspect zero-step partition artifacts" >&2
    return 1
  fi
  if test "${#captured_generated_files[@]}" -ne 0; then
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

require_bounded_project_failure() {
  local label=$1
  local status=$2
  local stdout_path=$3
  local stderr_path=$4
  local diagnostic=$5
  if test "${status}" -eq 0; then
    echo "${label} unexpectedly completed successfully" >&2
    return 1
  fi
  case "${status}" in
  124 | 137)
    echo "${label} reached the time limit with status ${status}" >&2
    return 1
    ;;
  esac
  local diagnostic_count
  set +e
  diagnostic_count=$(grep -Fxc -- "${diagnostic}" "${stderr_path}")
  local diagnostic_status=$?
  set -e
  case "${diagnostic_status}" in
  0) ;;
  1) diagnostic_count=0 ;;
  *)
    echo "unable to count ${label} diagnostic" >&2
    return 1
    ;;
  esac
  if test "${diagnostic_count}" -ne 1; then
    echo "${label} did not emit its stable diagnostic exactly once" >&2
    return 1
  fi
  local diagnostic_bytes
  if ! diagnostic_bytes=$(wc -c <"${stderr_path}"); then
    echo "unable to measure ${label} diagnostic output" >&2
    return 1
  fi
  if [[ ! ${diagnostic_bytes} =~ ^[0-9]+$ ]]; then
    echo "${label} diagnostic size is invalid" >&2
    return 1
  fi
  if test "${diagnostic_bytes}" -gt 4096; then
    echo "${label} diagnostic output is not bounded" >&2
    return 1
  fi
  if grep -E '^(VALID|STEP|FINISHED)( |$)' \
      "${stdout_path}" >/dev/null; then
    echo "${label} emitted a success record" >&2
    return 1
  else
    local grep_status=$?
    if test "${grep_status}" -ne 1; then
      echo "unable to inspect ${label} stdout" >&2
      return 1
    fi
  fi
}

require_no_restart_rejection_artifacts() {
  local case_directory=$1
  local output_name=$2
  local restart_name=$3
  local enumerator=$4
  local root
  for root in "${case_directory}/${output_name}" \
      "${case_directory}/${restart_name}"; do
    if ! test -e "${root}"; then
      continue
    fi
    if ! test -d "${root}"; then
      echo "restart rejection artifact root is not a directory: ${root}" >&2
      return 1
    fi
    if ! capture_generated_files "${enumerator}" "${root}" -type f \
        \( -name '*.vtk' -o -name 'restart.rank*.bin' -o \
        -name 'manifest.v1.bin' -o -name 'COMPLETED' -o -name '*.tmp' \) \
        -print; then
      echo "unable to inspect restart rejection artifacts" >&2
      return 1
    fi
    if test "${#captured_generated_files[@]}" -ne 0; then
      echo "restart rejection published a forbidden artifact" >&2
      return 1
    fi
  done
}

require_restart_rejection_result() {
  local label=$1
  local status=$2
  local stdout_path=$3
  local stderr_path=$4
  local expected_diagnostic=$5
  local case_directory=$6
  local output_name=$7
  local restart_name=$8
  local enumerator=$9

  case "${status}" in
  1) ;;
  0)
    echo "${label} unexpectedly completed successfully" >&2
    return 1
    ;;
  124)
    echo "${label} reached the local time limit" >&2
    return 1
    ;;
  137)
    echo "${label} required a forced kill" >&2
    return 1
    ;;
  *)
    echo "${label} returned unrelated status ${status}" >&2
    return 1
    ;;
  esac

  if ! test -f "${stdout_path}" || ! test -f "${stderr_path}"; then
    echo "${label} did not produce separate output captures" >&2
    return 1
  fi

  local diagnostic_count
  local diagnostic_status
  set +e
  diagnostic_count=$(grep -Fc -- "${expected_diagnostic}" "${stderr_path}")
  diagnostic_status=$?
  set -e
  case "${diagnostic_status}" in
  0) ;;
  1) diagnostic_count=0 ;;
  *)
    echo "unable to inspect ${label} diagnostic" >&2
    return 1
    ;;
  esac
  if [[ ! ${diagnostic_count} =~ ^[0-9]+$ ]] || \
      test "${diagnostic_count}" -ne 1; then
    echo "${label} did not emit its condition-specific diagnostic once" >&2
    return 1
  fi

  local diagnostic_bytes
  if ! diagnostic_bytes=$(wc -c <"${stderr_path}"); then
    echo "unable to measure ${label} diagnostic" >&2
    return 1
  fi
  if [[ ! ${diagnostic_bytes} =~ ^[0-9]+$ ]] || \
      test "${diagnostic_bytes}" -gt 4096; then
    echo "${label} diagnostic is not within the fixed byte bound" >&2
    return 1
  fi

  set +e
  grep -E '^(VALID|STEP|FINISHED)( |$)' "${stdout_path}" >/dev/null
  local success_record_status=$?
  set -e
  case "${success_record_status}" in
  0)
    echo "${label} emitted a success record" >&2
    return 1
    ;;
  1) ;;
  *)
    echo "unable to inspect ${label} stdout" >&2
    return 1
    ;;
  esac

  require_no_restart_rejection_artifacts \
    "${case_directory}" "${output_name}" "${restart_name}" "${enumerator}"
}

run_restart_rejection_case() {
  local label=$1
  local time_limit=$2
  local kill_after=$3
  local expected_diagnostic=$4
  local case_directory=$5
  local output_name=$6
  local restart_name=$7
  local stdout_path=$8
  local stderr_path=$9
  shift 9

  : >"${stdout_path}"
  : >"${stderr_path}"
  set +e
  timeout --kill-after="${kill_after}" "${time_limit}" "$@" \
    >"${stdout_path}" 2>"${stderr_path}"
  local status=$?
  set -e
  require_restart_rejection_result \
    "${label}" "${status}" "${stdout_path}" "${stderr_path}" \
    "${expected_diagnostic}" "${case_directory}" "${output_name}" \
    "${restart_name}" find
}

require_no_nonfinite_token() {
  local path=$1
  set +e
  grep -E -i '(^|[^[:alnum:]_])(inf|nan)([^[:alnum:]_]|$)' \
    "${path}" >/dev/null
  local status=$?
  set -e
  case "${status}" in
  0)
    echo "successful numerical output contains a nonfinite token" >&2
    return 1
    ;;
  1) ;;
  *)
    echo "unable to inspect successful numerical output" >&2
    return 1
    ;;
  esac
}

fixture_find_failure_empty() {
  return 23
}

fixture_find_failure_two() {
  printf '%s\n' 'rank-a.bin' 'rank-b.bin'
  return 23
}

fixture_find_failure_four() {
  printf '%s\n' 'field-a.vtk' 'field-b.vtk' 'field-c.vtk' 'field-d.vtk'
  return 23
}

fixture_find_success_empty() {
  return 0
}

fixture_find_success_two() {
  printf '%s\n' 'rank-a.bin' 'rank-b.bin'
}

fixture_find_success_four() {
  printf '%s\n' 'field-a.vtk' 'field-b.vtk' 'field-c.vtk' 'field-d.vtk'
}

fixture_vtk_enumerator_failure() {
  local directory=$1
  printf '%s\n' \
    "${directory}/scalar.step00000010.rank000000.vtk" \
    "${directory}/scalar.step00000020.rank000001.vtk" \
    "${directory}/scalar.step00000010.rank000001.vtk" \
    "${directory}/scalar.step00000020.rank000000.vtk"
  return 23
}

fixture_expect_generated_count() {
  local expected=$1
  local find_command=$2
  shift 2
  if ! capture_generated_files "${find_command}" "$@"; then
    return 1
  fi
  test "${#captured_generated_files[@]}" -eq "${expected}"
}

run_fast_shell_contract_fixtures() {
  local stdout_path="${work_root}/normal-log-fixture.stdout"
  local stderr_path="${work_root}/normal-log-fixture.stderr"
  local no_final_lf="${work_root}/normal-log-no-final-lf.stdout"
  local moved_stdout="${work_root}/normal-log-moved.stdout"
  local moved_stderr="${work_root}/normal-log-moved.stderr"
  local extra_stderr="${work_root}/normal-log-extra.stderr"

  printf '%s\n' \
    'HUNDUN-FLOW 0.0.0-stage1' \
    'CASE name=periodic_passive_scalar ranks=2 cells=64x8x8' \
    'STEP 10 time_s=0.1 mass=1 relative_mass_error=0' \
    'STEP 20 time_s=0.2 mass=1 relative_mass_error=0' \
    'FINISHED step=20 time_s=0.2' >"${stdout_path}"
  : >"${stderr_path}"
  check_normal_log "${stdout_path}" "${stderr_path}"

  printf '%s\n' \
    'HUNDUN-FLOW 0.0.0-stage1' \
    'CASE name=periodic_passive_scalar ranks=2 cells=64x8x8' \
    'STEP 10 time_s=0.1 mass=1 relative_mass_error=0' \
    'STEP 20 time_s=0.2 mass=1 relative_mass_error=0' >"${no_final_lf}"
  printf '%s' 'FINISHED step=20 time_s=0.2' >>"${no_final_lf}"
  if check_normal_log "${no_final_lf}" "${stderr_path}" >/dev/null 2>&1; then
    echo "normal-log contract accepted a missing final LF" >&2
    return 1
  fi

  printf '%s\n' \
    'HUNDUN-FLOW 0.0.0-stage1' \
    'CASE name=periodic_passive_scalar ranks=2 cells=64x8x8' \
    'STEP 20 time_s=0.2 mass=1 relative_mass_error=0' \
    'FINISHED step=20 time_s=0.2' >"${moved_stdout}"
  printf '%s\n' \
    'STEP 10 time_s=0.1 mass=1 relative_mass_error=0' >"${moved_stderr}"
  if check_normal_log "${moved_stdout}" "${moved_stderr}" \
      >/dev/null 2>&1; then
    echo "normal-log contract accepted a STEP record on stderr" >&2
    return 1
  fi

  printf '%s\n' 'arbitrary stderr record' >"${extra_stderr}"
  if check_normal_log "${stdout_path}" "${extra_stderr}" \
      >/dev/null 2>&1; then
    echo "normal-log contract accepted nonempty stderr" >&2
    return 1
  fi

  if fixture_expect_generated_count 0 fixture_find_failure_empty; then
    echo "generated-file capture accepted an empty failing traversal" >&2
    return 1
  fi
  if fixture_expect_generated_count 2 fixture_find_failure_two; then
    echo "two-file count accepted a partial failing traversal" >&2
    return 1
  fi
  if fixture_expect_generated_count 4 fixture_find_failure_four; then
    echo "four-file count accepted a partial failing traversal" >&2
    return 1
  fi
  if ! fixture_expect_generated_count 0 fixture_find_success_empty; then
    echo "expected-zero generated-file control was rejected" >&2
    return 1
  fi
  if ! fixture_expect_generated_count 2 fixture_find_success_two; then
    echo "exact-two generated-file control was rejected" >&2
    return 1
  fi
  if ! fixture_expect_generated_count 4 fixture_find_success_four; then
    echo "exact-four generated-file control was rejected" >&2
    return 1
  fi

  local rejection_root="${work_root}/restart-rejection-fixture"
  local rejection_stdout="${work_root}/restart-rejection-fixture.stdout"
  local rejection_stderr="${work_root}/restart-rejection-fixture.stderr"
  local rejection_diagnostic='fixture restart rejection'
  mkdir -p -- "${rejection_root}"

  if ! run_restart_rejection_case \
      'intended restart rejection fixture' 1 1 \
      "${rejection_diagnostic}" "${rejection_root}" output.resumed \
      Restart.resumed "${rejection_stdout}" "${rejection_stderr}" \
      sh -c 'printf "%s\n" "$1" >&2; exit 1' sh \
      "${rejection_diagnostic}"; then
    echo "intended restart rejection fixture was rejected" >&2
    return 1
  fi

  if run_restart_rejection_case \
      'unrelated restart rejection status fixture' 1 1 \
      "${rejection_diagnostic}" "${rejection_root}" output.resumed \
      Restart.resumed "${rejection_stdout}" "${rejection_stderr}" \
      sh -c 'printf "%s\n" "$1" >&2; exit 42' sh \
      "${rejection_diagnostic}" >/dev/null 2>&1; then
    echo "restart rejection fixture accepted unrelated status 42" >&2
    return 1
  fi

  if run_restart_rejection_case \
      'restart rejection timeout fixture' 0.2 0.1 \
      "${rejection_diagnostic}" "${rejection_root}" output.resumed \
      Restart.resumed "${rejection_stdout}" "${rejection_stderr}" \
      sh -c 'sleep 2' >/dev/null 2>&1; then
    echo "restart rejection fixture accepted a timeout" >&2
    return 1
  fi

  if run_restart_rejection_case \
      'restart rejection success-record fixture' 1 1 \
      "${rejection_diagnostic}" "${rejection_root}" output.resumed \
      Restart.resumed "${rejection_stdout}" "${rejection_stderr}" \
      sh -c 'printf "%s\n" "FINISHED step=20 time_s=1"; printf "%s\n" "$1" >&2; exit 1' \
      sh "${rejection_diagnostic}" >/dev/null 2>&1; then
    echo "restart rejection fixture accepted a success record" >&2
    return 1
  fi

  local -a forbidden_rejection_artifacts=(
    'output.resumed/stray.vtk'
    'Restart.resumed/restart.rank000000.bin'
    'Restart.resumed/manifest.v1.bin'
    'Restart.resumed/COMPLETED'
    'Restart.resumed/partial.tmp')
  local artifact
  for artifact in "${forbidden_rejection_artifacts[@]}"; do
    rm -rf -- "${rejection_root}/output.resumed" \
      "${rejection_root}/Restart.resumed"
    if run_restart_rejection_case \
        "restart rejection artifact fixture ${artifact}" 1 1 \
        "${rejection_diagnostic}" "${rejection_root}" output.resumed \
        Restart.resumed "${rejection_stdout}" "${rejection_stderr}" \
        sh -c 'mkdir -p -- "$(dirname -- "$1")"; : >"$1"; printf "%s\n" "$2" >&2; exit 1' \
        sh "${rejection_root}/${artifact}" "${rejection_diagnostic}" \
        >/dev/null 2>&1; then
      echo "restart rejection fixture accepted artifact ${artifact}" >&2
      return 1
    fi
  done
  rm -rf -- "${rejection_root}/output.resumed" \
    "${rejection_root}/Restart.resumed"

  mkdir -p -- "${rejection_root}/output.resumed"
  : >"${rejection_stdout}"
  printf '%s\n' "${rejection_diagnostic}" >"${rejection_stderr}"
  if require_restart_rejection_result \
      'restart rejection inspection-error fixture' 1 \
      "${rejection_stdout}" "${rejection_stderr}" \
      "${rejection_diagnostic}" "${rejection_root}" output.resumed \
      Restart.resumed fixture_find_failure_empty >/dev/null 2>&1; then
    echo "restart rejection fixture accepted an inspection error" >&2
    return 1
  fi

  local vtk_fixture_root="${work_root}/canonical-vtk-set-fixture"
  local -a expected_vtk_basenames=(
    'scalar.step00000010.rank000000.vtk'
    'scalar.step00000010.rank000001.vtk'
    'scalar.step00000020.rank000000.vtk'
    'scalar.step00000020.rank000001.vtk')
  local vtk_name
  rm -rf -- "${vtk_fixture_root}"
  mkdir -p -- "${vtk_fixture_root}"
  for vtk_name in "${expected_vtk_basenames[@]}"; do
    printf '%s\n' 'VTK fixture' >"${vtk_fixture_root}/${vtk_name}"
  done
  if ! require_exact_canonical_vtk_set "${vtk_fixture_root}" find; then
    echo "exact canonical VTK-set fixture was rejected" >&2
    return 1
  fi

  local replacement_index
  for replacement_index in 0 1 2 3; do
    rm -rf -- "${vtk_fixture_root}"
    mkdir -p -- "${vtk_fixture_root}"
    local expected_index
    for expected_index in 0 1 2 3; do
      if test "${expected_index}" -eq "${replacement_index}"; then
        vtk_name="scalar.step9999999${replacement_index}.rank99999${replacement_index}.vtk"
      else
        vtk_name=${expected_vtk_basenames[expected_index]}
      fi
      printf '%s\n' 'VTK fixture' >"${vtk_fixture_root}/${vtk_name}"
    done
    if require_exact_canonical_vtk_set "${vtk_fixture_root}" find \
        >/dev/null 2>&1; then
      echo "canonical VTK-set fixture accepted substitution ${replacement_index}" >&2
      return 1
    fi
  done

  rm -rf -- "${vtk_fixture_root}"
  mkdir -p -- "${vtk_fixture_root}"
  for vtk_name in "${expected_vtk_basenames[@]}"; do
    printf '%s\n' 'VTK fixture' >"${vtk_fixture_root}/${vtk_name}"
  done
  if require_exact_canonical_vtk_set \
      "${vtk_fixture_root}" fixture_vtk_enumerator_failure \
      >/dev/null 2>&1; then
    echo "canonical VTK-set fixture accepted a failing enumerator" >&2
    return 1
  fi
}

run_fast_shell_contract_fixtures

if test "${HUNDUN_ACCEPTANCE_FAST_FIXTURES_ONLY:-0}" = 1; then
  exit 0
fi

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
  "${minimum_valid_case_dir}" "${mpmd_root_a_dir}" \
  "${mpmd_root_b_dir}" "${time_overflow_case_dir}" \
  "${time_finite_case_dir}" "${state_nonfinite_case_dir}" \
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

cat >"${mpmd_root_a_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "rank_zero_authoritative_root"},
  "resources": {
    "expected_ranks": 2,
    "process_grid": [2, 1, 1]
  },
  "mesh": {
    "cells": [4, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 1.0, 1.0],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 0.125, "steps": 1},
  "transport": {
    "velocity_m_per_s": [0.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart"
  },
  "output": {
    "directory": "output",
    "write_interval": 1,
    "restart_interval": 2
  }
}
EOF
cp -- "${mpmd_root_a_dir}/case.json" "${mpmd_root_b_dir}/case.json"

cat >"${time_overflow_case_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "time_overflow"},
  "resources": {"expected_ranks": 1},
  "mesh": {
    "cells": [2, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 1.0, 1.0],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 1.7976931348623157e308, "steps": 2},
  "transport": {
    "velocity_m_per_s": [0.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart"
  },
  "output": {
    "directory": "output",
    "write_interval": 2,
    "restart_interval": 3
  }
}
EOF

cat >"${time_finite_case_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "time_finite_near_limit"},
  "resources": {"expected_ranks": 1},
  "mesh": {
    "cells": [2, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 1.0, 1.0],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 8.0e307, "steps": 2},
  "transport": {
    "velocity_m_per_s": [0.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart"
  },
  "output": {
    "directory": "output",
    "write_interval": 2,
    "restart_interval": 3
  }
}
EOF

cat >"${state_nonfinite_case_dir}/case.json" <<'EOF'
{
  "schema_version": 1,
  "case": {"name": "state_nonfinite"},
  "resources": {"expected_ranks": 1},
  "mesh": {
    "cells": [4, 2, 2],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 1.0, 1.0],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 1.0, "steps": 1},
  "transport": {
    "velocity_m_per_s": [1.7976931348623157e308, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {
    "read": false,
    "write_directory": "Restart"
  },
  "output": {
    "directory": "output",
    "write_interval": 2,
    "restart_interval": 1
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

mode_stdout="${work_root}/mode-mismatch.stdout"
mode_stderr="${work_root}/mode-mismatch.stderr"
set +e
timeout --kill-after=2 10 "${mpiexec_command}" \
  -n 1 "${build_dir}/hundun" "${mpmd_root_a_dir}/case.json" --validate \
  : -n 1 "${build_dir}/hundun" "${mpmd_root_a_dir}/case.json" \
  >"${mode_stdout}" 2>"${mode_stderr}"
mode_status=$?
set -e
require_bounded_project_failure \
  "communicator run-mode mismatch" "${mode_status}" \
  "${mode_stdout}" "${mode_stderr}" \
  "MPI run mode differs across communicator ranks"
set +e
grep -E '^(VALID|HUNDUN-FLOW|CASE|STEP|FINISHED)( |$)' \
  "${mode_stdout}" >/dev/null
mode_success_record_status=$?
set -e
case "${mode_success_record_status}" in
0)
  echo "communicator run-mode mismatch emitted a success record" >&2
  exit 1
  ;;
1) ;;
*)
  echo "unable to inspect communicator run-mode stdout" >&2
  exit 1
  ;;
esac
require_no_numerical_artifacts "${mpmd_root_a_dir}" output Restart

root_stdout="${work_root}/authoritative-root.stdout"
root_stderr="${work_root}/authoritative-root.stderr"
set +e
timeout --kill-after=2 20 "${mpiexec_command}" \
  -n 1 "${build_dir}/hundun" "${mpmd_root_a_dir}/case.json" \
  : -n 1 "${build_dir}/hundun" "${mpmd_root_b_dir}/case.json" \
  >"${root_stdout}" 2>"${root_stderr}"
root_status=$?
set -e
if test "${root_status}" -ne 0; then
  echo "authoritative-root MPMD run returned status ${root_status}" >&2
  exit 1
fi
if test -s "${root_stderr}"; then
  echo "authoritative-root MPMD run emitted stderr" >&2
  exit 1
fi
grep -Fx 'HUNDUN-FLOW 0.0.0-stage1' "${root_stdout}" >/dev/null
grep -Fx 'CASE name=rank_zero_authoritative_root ranks=2 cells=4x2x2' \
  "${root_stdout}" >/dev/null
grep -E '^STEP 1 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$' \
  "${root_stdout}" >/dev/null
grep -E '^FINISHED step=1 time_s=[0-9.eE+-]+$' \
  "${root_stdout}" >/dev/null
test -s "${mpmd_root_a_dir}/output/scalar.step00000001.rank000000.vtk"
test -s "${mpmd_root_a_dir}/output/scalar.step00000001.rank000001.vtk"
if ! capture_generated_files find "${mpmd_root_a_dir}/output" -maxdepth 1 \
    -type f -name 'scalar.step*.rank*.vtk' -print; then
  echo "unable to enumerate authoritative-root VTK files" >&2
  exit 1
fi
if test "${#captured_generated_files[@]}" -ne 2; then
  echo "authoritative root must contain exactly two VTK files" >&2
  exit 1
fi
if test -e "${mpmd_root_a_dir}/Restart"; then
  echo "one-step authoritative-root control created Restart output" >&2
  exit 1
fi
require_no_numerical_artifacts "${mpmd_root_b_dir}" output Restart

time_overflow_stdout="${work_root}/time-overflow.stdout"
time_overflow_stderr="${work_root}/time-overflow.stderr"
set +e
timeout --kill-after=2 10 "${build_dir}/hundun" \
  "${time_overflow_case_dir}/case.json" \
  >"${time_overflow_stdout}" 2>"${time_overflow_stderr}"
time_overflow_status=$?
set -e
require_bounded_project_failure \
  "nonfinite accumulated physical time" "${time_overflow_status}" \
  "${time_overflow_stdout}" "${time_overflow_stderr}" \
  "physical time is not finite"
require_no_numerical_artifacts "${time_overflow_case_dir}" output Restart

time_finite_stdout="${work_root}/time-finite.stdout"
time_finite_stderr="${work_root}/time-finite.stderr"
set +e
timeout --kill-after=2 10 "${build_dir}/hundun" \
  "${time_finite_case_dir}/case.json" \
  >"${time_finite_stdout}" 2>"${time_finite_stderr}"
time_finite_status=$?
set -e
if test "${time_finite_status}" -ne 0; then
  echo "finite near-limit time control returned status ${time_finite_status}" >&2
  exit 1
fi
if test -s "${time_finite_stderr}"; then
  echo "finite near-limit time control emitted stderr" >&2
  exit 1
fi
grep -E '^STEP 2 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$' \
  "${time_finite_stdout}" >/dev/null
grep -E '^FINISHED step=2 time_s=[0-9.eE+-]+$' \
  "${time_finite_stdout}" >/dev/null
require_no_nonfinite_token "${time_finite_stdout}"
test -s "${time_finite_case_dir}/output/scalar.step00000002.rank000000.vtk"
if test -e "${time_finite_case_dir}/Restart"; then
  echo "finite near-limit time control created Restart output" >&2
  exit 1
fi

state_stdout="${work_root}/state-nonfinite.stdout"
state_stderr="${work_root}/state-nonfinite.stderr"
set +e
timeout --kill-after=2 10 "${build_dir}/hundun" \
  "${state_nonfinite_case_dir}/case.json" \
  >"${state_stdout}" 2>"${state_stderr}"
state_status=$?
set -e
require_bounded_project_failure \
  "nonfinite passive-scalar state" "${state_status}" \
  "${state_stdout}" "${state_stderr}" \
  "passive-scalar state is not finite"
require_no_numerical_artifacts "${state_nonfinite_case_dir}" output Restart

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
  >"${case_dir}/output/run.log" 2>"${case_dir}/output/run.stderr"
check_normal_log "${case_dir}/output/run.log" \
  "${case_dir}/output/run.stderr"

check_checkpoint() {
  local directory=$1
  test -s "${directory}/manifest.v1.bin"
  test -s "${directory}/COMPLETED"
  if ! capture_generated_files find "${directory}" -maxdepth 1 -type f \
      -name 'restart.rank*.bin' -print; then
    echo "unable to enumerate rank checkpoint files" >&2
    return 1
  fi
  if test "${#captured_generated_files[@]}" -ne 2; then
    echo "rank checkpoint directory must contain exactly two rank files" >&2
    return 1
  fi
  test -s "${directory}/restart.rank000000.bin"
  test -s "${directory}/restart.rank000001.bin"
  local records
  records=$(od -An -tu4 -j 75 -N 4 "${directory}/manifest.v1.bin" | \
    tr -d '[:space:]')
  test "${records}" = "2"
}

check_checkpoint "${case_dir}/Restart/step00000010"
check_checkpoint "${case_dir}/Restart/step00000020"
require_exact_canonical_vtk_set "${case_dir}/output" find

check_missing_marker_rejection() {
  local marker="${case_dir}/Restart/step00000010/COMPLETED"
  local saved_marker="${marker}.saved"
  local rejection_status=0
  local restore_status=0

  mv -- "${marker}" "${saved_marker}"
  run_restart_rejection_case \
    'missing Restart completion marker' 20 2 \
    'Restart v1 completion marker read: unable to open Restart v1 file:' \
    "${case_dir}" output.resumed Restart.resumed \
    "${work_root}/missing-marker.stdout" \
    "${work_root}/missing-marker.stderr" \
    "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
    "${case_dir}/case_restart.json" || rejection_status=$?
  if ! mv -- "${saved_marker}" "${marker}"; then
    echo "unable to restore Restart completion marker fixture" >&2
    restore_status=1
  fi
  if test "${rejection_status}" -ne 0; then
    return "${rejection_status}"
  fi
  return "${restore_status}"
}

check_truncated_rank_rejection() {
  local rank_file=
  rank_file="${case_dir}/Restart/step00000010/restart.rank000001.bin"
  local saved_rank_file="${rank_file}.saved"
  local rejection_status=0
  local restore_status=0

  cp -p -- "${rank_file}" "${saved_rank_file}"
  : >"${rank_file}"
  run_restart_rejection_case \
    'truncated Restart rank file' 20 2 \
    'Restart v1 rank-file read: Restart v1 rank-file size or CRC does not match manifest' \
    "${case_dir}" output.resumed Restart.resumed \
    "${work_root}/truncated-rank.stdout" \
    "${work_root}/truncated-rank.stderr" \
    "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
    "${case_dir}/case_restart.json" || rejection_status=$?
  if ! mv -- "${saved_rank_file}" "${rank_file}"; then
    echo "unable to restore truncated Restart rank-file fixture" >&2
    restore_status=1
  fi
  if test "${rejection_status}" -ne 0; then
    return "${rejection_status}"
  fi
  return "${restore_status}"
}

check_missing_marker_rejection
check_truncated_rank_rejection

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
if ! capture_generated_files find "${case_dir}" -type f -name '*.tmp' \
    -print; then
  echo "unable to inspect generated case temporary files" >&2
  exit 1
fi
if test "${#captured_generated_files[@]}" -ne 0; then
  echo "generated case tree contains a residual temporary file" >&2
  exit 1
fi
