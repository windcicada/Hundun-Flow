#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -uo pipefail

fail() { printf 'stage4 acceptance: %s\n' "$*" >&2; exit 2; }
source_root=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)

if test "${1:-}" = --list; then
  cat <<'EOF'
debug-stage4	ctest -L stage4	Debug	1/2/4-rank compact Stage 4 matrix
release-focused	checkpoint/coupling/state/transport/boundary/IBM/pressure/driver/diagnostics	Release	focused
asan-focused	coupling/state/transport/boundary/driver/checkpoint/diagnostics	ASan	serial-small
ubsan-focused	coupling/state/transport/boundary/driver/checkpoint/diagnostics	UBSan	serial-small
cantera-conformance	backend/thermo/transport/chemistry/0D/PSR	GCC11-rootfs	low-cost
package-rpath	stage4_package_contract	Debug	package/RPATH
cli-v4	--validate/--print-resolved	Debug	schema-v4
tests-off	configure-only	tests-off	no tests registered
EOF
  exit 0
fi

test "${1:-}" = --run || fail 'use --list or --run'
shift
candidate=
debug_root=
release_root=
asan_root=
ubsan_root=
tests_off_root=
cantera_root=
mechanism=
evidence_root=
while test "$#" -gt 0; do
  case "$1" in
    --candidate-head) candidate=${2:-}; shift 2 ;;
    --debug-root) debug_root=${2:-}; shift 2 ;;
    --release-root) release_root=${2:-}; shift 2 ;;
    --asan-root) asan_root=${2:-}; shift 2 ;;
    --ubsan-root) ubsan_root=${2:-}; shift 2 ;;
    --tests-off-root) tests_off_root=${2:-}; shift 2 ;;
    --cantera-root) cantera_root=${2:-}; shift 2 ;;
    --mechanism) mechanism=${2:-}; shift 2 ;;
    --evidence-dir) evidence_root=${2:-}; shift 2 ;;
    *) fail "unknown argument: $1" ;;
  esac
done
[[ "${candidate}" =~ ^[0-9a-f]{40}$ ]] || fail 'candidate HEAD must be lowercase SHA-1'
for path in "${debug_root}" "${release_root}" "${asan_root}" \
            "${ubsan_root}" "${tests_off_root}" "${cantera_root}" \
            "${mechanism}" "${evidence_root}"; do
  [[ "${path}" = /* ]] || fail "required path is not absolute: ${path}"
done
test "$(git -C "${source_root}" rev-parse HEAD)" = "${candidate}" || fail 'HEAD mismatch'
test -z "$(git -C "${source_root}" status --porcelain=v1 --untracked-files=normal)" || fail 'dirty candidate'
test "$(<"${source_root}/VERSION")" = 0.3.0 || fail 'VERSION is not 0.3.0'
for root in "${debug_root}" "${release_root}" "${asan_root}" "${ubsan_root}" "${tests_off_root}"; do
  test -f "${root}/CMakeCache.txt" || fail "missing cache: ${root}"
done
grep -Fq 'HUNDUN_BUILD_TESTS:BOOL=OFF' "${tests_off_root}/CMakeCache.txt" || fail 'tests-off cache mismatch'
test -x "${debug_root}/src/hundun" || fail 'missing Debug hundun'
test -f "${mechanism}" || fail 'missing Cantera mechanism'
mkdir -p -- "${evidence_root}"

stage4_listing=$(ctest --test-dir "${debug_root}" -N -L stage4)
if grep -Eiq '48\^?3|48x48|96\^?3|96x96|flame|tpdf|spray|sanitizer-large' <<<"${stage4_listing}"; then
  fail 'forbidden Stage 4 selector is registered'
fi

run_logged() {
  local name=$1
  shift
  local log="${evidence_root}/${name}.log"
  local started ended status sha
  started=$(date --iso-8601=seconds)
  "$@" >"${log}" 2>&1
  status=$?
  ended=$(date --iso-8601=seconds)
  sha=$(sha256sum "${log}" | awk '{print $1}')
  printf '%s\t%s\t%s\t%s\t%s\n' "${name}" "${status}" "${sha}" "${started}" "${ended}" \
    >>"${evidence_root}/results.tsv"
  test "${status}" -eq 0 || fail "selector failed: ${name}"
}

: >"${evidence_root}/results.tsv"
run_logged debug-stage4 ctest --test-dir "${debug_root}" --output-on-failure -L stage4
focused='^(test_checkpoint_v4|test_reacting_(coupling|mms|state|transport|boundary|driver)|test_reacting_(decomposition|immersed|closed_pressure|diagnostics)_[12]_rank)$'
run_logged release-focused ctest --test-dir "${release_root}" --output-on-failure -R "${focused}"
serial='^(test_checkpoint_v4|test_reacting_(coupling|mms|state|transport|boundary|driver)|test_reacting_diagnostics_1_rank)$'
run_logged asan-focused ctest --test-dir "${asan_root}" --output-on-failure -R "${serial}"
run_logged ubsan-focused ctest --test-dir "${ubsan_root}" --output-on-failure -R "${serial}"
run_logged cli-validate "${debug_root}/src/hundun" \
  "${source_root}/tests/fixtures/reacting_driver_case_v4.json" --validate
run_logged cli-print-resolved "${debug_root}/src/hundun" \
  "${source_root}/tests/fixtures/reacting_driver_case_v4.json" --print-resolved

for binary in test_cantera_backend test_cantera_thermodynamics \
              test_cantera_transport test_cantera_chemistry_interval \
              test_reacting_zero_dimensional test_reacting_psr; do
  test -x "${cantera_root}/${binary}" || fail "missing Cantera binary: ${binary}"
  run_logged "cantera-${binary}" "${cantera_root}/${binary}" "${mechanism}"
done

tree=$(git -C "${source_root}" rev-parse 'HEAD^{tree}')
binary_sha=$(sha256sum "${debug_root}/src/hundun" | awk '{print $1}')
results_sha=$(sha256sum "${evidence_root}/results.tsv" | awk '{print $1}')
printf 'candidate_head=%s\ntree=%s\nbinary_sha256=%s\nresults_sha256=%s\n' \
  "${candidate}" "${tree}" "${binary_sha}" "${results_sha}" \
  >"${evidence_root}/acceptance.identity"
printf 'STAGE4_ACCEPTANCE_MATRIX_PASS candidate=%s tree=%s\n' "${candidate}" "${tree}"
