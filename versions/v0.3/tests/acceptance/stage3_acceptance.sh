#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -uo pipefail

fail() { printf 'stage3 acceptance: %s\n' "$*" >&2; exit 2; }
source_root=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
inventory="${source_root}/tests/acceptance/stage3_acceptance_inventory.tsv"
mode=
group=
candidate=
debug_root=
release_root=
asan_root=
ubsan_root=
tests_off_root=
while test "$#" -gt 0; do
  case "$1" in
    --list) mode=list; shift ;;
    --group) test "$#" -ge 2 || fail 'missing --group value'; mode=group; group=$2; shift 2 ;;
    --candidate-head) candidate=${2:-}; shift 2 ;;
    --debug-root) debug_root=${2:-}; shift 2 ;;
    --release-root) release_root=${2:-}; shift 2 ;;
    --asan-root) asan_root=${2:-}; shift 2 ;;
    --ubsan-root) ubsan_root=${2:-}; shift 2 ;;
    --tests-off-root) tests_off_root=${2:-}; shift 2 ;;
    *) fail "unknown argument: $1" ;;
  esac
done
test -n "${mode}" || fail 'use --list or --group'
if test "${mode}" = list; then
  LC_ALL=C sort -t $'\t' -k2,2 -k1,1 "${inventory}"
  exit 0
fi
case "${group}" in low-cost|scientific|performance|sanitizer|governance) ;; *) fail 'invalid group' ;; esac
[[ "${candidate}" =~ ^[0-9a-f]{40}$ ]] || fail 'candidate head must be 40 lowercase hex'
for root in "${debug_root}" "${release_root}" "${asan_root}" "${ubsan_root}" "${tests_off_root}"; do
  [[ "${root}" = /* ]] || fail 'all five build roots must be absolute'
  test -f "${root}/CMakeCache.txt" || fail "missing CMakeCache.txt: ${root}"
done
observed=$(git -C "${source_root}" rev-parse HEAD) || fail 'cannot read source HEAD'
test "${observed}" = "${candidate}" || fail 'candidate HEAD mismatch'
test -z "$(git -C "${source_root}" status --porcelain=v1 --untracked-files=normal)" || fail 'candidate worktree is not clean'
evidence_root=${HUNDUN_STAGE3_EVIDENCE_DIR:-}
[[ "${evidence_root}" = /* ]] || fail 'HUNDUN_STAGE3_EVIDENCE_DIR must be absolute'
tree=$(git -C "${source_root}" rev-parse 'HEAD^{tree}') || fail 'cannot read tree identity'
cpuset=$(sed -n 's/^Cpus_allowed_list:[[:space:]]*//p' /proc/self/status)
overall=0
while IFS=$'\t' read -r row row_group required resource ranks role producer executable argv timeout_s artifact_subdir; do
  test "${row}" != row_id || continue
  test "${row_group}" = "${group}" || continue
  case "${role}" in
    debug) root=${debug_root} ;;
    release) root=${release_root} ;;
    asan) root=${asan_root} ;;
    ubsan) root=${ubsan_root} ;;
    tests-off) root=${tests_off_root} ;;
    source-only) root=${debug_root} ;;
    *) fail "invalid build role in ${row}" ;;
  esac
  root=$(CDPATH= cd -- "${root}" && pwd -P) || fail "noncanonical build root: ${root}"
  cache=${root}/CMakeCache.txt
  cache_sha=$(sha256sum "${cache}" | awk '{print $1}')
  case "${role}" in
    tests-off) grep -Fq 'HUNDUN_BUILD_TESTS:BOOL=OFF' "${cache}" || fail "tests-off mismatch: ${row}" ;;
    source-only) ;;
    *) grep -Fq 'HUNDUN_BUILD_TESTS:BOOL=ON' "${cache}" || fail "tests-enabled mismatch: ${row}" ;;
  esac
  binary=${root}/tests/${producer}
  test -f "${binary}" || binary=${root}/src/${producer}
  test -f "${binary}" || binary=${cache}
  binary_sha=$(sha256sum "${binary}" | awk '{print $1}')
  binary_inode=$(stat -c %i "${binary}")
  row_dir=${evidence_root}/${artifact_subdir}
  mkdir -p -- "${row_dir}"
  log=${row_dir}/run.log
  time_log=${row_dir}/time.log
  started=$(date --iso-8601=seconds)
  start_epoch=$(date +%s)
  export HUNDUN_STAGE3_TREE_FINGERPRINT="git-tree:${tree}"
  export HUNDUN_STAGE3_BINARY_FINGERPRINT="sha256:${binary_sha}"
  export HUNDUN_STAGE3_CPUSET="${cpuset}"
  export HUNDUN_STAGE3_THREAD_BUDGET=1
  export HUNDUN_STAGE3_EVIDENCE_DIR="${row_dir}"
  row_environment=thread_budget=1
  command_environment=(env)
  if test "${role}" = asan; then
    command_environment+=(ASAN_OPTIONS=detect_leaks=0)
    row_environment+=';ASAN_OPTIONS=detect_leaks=0'
  fi
  if test "${executable}" = @ctest; then
    command=(ctest --test-dir "${root}" --output-on-failure -R "^${argv}$")
  else
    command=("${root}/${executable}" "${argv}")
  fi
  /usr/bin/time -v -o "${time_log}" timeout "${timeout_s}" \
    "${command_environment[@]}" "${command[@]}" </dev/null >"${log}" 2>&1
  status=$?
  ended=$(date --iso-8601=seconds)
  end_epoch=$(date +%s)
  duration=$((end_epoch - start_epoch))
  peak=$(sed -n 's/^[[:space:]]*Maximum resident set size (kbytes):[[:space:]]*//p' "${time_log}")
  peak=${peak:-0}
  log_sha=$(sha256sum "${log}" | awk '{print $1}')
  artifact_sha=$(sha256sum "${log}" | awk '{print $1}')
  manifest=${row_dir}/terminal-manifest.v1.json
  printf '{"schema_version":1,"row_id":"%s","candidate_head":"%s","tree_fingerprint":"git-tree:%s","diff_fingerprint":"clean","build_role":"%s","build_root":"%s","cache_sha256":"%s","binary_sha256":"%s","binary_inode":%s,"compiler_identity":"configured","libcxx_identity":"configured","mpi_identity":"configured","argv":"%s","environment":"%s","cpuset":"%s","ranks":%s,"started_at":"%s","ended_at":"%s","exit_status":%s,"duration_seconds":%s,"peak_rss_kib":%s,"log_sha256":"%s","artifact_sha256":["%s"]}\n' \
    "${row}" "${candidate}" "${tree}" "${role}" "${root}" "${cache_sha}" "${binary_sha}" "${binary_inode}" "${argv}" "${row_environment}" "${cpuset}" "${ranks}" "${started}" "${ended}" "${status}" "${duration}" "${peak}" "${log_sha}" "${artifact_sha}" >"${manifest}"
  printf 'STAGE3_ROW row=%s exit=%s manifest=%s\n' "${row}" "${status}" "${manifest}"
  if test "${required}" = 1 && test "${status}" -ne 0; then overall=1; fi
done <"${inventory}"
exit "${overall}"
