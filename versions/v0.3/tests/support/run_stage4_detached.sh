#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set -euo pipefail
test "$#" -ge 3 || { echo 'usage: run_stage4_detached.sh EVIDENCE_DIR BINARY ARGS...' >&2; exit 2; }
evidence=$1
binary=$2
shift 2
[[ "${evidence}" = /* && "${binary}" = /* ]] || { echo 'absolute paths required' >&2; exit 2; }
source_root=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
head=$(git -C "${source_root}" rev-parse HEAD)
test -z "$(git -C "${source_root}" status --porcelain=v1 --untracked-files=normal)" || { echo 'dirty source' >&2; exit 2; }
mkdir -p -- "${evidence}"
sha=$(sha256sum "${binary}" | awk '{print $1}')
printf 'candidate_head=%s\nbinary=%s\nbinary_sha256=%s\n' "${head}" "${binary}" "${sha}" >"${evidence}/identity"
setsid "${binary}" "$@" >"${evidence}/run.log" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "${pid}" >"${evidence}/pid"
printf 'STAGE4_DETACHED_STARTED pid=%s head=%s binary_sha256=%s\n' "${pid}" "${head}" "${sha}"
