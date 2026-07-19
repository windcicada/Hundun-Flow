#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if test "$#" -ne 1; then
  echo "usage: stage2_task2_dispatch.sh <hundun>" >&2
  exit 2
fi

hundun=$1
work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-stage2-dispatch.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

case_path="${work_root}/case.json"
resolved_path="${work_root}/resolved.json"
second_path="${work_root}/resolved-second.json"

cat >"${case_path}" <<'JSON'
{"schema_version":2,"case":{"name":"dispatch_only"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":10,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"no_slip_wall"},{"patch":"x_max","type":"no_slip_wall"},{"patch":"y_min","type":"symmetry"},{"patch":"y_max","type":"symmetry"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":10},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":true},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}}
JSON

validate_output=$("${hundun}" "${case_path}" --validate)
test "${validate_output}" = "VALID"

"${hundun}" "${case_path}" --print-resolved >"${resolved_path}"
test "$(LC_ALL=C tr -cd '\n' <"${resolved_path}" | wc -c)" -eq 1
grep -Fq '"schema_version":2' "${resolved_path}"
grep -Fq '"resources":{}' "${resolved_path}"
if grep -Fq 'HUNDUN-FLOW 0.0.0-stage1' "${resolved_path}"; then
  echo "schema-v2 resolved output printed the Stage 1 banner" >&2
  exit 1
fi
"${hundun}" "${resolved_path}" --print-resolved >"${second_path}"
cmp "${resolved_path}" "${second_path}"

if "${hundun}" "${case_path}" >"${work_root}/normal.stdout" \
    2>"${work_root}/normal.stderr"; then
  echo "schema-v2 normal execution unexpectedly succeeded" >&2
  exit 1
fi
test ! -s "${work_root}/normal.stdout"
test "$(cat "${work_root}/normal.stderr")" = \
  "Stage 2 variable-density flow driver is not implemented before Task 24"
