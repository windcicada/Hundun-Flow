#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set -euo pipefail

if test "$#" -ne 2; then
  echo "usage: stage2_task2_dispatch.sh <hundun> <mpiexec>" >&2
  exit 2
fi

hundun=$1
mpiexec=$2
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
if grep -Fq 'HUNDUN-FLOW 0.1.0' "${resolved_path}"; then
  echo "schema-v2 resolved output printed the Stage 1 banner" >&2
  exit 1
fi
"${hundun}" "${resolved_path}" --print-resolved >"${second_path}"
cmp "${resolved_path}" "${second_path}"

if ! "${hundun}" "${case_path}" >"${work_root}/normal.stdout" \
    2>"${work_root}/normal.stderr"; then
  echo "schema-v2 normal execution failed" >&2
  cat "${work_root}/normal.stdout" >&2
  cat "${work_root}/normal.stderr" >&2
  exit 1
fi
test ! -s "${work_root}/normal.stderr"
grep -Fxq 'HUNDUN-FLOW 0.1.0' "${work_root}/normal.stdout"
grep -Fxq \
  'CASE name=dispatch_only ranks=1 cells=8x8x4 density_model=constant' \
  "${work_root}/normal.stdout"
test "$(grep -c '^STEP ' "${work_root}/normal.stdout")" -eq 10
grep -Fq 'FINISHED step=10 ' "${work_root}/normal.stdout"
test -f "${work_root}/diagnostics/meshdiag.v2.rank-000000.bin"
test -f \
  "${work_root}/diagnostics/diagnostics.v1.rank-000000.step-00000000000000000010.jsonl"
test -f \
  "${work_root}/checkpoints/step00000000000000000010/COMPLETED"

mkdir -p "${work_root}/rank-zero" "${work_root}/rank-one"
cat >"${work_root}/rank-zero/invalid-v1.json" <<'JSON'
{"schema_version":1,"case":{"name":"invalid_cells"},"resources":{"expected_ranks":1},"mesh":{"cells":[0,4,2],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":2},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}}
JSON
cat >"${work_root}/rank-one/competing-v1.json" <<'JSON'
{"schema_version":1,"case":{"name":"non_authoritative"},"resources":{"expected_ranks":2},"mesh":{"cells":[8,4,2],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":2},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}}
JSON
if "${mpiexec}" -n 1 "${hundun}" \
      "${work_root}/rank-zero/invalid-v1.json" --validate : \
      -n 1 "${hundun}" \
      "${work_root}/rank-one/competing-v1.json" --validate \
      >"${work_root}/invalid-v1.stdout" \
      2>"${work_root}/invalid-v1.stderr"; then
  echo "invalid authoritative schema-v1 case unexpectedly succeeded" >&2
  exit 1
fi
if grep -Fxq 'VALID' "${work_root}/invalid-v1.stdout"; then
  echo "invalid authoritative schema-v1 case printed VALID" >&2
  exit 1
fi
if ! grep -Fxq '/mesh/cells/0: expected an integer of at least 1, got 0' \
    "${work_root}/invalid-v1.stderr"; then
  echo "schema-v1 invalid-case failure was not stable" >&2
  cat "${work_root}/invalid-v1.stdout" >&2
  cat "${work_root}/invalid-v1.stderr" >&2
  exit 1
fi

cat >"${work_root}/rank-zero/rank-mismatch-v1.json" <<'JSON'
{"schema_version":1,"case":{"name":"rank_mismatch"},"resources":{"expected_ranks":1},"mesh":{"cells":[8,4,2],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":2},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}}
JSON
if "${mpiexec}" -n 1 "${hundun}" \
      "${work_root}/rank-zero/rank-mismatch-v1.json" --validate : \
      -n 1 "${hundun}" \
      "${work_root}/rank-one/competing-v1.json" --validate \
      >"${work_root}/rank-mismatch-v1.stdout" \
      2>"${work_root}/rank-mismatch-v1.stderr"; then
  echo "schema-v1 expected-ranks mismatch unexpectedly succeeded" >&2
  exit 1
fi
if grep -Fxq 'VALID' "${work_root}/rank-mismatch-v1.stdout"; then
  echo "schema-v1 expected-ranks mismatch printed VALID" >&2
  exit 1
fi
if ! grep -Fxq 'expected MPI rank count 1, got 2' \
    "${work_root}/rank-mismatch-v1.stderr"; then
  echo "schema-v1 rank-mismatch failure was not stable" >&2
  cat "${work_root}/rank-mismatch-v1.stdout" >&2
  cat "${work_root}/rank-mismatch-v1.stderr" >&2
  exit 1
fi
