#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set -euo pipefail

if test "$#" -ne 4; then
  echo "usage: task19a_immersed_flow_dispatch.sh <hundun> <stl-writer> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
stl_writer=$2
mpiexec=$3
ranks=$4
run=("${mpiexec}" -n "${ranks}" "${hundun}")
process_grid="${ranks},1,1"
if test "${ranks}" -eq 4; then
  process_grid="2,2,1"
fi
work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-task19a-dispatch.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

case_path="${work_root}/case.json"
resolved_path="${work_root}/resolved.json"
second_path="${work_root}/resolved-second.json"
mkdir -p "${work_root}/geometry"
"${stl_writer}" "${work_root}/geometry/body.stl"

cat >"${case_path}" <<JSON
{"schema_version":3,"case":{"name":"task19a_dispatch"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"resources":{"expected_ranks":${ranks},"process_grid":[${process_grid}]},"mesh":{"cells":[12,12,12],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":1,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":1},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":false},"performance":{"enabled":false,"directory":"performance","warmup_steps":1,"measured_steps":1,"repetitions":1},"immersed_boundary":{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"outside"},"wall":{"velocity_m_per_s":[0.0,0.0,0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}},"les":{"model":"none"}}
JSON

validate_output=$("${run[@]}" "${case_path}" --validate)
test "${validate_output}" = "VALID"

"${run[@]}" "${case_path}" --print-resolved >"${resolved_path}"
test "$(LC_ALL=C tr -cd '\n' <"${resolved_path}" | wc -c)" -eq 1
grep -Fq '"schema_version":3' "${resolved_path}"
grep -Fq '"immersed_boundary":{"model":"local_flow_pattern_ghost_cell"' \
  "${resolved_path}"
grep -Fq '"les":{"model":"none"}' "${resolved_path}"

"${run[@]}" "${resolved_path}" --print-resolved >"${second_path}"
cmp "${resolved_path}" "${second_path}"

if ! "${run[@]}" "${case_path}" >"${work_root}/run.stdout" \
    2>"${work_root}/run.stderr"; then
  echo "Task 19A normal execution failed" >&2
  cat "${work_root}/run.stdout" >&2
  cat "${work_root}/run.stderr" >&2
  exit 1
fi
test ! -s "${work_root}/run.stderr"
grep -Fxq \
  "CASE name=task19a_dispatch ranks=${ranks} cells=12x12x12 immersed_boundary=local_flow_pattern_ghost_cell density_model=constant les=none" \
  "${work_root}/run.stdout"
grep -Eq '^STEP step=1 .*correctors=2 .*force_operator=' \
  "${work_root}/run.stdout"
grep -Eq '^FINISHED step=1 ' "${work_root}/run.stdout"

checkpoint_path="${work_root}/checkpoints/step00000000000000000001"
test -f "${checkpoint_path}/COMPLETED"
test -f "${checkpoint_path}/manifest.v3.bin"
test -f "${checkpoint_path}/rank-000000.v3.bin"
last_rank=$(printf '%06d' "$((ranks - 1))")
test -f "${checkpoint_path}/rank-${last_rank}.v3.bin"

restart_path="${work_root}/restart.json"
sed \
  -e 's/"steps":1/"steps":2/' \
  -e 's#"read":false,"write_directory":"checkpoints"#"read":true,"read_directory":"checkpoints/step00000000000000000001","write_directory":"resumed-checkpoints"#' \
  "${case_path}" >"${restart_path}"
if ! "${run[@]}" "${restart_path}" >"${work_root}/restart.stdout" \
    2>"${work_root}/restart.stderr"; then
  echo "Task 17A restart execution failed" >&2
  cat "${work_root}/restart.stdout" >&2
  cat "${work_root}/restart.stderr" >&2
  exit 1
fi
test ! -s "${work_root}/restart.stderr"
test "$(grep -c '^STEP ' "${work_root}/restart.stdout")" -eq 1
grep -Eq '^STEP step=2 .*correctors=2 .*force_operator=' \
  "${work_root}/restart.stdout"
grep -Eq '^FINISHED step=2 ' "${work_root}/restart.stdout"
test -f "${work_root}/resumed-checkpoints/step00000000000000000002/COMPLETED"
test -f "${work_root}/resumed-checkpoints/step00000000000000000002/rank-${last_rank}.v3.bin"

wale_body_path="${work_root}/wale-body.json"
cat >"${wale_body_path}" <<JSON
{"schema_version":3,"case":{"name":"task19b_wale_body"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"resources":{"expected_ranks":${ranks},"process_grid":[${process_grid}]},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":1,"initial_dt_s":0.0001,"min_dt_s":0.0000125,"max_dt_s":0.0001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"wale-checkpoints","write_interval":10},"diagnostics":{"directory":"wale-diagnostics","write_interval":1,"write_mesh":false},"performance":{"enabled":false,"directory":"wale-performance","warmup_steps":1,"measured_steps":1,"repetitions":1},"immersed_boundary":{"model":"none"},"les":{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7}}
JSON
if ! "${run[@]}" "${wale_body_path}" >"${work_root}/wale-body.stdout" \
    2>"${work_root}/wale-body.stderr"; then
  echo "Task 19B body-fitted WALE execution failed" >&2
  cat "${work_root}/wale-body.stdout" >&2
  cat "${work_root}/wale-body.stderr" >&2
  exit 1
fi
test ! -s "${work_root}/wale-body.stderr"
grep -Fxq \
  "CASE name=task19b_wale_body ranks=${ranks} cells=8x8x4 immersed_boundary=none density_model=constant les=wale" \
  "${work_root}/wale-body.stdout"
grep -Eq '^STEP step=1 .*correctors=2 .*wale_identity=[1-9][0-9]* .*wale_zero_count=' \
  "${work_root}/wale-body.stdout"
grep -Eq '^FINISHED step=1 ' "${work_root}/wale-body.stdout"

if test "${ranks}" -eq 1; then
  wale_path="${work_root}/wale.json"
  sed \
    's/"les":{"model":"none"}/"les":{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7}/' \
    "${case_path}" >"${wale_path}"
  if "${run[@]}" "${wale_path}" >"${work_root}/wale.stdout" \
      2>"${work_root}/wale.stderr"; then
    echo "Task 19B unexpectedly accepted combined IBM+WALE" >&2
    exit 1
  fi
  grep -Fxq 'immersed-flow driver does not yet support combined IBM+WALE' \
    "${work_root}/wale.stderr"
fi
