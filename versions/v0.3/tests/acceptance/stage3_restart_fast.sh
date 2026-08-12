#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if test "$#" -ne 4; then
  echo "usage: stage3_restart_fast.sh <hundun> <stl-writer> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
stl_writer=$2
mpiexec=$3
ranks=$4
if test "${ranks}" -ne 1 && test "${ranks}" -ne 2; then
  echo "stage3_restart_fast.sh supports only 1 or 2 ranks" >&2
  exit 2
fi

work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-stage3-restart-fast.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

mkdir -p "${work_root}/geometry"
"${stl_writer}" "${work_root}/geometry/body.stl"

write_case() {
  local path=$1
  local name=$2
  local density=$3
  local steps=$4
  local read_directory=$5
  local write_directory=$6
  local diagnostics_directory=$7
  local immersed=$8
  local wale=$9

  local physics='{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.01,"inlet_consistency_rtol":1e-12}'
  local scalars='[]'
  if test "${density}" = material; then
    scalars='[{"name":"alpha","diffusivity_m2_per_s":0.0}]'
  elif test "${density}" = ideal_gas; then
    physics='{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.01,"inlet_consistency_rtol":1e-12,"cp_J_per_kg_K":1000.0,"gas_constant_J_per_kg_K":287.0,"thermodynamic_pressure_pa":86100.0}'
    scalars='[{"name":"alpha","diffusivity_m2_per_s":0.0}]'
  fi

  local immersed_json='{"model":"none"}'
  if test "${immersed}" = yes; then
    immersed_json='{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"outside"},"wall":{"velocity_m_per_s":[0.0,0.0,0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}}'
  fi
  local les_json='{"model":"none"}'
  if test "${wale}" = yes; then
    les_json='{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7}'
  fi
  local restart_read='"read":false'
  if test "${read_directory}" != none; then
    restart_read="\"read\":true,\"read_directory\":\"${read_directory}\""
  fi
  local cells=8
  local origin='0.0,0.0,0.0'
  local length='1.0,1.0,1.0'
  if test "${immersed}" = yes && test "${density}" = constant; then
    cells=12
  elif test "${immersed}" = yes; then
    origin='-0.5,-0.5,-0.5'
    length='2.0,2.0,2.0'
  fi

  cat >"${path}" <<JSON
{"schema_version":3,"case":{"name":"${name}"},"simulation":{"type":"variable_density_flow","density_model":"${density}"},"resources":{"expected_ranks":${ranks},"process_grid":[${ranks},1,1]},"mesh":{"cells":[${cells},${cells},${cells}],"origin_m":[${origin}],"length_m":[${length}],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":${steps},"initial_dt_s":0.0001,"min_dt_s":0.0000125,"max_dt_s":0.0001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":${physics},"scalars":${scalars},"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{${restart_read},"write_directory":"${write_directory}","write_interval":1},"diagnostics":{"directory":"${diagnostics_directory}","write_interval":2,"write_mesh":false},"performance":{"enabled":false,"directory":"performance","warmup_steps":1,"measured_steps":1,"repetitions":1},"immersed_boundary":${immersed_json},"les":${les_json}}
JSON
}

run_case() {
  local case_path=$1
  local stdout_path=$2
  local stderr_path=$3
  if ! "${mpiexec}" -n "${ranks}" "${hundun}" "${case_path}" \
      >"${stdout_path}" 2>"${stderr_path}"; then
    sed -n '1,240p' "${stdout_path}" >&2
    sed -n '1,240p' "${stderr_path}" >&2
    return 1
  fi
  test ! -s "${stderr_path}"
}

check_profile() {
  local tag=$1
  local density=$2
  local immersed=$3
  local wale=$4
  local name="stage3_restart_fast_${tag}"
  local continuous_root="${tag}-continuous-checkpoints"
  local split_root="${tag}-split-checkpoints"

  write_case "${work_root}/${tag}-continuous.json" "${name}" "${density}" 2 \
    none "${continuous_root}" "${tag}-continuous-diagnostics" \
    "${immersed}" "${wale}"
  run_case "${work_root}/${tag}-continuous.json" \
    "${work_root}/${tag}-continuous.stdout" \
    "${work_root}/${tag}-continuous.stderr"
  test "$(grep -c '^STEP step=' "${work_root}/${tag}-continuous.stdout")" -eq 2
  grep -Eq '^FINISHED step=2 ' "${work_root}/${tag}-continuous.stdout"

  write_case "${work_root}/${tag}-first.json" "${name}" "${density}" 1 \
    none "${split_root}" "${tag}-first-diagnostics" "${immersed}" "${wale}"
  run_case "${work_root}/${tag}-first.json" \
    "${work_root}/${tag}-first.stdout" "${work_root}/${tag}-first.stderr"
  test -f \
    "${work_root}/${split_root}/step00000000000000000001/COMPLETED"

  write_case "${work_root}/${tag}-restart.json" "${name}" "${density}" 2 \
    "${split_root}/step00000000000000000001" "${split_root}" \
    "${tag}-restart-diagnostics" "${immersed}" "${wale}"
  run_case "${work_root}/${tag}-restart.json" \
    "${work_root}/${tag}-restart.stdout" \
    "${work_root}/${tag}-restart.stderr"
  test "$(grep -c '^STEP step=' "${work_root}/${tag}-restart.stdout")" -eq 1
  grep -Eq '^STEP step=2 .*correctors=2 ' "${work_root}/${tag}-restart.stdout"
  grep -Eq '^FINISHED step=2 ' "${work_root}/${tag}-restart.stdout"

  diff -r \
    "${work_root}/${continuous_root}/step00000000000000000002" \
    "${work_root}/${split_root}/step00000000000000000002"

  local expected_records=0
  if test "${immersed}" = yes; then
    expected_records=4
  fi
  if test "${wale}" = yes; then
    expected_records=$((expected_records + 1))
  fi
  local rank=0
  while test "${rank}" -lt "${ranks}"; do
    local rank_name
    rank_name=$(printf '%06d' "${rank}")
    local continuous_records="${work_root}/${tag}-continuous-diagnostics/diagnostics.v1.rank-${rank_name}.step-00000000000000000002.jsonl"
    local restart_records="${work_root}/${tag}-restart-diagnostics/diagnostics.v1.rank-${rank_name}.step-00000000000000000002.jsonl"
    test -f "${continuous_records}"
    test -f "${restart_records}"
    test "$(wc -l <"${continuous_records}")" -eq "${expected_records}"
    test "$(wc -l <"${restart_records}")" -eq "${expected_records}"
    diff -u \
      <(sed -E 's/.*"module_kind":"([^"]+)".*"state_fingerprint":\{"algorithm":"[^"]+","hex":"([^"]+)"\}.*/\1 \2/' \
          "${continuous_records}") \
      <(sed -E 's/.*"module_kind":"([^"]+)".*"state_fingerprint":\{"algorithm":"[^"]+","hex":"([^"]+)"\}.*/\1 \2/' \
          "${restart_records}")
    rank=$((rank + 1))
  done
}

check_profile profile_1_constant_ibm constant yes no
check_profile profile_2_constant_body_wale constant no yes
check_profile profile_3_constant_ibm_wale constant yes yes
check_profile profile_4_material_ibm material yes no
check_profile profile_5_material_body_wale material no yes
check_profile profile_6_material_ibm_wale material yes yes
check_profile profile_7_ideal_ibm ideal_gas yes no
check_profile profile_8_ideal_body_wale ideal_gas no yes
check_profile profile_9_ideal_ibm_wale ideal_gas yes yes
