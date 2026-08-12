#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if test "$#" -ne 4; then
  echo "usage: stage3_flow_models_fast.sh <hundun> <stl-writer> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
stl_writer=$2
mpiexec=$3
ranks=$4
if test "${ranks}" -ne 1 && test "${ranks}" -ne 2; then
  echo "stage3_flow_models_fast.sh supports only 1 or 2 ranks" >&2
  exit 2
fi

work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-stage3-flow-models-fast.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

mkdir -p "${work_root}/geometry"
"${stl_writer}" "${work_root}/geometry/body.stl"

write_case() {
  local profile=$1
  local density=$2
  local immersed=$3
  local wale=$4
  local physics='{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.01,"inlet_consistency_rtol":1e-12}'
  local scalars='[]'
  if test "${density}" = material; then
    scalars='[{"name":"alpha","diffusivity_m2_per_s":0.0}]'
  elif test "${density}" = ideal_gas; then
    physics='{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.01,"inlet_consistency_rtol":1e-12,"cp_J_per_kg_K":1000.0,"gas_constant_J_per_kg_K":287.0,"thermodynamic_pressure_pa":86100.0}'
    scalars='[{"name":"alpha","diffusivity_m2_per_s":0.0}]'
  fi

  local immersed_json='{"model":"none"}'
  local cells=8
  local origin='0.0,0.0,0.0'
  local length='1.0,1.0,1.0'
  if test "${immersed}" = yes; then
    immersed_json='{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"outside"},"wall":{"velocity_m_per_s":[0.0,0.0,0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}}'
    if test "${density}" = constant; then
      cells=12
    else
      origin='-0.5,-0.5,-0.5'
      length='2.0,2.0,2.0'
    fi
  fi
  local les_json='{"model":"none"}'
  if test "${wale}" = yes; then
    les_json='{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7}'
  fi

  cat >"${work_root}/profile-${profile}.json" <<JSON
{"schema_version":3,"case":{"name":"stage3_profile_${profile}_flow_fast"},"simulation":{"type":"variable_density_flow","density_model":"${density}"},"resources":{"expected_ranks":${ranks},"process_grid":[${ranks},1,1]},"mesh":{"cells":[${cells},${cells},${cells}],"origin_m":[${origin}],"length_m":[${length}],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":1,"initial_dt_s":0.0001,"min_dt_s":0.0000125,"max_dt_s":0.0001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":${physics},"scalars":${scalars},"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"profile-${profile}-checkpoints","write_interval":1},"diagnostics":{"directory":"profile-${profile}-diagnostics","write_interval":1,"write_mesh":false},"performance":{"enabled":false,"directory":"performance","warmup_steps":1,"measured_steps":1,"repetitions":1},"immersed_boundary":${immersed_json},"les":${les_json}}
JSON
}

require_diagnostics() {
  local profile=$1
  local immersed=$2
  local wale=$3
  local rank=0
  local expected=0
  if test "${immersed}" = yes; then
    expected=4
  fi
  if test "${wale}" = yes; then
    expected=$((expected + 1))
  fi
  while test "${rank}" -lt "${ranks}"; do
    local rank_name
    rank_name=$(printf '%06d' "${rank}")
    local records="${work_root}/profile-${profile}-diagnostics/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl"
    test -f "${records}"
    test "$(wc -l <"${records}")" -eq "${expected}"
    for kind in immersed_surface ghost_stencil local_flow_pattern wall_force; do
      local count=0
      if test "${immersed}" = yes; then
        count=1
      fi
      test "$(grep -Fc "\"module_kind\":\"${kind}\"" "${records}")" -eq "${count}"
    done
    local les_count=0
    if test "${wale}" = yes; then
      les_count=1
    fi
    test "$(grep -Fc '"module_kind":"les"' "${records}")" -eq "${les_count}"
    rank=$((rank + 1))
  done
}

run_profile() {
  local profile=$1
  local density=$2
  local immersed=$3
  local wale=$4
  write_case "${profile}" "${density}" "${immersed}" "${wale}"
  if ! "${mpiexec}" -n "${ranks}" "${hundun}" \
      "${work_root}/profile-${profile}.json" \
      >"${work_root}/profile-${profile}.stdout" \
      2>"${work_root}/profile-${profile}.stderr"; then
    echo "Stage 3 profile ${profile} fast execution failed" >&2
    sed -n '1,240p' "${work_root}/profile-${profile}.stdout" >&2
    sed -n '1,240p' "${work_root}/profile-${profile}.stderr" >&2
    exit 1
  fi
  test ! -s "${work_root}/profile-${profile}.stderr"
  grep -Fxq \
    "CASE name=stage3_profile_${profile}_flow_fast ranks=${ranks} cells=$(test "${immersed}" = yes && test "${density}" = constant && printf '12x12x12' || printf '8x8x8') immersed_boundary=$(test "${immersed}" = yes && printf 'local_flow_pattern_ghost_cell' || printf 'none') density_model=${density} les=$(test "${wale}" = yes && printf 'wale' || printf 'none')" \
    "${work_root}/profile-${profile}.stdout"
  grep -Eq '^STEP step=1 .*correctors=2 ' "${work_root}/profile-${profile}.stdout"
  if test "${immersed}" = yes; then
    grep -Eq '^STEP step=1 .*force_operator=[^ ]+ .*force_budget_reaction=[^ ]+ .*force_surface_traction=[^ ]+ .*force_consistency=[^ ]+' \
      "${work_root}/profile-${profile}.stdout"
  elif grep -Eq 'force_' "${work_root}/profile-${profile}.stdout"; then
    echo "Stage 3 profile ${profile} claimed absent IBM force" >&2
    exit 1
  fi
  if test "${wale}" = yes; then
    grep -Eq '^STEP step=1 .*wale_identity=[1-9][0-9]* .*wale_nu_t_min=[^ ]+ .*wale_nu_t_max=[^ ]+' \
      "${work_root}/profile-${profile}.stdout"
  elif grep -Eq 'wale_' "${work_root}/profile-${profile}.stdout"; then
    echo "Stage 3 profile ${profile} claimed absent WALE" >&2
    exit 1
  fi
  grep -Eq '^FINISHED step=1 ' "${work_root}/profile-${profile}.stdout"
  test -f "${work_root}/profile-${profile}-checkpoints/step00000000000000000001/COMPLETED"
  require_diagnostics "${profile}" "${immersed}" "${wale}"
}

run_profile 1 constant yes no
run_profile 2 constant no yes
run_profile 3 constant yes yes
run_profile 4 material yes no
run_profile 5 material no yes
run_profile 6 material yes yes
run_profile 7 ideal_gas yes no
run_profile 8 ideal_gas no yes
run_profile 9 ideal_gas yes yes
