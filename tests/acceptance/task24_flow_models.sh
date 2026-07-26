#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if test "$#" -ne 3; then
  echo "usage: task24_flow_models.sh <hundun> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
mpiexec=$2
ranks=$3
work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-task24-models.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

write_case() {
  model=$1
  physics=$2
  case_path=$3
  write_mesh=$4
  cat >"${case_path}" <<JSON
{"schema_version":2,"case":{"name":"task24_${model}"},"simulation":{"type":"variable_density_flow","density_model":"${model}"},"resources":{"expected_ranks":${ranks}},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":1,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":${physics},"scalars":[{"name":"alpha","diffusivity_m2_per_s":0.0}],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoint-${model}","write_interval":10},"diagnostics":{"directory":"diagnostics-${model}","write_interval":1,"write_mesh":${write_mesh}},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}}
JSON
}

common_physics='{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12}'
ideal_physics='{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12,"cp_J_per_kg_K":1000.0,"gas_constant_J_per_kg_K":287.0,"thermodynamic_pressure_pa":86100.0}'

write_case constant "${common_physics}" "${work_root}/constant.json" true
write_case material "${common_physics}" "${work_root}/material.json" false
write_case ideal_gas "${ideal_physics}" "${work_root}/ideal-gas.json" false
sed -e 's/task24_constant/task24_warped/' \
  -e 's/"mapping":"uniform_box"/"mapping":"analytic_warped_box","warp_amplitude":[0.01,-0.005,0.005]/' \
  -e 's/checkpoint-constant/checkpoint-warped/' \
  -e 's/diagnostics-constant/diagnostics-warped/' \
  -e 's/"write_mesh":true/"write_mesh":false/' \
  "${work_root}/constant.json" >"${work_root}/warped.json"
cat >"${work_root}/open-ideal.json" <<JSON
{"schema_version":2,"case":{"name":"task24_open_ideal"},"simulation":{"type":"variable_density_flow","density_model":"ideal_gas"},"resources":{"expected_ranks":${ranks}},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":1,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12,"cp_J_per_kg_K":1000.0,"gas_constant_J_per_kg_K":287.0,"thermodynamic_pressure_pa":86100.0},"scalars":[{"name":"alpha","diffusivity_m2_per_s":0.0}],"boundaries":[{"patch":"x_min","type":"velocity_inlet","velocity_m_per_s":[0.1,0.0,0.0],"thermal_authority":"temperature","temperature_K":300.0,"enthalpy_J_per_kg":300000.0,"density_kg_per_m3":1.0,"scalar_values":[{"name":"alpha","value":0.25}]},{"patch":"x_max","type":"pressure_outlet","pressure_perturbation_pa":0.0},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoint-open","write_interval":10},"diagnostics":{"directory":"diagnostics-open","write_interval":1,"write_mesh":false},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}}
JSON

for model in constant material ideal-gas; do
  "${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/${model}.json" \
    >"${work_root}/${model}.stdout" 2>"${work_root}/${model}.stderr"
  test ! -s "${work_root}/${model}.stderr"
  grep -Fxq 'HUNDUN-FLOW 0.0.0-stage2' "${work_root}/${model}.stdout"
  grep -Fq 'FINISHED step=1 time_s=0.001' "${work_root}/${model}.stdout"
done

if "${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/warped.json" \
    >"${work_root}/warped.stdout" 2>"${work_root}/warped.stderr"; then
  echo "pre-Task25 warped flow unexpectedly advanced" >&2
  exit 1
fi
if grep -Fq 'FINISHED ' "${work_root}/warped.stdout"; then
  echo "rejected warped flow printed FINISHED" >&2
  exit 1
fi
grep -Fq 'Task 18 fixed-step flow supports uniform geometry only' \
  "${work_root}/warped.stderr"

if ! "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/open-ideal.json" >"${work_root}/open.stdout" \
    2>"${work_root}/open.stderr"; then
  cat "${work_root}/open.stdout" >&2
  cat "${work_root}/open.stderr" >&2
  exit 1
fi
test ! -s "${work_root}/open.stderr"
grep -Fq 'FINISHED step=1 time_s=0.001' "${work_root}/open.stdout"

test -f \
  "${work_root}/diagnostics-constant/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"
test -f \
  "${work_root}/diagnostics-material/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"
test -f \
  "${work_root}/diagnostics-ideal_gas/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"

rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  test -f \
    "${work_root}/diagnostics-constant/meshdiag.v2.rank-${rank_name}.bin"
  test ! -e \
    "${work_root}/diagnostics-material/meshdiag.v2.rank-${rank_name}.bin"
  rank=$((rank + 1))
done

require_module_once() {
  file=$1
  module=$2
  test "$(grep -Fc "\"module_id\":\"${module}\"" "${file}")" -eq 1
}

constant_records="${work_root}/diagnostics-constant/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"
material_records="${work_root}/diagnostics-material/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"
ideal_records="${work_root}/diagnostics-ideal_gas/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"

for records in "${constant_records}" "${material_records}" "${ideal_records}"; do
  for module in \
    hundun.runtime.mpi_context \
    hundun.runtime.structured_decomposition \
    hundun.runtime.field_layout \
    hundun.runtime.halo \
    hundun.mesh.topology \
    hundun.mesh.geometry \
    hundun.execution.context \
    hundun.finite_volume.shared_flux \
    hundun.boundary.registry \
    hundun.application.flow_driver \
    hundun.flow.bdf2-retry-controller; do
    require_module_once "${records}" "${module}"
  done
done

require_module_once "${constant_records}" hundun.flow.constant_density_piso
require_module_once \
  "${material_records}" hundun.flow.fixed_step_material_density
require_module_once \
  "${ideal_records}" hundun.flow.fixed_step_material_density
require_module_once "${ideal_records}" flow.ideal-gas-closure

sed -e 's/"enabled":false/"enabled":true/' \
  -e 's/checkpoint-constant/checkpoint-performance-reject/' \
  -e 's/diagnostics-constant/diagnostics-performance-reject/' \
  "${work_root}/constant.json" >"${work_root}/performance-request.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/performance-request.json" \
    >"${work_root}/performance.stdout" \
    2>"${work_root}/performance.stderr"; then
  echo "unsupported Task 25 performance request unexpectedly succeeded" >&2
  exit 1
fi
if grep -Fq 'HUNDUN-FLOW 0.0.0-stage2' \
    "${work_root}/performance.stdout"; then
  echo "performance rejection printed the Stage 2 run banner" >&2
  exit 1
fi
grep -Fq \
  'Stage 2 performance artifacts are unavailable before Task 25' \
  "${work_root}/performance.stderr"
test ! -e "${work_root}/checkpoint-performance-reject"
test ! -e "${work_root}/diagnostics-performance-reject"
