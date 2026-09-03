#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set -euo pipefail

if test "$#" -ne 3; then
  echo "usage: task24_flow_models.sh <hundun> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
mpiexec=$2
ranks=$3
case "${ranks}" in
  1) warped_grid='1,1,1' ;;
  2) warped_grid='2,1,1' ;;
  4) warped_grid='2,2,1' ;;
  *)
    echo "unsupported Task 25 curved rank count: ${ranks}" >&2
    exit 2
    ;;
esac
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
for model in constant material ideal-gas; do
  source_name=${model}
  case_token=${model}
  if test "${model}" = "ideal-gas"; then
    source_name=ideal-gas
    case_token=ideal_gas
  fi
  sed -e "s/task24_${case_token}/task24_warped_${case_token}/" \
    -e "s/\"expected_ranks\":${ranks}/\"expected_ranks\":${ranks},\"process_grid\":[${warped_grid}]/" \
    -e 's/"mapping":"uniform_box"/"mapping":"analytic_warped_box","warp_amplitude":[0.02,-0.015,0.01]/' \
    -e "s/checkpoint-${case_token}/checkpoint-warped-${case_token}/" \
    -e "s/diagnostics-${case_token}/diagnostics-warped-${case_token}/" \
    -e 's/"write_mesh":true/"write_mesh":false/' \
    "${work_root}/${source_name}.json" \
    >"${work_root}/warped-${model}.json"
done
cat >"${work_root}/open-ideal.json" <<JSON
{"schema_version":2,"case":{"name":"task24_open_ideal"},"simulation":{"type":"variable_density_flow","density_model":"ideal_gas"},"resources":{"expected_ranks":${ranks}},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":1,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12,"cp_J_per_kg_K":1000.0,"gas_constant_J_per_kg_K":287.0,"thermodynamic_pressure_pa":86100.0},"scalars":[{"name":"alpha","diffusivity_m2_per_s":0.0}],"boundaries":[{"patch":"x_min","type":"velocity_inlet","velocity_m_per_s":[0.1,0.0,0.0],"thermal_authority":"temperature","temperature_K":300.0,"enthalpy_J_per_kg":300000.0,"density_kg_per_m3":1.0,"scalar_values":[{"name":"alpha","value":0.25}]},{"patch":"x_max","type":"pressure_outlet","pressure_perturbation_pa":0.0},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoint-open","write_interval":10},"diagnostics":{"directory":"diagnostics-open","write_interval":1,"write_mesh":false},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}}
JSON

for model in constant material ideal-gas; do
  if test "${model}" = constant; then
    run_model=(env
      "HUNDUN_DIAGNOSTIC_OBSERVER=${work_root}/observer-due"
      "${mpiexec}" -n "${ranks}" "${hundun}")
  else
    run_model=("${mpiexec}" -n "${ranks}" "${hundun}")
  fi
  if ! "${run_model[@]}" "${work_root}/${model}.json" \
      >"${work_root}/${model}.stdout" 2>"${work_root}/${model}.stderr"; then
    cat "${work_root}/${model}.stdout" >&2
    cat "${work_root}/${model}.stderr" >&2
    exit 1
  fi
  test ! -s "${work_root}/${model}.stderr"
  grep -Fxq 'HUNDUN-FLOW 0.1.0' "${work_root}/${model}.stdout"
  grep -Fq 'FINISHED step=1 time_s=0.001' "${work_root}/${model}.stdout"
done

rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  due_observer="${work_root}/observer-due.rank-${rank_name}.jsonl"
  test "$(wc -l <"${due_observer}")" -eq 1
  grep -Fq '"mode":"due"' "${due_observer}"
  grep -Fq '"allocation_equal":true,"runtime_halo_equal":true,"pressure_halo_equal":true,"fp64_equal":true,"solves_equal":true,"state_equal":true' \
    "${due_observer}"
  grep -Eq '"logical_bytes":[1-9][0-9]*}' "${due_observer}"
  rank=$((rank + 1))
done

for model in constant material ideal-gas; do
  "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/warped-${model}.json" \
    >"${work_root}/warped-${model}.stdout" \
    2>"${work_root}/warped-${model}.stderr"
  test ! -s "${work_root}/warped-${model}.stderr"
  grep -Fxq 'HUNDUN-FLOW 0.1.0' \
    "${work_root}/warped-${model}.stdout"
  grep -Fq 'FINISHED step=1 time_s=0.001' \
    "${work_root}/warped-${model}.stdout"
done

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

module_sequence() {
  sed -n 's/.*"module_id":"\([^"]*\)".*/\1/p' "$1"
}

require_solve_instances() {
  sed -n \
    '/"module_id":"hundun.linear.solve"/s/.*"instance_id":"\([^"]*\)".*/\1/p' \
    "$1" | diff -u - <(cat <<'EXPECTED'
momentum-x
momentum-y
momentum-z
pressure-corrector-1
pressure-corrector-2
EXPECTED
)
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

for records in "${constant_records}" "${material_records}" "${ideal_records}"; do
  test "$(grep -Fc '"module_id":"hundun.linear.solve"' "${records}")" -eq 5
  require_solve_instances "${records}"
done

require_module_once "${constant_records}" hundun.flow.constant_density_piso
require_module_once \
  "${material_records}" hundun.flow.fixed_step_material_density
require_module_once \
  "${ideal_records}" hundun.flow.fixed_step_material_density
require_module_once "${ideal_records}" flow.ideal-gas-closure

module_sequence "${constant_records}" | diff -u - <(cat <<'EXPECTED'
hundun.runtime.mpi_context
hundun.runtime.structured_decomposition
hundun.runtime.field_layout
hundun.runtime.halo
hundun.mesh.topology
hundun.mesh.geometry
hundun.execution.context
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.finite_volume.shared_flux
hundun.boundary.registry
hundun.flow.constant_density_piso
hundun.application.flow_driver
hundun.flow.bdf2-retry-controller
EXPECTED
)
module_sequence "${material_records}" | diff -u - <(cat <<'EXPECTED'
hundun.runtime.mpi_context
hundun.runtime.structured_decomposition
hundun.runtime.field_layout
hundun.runtime.halo
hundun.mesh.topology
hundun.mesh.geometry
hundun.execution.context
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.finite_volume.shared_flux
hundun.boundary.registry
hundun.application.flow_driver
hundun.flow.fixed_step_material_density
hundun.flow.bdf2-retry-controller
EXPECTED
)
module_sequence "${ideal_records}" | diff -u - <(cat <<'EXPECTED'
hundun.runtime.mpi_context
hundun.runtime.structured_decomposition
hundun.runtime.field_layout
hundun.runtime.halo
hundun.mesh.topology
hundun.mesh.geometry
hundun.execution.context
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.linear.solve
hundun.finite_volume.shared_flux
hundun.boundary.registry
hundun.application.flow_driver
hundun.flow.fixed_step_material_density
flow.ideal-gas-closure
hundun.flow.bdf2-retry-controller
EXPECTED
)

sed -e 's/"write_interval":1,"write_mesh":true/"write_interval":10,"write_mesh":false/' \
  -e 's/checkpoint-constant/checkpoint-null-session/' \
  -e 's/diagnostics-constant/diagnostics-null-session/' \
  "${work_root}/constant.json" >"${work_root}/null-session.json"
env "HUNDUN_DIAGNOSTIC_OBSERVER=${work_root}/observer-disabled" \
  "${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/null-session.json" \
  >"${work_root}/null-session.stdout" 2>"${work_root}/null-session.stderr"
test ! -s "${work_root}/null-session.stderr"
test ! -e "${work_root}/diagnostics-null-session"

sed -e 's/"write_interval":1,"write_mesh"/"write_interval":10,"write_mesh"/' \
  -e 's/checkpoint-constant/checkpoint-not-due/' \
  -e 's/diagnostics-constant/diagnostics-not-due/' \
  "${work_root}/constant.json" >"${work_root}/not-due.json"
env "HUNDUN_DIAGNOSTIC_OBSERVER=${work_root}/observer-not-due" \
  "${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/not-due.json" \
  >"${work_root}/not-due.stdout" 2>"${work_root}/not-due.stderr"
test ! -s "${work_root}/not-due.stderr"
test -f "${work_root}/diagnostics-not-due/meshdiag.v2.rank-000000.bin"
test ! -e \
  "${work_root}/diagnostics-not-due/diagnostics.v1.rank-000000.step-00000000000000000001.jsonl"

rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  for mode in disabled not-due; do
    observer_file="${work_root}/observer-${mode}.rank-${rank_name}.jsonl"
    test "$(wc -l <"${observer_file}")" -eq 1
    grep -Fq '"allocation_equal":true,"runtime_halo_equal":true,"pressure_halo_equal":true,"fp64_equal":true,"solves_equal":true,"state_equal":true,"files_equal":true,"logical_bytes":0}' \
      "${observer_file}"
  done
  grep -Fq '"mode":"disabled"' \
    "${work_root}/observer-disabled.rank-${rank_name}.jsonl"
  grep -Fq '"mode":"not_due"' \
    "${work_root}/observer-not-due.rank-${rank_name}.jsonl"
  rank=$((rank + 1))
done

sed -e 's/diagnostics-constant/diagnostics-sink-failure/' \
  "${work_root}/constant.json" >"${work_root}/sink-failure.json"
mkdir "${work_root}/diagnostics-sink-failure"
rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  mkdir "${work_root}/diagnostics-sink-failure/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl.tmp"
  rank=$((rank + 1))
done
if env "HUNDUN_DIAGNOSTIC_OBSERVER=${work_root}/observer-sink" \
    "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/sink-failure.json" \
    >"${work_root}/sink-failure.stdout" \
    2>"${work_root}/sink-failure.stderr"; then
  echo "diagnostic sink failure unexpectedly succeeded" >&2
  exit 1
fi
rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  sink_observer="${work_root}/observer-sink.rank-${rank_name}.jsonl"
  test "$(wc -l <"${sink_observer}")" -eq 1
  grep -Fq '"mode":"sink_failure"' "${sink_observer}"
  grep -Fq '"allocation_equal":true,"runtime_halo_equal":true,"pressure_halo_equal":true,"fp64_equal":true,"solves_equal":true,"state_equal":true,"files_equal":true,"logical_bytes":0}' \
    "${sink_observer}"
  obstruction="${work_root}/diagnostics-sink-failure/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl.tmp"
  final_record="${work_root}/diagnostics-sink-failure/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl"
  test -d "${obstruction}"
  test ! -e "${final_record}"
  rank=$((rank + 1))
done
test -d "${work_root}/diagnostics-sink-failure"

sed -e 's/diagnostics-constant/diagnostics-sink-file-failure/' \
  "${work_root}/constant.json" >"${work_root}/sink-file-failure.json"
mkdir "${work_root}/diagnostics-sink-file-failure"
rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  staging_file="${work_root}/diagnostics-sink-file-failure/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl.tmp"
  printf 'preexisting-staging-file-rank-%s\n' "${rank_name}" \
    >"${staging_file}"
  rank=$((rank + 1))
done
if env "HUNDUN_DIAGNOSTIC_OBSERVER=${work_root}/observer-sink-file" \
    "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/sink-file-failure.json" \
    >"${work_root}/sink-file-failure.stdout" \
    2>"${work_root}/sink-file-failure.stderr"; then
  echo "pre-existing diagnostic staging file unexpectedly succeeded" >&2
  exit 1
fi
rank=0
while test "${rank}" -lt "${ranks}"; do
  rank_name=$(printf '%06d' "${rank}")
  sink_observer="${work_root}/observer-sink-file.rank-${rank_name}.jsonl"
  test "$(wc -l <"${sink_observer}")" -eq 1
  grep -Fq '"mode":"sink_failure"' "${sink_observer}"
  grep -Fq '"allocation_equal":true,"runtime_halo_equal":true,"pressure_halo_equal":true,"fp64_equal":true,"solves_equal":true,"state_equal":true,"files_equal":true,"logical_bytes":0}' \
    "${sink_observer}"
  staging_file="${work_root}/diagnostics-sink-file-failure/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl.tmp"
  final_record="${work_root}/diagnostics-sink-file-failure/diagnostics.v1.rank-${rank_name}.step-00000000000000000001.jsonl"
  test -f "${staging_file}"
  expected_sha=$(printf 'preexisting-staging-file-rank-%s\n' "${rank_name}" |
    sha256sum | cut -d' ' -f1)
  test "$(sha256sum "${staging_file}" | cut -d' ' -f1)" = \
    "${expected_sha}"
  test ! -e "${final_record}"
  rank=$((rank + 1))
done

if "${mpiexec}" -n "${ranks}" sh -c \
    'exec 1>/dev/full; exec "$1" "$2"' _ "${hundun}" \
    "${work_root}/null-session.json" \
    2>"${work_root}/stdout-failure.stderr"; then
  echo "root output failure unexpectedly succeeded" >&2
  exit 1
fi
grep -Fq 'unable to write flow output' \
  "${work_root}/stdout-failure.stderr"

sed -e 's/"enabled":false/"enabled":true/' \
  -e 's/checkpoint-constant/checkpoint-performance-reject/' \
  -e 's/diagnostics-constant/diagnostics-performance-reject/' \
  "${work_root}/constant.json" >"${work_root}/performance-request.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/performance-request.json" \
    >"${work_root}/performance.stdout" \
    2>"${work_root}/performance.stderr"; then
  echo "mismatched Task 25 performance request unexpectedly succeeded" >&2
  exit 1
fi
if grep -Fq 'HUNDUN-FLOW 0.1.0' \
    "${work_root}/performance.stdout"; then
  echo "performance rejection printed the Stage 2 run banner" >&2
  exit 1
fi
grep -Fq \
  'performance time.steps must equal warmup_steps + measured_steps' \
  "${work_root}/performance.stderr"
test ! -e "${work_root}/checkpoint-performance-reject"
test ! -e "${work_root}/diagnostics-performance-reject"

if test -n "${HUNDUN_TASK25_OBSERVER_EVIDENCE_DIR:-}"; then
  mkdir -p "${HUNDUN_TASK25_OBSERVER_EVIDENCE_DIR}"
  for observer_file in "${work_root}"/observer-*.jsonl; do
    cp "${observer_file}" \
      "${HUNDUN_TASK25_OBSERVER_EVIDENCE_DIR}/${ranks}rank-$(basename "${observer_file}")"
  done
fi
