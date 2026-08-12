#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if test "$#" -ne 3; then
  echo "usage: task24_restart_continuation.sh <hundun> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
mpiexec=$2
ranks=$3
work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-task24-restart.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

cat >"${work_root}/continuous.json" <<JSON
{"schema_version":2,"case":{"name":"task24_restart"},"simulation":{"type":"variable_density_flow","density_model":"ideal_gas"},"resources":{"expected_ranks":${ranks}},"mesh":{"cells":[8,8,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":2,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12,"cp_J_per_kg_K":1000.0,"gas_constant_J_per_kg_K":287.0,"thermodynamic_pressure_pa":86100.0},"scalars":[{"name":"alpha","diffusivity_m2_per_s":0.0}],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoint","write_interval":1},"diagnostics":{"directory":"diagnostics","write_interval":2,"write_mesh":false},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5}}
JSON

"${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/continuous.json" \
  >"${work_root}/continuous.stdout" 2>"${work_root}/continuous.stderr"
test ! -s "${work_root}/continuous.stderr"
test -f \
  "${work_root}/checkpoint/step00000000000000000001/COMPLETED"
test -f \
  "${work_root}/checkpoint/step00000000000000000002/COMPLETED"
cp -a "${work_root}/checkpoint/step00000000000000000002" \
  "${work_root}/continuous-step2"

sed \
  -e 's|"read":false|"read":true,"read_directory":"missing-checkpoint"|' \
  -e 's|"write_directory":"checkpoint"|"write_directory":"checkpoint-missing-read"|' \
  -e 's|"directory":"diagnostics"|"directory":"diagnostics-missing-read"|' \
  "${work_root}/continuous.json" >"${work_root}/missing-read.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/missing-read.json" \
    >"${work_root}/missing-read.stdout" \
    2>"${work_root}/missing-read.stderr"; then
  echo "missing checkpoint read unexpectedly succeeded" >&2
  exit 1
fi
grep -Fq 'flow checkpoint read failed' \
  "${work_root}/missing-read.stderr"
test ! -e "${work_root}/checkpoint-missing-read"

printf 'not a directory\n' >"${work_root}/blocked-checkpoint"
sed \
  -e 's|"write_directory":"checkpoint"|"write_directory":"blocked-checkpoint"|' \
  -e 's|"directory":"diagnostics"|"directory":"diagnostics-blocked-write"|' \
  "${work_root}/continuous.json" >"${work_root}/blocked-write.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/blocked-write.json" \
    >"${work_root}/blocked-write.stdout" \
    2>"${work_root}/blocked-write.stderr"; then
  echo "blocked checkpoint write unexpectedly succeeded" >&2
  exit 1
fi
test ! -e "${work_root}/diagnostics-blocked-write"
if grep -Fq 'FINISHED ' "${work_root}/blocked-write.stdout"; then
  echo "failed checkpoint write printed FINISHED" >&2
  exit 1
fi

rm -rf -- "${work_root}/checkpoint/step00000000000000000002"

sed \
  's|"read":false|"read":true,"read_directory":"checkpoint/step00000000000000000001"|' \
  "${work_root}/continuous.json" >"${work_root}/restart.json"
"${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/restart.json" \
  >"${work_root}/restart.stdout" 2>"${work_root}/restart.stderr"
test ! -s "${work_root}/restart.stderr"
test "$(grep -c '^STEP ' "${work_root}/restart.stdout")" -eq 1
grep -Fq 'STEP 2 ' "${work_root}/restart.stdout"
grep -Fq 'FINISHED step=2 ' "${work_root}/restart.stdout"
diff -r "${work_root}/continuous-step2" \
  "${work_root}/checkpoint/step00000000000000000002"

records="${work_root}/diagnostics/diagnostics.v1.rank-000000.step-00000000000000000002.jsonl"
test -f "${records}"
test "$(grep -Fc '"module_id":"checkpoint-v2"' "${records}")" -eq 1
test "$(tail -n 1 "${records}" | grep -Fc '"module_id":"checkpoint-v2"')" -eq 1
