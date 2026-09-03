#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set -euo pipefail

if test "$#" -ne 3; then
  echo "usage: run_task25_performance_release.sh <hundun> <mpiexec> <output>" >&2
  exit 2
fi

hundun=$1
mpiexec=$2
output=$3
mkdir -p -- "${output}"
evidence_writer="$(dirname -- "${hundun}")/task25_performance_evidence"
test -x "${evidence_writer}"
evidence_manifest="${output}/performance-evidence-input.tsv"
: >"${evidence_manifest}"

write_case() {
  case_path=$1
  ranks=$2
  grid=$3
  cells=$4
  case_id=$5
  repetitions=$6
  interval=$7
  cat >"${case_path}" <<JSON
{"schema_version":2,"case":{"name":"${case_id}"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"resources":{"expected_ranks":${ranks},"process_grid":[${grid}]},"mesh":{"cells":[${cells}],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":25,"initial_dt_s":0.0001,"min_dt_s":0.0001,"max_dt_s":0.0001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.01,"inlet_consistency_rtol":1e-12},"scalars":[{"name":"alpha","diffusivity_m2_per_s":0.01}],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoint-${case_id}","write_interval":${interval}},"diagnostics":{"directory":"diagnostics-${case_id}","write_interval":${interval},"write_mesh":false},"performance":{"enabled":true,"directory":"artifact-${case_id}","warmup_steps":5,"measured_steps":20,"repetitions":${repetitions}}}
JSON
}

manifest_entries=
for ranks in 1 2 4; do
  case "${ranks}" in
    1)
      grid='1,1,1'
      weak_cells='32,32,32'
      ;;
    2)
      grid='2,1,1'
      weak_cells='64,32,32'
      ;;
    4)
      grid='2,2,1'
      weak_cells='64,64,32'
      ;;
  esac
  for scaling in strong weak; do
    case_id="task25-${scaling}-${ranks}rank"
    case_path="${output}/${case_id}.json"
    if test "${scaling}" = strong; then
      cells='64,64,64'
    else
      cells=${weak_cells}
    fi
    write_case "${case_path}" "${ranks}" "${grid}" "${cells}" \
      "${case_id}" 5 1000
    command_text="${mpiexec} -n ${ranks} ${hundun} ${case_path}"
    command_sha=$(printf '%s' "${command_text}" | sha256sum | cut -d' ' -f1)
    case_sha=$(sha256sum "${case_path}" | cut -d' ' -f1)
    time_file="${output}/${case_id}.time"
    if test -x /usr/bin/time; then
      /usr/bin/time -f '%e %M' -o "${time_file}" \
        "${mpiexec}" -n "${ranks}" "${hundun}" "${case_path}" \
        >"${output}/${case_id}.stdout" 2>"${output}/${case_id}.stderr"
      read -r wall_seconds max_rss_kb <"${time_file}"
    else
      wall_seconds=null
      max_rss_kb=null
      "${mpiexec}" -n "${ranks}" "${hundun}" "${case_path}" \
        >"${output}/${case_id}.stdout" 2>"${output}/${case_id}.stderr"
    fi
    test ! -s "${output}/${case_id}.stdout"
    test ! -s "${output}/${case_id}.stderr"
    artifact="${output}/artifact-${case_id}/performance.v1.json"
    artifact_sha=$(sha256sum "${artifact}" | cut -d' ' -f1)
    entry="{\"case_id\":\"${case_id}\",\"command_sha256\":\"${command_sha}\",\"case_sha256\":\"${case_sha}\",\"artifact_sha256\":\"${artifact_sha}\"}"
    if test -n "${manifest_entries}"; then
      manifest_entries="${manifest_entries},${entry}"
    else
      manifest_entries=${entry}
    fi
    evidence_input="${output}/artifact-${case_id}/performance-evidence-input.v1.json"
    test -s "${evidence_input}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${case_id}" "${command_sha}" "${case_sha}" "${artifact_sha}" \
      "${case_path}" "${artifact}" "${evidence_input}" \
      "${wall_seconds}" "${max_rss_kb}" >>"${evidence_manifest}"
  done

  io_case_id="task25-io-${ranks}rank"
  io_case_path="${output}/${io_case_id}.json"
  write_case "${io_case_path}" "${ranks}" "${grid}" "${weak_cells}" \
    "${io_case_id}" 1 2
  sed -i \
    -e 's/"steps":25/"steps":3/' \
    -e 's/"warmup_steps":5/"warmup_steps":1/' \
    -e 's/"measured_steps":20/"measured_steps":2/' \
    "${io_case_path}"
  command_text="${mpiexec} -n ${ranks} ${hundun} ${io_case_path}"
  command_sha=$(printf '%s' "${command_text}" | sha256sum | cut -d' ' -f1)
  case_sha=$(sha256sum "${io_case_path}" | cut -d' ' -f1)
  time_file="${output}/${io_case_id}.time"
  if test -x /usr/bin/time; then
    /usr/bin/time -f '%e %M' -o "${time_file}" \
      "${mpiexec}" -n "${ranks}" "${hundun}" "${io_case_path}" \
      >"${output}/${io_case_id}.stdout" 2>"${output}/${io_case_id}.stderr"
    read -r wall_seconds max_rss_kb <"${time_file}"
  else
    wall_seconds=null
    max_rss_kb=null
    "${mpiexec}" -n "${ranks}" "${hundun}" "${io_case_path}" \
      >"${output}/${io_case_id}.stdout" 2>"${output}/${io_case_id}.stderr"
  fi
  test ! -s "${output}/${io_case_id}.stdout"
  test ! -s "${output}/${io_case_id}.stderr"
  artifact="${output}/artifact-${io_case_id}/performance.v1.json"
  artifact_sha=$(sha256sum "${artifact}" | cut -d' ' -f1)
  entry="{\"case_id\":\"${io_case_id}\",\"command_sha256\":\"${command_sha}\",\"case_sha256\":\"${case_sha}\",\"artifact_sha256\":\"${artifact_sha}\"}"
  manifest_entries="${manifest_entries},${entry}"
  evidence_input="${output}/artifact-${io_case_id}/performance-evidence-input.v1.json"
  test -s "${evidence_input}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${io_case_id}" "${command_sha}" "${case_sha}" "${artifact_sha}" \
    "${io_case_path}" "${artifact}" "${evidence_input}" \
    "${wall_seconds}" "${max_rss_kb}" >>"${evidence_manifest}"
done

printf '%s\n' \
  "{\"schema_version\":1,\"collection_boundary\":\"manual-release\",\"entries\":[${manifest_entries}]}" \
  >"${output}/performance-manifest.v1.json"
"${evidence_writer}" "${evidence_manifest}" \
  "${output}/performance-evidence.v1.json"
