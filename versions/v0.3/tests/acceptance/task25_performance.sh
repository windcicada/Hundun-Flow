#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if test "$#" -ne 3; then
  echo "usage: task25_performance.sh <hundun> <mpiexec> <ranks>" >&2
  exit 2
fi

hundun=$1
mpiexec=$2
ranks=$3
case "${ranks}" in
  1)
    cells='8,8,4'
    grid='1,1,1'
    portable_allocated=67584
    portable_peak=87064
    portable_halo_bytes=2809856
    portable_halo_messages=3328
    positive_allocated=33792
    positive_peak=87064
    positive_halo_bytes=1404928
    positive_halo_messages=1664
    positive_checkpoint_bytes=86016
    positive_diagnostic_bytes=103097
    peak_ratio_bits=4075418000000000
    ;;
  2)
    cells='16,8,4'
    grid='2,1,1'
    portable_allocated=135168
    portable_peak=177200
    portable_halo_bytes=5644288
    portable_halo_messages=6704
    positive_allocated=67584
    positive_peak=177200
    positive_halo_bytes=2822144
    positive_halo_messages=3352
    positive_checkpoint_bytes=172032
    positive_diagnostic_bytes=206203
    peak_ratio_bits=4075a18000000000
    ;;
  4)
    cells='16,16,4'
    grid='2,2,1'
    portable_allocated=270336
    portable_peak=360544
    portable_halo_bytes=11337728
    portable_halo_messages=13504
    positive_allocated=135168
    positive_peak=360544
    positive_halo_bytes=5668864
    positive_halo_messages=6752
    positive_checkpoint_bytes=344064
    positive_diagnostic_bytes=412427
    peak_ratio_bits=4076018000000000
    ;;
  *)
    echo "unsupported Task 25 rank count: ${ranks}" >&2
    exit 2
    ;;
esac

work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-task25-performance.XXXXXX")
cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

write_case() {
  destination=$1
  repetitions=$2
  diagnostic_interval=$3
  checkpoint_interval=$4
  performance_directory=$5
  cat >"${destination}" <<JSON
{"schema_version":2,"case":{"name":"task25_performance"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"resources":{"expected_ranks":${ranks},"process_grid":[${grid}]},"mesh":{"cells":[${cells}],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":3,"initial_dt_s":0.0001,"min_dt_s":0.0001,"max_dt_s":0.0001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.01,"inlet_consistency_rtol":1e-12},"scalars":[{"name":"alpha","diffusivity_m2_per_s":0.01}],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoint-${performance_directory}","write_interval":${checkpoint_interval}},"diagnostics":{"directory":"diagnostics-${performance_directory}","write_interval":${diagnostic_interval},"write_mesh":false},"performance":{"enabled":true,"directory":"${performance_directory}","warmup_steps":1,"measured_steps":2,"repetitions":${repetitions}}}
JSON
}

write_case "${work_root}/portable.json" 2 100 100 performance-portable
"${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/portable.json" \
  >"${work_root}/portable.stdout" 2>"${work_root}/portable.stderr"
test ! -s "${work_root}/portable.stdout"
test ! -s "${work_root}/portable.stderr"
artifact="${work_root}/performance-portable/performance.v1.json"
evidence_input="${work_root}/performance-portable/performance-evidence-input.v1.json"
test -s "${artifact}"
test -s "${evidence_input}"
test ! -e "${artifact}.tmp"
test ! -e "${evidence_input}.tmp"
test ! -e "${evidence_input}.backup"
grep -Fq '"schema_version":1' "${artifact}"
grep -Fq '"correctness":{"passed":true' "${artifact}"
grep -Fq '"repetitions":2' "${artifact}"
grep -Fq \
  "allocation-bytes-per-owned-cell=4060800000000000;peak-allocation-bytes-per-owned-cell=${peak_ratio_bits}" \
  "${artifact}"
grep -Eq '"commit":"([0-9a-f]{40}|unavailable)"' "${artifact}"
grep -Eq \
  '"working_tree":\{"clean":true,"dirty_summary":""\}|"working_tree":\{"clean":false,"dirty_summary":"[^"]+"\}' \
  "${artifact}"
grep -Eq '"compiler":\{"identity":"[^"]+","version":"[^"]+","flags":"[^"]+"\}' \
  "${artifact}"
grep -Eq '"link_flags":"[^"]+"' "${artifact}"
grep -Fq '"hardware_identity":"unavailable"' "${artifact}"
grep -Fq '"cpu_affinity":"unavailable"' "${artifact}"
grep -Eq '"node_identity":"(unavailable|[^"]+)"' "${artifact}"
grep -Eq '"rank_placement":"0:[^"]+"' "${artifact}"
if grep -Fq '\u0000' "${artifact}"; then
  echo "MPI implementation identity contains a trailing NUL" >&2
  exit 1
fi
grep -Fq '"checkpoint":0,"fp64-reduction":' "${artifact}"
grep -Fq '"logical_io_bytes":{"checkpoint":0,"diagnostics":0}' "${artifact}"
grep -Fq \
  "\"allocated_bytes\":{\"execution.allocated\":${portable_allocated},\"execution.peak-live\":${portable_peak}}" \
  "${artifact}"
grep -Fq \
  "\"halo_payload_bytes\":{\"pack\":${portable_halo_bytes},\"receive\":${portable_halo_bytes},\"send\":${portable_halo_bytes},\"unpack\":${portable_halo_bytes}}" \
  "${artifact}"
grep -Fq \
  "\"halo_messages\":{\"receive\":${portable_halo_messages},\"send\":${portable_halo_messages}}" \
  "${artifact}"
grep -Fq \
  '"collectives":{"checkpoint":0,"fp64-reduction":632,"linear-reduction":440}' \
  "${artifact}"
grep -Fq \
  '"collective_logical_payload_bytes":{"fp64-reduction":9440}' \
  "${artifact}"
grep -Fq '"matvec":{"momentum":24,"pressure":16}' "${artifact}"
grep -Fq \
  '"preconditioner_applications":{"momentum":0,"pressure":0}' "${artifact}"
test "$(grep -o '"relative_rank":[0-9]*' "${artifact}" | wc -l)" \
  -eq "$((2 * ranks))"
test "$(grep -o '"elapsed_seconds":[0-9][0-9.eE+-]*' "${artifact}" |
  wc -l)" -eq "$((2 * ranks))"

performance_binary_directory=$(dirname -- "${hundun}")
performance_oracle="${performance_binary_directory}/task25_performance_oracle"
if test ! -x "${performance_oracle}"; then
  performance_oracle="${performance_binary_directory}/../tests/task25_performance_oracle"
fi
test -x "${performance_oracle}"
IFS=, read -r cells_x cells_y cells_z <<EOF
${cells}
EOF
IFS=, read -r grid_x grid_y grid_z <<EOF
${grid}
EOF
"${mpiexec}" -n "${ranks}" "${performance_oracle}" "${artifact}" \
  "${cells_x}" "${cells_y}" "${cells_z}" \
  "${grid_x}" "${grid_y}" "${grid_z}" 1

oracle_parse_error='unable to read or parse Task 25 artifact on rank zero'
expect_invalid_oracle_input() {
  input=$1
  label=$2
  if "${mpiexec}" -n "${ranks}" "${performance_oracle}" "${input}" \
      "${cells_x}" "${cells_y}" "${cells_z}" \
      "${grid_x}" "${grid_y}" "${grid_z}" 1 \
      >"${work_root}/oracle-${label}.stdout" \
      2>"${work_root}/oracle-${label}.stderr"; then
    echo "invalid ${label} artifact unexpectedly passed the C++/MPI oracle" >&2
    exit 1
  fi
  grep -Fq "${oracle_parse_error}" \
    "${work_root}/oracle-${label}.stderr"
}

expect_invalid_oracle_input \
  "${work_root}/missing-performance-artifact.json" missing
printf '{not-json\n' >"${work_root}/corrupt-performance-artifact.json"
expect_invalid_oracle_input \
  "${work_root}/corrupt-performance-artifact.json" corrupt

mutated_artifact="${work_root}/mutated-counter.json"
sed "0,/\"send\":${portable_halo_bytes}/s//\"send\":$((portable_halo_bytes + 1))/" \
  "${artifact}" >"${mutated_artifact}"
if "${mpiexec}" -n "${ranks}" "${performance_oracle}" \
    "${mutated_artifact}" "${cells_x}" "${cells_y}" "${cells_z}" \
    "${grid_x}" "${grid_y}" "${grid_z}" 1 \
    >"${work_root}/oracle-mutation.stdout" \
    2>"${work_root}/oracle-mutation.stderr"; then
  echo "mutated exact counter unexpectedly passed the C++/MPI oracle" >&2
  exit 1
fi
grep -Fq \
  'Task 25 artifact differs from independent exact-counter oracle' \
  "${work_root}/oracle-mutation.stderr"
sed "0,/\"receive\":${portable_halo_bytes}/s//\"receive\":$((portable_halo_bytes + 1))/" \
  "${artifact}" >"${mutated_artifact}"
if "${mpiexec}" -n "${ranks}" "${performance_oracle}" \
    "${mutated_artifact}" "${cells_x}" "${cells_y}" "${cells_z}" \
    "${grid_x}" "${grid_y}" "${grid_z}" 1 \
    >"${work_root}/oracle-receive-mutation.stdout" \
    2>"${work_root}/oracle-receive-mutation.stderr"; then
  echo "mutated receive counter unexpectedly passed the C++/MPI oracle" >&2
  exit 1
fi
grep -Fq \
  'Task 25 artifact differs from independent exact-counter oracle' \
  "${work_root}/oracle-receive-mutation.stderr"

missing_key_artifact="${work_root}/missing-counter-key.json"
sed 's/"execution.allocated":[0-9]*,//' \
  "${artifact}" >"${missing_key_artifact}"
expect_invalid_oracle_input "${missing_key_artifact}" missing-counter-key

extra_key_artifact="${work_root}/extra-counter-key.json"
sed 's/"execution.peak-live":\([0-9]*\)/"execution.peak-live":\1,"unapproved":0/' \
  "${artifact}" >"${extra_key_artifact}"
expect_invalid_oracle_input "${extra_key_artifact}" extra-counter-key

extra_map_artifact="${work_root}/extra-counter-map.json"
sed 's/"logical_io_bytes":/"unapproved_map":{},"logical_io_bytes":/' \
  "${artifact}" >"${extra_map_artifact}"
expect_invalid_oracle_input "${extra_map_artifact}" extra-counter-map

old_sha=$(sha256sum "${artifact}" | cut -d' ' -f1)
old_evidence_sha=$(sha256sum "${evidence_input}" | cut -d' ' -f1)

expect_injected_failure() {
  injection=$1
  expected=$2
  if env HUNDUN_PERFORMANCE_FAILURE="${injection}" \
      "${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/portable.json" \
      >"${work_root}/${injection}.stdout" \
      2>"${work_root}/${injection}.stderr"; then
    echo "injected performance ${injection} unexpectedly succeeded" >&2
    exit 1
  fi
  grep -Fq "${expected}" "${work_root}/${injection}.stderr"
  test "$(sha256sum "${artifact}" | cut -d' ' -f1)" = "${old_sha}"
  test "$(sha256sum "${evidence_input}" | cut -d' ' -f1)" = \
    "${old_evidence_sha}"
  test ! -e "${artifact}.tmp"
  test ! -e "${evidence_input}.tmp"
  test ! -e "${evidence_input}.backup"
}

expect_injected_failure root_overflow \
  'unable to serialize performance artifact'
expect_injected_failure missing_work \
  'unable to serialize performance artifact'
expect_injected_failure malformed_work \
  'unable to serialize performance artifact'
expect_injected_failure reordered_work \
  'unable to serialize performance artifact'
expect_injected_failure metadata_mismatch \
  'unable to serialize performance artifact'
expect_injected_failure missing_counter_key \
  'unable to serialize performance artifact'
expect_injected_failure extra_counter_key \
  'unable to serialize performance artifact'
expect_injected_failure stage \
  'unable to publish performance artifact'
expect_injected_failure evidence_stage \
  'unable to publish performance artifact'
expect_injected_failure evidence_rename \
  'unable to publish performance artifact'
expect_injected_failure rename \
  'unable to publish performance artifact'
if test "${ranks}" -gt 1; then
  expect_injected_failure nonroot_metadata \
    'performance metadata/controller differs between ranks'
  expect_injected_failure nonroot_work_key \
    'performance work reports differ between ranks'
  expect_injected_failure nonroot_work_count \
    'performance work report counts differ between ranks'
  expect_injected_failure nonroot_scalar_shape \
    'performance scalar collection shape differs between ranks'
  expect_injected_failure nonroot_path \
    'unable to derive performance artifact path'
  expect_injected_failure nonroot_diagnostic_path \
    'unable to derive performance diagnostics path'
  expect_injected_failure nonroot_checkpoint_path \
    'unable to derive performance checkpoint path'
  expect_injected_failure nonroot_nonfinite_io \
    'performance I/O duration contains an invalid rank sample'
  expect_injected_failure nonroot_negative_io \
    'performance I/O duration contains an invalid rank sample'
fi

sed 's/"steps":3/"steps":2/' "${work_root}/portable.json" \
  >"${work_root}/mismatched-steps.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/mismatched-steps.json" \
    >"${work_root}/mismatched.stdout" \
    2>"${work_root}/mismatched.stderr"; then
  echo "mismatched performance step count unexpectedly succeeded" >&2
  exit 1
fi
if grep -Fq 'HUNDUN-FLOW 0.1.0' \
    "${work_root}/mismatched.stdout"; then
  echo "performance preflight rejection printed the Stage 2 banner" >&2
  exit 1
fi
grep -Fq \
  'performance time.steps must equal warmup_steps + measured_steps' \
  "${work_root}/mismatched.stderr"
test "$(sha256sum "${artifact}" | cut -d' ' -f1)" = "${old_sha}"
test ! -e "${artifact}.tmp"

sed \
  -e 's/"read":false/"read":true,"read_directory":"missing-restart"/' \
  "${work_root}/portable.json" >"${work_root}/restart-request.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/restart-request.json" \
    >"${work_root}/restart.stdout" 2>"${work_root}/restart.stderr"; then
  echo "performance restart request unexpectedly succeeded" >&2
  exit 1
fi
grep -Fq 'performance mode does not accept restart input' \
  "${work_root}/restart.stderr"
test "$(sha256sum "${artifact}" | cut -d' ' -f1)" = "${old_sha}"

printf '%s\n' blocked >"${work_root}/performance-blocked"
sed 's/"directory":"performance-portable"/"directory":"performance-blocked"/' \
  "${work_root}/portable.json" >"${work_root}/blocked-output.json"
if "${mpiexec}" -n "${ranks}" "${hundun}" \
    "${work_root}/blocked-output.json" \
    >"${work_root}/blocked.stdout" 2>"${work_root}/blocked.stderr"; then
  echo "blocked performance publication unexpectedly succeeded" >&2
  exit 1
fi
grep -Fq 'unable to publish performance artifact' \
  "${work_root}/blocked.stderr"
test -f "${work_root}/performance-blocked"
test ! -e "${work_root}/performance-blocked/performance.v1.json.tmp"
test ! -e "${work_root}/performance-blocked/performance.v1.json"
test ! -e \
  "${work_root}/performance-blocked/performance-evidence-input.v1.json.tmp"
test ! -e \
  "${work_root}/performance-blocked/performance-evidence-input.v1.json"
test ! -e \
  "${work_root}/performance-blocked/performance-evidence-input.v1.json.backup"

write_case "${work_root}/positive-io.json" 1 2 2 performance-positive-io
"${mpiexec}" -n "${ranks}" "${hundun}" "${work_root}/positive-io.json" \
  >"${work_root}/positive.stdout" 2>"${work_root}/positive.stderr"
test ! -s "${work_root}/positive.stdout"
test ! -s "${work_root}/positive.stderr"
positive="${work_root}/performance-positive-io/performance.v1.json"
positive_evidence_input="${work_root}/performance-positive-io/performance-evidence-input.v1.json"
test -s "${positive}"
test -s "${positive_evidence_input}"
test ! -e "${positive}.tmp"
grep -Fq \
  "allocation-bytes-per-owned-cell=4060800000000000;peak-allocation-bytes-per-owned-cell=${peak_ratio_bits}" \
  "${positive}"
if grep -Fq \
    '"logical_io_bytes":{"checkpoint":0,"diagnostics":0}' "${positive}"; then
  echo "positive-I/O performance artifact recorded zero I/O" >&2
  exit 1
fi
grep -Eq \
  '"logical_io_bytes":\{"checkpoint":[1-9][0-9]*,"diagnostics":[1-9][0-9]*\}' \
  "${positive}"
grep -Eq '"collectives":\{"checkpoint":[1-9][0-9]*,' "${positive}"
grep -Fq \
  "\"allocated_bytes\":{\"execution.allocated\":${positive_allocated},\"execution.peak-live\":${positive_peak}}" \
  "${positive}"
grep -Fq \
  "\"halo_payload_bytes\":{\"pack\":${positive_halo_bytes},\"receive\":${positive_halo_bytes},\"send\":${positive_halo_bytes},\"unpack\":${positive_halo_bytes}}" \
  "${positive}"
grep -Fq \
  "\"halo_messages\":{\"receive\":${positive_halo_messages},\"send\":${positive_halo_messages}}" \
  "${positive}"
grep -Fq \
  '"collectives":{"checkpoint":32,"fp64-reduction":316,"linear-reduction":220}' \
  "${positive}"
grep -Fq \
  '"collective_logical_payload_bytes":{"fp64-reduction":4720}' "${positive}"
grep -Fq '"matvec":{"momentum":12,"pressure":8}' "${positive}"
grep -Fq \
  "\"logical_io_bytes\":{\"checkpoint\":${positive_checkpoint_bytes},\"diagnostics\":${positive_diagnostic_bytes}}" \
  "${positive}"

evidence_writer="${performance_binary_directory}/task25_performance_evidence"
if test ! -x "${evidence_writer}"; then
  evidence_writer="${performance_binary_directory}/../tests/task25_performance_evidence"
fi
test -x "${evidence_writer}"
portable_command_sha=$(printf '%s' \
  "${mpiexec} -n ${ranks} ${hundun} ${work_root}/portable.json" |
  sha256sum | cut -d' ' -f1)
portable_case_sha=$(sha256sum "${work_root}/portable.json" | cut -d' ' -f1)
positive_command_sha=$(printf '%s' \
  "${mpiexec} -n ${ranks} ${hundun} ${work_root}/positive-io.json" |
  sha256sum | cut -d' ' -f1)
positive_case_sha=$(sha256sum "${work_root}/positive-io.json" | cut -d' ' -f1)
positive_sha=$(sha256sum "${positive}" | cut -d' ' -f1)
evidence_manifest="${work_root}/performance-evidence-input.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tnull\tnull\n' \
  "task25-portable-${ranks}rank" "${portable_command_sha}" \
  "${portable_case_sha}" "${old_sha}" "${work_root}/portable.json" \
  "${artifact}" \
  "${work_root}/performance-portable/performance-evidence-input.v1.json" \
  >"${evidence_manifest}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tnull\tnull\n' \
  "task25-positive-io-${ranks}rank" "${positive_command_sha}" \
  "${positive_case_sha}" "${positive_sha}" \
  "${work_root}/positive-io.json" "${positive}" \
  "${positive_evidence_input}" >>"${evidence_manifest}"
"${evidence_writer}" "${evidence_manifest}" \
  "${work_root}/performance-evidence.v1.json"
test -s "${work_root}/performance-evidence.v1.json"
test ! -e "${work_root}/performance-evidence.v1.json.tmp"
grep -Fq '"family":"runtime"' \
  "${work_root}/performance-evidence.v1.json"
grep -Fq '"family":"pressure"' \
  "${work_root}/performance-evidence.v1.json"
grep -Fq '"family":"combined"' \
  "${work_root}/performance-evidence.v1.json"
grep -Fq '"throughput_bytes_per_second":{"unit":"byte/s","status":"available"' \
  "${work_root}/performance-evidence.v1.json"

if test -n "${HUNDUN_TASK25_EVIDENCE_DIRECTORY:-}"; then
  mkdir -p -- "${HUNDUN_TASK25_EVIDENCE_DIRECTORY}"
  cp -- "${artifact}" \
    "${HUNDUN_TASK25_EVIDENCE_DIRECTORY}/performance-portable-${ranks}rank.v1.json"
  cp -- "${positive}" \
    "${HUNDUN_TASK25_EVIDENCE_DIRECTORY}/performance-positive-io-${ranks}rank.v1.json"
fi
