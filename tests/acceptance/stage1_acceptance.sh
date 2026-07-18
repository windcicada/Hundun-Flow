#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
project_root=$(CDPATH= cd -- "${script_dir}/../.." && pwd -P)
cmake_command=${CMAKE_COMMAND:-cmake}
ctest_command=${CTEST_COMMAND:-ctest}
mpiexec_command=${MPIEXEC_COMMAND:-mpiexec}
work_root=$(mktemp -d "${TMPDIR:-/tmp}/hundun-stage1-acceptance.XXXXXX")
build_dir="${work_root}/build"
case_dir="${work_root}/case"

cleanup() {
  rm -rf -- "${work_root}"
}
trap cleanup EXIT

export OMPI_MCA_rmaps_base_oversubscribe=1

mkdir -p -- "${case_dir}" "${work_root}/launch-a" "${work_root}/launch-b"
cp -- "${project_root}/cases/passive_scalar/case.json" "${case_dir}/case.json"
cp -- "${project_root}/cases/passive_scalar/case_restart.json" \
  "${case_dir}/case_restart.json"

"${cmake_command}" -S "${project_root}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=ON
"${cmake_command}" --build "${build_dir}" -j 2

test "$("${build_dir}/hundun" --version)" = "HUNDUN-FLOW 0.0.0-stage1"
test "$("${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case.json" --validate)" = "VALID"

(
  cd -- "${work_root}/launch-a"
  "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
    "${case_dir}/case.json" --print-resolved
) >"${work_root}/resolved-a.json"
(
  cd -- "${work_root}/launch-b"
  "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
    "${case_dir}/case.json" --print-resolved
) >"${work_root}/resolved-b.json"
cmp --silent "${work_root}/resolved-a.json" "${work_root}/resolved-b.json"
grep -F '"process_grid":[2,1,1]' "${work_root}/resolved-a.json" >/dev/null
grep -F '"restart":{"read":false,"write_directory":"Restart"}' \
  "${work_root}/resolved-a.json" >/dev/null

mkdir -p -- "${case_dir}/output"
"${mpiexec_command}" -n 2 "${build_dir}/hundun" "${case_dir}/case.json" \
  >"${case_dir}/output/run.log" 2>&1
grep -E '^CASE name=periodic_passive_scalar ranks=2 cells=64x8x8$' \
  "${case_dir}/output/run.log" >/dev/null
grep -E '^STEP 10 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$' \
  "${case_dir}/output/run.log" >/dev/null
grep -E '^STEP 20 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$' \
  "${case_dir}/output/run.log" >/dev/null
grep -E '^FINISHED step=20 time_s=[0-9.eE+-]+$' \
  "${case_dir}/output/run.log" >/dev/null

check_checkpoint() {
  local directory=$1
  test -s "${directory}/manifest.v1.bin"
  test -s "${directory}/COMPLETED"
  test "$(find "${directory}" -maxdepth 1 -type f \
    -name 'restart.rank*.bin' | wc -l)" -eq 2
  test -s "${directory}/restart.rank000000.bin"
  test -s "${directory}/restart.rank000001.bin"
  local records
  records=$(od -An -tu4 -j 75 -N 4 "${directory}/manifest.v1.bin" | \
    tr -d '[:space:]')
  test "${records}" = "2"
}

check_checkpoint "${case_dir}/Restart/step00000010"
check_checkpoint "${case_dir}/Restart/step00000020"
test "$(find "${case_dir}/output" -maxdepth 1 -type f \
  -name 'scalar.step*.rank*.vtk' | wc -l)" -eq 4

mv -- "${case_dir}/Restart/step00000010/COMPLETED" \
  "${case_dir}/Restart/step00000010/COMPLETED.saved"
if "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case_restart.json" >"${work_root}/missing-marker.log" 2>&1; then
  echo "restart unexpectedly accepted a missing COMPLETED marker" >&2
  exit 1
fi
mv -- "${case_dir}/Restart/step00000010/COMPLETED.saved" \
  "${case_dir}/Restart/step00000010/COMPLETED"

cp -p -- "${case_dir}/Restart/step00000010/restart.rank000001.bin" \
  "${case_dir}/Restart/step00000010/restart.rank000001.bin.saved"
: >"${case_dir}/Restart/step00000010/restart.rank000001.bin"
if "${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case_restart.json" >"${work_root}/truncated-rank.log" 2>&1; then
  echo "restart unexpectedly accepted a truncated rank file" >&2
  exit 1
fi
mv -- "${case_dir}/Restart/step00000010/restart.rank000001.bin.saved" \
  "${case_dir}/Restart/step00000010/restart.rank000001.bin"

mkdir -p -- "${case_dir}/output.resumed"
"${mpiexec_command}" -n 2 "${build_dir}/hundun" \
  "${case_dir}/case_restart.json" >"${case_dir}/output.resumed/run.log" 2>&1
grep -E '^FINISHED step=20 time_s=[0-9.eE+-]+$' \
  "${case_dir}/output.resumed/run.log" >/dev/null
check_checkpoint "${case_dir}/Restart.resumed/step00000020"

for rank in 000000 000001; do
  "${build_dir}/compare_scalar_vtk" \
    "${case_dir}/output/scalar.step00000020.rank${rank}.vtk" \
    "${case_dir}/output.resumed/scalar.step00000020.rank${rank}.vtk" \
    1e-13
done
if "${build_dir}/compare_scalar_vtk" \
  "${case_dir}/output/scalar.step00000010.rank000000.vtk" \
  "${case_dir}/output/scalar.step00000020.rank000000.vtk" 1e-13; then
  echo "VTK comparator unexpectedly accepted distinct generated fields" >&2
  exit 1
fi

"${ctest_command}" --test-dir "${build_dir}" -j 1 --output-on-failure
"${ctest_command}" --test-dir "${build_dir}" -L provenance -j 1 \
  --output-on-failure

if ldd "${build_dir}/hundun" | \
    grep -E -i 'python|not found|coast|boffin'; then
  echo "hundun has forbidden or missing dynamic linkage" >&2
  exit 1
fi
if find "${case_dir}" -type f -name '*.tmp' -print -quit | grep -q .; then
  echo "generated case tree contains a residual temporary file" >&2
  exit 1
fi
