#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
set -euo pipefail
base=/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52
audit="$base/fgmres-observation-20260908"
source="$base/long-thermodynamic-35000-20260907"
physical="$base/method-frozen-thermodynamic-48493ed-20260907"
output="$base/pilot-fgmres-9500-9510-20260908"
service=hundun-re3900-thermodynamic-long-20260907.service
# The frozen long job retains its original lock and process state under SIGSTOP.
# This maintenance slot permits only that exact, fully stopped job to coexist.
# It neither releases the long-job lock nor sends SIGCONT at exit.
exec 9>"$base/.hundun-mpi-maintenance.lock"
flock -n 9 || { echo 'maintenance slot occupied' >&2; exit 3; }
if pgrep -x 'ctest|ninja|make|clang|clang-15|cc1plus' >/dev/null; then
  echo 'another test/compiler process exists' >&2; exit 3
fi
test "$(pgrep -xc v04_thin_domain)" = 128
for pid in $(pgrep -x 'mpirun|mpiexec|v04_thin_domain|hundun|coast|COAST'); do
  state=$(ps -p "$pid" -o stat=)
  case "$state" in T*) ;; *) echo "active competing process $pid $state" >&2; exit 3;; esac
  rg -Fq "$service" "/proc/$pid/cgroup" || { echo "unowned paused process $pid" >&2; exit 3; }
done
test "$(<"$source/Restart/current")" = generation-9500-139424060480232
test ! -e "$output"
test ! -e "$output.log"
test -f "$audit/LOCAL_ACCEPTED.json"
cd "$audit"
sha256sum --check SOURCE.sha256 >/dev/null
sha256sum --check FROZEN.sha256 >/dev/null
set -o noclobber
exec >"$output.log" 2>&1
date --iso-8601=seconds
export LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu:/home/wyf/.local/opt/hundun-toolchain/clang/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# One diagnostic configuration, ten steps. No method recovery or physics change.
# Same-method restart clears ephemeral warm authority by its existing contract;
# report the first step separately, not as a steady throughput measurement.
exec /usr/bin/time -p /usr/bin/mpirun -n 128 --bind-to core --map-by core \
  "$audit/frozen/v04_thin_domain_runner" --spec "$physical/statistics-long-20plus50D.d" \
  --case-root "$physical/case" --run-root "$output" --restart-root "$source/Restart" \
  --steps 10 --visit-interval 500 --observe-performance --observe-mg-cost --observe-fgmres-recovery
