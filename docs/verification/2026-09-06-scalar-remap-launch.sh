#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
set -euo pipefail
base=/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52
frozen="$base/method-frozen-scalar-20260906"
pilot="$base/scalar-recovery-pilot-20260906"
case "${1:-}" in
  pilot)
    output="$pilot"
    restart="$base/long-periodic-fresh-35000-20260906"
    steps=10
    visit=10
    policy=(--restart-method-recovery)
    ;;
  long)
    output="$base/long-scalar-repaired-35000-20260906"
    restart="$pilot"
    steps=33990
    visit=500
    policy=()
    ;;
  *) echo 'usage: launch.sh pilot|long' >&2; exit 2 ;;
esac
# Cooperating jobs share the lock. Also refuse any existing MPI launcher,
# including jobs not started with this script. Never stop another job here.
exec 9>"$base/.hundun-mpi-exclusive.lock"
flock -n 9 || { echo 'another HUNDUN MPI launch holds the lock' >&2; exit 3; }
if pgrep -x 'mpirun|mpiexec' >/dev/null; then
  echo 'an MPI job is already running; serialize before launch' >&2
  exit 3
fi
test ! -e "$output"
test ! -e "$output.log"
test -f "$restart/Restart/current"
cd "$frozen"
sha256sum --check FROZEN.sha256
set -o noclobber
exec >"$output.log" 2>&1
date --iso-8601=seconds
export LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu:/home/wyf/.local/opt/hundun-toolchain/clang/lib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# The pilot establishes the new method/statistics epoch at step 1000. Long
# continuation inherits that exact V3 history and ends at absolute step 35000.
exec /usr/bin/time -p /usr/bin/mpirun -n 128 --bind-to core --map-by core \
  "$frozen/v04_thin_domain_runner" \
  --spec "$frozen/statistics-long-20plus50D.d" \
  --case-root "$frozen/case" --run-root "$output" \
  --restart-root "$restart/Restart" "${policy[@]}" --steps "$steps" \
  --visit-interval "$visit" --observe-performance
