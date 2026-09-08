#!/usr/bin/env bash
set -euo pipefail
base=/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52
audit="$base/coupled-stage46-20260909"
frozen="$audit/frozen"
source="$base/integrated-stage46-20260908/source-10510"
output="$base/long-coupled-35000-20260909"
exec 9>"$base/.hundun-mpi-exclusive.lock"
flock -n 9 || exit 3
exec 8>"$base/.hundun-mpi-maintenance.lock"
flock -n 8 || exit 3
if pgrep -x 'mpirun|mpiexec|v04_thin_domain|hundun|coast|COAST|ctest|ninja|make|clang|clang-15|cc1plus' >/dev/null; then
  echo 'another solver/compiler/test process exists' >&2; exit 3
fi
test -f "$audit/PILOT_ACCEPTED.json"
test ! -e "$output"
test ! -e "$output.log"
cd "$frozen"
sha256sum --check "$audit/FROZEN.sha256" >/dev/null
cd "$source"
sha256sum --check "$audit/SOURCE-LONG.sha256" >/dev/null
cd "$frozen"
set -o noclobber
exec >"$output.log" 2>&1
date --iso-8601=seconds
echo 'Native coupled source=83db41df2ec9c351c71e3aebe812509ddfd7262b; exact restart=10510; target=35000; no method recovery'
export LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu:/home/wyf/.local/opt/hundun-toolchain/clang/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /usr/bin/time -p /usr/bin/mpirun -n 128 --bind-to core --map-by core \
 "$frozen/v04_thin_domain_runner" --spec "$frozen/statistics-long-20plus50D.d" \
 --case-root "$frozen/case" --run-root "$output" --restart-root "$source/Restart" \
 --steps 24490 --visit-interval 500 --observe-performance --observe-mg-cost --observe-fgmres-recovery
