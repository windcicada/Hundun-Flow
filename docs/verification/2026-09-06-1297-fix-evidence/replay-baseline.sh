#!/usr/bin/env bash
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
set -euo pipefail
root=/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52
frozen="$root/review-frozen-20260905"
restart="$root/long-review-749to35000-20260905/Restart"
output="$root/replay-baseline-1000to1296-20260906"
expected=3737f83420ce255c0268cf655855284a0690100a2af6456a85e075604952ac38
actual=$(sha256sum "$frozen/v04_thin_domain_runner")
[[ ${actual%% *} == "$expected" ]]
[[ ! -e "$output" ]]
read -r generation < "$restart/current"
[[ $generation == generation-1000-* ]]
export LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu:/home/wyf/.local/opt/hundun-toolchain/clang/lib
export OMP_NUM_THREADS=1
exec /usr/bin/time -p timeout --signal=TERM --kill-after=30s 9000s \
  mpirun -n 128 --bind-to core --map-by core \
  "$frozen/v04_thin_domain_runner" \
  --spec "$frozen/statistics-long-20plus50D.d" \
  --case-root "$frozen/case" --restart-root "$restart" \
  --run-root "$output" --steps 296 --visit-interval 0 --observe-performance
