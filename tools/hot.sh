#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Production cadence: compute every step, Visit every 100, Restart every 500.
# Positional arguments: executable case restart output steps [ranks].
set -euo pipefail
if (( $# < 5 || $# > 6 )); then
  echo 'usage: hot.sh executable case restart output steps [ranks]' >&2
  exit 2
fi
exec mpiexec -np "${6:-128}" --bind-to core "$1" run "$2" --restart "$3" \
  --output "$4" --steps "$5" --output-interval "${HF_PLOT_EVERY:-100}" \
  --restart-interval "${HF_SAVE_EVERY:-500}" --diagnostics-interval "${HF_DIAG_EVERY:-10}"
