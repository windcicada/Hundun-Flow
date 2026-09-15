#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Sample one MPI rank with a user-space profiler in a diagnostic run.
set -euo pipefail
if (( $# < 3 )); then
  echo 'usage: cpu.sh <libprofiler.so> <profile-output> <program> [args...]' >&2
  exit 2
fi
profile_library=$1
profile_output=$2
shift 2
profile_rank=${HUNDUN_PROFILE_RANK:-0}
process_rank=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}
if [[ "$process_rank" == "$profile_rank" ]]; then
  export LD_PRELOAD="$profile_library${LD_PRELOAD:+:$LD_PRELOAD}"
  export CPUPROFILE="$profile_output"
  export CPUPROFILE_FREQUENCY=200
fi
exec "$@"
