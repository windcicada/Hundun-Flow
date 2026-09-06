#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Exercise the public observer CLI, including old/new loop CSV contracts."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def write(path, records):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    script = sys.argv[1]
    with tempfile.TemporaryDirectory(prefix="hundun-observe-") as directory:
        root = Path(directory)
        base = dict(step=7, rank=0, attempt=1, corrector=1, refinement=0,
                    kind=0, dropped_loops=0, invoked=1, iterations=2,
                    A_calls=3, M_calls=2, A_ns=30, M_ns=20,
                    baseline_candidates=1, extrapolated_candidates=0,
                    ladder_candidates=0, incomplete_candidates=0,
                    linear_initial=1, linear_final=0.01, globalization_valid=0)
        total = dict(step=7, rank=0, dropped_loops=0, pressure_calls=1,
                     diagonal_calls=0, spatial_calls=0, A_apply_ns=30,
                     M_apply_ns=20)
        (root / "conservation.csv").write_text("step,rank\n7,0\n")

        def run(loops, performance, accepted):
            write(root / "solver-rank-0.csv", loops)
            write(root / "performance.csv", performance)
            result = subprocess.run([sys.executable, script, str(root)],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, timeout=10)
            assert (result.returncode == 0) == accepted, result.stderr or result.stdout
            return json.loads(result.stdout) if accepted else None

        old = run([base], [total], True)
        assert old["loops"][0]["scalar_coupling_sweep"] == 1
        assert old["rank_step_mean_ns"]["scalar_remap_ns"] == 0
        loops = [dict(base, scalar_coupling_sweep=i) for i in (1, 2)]
        summed = dict(total, pressure_calls=2, A_apply_ns=60, M_apply_ns=40,
                      scalar_remap_ns=9)
        new = run(loops, [summed], True)
        assert len(new["loops"]) == 2
        assert new["rank_step_mean_ns"]["scalar_remap_ns"] == 9
        run([loops[0], loops[0]], [summed], False)  # Duplicate logical loop.
        run([dict(base, dropped_loops=1)], [total], False)
        run([base], [dict(total, A_apply_ns=29)], False)
        run([base], [dict(total, pressure_calls=0, spatial_calls=1)], False)
        run([base, dict(base, step=8)], [total], False)  # Unaccounted step.
        run([base], [total, total], False)  # Duplicate performance row.
    print("solver observation CLI PASS (V1/V2, composition sweeps, invalid attribution)")


if __name__ == "__main__":
    main()
