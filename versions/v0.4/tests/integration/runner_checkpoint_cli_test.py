#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
# windcicada | Year.M: 2026.09
"""Real runner CLI, stationary IBM fixture, scoped real new failures."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--first-only", action="store_true")
    parser.add_argument("--visit", action="store_true")
    parser.add_argument("--only-index", type=int)
    parser.add_argument("--target-rank", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        args.output = Path(tempfile.mkdtemp(prefix="hundun-runner-checkpoint-")) / "audit"
    args.output.mkdir(parents=True, exist_ok=False)
    data = Path(__file__).resolve().parents[1] / "data"
    case = args.output / "case"
    case.mkdir()
    model = json.loads((data / "case_minimal_valid.json").read_text())
    model["mesh"]["domain"] = {"lower": [-2, -2, -2], "upper": [2, 2, 2]}
    model["mesh"]["exact_cells"] = [16, 16, 16]
    model["mesh"]["minimum_spacing"] = [.25, .25, .25]
    model["mesh"]["limits"]["max_memory_bytes_per_rank"] = 1073741824
    model["mesh"]["immersed_boundary"] = {"stl_file": "cylinder_ascii.stl", "fluid_side": "outside"}
    model["flow"]["pressure_reference"] = "closed_mass"
    for boundary in model["boundaries"].values():
        boundary["flow_kind"] = "periodic"
        boundary["velocity"] = [0, 0, 0]
        boundary["pressure"] = 0
    (case / "case.json").write_text(json.dumps(model))
    for name in ("thermophysics.d", "cylinder_ascii.stl"):
        shutil.copyfile(str(data / name), str(case / name))
    spec = args.output / "statistics.d"
    spec.write_text("HUNDUN_V04_LITERATURE_STATISTICS_V1\n"
                    "development_steps 10\ncollection_end_step 20\n"
                    "checkpoint_interval 1\nrho_ref 1\nu_ref 1\ndiameter 1\n"
                    "span 4\ncylinder_center_x 0\nstation_x_over_d 1.5\nend\n")
    results = []
    for target in ([args.target_rank] if args.target_rank is not None else sorted(set((0, args.ranks - 1)))):
        baseline_visit = {}
        def run(index):
            # Fixed width, including the baseline: path-length changes may
            # cross the standard library's small-string allocation threshold.
            path = args.output / ("rank-{}-index-{:07d}".format(target, index))
            env = dict(os.environ, LD_PRELOAD=str(args.probe.resolve()),
                       HUNDUN_TEST_ALLOC_RANK=str(target), HUNDUN_TEST_ALLOC_INDEX=str(index))
            if args.visit:
                env["HUNDUN_TEST_VISIT_ALLOC"] = "1"
            command = [args.mpi, "-n", str(args.ranks), str(args.binary.resolve()),
                       "--spec", str(spec.resolve()), "--case-root", str(case.resolve()),
                       "--run-root", str(path.resolve()), "--steps", "1", "--visit-interval", "1" if args.visit else "0"]
            if args.visit:
                command.append("--observe-performance")
            result = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=40,
                                    universal_newlines=True)
            (args.output / (path.name + ".log")).write_text(result.stdout)
            match = re.search(r"CHECKPOINT_ALLOC sites=(\d+) fired=(\d+) started=(\d+)", result.stdout)
            row = {"target": target, "index": index, "returncode": result.returncode,
                   "observation": match.groups() if match else None, "command": command}
            results.append(row)
            (args.output / "results.json").write_text(json.dumps(results, indent=2))
            assert match and match[3] == "1", result.stdout
            if index < 0:
                assert result.returncode == 0, result.stdout
            else:
                assert match[2] == "1", "selected real allocation did not fire"
                assert result.returncode in ((0, 6, 7) if args.visit else (0, 6)), result.stdout
                if result.returncode == 0 and not args.visit:
                    assert "CHECKPOINT_IO status=0/0" in result.stdout, result.stdout
                    assert list(path.glob("*.complete")), "cleanup warning lost committed marker"
            if args.visit:
                finalized = re.findall(r"VISIT_FINALIZE rank=(\d+) accepted_step=1", result.stdout)
                assert sorted(map(int, finalized)) == list(range(args.ranks)), result.stdout
                if "visit_status=" in result.stdout:
                    assert "accepted_step=1 committed_step=1" in result.stdout, result.stdout
                    times = re.search(r"accepted_time=(\S+) committed_time=(\S+)", result.stdout)
                    assert times and float(times[1]) == float(times[2]), result.stdout
            if result.returncode == 0:
                with (path / "health.csv").open() as stream:
                    health_rows = list(csv.DictReader(stream))
                assert len(health_rows) == 1
                health = health_rows[0]
                # The fixture STL is a cube [-1,1]^3 inside the 16^3 domain.
                assert int(health["solid_placeholder_count"]) == 8**3
                assert int(health["fluid_count"]) == 16**3 - 8**3
                for region in ("fluid", "solid_placeholder"):
                    for value in ("p_min", "p_max", "T_min", "T_max", "rho_min", "rho_max"):
                        assert math.isfinite(float(health[region + "_" + value]))
                        cell = int(health[region + "_" + value + "_cell"])
                        i, j, k = cell % 16, (cell // 16) % 16, cell // 256
                        solid = all(4 <= v < 12 for v in (i, j, k))
                        assert solid == (region == "solid_placeholder")
                if args.visit:
                    indices = list((path / "Visit").glob("*.visit"))
                    assert len(indices) == 1, "Visit was not exercised"
                    assert indices[0].read_text().startswith("!NBLOCKS {}\n".format(args.ranks))
                    payload = {p.name: p.read_bytes() for p in (path / "Visit").iterdir() if p.is_file()}
                    if index < 0:
                        baseline_visit.update(payload)
                    else:
                        assert payload == baseline_visit, "successful Visit changed accepted fields"
                with (path / "conservation.csv").open() as stream:
                    balances = list(csv.DictReader(stream))
                assert len(balances) == 1 and balances[0]["step"] == "1"
                balance = balances[0]
                assert balance["epoch_start_step"] == "0"
                assert all(math.isfinite(float(value)) for value in balance.values())
                assert float(balance["mass_kg"]) > 0.0
                assert abs(float(balance["cumulative_mass_defect_kg"])) < 1e-10
                assert abs(float(balance["cumulative_energy_defect_J"])) < 1e-6
                marker = path / "step-00000000000000000001.complete"
                assert marker.is_file() and marker.read_text().endswith("end\n")
                statistics = path / "step-00000000000000000001.statistics.json"
                assert json.loads(statistics.read_text())["snapshot_step"] == 1
                accumulator = path / "step-00000000000000000001.accumulator"
                assert accumulator.read_text().endswith("end\n")
                if index < 0 and target == 0 and not args.visit:
                    resumed = args.output / "method-recovery"
                    recovery_command = [args.mpi, "-n", str(args.ranks),
                        str(args.binary.resolve()), "--spec", str(spec.resolve()),
                        "--case-root", str(case.resolve()), "--run-root", str(resumed.resolve()),
                        "--restart-root", str((path / "Restart").resolve()),
                        "--restart-method-recovery", "--steps", "1", "--visit-interval", "0"]
                    recovery = subprocess.run(recovery_command, env=os.environ.copy(),
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=40,
                        universal_newlines=True)
                    (args.output / "method-recovery.log").write_text(recovery.stdout)
                    assert recovery.returncode == 0, recovery.stdout
                    with (resumed / "health.csv").open() as stream:
                        resumed_health = list(csv.DictReader(stream))
                    assert len(resumed_health) == 1
                    assert resumed_health[0]["step"] == "2"
                    assert resumed_health[0]["bdf_order"] == "1"
                    assert resumed_health[0]["attempts"] == "1"
                    assert resumed_health[0]["restart_recovery"] == "1"
                    assert "restart_method_recovery 1\n" in (resumed / "RUN.meta").read_text()
                    with (resumed / "conservation.csv").open() as stream:
                        resumed_balance = list(csv.DictReader(stream))
                    assert resumed_balance[0]["epoch_start_step"] == "1"
            return int(match[1])
        count = run(-1)
        assert count > 0
        for index in ([args.only_index] if args.only_index is not None else range(1 if args.first_only else count)):
            run(index)
        print("runner {} ranks={} target={} sites={} tested={}".format(
            "visit" if args.visit else "checkpoint", args.ranks, target, count,
            1 if args.first_only or args.only_index is not None else count), flush=True)


if __name__ == "__main__":
    main()
