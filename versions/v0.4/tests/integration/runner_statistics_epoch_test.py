#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Real CLI: nonzero statistics, exact continuation, method reset and resegmentation."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def metadata(path):
    lines = path.read_text().splitlines()[1:-1]
    return dict(line.split(" ", 1) for line in lines)


def inventory(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix="hundun-epoch-")) / "audit"
    root.mkdir(parents=True, exist_ok=False)
    data = Path(__file__).resolve().parents[1] / "data"
    case = root / "case"
    case.mkdir()
    model = json.loads((data / "case_minimal_valid.json").read_text())
    model["time"]["scheme"] = "backward_euler"
    model["mesh"]["domain"] = {"lower": [-2, -2, -2], "upper": [2, 2, 2]}
    model["mesh"]["exact_cells"] = [16, 16, 16]
    model["mesh"]["minimum_spacing"] = [.25, .25, .25]
    model["mesh"]["limits"]["max_memory_bytes_per_rank"] = 1073741824
    model["mesh"]["immersed_boundary"] = {"stl_file": "cylinder_ascii.stl", "fluid_side": "outside"}
    model["flow"]["pressure_reference"] = "closed_mass"
    for boundary in model["boundaries"].values():
        boundary.update(flow_kind="periodic", velocity=[0, 0, 0], pressure=0)
    (case / "case.json").write_text(json.dumps(model))
    for name in ("thermophysics.d", "cylinder_ascii.stl"):
        shutil.copyfile(str(data / name), str(case / name))
    spec = root / "statistics.d"
    spec.write_text("HUNDUN_V04_LITERATURE_STATISTICS_V1\n"
                    "development_steps 3\ncollection_end_step 20\ncheckpoint_interval 1\n"
                    "rho_ref 1\nu_ref 1\ndiameter 1\nspan 4\ncylinder_center_x 0\n"
                    "station_x_over_d 1.5\nend\n")
    def run(name, steps, source=None, method=False, development=None):
        path = root / name
        command = [args.mpi, "-n", str(args.ranks), str(args.binary.resolve()),
                   "--spec", str(spec.resolve()), "--case-root", str(case.resolve()),
                   "--run-root", str(path.resolve()), "--steps", str(steps), "--visit-interval", "0",
                   "--observe-performance"]
        if source:
            command += ["--restart-root", str((source / "Restart").resolve())]
        if method:
            command += ["--restart-method-recovery"]
        if development is not None:
            command += ["--restart-development-steps", str(development)]
        result = subprocess.run(command, env=os.environ.copy(), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=40, universal_newlines=True)
        (root / (name + ".log")).write_text(result.stdout)
        assert result.returncode == 0, result.stdout
        run_meta = metadata(path / "RUN.meta")
        assert run_meta["expected_ranks"] == str(args.ranks)
        assert run_meta["requested_steps"] == str(steps)
        assert run_meta["observation_schema"] == "4"
        source_hash = hashlib.sha256((path / "RUN.meta").read_bytes()).hexdigest()
        with (path / "performance.csv").open() as stream:
            performance = list(csv.DictReader(stream))
        with (path / "conservation.csv").open() as stream:
            balances = list(csv.DictReader(stream))
        for row in balances:
            assert int(row["ibm_adjacent_cells"]) > 0 and int(row["interior_cells"]) > 0
            assert row["normalization_valid"] == "0", "zero-velocity normalization invented a reference"
        for row in performance:
            assert row["source_meta_sha256"] == source_hash
            assert int(row["final_momentum_ns"]) > 0 and int(row["terminal_metrics_ns"]) > 0
            assert int(row["boundary_ledger_ns"]) > 0 and row["dropped_loops"] == "0"
            assert sum(int(row[k]) for k in ("final_momentum_ns", "terminal_metrics_ns", "boundary_ledger_ns")) <= int(row["advance_ns"])
        for rank in range(args.ranks):
            with (path / ("solver-rank-{}.csv".format(rank))).open() as stream:
                loops = list(csv.DictReader(stream))
            for row in loops:
                assert row["source_meta_sha256"] == source_hash
                assert row["attempt"] == "1" and row["corrector"] in ("1", "2")
                assert int(row["mg_copy_ns"]) <= int(row["mg_refill_ns"])
                assert row["dropped_loops"] == "0"
                assert row["linear_criterion_valid"] in ("0", "1")
                if row["invoked"] == "1":
                    assert row["linear_criterion_valid"] == "1", "successful invoked solve lost its criterion"
                rhs, atol, rtol, limit = (float(row[key]) for key in (
                    "linear_rhs_norm", "linear_atol", "linear_rtol", "linear_residual_limit"))
                assert all(math.isfinite(value) for value in (rhs, atol, rtol, limit))
                if row["linear_criterion_valid"] == "1":
                    assert rhs >= 0 and atol >= 0 and rtol >= 0 and (atol > 0 or rtol > 0)
                    assert math.isclose(limit, max(atol, rtol * rhs), rel_tol=8 * sys.float_info.epsilon)
                else:
                    assert (rhs, atol, rtol, limit) == (0, 0, 0, 0)
        observer = Path(__file__).resolve().parents[4] / "tools" / "v04_solver_observe.py"
        observed = subprocess.run([sys.executable, str(observer), str(path),
            "--output", str(root / (name + "-solver-observation.json"))],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=20)
        assert observed.returncode == 0, observed.stdout
        attribution = json.loads((root / (name + "-solver-observation.json")).read_text())
        assert attribution["complete"] and attribution["validated_steps"] == steps
        assert attribution["expected_ranks"] == list(range(args.ranks))
        assert attribution["schema"] == "HUNDUN_LOOP_OBSERVATION_V4"
        validator = Path(__file__).resolve().parents[4] / "tools" / "v04_evidence_validate.py"
        check = [sys.executable, str(validator), "runtime", str(path / "evidence.jsonl")]
        if source:
            generation = (source / "Restart" / "current").read_text().strip()
            check += ["--run-start-manifest", str(source / "Restart" / generation / "manifest.bin")]
        validated = subprocess.run(check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   universal_newlines=True, timeout=20)
        (root / (name + "-evidence.log")).write_text(validated.stdout)
        assert validated.returncode == 0, validated.stdout
        with (path / "health.csv").open() as stream:
            health = list(csv.DictReader(stream))
        stats = json.loads(sorted(path.glob("*.statistics.json"))[-1].read_text())
        return path, metadata(path / "RUN.meta"), health, stats
    source, source_meta, _, source_stats = run("source", 5)
    assert source_stats["sample_steps"] == 2, source_stats
    before = inventory(source)
    exact, exact_meta, exact_health, exact_stats = run("exact", 2, source)
    assert exact_meta["starting_sample_steps"] == "2"
    assert exact_stats["sample_steps"] == 4
    assert all(row["bdf_order"] == "1" for row in exact_health)
    exact_before = inventory(exact)
    recovered, recovered_meta, recovered_health, recovered_stats = run("recovered", 2, exact, True)
    print("epoch source_samples={} exact_samples={} method_start_samples={} method_end_samples={}".format(
        source_stats["sample_steps"], exact_stats["sample_steps"],
        recovered_meta["starting_sample_steps"], recovered_stats["sample_steps"]), flush=True)
    assert recovered_meta["starting_sample_steps"] == "0", "method recovery inherited old samples"
    assert recovered_stats["sample_steps"] == 0
    assert [row["bdf_order"] for row in recovered_health] == ["1", "1"]
    assert recovered_meta["statistics_epoch_start_step"] == "7"
    assert recovered_meta["statistics_sampling_start_step"] == "11"
    assert recovered_meta["statistics_reset_reason"] == "method_recovery"
    assert recovered_meta["statistics_discarded_samples"] == "4"
    assert int(recovered_meta["statistics_source_manifest_sha256"], 16) != 0
    resumed, resumed_meta, resumed_health, resumed_stats = run("resumed", 3, recovered)
    assert resumed_meta["statistics_epoch_start_step"] == "7"
    assert resumed_meta["statistics_sampling_start_step"] == "11"
    assert all(row["bdf_order"] == "1" for row in resumed_health)
    assert resumed_stats["sample_steps"] == 2, "exact segmentation lost the relative development window"
    validator = Path(__file__).resolve().parents[4] / "tools" / "v04_evidence_validate.py"
    generation = (exact / "Restart" / "current").read_text().strip()
    manifest = exact / "Restart" / generation / "manifest.bin"
    for key, value in (("policy", "require_compatible"), ("source_signature", 42),
                       ("source_format_version", 2)):
        rows = [json.loads(line) for line in (recovered / "evidence.jsonl").read_text().splitlines()]
        for row in rows:
            row["run_start"]["history"][key] = value
        damaged = root / ("bad-history-" + key + ".jsonl")
        damaged.write_text("".join(json.dumps(row) + "\n" for row in rows))
        rejected = subprocess.run([sys.executable, str(validator), "runtime", str(damaged),
                                   "--run-start-manifest", str(manifest)],
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  universal_newlines=True, timeout=20)
        assert rejected.returncode != 0, "tampered history anchor accepted: " + key
    _, override_meta, _, override_stats = run("override", 3, source, True, 1)
    assert override_meta["statistics_sampling_start_step"] == "7"
    assert override_stats["sample_steps"] == 2
    corrupt = root / "corrupt-source"
    shutil.copytree(str(source), str(corrupt))
    generation = (source / "Restart" / "current").read_text().strip()
    block = next((corrupt / "Restart" / generation).glob("rank-*"))
    payload = bytearray(block.read_bytes())
    payload[-1] ^= 1
    block.write_bytes(payload)
    rejected = subprocess.run([args.mpi, "-n", str(args.ranks), str(args.binary.resolve()),
        "--spec", str(spec.resolve()), "--case-root", str(case.resolve()),
        "--run-root", str(root / "corrupt-result"), "--steps", "1", "--visit-interval", "0",
        "--restart-root", str(corrupt / "Restart"), "--restart-method-recovery"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=40)
    (root / "corrupt-source.log").write_text(rejected.stdout)
    assert rejected.returncode != 0 and "initialize_status=" in rejected.stdout, \
        "method policy bypassed source integrity"
    assert inventory(source) == before, "source checkpoint or attachments were modified"
    assert inventory(exact) == exact_before, "exact-continuation source was modified by method recovery"
    print("PASS runner statistics epoch ranks={} source_readonly=true".format(args.ranks))


if __name__ == "__main__":
    main()
