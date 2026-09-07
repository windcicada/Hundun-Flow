#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Real public CLIs: local inversion exhaustion, common evidence and no commit.

The three-iteration fixture intentionally exhausts inversion. It is not the
production 200-iteration Re3900 configuration or a replay of step 7232.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

THERMO = """HUNDUN_THERMOPHYSICS_V1
temperature_bounds 273.15 6000
temperature_inversion 1e-12 3
closed_mass_newton 1e-12 32 0.2
species_count 1
species air
molecular_weight 28.850334
temperature_switch 1000
nasa7_low 3.5838100068 -7.2700635412e-4 1.67056387003e-6 -1.091801341e-10 -4.317787988e-13 -1050.5394088 3.1124135035
nasa7_high 3.1013370688 1.24138813631e-3 -4.1882038804e-7 6.641656204e-11 -3.9127843272e-15 -985.27467132 5.3560174057
transport_coast_native_air
end_species
end
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--app", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix="hundun-thermo-cli-")) / "audit"
    root.mkdir(parents=True, exist_ok=False)
    data = Path(__file__).resolve().parents[1] / "data"
    results = []
    for target, shift in ((0, -2), (1, 2)):
        case = root / ("case-{}".format(target))
        case.mkdir()
        model = json.loads((data / "case_minimal_valid.json").read_text())
        model["mesh"]["domain"] = {"lower": [-4, -2, -2], "upper": [4, 2, 2]}
        model["mesh"]["exact_cells"] = [16, 16, 16]
        model["mesh"]["minimum_spacing"] = [.5, .25, .25]
        # The public case contract restricts native-air transport to COAST axes.
        model["mesh"].update(kind="coast_runtime_axes_v1", axes_file="axes.dat",
                              base_spacing=[.5, .25, .25])
        axes = ["COAST_RUNTIME_AXES 1", "grid 16 16 16"]
        for axis, lower, spacing in (("x", -4, .5), ("y", -2, .25), ("z", -2, .25)):
            axes += [axis + " 17", " ".join(str(lower + i * spacing) for i in range(17))]
        (case / "axes.dat").write_text("\n".join(axes) + "\n")
        model["mesh"]["limits"]["max_memory_bytes_per_rank"] = 1073741824
        model["mesh"]["immersed_boundary"] = {"stl_file": "cube.stl", "fluid_side": "outside",
                                                "reconstruction_policy": "adaptive_order"}
        model["flow"]["pressure_reference"] = "closed_mass"
        for boundary in model["boundaries"].values():
            boundary["velocity"] = [0, 0, 0]
            boundary["pressure"] = 0
            for key in ("temperature", "backflow_temperature"):
                boundary[key] = 300.0
        for face in ("x_min", "x_max"):
            model["boundaries"][face]["flow_kind"] = "no_slip_wall"
        model["time"].update(control="fixed", initial_dt=.001, minimum_dt=.001, maximum_dt=.001)
        (case / "case.json").write_text(json.dumps(model))
        (case / "thermophysics.d").write_text(THERMO)
        lines = []
        for line in (data / "cylinder_ascii.stl").read_text().splitlines():
            parts = line.split()
            if parts and parts[0] == "vertex":
                line = "vertex {} {} {}".format(float(parts[1]) + shift, parts[2], parts[3])
            lines.append(line)
        (case / "cube.stl").write_text("\n".join(lines) + "\n")
        spec = root / ("statistics-{}.d".format(target))
        spec.write_text("HUNDUN_V04_LITERATURE_STATISTICS_V1\ndevelopment_steps 10\n"
                        "collection_end_step 20\ncheckpoint_interval 1\nrho_ref 1\nu_ref 10\n"
                        "diameter 2\nspan 4\ncylinder_center_x {}\nstation_x_over_d 0.6\nend\n".format(shift))
        source = root / ("seed-{}".format(target)) / "Restart"
        initialized = subprocess.run([args.mpi, "-n", "2", str(args.runner.resolve()),
            "--spec", str(spec), "--case-root", str(case), "--dry-plan"], env=os.environ.copy(),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=45)
        (root / ("uniform-{}.log".format(target))).write_text(initialized.stdout)
        assert initialized.returncode == 0, initialized.stdout
        if target == 0:
            # A closed, nonperiodic accepted step still has a coupled C2 and
            # exact component audits, even though it has no refinement loop.
            for name, binary in (("runner", args.runner), ("app", args.app)):
                normal = root / ("uniform-" + name)
                command = [args.mpi, "-n", "2", str(binary.resolve())]
                if name == "runner":
                    command += ["--spec", str(spec), "--case-root", str(case),
                                "--run-root", str(normal), "--steps", "1", "--visit-interval", "1"]
                else:
                    command += ["run", str(case), "--output", str(normal), "--steps", "1",
                                "--output-interval", "1", "--restart-interval", "1"]
                accepted = subprocess.run(command, env=os.environ.copy(), stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, universal_newlines=True, timeout=45)
                (root / ("uniform-" + name + ".log")).write_text(accepted.stdout)
                assert accepted.returncode == 0, accepted.stdout
                assert list(normal.rglob("current")), "accepted step did not persist its checkpoint"
                records = [json.loads(line) for line in (normal / "evidence.jsonl").read_text().splitlines()]
                assert len(records) == 1
                record = records[0]
                assert record["pressure_solve_contract"] == "continuity_energy_coupled", record
                assert record["pressure_energy_refinement_solve_calls"] == 0
                assert record["pressure_energy_refinement_termination"] == "component_residuals_converged"
                assert record["terminal_physical_audit"]["energy_tolerance"] > 0
                validator = Path(__file__).resolve().parents[4] / "tools" / "v04_evidence_validate.py"
                subprocess.run([sys.executable, str(validator), "runtime", str(normal / "evidence.jsonl")],
                               check=True, timeout=30)
        seeded = subprocess.run([args.mpi, "-n", "2", str(args.seed.resolve()),
            "--seed-cli-failure", str(case), str(source), str(target)],
            env=os.environ.copy(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True, timeout=45)
        (root / ("seed-{}.log".format(target))).write_text(seeded.stdout)
        assert seeded.returncode == 0, seeded.stdout
        # Manufacture zero-sample legacy attachments for this new V1 fixture.
        # Identities come from public dry-plan/seed outputs, not guessed hashes.
        dry_line = next(line for line in initialized.stdout.splitlines() if line.startswith("schema="))
        dry = dict(word.split("=", 1) for word in dry_line.split())
        seed_line = next(line for line in seeded.stdout.splitlines() if line.startswith("seed "))
        identity = dict(word.split("=", 1) for word in seed_line.split()[1:])
        assert identity["plan"] == dry["product"]
        stem = "step-00000000000000000001"
        header = ["spec_fingerprint " + dry["spec"], "plan " + identity["plan"],
                  "schema " + identity["schema"], "step 1", "time 0.001"]
        accumulator = ["HUNDUN_V04_THIN_DOMAIN_ACCUMULATOR_V1"] + header + [
            "sample_steps 0", "profile_size 96", "centerline_size 48"]
        accumulator += ["p {} 0".format(i) for i in range(96)]
        accumulator += ["c {} 0".format(i) for i in range(48)] + ["end"]
        (source.parent / (stem + ".accumulator")).write_text("\n".join(accumulator) + "\n")
        statistics = {"schema": "HUNDUN_V04_THIN_DOMAIN_STATISTICS_V1",
            "spec_fingerprint": int(dry["spec"]), "plan": int(identity["plan"]),
            "field_schema": int(identity["schema"]), "snapshot_step": 1,
            "snapshot_time": .001, "sample_steps": 0, "profiles": [], "centerline": []}
        (source.parent / (stem + ".statistics.json")).write_text(json.dumps(statistics))
        marker = ["HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1"] + header + [
            "restart_generation " + (source / "current").read_text().strip(),
            "statistics " + stem + ".statistics.json", "accumulator " + stem + ".accumulator", "end"]
        (source.parent / (stem + ".complete")).write_text("\n".join(marker) + "\n")
        original = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in source.parent.rglob("*") if p.is_file()}
        for name, binary in (("runner", args.runner), ("app", args.app)):
            run = root / ("{}-rank-{}".format(name, target))
            command = [args.mpi, "-n", "2", str(binary.resolve())]
            if name == "runner":
                command += ["--spec", str(spec), "--case-root", str(case), "--run-root", str(run),
                            "--restart-root", str(source), "--steps", "1", "--visit-interval", "1", "--observe-performance",
                            "--trace-cell", "{},7,{}".format(5 if target == 0 else 13, 7 if target == 0 else 8),
                            "--trace-first", "2", "--trace-last", "2"]
            else:
                command += ["run", str(case), "--output", str(run), "--steps", "1",
                            "--output-interval", "1", "--restart-interval", "1",
                            "--restart", str(source)]
            result = subprocess.run(command, env=os.environ.copy(), stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, universal_newlines=True, timeout=45)
            (root / (run.name + ".log")).write_text(result.stdout)
            results.append({"name": name, "target": target, "command": command,
                            "returncode": result.returncode})
            (root / "results.json").write_text(json.dumps(results, indent=2))
            assert result.returncode != 0 and "COMPLETED" not in result.stdout, result.stdout
            contexts = re.findall(r"^numerical_failure_context_v1 .*$", result.stdout, re.M)
            assert len(contexts) == 1, result.stdout
            context = dict(word.split("=", 1) for word in contexts[0].split()[1:])
            assert context["status"] == "5/804" and context["stage"] == "15", context
            assert int(context["rank"]) == target, context
            assert context["inversion_iterations"] == "3", context
            assert "step_failure_completion_v1" in result.stdout, result.stdout
            assert not list(run.rglob("*.complete")) and not list(run.rglob("current"))
            assert original == {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in source.parent.rglob("*") if p.is_file()}
            if name == "runner":
                with (run / ("cell-trace-rank-{}.csv".format(target))).open() as stream:
                    trace = list(csv.DictReader(stream))
                assert {11, 12}.issubset({int(row["stage"]) for row in trace}), trace
                for row in trace:
                    if int(row["stage"]) >= 12:
                        assert row["active"] == "0" and float(row["h"]) == float(row["h_n"]), row
            print("{} thermal failure rank={} returned consistently".format(name, target), flush=True)


if __name__ == "__main__":
    main()
