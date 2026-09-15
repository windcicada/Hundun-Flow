#!/usr/bin/env python3
"""Compare native evidence for two runs sharing one case, state and fixed dt."""
import argparse
import hashlib
import json
import math
import os
import re
import statistics


def modules(path, warmup, steps):
    log = path + ".log"
    if not os.path.exists(log):
        return {}
    fields = {
        "cold_momentum_complete": {"momentum_assembly": "source_assembly_s",
                                   "momentum_closure": "closure_s",
                                   "momentum_solve_setup": "solve_setup_s"},
        "cold_solve": {"pressure_assembly": "assembly", "pressure_setup": "setup",
                       "pressure_solve": "solve", "pressure_audit": "audit"},
        "cold_species": {"species": "solve_s"},
        "cold_energy_audit": {"energy_assembly": "assembly_audit_s"},
        "cold_enthalpy_EOS": {"enthalpy_correction": "total_s"},
        "cold_terminal_audit": {"terminal_audit": "seconds"},
        "cold_final_balance": {"boundary_balance": "seconds"},
        "cold_final_rates": {"endpoint_rates": "seconds"},
    }
    blocks = []
    with open(log) as stream:
        for line in stream:
            tag = line.split(" ", 1)[0]
            if tag == "cold_momentum_policy":
                blocks.append({})
            if tag not in fields or not blocks:
                continue
            values = dict(re.findall(r"(\w+)=([^\s]+)", line))
            for name, field in fields[tag].items():
                value = float(values[field])
                if not math.isfinite(value) or value < 0:
                    raise ValueError("{}: invalid module timer".format(log))
                blocks[-1][name] = blocks[-1].get(name, 0.0) + value
    if not blocks:
        return {}
    if len(blocks) < warmup + steps:
        raise ValueError("{}: incomplete module log".format(log))
    measured = blocks[warmup:warmup + steps]
    names = sorted(measured[0])
    if any(sorted(block) != names for block in measured):
        raise ValueError("{}: inconsistent module timers".format(log))
    # These are sums of each emitted substage timer, distinct from the
    # authoritative whole-step maximum-rank wall clock.
    return {name: statistics.mean(block[name] for block in measured) for name in names}


def read_run(path, warmup, steps, dt):
    evidence = os.path.join(path, "evidence.jsonl")
    with open(evidence, "rb") as stream:
        data = stream.read()
    rows = [json.loads(line) for line in data.splitlines() if line.strip()]
    if len(rows) < warmup + steps:
        raise ValueError("{}: need {} accepted steps, found {}".format(
            path, warmup + steps, len(rows)))
    rows = rows[:warmup + steps]
    first = rows[0]
    identity = {key: first[key] for key in ("case", "stl", "cpu_plan", "run_start")}
    for index, row in enumerate(rows):
        if any(row[key] != identity[key] for key in identity):
            raise ValueError("{}: run identity changed".format(path))
        actual_dt = row["time"] - row["previous_committed_time"]
        # Subtraction of absolute physical times carries their FP64 ULP.
        tolerance = max(abs(dt) * 1e-10,
                        8 * 2.220446049250313e-16 * abs(row["time"]))
        if abs(actual_dt - dt) > tolerance or row["retry"]:
            raise ValueError("{}: step {} changed dt or retried".format(path, row["step"]))
        if row["step"] != first["step"] + index:
            raise ValueError("{}: discontinuous step sequence".format(path))
        if index and row["previous_committed_time"] != rows[index - 1]["time"]:
            raise ValueError("{}: discontinuous physical time".format(path))
        audit = row["terminal_physical_audit"]
        if not audit["present"] or any(
                not math.isfinite(audit[key + "_residual"]) or
                audit[key + "_residual"] > audit[key + "_tolerance"]
                for key in ("eos", "continuity", "energy", "closed_mass", "gauge")):
            raise ValueError("{}: step {} physical audit".format(path, row["step"]))
    measured = rows[warmup:]
    seconds = [row["max_rank_step_ns"] / 1e9 for row in measured]
    stages = {}
    for row in measured:
        for stage in row["stages"]:
            stages.setdefault(str(stage["id"]), []).append(stage["max_ns"] / 1e9)
    return {"path": os.path.abspath(path), "identity": identity,
            "candidate": first["candidate_identity"],
            "evidence_sha256": hashlib.sha256(data).hexdigest(),
            "warmup": warmup, "measured_steps": steps, "dt": dt,
            "seconds": seconds, "mean_s": statistics.mean(seconds),
            "median_s": statistics.median(seconds),
            "minimum_s": min(seconds), "maximum_s": max(seconds),
            "stage_max_rank_mean_s": {key: statistics.mean(value)
                                      for key, value in stages.items()},
            "module_timer_mean_s": modules(path, warmup, steps),
            "outer_iterations": [row.get("cold", {}).get("outer_iterations")
                                 for row in measured],
            "linear_iterations": [row["linear_iterations"] for row in measured]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--warmup", type=int, default=10)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.warmup < 0 or args.steps < 1 or not math.isfinite(args.dt) or args.dt <= 0:
        parser.error("warmup >= 0, steps >= 1 and finite dt > 0 required")
    try:
        runs = [read_run(path, args.warmup, args.steps, args.dt)
                for path in (args.baseline, args.candidate)]
        if runs[0]["identity"] != runs[1]["identity"]:
            raise ValueError("case, MPI plan or starting state differs")
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, str(error) + "\n")
    ratio = runs[1]["mean_s"] / runs[0]["mean_s"]
    report = {"schema": "hundun_timing_v1", "runs": runs, "ratio": ratio,
              "within_five_percent": ratio <= 1.05,
              "formal_window": args.warmup >= 10 and args.steps >= 100}
    with open(args.output, "w") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print("{:.6f} -> {:.6f} s/step, ratio {:.6f}, formal_window={}".format(
        runs[0]["mean_s"], runs[1]["mean_s"], ratio, report["formal_window"]))


if __name__ == "__main__":
    main()
