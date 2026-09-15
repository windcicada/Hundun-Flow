#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Extract the fixed-dt, single-field COAST cylinder reference window.

COAST timing.dat stores cumulative CPU means per call. Convert those means
using the actual legacy-driver call multiplicity; keep root wall time separate.
The final output step is outside the measured interval. Python 3.6 compatible.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report(case, log, program, dt, warmup, steps, timing=None):
    control = (case / 'input.d').read_text()
    iterations = int(re.search(r'^\s*(\d+)\s*/itnum', control, re.M).group(1))
    fields = int(re.search(r'^\s*(\d+)\s*/number of stochastic fields',
                           control, re.M).group(1))
    if fields != 1:
        raise ValueError('this reference entry requires one mean field')
    screen = (case / 'screen').read_text()
    trace = log.read_text()
    if 'finished' not in screen or any(v in trace for v in
            ('CORE BAD VALUE', 'MPI_ABORT', 'VELOCITY MACH LIMIT', 'STOP ')):
        raise ValueError('reference requires normal completion and finite states')
    timing = timing or case / 'monitor/timing.dat'
    rows = [[float(v) for v in line.split()]
            for line in timing.read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    wall = [float(v) for v in re.findall(r'wall_step\[s\]\s+([0-9.]+)', screen)]
    end = warmup + steps
    if len(rows) != end + 1 or len(wall) != len(rows):
        raise ValueError('expected warmup, measurement and one final output step')
    if any(len(r) != 15 or not all(math.isfinite(v) for v in r) for r in rows):
        raise ValueError('invalid timing records')
    if any(r[0] != rows[0][0]+i for i, r in enumerate(rows)):
        raise ValueError('step records must be continuous with one run per log')
    if any(abs(r[2]/dt-1) > 1e-6 for r in rows):
        raise ValueError('reference changed dt beyond the FP32 reporting tolerance')
    # Input boundary once; post-flow NSCBC and wall stages once each.
    # TPDF includes scalar transport. These two quantities are nested.
    columns = {'momentum': (4, iterations), 'continuity': (5, iterations),
               'boundary': (6, 3), 'statistics': (7, 1), 'scalar': (8, 1),
               'mix': (9, 1), 'reactor': (10, 1), 'tpdf': (11, 1), 'spray': (12, 1)}
    cpu = {}
    for name, (column, calls) in columns.items():
        start_sum = rows[warmup-1][column]*warmup if warmup else 0.0
        cpu[name] = (rows[end-1][column]*end-start_sum)*calls/steps
    balances = [float(v) for v in re.findall(
        r'(?<!_)relative_residual=\s*([\d.E+\-]+)', screen)]
    if len(balances) != len(rows):
        raise ValueError('one domain mass balance record per step is required')
    solves = re.findall(r'压力后端\s+(\w+)\s+iter\s+(\d+)\s+res\s+([\d.E+\-]+)', screen)
    fallbacks = re.findall(r'pressure_time\[s\].*?fallback\s+([TF])\s+(\S+)', screen)
    if len(solves) != len(rows) or len(fallbacks) != len(rows):
        raise ValueError('one last-pressure-solve backend record per step is required')
    result = {
        'status': 'continuous fixed-dt reference; accuracy acceptance recorded separately',
        'warmup_steps': warmup, 'timed_steps': steps,
        'final_output_step_excluded': True,
        'mean_root_wall_s': sum(wall[warmup:end])/steps,
        'wall_display_resolution_s': 0.01,
        'root_wall_s': wall, 'root_cpu_module_s_per_step': cpu,
        'module_calls_per_step': {k: v[1] for k, v in columns.items()},
        'cpu_scope': 'root cumulative means converted to per-step CPU; TPDF includes scalar',
        'dt_min': min(r[2] for r in rows), 'dt_max': max(r[2] for r in rows),
        'cfl_max': max(r[-1] for r in rows),
        'mass_relative_max': max(balances), 'rows': rows,
        'mass_audit_phase': 'input boundary after PDF, before current flow correction',
        'last_pressure_solve_backend_counts': dict(Counter(s[0] for s in solves)),
        'last_pressure_solve_iterations': [int(s[1]) for s in solves],
        'last_pressure_solve_fallback_counts': dict(Counter(' '.join(s) for s in fallbacks)),
        'input_sha256': digest(case / 'input.d'),
        'boundary_sha256': digest(case / 'boundary_conditions.d'),
        'program_sha256': digest(program), 'trace_sha256': digest(log),
        'timing_sha256': digest(timing),
    }
    observed = [[float(v) for v in line.split()[1:]]
                for line in trace.splitlines() if line.startswith('COAST_COST ')]
    if observed:
        from ccost import LABELS
        if len(observed) != len(rows) or any(len(r) != 26 for r in observed):
            raise ValueError('one complete endpoint observation per step is required')
        if any(r[0] != step[0] or not all(math.isfinite(v) for v in r)
               for r,step in zip(observed,rows)):
            raise ValueError('endpoint observations must match timing steps')
        measured = observed[warmup:end]
        if any(r[1] <= 0 for r in measured):
            raise ValueError('observed wall window starts after the first output call')
        result['endpoint_observation'] = {
            'scope': 'max-rank interval between output calls, excluding output and observation; TPDF nests scalar',
            'mean_max_rank_wall_s': sum(r[1] for r in measured)/steps,
            'kernel_max_rank_mean_s': {name: sum(r[2+i] for r in measured)/steps
                                       for i,name in enumerate(LABELS)},
            'kernel_calls_per_step': {name: sorted(set(int(r[12+i]) for r in measured))
                                       for i,name in enumerate(LABELS)},
            'mass_relative_max': max(r[24] for r in measured),
            'cell_rho_dt_relative_max': max(r[25] for r in measured),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    parser.add_argument('--log', required=True, type=Path)
    parser.add_argument('--program', required=True, type=Path)
    parser.add_argument('--timing', type=Path,
                        help='captured timing table for exactly one run')
    parser.add_argument('--dt', required=True, type=float)
    parser.add_argument('--warmup', default=10, type=int)
    parser.add_argument('--steps', default=100, type=int)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.dt) or args.dt <= 0 or args.warmup < 0 or args.steps <= 0:
        parser.error('positive dt and measured steps, and nonnegative warmup are required')
    result = report(args.case, args.log, args.program, args.dt, args.warmup, args.steps,
                    args.timing)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('{:.6f} s/step (root wall), mass relative max {:.6g}'.format(
        result['mean_root_wall_s'], result['mass_relative_max']))


if __name__ == '__main__':
    main()
