#!/usr/bin/env python3
"""Audit a same-binary, same-state pressure-solver diagnostic comparison."""
import argparse
import copy
import hashlib
import json
import math
import os
import re

from timing import read_run


def digest(path):
    result = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def case_snapshot(path):
    with open(os.path.join(path, 'case.json')) as stream:
        original = json.load(stream)
    normalized = copy.deepcopy(original)
    pressure = normalized['solver']['pressure_linear']
    method = {key: pressure.pop(key) for key in
              ('algorithm', 'krylov_restart', 'mg_correction_scaling')
              if key in pressure}
    assets = {}
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for name in sorted(files):
            filename = os.path.join(root, name)
            relative = os.path.relpath(filename, path)
            if relative != 'case.json':
                assets[relative] = digest(filename)
    return normalized, {'path': os.path.abspath(path), 'method': method,
                        'case_json_sha256': digest(os.path.join(path, 'case.json')),
                        'assets_sha256': assets}


def pressure_trace(path, steps):
    matrices, solves = [], []
    blocks = 0
    with open(path + '.log') as stream:
        for line in stream:
            tag = line.split(' ', 1)[0]
            if tag == 'cold_momentum_policy':
                blocks += 1
            if blocks > steps:
                break
            fields = dict(re.findall(r'(\w+)=([^\s]+)', line))
            if tag == 'cold_matrix':
                matrices.append(fields)
            if tag == 'cold_solve':
                if fields['status'] != '0/0':
                    raise ValueError('pressure solve failed: ' + path)
                solves.append(fields)
    if blocks < steps or not solves or len(solves) != len(matrices):
        raise ValueError('incomplete matrix/solver diagnostic: ' + path)
    if matrices[0]['outer'] != '0':
        raise ValueError('first pressure matrix is absent')
    names = ('solve', 'operator_s', 'preconditioner_s', 'arnoldi_dot_s',
             'arnoldi_reduce_s', 'arnoldi_update_s', 'operator_calls',
             'preconditioner_calls', 'reduction_calls', 'iterations')
    values = {key: [float(row[key]) for row in solves] for key in names}
    if any(not math.isfinite(v) or v < 0 for series in values.values() for v in series):
        raise ValueError('invalid solver counter')
    return {'log_sha256': digest(path + '.log'), 'first_matrix': matrices[0],
            'first_solve': solves[0], 'solve_count': len(solves),
            'per_step': {key: sum(series) / steps for key, series in values.items()},
            'maximum_continuity_residual': max(float(row['continuity_max'])
                                               for row in solves)}


def compare(paths, cases, steps, dt):
    snapshots = [case_snapshot(case) for case in cases]
    if snapshots[0][0] != snapshots[1][0]:
        raise ValueError('case physics, tolerances or non-admitted solver settings differ')
    if snapshots[0][1]['assets_sha256'] != snapshots[1][1]['assets_sha256']:
        raise ValueError('case asset bytes differ')
    runs = [read_run(path, 0, steps, dt) for path in paths]
    if runs[0]['candidate'] != runs[1]['candidate']:
        raise ValueError('runtime executable/build identities differ')
    identities = [copy.deepcopy(run['identity']) for run in runs]
    for identity in identities:
        identity.pop('case')
        history = identity['run_start'].get('history', {})
        source = history.pop('transport_source_case', None)
        if source is not None and source != runs[0]['identity']['case']:
            raise ValueError('solver migration refers to another source case')
    if identities[0] != identities[1]:
        raise ValueError('MPI plan, geometry, initial state or method history differs')
    traces = [pressure_trace(path, steps) for path in paths]
    for key in ('matrix', 'rhs'):
        if traces[0]['first_matrix'][key] != traces[1]['first_matrix'][key]:
            raise ValueError('first pressure ' + key + ' differs')
    return {'schema': 'hundun_solver_comparison_v1',
            'scope': 'diagnostic: first matrix/RHS identity; accepted fixed-dt steps',
            'case_snapshot_scope': 'current files; runtime source-case migration binds plan compatibility',
            'timer_scope': 'substage max-rank sums; whole-step wall clock reported separately',
            'formal_window': False, 'profiling_enabled': True,
            'cases': [snapshot[1] for snapshot in snapshots], 'runs': runs,
            'pressure': traces, 'ratio': runs[1]['mean_s'] / runs[0]['mean_s']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline')
    parser.add_argument('candidate')
    parser.add_argument('--cases', nargs=2, required=True)
    parser.add_argument('--steps', type=int, default=3)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.steps < 1 or not math.isfinite(args.dt) or args.dt <= 0:
        parser.error('steps >= 1 and finite dt > 0 required')
    try:
        report = compare((args.baseline, args.candidate), args.cases, args.steps, args.dt)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, str(error) + '\n')
    with open(args.output, 'w') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('diagnostic {:.6f} -> {:.6f} s/step; ratio {:.6f}'.format(
        report['runs'][0]['mean_s'], report['runs'][1]['mean_s'], report['ratio']))


if __name__ == '__main__':
    main()
