#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Real CLI MG observation: independent layout, step/loop accounting, same physics."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


PHASES = ('pre_smooth', 'residual', 'restriction', 'prolongation',
          'post_smooth', 'terminal', 'direct_mpi')
TOTALS = ('attempts', 'successes', 'failures', 'apply_ns', 'reduction_ns')
LEVELS = tuple(p + s for p in PHASES for s in ('_calls', '_ns')) + (
    'halo_wait_ns', 'halo_control_ns', 'halo_control_calls')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--mpi', default='mpirun')
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--fgmres-recovery', action='store_true')
    args = parser.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix='hundun-mg-cli-')) / 'audit'
    root.mkdir(parents=True, exist_ok=False)
    case = root / 'case'
    case.mkdir()
    data = Path(__file__).resolve().parents[1] / 'data'
    model = json.loads((data / 'case_minimal_valid.json').read_text())
    model['mesh']['domain'] = {'lower': [-2, -2, -2], 'upper': [2, 2, 2]}
    model['mesh']['exact_cells'] = [16, 16, 16]
    model['mesh']['minimum_spacing'] = [.25, .25, .25]
    model['mesh']['limits']['max_memory_bytes_per_rank'] = 1073741824
    model['mesh']['immersed_boundary'] = {'stl_file': 'cylinder_ascii.stl', 'fluid_side': 'outside'}
    model['flow']['pressure_reference'] = 'closed_mass'
    for boundary in model['boundaries'].values():
        boundary.update(flow_kind='periodic', velocity=[0, 0, 0], pressure=0)
    (case / 'case.json').write_text(json.dumps(model))
    for name in ('thermophysics.d', 'cylinder_ascii.stl'):
        shutil.copyfile(str(data / name), str(case / name))
    spec = root / 'statistics.d'
    spec.write_text('HUNDUN_V04_LITERATURE_STATISTICS_V1\n'
                    'development_steps 10\ncollection_end_step 20\ncheckpoint_interval 1\n'
                    'rho_ref 1\nu_ref 1\ndiameter 1\nspan 4\ncylinder_center_x 0\n'
                    'station_x_over_d 1.5\nend\n')

    def run(name, enabled, recovery=False, case_root=case):
        path = root / name
        command = [args.mpi, '-n', str(args.ranks), str(args.binary.resolve()),
                   '--spec', str(spec), '--case-root', str(case_root), '--run-root', str(path),
                   '--steps', '2', '--visit-interval', '0', '--observe-performance']
        if enabled:
            command.append('--observe-mg-cost')
        if recovery:
            command.append('--observe-fgmres-recovery')
        result = subprocess.run(command, env=os.environ.copy(), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True, timeout=60)
        (root / (name + '.log')).write_text(result.stdout)
        assert result.returncode == 0, result.stdout
        assert 'COMPLETED steps=2 final_step=2' in result.stdout, result.stdout
        return path

    off = run('off', False)
    on = run('on', True, args.fgmres_recovery)
    metadata = dict(line.split(' ', 1) for line in (on / 'RUN.meta').read_text().splitlines()[1:-1])
    assert metadata['observation_schema'] == ('6' if args.fgmres_recovery else '5')
    assert metadata['observe_mg_cost'] == '1'
    assert metadata['expected_ranks'] == str(args.ranks)
    assert metadata['mg_layout_sha256'] == digest(on / 'mg-layout.meta')
    assert not (off / 'mg-layout.meta').exists() and not list(off.glob('mg-rank-*.csv'))
    off_metadata = dict(line.split(' ', 1) for line in (off / 'RUN.meta').read_text().splitlines()[1:-1])
    assert off_metadata['observation_schema'] == '4'
    layouts = (on / 'mg-layout.meta').read_text().splitlines()
    assert layouts[0] == 'HUNDUN_MG_LAYOUT_V1' and layouts[-1] == 'end'
    assert layouts[1] == 'expected_ranks {}'.format(args.ranks)
    assert layouts[2] == 'maximum_levels 32'
    counts, geometry = {}, {}
    for line in layouts[3:-1]:
        kind, *words = line.split()
        values = list(map(int, words))
        if kind == 'rank':
            rank, count = values
            assert rank not in counts and 0 < count <= 32
            counts[rank] = count
        else:
            assert kind == 'level' and len(values) == 10
            rank, level = values[:2]
            assert (rank, level) not in geometry and all(v > 0 for v in values[2:8])
            geometry[rank, level] = values[2:]
    assert set(counts) == set(range(args.ranks))
    assert set(geometry) == {(r, l) for r in counts for l in range(counts[r])}
    assert all(geometry[r, 0][:3] == [16, 16, 16] for r in counts)
    assert sum(geometry[r, 0][3] * geometry[r, 0][4] * geometry[r, 0][5] for r in counts) == 16**3
    source = digest(on / 'RUN.meta')
    for rank in range(args.ranks):
        loops = rows(on / 'solver-rank-{}.csv'.format(rank))
        detail = rows(on / 'mg-rank-{}.csv'.format(rank))
        assert all(r['source_meta_sha256'] == source for r in loops + detail)
        keys = [(int(r['step']), int(r['rank']), int(r['level'])) for r in detail]
        assert len(keys) == len(set(keys))
        assert set(keys) == {(s, rank, l) for s in (1, 2) for l in range(-1, counts[rank])}
        assert all(r['complete'] == '1' and r['initialized'] == '1' for r in detail)
        for step in (1, 2):
            step_loops = [r for r in loops if int(r['step']) == step]
            step_detail = [r for r in detail if int(r['step']) == step]
            total = next(r for r in step_detail if r['level'] == '-1')
            levels = [r for r in step_detail if r['level'] != '-1']
            assert step_loops and all(r['mg_enabled'] == '1' and r['mg_complete'] == '1' for r in step_loops)
            assert int(total['attempts']) == int(total['successes']) + int(total['failures'])
            for key in TOTALS:
                assert sum(int(r['mg_' + key]) for r in step_loops) == int(total[key])
                assert all(int(r[key]) == 0 for r in levels)
            for key in LEVELS:
                assert int(total[key]) == 0
                assert sum(int(r['mg_' + key]) for r in step_loops) == sum(int(r[key]) for r in levels)
            assert int(total['visits']) == 0
            for phase in ('pre_smooth', 'post_smooth'):
                finest = next(r for r in levels if r['level'] == '0')
                assert sum(int(r['mg_finest_' + phase + '_ns']) for r in step_loops) == int(finest[phase + '_ns'])
    def checkpoint(path):
        generation = (path / 'Restart/current').read_text().strip()
        files = sorted((path / 'Restart' / generation).glob('rank-*.bin'))
        assert len(files) == args.ranks
        return {p.name: digest(p) for p in files}
    assert checkpoint(off) == checkpoint(on), 'observation changed committed checkpoint bytes'
    observer = Path(__file__).resolve().parents[4] / 'tools/v04_solver_observe.py'
    result = subprocess.run([sys.executable, str(observer), str(on), '--output', str(root / 'observed.json'),
                            '--details-output', str(root / 'details.jsonl')],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=15)
    assert result.returncode == 0, result.stdout
    observed = json.loads((root / 'observed.json').read_text())
    assert observed['complete'] and observed['validated_steps'] == 2
    assert len(observed['mg']['levels']) == sum(counts.values())
    attempts = [r['attempts'] for r in observed['mg']['totals_by_rank']]
    if args.fgmres_recovery:
        alone = run('recovery-only', False, True)
        assert checkpoint(alone) == checkpoint(off)
        assert not (alone / 'mg-layout.meta').exists() and not list(alone.glob('mg-rank-*.csv'))
        a_keys = ('initial_residual_applies', 'arnoldi_applies', 'unsafe_residual_applies',
                  'interior_residual_applies', 'cycle_residual_applies')
        for path in (on, alone):
            meta = dict(line.split(' ', 1) for line in (path / 'RUN.meta').read_text().splitlines()[1:-1])
            assert meta['observe_fgmres_recovery'] == '1'
            assert meta['pressure_linear_algorithm'] == 'fgmres' and meta['observation_schema'] == '6'
            seen = 0
            for rank in range(args.ranks):
                for row in rows(path / 'solver-rank-{}.csv'.format(rank)):
                    assert int(row['fgmres_available']) == int(row['invoked'])
                    assert sum(int(row['fgmres_' + k]) for k in a_keys) == int(row['A_calls'])
                    assert int(row['fgmres_unsafe_restarts']) + int(row['fgmres_happy_restarts']) == int(row['fgmres_norm_restarts'])
                    seen += int(row['fgmres_arnoldi_applies'])
            assert seen > 0, 'fixture must execute actual Arnoldi work'
            result = subprocess.run([sys.executable, str(observer), str(path)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=15)
            assert result.returncode == 0, result.stderr
            report = json.loads(result.stdout)
            assert report['complete'] and report['schema'] == 'HUNDUN_LOOP_OBSERVATION_V6'
            assert ('mg' in report) == (path == on)
            assert all(r['fgmres_recovery']['available'] for r in report['loops'])
            print('FGMRES CLI {} ranks={} Arnoldi={} complete/checkpoint parity PASS'.format(path.name, args.ranks, seen))
        assert all(not any(k.startswith('fgmres_') for k in row)
                   for row in rows(off / 'solver-rank-0.csv'))
        # Exercise the actual non-FGMRES producer, not just a synthetic CSV.
        # The public case parser selects the matching unit-linear MG contract.
        bicg_case = root / 'bicg-case'
        shutil.copytree(str(case), str(bicg_case))
        bicg_model = json.loads((bicg_case / 'case.json').read_text())
        bicg_model['solver']['pressure_linear']['algorithm'] = 'bicgstab'
        bicg_model['solver']['pressure_linear']['krylov_restart'] = 0
        (bicg_case / 'case.json').write_text(json.dumps(bicg_model))
        bicg_off = run('bicg-off', False, case_root=bicg_case)
        bicg_on = run('bicg-on', False, True, bicg_case)
        assert checkpoint(bicg_off) == checkpoint(bicg_on)
        work = 0
        for rank in range(args.ranks):
            for row in rows(bicg_on / 'solver-rank-{}.csv'.format(rank)):
                assert all(int(v) == 0 for k, v in row.items() if k.startswith('fgmres_'))
                work += int(row['A_calls']) + int(row['M_calls'])
        assert work > 0, 'non-FGMRES fixture must execute actual linear work'
        result = subprocess.run([sys.executable, str(observer), str(bicg_on)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=15)
        assert result.returncode == 0, result.stderr
        report = json.loads(result.stdout)
        assert report['complete'] and report['pressure_linear_algorithm'] == 'bicgstab'
        assert all(r['fgmres_recovery'] == {'available': False, 'counts_by_rank': None}
                   for r in report['loops'])
        print('BiCGStab CLI ranks={} A+M={} unavailable recovery/checkpoint parity PASS'.format(args.ranks, work))
    # Both flags alter communication branches; inconsistent flags must reject
    # at the cold entrance, before creating output or running a time step.
    base = [str(args.binary.resolve()), '--spec', str(spec), '--case-root', str(case), '--steps', '1']
    contracts = [
            ('needs-performance', ['--observe-mg-cost'], None),
            ('duplicate-mg', ['--observe-performance', '--observe-mg-cost', '--observe-mg-cost'], None),
            ('mixed-mg', ['--observe-performance', '--observe-mg-cost'], ['--observe-performance']),
            ('mixed-performance', ['--observe-performance'], [])]
    if args.fgmres_recovery:
        contracts += [
            ('needs-performance-fgmres', ['--observe-fgmres-recovery'], None),
            ('duplicate-fgmres', ['--observe-performance', '--observe-fgmres-recovery', '--observe-fgmres-recovery'], None),
            ('mixed-fgmres', ['--observe-performance', '--observe-fgmres-recovery'], ['--observe-performance'])]
    for name, left, right in contracts:
        if right is not None and args.ranks < 2:
            continue
        path = root / name
        common = base + ['--run-root', str(path)]
        command = [args.mpi, '-n', str(1 if right is not None else args.ranks)] + common + left
        if right is not None:
            command += [':', '-n', str(args.ranks - 1)] + common + right
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=15)
        (root / (name + '.log')).write_text(result.stdout)
        assert result.returncode == (2 if right is None else 3), result.stdout
        assert not path.exists() and 'COMPLETED' not in result.stdout, result.stdout
        if right is not None:
            assert 'spec_consensus_failure' in result.stdout, result.stdout
    print('MG CLI ranks={} steps=2 apply_attempts={} layout/loop/level identities, observer, '
          'checkpoint parity and cold flag contract PASS output={}'.format(args.ranks, attempts, root))


if __name__ == '__main__':
    main()
