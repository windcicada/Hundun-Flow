#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Real CLI end-of-run libc failures, after a durable accepted checkpoint."""
import argparse
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--mpi', default='mpirun')
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--first-only', action='store_true')
    parser.add_argument('--mg-only', action='store_true')
    parser.add_argument('--fgmres-only', action='store_true')
    args = parser.parse_args()
    output = args.output or Path(tempfile.mkdtemp(prefix='hundun-log-close-')) / 'audit'
    output.mkdir(parents=True, exist_ok=False)
    case = output / 'case'
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
        boundary['flow_kind'] = 'periodic'
        boundary['velocity'] = [0, 0, 0]
        boundary['pressure'] = 0
    (case / 'case.json').write_text(json.dumps(model))
    for name in ('thermophysics.d', 'cylinder_ascii.stl'):
        shutil.copyfile(str(data / name), str(case / name))
    spec = output / 'statistics.d'
    spec.write_text('HUNDUN_V04_LITERATURE_STATISTICS_V1\n'
                    'development_steps 10\ncollection_end_step 20\n'
                    'checkpoint_interval 1\nrho_ref 1\nu_ref 1\ndiameter 1\n'
                    'span 4\ncylinder_center_x 0\nstation_x_over_d 1.5\nend\n')
    results = []
    baseline_payload = None

    def run(name, stream, target, observe, flush_error=0, close_error=0, flush_at=1):
        nonlocal baseline_payload
        root = output / name
        prefix = output / (name + '-exit')
        env = dict(os.environ, LD_PRELOAD=str(args.probe.resolve()),
                   HUNDUN_LOG_TARGET=stream, HUNDUN_LOG_RANK=str(target),
                   HUNDUN_LOG_ROOT=str(root.resolve()),
                   HUNDUN_LOG_EXIT_PREFIX=str(prefix.resolve()))
        if flush_error:
            env.update(HUNDUN_LOG_FLUSH_ERRNO=str(flush_error), HUNDUN_LOG_FLUSH_AT=str(flush_at))
        if close_error:
            env['HUNDUN_LOG_CLOSE_ERRNO'] = str(close_error)
        command = [args.mpi, '-n', str(args.ranks), str(args.binary.resolve()),
                   '--spec', str(spec.resolve()), '--case-root', str(case.resolve()),
                   '--run-root', str(root.resolve()), '--steps', '1', '--visit-interval', '0']
        if observe:
            command.append('--observe-performance')
            if args.mg_only:
                command.append('--observe-mg-cost')
            if args.fgmres_only:
                command.append('--observe-fgmres-recovery')
        result = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True, timeout=45)
        (output / (name + '.log')).write_text(result.stdout)
        exits = [list(map(int, Path(str(prefix) + '-{}.txt'.format(rank)).read_text().split()))
                 for rank in range(args.ranks)]
        row = {'name': name, 'returncode': result.returncode, 'exits': exits, 'command': command}
        results.append(row)
        (output / 'results.json').write_text(json.dumps(results, indent=2))
        generation = (root / 'Restart/current').read_text().strip()
        payloads = sorted((root / 'Restart' / generation).glob('rank-*.bin'))
        assert len(payloads) == args.ranks, 'missing committed rank payload'
        payload = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in payloads}
        if baseline_payload is None:
            baseline_payload = payload
        assert payload == baseline_payload, 'I/O completion changed accepted checkpoint bytes'
        assert (root / 'step-00000000000000000001.complete').is_file()
        observation_stream = stream == 'performance.csv' or stream.startswith(('solver-rank-', 'mg-rank-'))
        active_failure = (flush_error or close_error) and (observe or not observation_stream)
        if active_failure:
            assert exits[target][3] == bool(flush_error), row
            assert exits[target][4] == bool(close_error), row
            assert result.returncode == 6 and all(value[0] == 6 for value in exits), row
            assert 'COMPLETED steps=' not in result.stdout, result.stdout
            assert 'log_completion_failure' in result.stdout, result.stdout
            if flush_error and close_error:
                assert 'operation=flush errno={}'.format(flush_error) in result.stdout, result.stdout
        else:
            assert result.returncode == 0 and all(value[0] == 0 for value in exits), row
            assert 'COMPLETED steps=1 final_step=1' in result.stdout, result.stdout
            assert all(value[3:] == [0, 0] for value in exits), row
        print('{} returncode={} rank_exits={}'.format(name, result.returncode, exits), flush=True)
        return exits[target][1]

    run('baseline-off', 'force.csv', 0, False)
    if not args.mg_only and not args.fgmres_only:
        run('force-close-eio', 'force.csv', 0, False, close_error=errno.EIO)
    if args.first_only:
        return
    disabled = ('solver-rank-0.csv' if args.fgmres_only else
                'mg-rank-0.csv' if args.mg_only else 'performance.csv')
    run('disabled-observer', disabled, 0, False, close_error=errno.ENOSPC)
    streams = (('solver-rank-0.csv', 'solver-rank-{}.csv'.format(args.ranks - 1)) if args.fgmres_only else
               ('mg-rank-0.csv', 'mg-rank-{}.csv'.format(args.ranks - 1)) if args.mg_only else
               ('force.csv', 'health.csv', 'conservation.csv', 'probe.csv',
                'performance.csv', 'solver-rank-{}.csv'.format(args.ranks - 1)))
    for index, stream in enumerate(streams):
        target = int(stream.split('-')[-1].split('.')[0]) if '-rank-' in stream else 0
        count = run('baseline-{}'.format(index), stream, target, True)
        # The baseline's last libc flush belongs to filebuf close. The repaired
        # code may add an earlier explicit flush; select that final-stage call.
        assert count >= 1
        at = max(1, count - 1)
        for error in (errno.EIO, errno.ENOSPC):
            run('close-{}-{}'.format(index, error), stream, target, True, close_error=error)
            run('flush-{}-{}'.format(index, error), stream, target, True, flush_error=error, flush_at=at)
        run('first-error-{}'.format(index), stream, target, True,
            flush_error=errno.EIO, close_error=errno.ENOSPC, flush_at=at)


if __name__ == '__main__':
    main()
