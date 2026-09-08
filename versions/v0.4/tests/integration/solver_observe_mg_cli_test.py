#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Public observer CLI: worked MG ledger (synthetic ns, not CFD measurements)."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def write(path, records):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def fixture(root):
    layout = ('HUNDUN_MG_LAYOUT_V1\nexpected_ranks 2\nmaximum_levels 32\n'
              'rank 0 2\nlevel 0 0 4 4 4 2 4 4 0 0\nlevel 0 1 2 2 2 2 2 2 0 0\n'
              'rank 1 2\nlevel 1 0 4 4 4 2 4 4 0 0\nlevel 1 1 2 2 2 2 2 2 0 0\nend\n')
    (root / 'mg-layout.meta').write_text(layout)
    metadata = ('HUNDUN_V04_THIN_DOMAIN_RUN_V1\nstarting_step 6\nrequested_steps 2\n'
                'expected_ranks 2\nobservation_schema 5\nobserve_performance 1\n'
                'observe_mg_cost 1\nmg_layout_sha256 {}\nend\n').format(
                    hashlib.sha256(layout.encode()).hexdigest())
    (root / 'RUN.meta').write_text(metadata)
    source = hashlib.sha256(metadata.encode()).hexdigest()
    loop = dict(attempt=1, scalar_coupling_sweep=1, corrector=2, refinement=1,
                kind=1, dropped_loops=0, invoked=1, iterations=2, A_calls=3,
                M_calls=2, A_ns=30, M_ns=120, baseline_candidates=1,
                extrapolated_candidates=0, ladder_candidates=0, incomplete_candidates=0,
                linear_initial=1, linear_final=0.01, globalization_valid=0,
                linear_criterion_valid=1, linear_rhs_norm=1, linear_atol=0.1,
                linear_rtol=0, linear_residual_limit=0.1,
                mg_enabled=1, mg_complete=1, mg_attempts=2, mg_successes=2, mg_failures=0,
                mg_apply_ns=100, mg_reduction_ns=5, mg_halo_wait_ns=5,
                mg_halo_control_ns=5, mg_halo_control_calls=5,
                mg_pre_smooth_calls=2, mg_pre_smooth_ns=20, mg_residual_calls=2, mg_residual_ns=10,
                mg_restriction_calls=2, mg_restriction_ns=10, mg_prolongation_calls=2, mg_prolongation_ns=15,
                mg_post_smooth_calls=2, mg_post_smooth_ns=15, mg_terminal_calls=2, mg_terminal_ns=20,
                mg_direct_mpi_calls=2, mg_direct_mpi_ns=4,
                mg_finest_pre_smooth_ns=20, mg_finest_post_smooth_ns=15,
                source_meta_sha256=source)
    empty = dict(initialized=1, complete=1, attempts=0, successes=0, failures=0,
                 apply_ns=0, reduction_ns=0, visits=0, pre_smooth_calls=0, pre_smooth_ns=0,
                 residual_calls=0, residual_ns=0, restriction_calls=0, restriction_ns=0,
                 prolongation_calls=0, prolongation_ns=0, post_smooth_calls=0, post_smooth_ns=0,
                 terminal_calls=0, terminal_ns=0, direct_mpi_calls=0, direct_mpi_ns=0,
                 halo_wait_ns=0, halo_control_ns=0, halo_control_calls=0,
                 source_meta_sha256=source)
    records = [dict(empty, level=-1, attempts=2, successes=2, apply_ns=100, reduction_ns=5),
               dict(empty, level=0, visits=2, pre_smooth_calls=2, pre_smooth_ns=20,
                    residual_calls=2, residual_ns=10, restriction_calls=2, restriction_ns=10,
                    prolongation_calls=2, prolongation_ns=15, post_smooth_calls=2, post_smooth_ns=15,
                    halo_wait_ns=3, halo_control_ns=2, halo_control_calls=2),
               dict(empty, level=1, visits=2, terminal_calls=2, terminal_ns=20,
                    direct_mpi_calls=2, direct_mpi_ns=4, halo_wait_ns=2,
                    halo_control_ns=3, halo_control_calls=3)]
    for rank in (0, 1):
        write(root / 'solver-rank-{}.csv'.format(rank),
              [dict(loop, step=step, rank=rank) for step in (7, 8)])
        write(root / 'mg-rank-{}.csv'.format(rank),
              [dict(row, step=step, rank=rank) for step in (7, 8) for row in records])
    write(root / 'performance.csv', [dict(step=s, rank=r, dropped_loops=0,
          pressure_calls=0, diagonal_calls=1, spatial_calls=0, A_apply_ns=30,
          M_apply_ns=120, source_meta_sha256=source) for s in (7, 8) for r in (0, 1)])
    (root / 'conservation.csv').write_text('step,rank\n7,0\n8,0\n')


def main():
    script = sys.argv[1]
    with tempfile.TemporaryDirectory(prefix='hundun-observe-mg-') as directory:
        root = Path(directory)
        fixture(root)
        def run(*flags):
            return subprocess.run([sys.executable, script, str(root)] + list(flags),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True, timeout=10)
        result = run()
        assert result.returncode == 0, result.stderr
        payload = json.loads(result.stdout)
        assert payload['schema'] == 'HUNDUN_LOOP_OBSERVATION_V5' and payload['complete']
        assert payload['validated_steps'] == 2 and payload['loop_count'] == 2
        assert payload['mg']['totals_by_rank'][0]['attempts'] == 4
        assert payload['mg']['totals_by_rank'][1]['apply_ns'] == 200
        levels = payload['mg']['levels']
        assert len(levels) == 4
        assert [(r['rank'], r['level'], r['pre_smooth_ns'], r['terminal_ns']) for r in levels] == [
            (0, 0, 40, 0), (0, 1, 0, 40), (1, 0, 40, 0), (1, 1, 0, 40)]
        details = root / 'details.jsonl'
        streamed = run('--details-output', str(details))
        assert streamed.returncode == 0, streamed.stderr
        streamed = json.loads(streamed.stdout)
        assert streamed['mg'] == payload['mg'] and streamed['loops'] == []
        assert [json.loads(line) for line in details.read_text().splitlines()] == payload['loops']
        baseline = {p.name: p.read_bytes() for p in root.iterdir() if p.name != details.name}
        checks = []

        def reset():
            for p in root.iterdir():
                p.unlink()
            for name, value in baseline.items():
                (root / name).write_bytes(value)

        def change(name, modify):
            with (root / name).open() as stream:
                rows = list(csv.DictReader(stream))
            modify(rows)
            write(root / name, rows)

        def reject(label, reason=None):
            before = {p.name: p.read_bytes() for p in root.iterdir()}
            result = run('--output', str(root / 'strict.json'), '--details-output', str(details))
            assert result.returncode == 2, (label, result.stdout, result.stderr)
            if reason:
                assert reason in result.stderr, (label, result.stderr)
            assert not (root / 'strict.json').exists() and not details.exists(), label
            assert {p.name: p.read_bytes() for p in root.iterdir()} == before, label
            checks.append(label)

        for name in ('mg-rank-1.csv', 'solver-rank-1.csv'):
            reset()
            (root / name).unlink()
            reject('missing ' + name, 'rank set')
        reset()
        for name in ('mg-rank-1.csv', 'solver-rank-1.csv'):
            (root / name).write_text(baseline[name].decode().splitlines()[0] + '\n')
        change('performance.csv', lambda rows: rows.__setitem__(slice(None), [r for r in rows if r['rank'] == '0']))
        reject('same rank absent from all ledgers')

        for key, value in (('complete', 0), ('initialized', 0), ('rank', 0),
                           ('level', 2), ('visits', -1), ('apply_ns', 'nan'),
                           ('attempts', 2**64 - 1), ('successes', 1),
                           ('pre_smooth_ns', 1), ('source_meta_sha256', 'f' * 64)):
            reset()
            change('mg-rank-1.csv', lambda rows: rows[0].update({key: value}))
            reject('MG invalid ' + key)
        for key, value in (('mg_enabled', 0), ('mg_complete', 0), ('mg_attempts', 3),
                           ('mg_apply_ns', 2**64 - 1), ('mg_apply_ns', 121),
                           ('mg_pre_smooth_ns', 101), ('mg_finest_pre_smooth_ns', 19),
                           ('mg_halo_wait_ns', 6)):
            reset()
            change('solver-rank-1.csv', lambda rows: rows[0].update({key: value}))
            reject('loop invalid ' + key)
        reset()
        change('solver-rank-1.csv', lambda rows: [r.pop('mg_terminal_calls') for r in rows])
        reject('missing MG loop column', 'MG loop columns')
        reset()
        change('mg-rank-1.csv', lambda rows: rows.pop(2))
        reject('missing finest/terminal level', 'level coverage')
        reset()
        change('mg-rank-1.csv', lambda rows: rows.insert(1, dict(rows[0])))
        reject('duplicate MG row', 'capacity')
        reset()
        change('mg-rank-1.csv', lambda rows: rows[2].update(terminal_ns=19))
        reject('level cost lost', 'loop/level attribution')
        reset()
        path = root / 'mg-rank-1.csv'
        path.write_bytes(path.read_bytes().rstrip(b'\r\n'))
        reject('truncated MG newline', 'truncated MG CSV')
        reset()
        (root / 'mg-rank-2.csv').write_bytes(baseline['mg-rank-1.csv'])
        reject('unexpected MG rank file', 'rank set')
        reset()
        change('mg-rank-1.csv', lambda rows: rows.extend(dict(r, step=9) for r in rows[-3:]))
        reject('MG rows beyond window', 'beyond frozen')
        reset()
        path = root / 'mg-layout.meta'
        path.write_text(path.read_text().replace('level 0 0 4', 'level 0 0 3'))
        reject('unbound layout change', 'layout source identity')

        # Correctly re-bind an invalid producer layout: the independent rank /
        # level contract must reject it, not only the old content hash.
        for original, replacement in (('rank 1 2\n', ''), ('rank 0 2', 'rank 0 33'),
                                       ('level 0 0 4', 'level 0 0 0')):
            reset()
            layout = baseline['mg-layout.meta'].decode().replace(original, replacement)
            (root / 'mg-layout.meta').write_text(layout)
            meta = baseline['RUN.meta'].decode().replace(
                hashlib.sha256(baseline['mg-layout.meta']).hexdigest(), hashlib.sha256(layout.encode()).hexdigest())
            (root / 'RUN.meta').write_text(meta)
            identity = hashlib.sha256(meta.encode()).hexdigest()
            for name in ('solver-rank-0.csv', 'solver-rank-1.csv', 'mg-rank-0.csv', 'mg-rank-1.csv', 'performance.csv'):
                change(name, lambda rows: [r.update(source_meta_sha256=identity) for r in rows])
            reject('invalid rebound layout', 'MG layout')
        reset()
        change('mg-rank-1.csv', lambda rows: rows.__setitem__(slice(None), rows[:3]))
        reject('failed tail step', 'MG step 8')
        partial = run('--allow-partial', '--details-output', str(details))
        assert partial.returncode == 0, partial.stderr
        partial = json.loads(partial.stdout)
        assert not partial['complete'] and partial['validated_steps'] == 1 and partial['loop_count'] == 1
        assert partial['mg']['totals_by_rank'][1]['attempts'] == 2
        assert not partial['details']['complete'] and len(details.read_text().splitlines()) == 1
        reset()
        (root / 'RUN.meta').write_text('damaged\n')
        partial = run('--allow-partial')
        assert partial.returncode == 0, partial.stderr
        assert not json.loads(partial.stdout)['complete']
        print('MG observer synthetic ledger, details, partial and {} rejection checks PASS'.format(len(checks)))


if __name__ == '__main__':
    main()
