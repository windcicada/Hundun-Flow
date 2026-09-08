#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Public V6 observer CLI: synthetic work ledger, never CFD timing evidence."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from solver_observe_mg_cli_test import fixture as mg_fixture, write


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def bind(root, text):
    (root / 'RUN.meta').write_text(text)
    digest = hashlib.sha256(text.encode()).hexdigest()
    for path in list(root.glob('solver-rank-*.csv')) + list(root.glob('mg-rank-*.csv')) + [root / 'performance.csv']:
        rows = read(path)
        for row in rows:
            row['source_meta_sha256'] = digest
        write(path, rows)


def fixture(root, mg=True):
    mg_fixture(root)
    text = (root / 'RUN.meta').read_text().replace('observation_schema 5', 'observation_schema 6')
    text = text.replace('end\n', 'observe_fgmres_recovery 1\npressure_linear_algorithm fgmres\nend\n')
    if not mg:
        text = '\n'.join(line for line in text.splitlines() if not line.startswith('mg_layout_sha256')) + '\n'
        text = text.replace('observe_mg_cost 1', 'observe_mg_cost 0')
        (root / 'mg-layout.meta').unlink()
        for path in root.glob('mg-rank-*.csv'):
            path.unlink()
    # Worked counts: 3 Arnoldi A + initial A + discarded-column residual A +
    # interior candidate residual A = 6, with one unsafe restart. Synthetic ns
    # stay those of the independent MG fixture and are not performance samples.
    counters = dict(available=1, norm_restarts=1, unsafe_norms=1, discarded_columns=1,
        explicit_reorthogonalizations=0, unsafe_restarts=1, happy_restarts=0, length_restarts=0,
        initial_residual_applies=1, arnoldi_applies=3, unsafe_residual_applies=1,
        interior_residual_applies=1, cycle_residual_applies=0)
    for path in root.glob('solver-rank-*.csv'):
        rows = read(path)
        for row in rows:
            row.update(iterations=3, A_calls=6, M_calls=3)
            row.update({'fgmres_' + k: v for k, v in counters.items()})
            if mg:
                row.update(mg_attempts=3, mg_successes=3)
            else:
                for key in list(row):
                    if key.startswith('mg_'):
                        del row[key]
        write(path, rows)
    for path in root.glob('mg-rank-*.csv'):
        rows = read(path)
        for row in rows:
            if row['level'] == '-1':
                row.update(attempts=3, successes=3)
        write(path, rows)
    bind(root, text)


def main():
    script = sys.argv[1]
    checks = 0
    with tempfile.TemporaryDirectory(prefix='hundun-observe-recovery-') as directory:
        root = Path(directory)
        fixture(root)
        baseline = {p.name: p.read_bytes() for p in root.iterdir()}
        def restore():
            for name, data in baseline.items():
                (root / name).write_bytes(data)
        def run(*flags):
            return subprocess.run([sys.executable, script, str(root)] + list(flags),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=10)
        def good():
            result = run()
            assert result.returncode == 0, result.stderr
            payload = json.loads(result.stdout)
            assert payload['complete'] and payload['schema'] == 'HUNDUN_LOOP_OBSERVATION_V6'
            return payload
        payload = good()
        assert payload['validated_steps'] == 2 and len(payload['loops']) == 2
        for row in payload['loops']:
            profile = row['fgmres_recovery']
            assert profile['available']
            assert profile['counts_by_rank']['unsafe_restarts'] == [1]
            assert profile['counts_by_rank']['happy_restarts'] == [0]
            assert profile['counts_by_rank']['arnoldi_applies'] == [3]
        details = root / 'details.jsonl'
        result = run('--details-output', str(details))
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)['loops'] == []
        assert [json.loads(line) for line in details.read_text().splitlines()] == payload['loops']
        assert all((root / name).read_bytes() == data for name, data in baseline.items())

        def bad(name, mutate):
            nonlocal checks
            restore()
            mutate()
            result = run()
            assert result.returncode != 0, name + ': invalid evidence accepted'
            checks += 1
        path = root / 'solver-rank-1.csv'
        def change(key, value):
            rows = read(path)
            rows[0][key] = value
            write(path, rows)
        def omit(key):
            rows = read(path)
            for row in rows:
                del row[key]
            write(path, rows)
        def duplicate_header():
            text = path.read_text().replace('fgmres_happy_restarts', 'fgmres_unsafe_restarts', 1)
            path.write_text(text)
        for key, value in (('fgmres_available', 2), ('fgmres_arnoldi_applies', -1),
                ('fgmres_arnoldi_applies', 2**64 - 1), ('fgmres_arnoldi_applies', 2**64),
                ('fgmres_arnoldi_applies', 'nan'), ('fgmres_arnoldi_applies', '1.5'),
                ('fgmres_initial_residual_applies', 2), ('fgmres_norm_restarts', 2),
                ('fgmres_available', 0), ('fgmres_discarded_columns', 0),
                ('fgmres_explicit_reorthogonalizations', 3), ('M_calls', 2),
                ('iterations', 4), ('invoked', 0), ('source_meta_sha256', '0' * 64),
                ('dropped_loops', 1)):
            bad(key + '=' + str(value), lambda k=key, v=value: change(k, v))
        bad('missing column', lambda: omit('fgmres_unsafe_restarts'))
        bad('duplicate column', duplicate_header)
        bad('duplicate loop', lambda: write(path, [read(path)[0]] + read(path)))
        bad('missing loop', lambda: write(path, read(path)[1:]))
        # DictWriter emits CRLF; removing only LF leaves a valid CR newline
        # under Python's universal-newline reader, not a truncated record.
        bad('truncated tail', lambda: path.write_bytes(path.read_bytes().rstrip(b'\r\n')))
        def missing_rank():
            path.unlink()
            write(root / 'performance.csv', [r for r in read(root / 'performance.csv') if r['rank'] == '0'])
            (root / 'mg-rank-1.csv').unlink()
        bad('whole rank missing in both ledgers', missing_rank)
        for old, new in (('observe_fgmres_recovery 1', 'observe_fgmres_recovery 0'),
                ('observation_schema 6', 'observation_schema 5'),
                ('pressure_linear_algorithm fgmres', 'pressure_linear_algorithm unknown'),
                ('pressure_linear_algorithm fgmres', 'pressure_linear_algorithm bicgstab')):
            bad('metadata ' + new, lambda a=old, b=new: bind(root, (root / 'RUN.meta').read_text().replace(a, b)))
        def rank_disagreement():
            rows = read(path)
            rows[0].update(fgmres_unsafe_norms=0, fgmres_discarded_columns=0,
                fgmres_explicit_reorthogonalizations=1, fgmres_unsafe_restarts=0,
                fgmres_happy_restarts=1, fgmres_unsafe_residual_applies=0,
                fgmres_cycle_residual_applies=1)
            write(path, rows)
        bad('locally consistent but different recovery reasons across ranks', rank_disagreement)
        restore()
        write(path, read(path)[:1])
        partial = run('--allow-partial')
        assert partial.returncode == 0, partial.stderr
        partial = json.loads(partial.stdout)
        assert not partial['complete'] and partial['validated_steps'] == 1 and partial['loop_count'] == 1
        assert partial['loops'][0]['fgmres_recovery']['counts_by_rank']['unsafe_restarts'] == [1]
        restore()

        # An explicit non-FGMRES run may retain ordinary A/M work, but must
        # report unavailable recovery rather than fabricating zero samples.
        for name in ('solver-rank-0.csv', 'solver-rank-1.csv'):
            rows = read(root / name)
            for row in rows:
                for key in row:
                    if key.startswith('fgmres_'):
                        row[key] = 0
            write(root / name, rows)
        bind(root, (root / 'RUN.meta').read_text().replace('pressure_linear_algorithm fgmres',
                                                       'pressure_linear_algorithm bicgstab'))
        unsupported = good()
        assert all(r['fgmres_recovery'] == {'available': False, 'counts_by_rank': None}
                   for r in unsupported['loops'])
        bind(root, (root / 'RUN.meta').read_text().replace('pressure_linear_algorithm bicgstab',
                                                       'pressure_linear_algorithm fgmres'))
        assert run().returncode != 0, 'unobserved FGMRES work must not be complete'
        checks += 1

    with tempfile.TemporaryDirectory(prefix='hundun-recovery-without-mg-') as directory:
        root = Path(directory)
        fixture(root, mg=False)
        result = subprocess.run([sys.executable, script, str(root)], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True, timeout=10)
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)['complete'] and 'mg' not in json.loads(result.stdout)
    print('V6 FGMRES observer: {} rejection checks, independent work sums, optional MG, '
          'unsupported algorithm, partial whole-step output, streaming and source immutability PASS'.format(checks))


if __name__ == '__main__':
    main()
