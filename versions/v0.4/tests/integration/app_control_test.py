#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exercise accepted-step monitor, one-shot output and save/stop across ranks."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary, mpi = map(lambda p: str(Path(p).resolve()), sys.argv[1:3])

def command(n, *args):
    result = subprocess.run([mpi, '--oversubscribe', '--bind-to', 'none', '-n', str(n), binary,
                             *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=45)
    assert result.returncode == 0, result.stdout[-6000:]
    return result.stdout

with tempfile.TemporaryDirectory(prefix='hf-control-') as tmp:
    root = Path(tmp)
    case, run = root/'c', root/'r'
    command(1, 'init-case', '--output', case)
    run.mkdir()
    (run/'stop').touch()
    (run/'output').touch()
    text = command(2, 'run', case, '--output', run, '--steps', 5,
                   '--max-dt', '1e-5', '--output-interval', 0, '--restart-interval', 0)
    assert 'COMPLETED steps=1 ' in text and 'step=1 time=' in text, text
    assert json.loads(command(1, 'status', run))['phase'] == 'stopped'
    receipts = [json.loads(s) for s in (run/'control.jsonl').read_text().splitlines()]
    assert [(r['request'], r['step']) for r in receipts] == [('stop', 1), ('output', 1)]
    assert (run/'Restart/current').exists()
    assert list((run/'Visit').glob('*.visit'))
    assert not any((run/p).exists() for p in ('stop','output','stop.pending','output.pending'))
    rows = [json.loads(s) for s in (run/'monitor.jsonl').read_text().splitlines()]
    assert len(rows) == 1 and rows[0]['payload']['seconds'] > 0
    measured = rows[0]['payload']
    phases = measured['cn_phases']
    assert phases['scope'] == 'max_rank_all_attempts_inclusive_mpi'
    assert len(phases['seconds']) == 10
    assert all(0 <= value <= measured['seconds'] for value in phases['seconds'].values())
    if measured['outer_iterations']:
        assert phases['seconds']['momentum_solve'] > 0
        assert phases['seconds']['pressure_solve'] > 0
    resumed = root/'s'
    command(4, 'run', case, '--output', resumed, '--steps', 3, '--max-dt','1e-5',
            '--restart', run/'Restart', '--output-interval',0, '--restart-interval',0,
            '--monitor-interval',2)
    assert not (resumed/'Visit').exists()
    rows = [json.loads(s) for s in (resumed/'monitor.jsonl').read_text().splitlines()]
    assert [r['step'] for r in rows] == [2,4]
    assert json.loads(command(1,'status',resumed))['phase'] == 'completed'
    model = json.loads((case/'case.json').read_text())
    model['time']['scheme'] = 'cn_be'
    (case/'case.json').write_text(json.dumps(model))
    cn = root/'cn'
    command(2, 'run', case, '--output', cn, '--steps', 1, '--max-dt', '1e-5',
            '--output-interval', 0, '--restart-interval', 0)
    measured = json.loads((cn/'monitor.jsonl').read_text())['payload']
    phases = measured['cn_phases']['seconds']
    assert measured['outer_iterations'] > 0
    assert phases['momentum_solve'] > 0 and phases['pressure_solve'] > 0
    assert all(0 <= value <= measured['seconds'] for value in phases.values())
    print('monitor independent; output receipt; save/stop; 2-to-4 rank restart; CN timing: PASS')
