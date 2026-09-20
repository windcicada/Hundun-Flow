#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""GTMC JL4/4-field dynamic TCR through native run and partitioned Restart."""
from flow_budget import check_flow, check_omission_counterexamples
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

binary, fixture, mpi, work, validator, compare = map(Path, sys.argv[1:])
if work.exists():
    shutil.rmtree(work)
shutil.copytree(fixture, work)
dt = json.loads((work/'case.json').read_text())['time']['initial_dt']

def call(args, log):
    with log.open('w') as stream:
        result = subprocess.run(list(map(str, args)), stdout=stream, stderr=stream)
    assert result.returncode == 0, log.read_text()[-16000:]

call([binary, 'check', work], work/'check.log')
text = (work/'check.log').read_text()
assert 'sgs=vreman' in text and 'tcr_model=cdphyso_dynamic_v1 tcr_mode=experimental' in text

def run(label, ranks, steps, restart=None):
    output = work.with_name(work.name+label)
    if output.exists():
        shutil.rmtree(output)
    args = [mpi, '--oversubscribe', '--bind-to', 'none', '-n', ranks, binary,
            'run', work, '--output', output, '--steps', steps, '--max-dt', dt,
            '--output-interval', 0, '--restart-interval', 1, '--diagnostics-interval', 1]
    if restart:
        args += ['--restart', restart/'Restart']
    else:
        args += ['--initial-state', '100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001']
    call(args, work/(label+'.log'))
    checks = [sys.executable, validator, 'runtime', output/'evidence.jsonl']
    if restart:
        generation = (restart/'Restart/current').read_text().strip()
        checks += ['--run-start-manifest', restart/'Restart'/generation/'manifest.bin']
    call(checks, work/(label+'-check.log'))
    rows = [json.loads(line) for line in (output/'diagnostics.jsonl').read_text().splitlines()]
    for row in rows:
        p = row['payload']
        check_flow(p)
        check_omission_counterexamples(p)
    return output

seed = run('1', 1, 3)
two = run('2', 2, 2, seed)
one = run('r', 1, 2, seed)
call([mpi, '--oversubscribe', '--bind-to', 'none', '-n', 2, compare,
      work, one/'Restart', two/'Restart'], work/'compare.log')
(work/'result.json').write_text(json.dumps(dict(
    scope='JL4 4-field ESF, Vreman, dynamic TCR; native steps 1-5 and 1/2-rank Restart',
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    comparison=(work/'compare.log').read_text(), passed=True), indent=2)+'\n')
print('JL4 dynamic TCR: native chemistry, history and 1/2-rank Restart passed')
