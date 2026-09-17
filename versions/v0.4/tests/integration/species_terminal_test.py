# SPDX-License-Identifier: Apache-2.0
"""The default local-time species audit drives refinement across outer solves."""
import json
from timing_check import check_timing
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

program, launcher, validator = sys.argv[1:]
fixture = Path(__file__).resolve().parents[1] / 'fixtures' / 'esf-cn'
with tempfile.TemporaryDirectory(prefix='species-') as directory:
    root = Path(directory)
    case = root / 'case'
    shutil.copytree(fixture, case)
    model = json.loads((case / 'case.json').read_text())
    reaction = model['reaction']
    reaction['model'] = 'finite_rate_mean'
    reaction.pop('ensemble')
    reaction.pop('mixing')
    reaction['chemistry_solver'].update(relative_tolerance=1e-10, absolute_tolerance=1e-14)
    (case / 'case.json').write_text(json.dumps(model))
    for ranks in (1, 2, 4):
        output = root / ('run' + str(ranks))
        args = [launcher, '-n', str(ranks), program, 'run', str(case), '--output', str(output),
                '--steps', '2', '--initial-state', '101325,300,0,0,0,0.25',
                '--output-interval', '0', '--restart-interval', '1', '--diagnostics-interval', '0']
        result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=60)
        if result.returncode:
            raise RuntimeError(' '.join(args) + '\n' + result.stdout[-7000:])
        check_timing(output/'monitor.jsonl', {'reaction_sources', 'mean_reaction'})
        records = [json.loads(line) for line in (output / 'evidence.jsonl').read_text().splitlines()]
        assert [row['step'] for row in records] == [1, 2]
        for row in records:
            assert row['cold']['normalization'] == 'local_time'
            assert row['cold']['species_residual'] < 128 * sys.float_info.epsilon
        subprocess.run([sys.executable, validator, 'runtime', str(output / 'evidence.jsonl')], check=True)
        print('ranks={} second_species_residual={:.17g} outer={}'.format(
            ranks, records[-1]['cold']['species_residual'], records[-1]['cold']['outer_iterations']))
