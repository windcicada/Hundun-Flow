"""Physical ESF species/element budgets through native MPI and Restart CLI."""
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

product, case, launcher, rank_flag, destination, validator = sys.argv[1:]
root = Path(destination)
if root.exists():
    shutil.rmtree(str(root))
root.mkdir(parents=True)


def run(name, ranks, steps, restart=None):
    output = root / name
    command = [launcher, rank_flag, str(ranks), product, 'run', case,
               '--output', str(output), '--steps', str(steps),
               '--output-interval', '0', '--restart-interval', '1',
               '--diagnostics-interval', '1']
    if restart is not None:
        command += ['--restart', str(restart)]
    with (root / (name + '.log')).open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=40)
    validation = [sys.executable, validator, 'runtime', str(output / 'evidence.jsonl')]
    if restart is not None:
        generation = (restart / 'current').read_text().strip()
        validation += ['--run-start-manifest', str(restart / generation / 'manifest.bin')]
    subprocess.run(validation, check=True, timeout=10)
    rows = [json.loads(line) for line in (output / 'diagnostics.jsonl').read_text().splitlines()]
    assert len(rows) == steps
    for row in rows:
        budget = row['payload']['composition_balance']
        assert budget['scope'] == 'gas_transport_reaction'
        assert budget['density'] == 'field0' and budget['composition'] == 'physical_ensemble_mean'
        assert budget['revision'] > 0 and budget['duration_s'] == row['payload']['dt']
        assert not budget['after_parcel_exchange']
        assert {r['name'] for r in budget['species']} == {'A', 'B'}
        assert {r['name'] for r in budget['elements']} == {'N'}
        for group in ('species', 'elements'):
            for entry in budget[group]:
                assert all(math.isfinite(v) for k, v in entry.items() if k != 'name')
                residual = (entry['temporal_rate'] + entry['transport_outflow'] + entry['pressure_outflow']
                            - entry['noise_source'] - entry['mixing_source'] - entry['chemistry_source'])
                assert math.isclose(residual, entry['defect'], rel_tol=1e-10, abs_tol=1e-12)
                assert entry['relative_defect'] < 1e-6, (name, row['step'], entry)
                rate = (entry['current_inventory'] - entry['accepted_inventory']) / budget['duration_s']
                assert math.isclose(rate, entry['temporal_rate'], rel_tol=1e-9, abs_tol=1e-10)
        mass = sum(r['current_inventory'] for r in budget['species'])
        assert math.isclose(mass, row['payload']['mass_kg'], rel_tol=1e-12, abs_tol=1e-12)
        assert abs(sum(r['chemistry_source'] for r in budget['species'])) < 1e-12
        assert abs(budget['elements'][0]['chemistry_source']) < 1e-12
    for previous, current in zip(rows, rows[1:]):
        for group in ('species', 'elements'):
            for a, b in zip(previous['payload']['composition_balance'][group],
                            current['payload']['composition_balance'][group]):
                assert math.isclose(a['current_inventory'], b['accepted_inventory'], rel_tol=1e-12, abs_tol=1e-12)
    return rows


def compare(a, b):
    assert a['step'] == b['step'] and a['time'] == b['time']
    for group in ('species', 'elements'):
        for x, y in zip(a['payload']['composition_balance'][group], b['payload']['composition_balance'][group]):
            assert x['name'] == y['name']
            for key in x:
                if key != 'name':
                    # Use the same max(1, |a|, |b|) normalization as the
                    # native cross-partition field comparison.
                    scale = max(1., abs(x[key]), abs(y[key]))
                    assert abs(x[key] - y[key]) / scale < 1e-9, (group, key, x, y)


reference = run('r1', 1, 3)
two = run('r2', 2, 2)
four = run('r4', 4, 2)
resumed = run('resume', 4, 1, root / 'r2' / 'Restart')
for i in range(2):
    compare(reference[i], two[i])
    compare(reference[i], four[i])
compare(reference[2], resumed[0])
for group in ('species', 'elements'):
    for a, b in zip(two[-1]['payload']['composition_balance'][group],
                    resumed[0]['payload']['composition_balance'][group]):
        assert math.isclose(a['current_inventory'], b['accepted_inventory'], rel_tol=1e-12, abs_tol=1e-12)
print('composition CLI: 1/2/4 ranks, 2-to-4 Restart, inventories, species/element terms and gates passed')
