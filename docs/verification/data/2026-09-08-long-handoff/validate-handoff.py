#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Verify the completed checkpoint bundle and accepted bridge prefix, read-only."""
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

AUDIT = Path(__file__).resolve().parent
SOURCE = AUDIT / 'source-10000'
OLD = AUDIT.parent / 'long-thermodynamic-35000-20260907'
module_spec = importlib.util.spec_from_file_location('evidence_validator', str(AUDIT / 'frozen/v04_evidence_validate.py'))
validator = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(validator)
completion = (SOURCE / 'step-00000000000000010000.complete').read_text().splitlines()
assert completion[0] == 'HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1' and completion[-1] == 'end'
done = dict(line.split(' ', 1) for line in completion[1:-1])
assert done['step'] == '10000'
generation = (SOURCE / 'Restart/current').read_text().strip()
assert generation == done['restart_generation'] and generation.startswith('generation-10000-')
manifest = validator.load_v04_restart_manifest(SOURCE / 'Restart' / generation / 'manifest.bin')
assert manifest['step'] == 10000 and manifest['source_format_version'] == 3
assert manifest['method_history_signature'] == 12213963202598979269
assert not manifest['backward_euler_recovery']
files = sorted((SOURCE / 'Restart' / generation).glob('rank-*.bin'))
assert [p.name for p in files] == ['rank-{:08d}.bin'.format(r) for r in range(128)]
assert all(p.stat().st_size > 0 for p in files)
stats = json.loads((SOURCE / done['statistics']).read_text())
assert stats['snapshot_step'] == 10000 and stats['sample_steps'] == 0
assert stats['statistics_epoch'] == dict(start_step=7000, development_steps=10000,
                                        sampling_start_step=17001, reset_reason='method_recovery')
assert (SOURCE / done['accumulator']).stat().st_size > 0
rows = []
with (OLD / 'health.csv').open(newline='') as stream:
    for row in csv.DictReader(stream):
        if None not in row and all(v is not None for v in row.values()):
            if 9951 <= int(row['step']) <= 10000:
                rows.append(row)
assert [int(r['step']) for r in rows] == list(range(9951, 10001))
maxima = dict.fromkeys(('eos', 'continuity', 'energy', 'closed_mass', 'gauge'), 0.0)
solid = None
for row in rows:
    assert int(row['attempts']) == 1 and int(row['retry']) == 0
    assert int(row['bdf_order']) == 2 and int(row['restart_recovery']) == 0
    assert int(row['terminal_audit_present']) == 1
    assert int(row['fluid_count']) == 5994352 and int(row['solid_placeholder_count']) == 75920
    for key in maxima:
        value = float(row[key])
        assert math.isfinite(value) and 0 <= value <= float(row[key + '_tolerance'])
        maxima[key] = max(maxima[key], value)
    current_solid = {k: v for k, v in row.items() if k.startswith('solid_placeholder_')}
    if solid is None:
        solid = current_solid
    assert current_solid == solid, 'solid regional extrema drift'
    for region in ('fluid', 'solid_placeholder'):
        for field in ('T', 'p', 'rho'):
            low, high = [float(row[region + '_' + field + '_' + end]) for end in ('min', 'max')]
            assert 0 < low <= high and math.isfinite(high)
visit = sorted((OLD / 'Visit').glob('step-00000000000000010000-rank-*.vtr'))
assert len(visit) == 128 and all(p.stat().st_size > 0 for p in visit)
assert (OLD / 'Visit/step-00000000000000010000.visit').is_file()
result = dict(schema='HUNDUN_RE3900_CHECKPOINT_HANDOFF_V1', checkpoint_bundle_complete=True,
              original_run_complete=False, bridge_steps=[9951, 10000], bridge_steps_accepted=50,
              source=str(SOURCE), generation=generation, manifest=manifest,
              terminal_maxima=maxima, bdf2_steps=50, retries=0,
              statistics_epoch=stats['statistics_epoch'], samples=0, rank_files=128,
              visit_files=128, solid_regional_extrema_unchanged=True,
              runner_sha256=hashlib.sha256((AUDIT / 'frozen/v04_thin_domain_runner').read_bytes()).hexdigest(),
              production_source='7c03a54a9f1061d1de0e0922ed38a3bcc8f2412f',
              repository_head='d706a41f37801553de29bc52a5ce85209a10aeee',
              notes=['Step 9951 timing includes SIGSTOP interval; exclude it from throughput.',
                     'SHA-256 copy check and manifest integrity verified separately; native Reader checks all rank contents on restart.',
                     'Stopped old output prefix is not a successful final stream-close or COMPLETED claim.'])
assert result['runner_sha256'] == '345cca4802286a7ae1b2bd39b7c25afceeedb8a13fc60b5fdf71f80d92f0f67a'
with (AUDIT / 'HANDOFF_VALIDATED.json').open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(result, indent=2))
