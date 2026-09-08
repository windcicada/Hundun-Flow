#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
import csv
import hashlib
import json
import math
from pathlib import Path

audit = Path(__file__).resolve().parent
prefix = audit / 'startup-prefix'
meta = dict(line.split(' ', 1) for line in (prefix / 'RUN.meta').read_text().splitlines() if ' ' in line)
for key, expected in dict(starting_step='10000', requested_steps='25000', expected_ranks='128',
                          starting_sample_steps='0', statistics_epoch_start_step='7000',
                          statistics_sampling_start_step='17001', restart_requires_recovery='0',
                          restart_method_recovery='0', observation_schema='6', observe_mg_cost='1',
                          observe_fgmres_recovery='1', trace_cell_count='0').items():
    assert meta[key] == expected, key
assert meta['executable_sha256'] == '345cca4802286a7ae1b2bd39b7c25afceeedb8a13fc60b5fdf71f80d92f0f67a'
with (prefix / 'health.csv').open(newline='') as stream:
    rows = [r for r in csv.DictReader(stream) if None not in r and all(v is not None for v in r.values())]
assert len(rows) >= 3
assert [int(r['step']) for r in rows] == list(range(10001, 10001 + len(rows)))
maxima = dict.fromkeys(('eos', 'continuity', 'energy', 'closed_mass', 'gauge'), 0.0)
solid = None
for row in rows:
    assert row['bdf_order'] == '2' and row['restart_recovery'] == '0'
    assert row['attempts'] == '1' and row['retry'] == '0'
    assert row['terminal_audit_present'] == '1'
    for key in maxima:
        value = float(row[key])
        assert math.isfinite(value) and 0 <= value <= float(row[key + '_tolerance'])
        maxima[key] = max(maxima[key], value)
    current = {k:v for k,v in row.items() if k.startswith('solid_placeholder_')}
    if solid is None:
        solid = current
    assert current == solid
with (prefix / 'evidence.jsonl').open() as stream:
    evidence = [json.loads(line) for line in stream]
anchor = evidence[0]['run_start']
assert anchor['previous_step'] == 10000
assert anchor['history'] == dict(source_format_version=3, source_signature=12213963202598979269,
                                  target_signature=12213963202598979269, policy='require_compatible')
assert anchor['restart_manifest_sha256'] == 'b5bafe56b0c600016b7ebadda82edd896fe83fe6a3c6f17cef26b895f57310a2'
observer = json.loads((audit / 'startup-observation.json').read_text())
assert observer['complete'] is False and observer['validated_steps'] >= len(rows)
assert observer['expected_rank_count'] == 128 and observer['loop_count'] > 0
assert observer['issues'] == ['missing/truncated expected step {}'.format(10001 + observer['validated_steps'])]
artifacts = {}
for name in ('HANDOFF_VALIDATED.json', 'SOURCE-10000.sha256', 'FROZEN.sha256', 'bridge.log',
             'launch-long.sh', 'GITHUB_MAIN_VERIFIED.txt', 'startup-observation.json',
             'startup-loops.jsonl', 'startup-runtime-validation.log', 'source-unchanged-after-start.log'):
    artifacts[name] = hashlib.sha256((audit / name).read_bytes()).hexdigest()
result = dict(schema='HUNDUN_RE3900_LONG_STARTUP_V1', startup_verified=True, long_run_complete=False,
              run_root=str(audit.parent / 'long-observed-35000-20260908'),
              service='hundun-re3900-observed-long-20260908.service',
              health_step_range=[10001, int(rows[-1]['step'])], accepted_health_rows=len(rows),
              terminal_maxima=maxima, all_bdf2=True, retries=0, restart_anchor=anchor,
              observer_validated_steps=observer['validated_steps'], observer_loops=observer['loop_count'],
              observer_complete=False, observer_issues=observer['issues'],
              evidence_records=len(evidence), runtime_validator_exit=0,
              epoch_start=7000, sampling_start=17001, sample_steps=0,
              artifacts_sha256=artifacts,
              limitations=['Live output copied without stopping the run; stream prefixes can end on different steps.',
                           'Startup/partial observation is not 35000-step stability or COAST replacement acceptance.',
                           'Next new-run checkpoint and Visit are due at step 10500; source 10000 remains durable.',
                           'Statistics reset_reason is inherited epoch provenance, not a new reset.'])
with (audit / 'STARTUP_VERIFIED.json').open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(result, indent=2))
