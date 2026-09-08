#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Compare the two already completed, single-round windows; never launch CFD."""
import hashlib
import json
from pathlib import Path

audit = Path(__file__).resolve().parent
baseline = audit.parent / 'mg-observation-20260908'
output = audit / 'PERFORMANCE_DECISION.json'
assert not output.exists()

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()

def read(root):
    receipt = json.loads((root / 'PILOT_ACCEPTED.json').read_text())
    path = root / 'cost-summary.json'
    assert receipt['artifacts_sha256'][str(path.relative_to(root.parent))] == sha(path)
    cost = json.loads(path.read_text())
    assert cost['candidate_identity']['head'] == receipt['source_commit']
    window = cost['windows']['subsequent_steps']
    assert window['step_range'] == [9502, 9510] and window['logical_loops'] == 63
    return receipt, window

old_receipt, old = read(baseline)
new_receipt, new = read(audit)
assert old_receipt['source_commit'] == '461d7631aa251969ddfc4e639bcadb74c640070c'
assert new_receipt['source_commit'] == 'ed9b09ab2a019a7957030844e5437421351fd044'
assert new_receipt['payload_equivalence']['solver_nontiming_columns']['all_equal']
assert old['counts'] == new['counts']
metrics = {}
for name, receipt, window in [('baseline', old_receipt, old), ('candidate', new_receipt, new)]:
    p = window['rank_step_mean_seconds']
    q = window['rank_step_mean_loop_seconds']
    metrics[name] = {
        'process_wall_seconds': receipt['process_wall_seconds'],
        'mean_max_rank_advance_seconds': window['max_rank_advance_seconds']['mean'],
        'rank_mean_krylov_seconds_per_step': p['krylov'],
        'rank_mean_candidate_seconds_per_step': p['candidate'],
        'rank_mean_mg_seconds_per_step': q['mg_apply'],
        'rank_mean_mg_finest_pre_seconds_per_step': q['mg_finest_pre_smooth'],
        'mg': window['mg'],
        'counts': window['counts'],
    }
# The agreed rule rejects a candidate without a total-cost gain; do not
# discard a slow loop or rerun the same configuration to manufacture one.
no_total_gain = metrics['candidate']['process_wall_seconds'] >= metrics['baseline']['process_wall_seconds']
assert no_total_gain
result = {
    'schema': 'HUNDUN_MG_SINGLE_STAGE_PERFORMANCE_DECISION_V1',
    'decision': 'withdraw_algorithm_experiment',
    'reason': 'No total wall-time gain in the only permitted candidate round.',
    'scientific_and_payload_checks_passed': True,
    'rounds_per_configuration': 1,
    'metrics': metrics,
    'limitations': ['Not a statistical estimate of long-term slowdown or speedup.',
                    'Finest pre decreased, but inclusive MG/advance/process costs did not.',
                    'Communication is nested; increased wait has no established causal attribution.',
                    'Original frozen long run remains stopped; no COAST comparison or two-phase integration.'],
    'source_sha256': {str(p): sha(p) for p in [Path(__file__),
        baseline / 'PILOT_ACCEPTED.json', baseline / 'cost-summary.json',
        audit / 'PILOT_ACCEPTED.json', audit / 'cost-summary.json']},
}
with output.open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print('WITHDRAW: process %.2f -> %.2f s; advance %.6f -> %.6f s/step; numerical payload unchanged.' % (
    metrics['baseline']['process_wall_seconds'], metrics['candidate']['process_wall_seconds'],
    metrics['baseline']['mean_max_rank_advance_seconds'], metrics['candidate']['mean_max_rank_advance_seconds']))
