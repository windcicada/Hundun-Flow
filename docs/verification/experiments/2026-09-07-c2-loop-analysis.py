#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Read-only analysis of the independently revalidated 5501-5510 window."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

p = argparse.ArgumentParser()
p.add_argument('--run-root', type=Path, required=True)
p.add_argument('--observation', type=Path, required=True)
p.add_argument('--details', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


observation = json.loads(a.observation.read_text())
assert observation['complete'] and not observation['issues']
assert observation['expected_step_range'] == [5501, 5510]
assert observation['expected_rank_count'] == 128 and observation['loop_count'] == 61
assert observation['details']['complete']
assert digest(a.details) == observation['details']['sha256']
for source in ('solver-rank-0.csv', 'conservation.csv', 'RUN.meta'):
    assert digest(a.run_root / source) == observation['files_sha256'][source], source
with (a.run_root / 'solver-rank-0.csv').open() as stream:
    rows = {(int(r['step']), int(r['corrector']), int(r['refinement'])): r for r in csv.DictReader(stream)}
assert len(rows) == 61 and all(int(r['attempt']) == int(r['scalar_coupling_sweep']) == 1 for r in rows.values())


def merit(row):
    return math.hypot(float(row['selected_continuity']), float(row['selected_energy']))


def forcing(key):
    # Source-rule reconstruction, NOT a runtime tolerance column. This frozen
    # zero-scalar case has base=terminal=1e-6, guarded ceiling=1e-4.
    step, corrector, refinement = key
    if refinement == 0:
        return 1e-4, 'mandatory_inexact'
    current = rows[(step, 2, refinement - 1)]
    previous = rows[(step, 1, 0)] if refinement == 1 else rows[(step, 2, refinement - 2)]
    if max(float(current['selected_continuity']), float(current['selected_energy'])) <= 1e-5:
        return 1e-6, 'terminal_band'
    if merit(current) >= 0.9 * merit(previous):
        return 1e-6, 'nondecrease_vs_C1' if refinement == 1 else 'nondecrease_same_C2'
    return min(1e-4, max(1e-6, 0.1 * merit(current))), 'residual_aware'


groups = defaultdict(list)
with a.details.open() as stream:
    for line in stream:
        r = json.loads(line)
        key = (r['step'], r['corrector'], r['refinement'])
        assert key in rows and r['ranks'] == 128
        assert all(len(v) == 1 for v in r['counts_by_rank'].values())
        groups[(r['corrector'], r['refinement'], r['kind'])].append((key, r))
table = []
for group, entries in sorted(groups.items()):
    timings = {k: sum(v['rank_mean_ns'][k] for _, v in entries) / 1e9
               for k in entries[0][1]['rank_mean_ns']}
    counts = {k: sum(v['counts_by_rank'][k][0] for _, v in entries)
              for k in entries[0][1]['counts_by_rank']}
    row = {'corrector': group[0], 'refinement': group[1], 'kind': group[2],
           'loops': len(entries), 'counts': counts, 'rank_mean_window_s': timings,
           'M_ms_per_apply': timings['M_ns'] * 1e3 / counts['M_calls'],
           'A_ms_per_apply': timings['A_ns'] * 1e3 / counts['A_calls'],
           'linear_contraction_median_across_loops': statistics.median(v['linear_contraction'][0] for _, v in entries),
           'linear_initial_range': [min(float(rows[k]['linear_initial']) for k, _ in entries),
                                    max(float(rows[k]['linear_initial']) for k, _ in entries)],
           'linear_final_range': [min(float(rows[k]['linear_final']) for k, _ in entries),
                                  max(float(rows[k]['linear_final']) for k, _ in entries)],
           'reconstructed_forcing_reasons': dict(Counter(forcing(k)[1] for k, _ in entries)),
           'reconstructed_relative_tolerances': sorted(set(forcing(k)[0] for k, _ in entries))}
    table.append(row)
    print('C{} r{} {} loops={} iterations={} M/solve={:.6f}/{:.6f}s M_ms/apply={:.6f} forcing={}'.format(
        *group, len(entries), counts['iterations'], timings['M_ns'], timings['solve_ns'],
        row['M_ms_per_apply'], row['reconstructed_forcing_reasons']))
with (a.run_root / 'conservation.csv').open() as stream:
    balances = list(csv.DictReader(stream))
assert [int(r['step']) for r in balances] == list(range(5501, 5511))
dt = {float(r['dt']) for r in rows.values()}
assert len(dt) == 1 and all(int(r['bdf_order']) == 2 for r in balances)
momentum = {key: {'minimum': min(float(r[key]) for r in balances),
                  'maximum': max(float(r[key]) for r in balances),
                  'last': float(balances[-1][key])}
            for key in ('momentum_x_linf_N', 'momentum_y_linf_N', 'momentum_z_linf_N',
                        'U_rms_m_s', 'momentum_x_normalized', 'momentum_y_normalized', 'momentum_z_normalized')}
result = {'observation': str(a.observation), 'observation_sha256': digest(a.observation),
          'details': str(a.details), 'details_sha256': observation['details']['sha256'],
          'script_sha256': digest(Path(__file__)),
          'complete_window': True, 'formal_speed_comparison': False, 'groups': table,
          'C2_r0_selected_merit_over_C1': [merit(rows[(s, 2, 0)]) / merit(rows[(s, 1, 0)]) for s in range(5501, 5511)],
          'dt': next(iter(dt)), 'a0': 1.5 / next(iter(dt)), 'momentum': momentum,
          'limits': ['No new solver run or kernel optimization.',
                     'Forcing is reconstructed from current unchanged source rules, not a measured CSV field.',
                     'M includes its communication. No level-specific smooth/coarse/wait attribution exists here.',
                     'Only one dt; raw momentum term-scale and equal-physical-time refinement comparisons remain open.']}
with a.output.open('x') as stream:
    json.dump(result, stream, indent=2, allow_nan=False)
    stream.write('\n')
