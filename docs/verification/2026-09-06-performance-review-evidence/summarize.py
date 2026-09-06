#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Read retained observations; do not launch solvers or rewrite old evidence."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
BENCH = Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903')
RUN = BENCH / 'trial-D0p02-zpi2-52/long-review-749to35000-20260905'
sys.path.insert(0, str(BENCH))
from analyze_strict_results import validate_coast


def mean(items):
    return statistics.mean(items)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--observation', type=Path, default=HERE / 'hundun-long-observation.json')
    parser.add_argument('--output', type=Path, default=HERE / 'module-summary.json')
    arguments = parser.parse_args()
    observation = json.loads(arguments.observation.read_text())
    steps = observation['steps']
    with (RUN / 'evidence.jsonl').open() as stream:
        evidence = [json.loads(line) for line in stream]
    assert [s['step'] for s in steps] == [e['step'] for e in evidence]
    assert [e['step'] for e in evidence] == list(range(750, 1297))
    assert all(s['ranks'] == 128 for s in steps)
    for e in evidence:
        assert e['candidate_identity']['executable_sha256'] == observation['binary_sha256']
        assert e['bdf_order'] == 2 and not e['retry'] and not e['restart_recovery']
        assert not e['temporal_method_fallback']
        audit = e['terminal_physical_audit']
        assert audit['present']
        for name in ('eos', 'continuity', 'energy', 'closed_mass', 'gauge'):
            assert audit[name + '_residual'] <= audit[name + '_tolerance']
    h = {'run': str(RUN), 'first_step': 750, 'last_step': 1296,
         'accepted_steps': len(evidence), 'binary_sha256': observation['binary_sha256']}
    h['mean_max_rank_advance_seconds'] = mean(s['maximum_advance_ns'] for s in steps) / 1e9
    h['mean_rank_advance_seconds'] = mean(s['mean_advance_ns'] for s in steps) / 1e9
    h['mean_max_rank_full_step_seconds'] = mean(s['maximum_full_step_ns'] for s in steps) / 1e9
    h['sum_max_rank_full_step_seconds'] = sum(s['maximum_full_step_ns'] for s in steps) / 1e9
    h['critical_rank_step_phase_mean_seconds'] = {
        k: mean(s['critical_rank_phases'][k] for s in steps) / 1e9
        for k in steps[0]['critical_rank_phases']}
    assert math.isclose(sum(h['critical_rank_step_phase_mean_seconds'].values()),
                        h['mean_max_rank_full_step_seconds'], abs_tol=1e-9)
    h['rank_mean_stage_seconds'] = {
        str(stage['id']): mean(next(t['mean_ns'] for t in e['stages']
                                   if t['id'] == stage['id']) for e in evidence) / 1e9
        for stage in evidence[0]['stages']}
    h['rank_mean_stage_unaccounted_seconds'] = h['mean_rank_advance_seconds'] - sum(h['rank_mean_stage_seconds'].values())
    for output, source in [('rank_mean_candidate_and_solve_seconds', 'rank_mean_candidate_phases'),
                           ('rank_mean_solve_detail_seconds', 'rank_mean_solve_detail_ns')]:
        h[output] = {k: mean(s[source][k] for s in steps) / 1e9 for k in steps[0][source]}
    h['mean_calls_by_kind'] = {}
    h['mean_candidate_evaluations'] = {}
    for s in steps:
        assert all(len(values) == 1 for values in s['solve_calls_by_kind'].values())
        assert all(len(values) == 1 for values in s['work_counts'].values())
    for kind in steps[0]['solve_calls_by_kind']:
        h['mean_calls_by_kind'][kind] = mean(s['solve_calls_by_kind'][kind][0] for s in steps)
    for kind in steps[0]['work_counts']:
        h['mean_candidate_evaluations'][kind] = mean(s['work_counts'][kind][0] for s in steps)
    counters = ('linear_iterations', 'pressure_energy_refinement_solve_calls',
                'blocking_collectives', 'structured_messages', 'structured_bytes',
                'ibm_messages', 'ibm_bytes', 'exact_numeric_refills', 'coarse_numeric_refills')
    h['mean_work_counts'] = {key: mean(e[key] for e in evidence) for key in counters}
    h['mean_reduction_engine_seconds'] = mean(e['reduction_ns'] for e in evidence) / 1e9
    with (RUN / 'performance.csv').open() as stream:
        performance = list(csv.DictReader(stream))
    covered = collections.defaultdict(set)
    for row in performance:
        covered[int(row['step'])].add(int(row['rank']))
    assert all(ranks == set(range(128)) for ranks in covered.values())
    rank_zero = [row for row in performance if row['rank'] == '0']
    assert len(rank_zero) == len(evidence)
    h['rank0_mean_control_collectives'] = {
        key: mean(int(row[key]) for row in rank_zero)
        for key in ('reported_collective_subtotal', 'structured_control', 'ibm_control')}
    h['mean_pressure_iterations'] = mean(sum(t['iterations'] for t in e['pressure']) for e in evidence)
    h['mean_refinement_iterations'] = mean(sum(t['iterations'] for t in e['pressure_energy_refinement']) for e in evidence)
    h['mean_all_momentum_iterations'] = h['mean_work_counts']['linear_iterations'] - h['mean_pressure_iterations'] - h['mean_refinement_iterations']
    h['mean_pressure_refinement_norm_breakdown_restarts'] = mean(
        sum(t['norm_breakdown_restarts'] for t in e['pressure'] + e['pressure_energy_refinement'])
        for e in evidence)
    h['refinement_histogram'] = dict(sorted(collections.Counter(
        e['pressure_energy_refinement_solve_calls'] for e in evidence).items()))
    h['time_bands'] = []
    for lo, hi in [(750, 849), (850, 949), (950, 1049), (1050, 1149), (1150, 1249), (1250, 1296)]:
        group = [e for e in evidence if lo <= e['step'] <= hi]
        h['time_bands'].append({'first': lo, 'last': hi, 'count': len(group),
            'mean_advance_seconds': mean(e['max_rank_step_ns'] for e in group) / 1e9,
            'mean_refinements': mean(e['pressure_energy_refinement_solve_calls'] for e in group),
            'mean_linear_iterations': mean(e['linear_iterations'] for e in group)})
    h['max_rank_rss_first_last_peak_bytes'] = [evidence[0]['max_rank_rss_bytes'],
        evidence[-1]['max_rank_rss_bytes'], max(e['max_rank_rss_bytes'] for e in evidence)]
    h['rss_around_output'] = [{k: e[k] for k in ('step', 'max_rank_rss_bytes', 'max_node_rss_bytes')}
                             for e in evidence if e['step'] in (999, 1000, 1001, 1002, 1296)]
    h['checkpoint_phases'] = [{'step': s['step'],
         'critical_rank_checkpoint_seconds': s['critical_rank_phases']['checkpoint_ns'] / 1e9,
         'critical_rank_visit_seconds': s['critical_rank_phases']['visit_ns'] / 1e9}
        for s in steps if s['step'] % 500 == 0]
    c = validate_coast()
    c['grid_cells'] = [456, 256, 104]
    c['domain_extent_m'] = [0.4, 0.2, 0.06]
    with Path(c['evidence']).open() as stream:
        modules = [json.loads(line) for line in stream]
    detailed = collections.defaultdict(float)
    calls = collections.defaultdict(float)
    for row in modules:
        if row.get('kind') != 'module' or row['step'] not in range(2, 7):
            continue
        key = row['module'] + '@' + row['route']
        detailed[key] += row['critical_rank_self_seconds'] / 5
        calls[key] += row['critical_rank_call_count'] / 5
    c['critical_rank_self_module_mean_seconds'] = dict(sorted(detailed.items(), key=lambda kv: -kv[1]))
    c['critical_rank_module_mean_calls'] = dict(calls)
    hashes = {str(p): digest(p) for p in [RUN / 'performance.csv', RUN / 'evidence.jsonl',
        RUN / 'health.csv', RUN / 'RUN.meta', Path(str(RUN) + '.launch.log'),
        Path(c['evidence']), Path(c['evidence']).parent / 'monitor/timing.dat',
        Path(c['evidence']).parent / 'input.d',
        Path(__file__), BENCH / 'analyze_results.py', BENCH / 'analyze_strict_results.py',
        HERE.parents[2] / 'tools/v04_performance_observe.py']}
    result = {'schema': 'HUNDUN_COAST_DIAGNOSTIC_MODULE_REVIEW_V1',
              'same_case_speed_ratio_permitted': False,
              'reason': 'Different Z grid/domain, different trajectory windows; module scopes also differ.',
              'hundun': h, 'ordinary_compressible_coast_historical': c, 'source_hashes': hashes}
    with arguments.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'hundun': {k: v for k, v in h.items() if k != 'checkpoint_phases'}, 'coast': {k: c[k] for k in (
        'statistics', 'category_mean_seconds', 'mean_halo_messages', 'mean_halo_endpoint_bytes',
        'mean_critical_rank_collective_calls', 'reported_final_pressure_iterations',
        'maximum_timed_mass_balance_relative_residual')}}, indent=2))


if __name__ == '__main__':
    main()
