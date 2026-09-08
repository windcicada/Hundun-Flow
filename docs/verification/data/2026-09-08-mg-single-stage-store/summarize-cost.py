#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Post-completion MG cost attribution for the one frozen 9500-9510 pilot; not a speed gate."""
import csv,json,hashlib,math
from pathlib import Path
from collections import defaultdict
import argparse
from datetime import datetime, timezone

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--run-root", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--observation", type=Path, required=True)
args = parser.parse_args()
if args.output.exists():
    parser.error("refusing to overwrite existing output: " + str(args.output))
root = args.run_root.resolve()

paths=[root/'RUN.meta',root/'performance.csv',root/'evidence.jsonl']+[root/('solver-rank-%d.csv'%r) for r in range(128)]
stamps={str(p):(p.stat().st_size,p.stat().st_mtime_ns) for p in paths}
def sha(p):
 d=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):d.update(b)
 return d.hexdigest()
hashes={p.name:sha(p) for p in paths}
meta={}
for line in (root/'RUN.meta').read_text().splitlines()[1:-1]:
 k,v=line.split(' ',1);meta[k]=v
assert meta['candidate_head']=='ed9b09ab2a019a7957030844e5437421351fd044'
assert meta['observation_schema']=='5' and meta['expected_ranks']=='128'
floats={'dt','linear_initial','linear_final','linear_rhs_norm','linear_atol','linear_rtol','linear_residual_limit','baseline_continuity','baseline_energy','selected_continuity','selected_energy','selected_alpha'}
def read_csv(p):
 with p.open() as f:
  for raw in csv.DictReader(f):
   assert raw['source_meta_sha256']==hashes['RUN.meta']
   yield {k:(v if k=='source_meta_sha256' else float(v) if k in floats else int(v)) for k,v in raw.items()}
performance=list(read_csv(root/'performance.csv'))
pby=defaultdict(list)
for r in performance:
 pby[r['step']].append(r)
 assert r['dropped_loops']==0
 assert sum(r[k] for k in ('advance_ns','observables_ns','visit_ns','evidence_resources_ns','csv_ns','checkpoint_ns'))==r['full_step_ns']
assert sorted(pby)==list(range(9501,9511))
for v in pby.values():assert len(v)==128 and {r['rank'] for r in v}==set(range(128))
with (root/'evidence.jsonl').open() as f:evidence=[json.loads(l) for l in f]
assert len(evidence)==10
eby={r['step']:r for r in evidence}
assert set(eby)==set(pby)
for e in evidence:
 assert e['coupling']=='SIMPLE' and not e['retry'] and not e['temporal_method_fallback']
 assert e['candidate_identity']['head']==meta['candidate_head']
lby=defaultdict(list)
for rank in range(128):
 for r in read_csv(root/('solver-rank-%d.csv'%rank)):
  assert r['rank']==rank and r['attempt']==1 and r['scalar_coupling_sweep']==1
  assert r['attempt_status']==0 and r['dropped_loops']==0
  key=tuple(r[k] for k in ('step','attempt','scalar_coupling_sweep','corrector','refinement','kind'))
  lby[key].append(r)
countnames=('invoked','iterations','A_calls','M_calls','baseline_candidates','extrapolated_candidates','ladder_candidates','incomplete_candidates')
jsoncounts=('norm_breakdown_restarts','convergence_audits','convergence_rejections','reduction_calls','recycle_offered_directions','recycle_retained_directions','recycle_operator_applies','recycle_reduction_calls','recycle_cycle_corrections','recycle_capture_vector_passes','recycle_capture_cycle_attempts','recycle_capture_reduction_calls','recycle_capture_blocking_operations')
jsonflags=('recycle_projection_attempted','recycle_projection_accepted')
loops=[]
for key,rows in sorted(lby.items()):
 assert len(rows)==128 and {r['rank'] for r in rows}==set(range(128))
 row=rows[0]
 invariant=list(countnames)+list(floats)+['linear_criterion_valid','globalization_valid','attempt_status']
 for k in invariant:assert all(r[k]==row[k] for r in rows),(key,k)
 assert row['invoked']==1 and row['linear_criterion_valid']==1 and row['globalization_valid']==1
 assert row['linear_residual_limit']==max(row['linear_atol'],row['linear_rtol']*row['linear_rhs_norm'])
 e=eby[row['step']]
 if row['refinement']==0:
  matches=[x for x in e['pressure'] if x['corrector']==row['corrector']]
  jpath='pressure[corrector=%d]'%row['corrector']
 else:
  assert row['corrector']==2
  matches=[x for x in e['pressure_energy_refinement'] if x['ordinal']==row['refinement']]
  jpath='pressure_energy_refinement[ordinal=%d]'%row['refinement']
 assert len(matches)==1
 solve=matches[0]
 assert solve['status_code']==0 and solve['termination']=='converged'
 for ck,jk in [('iterations','iterations'),('A_calls','operator_applies'),('M_calls','preconditioner_applies'),('linear_initial','initial_true_residual'),('linear_final','final_true_residual')]:
  assert row[ck]==solve[jk],(key,ck)
 kind=('pressure','diagonal','spatial')[row['kind']]
 loop={'step':row['step'],'attempt':row['attempt'],'scalar_coupling_sweep':row['scalar_coupling_sweep'],'corrector':row['corrector'],'refinement':row['refinement'],'kind':kind,
       'group':'C%d-r%d-%s'%(row['corrector'],row['refinement'],kind),
       'rank_count':128,'counts':{k:row[k] for k in countnames},
       'rank_mean_seconds':{k[:-3]:math.fsum(r[k] for r in rows)/128/1e9 for k in row if k.endswith('_ns')},
       'rank_max_solve_seconds':max(r['solve_ns'] for r in rows)/1e9,
       'linear':{k:row[k] for k in ('linear_initial','linear_final','linear_rhs_norm','linear_atol','linear_rtol','linear_residual_limit')},
       'initial_over_rhs':row['linear_initial']/row['linear_rhs_norm'],
       'linear_contraction':row['linear_final']/row['linear_initial'],
       'final_over_limit':row['linear_final']/row['linear_residual_limit'],
       'absolute_floor_dominates':row['linear_atol']>=row['linear_rtol']*row['linear_rhs_norm'],
       'physical':{k:row[k] for k in ('baseline_continuity','baseline_energy','selected_continuity','selected_energy','selected_alpha')},
       'continuity_contraction':row['selected_continuity']/row['baseline_continuity'],
       'energy_contraction':row['selected_energy']/row['baseline_energy'],
       'evidence_path':jpath,'evidence_counts':{k:solve[k] for k in jsoncounts},
       'evidence_flags':{k:solve[k] for k in jsonflags}}
 for k in ('target_generation','collective_lineage'):
  if k in solve:loop[k]=solve[k]
 loops.append(loop)
for s,ps in pby.items():
 ls=[l for l in loops if l['step']==s]
 assert len(ls)==eby[s]['pressure_solve_calls']+eby[s]['pressure_energy_refinement_solve_calls']
 for p in ps:
  rank=p['rank'];rr=[r for key,rs in lby.items() if key[0]==s for r in rs if r['rank']==rank]
  assert sum(r['A_ns'] for r in rr)==p['A_apply_ns']
  assert sum(r['M_ns'] for r in rr)==p['M_apply_ns']
  for kind,name in enumerate(('pressure','diagonal','spatial')):
   assert sum(r['invoked'] for r in rr if r['kind']==kind)==p[name+'_calls']
def stats(v):
 v=list(v)
 return {'min':min(v),'mean':math.fsum(v)/len(v),'max':max(v)} if v else None
def group_summary(ls):
 ns=len(ls)
 timenames=ls[0]['rank_mean_seconds']
 counts={k:sum(l['counts'][k] for l in ls) for k in countnames}
 times={k:math.fsum(l['rank_mean_seconds'][k] for l in ls) for k in timenames}
 return {'logical_loops':ns,'counts':counts,'iterations_per_loop':counts['iterations']/ns,
         'rank_mean_seconds_total':times,
         'rank_mean_solve_seconds_per_loop':times['solve']/ns,
         'rank_mean_A_seconds_per_apply':times['A']/counts['A_calls'] if counts['A_calls'] else None,
         'rank_mean_M_seconds_per_apply':times['M']/counts['M_calls'] if counts['M_calls'] else None,
         'rank_mean_candidate_seconds_per_loop':times['candidate']/ns,
         'linear_ranges':{k:stats(l['linear'][k] for l in ls) for k in ls[0]['linear']},
         'initial_over_rhs':stats(l['initial_over_rhs'] for l in ls),
         'initial_equals_rhs_loops':sum(l['linear']['linear_initial']==l['linear']['linear_rhs_norm'] for l in ls),
         'absolute_floor_dominates_loops':sum(l['absolute_floor_dominates'] for l in ls),
         'linear_contraction':stats(l['linear_contraction'] for l in ls),
         'final_over_limit':stats(l['final_over_limit'] for l in ls),
         'continuity_contraction':stats(l['continuity_contraction'] for l in ls),
         'energy_contraction':stats(l['energy_contraction'] for l in ls),
         'evidence_counts':{k:sum(l['evidence_counts'][k] for l in ls) for k in jsoncounts},
         'norm_breakdown_restarts_per_loop':math.fsum(l['evidence_counts']['norm_breakdown_restarts'] for l in ls)/ns,
         'projection_attempted_loops':sum(l['evidence_flags']['recycle_projection_attempted'] for l in ls),
         'projection_accepted_loops':sum(l['evidence_flags']['recycle_projection_accepted'] for l in ls)}
def window(first,last):
 ss=list(range(first,last+1))
 ls=[l for l in loops if first<=l['step']<=last]
 ps=[p for p in performance if first<=p['step']<=last]
 timings={k[:-3]:math.fsum(p[k] for p in ps)/len(ps)/1e9 for k in ps[0] if k.endswith('_ns')}
 groups={g:group_summary([l for l in ls if l['group']==g]) for g in sorted({l['group'] for l in ls})}
 return {'step_range':[first,last],'steps':len(ss),'rank_count':128,'logical_loops':len(ls),
         'max_rank_full_step_seconds':stats(max(p['full_step_ns'] for p in pby[s])/1e9 for s in ss),
         'max_rank_advance_seconds':stats(max(p['advance_ns'] for p in pby[s])/1e9 for s in ss),
         'rank_step_mean_seconds':timings,
         'rank_step_mean_loop_seconds':{k:math.fsum(l['rank_mean_seconds'][k] for l in ls)/len(ss) for k in ls[0]['rank_mean_seconds']},
         'rank_mean_advance_fraction':{k:timings[k]/timings['advance'] for k in ('krylov','candidate','A_apply','M_apply','schur_prepare','final_momentum','terminal_metrics','boundary_ledger')},
         'counts':{k:sum(l['counts'][k] for l in ls) for k in countnames},
         'groups':groups,
         'evidence_counts':{k:sum(l['evidence_counts'][k] for l in ls) for k in jsoncounts}}
perstep=[]
for s in sorted(pby):
 ps=pby[s];ls=[l for l in loops if l['step']==s]
 perstep.append({'step':s,'logical_loops':len(ls),'max_rank_full_step_seconds':max(p['full_step_ns'] for p in ps)/1e9,
                 'max_rank_advance_seconds':max(p['advance_ns'] for p in ps)/1e9,
                 'rank_mean_krylov_seconds':math.fsum(p['krylov_ns'] for p in ps)/128/1e9,
                 'rank_mean_candidate_seconds':math.fsum(p['candidate_ns'] for p in ps)/128/1e9,
                 'norm_breakdown_restarts':sum(l['evidence_counts']['norm_breakdown_restarts'] for l in ls)})
assert all((p.stat().st_size,p.stat().st_mtime_ns)==stamps[str(p)] for p in paths)
out={'schema':'HUNDUN_MG_COST_SUMMARY_V1','scope':'Independent cost aggregation only; formal observation/runtime/checkpoint/Visit validation belongs to root agent.',
     'source_root':str(root),'source_sha256':hashes,'candidate_identity':evidence[0]['candidate_identity'],
     'run_start':evidence[0]['run_start'],'observation_schema':5,'expected_rank_count':128,
     'analysis_range':[9501,9510],'cost_input_consistency_checks_passed':True,
     'formal_speed_comparison':False,'scientific_acceptance':False,
     'methods':['Logical loop key=(step,attempt,scalar_coupling_sweep,corrector,refinement,kind). Every key has 128 unique ranks; integer counts and global residuals counted once, not 128 times.',
                'Time per logical loop=mean of 128 rank-local inclusive timings. Per-apply cost=sum of rank-mean A/M seconds divided by once-per-logical-loop A/M calls.',
                'Full/advance statistic=max rank separately for each step, then arithmetic mean across selected steps; independent maxima are never summed.',
                'All reported A/M, MG, structured communication and Arnoldi costs are nested where recorded, not additive external wall-time components. MG copy is within refill.',
                'Evidence join: r0 uses pressure.corrector; refinement>0 uses pressure_energy_refinement.ordinal. Join is valid here because every step has exactly one attempt and one scalar sweep; counts and initial/final residuals match every rank.',
                '9501 is the restart first step; 9502-9510 are subsequent steps, not a claim of statistically steady flow.',
                'Restart clears ephemeral warm-start authority. Same step numbers in old long run are not a strict equal-initial-guess work baseline.',
                'No CFD/MPI, solver change, repeated configuration, or speed/physics acceptance experiment was performed by this cost analysis.'],
     'windows':{'restart_first_step':window(9501,9501),'subsequent_steps':window(9502,9510),'all10':window(9501,9510)},
     'per_step':perstep,'logical_loop_evidence_joins':loops}

# Read the independently validated V5 layout/ledger. Bind every consumed file
# to its recorded digest before deriving any MG cost fraction.
observed=json.loads(args.observation.read_text())
assert observed['complete'] and observed['schema']=='HUNDUN_LOOP_OBSERVATION_V5'
assert observed['validated_steps']==10 and observed['loop_count']==len(loops)
for name,digest in observed['files_sha256'].items():
 assert sha(root/name)==digest,name
geometry={(r['rank'],r['level']):{k:r[k] for k in ('global_shape','local_shape','coarsening','line_axis_mask')}
          for r in observed['mg']['levels']}
mg_rows=[]
for rank in range(128):
 for r in read_csv(root/('mg-rank-%d.csv'%rank)):
  assert r['rank']==rank and r['complete']==1 and r['initialized']==1
  mg_rows.append(r)
for label,w in out['windows'].items():
 first,last=w['step_range'];steps=w['steps']
 selected=[r for r in mg_rows if first<=r['step']<=last]
 totals=[r for r in selected if r['level']==-1]
 denom=128*steps
 aggregate={k:math.fsum(r[k] for r in totals)/denom/1e9 for k in ('apply_ns','reduction_ns')}
 levels=[]
 for level in sorted({r['level'] for r in selected if r['level']>=0}):
  rr=[r for r in selected if r['level']==level]
  assert len(rr)==denom
  times={k:math.fsum(r[k] for r in rr)/denom/1e9 for k in rr[0] if k.endswith('_ns')}
  calls={k:sum(r[k] for r in rr)/denom for k in rr[0] if k.endswith('_calls') or k=='visits'}
  shapes=sorted({tuple(geometry[r['rank'],level]['local_shape']) for r in rr})
  levels.append(dict(level=level,global_shape=geometry[0,level]['global_shape'],
     local_shapes=shapes,rank_step_mean_seconds=times,rank_step_mean_calls=calls))
 phase_names=('pre_smooth_ns','residual_ns','restriction_ns','prolongation_ns','post_smooth_ns','terminal_ns')
 phase={k:math.fsum(l['rank_step_mean_seconds'][k] for l in levels) for k in phase_names}
 nested={k:math.fsum(l['rank_step_mean_seconds'][k] for l in levels)
         for k in ('direct_mpi_ns','halo_wait_ns','halo_control_ns')}
 assert sum(phase.values())<=aggregate['apply_ns']
 w['mg']={'rank_step_mean_seconds':aggregate,'phases_seconds':phase,
          'phases_fraction_of_mg_apply':{k:v/aggregate['apply_ns'] for k,v in phase.items()},
          'nested_communication_seconds':nested,'levels':levels,
          'outside_six_phases_seconds':aggregate['apply_ns']-sum(phase.values()),
          'scope':'Native pressure M only. Six phases are disjoint; communication is nested. No per-loop full hierarchy reconstruction.'}
out['generation_metadata']={'script_path':str(Path(__file__).resolve()),'script_sha256':sha(Path(__file__).resolve()),
 'generated_utc':datetime.now(timezone.utc).isoformat(),
 'derived_from':'criterion-observation-20260908/summarize-cost.py; same once-only logical counts and rank-time means',
 'observation_sha256':sha(args.observation),
 'generation_method':'Python JSON preserves uint64 identities; post-hoc analysis only, no CFD/MPI or input writes'}
with args.output.open('x') as stream:
 json.dump(out,stream,indent=2,allow_nan=False);stream.write('\n')
w=out['windows']['subsequent_steps']
print(json.dumps({'window':w['step_range'],'advance_max_rank_seconds':w['max_rank_advance_seconds'],
 'rank_step_mean_seconds':w['rank_step_mean_seconds'],
 'mg':w['mg'],
 'groups':{g:{k:v[k] for k in ('logical_loops','iterations_per_loop','rank_mean_M_seconds_per_apply')} for g,v in w['groups'].items()},
 'script_sha256':out['generation_metadata']['script_sha256'],
 'output_sha256':sha(args.output)},indent=2))
