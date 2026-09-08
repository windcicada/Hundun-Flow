import csv,hashlib,json,math,os,shutil,subprocess
from pathlib import Path
base=Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52');a=base/'coupled-stage46-20260909';run=base/'long-coupled-35000-20260909';snap=a/'startup-prefix'
svc='hundun-re3900-coupled-long-20260909.service'
state=subprocess.check_output(['systemctl','--user','show',svc,'-p','MainPID','-p','ActiveState','-p','SubState'],universal_newlines=True)
assert 'ActiveState=active' in state and 'SubState=running' in state and 'MainPID=0' not in state
snap.mkdir()
for p in run.iterdir():
 if p.is_file() and (p.suffix in ('.csv','.jsonl') or p.name in ('RUN.meta','mg-layout.meta')):shutil.copy2(str(p),str(snap/p.name))
raw=(snap/'evidence.jsonl').read_bytes();shutil.copy2(str(snap/'evidence.jsonl'),str(snap/'evidence.raw.jsonl'))
end=raw.rfind(b'\n')+1;assert end>0;(snap/'evidence.jsonl').write_bytes(raw[:end])
meta=dict(l.split(' ',1) for l in (snap/'RUN.meta').read_text().splitlines() if ' ' in l)
expected=dict(starting_step='10510',requested_steps='24490',expected_ranks='128',starting_sample_steps='0',statistics_epoch_start_step='7000',statistics_sampling_start_step='17001',restart_requires_recovery='0',restart_method_recovery='0',observation_schema='6',observe_mg_cost='1',observe_fgmres_recovery='1',trace_cell_count='0',candidate_head='83db41df2ec9c351c71e3aebe812509ddfd7262b')
for k,v in expected.items():assert meta[k]==v,(k,meta[k],v)
assert meta['executable_sha256']==hashlib.sha256((a/'frozen/v04_thin_domain_runner').read_bytes()).hexdigest()
rows=[r for r in csv.DictReader((snap/'health.csv').read_text().splitlines()) if None not in r and all(v is not None for v in r.values())]
assert len(rows)>=3
assert [int(r['step']) for r in rows]==list(range(10511,10511+len(rows)))
maxima=dict.fromkeys(['eos','continuity','energy','closed_mass','gauge'],0.)
for r in rows:
 assert int(r['fluid_count'])==5994352 and int(r['solid_placeholder_count'])==75920
 assert math.isfinite(float(r['committed_cfl_out'])) and float(r['committed_cfl_out'])<=float(r['committed_cfl_limit'])
 assert r['bdf_order']=='2' and r['attempts']=='1' and r['retry']=='0' and r['restart_recovery']=='0' and r['terminal_audit_present']=='1'
 for k in maxima:
  x=float(r[k]);assert math.isfinite(x) and 0<=x<=float(r[k+'_tolerance']);maxima[k]=max(maxima[k],x)
evidence=[json.loads(l) for l in (snap/'evidence.jsonl').read_text().splitlines()]
assert all(e['candidate_identity']['head']==expected['candidate_head'] and e['heap_allocations']==0 and not e['retry'] for e in evidence)
anchor=evidence[0]['run_start'];assert anchor['previous_step']==10510
assert anchor['history']==dict(source_format_version=3,source_signature=12213963202598979269,target_signature=12213963202598979269,policy='require_compatible')
g=(base/'integrated-stage46-20260908/source-10510/Restart/current').read_text().strip()
with (a/'startup-runtime.log').open('x') as f:subprocess.check_call(['python3',str(a/'frozen/v04_evidence_validate.py'),'runtime',str(snap/'evidence.jsonl'),'--run-start-manifest',str(base/'integrated-stage46-20260908/source-10510/Restart'/g/'manifest.bin')],stdout=f,stderr=subprocess.STDOUT)
with (a/'startup-observer.log').open('x') as f:subprocess.check_call(['python3',str(a/'frozen/v04_solver_observe.py'),str(snap),'--allow-partial','--output',str(a/'startup-observation.json')],stdout=f,stderr=subprocess.STDOUT)
o=json.loads((a/'startup-observation.json').read_text());assert not o['complete'] and o['validated_steps']>=3 and o['expected_rank_count']==128
assert o['issues']==['missing/truncated expected step {}'.format(10511+o['validated_steps'])],o['issues']
subprocess.check_call(['sha256sum','--check',str(a/'SOURCE-LONG.sha256')],cwd=str(base/'integrated-stage46-20260908/source-10510'),stdout=subprocess.DEVNULL)
subprocess.check_call(['sha256sum','--check',str(a/'FROZEN.sha256')],cwd=str(a/'frozen'),stdout=subprocess.DEVNULL)
ranks=[];rss=0
for l in subprocess.check_output(['ps','-eo','pid=,stat=,comm=,rss='],universal_newlines=True).splitlines():
 pid,st,comm,mem=l.split()
 if comm in ('v04_thin_domain','mpirun','mpiexec'):
  assert svc in Path('/proc/'+pid+'/cgroup').read_text(),pid
  assert not st.startswith(('T','Z')),st
  if comm=='v04_thin_domain':
   assert Path('/proc/'+pid+'/exe').resolve()==(a/'frozen/v04_thin_domain_runner').resolve(),pid
   ranks.append(pid);rss+=int(mem)
assert len(ranks)==128
result=dict(service=svc,service_state=state,source_head=expected['candidate_head'],source_step=10510,target_step=35000,run_root=str(run),startup_verified=True,long_run_complete=False,health_step_range=[10511,int(rows[-1]['step'])],bdf2_steps=len(rows),retries=0,terminal_maxima=maxima,restart_anchor=anchor,observer_validated_steps=o['validated_steps'],observer_issues=o['issues'],ranks=128,all_rank_executables_match_frozen=True,ranks_rss_sum_kib=rss,next_checkpoint_and_visit_step=11000,source_and_frozen_hashes_unchanged=True,evidence_partial_tail_bytes_excluded=len(raw)-end,checksums={n:hashlib.sha256((a/n).read_bytes()).hexdigest() for n in ['PILOT_ACCEPTED.json','FROZEN.sha256','SOURCE-LONG.sha256','startup-observation.json']},limitations=['Read-only asynchronous snapshot; CSV streams can stop on different steps.','Raw evidence retained; only complete JSONL records submitted to runtime validator.','RSS is a single process-sum snapshot, not leak or memory-trend acceptance.'])
(a/'STARTUP_VERIFIED.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
