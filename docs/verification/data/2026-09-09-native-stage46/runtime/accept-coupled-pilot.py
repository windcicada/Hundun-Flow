import datetime,hashlib,json,subprocess
from pathlib import Path
b=Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52');a=b/'coupled-stage46-20260909';r=b/'pilot-coupled-10500-10510-20260909';f=a/'frozen';source=b/'integrated-stage46-20260908/source-10500';svc='hundun-re3900-coupled-pilot-20260909.service'
state=subprocess.check_output(['systemctl','--user','show',svc,'-p','MainPID','-p','ActiveState','-p','Result','-p','ExecMainStatus'],universal_newlines=True)
assert 'MainPID=0' in state and 'Result=success' in state and 'ExecMainStatus=0' in state,state
assert 'COMPLETED steps=10 final_step=10510 samples=0' in Path(str(r)+'.log').read_text()
def run(name,args):
 with (a/name).open('x') as out:subprocess.check_call(args,stdout=out,stderr=subprocess.STDOUT)
run('pilot-audit.log',['python3',str(a/'verify-pilot.py'),str(r),str(a/'PILOT_AUDIT.json')])
g=(source/'Restart/current').read_text().strip()
run('pilot-runtime.log',['python3',str(f/'v04_evidence_validate.py'),'runtime',str(r/'evidence.jsonl'),'--run-start-manifest',str(source/'Restart'/g/'manifest.bin')])
run('pilot-observer.log',['python3',str(f/'v04_solver_observe.py'),str(r),'--output',str(a/'pilot-observation.json')])
run('pilot-control.log',['python3',str(a/'compare-control.py')])
o=json.loads((a/'pilot-observation.json').read_text());assert o['complete'] and o['validated_steps']==10 and o['expected_rank_count']==128 and not o['issues'],o
m=dict(l.split(' ',1) for l in (r/'RUN.meta').read_text().splitlines() if ' ' in l);local=json.loads((a/'LOCAL_ACCEPTED.json').read_text());assert m['candidate_head']==local['source_head'] and m['executable_sha256']==local['runner_sha256']
for root,name in [(f,'FROZEN.sha256'),(source,'SOURCE.sha256')]:subprocess.check_call(['sha256sum','--check',str(a/name)],cwd=str(root),stdout=subprocess.DEVNULL)
c=json.loads((a/'CONTROL_COMPARISON.json').read_text());assert c['checkpoint_and_visit_payloads_byte_identical'] and c['files_compared']==260
result=dict(accepted_at=datetime.datetime.now().isoformat(),source_head=local['source_head'],executable_sha256=local['runner_sha256'],service=svc,service_state=state,complete=True,step_range=[10501,10510],ranks=128,retries=0,bdf2_steps=10,control_payloads_byte_identical=260,physical_and_stable_revision_csv_equal=True,runtime_validator_passed=True,observer_complete=True,source_and_frozen_hashes_unchanged=True,long_run_complete=False,checksums={n:hashlib.sha256((a/n).read_bytes()).hexdigest() for n in ['PILOT_AUDIT.json','CONTROL_COMPARISON.json','pilot-observation.json','FROZEN.sha256','LOCAL_ACCEPTED.json']})
(a/'PILOT_ACCEPTED.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
