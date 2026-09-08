import csv, hashlib, importlib.util, json, math, os, shutil, subprocess
from pathlib import Path
base=Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52')
audit=base/'integrated-stage46-20260908'; old=base/'long-observed-35000-20260908'; dest=audit/'source-10500'
service='hundun-re3900-observed-long-20260908.service'
assert subprocess.check_output(['systemctl','--user','show',service,'-p','MainPID','--value'],universal_newlines=True).strip()=='142362'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
 return h.hexdigest()
dest.mkdir(); (dest/'Restart').mkdir()
generation=(old/'Restart/current').read_text().strip(); assert generation.startswith('generation-10500-')
donefile='step-00000000000000010500.complete'
lines=(old/donefile).read_text().splitlines(); assert lines[0]=='HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1' and lines[-1]=='end'
done=dict(l.split(' ',1) for l in lines[1:-1]); assert done['step']=='10500' and done['restart_generation']==generation
shutil.copytree(str(old/'Restart'/generation),str(dest/'Restart'/generation))
for rel in ['Restart/current',donefile,done['statistics'],done['accumulator']]: shutil.copy2(str(old/rel),str(dest/rel))
hashes={}
for p in sorted(dest.rglob('*')):
 if p.is_file():
  rel=str(p.relative_to(dest)); hashes[rel]=sha(p); assert hashes[rel]==sha(old/rel),rel
(audit/'SOURCE.sha256').write_text(''.join(h+'  '+r+'\n' for r,h in hashes.items()))
# Preserve stopped logs verbatim; no truncation of original partial trailing records.
(audit/'old-prefix').mkdir()
for name in ['health.csv','force.csv','conservation.csv','probe.csv','performance.csv','RUN.meta','evidence.jsonl']: shutil.copy2(str(old/name),str(audit/'old-prefix'/name))
module_spec=importlib.util.spec_from_file_location('ev',str(base/'long-handoff-20260908/frozen/v04_evidence_validate.py'))
v=importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(v)
m=v.load_v04_restart_manifest(dest/'Restart'/generation/'manifest.bin')
assert m['step']==10500 and m['source_format_version']==3 and m['method_history_signature']==12213963202598979269 and not m['backward_euler_recovery']
ranks=list((dest/'Restart'/generation).glob('rank-*.bin')); assert len(ranks)==128
stats=json.loads((dest/done['statistics']).read_text()); assert stats['snapshot_step']==10500 and stats['sample_steps']==0
assert stats['statistics_epoch']==dict(start_step=7000,development_steps=10000,sampling_start_step=17001,reset_reason='method_recovery')
rows=[r for r in csv.DictReader((audit/'old-prefix/health.csv').read_text().splitlines()) if None not in r and all(x is not None for x in r.values())]
assert [int(r['step']) for r in rows]==list(range(10001,int(rows[-1]['step'])+1))
maxima=dict.fromkeys(['eos','continuity','energy','closed_mass','gauge'],0.)
for r in rows:
 assert r['attempts']=='1' and r['retry']=='0' and r['bdf_order']=='2' and r['restart_recovery']=='0' and r['terminal_audit_present']=='1'
 for k in maxima:
  x=float(r[k]); assert math.isfinite(x) and 0<=x<=float(r[k+'_tolerance']); maxima[k]=max(maxima[k],x)
visits=list((old/'Visit').glob('step-00000000000000010500-rank-*.vtr')); assert len(visits)==128 and all(x.stat().st_size>0 for x in visits)
assert (old/'Visit/step-00000000000000010500.visit').is_file()
result=dict(source=str(dest),generation=generation,manifest=m,files=len(hashes),rank_files=128,visit_files=128,statistics_epoch=stats['statistics_epoch'],last_complete_health_step=int(rows[-1]['step']),terminal_maxima=maxima,old_service=service,old_main_pid=142362,replay_from=10500,old_run_complete=False)
(audit/'SOURCE_VALIDATED.json').write_text(json.dumps(result,indent=2)+'\n'); os.sync(); print(json.dumps(result,indent=2))
