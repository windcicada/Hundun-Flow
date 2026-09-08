import csv,hashlib,json
from pathlib import Path
base=Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52'); audit=base/'coupled-stage46-20260909'
new=base/'pilot-coupled-10500-10510-20260909'; old=base/'control-mainline-10500-10510-20260908'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def done(root):
 assert 'COMPLETED steps=10 final_step=10510 samples=0' in Path(str(root)+'.log').read_text()
 lines=(root/'step-00000000000000010510.complete').read_text().splitlines()
 assert lines[0]=='HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1' and lines[-1]=='end'
 d=dict(l.split(' ',1) for l in lines[1:-1]); assert d['restart_generation']==(root/'Restart/current').read_text().strip();return d
nd,od=done(new),done(old)
assert {k:v for k,v in nd.items() if k!='restart_generation'}=={k:v for k,v in od.items() if k!='restart_generation'}
ng,og=[root/'Restart'/d['restart_generation'] for root,d in [(new,nd),(old,od)]]
names=sorted(p.name for p in ng.iterdir() if p.is_file());assert names==sorted(p.name for p in og.iterdir() if p.is_file());assert len(names)==129
files={}
for name in names:
 x,y=sha(ng/name),sha(og/name);assert x==y,name;files['Restart/generation/'+name]=x
for rel in [nd['statistics'],nd['accumulator']]+[str(p.relative_to(new)) for p in sorted((new/'Visit').glob('*')) if p.is_file()]:
 x,y=sha(new/rel),sha(old/rel);assert x==y,rel;files[rel]=x
csv_rows={}
process_bindings={}
for name in ['health.csv','force.csv','probe.csv','conservation.csv']:
 with (new/name).open(newline='') as f:nr=list(csv.DictReader(f))
 with (old/name).open(newline='') as f:orr=list(csv.DictReader(f))
 assert len(nr)==len(orr)>0
 for x,y in zip(nr,orr):
  assert None not in x and None not in y and all(v is not None for v in x.values()) and all(v is not None for v in y.values())
  ignored={'certificate_terminal_state','certificate_state'} if name=='force.csv' else set()
  if ignored:
   process_bindings[x['step']]={k:[x[k],y[k]] for k in sorted(ignored)}
   assert all(int(x[k])>0 and int(y[k])>0 for k in ignored)
  assert {k:v for k,v in x.items() if not k.endswith('_ns') and k not in ignored}=={k:v for k,v in y.items() if not k.endswith('_ns') and k not in ignored},(name,x['step'])
 csv_rows[name]=len(nr)
result=dict(same_restart_source=10500,final_step=10510,checkpoint_and_visit_payloads_byte_identical=True,files_compared=len(files),physical_and_stable_revision_csv_rows_equal=csv_rows,process_bound_certificate_values=process_bindings,completion_equal_except_generation=True,sha256=files,notes=['One run per program from the same immutable source; this is compatibility evidence, not a performance speedup claim.','Generation names and candidate provenance remain distinct.', 'Force certificate_terminal_state hashes rank-local field bindings including base addresses (solver_piso.cpp:451,642,8544,9063); certificate_state incorporates it (solver_ibm_force.cpp:132). They are process-bound authority, not portable physical checksums; both retained, not rewritten.'])
(audit/'CONTROL_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))
