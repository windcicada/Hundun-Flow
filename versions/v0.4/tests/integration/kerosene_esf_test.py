#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""624CF gas assets with two ESF fields, Vreman and native rank-change history."""
from pathlib import Path
import json,shutil,subprocess,sys,tempfile
program,checker,launcher,validator=sys.argv[1:]
fixture=Path(__file__).resolve().parents[1]/'fixtures'/'ker-flow'
def call(args):
 result=subprocess.run(list(map(str,args)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                       universal_newlines=True,timeout=60)
 if result.returncode:raise RuntimeError(' '.join(map(str,args))+'\n'+result.stdout[-7000:])
 return result.stdout
def run(ranks,case,out,steps,*extra):
 call([launcher,'-n',ranks,program,'run',case,'--output',out,'--steps',steps,
       '--output-interval',0,'--restart-interval',steps,'--diagnostics-interval',1,*extra])
with tempfile.TemporaryDirectory(prefix='ke-') as folder:
 root=Path(folder);case=root/'case';shutil.copytree(fixture,case)
 model=json.loads((case/'case.json').read_text())
 model['reaction'].update(model='esf_tpdf',
   ensemble=dict(fields=2,seed=1234,
     initial_species_offsets=[1e-4,0,0,0,0,0,-1e-4,0,0,0,0,0],tcr=dict(mode='off')),
   mixing=dict(c_z=.001,turbulent_schmidt=.7))
 model['turbulence']=dict(model='vreman')
 (case/'case.json').write_text(json.dumps(model))
 seed=root/'seed'
 run(1,case,seed,1,'--initial-state','100000,1200,0,0,0,.001,.05,.005,.01,.2,.02')
 manifest=next((seed/'Restart').glob('generation-*/manifest.bin'))
 call([sys.executable,validator,'runtime',seed/'evidence.jsonl'])
 for ranks in (1,2,4):
  out=root/str(ranks);run(ranks,case,out,2,'--restart',seed/'Restart')
  call([sys.executable,validator,'runtime',out/'evidence.jsonl','--run-start-manifest',manifest])
  rows=[json.loads(line) for line in (out/'diagnostics.jsonl').read_text().splitlines()]
  assert [r['step'] for r in rows]==[2,3]
  for row in rows:
   assert row['payload']['dt']==1e-5
   balance=row['payload']['composition_balance']
   assert len(balance['species'])==7 and len(balance['elements'])==4
   assert all(x['relative_defect']<1e-6 for k in ('species','elements') for x in balance[k])
  if ranks>1:
   comparison=call([launcher,'-n',4,checker,case,root/'1/Restart',out/'Restart'])
   summary=[x for x in comparison.splitlines() if x.startswith('kerosene_esf_compare ')]
   assert len(summary)==1 and 'passed=1' in summary[0]
   print(summary[0])
 print('two-field kerosene ESF: full-state 1/2/4 continuation and physical/element budgets passed')
