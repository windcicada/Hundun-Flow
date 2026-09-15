# SPDX-License-Identifier: Apache-2.0
"""Refine only chemistry controls; preserve native histories across ranks."""
from pathlib import Path
import copy,hashlib,json,shutil,subprocess,sys,tempfile

program,checker,launcher,validator=sys.argv[1:]
fixture=Path(__file__).resolve().parents[1]/'fixtures'/'esf-cn'
def call(args,success=True):
 result=subprocess.run([str(x) for x in args],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                       universal_newlines=True,timeout=60)
 if (result.returncode==0)!=success:
  raise RuntimeError(' '.join(map(str,args))+'\n'+result.stdout[-6000:])
 return result.stdout
def run(ranks,case,output,*extra,success=True):
 return call([launcher,'-n',ranks,program,'run',case,'--output',output,'--steps',1,
              '--output-interval',0,'--restart-interval',1,'--diagnostics-interval',0,*extra],success)
def hashes(root):
 return {str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest()
         for f in root.rglob('*') if f.is_file()}
with tempfile.TemporaryDirectory(prefix='chem-') as directory:
 root=Path(directory)
 for name,mode in [('mean','finite_rate_mean'),('pdf','esf_tpdf')]:
  source=root/(name+'0');target=root/(name+'1')
  shutil.copytree(fixture,source)
  model=json.loads((source/'case.json').read_text());model['reaction']['model']=mode
  if mode=='finite_rate_mean':
   model['reaction'].pop('ensemble');model['reaction'].pop('mixing')
  # Use the production reference-scaled equation audit in both cases.
  model['solver']['cold_stopping']=dict(reference_time=1.,momentum=1e-7,enthalpy=1e-10,species=1e-10)
  model['reaction']['chemistry_solver'].update(relative_tolerance=1e-10,absolute_tolerance=1e-14)
  (source/'case.json').write_text(json.dumps(model))
  shutil.copytree(source,target);fine=copy.deepcopy(model)
  fine['reaction']['chemistry_solver'].update(relative_tolerance=1e-12,absolute_tolerance=1e-18)
  (target/'case.json').write_text(json.dumps(fine))
  seed=root/(name+'s');run(1,source,seed,'--initial-state','101325,300,0,0,0,0.25')
  before=hashes(seed/'Restart');manifest=next((seed/'Restart').glob('generation-*/manifest.bin'))
  # Ordinary strict Restart continues to bind the complete source identity.
  run(1,target,root/(name+'strict'),'--restart',seed/'Restart',success=False)
  for ranks in [1,4]:
   output=root/(name+str(ranks)+'out')
   run(ranks,target,output,'--restart',seed/'Restart','--restart-refine-chemistry',source)
   rows=[json.loads(l) for l in (output/'evidence.jsonl').read_text().splitlines()]
   assert len(rows)==1 and rows[0]['step']==2 and not rows[0]['restart_recovery']
   history=rows[0]['run_start']['history']
   assert history['policy']=='refine_chemistry' and history['source_signature']==history['target_signature']
   assert history['chemistry_source_case']>0
   call([sys.executable,validator,'runtime',output/'evidence.jsonl','--run-start-manifest',manifest])
   for mutation in ['missing-source','wrong-policy','wrong-signature','recovery']:
    altered=copy.deepcopy(rows[0]);anchor=altered['run_start']['history']
    if mutation=='missing-source':anchor.pop('chemistry_source_case')
    elif mutation=='wrong-policy':anchor['policy']='require_compatible'
    elif mutation=='wrong-signature':anchor['target_signature']+=1
    else:altered['restart_recovery']=True
    evidence=root/'altered.jsonl';evidence.write_text(json.dumps(altered)+'\n')
    call([sys.executable,validator,'runtime',evidence,'--run-start-manifest',manifest],success=False)
   if checker!='-':
    text=call([launcher,'-n',ranks,checker,source,target,seed/'Restart'])
    assert 'chemistry_refinement exact=1 values=' in text
   # The refined case can subsequently use ordinary native Restart.
   resumed=root/(name+str(ranks)+'next')
   run(ranks,target,resumed,'--restart',output/'Restart')
   refined_manifest=next((output/'Restart').glob('generation-*/manifest.bin'))
   call([sys.executable,validator,'runtime',resumed/'evidence.jsonl','--run-start-manifest',refined_manifest])
  assert before==hashes(seed/'Restart')
  for index,change in enumerate(['looser','same','budget','time','transport','asset']):
   bad=root/(name+'b'+str(index));shutil.copytree(target,bad);candidate=copy.deepcopy(fine)
   controls=candidate['reaction']['chemistry_solver']
   if change=='looser':controls['relative_tolerance']=1e-8
   elif change=='same':controls.update(model['reaction']['chemistry_solver'])
   elif change=='budget':controls['maximum_internal_steps']=1000
   elif change=='time':candidate['time']['maximum_dt']=.002
   elif change=='transport':candidate['transported_scalars'][0]['molecular_schmidt']=1.1
   else:
    asset=bad/'mechanism.yaml';asset.write_text(asset.read_text()+'\n# changed asset identity\n')
    candidate['reaction']['mechanism_sha256']=hashlib.sha256(asset.read_bytes()).hexdigest()
   (bad/'case.json').write_text(json.dumps(candidate))
   run(4,bad,root/(name+'reject'+str(index)),'--restart',seed/'Restart',
       '--restart-refine-chemistry',source,success=False)
 print('mean and ESF:1/4-rank exact histories, native continuation,12 incompatible refinements rejected')
