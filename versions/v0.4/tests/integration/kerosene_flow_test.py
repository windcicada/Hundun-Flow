# SPDX-License-Identifier: Apache-2.0
"""Native kerosene mixture transport and original-history continuation."""
from pathlib import Path
import shutil,subprocess,sys,tempfile,json,re,hashlib
program,launcher,validator=sys.argv[1:]
fixture=Path(__file__).resolve().parents[1]/'fixtures'/'ker-flow'
def call(args,success=True):
 result=subprocess.run(list(map(str,args)),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                       universal_newlines=True,timeout=60)
 if (result.returncode==0)!=success:raise RuntimeError(' '.join(map(str,args))+'\n'+result.stdout[-7000:])
 return result.stdout
def native(n,case,out,steps,*args,success=True):
 return call([launcher,'-n',n,program,'run',case,'--output',out,'--steps',steps,
              '--output-interval',0,'--restart-interval',1,'--diagnostics-interval',0,*args],success)
with tempfile.TemporaryDirectory(prefix='ker-') as directory:
 root=Path(directory);case=root/'case';shutil.copytree(fixture,case)
 # Exercise the mixed-family dispatch with kerosene first in the property catalog.
 props=(case/'thermophysics.d').read_text()
 blocks=re.findall(r'^species .*?^end_species\n',props,re.M|re.S)
 assert len(blocks)==7 and blocks[-1].startswith('species C12H23\n')
 prefix=props[:props.index('species H2\n')]
 (case/'thermophysics.d').write_text(prefix+blocks[-1]+''.join(blocks[:-1])+'end\n')
 # The chemistry asset and property catalog share an ordered species identity.
 gas=(case/'gas.yaml').read_text()
 old='species: [H2, H2O, CO, CO2, O2, N2, C12H23]'
 assert gas.count(old)==1
 gas=gas.replace(old,'species: [C12H23, H2, H2O, CO, CO2, O2, N2]')
 (case/'gas.yaml').write_text(gas)
 config=json.loads((case/'case.json').read_text())
 config['reaction']['mechanism_sha256']=hashlib.sha256(gas.encode()).hexdigest()
 (case/'case.json').write_text(json.dumps(config))
 first=root/'one';native(1,case,first,1,'--initial-state','100000,1200,0,0,0,.001,.05,.005,.01,.2,.02')
 manifest=next((first/'Restart').glob('generation-*/manifest.bin'))
 call([sys.executable,validator,'runtime',first/'evidence.jsonl'])
 continued=root/'four';native(4,case,continued,2,'--restart',first/'Restart')
 rows=[json.loads(l) for l in (continued/'evidence.jsonl').read_text().splitlines()]
 assert [r['step'] for r in rows]==[2,3] and all(not r['restart_recovery'] for r in rows)
 call([sys.executable,validator,'runtime',continued/'evidence.jsonl','--run-start-manifest',manifest])
 original=(case/'thermophysics.d').read_text()
 for index,altered in enumerate([
  original.replace('transport_perry 33.3 12.8','transport_kerosene'),
  original.replace('transport_perry 33.3 12.8','transport_sutherland 1e-5 300 0 0.7'),
  original.replace('transport_kerosene','transport_kerosene 300'),
  original.replace('transport_kerosene','transport_kerosene_unknown')]):
  assert altered!=original
  (case/'thermophysics.d').write_text(altered)
  native(4,case,root/('bad'+str(index)),1,'--initial-state','100000,1200,0,0,0,.001,.05,.005,.01,.2,.02',success=False)
 (case/'thermophysics.d').write_text(original.replace('transport_kerosene','transport_perry 369.85 41.92'))
 native(4,case,root/'identity',1,'--restart',first/'Restart',success=False)
 print('kerosene mixture: native three steps, 1->4 Restart, four invalid inputs and changed physical identity checked')
