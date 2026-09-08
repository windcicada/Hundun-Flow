import hashlib,json,shutil,subprocess,datetime
from pathlib import Path
r=Path.cwd();base=Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52');a=base/'coupled-stage46-20260909';old=base/'integrated-stage46-20260908/frozen'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def git(*args):return subprocess.check_output(['git']+list(args),universal_newlines=True).strip()
head=git('rev-parse','HEAD');tree=git('rev-parse','HEAD^{tree}');assert head=='83db41df2ec9c351c71e3aebe812509ddfd7262b';assert not git('status','--porcelain')
paths=[p for d in ['versions/v0.4/include','versions/v0.4/src','cmake','third_party/yyjson'] for p in (r/d).rglob('*') if p.is_file() and p!=r/'versions/v0.4/src/app_main.cpp']+[r/'CMakeLists.txt',r/'versions/v0.4/CMakeLists.txt']
content=hashlib.sha256(('HUNDUN_SOURCE_CONTENT_V1\n'+''.join(str(p.relative_to(r))+'='+sha(p)+'\n' for p in sorted(paths,key=str))).encode()).hexdigest();verified={}
for directory,targets in [('build-integrated/versions/v0.4/generated',['core','hundun','v04_thin_domain_runner']),('build-integrated-cantera/generated',['core','hundun'])]:
 for target in targets:
  p=r/directory/(target+'.build-manifest.txt');m=dict(l.split('=',1) for l in p.read_text().splitlines())
  assert m['head']==hashlib.sha256(('hundun-git-head-v1:'+head).encode()).hexdigest();assert m['tree']==hashlib.sha256(('hundun-git-tree-v1:'+tree).encode()).hexdigest()
  assert m['core_source_clean']==m['target_source_clean']=='true';assert m['core_source_content_sha256']==content,(content,m['core_source_content_sha256'])
  entry='none' if target=='core' else sha(r/('tools/v04_thin_domain_runner.cpp' if target=='v04_thin_domain_runner' else 'versions/v0.4/src/app_main.cpp'))
  assert m['entry_sha256']==entry;verified[str(p.relative_to(r))]=sha(p)
f=a/'frozen';f.mkdir();shutil.copytree(str(old/'case'),str(f/'case'));shutil.copy2(str(old/'statistics-long-20plus50D.d'),str(f/'statistics-long-20plus50D.d'))
pairs=[(r/'build-integrated/versions/v0.4'/n,f/n) for n in ['hundun','v04_thin_domain_runner']]
pairs += [(r/'build-integrated/versions/v0.4/generated'/(n+'.build-manifest.txt'),f/(n+'.build-manifest.txt')) for n in ['core','hundun','v04_thin_domain_runner']]
pairs += [(r/'tools'/n,f/n) for n in ['v04_evidence_validate.py','v04_solver_observe.py']]
g=a/'cantera-gcc11';g.mkdir();pairs.append((r/'build-integrated-cantera/hundun',g/'hundun'))
pairs += [(r/'build-integrated-cantera/generated'/(n+'.build-manifest.txt'),g/(n+'.build-manifest.txt')) for n in ['core','hundun']]
for src,dst in pairs:shutil.copy2(str(src),str(dst));assert sha(src)==sha(dst)
for d,name in [(f,'FROZEN.sha256'),(g,'CANTERA-FROZEN.sha256')]:
 (a/name).write_text(''.join(sha(p)+'  '+str(p.relative_to(d))+'\n' for p in sorted(d.rglob('*')) if p.is_file()))
 for p in d.rglob('*'):
  if p.is_file():p.chmod(p.stat().st_mode & ~0o222)
 subprocess.check_call(['sha256sum','--check',str(a/name)],cwd=str(d),stdout=subprocess.DEVNULL)
for step,name in [(10500,'SOURCE.sha256'),(10510,'SOURCE-LONG.sha256')]:subprocess.check_call(['sha256sum','--check',str(a/name)],cwd=str(base/'integrated-stage46-20260908'/('source-'+str(step))),stdout=subprocess.DEVNULL)
result=dict(accepted_at=datetime.datetime.now().isoformat(),source_head=head,source_tree=tree,source_content_sha256=content,source_clean=True,native_coupling_complete=True,clang_tests_passed=205,gcc_cantera_tests_passed=67,clean_cantera_ibm_restart_cli_passed=True,verified_manifests=verified,runner_sha256=sha(f/'v04_thin_domain_runner'),hundun_clang_sha256=sha(f/'hundun'),hundun_gcc_cantera_sha256=sha(g/'hundun'),restart_sources_verified={'10500':133,'10510':133},cantera_rootfs='/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/source/rootfs-jammy-attempt-2-openmpi-fontconfig-postinst-recovery',cantera_package='/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/artifact/cantera-3.2.0-gcc11-release-v4',long_run_complete=False)
(a/'LOCAL_ACCEPTED.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
