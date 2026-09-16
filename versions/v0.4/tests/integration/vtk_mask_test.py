#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""VTK mask import: exact owned cells, ghost isolation and atomic publication."""
from pathlib import Path
import struct,tempfile,subprocess,sys,tarfile,json,hashlib
root=Path(__file__).resolve().parents[4]
script=root/'tools/v04_vtk_mask.py'

def vtk(path, changed=None):
 dims=(4,5,6);n=120;points=[];mask=[];owned=[];xyz=[]
 for k in range(6):
  for j in range(5):
   for i in range(4):
    p=(float(i),float(j),float(k));points.extend(p)
    inside=0<i<3 and 0<j<4 and 0<k<5
    v=float((i+j+k)%2) if inside else .5
    if inside:owned.append(int(v));xyz.extend(p)
    mask.append(v)
 if changed=='mask':mask[(1*5+1)*4+1]=1-mask[(1*5+1)*4+1]
 if changed=='fraction':mask[(1*5+1)*4+1]=.5
 if changed=='point':points[3*((1*5+1)*4+1)]+=.25
 if changed=='nan':points[3*((1*5+1)*4+1)]=float('nan')
 data=b'# vtk DataFile Version 2.0\nfixture\nBINARY\nDATASET STRUCTURED_GRID\nDIMENSIONS 4 5 6\nPOINTS 120 float\n'
 data+=struct.pack('>360f',*points)+b'\nPOINT_DATA 120\nSCALARS avtGhostNodes unsigned_char 1\nLOOKUP_TABLE default\n'+bytes(n)
 data+=b'\nVECTORS Velocity float\n'+bytes(n*12)
 data+=b'\nSCALARS 07_IBM_cell_type float 1\nLOOKUP_TABLE default\n'+struct.pack('>120f',*mask)+b'\n'
 if changed=='truncate':data=data[:-9]
 path.write_bytes(data)
 return bytes(owned),hashlib.sha256(struct.pack('<%df'%len(xyz),*xyz)).hexdigest()

for optimized in (False,True):
 with tempfile.TemporaryDirectory(prefix='mask-') as folder:
  folder=Path(folder);a=folder/'solution.00000010.domain.000.vtk';b=folder/'solution.00000020.domain.000.vtk'
  expected,coord=vtk(a);vtk(b);output=folder/'mask.tgz'
  cmd=[sys.executable]+(['-O'] if optimized else [])+[str(script),'--visit',str(folder),'--steps','10','20','--ranks','1','--dims','4','5','6','--output',str(output)]
  def run():return subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  r=run();assert r.returncode==0,r.stdout
  with tarfile.open(str(output)) as t:
   assert t.extractfile('000.bin').read()==expected
   m=json.load(t.extractfile('manifest.json'));s=m['ranks'][0]['before']
   assert s['owned_coordinates_sha256']==coord and s['marker_fractional_storage_values']==96
  saved=output.read_bytes();assert run().returncode!=0 and output.read_bytes()==saved
  output.unlink()
  for change in ['mask','fraction','point','nan','truncate']:
   vtk(b,change);r=run();assert r.returncode!=0,(change,r.stdout)
   assert not output.exists() and not list(folder.glob('.mask-*'))
print('VTK static masks: exact owned data, ghost isolation, source changes and atomic output passed with/without -O')
