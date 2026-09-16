#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Analytic cell volumes and source-index alignment for legacy grid audits."""
import importlib.util
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('grid',root/'tools/v04_grid_inventory.py')
grid = importlib.util.module_from_spec(spec);spec.loader.exec_module(grid)
z,y,x = np.meshgrid(np.arange(7.),np.arange(6.),np.arange(5.),indexing='ij')
# Affine tilted cell: determinant is 2*3*4 even with off-diagonal shear.
a = np.array([2*x+.25*y,3*y,4*z])
r = grid.metrics(a)
assert r['cells']==24 and r['source_abs_jacobian_sum_m3']==576
assert abs(r['trilinear_abs_volume_sum_m3']-576)<1e-12
assert r['quadrature_vs_midpoint_relative_max']<1e-14
assert r['scan_xz_drift_max_m']==[.5,0]
assert r['center_bounds_m']==[[3.375,5.875],[4.5,10.5],[6.,18.]]
# Trilinear warped map X=x(1+a*z), Y=y(1+b*z), Z=z.
# J=(1+a*z)(1+b*z); exact volume differs from midpoint by a*b/12.
a,b = .2,.3
xyz=np.array([x*(1+a*z), y*(1+b*z), z])
r=grid.metrics(xyz)
exact=6*sum(1+(a+b)*(k+.5)+a*b*(k*k+k+1/3) for k in range(1,5))
assert abs(r['trilinear_abs_volume_sum_m3']-exact)<1e-12
assert abs(r['trilinear_abs_volume_sum_m3']-r['source_abs_jacobian_sum_m3']-24*a*b/12)<1e-12
# Reversing a coordinate preserves absolute inventory and reports orientation.
xyz[0]*=-1
assert grid.metrics(xyz)['negative_jacobians']==24
# Collapsed source cell is rejected before publishing an inventory.
xyz[0]=0
try:grid.metrics(xyz)
except ValueError as e:assert 'zero source' in str(e)
else:raise AssertionError('collapsed cell admitted')
print('legacy geometry: affine/warped volume, source indexing, orientation and collapse passed')

# The CLI verifies original bytes and preserves both inputs and an existing report.
import hashlib,json,subprocess,sys,tempfile
with tempfile.TemporaryDirectory(prefix='grid-') as folder:
    folder=Path(folder);case=folder/'case';(case/'Decomp').mkdir(parents=True)
    path=case/'Decomp/grid_vv.000'
    xyz=np.array([2*x+.25*y,3*y,4*z])
    path.write_text('grid\n4 5 6\n1 1\n'+' '.join(map(str,xyz.reshape(-1)))+'\n')
    original=path.read_bytes()
    sha=hashlib.sha256(original).hexdigest()
    manifest=folder/'manifest.json'
    manifest.write_text(json.dumps(dict(files=[dict(path='Decomp/grid_vv.000',size=len(original),sha256=sha)])))
    output=folder/'grid.json'
    command=[sys.executable,str(root/'tools/v04_grid_inventory.py'),'--case',str(case),
             '--manifest',str(manifest),'--output',str(output)]
    def call():return subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    r=call();assert r.returncode==0,r.stdout
    assert json.loads(output.read_text())['cells']==24 and path.read_bytes()==original
    saved=output.read_bytes()
    r=call();assert r.returncode!=0 and b'output already exists' in r.stdout
    assert output.read_bytes()==saved
    output.unlink();path.write_bytes(original+b' ')
    r=call();assert r.returncode!=0 and b'input identity' in r.stdout and not output.exists()
    path.write_bytes(original[:-10])
    try:grid.read_grid(path)
    except ValueError as e:assert 'payload' in str(e)
    else:raise AssertionError('truncated coordinates admitted')
print('legacy geometry CLI: source identity, output preservation and truncation passed')
