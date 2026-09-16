#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
"""Audit legacy cell geometry before conservative transfer; source files stay intact."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_grid(path):
    with path.open() as stream:
        stream.readline()
        dims = tuple(map(int, stream.readline().split()))
        scales = tuple(map(float, stream.readline().split()))
        if len(dims) != 3 or min(dims) < 3 or len(scales) != 2:
            raise ValueError('invalid legacy grid header')
        count = math.prod(n + 1 for n in dims)
        xyz = np.fromfile(stream, sep=' ', count=3 * count)
    if xyz.size != 3 * count or not np.isfinite(xyz).all():
        raise ValueError('invalid coordinate payload')
    return xyz.reshape((3,) + tuple(n + 1 for n in dims[::-1])), scales


def corners(xyz):
    """Cell i=2:l has vertices i-1,i; k,j,i storage, component first."""
    if xyz.ndim != 4 or xyz.shape[0] != 3 or min(xyz.shape[1:]) < 4:
        raise ValueError('coordinate shape')
    return {(i,j,k): xyz[:, slice(1+k,-2+k), slice(1+j,-2+j), slice(1+i,-2+i)]
            for k in (0,1) for j in (0,1) for i in (0,1)}


def jacobian(c, point=(.5,.5,.5)):
    """Determinant of the trilinear reference-cell map; unit cube coordinates."""
    derivatives = []
    for axis in range(3):
        d = np.zeros_like(c[0,0,0])
        for corner, value in c.items():
            if corner[axis]:
                continue
            opposite = list(corner)
            opposite[axis] = 1
            weight = 1.
            for other in range(3):
                if other != axis:
                    weight *= point[other] if corner[other] else 1.-point[other]
            d += weight * (c[tuple(opposite)] - value)
        derivatives.append(d)
    a,b,c = derivatives
    return a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+c[0]*(a[1]*b[2]-a[2]*b[1])


def metrics(xyz):
    c = corners(xyz)
    centers = sum(c.values()) / 8.
    # This is the source center-Jacobian volume identity, evaluated in FP64.
    midpoint = jacobian(c)
    integrated = np.zeros_like(midpoint)
    gauss = (.5-.5/math.sqrt(3), .5+.5/math.sqrt(3))
    gauss_min = np.full_like(midpoint, np.inf)
    gauss_max = np.full_like(midpoint, -np.inf)
    for z in gauss:
        for y in gauss:
            for x in gauss:
                sample = jacobian(c, (x,y,z))
                integrated += sample / 8.
                gauss_min = np.minimum(gauss_min, sample)
                gauss_max = np.maximum(gauss_max, sample)
    if not all(np.isfinite(v).all() for v in (centers,midpoint,integrated)):
        raise ValueError('nonfinite cell geometry')
    volume = np.abs(midpoint)
    if (volume == 0).any():
        raise ValueError('zero source center Jacobian')
    # A y-directed IBM ray uses one x,z pair for each fixed logical i,k.
    # Peak-to-peak over INTERNAL j rows proves the required line invariance
    # independently of boundary/halo treatment of the source j=1 anchor.
    drift = np.ptp(centers[[0,2]], axis=2)
    return dict(cells=int(midpoint.size),
        center_bounds_m=[[float(v.min()),float(v.max())] for v in centers],
        jacobian_range_m3=[float(midpoint.min()),float(midpoint.max())],
        negative_jacobians=int((midpoint<0).sum()),
        gauss_orientation_changes=int(((gauss_min<=0)&(gauss_max>=0)).sum()),
        source_abs_jacobian_sum_m3=float(volume.sum()),
        trilinear_abs_volume_sum_m3=float(np.abs(integrated).sum()),
        quadrature_vs_midpoint_relative_max=float((np.abs(integrated-midpoint)/volume).max()),
        scan_xz_drift_max_m=[float(v.max()) for v in drift],
        logical_ik_lines=int(drift[0].size),
        lines_x_drift_above_1um=int((drift[0]>1e-6).sum()),
        lines_z_drift_above_1um=int((drift[1]>1e-6).sum()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists')
    manifest = json.loads(args.manifest.read_text())
    entries = [x for x in manifest['files'] if x['path'].startswith('Decomp/grid_vv.')]
    if not entries:
        parser.error('manifest has no legacy grid entries')
    rows = []
    for entry in sorted(entries,key=lambda x:x['path']):
        path = args.case/entry['path']
        if path.stat().st_size != entry['size'] or digest(path) != entry['sha256']:
            raise ValueError('input identity: '+entry['path'])
        xyz,scales = read_grid(path)
        row = metrics(xyz)
        row.update(path=entry['path'], sha256=entry['sha256'], scales=scales)
        if digest(path) != entry['sha256']:
            raise ValueError('input changed: '+entry['path'])
        rows.append(row)
        if len(rows)%16==0: print('geometry audited',len(rows),flush=True)
    result = dict(schema='hundun_legacy_geometry_inventory_v1',
        scope='All interior logical cells including solids; source FP64 center Jacobian and trilinear quadrature, before IBM classification.',
        manifest_sha256=digest(args.manifest), script_sha256=digest(Path(__file__)),
        cells=sum(r['cells'] for r in rows),
        source_abs_jacobian_sum_m3=math.fsum(r['source_abs_jacobian_sum_m3'] for r in rows),
        trilinear_abs_volume_sum_m3=math.fsum(r['trilinear_abs_volume_sum_m3'] for r in rows),
        quadrature_vs_midpoint_relative_max=max(r['quadrature_vs_midpoint_relative_max'] for r in rows),
        scan_xz_drift_max_m=[max(r['scan_xz_drift_max_m'][i] for r in rows) for i in range(2)],
        logical_ik_lines=sum(r['logical_ik_lines'] for r in rows),
        lines_x_drift_above_1um=sum(r['lines_x_drift_above_1um'] for r in rows),
        lines_z_drift_above_1um=sum(r['lines_z_drift_above_1um'] for r in rows),
        gauss_orientation_changes=sum(r['gauss_orientation_changes'] for r in rows),rows=rows)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':
    main()
