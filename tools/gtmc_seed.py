#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
"""Prepare the frozen Cartesian GTMC reference for explicit dyn711 import."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import numpy as np


def digest(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''):value.update(chunk)
    return value.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('base','transfer','source','manifest','case','output'):
        parser.add_argument('--'+key,required=True,type=Path)
    args=parser.parse_args()
    for out in (args.case,args.output):
        if out.exists():raise ValueError('output already exists: '+str(out))
    if args.case.resolve()==args.output.resolve():raise ValueError('outputs require separate paths')
    inputs=[p.resolve() for p in (args.base,args.transfer,args.source,args.manifest.parent)]
    outputs=[args.case.resolve(),args.output.resolve()]
    for out in outputs:
        if any(out==source or source in out.parents or out in source.parents for source in inputs):
            raise ValueError('output overlaps an input directory')
    if outputs[0] in outputs[1].parents or outputs[1] in outputs[0].parents:
        raise ValueError('outputs require independent directories')
    pack=json.loads((args.transfer/'pack.json').read_text())
    authority=json.loads(args.manifest.read_text())['files']
    header=(args.transfer/'state.txt').read_text().split()
    if header[:2]!=['HUNDUN_PDF_TRANSFER','1']:raise ValueError('expected frozen V1 transfer')
    nx,ny,nz=map(int,header[2:5]);step=int(header[5]);nf,ns=map(int,header[9:11])
    if nf!=4 or ns!=7 or len(header)!=18:raise ValueError('GTMC four-field seven-species reference required')
    model=json.loads((args.base/'case.json').read_text())
    if model['mesh']['exact_cells']!=[nx,ny,nz]:raise ValueError('grid identity mismatch')
    for name,entry in pack['artifacts'].items():
        path=args.transfer/name
        if path.stat().st_size!=entry['bytes'] or digest(path)!=entry['sha256']:
            raise ValueError('transfer identity mismatch: '+name)
    args.case.mkdir(parents=True);args.output.mkdir(parents=True)
    for path in args.base.iterdir():
        if not path.is_file():raise ValueError('expected flat native case assets')
        if path.name!='case.json':shutil.copyfile(path,args.case/path.name)
    for name in pack['artifacts']:
        if name!='state.txt':
            shutil.copyfile(args.transfer/name,args.output/name)
            if digest(args.output/name)!=pack['artifacts'][name]['sha256']:
                raise ValueError('copied transfer identity mismatch: '+name)
    shape=(nz,ny,nx)
    values=np.memmap(str(args.output/'passive.f64'),dtype='<f8',mode='w+',shape=shape)
    coverage=np.zeros(shape,dtype='u1')
    source_entries={}
    def read(name):
        data=(args.source/name).read_bytes();identity=authority[name]
        if len(data)!=identity['bytes'] or hashlib.sha256(data).hexdigest()!=identity['sha256']:
            raise ValueError('source identity mismatch: '+name)
        source_entries[name]=identity
        return data
    def records(data):
        offset=0;result=[]
        while offset<len(data):
            if offset+8>len(data):raise ValueError('truncated record')
            size,=struct.unpack_from('<i',data,offset)
            if size<0 or offset+size+8>len(data) or struct.unpack_from('<i',data,offset+size+4)[0]!=size:
                raise ValueError('invalid Fortran record')
            result.append(memoryview(data)[offset+4:offset+size+4]);offset+=size+8
        return result
    for rank in range(128):
        lines=read('runtime_mesh/block_%03d.dat'%rank).decode().splitlines()
        dims=tuple(map(int,lines[2].split()[1:]));bounds=tuple(map(int,lines[3].split()[1:]))
        region=np.s_[bounds[4]:bounds[5]+1,bounds[2]:bounds[3]+1,bounds[0]:bounds[1]+1]
        flow=records(read('restart.%03d'%rank))
        source_step,source_time,source_dt=struct.unpack('<iff',flow[0])
        if len(flow)!=9 or (source_step,source_time,source_dt)!=(step,float(header[6]),float(header[7])):
            raise ValueError('source epoch mismatch')
        if flow[5]!=flow[8]:raise ValueError('mixture fraction record alias mismatch')
        storage=tuple(n+2 for n in dims[::-1])
        field=np.frombuffer(flow[5],dtype='<f4').reshape(storage)[2:-2,2:-2,2:-2]
        if not np.isfinite(field).all():raise ValueError('nonfinite source mixture fraction')
        values[region]=field;coverage[region]+=1
    if not np.all(coverage==1):raise ValueError('source blocks must cover each target cell once')
    values.flush()
    minimum,maximum=float(values.min()),float(values.max())
    names=header[11:]
    oxygen=.23290921795842306
    model['reaction']['ensemble']['tcr']=dict(model='dyn711_v1',mode='experimental',fuel='CH4',
        weak_rate_threshold=1e-30,mixture_fraction='Z',oxidizer_oxygen_mass_fraction=oxygen)
    model['transported_scalars'].append(dict(stable_name='Z',role='passive_scalar',
        molecular_schmidt=.7,turbulent_schmidt=.7))
    boundaries=list(model['boundaries'].values())+[p['boundary'] for p in model['patch_inlets']['patches']]
    for boundary in boundaries:
        scalar={v['stable_name']:v for v in boundary['scalars']}
        if not scalar:continue
        for key in ('value','backflow_value'):
            fuel=scalar['CH4'][key]
            expected={name:(fuel if name=='CH4' else (1-fuel)*oxygen if name=='O2' else 0.) for name in scalar}
            if fuel not in (0.,1.) or any(scalar[name][key]!=value for name,value in expected.items()):
                raise ValueError('GTMC reference requires pure fuel or air boundary streams')
        prototype=scalar['CH4']
        boundary['scalars'].append(dict(stable_name='Z',kind=prototype['kind'],value=prototype['value'],
            backflow_kind=prototype['backflow_kind'],backflow_value=prototype['backflow_value']))
    model['time']['convective_cfl_definition']='directional_max'
    model['time']['convective_cfl']=.30;model['time']['convective_cfl_margin']=.05
    (args.case/'case.json').write_text(json.dumps(model,indent=2)+'\n')
    (args.output/'state.txt').write_text('HUNDUN_PDF_TRANSFER 2\n'+' '.join(header[2:11])+' 1\n'+' '.join(names)+'\nZ\n')
    report=dict(schema='gtmc_dyn711_seed_v1',source_pack_sha256=digest(args.transfer/'pack.json'),
        source_manifest_sha256=digest(args.manifest),source_files=source_entries,
        source_step=step,source_time=source_time,source_dt=source_dt,grid=[nx,ny,nz],
        passive=dict(name='Z',source='restart variable 5, checked against duplicate nvf record',
                     minimum=minimum,maximum=maximum,sha256=digest(args.output/'passive.f64')),
        model_history='explicit target initialization required',case_sha256=digest(args.case/'case.json'),
        geometry='same frozen Cartesian cells and IBM; coordinate mapping is identity')
    (args.output/'prepare.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_files'},indent=2))

if __name__=='__main__':main()
