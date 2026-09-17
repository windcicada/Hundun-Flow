#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Real JL4 output, independent VTK payload reads and unchanged Restart states."""
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

binary, probe, fixture, mpi, work = map(Path, sys.argv[1:])
if work.exists():
    shutil.rmtree(work)
shutil.copytree(fixture, work)
aliases = dict(U='Velocity', pi='PressureGauge', rho='Density', T='Temperature', h='Enthalpy',
    rho_mean_eos='DensityMeanEOS', rho_pdf_mean='DensityStatistical', rho_field0='DensityField0',
    T_field0='TemperatureField0', h_field0='EnthalpyField0', Y_field0='MassFractionsField0')
required = set(list(aliases.values())[5:])

def call(args, name):
    with (work/(name+'.log')).open('w') as log:
        result = subprocess.run(list(map(str,args)), stdout=log, stderr=log)
    assert result.returncode == 0, (name, (work/(name+'.log')).read_text()[-10000:])

def run(name, ranks, steps, output=1, restart=None, fmt='legacy'):
    root=work.with_name(work.name+'-'+name)
    if root.exists():
        shutil.rmtree(root)
    args=[mpi,'--oversubscribe','--bind-to','none','-n',ranks,binary,'run',work,
          '--output',root,'--steps',steps,'--output-interval',output,'--restart-interval',steps]
    if fmt=='xml':
        args+=['--visit-format','xml']
    args+=['--restart',restart/'Restart'] if restart else \
        ['--initial-state','100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001']
    call(args,name)
    return root

def checkpoint(root):
    path=root/'Restart'
    return path/(path/'current').read_text().strip()

def identical(a,b):
    left,right=checkpoint(a),checkpoint(b)
    files={p.name:p.read_bytes() for p in left.iterdir() if p.is_file()}
    other={p.name:p.read_bytes() for p in right.iterdir() if p.is_file()}
    assert files==other, 'output changed a complete native checkpoint'

def legacy(path):
    stream=io.BytesIO(path.read_bytes())
    def line():
        return stream.readline().decode().strip()
    assert line()=='# vtk DataFile Version 3.0'
    line()
    assert line()=='BINARY' and line()=='DATASET STRUCTURED_GRID'
    assert line().startswith('DIMENSIONS ')
    words=line().split();assert words[0]=='POINTS' and words[2]=='double'
    count=int(words[1])
    xyz=list(struct.iter_unpack('>3d',stream.read(count*24)))
    assert line()=='' and line()=='POINT_DATA '+str(count)
    fields={}
    while stream.tell()<len(stream.getbuffer()):
        words=line().split()
        if not words:continue
        if words[0]=='FIELD':
            assert words[2]=='1'
            name,width,points,dtype=line().split();assert int(points)==count
        elif words[0]=='VECTORS':
            _,name,dtype=words;width=3
        else:
            assert words[0]=='SCALARS'
            _,name,dtype,width=words
            assert line()=='LOOKUP_TABLE default'
        width=int(width);code='d' if dtype=='double' else 'B'
        fields[name]=list(struct.iter_unpack('>'+code*width,stream.read(count*width*struct.calcsize(code))))
        assert line()==''
    assert required <= fields.keys()
    assert len(fields['MassFractionsField0'][0])==7
    assert all(math.isfinite(v) for rows in fields.values() for row in rows for v in row)
    return xyz,fields

def legacy_frame(root, step):
    owned={};duplicates=[]
    for path in sorted((root/'Visit').glob('s%dr*.vtk'%step)):
        points,fields=legacy(path)
        for i,xyz in enumerate(points):
            values={name:rows[i] for name,rows in fields.items() if name!='vtkGhostType'}
            if fields['vtkGhostType'][i][0]:
                duplicates.append((xyz,values))
            else:
                assert xyz not in owned
                owned[xyz]=values
    assert len(owned)==512 and duplicates
    for xyz,values in duplicates:
        assert values==owned[xyz], 'halo output differs from its owner'
    return owned

def xml_frame(root):
    result={}
    for path in (root/'Visit').glob('*.vti'):
        data=path.read_bytes();split=data.index(b'<AppendedData')
        tree=ET.fromstring(data[:split]+b'</VTKFile>')
        piece=tree.find('./ImageData/Piece')
        extent=list(map(int,piece.attrib['Extent'].split()));begin=extent[::2]
        shape=[extent[j+1]-extent[j] for j in (0,2,4)];count=math.prod(shape)
        base=data.index(b'_',data.index(b'>',split))+1
        fields={}
        for node in piece.find('CellData'):
            width=int(node.attrib['NumberOfComponents']);offset=base+int(node.attrib['offset'])
            size,=struct.unpack_from('<Q',data,offset);assert size==count*width*8
            fields[aliases.get(node.attrib['Name'],node.attrib['Name'])]=list(
                struct.iter_unpack('<'+'d'*width,data[offset+8:offset+8+size]))
        assert required <= fields.keys()
        for i in range(count):
            coords=(i%shape[0],i//shape[0]%shape[1],i//(shape[0]*shape[1]))
            xyz=tuple((b+x+.5)/8 for b,x in zip(begin,coords))
            assert xyz not in result
            result[xyz]={name:rows[i] for name,rows in fields.items()}
    assert len(result)==512
    return result

fresh=run('fresh',2,2)
quiet=run('quiet',2,2,output=0)
identical(fresh,quiet)
resumed=run('resumed',4,1,restart=fresh)
xml=run('xml',4,1,restart=fresh,fmt='xml')
identical(resumed,xml)
for ranks,root in ((1,fresh),(2,fresh),(4,resumed)):
    call([mpi,'--oversubscribe','--bind-to','none','-n',ranks,probe,work,root/'Restart'],'probe'+str(ranks))
legacy_frame(fresh,1)
legacy_frame(fresh,2)
actual=legacy_frame(resumed,3)
for xyz,fields in xml_frame(xml).items():
    assert all(actual[xyz][name]==values for name,values in fields.items())
# The ESF-only path uses the same diagnostics without a velocity-gradient pass.
original=(work/'case.json').read_text()
model=json.loads(original);model['turbulence']={'model':'none'}
(work/'case.json').write_text(json.dumps(model))
plain=run('plain',2,1)
legacy_frame(plain,1)
call([mpi,'--oversubscribe','--bind-to','none','-n',2,probe,work,plain/'Restart'],'plain-probe')
(work/'case.json').write_text(original)
report=dict(passed=True, cells=512, source_ranks=2, restart_ranks=4,
    closure_checks=[(work/('probe'+str(r)+'.log')).read_text() for r in (1,2,4)],
    output_quiet_restart='byte-identical', xml_legacy_restart='byte-identical',
    vtk_payloads='identical by global coordinate, including overlap',
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
(work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
