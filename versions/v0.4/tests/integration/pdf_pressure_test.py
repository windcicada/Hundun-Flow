#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""PDF import separates signed mechanical pressure from fixed thermodynamic p0."""
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

binary,application,fixture,mpi,work=map(Path,sys.argv[1:])
if work.exists():shutil.rmtree(work)
shutil.copytree(fixture,work)
model=json.loads((work/'case.json').read_text())
model['reaction']['model']='esf_tpdf'
model['reaction']['ensemble']=dict(fields=2,seed=1234,initial_species_offsets=[0.]*12,tcr=dict(mode='off'))
model['reaction']['mixing']=dict(c_z=.25,turbulent_schmidt=.7)
p0=790216.58
model['flow']['thermodynamic_pressure_pa']=p0
(work/'case.json').write_text(json.dumps(model))
thermal=(work/'thermophysics.d').read_text()
species=re.findall(r'^species (\S+)$',thermal,re.M)
assert 'N2' in species and len(species)==7
fractions=[.232 if name=='O2' else .768 if name=='N2' else 0. for name in species]
h=gas=0.
temperature=300.
for name,y in zip(species,fractions):
    block=thermal.split('species '+name+'\n')[1].split('end_species')[0]
    mw=float(re.search(r'molecular_weight (\S+)',block).group(1))
    a=list(map(float,re.search(r'nasa7_low (.*)',block).group(1).split()))
    gas+=8314.46261815324*y/mw
    h+=8314.46261815324*y/mw*(sum(a[j]*temperature**(j+1)/(j+1) for j in range(5))+a[5])

def pack(path,values):
    path.write_bytes(struct.pack('<'+'d'*len(values),*values))

def transfer(name,pressure,eos_pressure):
    root=work/name;root.mkdir()
    (root/'state.txt').write_text('HUNDUN_PDF_TRANSFER 1\n8 8 8 17 .001 1e-5 100000 2 7\n'+' '.join(species)+'\n')
    pack(root/'flow.f64',[value for i in range(512) for value in (0.,0.,0.,pressure(i))])
    pack(root/'rho_ref.f64',[eos_pressure(i)/(gas*temperature) for i in range(512)])
    (root/'fluid.u8').write_bytes(bytes([1])*512)
    for f in range(2):pack(root/('pdf%d.f64'%f),(fractions+[h])*512)
    return root

def run(source,label,ranks,success=True):
    output=work.with_name(work.name+'-'+label)
    if output.exists():shutil.rmtree(output)
    with (work/(label+'.log')).open('w') as log:
        result=subprocess.run(list(map(str,[mpi,'--oversubscribe','--bind-to','none','-n',ranks,binary,work,source,output])),stdout=log,stderr=log)
    text=(work/(label+'.log')).read_text()
    assert (result.returncode==0)==success,(label,text)
    if not success:
        assert 'reconstruct status=1/24106' in text,text
        assert not (output/'current').exists()
        return
    assert 'physical_readback=exact' in text and 'native_restart status=0/0' in text,text
    return output,json.loads((source/'native.json').read_text())

def pressure_records(root,expected):
    generation=root/(root/'current').read_text().strip()
    owned=set()
    for path in generation.iterdir():
        data=path.read_bytes()
        if data[:8]!=b'H4RANK01':continue
        offset=0
        def take(fmt):
            nonlocal offset
            result=struct.unpack_from('<'+fmt,data,offset);offset+=struct.calcsize('<'+fmt)
            return result
        magic,version,ranks,rank=take('8sIII');assert version==1
        common=take('iiiQQQdddQQI');assert common[:3]==(8,8,8)
        reference=common[8];fields=[take('BHB') for _ in range(common[-1])]
        box=take('iiiiii');begin=box[:3];shape=box[3:];count=math.prod(shape)
        scalar=0
        for role,identity,width in fields:
            size,=take('Q');assert size==width*count
            values=take('d'*size)
            if role==4:
                name=model['transported_scalars'][scalar]['stable_name'];scalar+=1
                assert width==1 and all(v==fractions[species.index(name)] for v in values)
            if role!=1:continue
            for i,value in enumerate(values):
                xyz=(begin[0]+i%shape[0],begin[1]+i//shape[0]%shape[1],begin[2]+i//(shape[0]*shape[1]))
                assert xyz not in owned;owned.add(xyz)
                index=xyz[0]+8*(xyz[1]+8*xyz[2])
                assert value+reference==expected(index),(xyz,value,reference,expected(index))
        assert scalar==6
    assert len(owned)==512

signed=lambda i: -17942736.+1000*(i%8)
positive=lambda i: 100000.+1000*(i%8)
negative=transfer('negative',signed,lambda i:p0)
plus=transfer('positive',positive,lambda i:p0)
two,audit=run(negative,'two',2)
four,again=run(negative,'four',4)
pos,other=run(plus,'positive',2)
for root,expected in ((two,signed),(four,signed),(pos,positive)):pressure_records(root,expected)
for report in (audit,again,other):
    assert report['pressure_model']=='fixed_thermodynamic'
    assert report['thermodynamic_pressure_pa']==p0
    assert report['mechanical_pressure_reference_pa']==100000.
    assert report['species_order']==species
    assert report['independent_species_order']==[s['stable_name'] for s in model['transported_scalars']]
    assert report['history']=='current_state_v1_rebuild'
    assert abs(report['native_mass_kg']-p0/(gas*temperature))<1e-11
    assert abs(report['native_energy_J']-(p0/(gas*temperature)*h-p0))<1e-8
    assert abs(report['relative_mass_change'])<1e-12
    assert abs(report['relative_energy_change'])<1e-12
for key in ('native_mass_kg','native_energy_J','reference_mass_kg','reference_energy_J'):
    assert audit[key]==other[key],(key,audit[key],other[key])
# A boundary-matched uniform field checks the first native recovery step.
equilibrium=transfer('equilibrium',lambda i:100000.,lambda i:p0)
seed,equilibrium_audit=run(equilibrium,'seed',2)
resumed=work.with_name(work.name+'-resumed')
if resumed.exists():shutil.rmtree(resumed)
with (work/'resume.log').open('w') as log:
    result=subprocess.run(list(map(str,[mpi,'--oversubscribe','--bind-to','none','-n',4,
        application,'run',work,'--output',resumed,'--steps',1,'--restart',seed,
        '--output-interval',0,'--restart-interval',1,'--diagnostics-interval',1])),stdout=log,stderr=log)
assert result.returncode==0,(work/'resume.log').read_text()[-12000:]
generation=resumed/'Restart'/(resumed/'Restart/current').read_text().strip()
manifest=(generation/'manifest.bin').read_bytes()
header=struct.unpack_from('<8sIIiiiQQQdddQQI',manifest)
assert header[12]==18 and abs(header[9]-.00101)<1e-16
# Coupled EOS continues to require a positive absolute pressure.
model['flow'].pop('thermodynamic_pressure_pa');(work/'case.json').write_text(json.dumps(model))
run(negative,'rejected',2,False)
coupled=transfer('coupled',positive,positive)
root,coupled_audit=run(coupled,'coupled',2)
pressure_records(root,positive)
assert coupled_audit['pressure_model']=='coupled_eos'
assert coupled_audit['thermodynamic_pressure_pa'] is None
assert abs(coupled_audit['relative_mass_change'])<1e-12
report=dict(passed=True,source_pressure_range=[signed(0),signed(7)],
    fixed=audit,coupled=coupled_audit,source_ranks=[2,4],pressure_readback='exact',
    native_recovery=dict(source_step=17,accepted_step=header[12],accepted_time=header[9],source_ranks=2,restart_ranks=4),
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
(work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
