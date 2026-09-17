#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Native dyn711: partial-window Restart, terminal statistics and physical ledgers."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

binary, fixture, mpi, work, validator, compare = map(Path, sys.argv[1:7])
fields = int(sys.argv[7])
cold = len(sys.argv)>8 and sys.argv[8]=="cold"
if work.exists():
    shutil.rmtree(work)
shutil.copytree(fixture, work)
model = json.loads((work/'case.json').read_text())
if cold:
    # An O2/N2 stream has exactly zero chemical heat. Its nonzero
    # oxygen variance exercises the source's zero-heat remix branch directly.
    model['boundaries']['x_min'].update(temperature=1100,flow_kind='velocity_inlet',
        velocity=[1,0,0],mass_flow_rate=0)
    model['boundaries']['x_max']['backflow_temperature']=1100
    for side in ('x_min','x_max'):
        for scalar in model['boundaries'][side]['scalars']:
            value=(.23 if side=='x_min' else .232) if scalar['stable_name']=='O2' else 0.
            scalar['value']=value if side=='x_min' else 0.
            scalar['backflow_value']=value
ensemble = model['reaction']['ensemble']
ensemble['fields'] = fields
ensemble['initial_species_offsets'] = [
    (.0001 if f % 2 == 0 else -.0001) if q==(1 if cold else 0) else 0.
    for f in range(fields) for q in range(6)]
oxygen = .232
z = (4*.05/16.043 - 2*.2/31.998 + .005/28.01 + .001/2.016 + 2*oxygen/31.998)/(4/16.043 + 2*oxygen/31.998)
if cold:
    z=0.
ensemble['tcr'] = dict(model='dyn711_v1', mode='experimental', fuel='CH4',
    weak_rate_threshold=1e-30, mixture_fraction='Z', oxidizer_oxygen_mass_fraction=oxygen)
model['transported_scalars'].append(dict(stable_name='Z', role='passive_scalar',
    molecular_schmidt=.7, turbulent_schmidt=.7))
for side in ('x_min', 'x_max'):
    boundary_z=(2*(oxygen-.23)/31.998)/(4/16.043+2*oxygen/31.998) if cold and side=='x_min' else z
    model['boundaries'][side]['scalars'].append(dict(stable_name='Z',
        kind='dirichlet' if side == 'x_min' else 'zero_gradient',
        value=boundary_z if side == 'x_min' else 0, backflow_kind='dirichlet', backflow_value=z))
(work/'case.json').write_text(json.dumps(model, indent=2)+'\n')

def call(args, path):
    with path.open('w') as log:
        result = subprocess.run(list(map(str,args)), stdout=log, stderr=log)
    assert result.returncode == 0, path.read_text()[-16000:]

call([binary, 'check', work], work/'check.log')
check = (work/'check.log').read_text()
assert 'tcr_model=dyn711_v1 tcr_mode=experimental' in check and 'passive_scheme=CN' in check
max_mass = max_energy = max_element = 0.

def run(label, ranks, steps, restart=None):
    global max_mass, max_energy, max_element
    output = work.with_name(work.name+label)
    if output.exists():
        shutil.rmtree(output)
    args = [mpi, '--oversubscribe', '--bind-to', 'none', '-n', ranks, binary,
        'run', work, '--output', output, '--steps', steps,
        '--output-interval', 0, '--restart-interval', steps, '--diagnostics-interval', 1]
    if restart:
        args += ['--restart', restart/'Restart']
    else:
        initial='100000,1100,1,0,0,0,0.232,0,0,0,0,0' if cold else \
            '100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001,'+str(z)
        args += ['--initial-state', initial]
    call(args, work/(label+'.log'))
    args = [sys.executable, validator, 'runtime', output/'evidence.jsonl']
    if restart:
        generation = (restart/'Restart/current').read_text().strip()
        args += ['--run-start-manifest', restart/'Restart'/generation/'manifest.bin']
    call(args, work/(label+'-audit.log'))
    for row in map(json.loads, (output/'diagnostics.jsonl').read_text().splitlines()):
        p = row['payload']
        mass = abs(p['mass_balance_defect_kg_s']*p['dt'])/p['mass_kg']
        energy = abs(p['total_energy_balance_defect_W']*p['dt'])/max(1.,abs(p['internal_energy_J'])+p['kinetic_energy_J'])
        element = max(v['relative_defect'] for v in p['composition_balance']['elements'])
        assert mass < 1e-12 and energy < 1e-12 and element < 1e-6, (mass,energy,element)
        max_mass=max(max_mass,mass);max_energy=max(max_energy,energy);max_element=max(max_element,element)
    return output

# The restart carries eight accepted rate intervals into the ninth evaluation.
if cold:
    two=run('s',2,1)
    four=two
else:
    seed = run('s',2,8)
    two = run('r',2,1,seed)
    four = run('x',4,1,seed)
args=[mpi, '--oversubscribe', '--bind-to', 'none', '-n', 4, compare,
    work, two/'Restart', four/'Restart']
if cold:
    args += ['--collapsed']
call(args, work/'compare.log')
comparison = (work/'compare.log').read_text()
assert 'dyn711_history clocks=restored' in comparison and 'passed=1' in comparison
report = dict(fields=fields, steps=1 if cold else 9, source_ranks=2,
    restart_ranks=None if cold else 4, reader_ranks=4, cold_remix=cold,
    mass_defect=max_mass, energy_defect=max_energy, element_defect=max_element,
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    comparison=comparison, passed=True)
(work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print('dyn711 fields=%d cold_remix=%s mass=%.5g energy=%.5g elements=%.5g passed=1' %
      (fields,cold,max_mass,max_energy,max_element))
