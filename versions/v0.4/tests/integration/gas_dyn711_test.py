#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Native dyn711: partial-window Restart, terminal statistics and physical ledgers."""
import hashlib
import math
import json
from timing_check import check_timing
from pathlib import Path
import shutil
import subprocess
import sys

binary, fixture, mpi, work, validator, compare = map(Path, sys.argv[1:7])
fields = int(sys.argv[7])
cold = len(sys.argv)>8 and sys.argv[8] in ("cold","cold300","quiet300")
uniform=cold and sys.argv[8]=="quiet300"
temperature=300 if cold and sys.argv[8]!="cold" else 1100
if work.exists():
    shutil.rmtree(work)
shutil.copytree(fixture, work)
model = json.loads((work/'case.json').read_text())
if cold:
    # An O2/N2 stream has exactly zero chemical heat. Its nonzero
    # oxygen variance exercises the source's zero-heat remix branch directly.
    model['boundaries']['x_min'].update(temperature=temperature,flow_kind='velocity_inlet',
        velocity=[1,0,0],mass_flow_rate=0)
    model['boundaries']['x_max']['backflow_temperature']=temperature
    for side in ('x_min','x_max'):
        for scalar in model['boundaries'][side]['scalars']:
            value=(.23 if side=='x_min' and not uniform else .232) if scalar['stable_name']=='O2' else 0.
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
    boundary_z=(2*(oxygen-.23)/31.998)/(4/16.043+2*oxygen/31.998) if cold and not uniform and side=='x_min' else z
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
roundoff_applications=0

def run(label, ranks, steps, restart=None):
    global max_mass, max_energy, max_element, roundoff_applications
    output = work.with_name(work.name+label)
    if output.exists():
        shutil.rmtree(output)
    args = [mpi, '--oversubscribe', '--bind-to', 'none', '-n', ranks, binary,
        'run', work, '--output', output, '--steps', steps,
        '--output-interval', 0, '--restart-interval', steps, '--diagnostics-interval', 1]
    if restart:
        args += ['--restart', restart/'Restart']
    else:
        initial=f'100000,{temperature},1,0,0,0,0.232,0,0,0,0,0' if cold else \
            '100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001,'+str(z)
        args += ['--initial-state', initial]
    call(args, work/(label+'.log'))
    args = [sys.executable, validator, 'runtime', output/'evidence.jsonl']
    if restart:
        generation = (restart/'Restart/current').read_text().strip()
        args += ['--run-start-manifest', restart/'Restart'/generation/'manifest.bin']
    call(args, work/(label+'-audit.log'))
    check_timing(output/'monitor.jsonl', {'esf_reaction', 'tcr_statistics'})
    evidence=[json.loads(line) for line in (output/'evidence.jsonl').read_text().splitlines()]
    for row in evidence:
        flow=row['cold']
        assert flow['pdf_before_flow'] is True
        assert flow['outer_iterations']==(model['solver'].get('reference_outer_iterations',0) or 2)
        assert flow['species_solve_calls']==0 and flow['species_endpoint_solve_calls']==0
        assert flow['species_iterations']==0
        assert flow['enthalpy_solve_calls']==0 and flow['enthalpy_iterations']==0
        assert flow['momentum_iterations']<=3*50*flow['outer_iterations']
        assert flow['pressure_iterations']<=500*flow['outer_iterations']
        assert flow.get('passive_iterations',0)<=20*flow.get('passive_solve_calls',0)

    assert 'cold_esf_flux' not in (work/(label+'.log')).read_text()
    assert 'esf_transport_bound' not in (work/(label+'.log')).read_text()
    # The observed-balance contract still rejects fabricated solve counts.
    forged=json.loads(json.dumps(evidence))
    forged[-1]['cold']['species_solve_calls']=1
    bad=work/(label+'-schedule-bad.jsonl')
    bad.write_text(''.join(json.dumps(row)+'\n' for row in forged))
    denied=subprocess.run([sys.executable,str(validator),'runtime',str(bad)],
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    assert denied.returncode!=0,'fabricated PDF flow solve was accepted'


    for row in map(json.loads, (output/'diagnostics.jsonl').read_text().splitlines()):
        p = row['payload']
        mass = abs(p['mass_balance_defect_kg_s']*p['dt'])/p['mass_kg']
        energy = abs(p['total_energy_balance_defect_W']*p['dt'])/max(1.,abs(p['internal_energy_J'])+p['kinetic_energy_J'])
        element = max(v['relative_defect'] for v in p['composition_balance']['elements'])
        assert mass < 1e-6 and math.isfinite(energy), (mass,energy)
        budget=p['composition_balance']
        assert budget['roundoff_rule']=='fp64-local-storage-v1'
        for group in ('species','elements'):
            for entry in budget[group]:
                if entry['relative_defect']>=1e-6:
                    assert math.isfinite(entry['defect']) and math.isfinite(entry['relative_defect']),entry
                    if entry['roundoff_applied']:
                        assert abs(entry['defect'])<=entry['storage_roundoff_bound'],entry
                        roundoff_applications+=1
                else:
                    assert not entry['roundoff_applied'],entry
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
else:
    args += ['--reference-linear']
call(args, work/'compare.log')
comparison = (work/'compare.log').read_text()
assert 'dyn711_history clocks=restored' in comparison and 'passed=1' in comparison
assert not uniform or roundoff_applications>0
current_history=None
if not cold and len(sys.argv)>8:
    migrated=work.with_name(work.name+'m')
    if migrated.exists():shutil.rmtree(migrated)
    call([mpi,'--oversubscribe','--bind-to','none','-n',2,sys.argv[8],
          work,seed/'Restart',migrated/'Restart'],work/'current-history.log')
    assert 'model_records=exact' in (work/'current-history.log').read_text()
    resumed=run('mr',4,1,migrated)
    evidence=[json.loads(line) for line in (resumed/'evidence.jsonl').read_text().splitlines()]
    for row in evidence:
        if 'run_start' in row:
            row['run_start']['history']['source_signature']=1
    tampered=work/'current-bad.jsonl'
    tampered.write_text(''.join(json.dumps(row)+'\n' for row in evidence))
    denied=subprocess.run([sys.executable,str(validator),'runtime',str(tampered)],
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    assert denied.returncode!=0,'V6 fabricated time signature was accepted'
    current_history=dict(version=6,source_ranks=2,restart_ranks=4,model_records='exact',
                         time_history='rebuild',resumed=str(resumed))
report = dict(fields=fields, steps=1 if cold else 9, source_ranks=2,
    restart_ranks=None if cold else 4, reader_ranks=4, cold_remix=cold,
    temperature_K=temperature,
    uniform=uniform, roundoff_applications=roundoff_applications,
    mass_defect=max_mass, energy_defect=max_energy, element_defect=max_element,
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    comparison=comparison, current_history=current_history, passed=True)
(work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print('dyn711 fields=%d cold_remix=%s mass=%.5g energy=%.5g elements=%.5g passed=1' %
      (fields,cold,max_mass,max_energy,max_element))
