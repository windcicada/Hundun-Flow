#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""624CF approved liquid transfer: source mass, momentum and reconstructed h.

The frozen source defines h_liquid(T)=h_NASA,vapor(T)-L(T), with its own R.
The native executable inverts the admitted target material. Output is a
thermodynamic transfer inventory; target ownership and SGS history admission
are separate restart operations. Source records remain intact.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

SOURCE_GAS_SHA = '0b8ff99b7c8c7937b39019d946f543210e694f2241d330d893cf147851b1c604'
SOURCE_PROPERTIES_SHA = 'dfcd4fb17f03fceb6ec2208a948c7b2d613506868ca9f02adf5ea5bdc27be383'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def source_state(p, gas):
    """The checkpoint carries T; absolute h is the approved reconstruction."""
    t, d, n = p['temperature_k'], p['droplet_diameter_m'], p['multiplicity']
    require(all(math.isfinite(v) and v > 0 for v in (t, d, n)), 'positive finite primitives')
    require(len(p['position_m']) == len(p['velocity_m_per_s']) == 3, '3D parcel')
    require(all(math.isfinite(v) for v in p['position_m']+p['velocity_m_per_s']+p['legacy_state']),
            'finite vectors and history')
    require(p['kind'] in ('retained', 'daughter'), 'source kind')
    require(len(p['legacy_state']) == (22 if p['kind'] == 'retained' else 10), 'source history size')
    require(p['source_rank'] >= 0 and p['source_ordinal'] > 0, 'source identity')
    require(len(p['source_sha256']) == 64 and
            all(c in '0123456789abcdef' for c in p['source_sha256']), 'source SHA256')
    fuel = next(x for x in gas['species'] if x['name'] == 'C12H23')
    require(fuel['meta'][1] <= t < 684.26, 'source property domain')
    a = fuel['low'] if t <= fuel['meta'][3] else fuel['high']
    vapor = gas['reference_constants'][0]/fuel['meta'][0]*(
        sum(a[i]*t**(i+1)/(i+1) for i in range(5))+a[5])
    h = vapor-250183*((684.26-t)/(684.26-483.15))**.38
    # Literal old Fortran 9.255437+03: keep original mass before correction.
    rho = 1037.096-.7233865*t-9.255437+3/(733-t)
    mass = rho*math.pi/6*d**3
    require(math.isfinite(mass) and mass > 0 and math.isfinite(mass*n) and mass*n > 0
            and math.isfinite(h), 'source mass and enthalpy')
    return mass, h


def convert(args):
    output = Path(args.output)
    require(not output.exists(), 'output directory already exists')
    inputs = {name: Path(getattr(args, name)).resolve()
              for name in ('inventory', 'source_gas', 'gas', 'liquid', 'native')}
    hashes = {name: sha(path) for name, path in inputs.items()}
    require(hashes['inventory'] == args.inventory_sha, 'inventory SHA256 mismatch')
    require(hashes['source_gas'] == SOURCE_GAS_SHA, 'frozen 624CF thermodynamics identity')
    gas = json.loads(inputs['source_gas'].read_text())
    records = [json.loads(line) for line in inputs['inventory'].read_text().splitlines()]
    require(records, 'empty inventory')
    keys = [(p['source_sha256'], p['source_rank'], p['kind'], p['source_ordinal']) for p in records]
    require(len(set(keys)) == len(keys), 'duplicate source record identity')
    invariants = [source_state(p, gas) for p in records]
    numeric = ''.join('{:.17g} {:.17g}\n'.format(m, h) for m, h in invariants)
    fnv = 14695981039346656037
    for byte in inputs['liquid'].read_bytes():
        fnv = ((fnv ^ byte)*1099511628211) & ((1 << 64)-1)
    command = [str(inputs['native']), str(inputs['gas']), hashes['gas'], args.phase,
               str(inputs['liquid']), str(fnv), args.species]
    result = subprocess.run(command, input=numeric, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True, timeout=180)
    require(result.returncode == 0, 'native conversion failed: '+result.stderr[-2000:])
    lines = result.stdout.splitlines()
    require(len(lines) == len(records), 'native row coverage')
    mapped, source_m, target_m, source_h, target_h, kinetics = [], [], [], [], [], []
    shifts, ratios, h_errors, m_errors, target_kinetics = [], [], [], [], []
    source_p, target_p = [[], [], []], [[], [], []]
    material = None
    for p, (mass, h), line in zip(records, invariants, lines):
        fields = line.split()
        require(len(fields) == 7, 'native row width')
        t, d, rho, hn, residual, geometry = map(float, fields[:6])
        fingerprint = int(fields[6])
        require(all(math.isfinite(v) for v in (t, d, rho, hn, residual, geometry)), 'finite result')
        require(t > 0 and d > 0 and rho > 0 and 0 < fingerprint < 2**64, 'native state')
        require(material is None or material == fingerprint, 'single material identity')
        material = fingerprint
        reconstructed_mass = rho*math.pi/6*d**3
        rel_m = abs(reconstructed_mass-mass)/mass
        rel_h = abs(hn-h)/max(1., abs(h), abs(hn))
        require(rel_m < 1e-11 and rel_h < 1e-11, 'per-parcel conservation gate')
        require(abs(residual-(hn-h)) <= 1e-12*max(1.,abs(h)), 'enthalpy audit consistency')
        require(abs(geometry) < 1e-11, 'native geometry gate')
        weight, velocity = p['multiplicity'], p['velocity_m_per_s']
        source_m.append(mass*weight)
        target_m.append(reconstructed_mass*weight)
        source_h.append(mass*weight*h)
        target_h.append(mass*weight*hn)
        kinetics.append(mass*weight*.5*sum(v*v for v in velocity))
        target_kinetics.append(reconstructed_mass*weight*.5*sum(v*v for v in velocity))
        for axis in range(3):
            source_p[axis].append(mass*weight*velocity[axis])
            target_p[axis].append(reconstructed_mass*weight*velocity[axis])
        shifts.append(t-p['temperature_k'])
        ratios.append(d/p['droplet_diameter_m'])
        h_errors.append(abs(hn-h))
        m_errors.append(rel_m)
        mapped.append(dict(source=p, droplet_mass_kg=mass, source_specific_enthalpy_j_per_kg=h,
                           target_thermodynamics=dict(temperature_k=t, droplet_diameter_m=d,
                             density_kg_per_m3=rho, specific_enthalpy_j_per_kg=hn,
                             liquid_material_fingerprint=fingerprint)))
    total_m = math.fsum(source_m)
    enthalpy_error = abs(math.fsum(target_h)-math.fsum(source_h))/math.fsum(abs(x) for x in source_h)
    require(enthalpy_error < 1e-11, 'global enthalpy gate')
    mass_error = abs(math.fsum(target_m)-total_m)/total_m
    momentum_errors = [abs(math.fsum(b)-math.fsum(a))/max(1e-300,math.fsum(abs(x) for x in a))
                       for a,b in zip(source_p,target_p)]
    kinetic_error = abs(math.fsum(target_kinetics)-math.fsum(kinetics))/max(1e-300,math.fsum(kinetics))
    require(max([mass_error,kinetic_error]+momentum_errors) < 1e-11, 'global mechanical gate')
    for name, path in inputs.items():
        require(sha(path) == hashes[name], 'input changed during conversion: '+name)
    report = dict(format='hundun_liquid_transfer_v1', scope='thermodynamic transfer inventory',
        policy='preserve source mass, weight, velocity and reconstructed absolute enthalpy',
        source_enthalpy='NASA vapor h(T; source R and MW) minus latent heat at checkpoint T',
        target_enthalpy='native liquid asset reference plus integrated cp',
        source_density='1037.096 - 0.7233865*T - 9.255437 + 3/(733-T)',
        source_properties_sha256=SOURCE_PROPERTIES_SHA,
        pending=['target mesh ownership', 'legacy SGS history physical admission', 'native Restart assembly'],
        inputs={name: dict(path=str(path), sha256=hashes[name]) for name, path in inputs.items()},
        parcels=len(records), retained=sum(p['kind']=='retained' for p in records),
        daughters=sum(p['kind']=='daughter' for p in records),
        halo=sum(p['source_bin']=='halo' for p in records),
        represented_particles=math.fsum(p['multiplicity'] for p in records),
        source_liquid_mass_kg=total_m, target_geometry_mass_kg=math.fsum(target_m),
        source_momentum_kg_m_s=[math.fsum(x) for x in source_p],
        target_geometry_momentum_kg_m_s=[math.fsum(x) for x in target_p],
        source_enthalpy_j=math.fsum(source_h), target_enthalpy_j=math.fsum(target_h),
        kinetic_energy_j=math.fsum(kinetics), target_geometry_kinetic_energy_j=math.fsum(target_kinetics),
        relative_global_mass_error=mass_error, relative_global_momentum_errors=momentum_errors,
        relative_global_kinetic_error=kinetic_error, relative_global_enthalpy_error=enthalpy_error,
        temperature_shift_range_k=[min(shifts),max(shifts)], diameter_ratio_range=[min(ratios),max(ratios)],
        max_absolute_specific_enthalpy_error_j_kg=max(h_errors), max_relative_geometry_mass_error=max(m_errors))
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.remap-', dir=str(output.parent)))
    try:
        parcel_path = stage/'parcels.jsonl'
        with parcel_path.open('w') as stream:
            for p in mapped:
                stream.write(json.dumps(p, sort_keys=True, allow_nan=False, separators=(',', ':'))+'\n')
        report['mapped_sha256'] = sha(parcel_path)
        (stage/'report.json').write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
        # All source, native, per-parcel and global checks precede publication.
        require(not output.exists(), 'output directory appeared during conversion')
        os.rename(str(stage), str(output))
    finally:
        if stage.exists():
            shutil.rmtree(str(stage))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('inventory', 'inventory-sha', 'source-gas', 'gas', 'liquid', 'native', 'output'):
        p.add_argument('--'+name, required=True)
    p.add_argument('--phase', default='kerosene-thermo')
    p.add_argument('--species', default='H2,H2O,CO,CO2,O2,N2,C12H23')
    args = p.parse_args()
    result = convert(args)
    print(json.dumps({k: result[k] for k in ('parcels', 'source_liquid_mass_kg',
                    'source_enthalpy_j', 'temperature_shift_range_k', 'diameter_ratio_range',
                    'max_relative_geometry_mass_error', 'relative_global_enthalpy_error')}, indent=2))


if __name__ == '__main__':
    main()
