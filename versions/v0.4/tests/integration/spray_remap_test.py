# SPDX-License-Identifier: Apache-2.0
"""Exercise the real native remapper and transactional inventory publication."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location('remap', sys.argv[1])
remap = importlib.util.module_from_spec(spec)
spec.loader.exec_module(remap)
fixtures, native = Path(sys.argv[3]), Path(sys.argv[2])
with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    args = SimpleNamespace(inventory=fixtures/'spray-parcels.jsonl',
        inventory_sha=remap.sha(fixtures/'spray-parcels.jsonl'),
        source_gas=fixtures/'spray-source.json', gas=fixtures/'spray-gas.yaml',
        liquid=fixtures/'spray-liquid.asset', native=native, output=root/'map',
        phase='kerosene-thermo', species='H2,H2O,CO,CO2,O2,N2,C12H23',
        injector_id=1, injector_next_ordinal=77, injector_residual_mass=1e-12)
    report = remap.convert(args)
    first = (root/'map/parcels.jsonl').read_bytes()
    other = copy.copy(args)
    other.output = root/'again'
    remap.convert(other)
    assert first == (root/'again/parcels.jsonl').read_bytes()
    source = [json.loads(s) for s in args.inventory.read_text().splitlines()]
    mapped = [json.loads(s) for s in first.splitlines()]
    assert [p['source'] for p in mapped] == source
    assert report['halo'] == 1 and report['daughters'] == 1
    restart = json.loads((root/'map/restart.json').read_text())
    assert restart['format'] == 'hundun_spray_restart_import_v1'
    assert restart['parcels'] == len(mapped)
    assert restart['injectors'] == [dict(id=1, next_ordinal=77,
                                         residual_mass_kg=1e-12)]
    for p in mapped:
        assert p['target_thermodynamics']['temperature_k'] > p['source']['temperature_k']+.1
        assert p['target_thermodynamics']['droplet_diameter_m'] > p['source']['droplet_diameter_m']*1.006
    failures = 0
    def rejects(changed):
        global failures
        try:
            remap.convert(changed)
        except (ValueError, subprocess.SubprocessError):
            failures += 1
        else:
            raise AssertionError('invalid request published')
        assert (root/'map/parcels.jsonl').read_bytes() == first
        if Path(changed.output) != args.output:
            assert not Path(changed.output).exists()
    rejects(args)  # completed inventory stays intact
    other.output = root/'bad'
    other.inventory_sha = '0'*64
    rejects(other)
    other = copy.copy(args)
    other.output = root/'bad'
    other.source_gas = args.gas
    rejects(other)
    other = copy.copy(args)
    other.output = root/'bad'
    corrupt = root/'duplicate.jsonl'
    corrupt.write_text(args.inventory.read_text()+args.inventory.read_text().splitlines()[0]+'\n')
    other.inventory, other.inventory_sha = corrupt, remap.sha(corrupt)
    rejects(other)
    invalid = copy.deepcopy(source)
    invalid[0]['droplet_diameter_m'] = 0
    corrupt.write_text(''.join(json.dumps(p)+'\n' for p in invalid))
    other.inventory_sha = remap.sha(corrupt)
    rejects(other)
    other = copy.copy(args)
    other.output = root/'bad'
    liquid = root/'bad.asset'
    liquid.write_text(args.liquid.read_text().replace('gas_sha 2','gas_sha 3'))
    other.liquid = liquid
    rejects(other)  # full native gas/material identity validation
    fnv = 14695981039346656037
    for byte in args.liquid.read_bytes():
        fnv = ((fnv ^ byte)*1099511628211) & ((1 << 64)-1)
    command = [str(native), str(args.gas), remap.sha(args.gas), args.phase,
               str(args.liquid), str(fnv), args.species]
    for text in ('', '1e-12 1e90\n', '1e-12 -1e6 trailing\n', '0 -1e6\n'):
        r = subprocess.run(command,input=text,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                           universal_newlines=True,timeout=30)
        assert r.returncode != 0 and r.stdout == ''
    assert not list(root.glob('.remap-*'))
    print('real native conversion, byte repeat, full source retention; {} publication failures and 4 native failures passed'.format(failures))
