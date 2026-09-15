#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Audit native SGS output for the constant-cp Smagorinsky/IBM spray fixture."""
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET


def product(values):
    result = 1
    for value in values:
        result *= value
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(run, case):
    model = json.loads((case / 'case.json').read_text())
    assert model['turbulence']['model'] == 'smagorinsky'
    assert model['mesh']['kind'] == 'uniform'
    cs = model['turbulence']['coefficient']
    shape = model['mesh']['exact_cells']
    lower, upper = (model['mesh']['domain'][k] for k in ('lower', 'upper'))
    spacing = [(hi-lo)/n for lo, hi, n in zip(lower, upper, shape)]
    delta = product(spacing) ** (1/3)
    assert cs > 0
    # This test intentionally uses the fixture's independent constant-cp EOS
    # and Sutherland law. Keep that physical reference explicit.
    thermal = (case / 'thermophysics.d').read_text()
    assert thermal.count('molecular_weight 28\n') == 2
    assert thermal.count('transport_sutherland 1e-5 300 0 0.7\n') == 2
    assert thermal.count('nasa7_low 3.3676259411963292') == 2
    stl = case / model['mesh']['immersed_boundary']['stl_file']
    vertices = [tuple(map(float, line.split()[1:])) for line in stl.read_text().splitlines()
                if line.strip().startswith('vertex ')]
    box = [(min(v[a] for v in vertices), max(v[a] for v in vertices)) for a in range(3)]
    assert len(set(vertices)) == 8
    names = ('nu_sgs', 'k_sgs', 'eps_sgs_volume', 'eps_sgs_specific')
    errors = {'kinetic_energy': 0., 'dissipation_specific': 0., 'dissipation_volume': 0.}
    coverage, files, solid_cells, active_sgs = {}, {}, 0, 0
    for segment, ranks in (('fresh', 1), ('resumed', 4)):
        root = run / segment
        restart = root / 'Restart'
        generation = restart / (restart / 'current').read_text().strip()
        header = struct.unpack_from('<8sIIiiiQQQdddQQI', (generation / 'manifest.bin').read_bytes())
        assert header[0] == b'H4MANI01' and header[1] == 5 and header[2] == ranks
        final_step, pressure_reference = header[12], header[11]
        for path in sorted((root / 'Visit').glob('*.vti')):
            step = int(path.name.split('-')[1])
            key = '{}:{}'.format(segment, step)
            occupied = coverage.setdefault(key, set())
            data = path.read_bytes()
            split = data.index(b'<AppendedData')
            tree = ET.fromstring(data[:split] + b'</VTKFile>')
            piece = tree.find('./ImageData/Piece')
            extent = list(map(int, piece.attrib['Extent'].split()))
            begin = extent[::2]
            counts = [extent[a+1]-extent[a] for a in (0, 2, 4)]
            n = product(counts)
            base = data.index(b'_', data.index(b'>', split)) + 1
            fields = {}
            for node in piece.find('CellData'):
                offset = base + int(node.attrib['offset'])
                length, = struct.unpack_from('<Q', data, offset)
                components = int(node.attrib['NumberOfComponents'])
                assert node.attrib['type'] == 'Float64' and length == n*components*8
                values = struct.unpack_from('<{}d'.format(n*components), data, offset+8)
                assert all(math.isfinite(value) for value in values)
                fields[node.attrib['Name']] = values
            assert set(fields) == set(names) | {'U', 'pi', 'h', 'A'}
            for cell in range(n):
                local = (cell % counts[0], (cell // counts[0]) % counts[1], cell // (counts[0]*counts[1]))
                global_cell = tuple(b+i for b, i in zip(begin, local))
                assert global_cell not in occupied
                occupied.add(global_cell)
                center = [lo+(i+.5)*dx for lo, i, dx in zip(lower, global_cell, spacing)]
                solid = all(lo < x < hi for x, (lo, hi) in zip(center, box))
                values = [fields[name][cell] for name in names]
                assert all(v >= 0 for v in values)
                if solid:
                    assert values == [0., 0., 0., 0.]
                    solid_cells += 1
                    continue
                nu, k, volume, specific = values
                active_sgs += nu > 0
                expected_k = (nu / delta)**2 / cs**(8/3)
                errors['kinetic_energy'] = max(errors['kinetic_energy'], abs(k-expected_k)/max(1e-300, k, expected_k))
                if step == final_step:
                    temperature = 298.15+(fields['h'][cell]-1e5*fields['A'][cell])/1000
                    rho = (pressure_reference+fields['pi'][cell])*28/(8314.46261815324*temperature)
                    mu = 1e-5*math.sqrt(temperature/300)
                    strain = nu/(cs*delta)**2
                    expected_specific = (nu+mu/rho)*strain**2
                    for name, actual, expected in (
                            ('dissipation_specific', specific, expected_specific),
                            ('dissipation_volume', volume, rho*expected_specific)):
                        errors[name] = max(errors[name], abs(actual-expected)/max(1e-300, actual, expected))
            files[str(path.relative_to(run))] = sha(path)
    assert len(coverage) == 4 and len(files) == 10
    assert all(len(cells) == product(shape) for cells in coverage.values())
    assert solid_cells > 0 and active_sgs > 0
    assert max(errors.values()) <= 1e-11, errors
    report = {'schema': 'hundun_sgs_output_audit_v1', 'passed': True,
              'scope': 'native Smagorinsky IBM spray ESF fixture; fresh and 1-to-4 restart output',
              'frames': len(coverage), 'files': len(files), 'solid_samples': solid_cells,
              'active_sgs_samples': active_sgs, 'relative_errors': errors,
              'artifact_sha256': files, 'case_sha256': sha(case / 'case.json')}
    (run / 'sgs.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'artifact_sha256'}, indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
