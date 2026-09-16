#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Freeze a read-only original gas source snapshot and its call-site evidence."""
import argparse
import hashlib
import json
import re
from pathlib import Path


def inventory(root):
    paths = sorted(p for p in root.rglob('*') if p.is_file() and
                   (p.suffix.lower() in ('.f', '.f90', '.inc', '.h') or
                    p.name in ('Makefile', 'VERSION')))
    records = [{'path': str(p.relative_to(root)), 'bytes': p.stat().st_size,
                'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]
    identity = hashlib.sha256(json.dumps(records, sort_keys=True,
                                         separators=(',', ':')).encode()).hexdigest()
    calls = []
    for name in ('boffin.F90', 'fieldpdf.F90', 'mixer.F90', 'reactor.F90',
                 'statistics.F90', 'compress.F90', 'update.F90'):
        for n, line in enumerate((root/name).read_text(errors='replace').splitlines(), 1):
            code = line.split('!')[0].strip()
            if re.search(r'\bcall\s+|\bdo\s+(iter|jstep)\s*=', code, re.I):
                calls.append({'file': name, 'line': n, 'statement': code})
    return dict(schema='hundun.gas.reference.v1', source_root=str(root),
                source_sha256=identity, files=records, call_sites=calls,
                build=dict(declared_compiler='mpif90', declared_flags='-O3 -IINCLIB',
                           declared_real='default REAL / MPI_REAL',
                           executable_provenance='pending matched build and runtime receipt'),
                scope='source identity; runtime model and compiler precision require case receipts')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    record = inventory(args.source.resolve())
    if args.verify:
        old = json.loads(args.output.read_text())
        if record != old:
            raise SystemExit('reference changed')
    else:
        args.output.write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n')
    print(record['source_sha256'], len(record['files']), 'files')


if __name__ == '__main__':
    main()
