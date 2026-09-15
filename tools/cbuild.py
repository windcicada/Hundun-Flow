#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Build an unmodified COAST source reference in a short external directory.

Source files are symlinked read-only by convention; compiler outputs and modules
belong to the external directory. The original Makefile and serial module order
are retained. Run tools/cref.py --build <directory> for the periodic probe link.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coast', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    root, out = args.coast.resolve(), args.output.resolve()
    if out == root or root in out.parents:
        parser.error('choose a build directory outside the source checkout')
    src = root / 'SRC.Coast'
    selected = {}
    suffixes = {'.F90', '.F', '.f', '.cpp', '.hpp', '.h', '.inc'}
    for directory, dirs, files in os.walk(str(src)):
        dirs[:] = [d for d in dirs if d not in
                   {'x86_64', 'tests', 'ai', 'nanobot', '.codegraphf', '__pycache__'}]
        for name in files:
            p = Path(directory) / name
            if p.suffix in suffixes or p == src / 'Makefile':
                selected[str(p.relative_to(src))] = hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = out / 'build.json'
    if out.exists():
        if not manifest.is_file():
            parser.error('existing output requires its source build manifest')
        old = json.loads(manifest.read_text())
        if old.get('source') != str(src) or old.get('files') != selected:
            parser.error('existing build requires the same source identity')
    out.mkdir(parents=True, exist_ok=True)
    for name in selected:
        target = out / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.symlink_to(src / name)
        elif not target.is_symlink() or target.resolve() != (src / name).resolve():
            parser.error('source link identity mismatch: '+name)
    command = ['make', '-j1', 'check', 'coast', 'EXEDIR='+str(out),
               'OCE_ROOT='+str(root / 'third_party/oce_sysroot/usr')]
    report = {'source': str(src), 'files': selected, 'command': command,
              'status': 'prepared',
              'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                      cwd=str(root), universal_newlines=True).strip(),
              'compiler': subprocess.check_output(['mpif90', '--version'],
                          universal_newlines=True).splitlines()[0]}
    manifest.write_text(json.dumps(report, indent=2)+'\n')
    if args.prepare_only:
        print('prepared {} source links in {}'.format(len(selected), out))
        return
    with (out / 'build.log').open('w') as log:
        status = subprocess.call(command, cwd=str(out), stdout=log, stderr=subprocess.STDOUT)
    report['status'] = 'built' if status == 0 else 'build failure'
    if status == 0:
        report['executable_sha256'] = hashlib.sha256((out / 'coast').read_bytes()).hexdigest()
    changed = [name for name, expected in selected.items()
               if hashlib.sha256((src / name).read_bytes()).hexdigest() != expected]
    report['source_changes'] = changed
    manifest.write_text(json.dumps(report, indent=2)+'\n')
    print(report['status']+': '+str(out / 'build.log'))
    if status or changed:
        raise SystemExit(status or 1)


if __name__ == '__main__':
    main()
