#!/usr/bin/env python3
"""Internal release-content gate: paths and bytes, including binary strings."""
import argparse
from pathlib import Path
import sys

FORBIDDEN = (b'coast', b'codex')


def scan(root):
    failures = []
    root = Path(root)
    for path in sorted(root.rglob('*')):
        relative = str(path.relative_to(root))
        if any(word in relative.lower().encode('utf-8') for word in FORBIDDEN):
            failures.append((relative, 'path'))
        if path.is_symlink():
            import os
            raw = os.readlink(str(path)).lower().encode('utf-8')
            if any(word in raw for word in FORBIDDEN):
                failures.append((relative, 'link'))
        elif path.is_file():
            # Chunk overlap also checks strings crossing a read boundary.
            tail = b''
            with path.open('rb') as stream:
                while True:
                    block = stream.read(1024*1024)
                    if not block:
                        break
                    raw = tail + block.lower()
                    if any(word in raw for word in FORBIDDEN):
                        failures.append((relative, 'contents'))
                        break
                    tail = raw[-4:]
    return failures


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    if not args.directory.is_dir():
        parser.error('expected an extracted source or runtime directory')
    failures = scan(args.directory)
    for path, kind in failures:
        print('{}: {}'.format(kind, path))
    print('release content scan: {} issue(s)'.format(len(failures)))
    sys.exit(bool(failures))
