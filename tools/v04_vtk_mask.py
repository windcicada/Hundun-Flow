#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
"""Extract static owned-cell masks from two COAST legacy VTK snapshots."""
from pathlib import Path
import hashlib, json, io, tarfile, sys, argparse, math, os, tempfile
from array import array

def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(p, expected_dims):
    before = p.stat()
    result = {}
    fields = []
    with p.open('rb') as f:
        require(f.readline().startswith(b'# vtk DataFile Version'), 'invalid VTK input')
        result['title'] = f.readline().decode('ascii').strip()
        require(f.readline().strip() == b'BINARY', 'invalid VTK input')
        require(f.readline().strip() == b'DATASET STRUCTURED_GRID', 'invalid VTK input')
        parts = f.readline().split()
        require(parts[0] == b'DIMENSIONS', 'invalid VTK input')
        dims = tuple(map(int, parts[1:]))
        require(dims == tuple(expected_dims), 'invalid VTK input')
        (nx, ny, nz) = dims
        require(min(dims) >= 3, 'invalid VTK input')
        parts = f.readline().split()
        require(parts == [b'POINTS', str(nx * ny * nz).encode(), b'float'], 'invalid VTK input')
        n = int(parts[1])
        raw = f.read(n * 12)
        require(len(raw) == n * 12, 'invalid VTK input')
        xyz = bytearray()
        for k in range(1, nz - 1):
            for j in range(1, ny - 1):
                start = ((k * ny + j) * nx + 1) * 12
                xyz.extend(raw[start:start + (nx - 2) * 12])
        native = array('f')
        native.frombytes(xyz)
        require(sys.byteorder == 'little', 'invalid VTK input')
        native.byteswap()
        require(all((math.isfinite(v) for v in native)), 'nonfinite owned coordinates')
        result['owned_coordinates_sha256'] = hashlib.sha256(native.tobytes()).hexdigest()
        result['points_payload_sha256'] = hashlib.sha256(raw).hexdigest()
        marker = None
        while True:
            line = f.readline()
            if not line:
                break
            a = line.split()
            if not a:
                continue
            if a[0] == b'POINT_DATA':
                require(a == [b'POINT_DATA', str(n).encode()], 'invalid VTK input')
                continue
            if a[0] == b'VECTORS':
                count = 3
            elif a[0] == b'SCALARS':
                count = int(a[3]) if len(a) > 3 else 1
                require(f.readline().split() == [b'LOOKUP_TABLE', b'default'], 'invalid VTK input')
            else:
                raise ValueError('unexpected VTK record: ' + repr(line))
            require(a[2] in (b'float', b'unsigned_char'), a)
            fields.append(a[1].decode('ascii'))
            size = n * count * (4 if a[2] == b'float' else 1)
            if a[1] == b'07_IBM_cell_type':
                require(count == 1 and marker is None, 'invalid VTK input')
                raw = f.read(size)
                require(len(raw) == size, 'invalid VTK input')
                values = array('f')
                values.frombytes(raw)
                values.byteswap()
                result['marker_fractional_storage_values'] = sum((v not in (0.0, 1.0) for v in values))
                parts = []
                for k in range(1, nz - 1):
                    for j in range(1, ny - 1):
                        start = (k * ny + j) * nx + 1
                        row = values[start:start + nx - 2]
                        require(set(row) <= {0.0, 1.0}, ('fractional owned marker', k, j, sorted(set(row))[:25]))
                        parts.append(bytes((int(v) for v in row)))
                marker = b''.join(parts)
                result['marker_payload_sha256'] = hashlib.sha256(raw).hexdigest()
            else:
                require(f.tell() + size <= before.st_size, 'truncated scalar payload')
                f.seek(size, 1)
    after = p.stat()
    require((before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino), 'invalid VTK input')
    require(marker is not None, 'invalid VTK input')
    result.update(file=str(p), bytes=before.st_size, mtime_ns=before.st_mtime_ns, dims=dims, fields=fields, owned_cells=len(marker), fluid_cells=sum(marker), owned_marker_sha256=hashlib.sha256(marker).hexdigest())
    return (result, marker)

def extract(root, steps, ranks, dims, stream):
    with tarfile.open(fileobj=stream, mode='w|gz') as archive:

        def add(name, data):
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            entry.mtime = 0
            archive.addfile(entry, io.BytesIO(data))
        rows = []
        for rank in range(ranks):
            (a, mask) = read(root / ('solution.%08d.domain.%03d.vtk' % (steps[0], rank)), dims)
            (b, other) = read(root / ('solution.%08d.domain.%03d.vtk' % (steps[1], rank)), dims)
            require(mask == other, 'static marker changed')
            require(a['owned_coordinates_sha256'] == b['owned_coordinates_sha256'], 'grid coordinates changed')
            rows.append(dict(rank=rank, before=a, after=b))
            add('%03d.bin' % rank, mask)
            if (rank + 1) % 16 == 0:
                print('extracted static masks', rank + 1, file=sys.stderr, flush=True)
        add('manifest.json', (json.dumps(dict(schema='hundun_vtk_static_mask_v1', steps=steps, ranks=rows), indent=2) + '\n').encode())

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--visit', type=Path, required=True)
    parser.add_argument('--steps', type=int, nargs=2, required=True)
    parser.add_argument('--ranks', type=int, required=True)
    parser.add_argument('--dims', type=int, nargs=3, required=True)
    parser.add_argument('--output', default='-')
    args = parser.parse_args()
    if args.ranks < 1 or min(args.dims) < 3 or min(args.steps) < 0 or (args.steps[0] >= args.steps[1]):
        parser.error('invalid grid, rank or snapshot interval')
    if args.output == '-':
        extract(args.visit, args.steps, args.ranks, args.dims, sys.stdout.buffer)
        return
    output = Path(args.output)
    if output.exists():
        parser.error('output already exists')
    name = None
    try:
        with tempfile.NamedTemporaryFile(dir=str(output.parent), prefix='.mask-', delete=False) as f:
            name = f.name
            extract(args.visit, args.steps, args.ranks, args.dims, f)
            f.flush()
            os.fsync(f.fileno())
        os.link(name, str(output))
    finally:
        if name is not None:
            os.unlink(name)
if __name__ == '__main__':
    main()
