#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
"""Coarsen a cylinder mesh outside its protected cylinder/wake interval."""
import argparse
import bisect
import json
from pathlib import Path


def read_axes(path):
    tokens = iter(path.read_text().split())
    assert [next(tokens), next(tokens), next(tokens)] == [
        'COAST_RUNTIME_AXES', '1', 'grid']
    counts = [int(next(tokens)) for _ in range(3)]
    axes = []
    for label, count in zip('xyz', counts):
        assert next(tokens) == label
        assert int(next(tokens)) == count + 1
        axes.append([float(next(tokens)) for _ in range(count + 1)])
    return axes


def tail(length, width, count):
    low, high = 1.0, 2.0
    while sum(width * high ** k for k in range(1, count + 1)) < length:
        high *= 2.0
    for _ in range(80):
        ratio = .5 * (low + high)
        if sum(width * ratio ** k for k in range(1, count + 1)) < length:
            low = ratio
        else:
            high = ratio
    ratio = .5 * (low + high)
    widths = [width * ratio ** k for k in range(1, count + 1)]
    assert abs(sum(widths) - length) < 1e-12
    return widths, ratio


def coarsen(axis, lower, upper, count):
    left = bisect.bisect_right(axis, lower) - 1
    right = bisect.bisect_left(axis, upper)
    core = axis[left:right + 1]
    remaining = count - (right - left)
    candidates = []
    for nleft in range(2, remaining - 1):
        wl, rl = tail(core[0] - axis[0], core[1] - core[0], nleft)
        wr, rr = tail(axis[-1] - core[-1], core[-1] - core[-2], remaining - nleft)
        candidates.append((max(rl, rr), wl, wr))
    growth, wl, wr = min(candidates, key=lambda item: item[0])
    left_nodes, position = [], core[0]
    for width in wl:
        position -= width
        left_nodes.append(position)
    right_nodes, position = [], core[-1]
    for width in wr:
        position += width
        right_nodes.append(position)
    result = list(reversed(left_nodes)) + core + right_nodes
    result[0], result[-1] = axis[0], axis[-1]
    assert len(result) == count + 1
    assert all(b > a for a, b in zip(result, result[1:]))
    assert result[len(wl):len(wl) + len(core)] == core
    return result, {'protected': [core[0], core[-1]],
                    'preserved_cells': len(core) - 1,
                    'tail_growth': growth}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    old = read_axes(args.source)
    x, mx = coarsen(old[0], -.012, .09, 448)
    y, my = coarsen(old[1], -.03, .03, 280)
    assert len(old[2]) == 81
    axes = [x, y, old[2][::2]]
    counts = [len(a) - 1 for a in axes]
    lines = ['COAST_RUNTIME_AXES 1', 'grid ' + ' '.join(map(str, counts))]
    for label, axis in zip('xyz', axes):
        lines.append('{} {}'.format(label, len(axis)))
        lines.extend(' '.join('{:.17g}'.format(v) for v in axis[i:i + 4])
                     for i in range(0, len(axis), 4))
    args.output.write_text('\n'.join(lines) + '\n')
    print(json.dumps({'counts': counts, 'cells': counts[0] * counts[1] * counts[2],
                      'x': mx, 'y': my, 'spanwise_stride': 2},
                     indent=2))


if __name__ == '__main__':
    main()
