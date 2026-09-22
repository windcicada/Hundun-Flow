#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Create the audited watertight SI STL used by the native 624CF case."""

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import struct


SOURCE_SHA256 = "307879356b66179c1129dedd6f25213ecf240b5093d994cd0e9d361b505df65b"
TARGET_SHA256 = "c1d3c7dfaa402931dae23e73677b7d4f42206a6977d02f6f308b841fc280b10f"
SOURCE_TRIANGLES = 280_974
TARGET_TRIANGLES = 280_968
TARGET_EDGES = 421_452
REMOVED_TRIANGLES = (143935, 144161, 145448, 157028, 157226, 158867)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def digest_file(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def subtract(left, right):
    return tuple(a - b for a, b in zip(left, right))


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def signed_volume(points):
    return dot(points[0], cross(subtract(points[1], points[0]),
                                subtract(points[2], points[0]))) / 6.0


def read_source(path):
    raw = path.read_bytes()
    require(digest_bytes(raw) == SOURCE_SHA256, "624CF source STL SHA-256 mismatch")
    require(len(raw) >= 84, "truncated binary STL")
    triangle_count, = struct.unpack_from("<I", raw, 80)
    require(triangle_count == SOURCE_TRIANGLES, "unexpected source triangle count")
    require(len(raw) == 84 + 50 * triangle_count, "invalid binary STL size")

    triangles = []
    for index in range(triangle_count):
        values = struct.unpack_from("<9f", raw, 96 + 50 * index)
        require(all(math.isfinite(value) for value in values),
                "non-finite source vertex")
        triangles.append(tuple(tuple(value * 0.001 for value in values[offset:offset + 3])
                               for offset in (0, 3, 6)))
    return raw, triangles


def inspect_topology(triangles):
    edges = defaultdict(list)
    duplicate_faces = defaultdict(list)
    for index, points in enumerate(triangles):
        duplicate_faces[tuple(sorted(points))].append(index)
        for left, right in zip(points, points[1:] + points[:1]):
            edges[tuple(sorted((left, right)))].append(index)

    duplicates = [indices for indices in duplicate_faces.values() if len(indices) > 1]
    require(tuple(sorted(index for group in duplicates for index in group)) ==
            REMOVED_TRIANGLES, "unexpected duplicate-triangle topology")
    require(all(len(group) == 2 for group in duplicates),
            "duplicate-triangle group is not a pair")

    adjacency = defaultdict(list)
    for edge, owners in edges.items():
        if len(owners) == 1:
            left, right = edge
            adjacency[left].append(right)
            adjacency[right].append(left)
        else:
            require(len(owners) in (2, 4), "unexpected source edge incidence")
    require(len(adjacency) == 12 and all(len(neighbours) == 2
                                         for neighbours in adjacency.values()),
            "unexpected source boundary loops")

    loops = []
    visited = set()
    for start in adjacency:
        if start in visited:
            continue
        loop = []
        previous = None
        point = start
        while point not in visited:
            visited.add(point)
            loop.append(point)
            following = next(value for value in adjacency[point] if value != previous)
            previous, point = point, following
        require(point == start and len(loop) == 6, "invalid source boundary loop")
        loops.append(loop)
    require(len(loops) == 2, "unexpected boundary-loop count")
    return edges, duplicates, loops


def repair(triangles, edges, duplicates, loops):
    removed = set()
    for indices in duplicates:
        left, right = (triangles[index] for index in indices)
        require(any(left == (right[offset], right[(offset + 2) % 3],
                              right[(offset + 1) % 3]) for offset in range(3)),
                "duplicate faces do not have opposite winding")
        require(all(len(edges[tuple(sorted((a, b)))]) in (2, 4)
                    for a, b in zip(left, left[1:] + left[:1])),
                "duplicate face is not an attached zero-thickness pair")
        removed.update(indices)
    require(tuple(sorted(removed)) == REMOVED_TRIANGLES, "repair removal identity")

    mapping = {}
    seam_pairs = []
    for loop in loops:
        for left_index, right_index in ((0, 4), (1, 3)):
            left, right = loop[left_index], loop[right_index]
            midpoint = tuple((a + b) / 2.0 for a, b in zip(left, right))
            gap = math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
            require(gap < 0.00018, "624CF seam exceeds bounded repair distance")
            require(left not in mapping and right not in mapping,
                    "source seam vertex selected more than once")
            mapping[left] = mapping[right] = midpoint
            seam_pairs.append({"a": left, "b": right, "target": midpoint, "gap_m": gap})

    repaired = []
    changed = []
    for index, points in enumerate(triangles):
        if index in removed:
            continue
        candidate = tuple(mapping.get(point, point) for point in points)
        require(len(set(candidate)) == 3, "repair created a degenerate triangle")
        old_normal = cross(subtract(points[1], points[0]),
                           subtract(points[2], points[0]))
        new_normal = cross(subtract(candidate[1], candidate[0]),
                           subtract(candidate[2], candidate[0]))
        require(dot(old_normal, new_normal) > 0, "repair changed face orientation")
        if points != candidate:
            changed.append(index)
        repaired.append(candidate)
    require(len(repaired) == TARGET_TRIANGLES, "repaired triangle count")

    repaired_edges = defaultdict(list)
    for index, points in enumerate(repaired):
        for left, right in zip(points, points[1:] + points[:1]):
            repaired_edges[tuple(sorted((left, right)))].append((index, left < right))
    require(len(repaired_edges) == TARGET_EDGES, "repaired edge count")
    require(all(len(owners) == 2 and owners[0][1] != owners[1][1]
                for owners in repaired_edges.values()),
            "repaired surface is not closed with opposite edge winding")
    require(len(set(tuple(sorted(points)) for points in repaired)) == len(repaired),
            "repair retained a duplicate face")
    return repaired, seam_pairs, changed


def write_ascii_stl(path, triangles):
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("solid cf_native\n")
        for points in triangles:
            stream.write("facet normal 0 0 0\nouter loop\n")
            for point in points:
                stream.write("vertex " + " ".join(format(value, ".17g")
                                                    for value in point) + "\n")
            stream.write("endloop\nendfacet\n")
        stream.write("endsolid cf_native\n")


def run(args):
    source = args.source.resolve()
    target = args.output.resolve()
    report_path = args.report.resolve()
    require(source.is_file(), "source STL does not exist")
    require(not target.exists(), "output STL already exists")
    require(not report_path.exists(), "report already exists")
    target.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    raw, triangles = read_source(source)
    edges, duplicates, loops = inspect_topology(triangles)
    repaired, seam_pairs, changed = repair(triangles, edges, duplicates, loops)
    write_ascii_stl(target, repaired)
    target_sha = digest_file(target)
    require(target_sha == TARGET_SHA256, "repaired STL SHA-256 mismatch")
    require(digest_file(source) == SOURCE_SHA256, "source STL changed during repair")

    original_volume = math.fsum(signed_volume(points) for points in triangles)
    repaired_volume = math.fsum(signed_volume(points) for points in repaired)
    report = {
        "schema": "hundun.624cf.stl_repair.v1",
        "source": str(source),
        "output": str(target),
        "source_sha256": digest_bytes(raw),
        "target_sha256": target_sha,
        "removed_opposite_duplicate_triangles": list(REMOVED_TRIANGLES),
        "seam_pairs": seam_pairs,
        "changed_triangles": changed,
        "triangles": len(repaired),
        "edges": TARGET_EDGES,
        "original_signed_volume_m3": original_volume,
        "repaired_signed_volume_m3": repaired_volume,
        "volume_change_m3": repaired_volume - original_volume,
        "relative_volume_change": abs((repaired_volume - original_volume) /
                                      original_volume),
        "max_vertex_movement_m": max(pair["gap_m"] / 2.0 for pair in seam_pairs),
        "scope": "exact edge and winding closure; all retained face orientations preserved",
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path,
                        help="original binary cf2222.STL in millimetres")
    parser.add_argument("--output", required=True, type=Path,
                        help="new watertight ASCII STL in metres")
    parser.add_argument("--report", required=True, type=Path,
                        help="new JSON identity and repair audit")
    args = parser.parse_args()
    try:
        run(args)
    except (OSError, ValueError, struct.error) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
