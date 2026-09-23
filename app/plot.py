# SPDX-License-Identifier: Apache-2.0
"""Read-only VTK field discovery and bounded scientific slice rendering.

Synchronous functions: call render_slice in the server's worker thread. One
render runs at a time; source files are never modified. Native VTK XML and
legacy binary structured grids retain their coordinates and IBM masks.
"""

from __future__ import annotations
import hashlib
import json
import math
import re
import threading
import time
from functools import lru_cache
from pathlib import Path

MAX_FILES = 1024
MAX_FRAME_BYTES = 3 * 1024**3
MAX_FILE_BYTES = 256 * 1024**2
MAX_POINTS = 12_000_000
MAX_SLICE_CELLS = 800_000
LOCK = threading.Lock()
SUFFIXES = {
    ".vtk",
    ".vtr",
    ".vti",
    ".vtu",
    ".vts",
    ".vtm",
    ".pvtr",
    ".pvti",
    ".pvtu",
    ".pvts",
}
ALIASES = {
    "speed": ("velocity", "Velocity", "U", "01_Velocity"),
    "pressure": (
        "pressure_perturbation",
        "PressureGauge",
        "pi",
        "pressure",
        "Pressure",
        "02_Static_pressure",
    ),
    "temperature": ("temperature", "Temperature", "T", "02_Temperature"),
}
MASKS = ("07_IBM_cell_type", "cell_activity", "fluid_mask", "ibm_fluid")
GHOSTS = ("vtkGhostType", "avtGhostNodes", "avtGhostZones")


class PlotError(RuntimeError):
    def __init__(self, code, message):
        self.code = code
        self.message = message
        super().__init__(message)


def _root(run_path):
    root = Path(run_path).resolve()
    if not root.is_dir():
        raise PlotError("missing_run", "Run directory is unavailable")
    return root


def _inside(path, root):
    path = path.resolve()
    if not path.is_relative_to(root):
        raise PlotError("outside_run", "Field reference leaves the registered run")
    return path


def _identity(path):
    stat = path.stat()
    return str(path), stat.st_size, stat.st_mtime_ns, stat.st_ino


def _step(name):
    match = re.search(r"(?:solution[.]|step[-_.]?)(\d+)", name, re.I) or re.fullmatch(
        r"s(\d+)r\d+\.vtk", name, re.I
    )
    return int(match.group(1)) if match else None


def _units(variable, field):
    if variable == "speed":
        return "m/s", "Velocity magnitude"
    if variable == "temperature":
        return "K", "Temperature"
    if variable == "pressure":
        perturbation = field.lower() in ("pi", "pressure_perturbation", "pressuregauge")
        return "Pa", "Pressure perturbation" if perturbation else "Static pressure"
    native = {
        "Density": ("kg/m³", "Density"),
        "Enthalpy": ("J/kg", "Specific enthalpy"),
        "rho": ("kg/m³", "Density"),
        "h": ("J/kg", "Specific enthalpy"),
    }
    return native.get(field, ("unspecified", field))


def _hidden_field(name):
    # Domain markers and identifiers carry topology, not scalar physics.
    return name in (*MASKS, *GHOSTS) or bool(
        re.search(
            r"ghost|solid|(?:^|[_\W])(?:id|ids)(?:$|[_\W])|(?:node|cell|point|global|material|region)ids?$",
            name,
            re.I,
        )
    )


def _variables(metadata):
    components = metadata["components"]
    result, used = [], set()
    for name, candidates in ALIASES.items():
        field = next(
            (
                f
                for f in candidates
                if components.get(f) == (3 if name == "speed" else 1)
            ),
            None,
        )
        if field:
            units, label = _units(name, field)
            result.append(
                {"name": name, "label": label, "units": units, "field": field}
            )
            used.add(field)
    for field in metadata["fields"]:
        if field in used or components.get(field) != 1 or _hidden_field(field):
            continue
        # Prefix raw fields to keep arbitrary file names distinct from aliases.
        units, label = _units("field:" + field, field)
        result.append(
            {"name": "field:" + field, "label": label, "units": units, "field": field}
        )
    return result


def _field_for(variable, metadata):
    item = next((v for v in _variables(metadata) if v["name"] == variable), None)
    if item is None:
        raise PlotError(
            "invalid_variable",
            "Requested scalar field is absent or unavailable for rendering",
        )
    return item["field"]


def _legacy(path, wanted=()):
    """Index binary legacy blocks with seeks; read only selected arrays/geometry."""
    import numpy as np

    types = {
        "float": ">f4",
        "double": ">f8",
        "unsigned_char": "u1",
        "int": ">i4",
        "unsigned_int": ">u4",
    }
    arrays, fields, component_counts, points, dims = {}, [], {}, None, None
    association, count = "point", 0
    with path.open("rb") as stream:
        if not stream.readline().startswith(b"# vtk DataFile"):
            return None
        stream.readline()
        if stream.readline().strip() != b"BINARY":
            return None
        if stream.readline().strip() != b"DATASET STRUCTURED_GRID":
            return None
        while line := stream.readline():
            parts = line.decode("ascii", errors="strict").split()
            if not parts:
                continue
            key = parts[0]
            if key == "DIMENSIONS":
                dims = tuple(map(int, parts[1:4]))
                if math.prod(dims) > MAX_POINTS:
                    raise PlotError(
                        "resource_limit", "A partition exceeds the point limit"
                    )
            elif key == "POINTS":
                n, kind = int(parts[1]), parts[2]
                size = n * 3 * np.dtype(types[kind]).itemsize
                if stream.tell() + size > path.stat().st_size:
                    raise PlotError("incomplete_frame", "Truncated VTK coordinates")
                if "__geometry__" in wanted:
                    points = (
                        np.frombuffer(stream.read(size), dtype=types[kind])
                        .reshape(-1, 3)
                        .astype(float)
                    )
                else:
                    stream.seek(size, 1)
            elif key in ("POINT_DATA", "CELL_DATA"):
                association, count = (
                    ("point" if key == "POINT_DATA" else "cell"),
                    int(parts[1]),
                )
            elif key in ("SCALARS", "VECTORS"):
                name, kind = parts[1:3]
                components = (
                    3 if key == "VECTORS" else (int(parts[3]) if len(parts) > 3 else 1)
                )
                if key == "SCALARS":
                    if not stream.readline().startswith(b"LOOKUP_TABLE"):
                        raise PlotError("unsupported_data", "Invalid VTK lookup record")
                if kind not in types:
                    raise PlotError(
                        "unsupported_data", "Unsupported VTK scalar encoding"
                    )
                size = count * components * np.dtype(types[kind]).itemsize
                if stream.tell() + size > path.stat().st_size:
                    raise PlotError("incomplete_frame", "Truncated VTK field: " + name)
                fields.append(name)
                component_counts[name] = components
                if name in wanted:
                    value = np.frombuffer(stream.read(size), dtype=types[kind]).astype(
                        float
                    )
                    arrays[name] = (
                        association,
                        value.reshape(-1, components) if components > 1 else value,
                    )
                else:
                    stream.seek(size, 1)
            else:
                return None
    return {
        "points": points,
        "dimensions": dims,
        "fields": fields,
        "components": component_counts,
        "arrays": arrays,
    }


def _xml_components(path):
    """Read numeric array declarations without loading appended field payloads."""
    import xml.etree.ElementTree as ET

    with path.open("rb") as stream:
        header = stream.read(128 * 1024).split(b"<AppendedData", 1)[0]
    result = {}
    for section in re.findall(
        rb"<(?:PointData|CellData)\b[^>]*>(.*?)</(?:PointData|CellData)>", header, re.S
    ):
        for tag in re.findall(rb"<DataArray\b[^>]*>", section):
            attributes = ET.fromstring(tag.rstrip(b">").rstrip(b"/") + b"/>").attrib
            if (
                re.fullmatch(r"(?:Float|U?Int)\d+", attributes.get("type", ""))
                and "Name" in attributes
            ):
                result[attributes["Name"]] = int(
                    attributes.get("NumberOfComponents", "1")
                )
    return result


def _read(path, wanted=()):
    import pyvista as pv

    old = _legacy(path, ("__geometry__", *wanted)) if path.suffix == ".vtk" else None
    if old is not None:
        grid = pv.StructuredGrid()
        grid.points = old["points"]
        grid.dimensions = old["dimensions"]
        for name, (association, values) in old["arrays"].items():
            (grid.point_data if association == "point" else grid.cell_data)[name] = (
                values
            )
        return grid
    if path.suffix in (".vtr", ".vti", ".vts"):
        with path.open("rb") as stream:
            header = stream.read(128 * 1024)
        extent = re.search(rb"<Piece[^>]*Extent=[\"\']([^\"\']+)", header)
        if extent:
            values = list(map(int, extent.group(1).split()))
            if (
                len(values) != 6
                or math.prod(values[i + 1] - values[i] + 1 for i in (0, 2, 4))
                > MAX_POINTS
            ):
                raise PlotError("resource_limit", "VTK extent exceeds the point budget")
    reader = pv.get_reader(path)
    for association in ("point", "cell"):
        disable = getattr(reader, "disable_all_" + association + "_arrays", None)
        if disable:
            disable()
            names = getattr(reader, association + "_array_names", [])
            for name in names:
                if name in wanted:
                    getattr(reader, "enable_" + association + "_array")(name)
    grid = reader.read()
    if isinstance(grid, pv.MultiBlock):
        raise PlotError(
            "unsupported_data", "Resolve multiblock files to leaf partitions first"
        )
    if grid.n_points > MAX_POINTS:
        raise PlotError("resource_limit", "A partition exceeds the point limit")
    return grid


@lru_cache(maxsize=2048)
def _metadata(identity):
    path = Path(identity[0])
    if identity[1] > MAX_FILE_BYTES:
        raise PlotError("resource_limit", "A partition exceeds 256 MiB")
    legacy = _legacy(path) if path.suffix == ".vtk" else None
    if legacy is not None:
        fields = legacy["fields"]
        components = legacy["components"]
    else:
        import pyvista as pv

        reader = pv.get_reader(path)
        fields = list(
            dict.fromkeys(
                [
                    *getattr(reader, "point_array_names", []),
                    *getattr(reader, "cell_array_names", []),
                ]
            )
        )
        components = _xml_components(path) if path.suffix != ".vtk" else {}
    geometry = _read(path)
    times = {}
    for name in ("TimeValue", "TIME", "time", "Time"):
        if name in geometry.field_data and len(geometry.field_data[name]):
            value = float(geometry.field_data[name][0])
            if math.isfinite(value):
                times["time"] = value
    if not fields:
        # Legacy ASCII reader metadata is only populated after Update.
        geometry = _read(
            path, (*sum((list(v) for v in ALIASES.values()), []), *MASKS, *GHOSTS)
        )
        fields = list(
            dict.fromkeys([*geometry.point_data.keys(), *geometry.cell_data.keys()])
        )
    for association in (geometry.point_data, geometry.cell_data):
        for name in association:
            array = association[name]
            if array.dtype.kind in "biuf":
                components[name] = array.shape[1] if array.ndim == 2 else 1
    return {
        "fields": fields,
        "components": components,
        "bounds": list(geometry.bounds),
        "points": geometry.n_points,
        **times,
    }


def _leaves(path, root, depth=0):
    if depth > 6:
        raise PlotError("resource_limit", "Nested VTK references exceed six levels")
    path = _inside(path, root)
    if not path.is_file():
        raise PlotError("incomplete_frame", "Missing field partition: " + path.name)
    if path.suffix in (".vtm", ".pvtr", ".pvti", ".pvtu", ".pvts"):
        import xml.etree.ElementTree as ET

        if path.stat().st_size > 4 * 1024**2:
            raise PlotError("resource_limit", "VTK index is too large")
        result = []
        for node in ET.parse(path).iter():
            reference = node.attrib.get("file") or node.attrib.get("Source")
            if reference:
                result.extend(_leaves(path.parent / reference, root, depth + 1))
                if len(result) > MAX_FILES:
                    raise PlotError("resource_limit", "Too many VTK partitions")
        return result
    return [path]


def _times(root):
    result = {}
    # Bounded read of explicit step/time records; never derive time from dt.
    for name in ("diagnostics.jsonl", "monitor.jsonl"):
        path = root / name
        if not path.is_file():
            continue
        with path.open("rb") as stream:
            if path.stat().st_size > 16 * 1024**2:
                stream.seek(-16 * 1024**2, 2)
                stream.readline()
            for line in stream:
                try:
                    record = json.loads(line)
                    step, value = int(record["step"]), float(record["time"])
                    if math.isfinite(value):
                        result[step] = (value, str(path))
                except (KeyError, TypeError, ValueError):
                    continue
    legacy = root / "monitor" / "chemistry_backend.dat"
    if legacy.is_file():
        with legacy.open("rb") as stream:
            header = stream.readline().decode("ascii", errors="replace").split()
            if header[:3] == ["#", "step", "time"]:
                if legacy.stat().st_size > 16 * 1024**2:
                    stream.seek(-16 * 1024**2, 2)
                    stream.readline()
                for line in stream:
                    try:
                        parts = line.split()
                        step, value = int(parts[0]), float(parts[1].replace(b"D", b"E"))
                        if math.isfinite(value) and step not in result:
                            result[step] = (value, str(legacy))
                    except (IndexError, ValueError):
                        continue
    return result


def _frames(root):
    visit = root / "Visit" if (root / "Visit").is_dir() else root
    groups, counts, published_times = {}, {}, {}
    expected = None
    files = list(visit.iterdir())
    if len(files) > 100_000:
        raise PlotError("resource_limit", "Field directory exceeds listing limit")
    indices = sorted(p for p in files if p.suffix == ".visit")
    for index in indices:
        if index.stat().st_size > 8 * 1024**2:
            raise PlotError("resource_limit", "Visit index is too large")
        lines = index.read_text(encoding="utf-8").splitlines()
        blocks = (
            int(lines[0].split()[1])
            if lines and lines[0].startswith("!NBLOCKS")
            else None
        )
        if index.name == "solution.visit":
            expected = blocks
        index_time = None
        for line in lines:
            if line.startswith("!TIME "):
                try:
                    index_time = float(line.split()[1])
                except (ValueError, IndexError):
                    index_time = None
            if line.strip() and not line.startswith("!"):
                file = _inside(visit / line.strip(), root)
                key = (
                    str(_step(file.name)) if _step(file.name) is not None else file.stem
                )
                groups.setdefault(key, set()).add(file)
                counts[key] = blocks
                if index_time is not None and math.isfinite(index_time):
                    published_times[key] = (index_time, str(index))
    masters = [
        p for p in files if p.suffix in (".vtm", ".pvtr", ".pvti", ".pvtu", ".pvts")
    ]
    for file in masters or [p for p in files if p.suffix in SUFFIXES]:
        step = _step(file.name)
        key = str(step) if step is not None else file.stem
        # Native rank dumps become visible only after their .visit publication.
        if (
            ("-rank-" in file.name or re.fullmatch(r"s\d+r\d+\.vtk", file.name, re.I))
            and indices
            and file not in groups.get(key, set())
        ):
            continue
        groups.setdefault(key, set()).add(_inside(file, root))
    frames = []
    times = _times(root)
    for key, sources in groups.items():
        step = _step(next(iter(sources)).name)
        physical_time, time_source = times.get(
            step, published_times.get(key, (None, "unavailable"))
        )
        try:
            leaves = sorted(set(p for source in sources for p in _leaves(source, root)))
            if (
                len(leaves) > MAX_FILES
                or sum(p.stat().st_size for p in leaves) > MAX_FRAME_BYTES
            ):
                raise PlotError(
                    "resource_limit", "Frame exceeds partition or byte budget"
                )
            count = counts.get(key, expected)
            complete = not count or len(leaves) == count
            frames.append(
                {
                    "id": key,
                    "step": step,
                    "time": physical_time,
                    "time_source": time_source,
                    "files": [str(p) for p in leaves],
                    "source": [str(p) for p in sorted(sources)],
                    "status": "ready" if complete else "incomplete",
                    "partitions": len(leaves),
                }
            )
        except PlotError as error:
            frames.append(
                {
                    "id": key,
                    "step": step,
                    "time": physical_time,
                    "files": [],
                    "source": [],
                    "status": error.code,
                    "error": error.message,
                }
            )
    return sorted(
        frames, key=lambda f: (f["step"] is not None, f["step"] or 0, f["id"])
    )


def discover_fields(run_path):
    root = _root(run_path)
    frames = _frames(root)
    variables = []
    sample = next(
        (f for f in reversed(frames) if f["status"] == "ready" and f["files"]), None
    )
    if sample:
        metadata = _metadata(_identity(Path(sample["files"][0])))
        sample["time"] = metadata.get("time", sample["time"])
        variables = _variables(metadata)
    # Existing assets are source assets only; backend decides which to serve.
    images = [
        {"path": str(p), "source": "existing run asset"}
        for p in root.iterdir()
        if p.suffix.lower() in (".png", ".jpg", ".jpeg") and p.is_file()
    ][:40]
    return {
        "status": "ready" if sample else "no_fields",
        "frames": frames,
        "variables": variables,
        "images": images,
        "limits": {
            "max_frame_bytes": MAX_FRAME_BYTES,
            "max_partition_bytes": MAX_FILE_BYTES,
            "max_points": MAX_POINTS,
            "max_slice_cells": MAX_SLICE_CELLS,
        },
    }


def render_slice(run_path, options, out_dir):
    if not LOCK.acquire(blocking=False):
        raise PlotError("busy", "A field image is already being generated")
    try:
        return _render(run_path, options, out_dir)
    except PlotError:
        raise
    except ImportError as error:
        raise PlotError(
            "dependencies_unavailable",
            "Field plotting requires numpy, matplotlib and pyvista/VTK",
        ) from error
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        raise PlotError("invalid_field", str(error)) from error
    finally:
        LOCK.release()


def _render(run_path, options, out_dir):
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize

    started = time.monotonic()
    root = _root(run_path)
    output = Path(out_dir).resolve()
    if output.is_relative_to(root):
        raise PlotError("invalid_output", "Plot output must be outside the source run")
    frames = _frames(root)
    selector = options.get("frame", options.get("step"))
    frame = next(
        (
            f
            for f in reversed(frames)
            if (selector is None or str(selector) in (f["id"], str(f["step"])))
            and f["status"] == "ready"
        ),
        None,
    )
    if frame is None:
        raise PlotError("missing_frame", "Selected frame is unavailable or incomplete")
    variable = str(options.get("variable", "speed"))
    normal = str(options.get("normal", "z")).lower()
    if normal not in ("x", "y", "z"):
        raise PlotError("invalid_slice", "Normal must be x, y or z")
    axis = "xyz".index(normal)
    files = [Path(p) for p in frame["files"]]
    identities = [_identity(p) for p in files]
    key = hashlib.sha256(
        json.dumps(
            [
                identities,
                variable,
                normal,
                options.get("coordinate", "center"),
                frame["time"],
                4,
            ],
            sort_keys=True,
        ).encode()
    ).hexdigest()[:24]
    png, manifest = output / (key + ".png"), output / (key + ".json")
    if png.is_file() and manifest.is_file():
        result = json.loads(manifest.read_text(encoding="utf-8"))
        result["cache_hit"] = True
        return result
    metadata = []
    for identity in identities:
        if time.monotonic() - started > 120:
            raise PlotError(
                "resource_limit", "Field metadata scan exceeded 120 seconds"
            )
        metadata.append(_metadata(identity))
    if sum(m["points"] for m in metadata) > MAX_POINTS:
        raise PlotError("resource_limit", "Frame exceeds 12 million points")
    bounds = (
        [min(m["bounds"][2 * a] for m in metadata) for a in range(3)],
        [max(m["bounds"][2 * a + 1] for m in metadata) for a in range(3)],
    )
    coordinate = options.get("coordinate", "center")
    coordinate = (
        (bounds[0][axis] + bounds[1][axis]) / 2
        if coordinate in (None, "center")
        else float(coordinate)
    )
    if (
        not math.isfinite(coordinate)
        or not bounds[0][axis] <= coordinate <= bounds[1][axis]
    ):
        raise PlotError(
            "invalid_slice", "Slice coordinate is outside the physical domain"
        )
    direction = np.eye(3)[axis]
    origin = np.zeros(3)
    origin[axis] = coordinate
    projected = [a for a in range(3) if a != axis]
    polygons, values, selected, fields, masks, times = [], [], [], set(), set(), set()
    for path, meta in zip(files, metadata):
        if not meta["bounds"][2 * axis] <= coordinate <= meta["bounds"][2 * axis + 1]:
            continue
        if time.monotonic() - started > 120:
            raise PlotError("resource_limit", "Slice generation exceeded 120 seconds")
        field = _field_for(variable, meta)
        grid = _read(path, (field, *MASKS, *GHOSTS))
        if field not in grid.array_names:
            raise PlotError("missing_variable", "Reader did not return " + field)
        # Mask before cutting: no interpolation across solid/ghost point values.
        keep = np.ones(grid.n_cells, dtype=bool)
        for ghost in GHOSTS:
            if ghost in grid.cell_data:
                keep &= np.asarray(grid.cell_data[ghost]) == 0
        if grid.point_data:
            for name in (*MASKS, *GHOSTS):
                if name in grid.point_data:
                    raw = np.asarray(grid.point_data[name])
                    point_keep = raw > 0.5 if name in MASKS else raw == 0
                    if name == "vtkGhostType":
                        # Duplicate points provide interpolation across MPI interfaces.
                        # Hidden points remain excluded; wholly duplicate cells are removed.
                        point_keep = (raw.astype(np.uint8) & 254) == 0
                        grid.point_data["_owned_point"] = (
                            (raw.astype(np.uint8) & 1) == 0
                        ).astype(float)
                    grid.point_data["_keep_" + name] = point_keep.astype(float)
            averaged = grid.point_data_to_cell_data(pass_point_data=True)
            if "_owned_point" in averaged.cell_data:
                keep &= np.asarray(averaged.cell_data["_owned_point"]) > 0
            for name in (*MASKS, *GHOSTS):
                keyname = "_keep_" + name
                if keyname in averaged.cell_data:
                    keep &= np.asarray(averaged.cell_data[keyname]) >= 1 - 1e-12
        for mask in MASKS:
            if mask in grid.array_names:
                masks.add(mask)
            if mask in grid.cell_data:
                keep &= np.asarray(grid.cell_data[mask]) > 0.5
        if not keep.any():
            continue
        if not keep.all():
            grid = grid.extract_cells(np.flatnonzero(keep))
        cut = grid.slice(normal=direction, origin=origin).triangulate()
        if cut.n_cells == 0:
            continue
        if field in cut.point_data:
            data = np.asarray(cut.point_data[field])
            if variable == "speed":
                if data.ndim != 2 or data.shape[1] != 3:
                    raise PlotError(
                        "invalid_variable", "Velocity is not a three-component vector"
                    )
                data = np.linalg.norm(data, axis=1)
            cut.point_data["_plot_value"] = data
            cut = cut.point_data_to_cell_data()
            data = np.asarray(cut.cell_data["_plot_value"])
        else:
            data = np.asarray(cut.cell_data[field])
            if variable == "speed":
                if data.ndim != 2 or data.shape[1] != 3:
                    raise PlotError(
                        "invalid_variable", "Velocity is not a three-component vector"
                    )
                data = np.linalg.norm(data, axis=1)
        if data.ndim == 2 and data.shape[1] == 1:
            data = data[:, 0]
        if data.ndim != 1 or data.dtype.kind not in "biuf":
            raise PlotError(
                "invalid_variable", "Selected field is not a numeric scalar"
            )
        faces = np.asarray(cut.faces).reshape(-1, 4)[:, 1:]
        polygons.extend(np.asarray(cut.points)[faces][:, :, projected])
        values.extend(data.tolist())
        if len(values) > MAX_SLICE_CELLS:
            raise PlotError("resource_limit", "Slice exceeds polygon budget")
        selected.append(str(path))
        fields.add(field)
        if meta.get("time") is not None:
            times.add(meta["time"])
    if not values:
        raise PlotError("empty_slice", "Slice contains no owned fluid cells")
    data = np.asarray(values)
    if not np.all(np.isfinite(data)):
        raise PlotError("nonfinite_field", "Slice contains non-finite physical values")
    if len(fields) != 1 or len(times) > 1:
        raise PlotError(
            "inconsistent_frame", "Partition field identities or times differ"
        )
    if identities != [_identity(p) for p in files]:
        raise PlotError(
            "source_changed",
            "Field output changed while rendering; retry after publication",
        )
    field = next(iter(fields))
    units, label = _units(variable, field)
    physical_time = next(iter(times)) if times else frame["time"]
    figure = Figure(figsize=(9, 6), dpi=160, layout="constrained")
    FigureCanvasAgg(figure)
    ax = figure.subplots()
    vmin, vmax = float(data.min()), float(data.max())
    artist = PolyCollection(
        polygons,
        array=data,
        cmap="coolwarm" if variable == "pressure" else "viridis",
        norm=Normalize(
            vmin=vmin, vmax=vmax if vmax > vmin else vmin + max(abs(vmin) * 1e-9, 1e-12)
        ),
        edgecolors="none",
        antialiaseds=False,
        rasterized=True,
    )
    ax.add_collection(artist)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_xlabel("xyz"[projected[0]] + " [m]")
    ax.set_ylabel("xyz"[projected[1]] + " [m]")
    timestamp = (
        f"t = {physical_time:.9g} s"
        if physical_time is not None
        else "time metadata unavailable"
    )
    ax.set_title(
        f"{label} | step {frame['step'] if frame['step'] is not None else frame['id']}\n{normal} = {coordinate:.7g} m | {timestamp}",
        fontsize=11,
    )
    figure.colorbar(
        artist,
        ax=ax,
        label=f"{label} [{'original units' if units == 'unspecified' else units}]",
        shrink=0.82,
    )
    figure.text(
        0.01,
        0.005,
        f"{field}; IBM mask: {', '.join(sorted(masks)) or 'not supplied'}; source: {root.name}",
        fontsize=7,
        color="#555555",
    )
    output.mkdir(parents=True, exist_ok=True)
    temporary = output / (key + ".tmp.png")
    figure.savefig(temporary, dpi=160)
    temporary.replace(png)
    result = {
        "path": str(png),
        "step": frame["step"],
        "frame": frame["id"],
        "time": physical_time,
        "variable": variable,
        "field": field,
        "units": units,
        "normal": normal,
        "coordinate": coordinate,
        "minimum": vmin,
        "maximum": vmax,
        "cache_hit": False,
        "provenance": {
            "run": str(root),
            "sources": selected,
            "identities": identities,
            "mask_fields": sorted(masks),
            "mask_policy": "all cell vertices fluid and owned before slicing",
            "interpolation": "VTK physical plane cut; point values averaged per triangle",
            "time_source": "VTK field data"
            if times
            else frame.get("time_source", "unavailable"),
            "polygons": len(values),
        },
    }
    manifest.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
