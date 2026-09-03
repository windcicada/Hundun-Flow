#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Extract Parnaudeau Fig. 11--15 filled PIV marker objects reproducibly."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


SOURCE_SHA256 = (
    "b8a775e5a5078e19fc47d9c5f47e95b81b4a68d4449e4444459088ed9befcdd4")
SOURCE_BYTES = 1892798
SVG_HEIGHT_PT = 792.0
STATIONS = (1.06, 1.54, 2.02)
CANONICAL_INVENTORY_SHA256 = (
    "9972aa40e1bcb57513dc739f203257e9e020783622ecf90af055aefad107c537")
EXPECTED_SEQUENCE_COUNTS = {
    11: (130, 101, 158),
    12: (158, 157, 158),
    13: (128, 91, 118),
    14: (158, 139, 97),
    15: (152, 111, 102),
}
EXPECTED_TICK_CENTERS = {
    11: {
        "x": (99.882812, 121.453125, 143.027344, 164.601562,
              186.175781, 207.75, 229.320312, 250.894531, 272.46875),
        "y": (226.503906, 204.933594, 183.359375, 161.785156,
              140.210938, 118.636719, 97.0625, 75.488281, 53.917969),
    },
    12: {
        "x": (363.300781, 384.871094, 406.445312, 428.019531,
              449.59375, 471.167969, 492.738281, 514.3125, 535.886719),
        "y": (226.503906, 197.734375, 168.980469, 140.210938,
              111.441406, 82.6875, 53.917969),
    },
    13: {
        "x": (100.46875, 121.96875, 143.472656, 164.972656,
              186.472656, 207.976562, 229.476562, 250.976562, 272.476562),
        "y": (225.601562, 206.492188, 187.378906, 168.269531,
              149.160156, 130.035156, 110.925781, 91.8125, 72.703125,
              53.59375),
    },
    14: {
        "x": (99.652344, 121.296875, 142.9375, 164.578125,
              186.222656, 207.863281, 229.507812, 251.148438, 272.792969),
        "y": (697.34375, 668.480469, 639.632812, 610.773438,
              581.910156, 553.0625, 524.203125),
    },
    15: {
        "x": (363.886719, 385.386719, 406.890625, 428.390625,
              449.890625, 471.394531, 492.894531, 514.394531, 535.894531),
        "y": (225.601562, 201.027344, 176.453125, 151.878906,
              127.316406, 102.742188, 78.167969, 53.59375),
    },
}
EXPECTED_AXIS_ELEMENT_SHA256 = {
    11: "732caef37073cba2634703a9e0a02094010c553fd317c4657b06e503d55c929d",
    12: "e2112c92b0097e04260035bbcc8884c24c323f618e0f976524c881e0b46a1aef",
    13: "ec823c99b8e1dff1f09d3af240a98fdaed61959657820b785fa0e8f3a5842307",
    14: "45e7abb66696cffd698b36a8d5166b792fddb7d38b8b1a6fd67c80ee83c75dfb",
    15: "839090d73e9c1c669c1cc9eb1dd937563d031963d0dcacf7b948c82aea683e9c",
}
BLACK_CIRCLE_STYLE = (
    " stroke:none;fill-rule:nonzero;fill:rgb(0%,0%,0%);fill-opacity:1;")
PATH_RE = re.compile(
    r'<path style="([^"]+)" d="M ([0-9.]+) ([0-9.]+) C ([^"]+)"/>')
PATH_TAG_RE = re.compile(r'<path\b([^>]*)/>')
ATTRIBUTE_RE = re.compile(r'([A-Za-z_:][-A-Za-z0-9_.:]*)="([^"]*)"')
NUMBER_RE = re.compile(r"-?[0-9]+(?:\.[0-9]+)?")
PATH_TOKEN_RE = re.compile(r"[MLZ]|-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)")


PANELS = (
    {
        "figure": 11, "page": 11, "quantity": "mean_u_over_uc",
        "box": (99.882812, 272.468750, 53.917969, 226.503906),
        "x_ticks": (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0),
        "y_ticks": (-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5),
        "offsets": (0.0, 1.0, 2.0),
    },
    {
        "figure": 12, "page": 11, "quantity": "mean_v_over_uc",
        "box": (363.300781, 535.886719, 53.917969, 226.503906),
        "x_ticks": (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0),
        "y_ticks": (-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5),
        "offsets": (0.0, 1.0, 2.0),
    },
    {
        "figure": 13, "page": 12, "quantity": "uu_over_uc2",
        "box": (100.468750, 272.476562, 53.593750, 225.601562),
        "x_ticks": (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0),
        "y_ticks": (-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3),
        "offsets": (0.0, 0.3, 0.6),
    },
    {
        "figure": 14, "page": 12, "quantity": "vv_over_uc2",
        "box": (99.652344, 272.792969, 524.203125, 697.343750),
        "x_ticks": (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0),
        "y_ticks": (-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2),
        "offsets": (0.0, 0.4, 0.8),
    },
    {
        "figure": 15, "page": 12, "quantity": "uv_over_uc2",
        "box": (363.886719, 535.894531, 53.593750, 225.601562),
        "x_ticks": (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0),
        "y_ticks": (-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1),
        "offsets": (0.0, 0.2, 0.4),
    },
)


class ExtractionError(RuntimeError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    with open(str(path), "w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def run(command, log_path):
    completed = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True)
    with open(str(log_path), "w", encoding="utf-8") as stream:
        stream.write("command: {}\n".format(json.dumps(command)))
        stream.write("exit_code: {}\n".format(completed.returncode))
        stream.write(completed.stdout)
    if completed.returncode != 0:
        raise ExtractionError("command failed; see {}".format(log_path))


def affine_fit(coordinates, values):
    if len(coordinates) != len(values) or len(values) < 3:
        raise ExtractionError("affine calibration requires at least three ticks")
    mean_p = sum(coordinates) / len(coordinates)
    mean_z = sum(values) / len(values)
    denominator = sum((p - mean_p) ** 2 for p in coordinates)
    if denominator == 0.0:
        raise ExtractionError("degenerate tick coordinates")
    slope = sum((p - mean_p) * (z - mean_z)
                for p, z in zip(coordinates, values)) / denominator
    intercept = mean_z - slope * mean_p
    residual = max(abs(slope * p + intercept - z)
                   for p, z in zip(coordinates, values))
    return slope, intercept, residual


def parse_circles(svg_path):
    text = Path(svg_path).read_text(encoding="utf-8")
    circles = []
    for source_index, match in enumerate(PATH_RE.finditer(text)):
        if match.group(1) != BLACK_CIRCLE_STYLE:
            continue
        x0 = float(match.group(2))
        y0 = float(match.group(3))
        tail = [float(value) for value in NUMBER_RE.findall(match.group(4))]
        if len(tail) < 12 or len(tail) % 2:
            continue
        xs = [x0] + tail[0::2]
        ys = [y0] + tail[1::2]
        width = max(xs) - min(xs)
        height = max(ys) - min(ys)
        if not (1.0 <= width <= 1.15 and 1.0 <= height <= 1.15 and
                abs(width - height) <= 0.02):
            continue
        circles.append({
            "source_index": source_index,
            "x_pt": (min(xs) + max(xs)) / 2.0,
            "y_pt": (min(ys) + max(ys)) / 2.0,
            "width_pt": width,
            "height_pt": height,
        })
    return circles


def _path_subpaths(path_d, transform):
    """Parse the M/L/Z-only axis paths emitted by pdftocairo."""
    tokens = PATH_TOKEN_RE.findall(path_d)
    remainder = PATH_TOKEN_RE.sub("", path_d)
    if remainder.strip(" ,\t\r\n"):
        raise ExtractionError("unsupported command in axis path")
    subpaths = []
    current = None
    start = None
    index = 0
    while index < len(tokens):
        command = tokens[index]
        if command == "M":
            if index + 2 >= len(tokens):
                raise ExtractionError("truncated axis move")
            point = (float(tokens[index + 1]), float(tokens[index + 2]))
            current = point
            start = point
            subpaths.append([point])
            index += 3
        elif command == "L":
            if current is None or index + 2 >= len(tokens):
                raise ExtractionError("truncated axis line")
            point = (float(tokens[index + 1]), float(tokens[index + 2]))
            current = point
            subpaths[-1].append(point)
            index += 3
        elif command == "Z":
            if current is None:
                raise ExtractionError("axis close before move")
            if current != start:
                subpaths[-1].append(start)
            current = start
            index += 1
        else:
            raise ExtractionError("unsupported command in axis path")
    if transform is None:
        return subpaths
    if transform != "matrix(1,0,0,-1,0,792)":
        raise ExtractionError("unsupported axis transform")
    return [[(x, round(SVG_HEIGHT_PT - y, 12)) for x, y in subpath]
            for subpath in subpaths]


def _unique_coordinates(values):
    unique = []
    for value in sorted(values):
        if not unique or abs(value - unique[-1]) > 1.0e-7:
            unique.append(value)
    return unique


def _axis_frame_matches(points, panel):
    x0, x1, y0, y1 = panel["box"]
    corners = ((x0, y0), (x0, y1), (x1, y0), (x1, y1))
    return all(any(abs(x - cx) <= 1.0e-7 and abs(y - cy) <= 1.0e-7
                   for x, y in points) for cx, cy in corners)


def _axis_path_ticks(path_d, transform, panel):
    subpaths = _path_subpaths(path_d, transform)
    points = [point for subpath in subpaths for point in subpath]
    if not _axis_frame_matches(points, panel):
        return None
    vertical = []
    horizontal = []
    for subpath in subpaths:
        for left, right in zip(subpath, subpath[1:]):
            dx = abs(left[0] - right[0])
            dy = abs(left[1] - right[1])
            if dx <= 1.0e-7 and 3.0 <= dy <= 5.0:
                vertical.append(left[0])
            elif dy <= 1.0e-7 and 3.0 <= dx <= 5.0:
                horizontal.append(left[1])
    x_coordinates = _unique_coordinates(vertical)
    y_coordinates = list(reversed(_unique_coordinates(horizontal)))
    if (len(x_coordinates) != len(panel["x_ticks"]) or
            len(y_coordinates) != len(panel["y_ticks"])):
        return None
    x0, x1, y0, y1 = panel["box"]
    if (not all(x0 - 1.0e-7 <= x <= x1 + 1.0e-7
                for x in x_coordinates) or
            not all(y0 - 1.0e-7 <= y <= y1 + 1.0e-7
                    for y in y_coordinates)):
        return None
    return x_coordinates, y_coordinates


def axis_tick_inventory(svg_by_page):
    """Recover every major tick centre from each source axis/frame path."""
    paths_by_page = {}
    for page, svg_path in svg_by_page.items():
        text = Path(svg_path).read_text(encoding="utf-8")
        paths_by_page[page] = [
            (index, match, dict(ATTRIBUTE_RE.findall(match.group(1))))
            for index, match in enumerate(PATH_TAG_RE.finditer(text))]
    inventory = []
    for panel in PANELS:
        candidates = []
        for index, match, attributes in paths_by_page[panel["page"]]:
            style = attributes.get("style", "")
            if ("fill:none" not in style or
                    "stroke-miterlimit:10" not in style or
                    "stroke-dasharray" in style):
                continue
            path_d = attributes.get("d")
            if path_d is None:
                continue
            try:
                ticks = _axis_path_ticks(
                    path_d, attributes.get("transform"), panel)
            except ExtractionError:
                continue
            if ticks is None:
                continue
            x_coordinates, y_coordinates = ticks
            element = match.group(0).encode("utf-8")
            element_sha256 = hashlib.sha256(element).hexdigest()
            expected_ticks = EXPECTED_TICK_CENTERS[panel["figure"]]
            if (x_coordinates != list(expected_ticks["x"]) or
                    y_coordinates != list(expected_ticks["y"])):
                raise ExtractionError(
                    "Fig. {} axis tick centers differ from frozen baseline"
                    .format(panel["figure"]))
            if element_sha256 != EXPECTED_AXIS_ELEMENT_SHA256[panel["figure"]]:
                raise ExtractionError(
                    "Fig. {} axis/frame element differs from frozen baseline"
                    .format(panel["figure"]))
            candidates.append({
                "source_identity": {
                    "page": panel["page"],
                    "svg_path_ordinal": index,
                    "element_sha256": element_sha256,
                    "d_sha256": hashlib.sha256(
                        path_d.encode("utf-8")).hexdigest(),
                    "style": style,
                    "transform": attributes.get("transform"),
                },
                "x_tick_centers_pt": x_coordinates,
                "y_tick_centers_pt": y_coordinates,
            })
        if len(candidates) != 1:
            raise ExtractionError(
                "Fig. {} did not yield one complete axis/frame path: {}"
                .format(panel["figure"], len(candidates)))
        tick_record = candidates[0]
        tick_record.update({
            "figure": panel["figure"],
            "page": panel["page"],
            "quantity": panel["quantity"],
            "x_tick_coordinate_value": [
                list(pair) for pair in zip(
                    tick_record["x_tick_centers_pt"], panel["x_ticks"])],
            "y_tick_coordinate_value": [
                list(pair) for pair in zip(
                    tick_record["y_tick_centers_pt"], panel["y_ticks"])],
        })
        inventory.append(tick_record)
    return inventory


def panel_sequences(circles, panel):
    x0, x1, y0, y1 = panel["box"]
    selected = []
    for circle in circles:
        half_w = circle["width_pt"] / 2.0
        half_h = circle["height_pt"] / 2.0
        if (circle["x_pt"] - half_w >= x0 and
                circle["x_pt"] + half_w <= x1 and
                circle["y_pt"] - half_h >= y0 and
                circle["y_pt"] + half_h <= y1):
            selected.append(circle)
    groups = [[]]
    reset_threshold = (x1 - x0) / 4.0
    for circle in selected:
        if (groups[-1] and
                circle["x_pt"] - groups[-1][-1]["x_pt"] > reset_threshold):
            groups.append([])
        groups[-1].append(circle)
    if len(groups) != 3 or any(len(group) < 8 for group in groups):
        raise ExtractionError(
            "Fig. {} did not yield three sufficient marker sequences: {}"
            .format(panel["figure"], [len(group) for group in groups]))
    mean_y = [sum(item["y_pt"] for item in group) / len(group)
              for group in groups]
    if not (mean_y[0] < mean_y[1] < mean_y[2]):
        raise ExtractionError(
            "Fig. {} station sequence order is not top-to-bottom".format(
                panel["figure"]))
    for group in groups:
        prior = None
        for item in group:
            if prior is not None and item["x_pt"] > prior + 1.0e-9:
                raise ExtractionError(
                    "Fig. {} marker sequence is not decreasing".format(
                        panel["figure"]))
            prior = item["x_pt"]
    return groups


def canonical_inventory(svg_by_page):
    circles_by_page = {page: parse_circles(path)
                       for page, path in svg_by_page.items()}
    inventory = []
    for panel in PANELS:
        groups = panel_sequences(circles_by_page[panel["page"]], panel)
        inventory.append({
            "figure": panel["figure"],
            "page": panel["page"],
            "quantity": panel["quantity"],
            "sequences": groups,
        })
    counts = {
        item["figure"]: tuple(len(group) for group in item["sequences"])
        for item in inventory}
    if counts != EXPECTED_SEQUENCE_COUNTS:
        raise ExtractionError(
            "canonical sequence cardinality mismatch: {}".format(counts))
    digest = hashlib.sha256(canonical_bytes(inventory)).hexdigest()
    if digest != CANONICAL_INVENTORY_SHA256:
        raise ExtractionError(
            "canonical marker inventory hash mismatch: {}".format(digest))
    return inventory


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def make_calibration(panel, tick_inventory):
    x_coordinates = tick_inventory["x_tick_centers_pt"]
    y_coordinates = tick_inventory["y_tick_centers_pt"]
    ax, bx, rx = affine_fit(x_coordinates, panel["x_ticks"])
    ay, by, ry = affine_fit(y_coordinates, panel["y_ticks"])
    return {
        "tick_source_identity": tick_inventory["source_identity"],
        "x_axis": {
            "tick_coordinate_value": list(map(list, zip(x_coordinates, panel["x_ticks"]))),
            "slope": ax, "intercept": bx, "max_tick_residual": rx,
        },
        "y_axis": {
            "tick_coordinate_value": list(map(list, zip(y_coordinates, panel["y_ticks"]))),
            "slope": ay, "intercept": by, "max_tick_residual": ry,
        },
    }


def map_inventory(inventory, tick_inventory):
    profiles = []
    calibration_panels = []
    raw_points = []
    by_figure = {item["figure"]: item for item in inventory}
    ticks_by_figure = {item["figure"]: item for item in tick_inventory}
    for panel in PANELS:
        panel_ticks = ticks_by_figure[panel["figure"]]
        calibration = make_calibration(panel, panel_ticks)
        ax = calibration["x_axis"]["slope"]
        bx = calibration["x_axis"]["intercept"]
        rx = calibration["x_axis"]["max_tick_residual"]
        ay = calibration["y_axis"]["slope"]
        by = calibration["y_axis"]["intercept"]
        ry = calibration["y_axis"]["max_tick_residual"]
        groups = by_figure[panel["figure"]]["sequences"]
        panel_bounds = []
        for station_index, (station, offset, group) in enumerate(
                zip(STATIONS, panel["offsets"], groups)):
            data = []
            display_data = []
            coordinate_bounds = []
            value_bounds = []
            mapped_objects = []
            for circle in reversed(group):
                coordinate = ax * circle["x_pt"] + bx
                displayed = ay * circle["y_pt"] + by
                physical = displayed + offset
                q_pt = 0.5e-6
                coordinate_bound = rx + abs(ax) * (
                    q_pt + circle["width_pt"] / 2.0)
                value_bound = ry + abs(ay) * (
                    q_pt + circle["height_pt"] / 2.0)
                data.append([coordinate, physical])
                display_data.append([coordinate, displayed])
                coordinate_bounds.append(coordinate_bound)
                value_bounds.append(value_bound)
                mapped_objects.append({
                    "source_index": circle["source_index"],
                    "x_pt": circle["x_pt"], "y_pt": circle["y_pt"],
                    "width_pt": circle["width_pt"],
                    "height_pt": circle["height_pt"],
                    "y_over_d": coordinate,
                    "printed_ordinate": displayed,
                    "physical_value": physical,
                    "coordinate_bound": coordinate_bound,
                    "value_bound": value_bound,
                })
            for left, right in zip(data, data[1:]):
                if not left[0] < right[0]:
                    raise ExtractionError("mapped coordinate is not increasing")
            coordinate_max = max(coordinate_bounds)
            value_max = max(value_bounds)
            panel_bounds.append({
                "station_x_over_d": station,
                "coordinate_absolute": coordinate_max,
                "value_absolute": value_max,
            })
            profiles.append({
                "station_x_over_d": station,
                "quantity": panel["quantity"],
                "figure": panel["figure"],
                "legend_identity":
                    "filled circles: present PIV; article scope fixes case 2",
                "source_locator":
                    "Parnaudeau et al. (2008), PDF page {}, Fig. {}, x/D={}"
                    .format(panel["page"], panel["figure"], station),
                "extraction_method":
                    "exact filled-circle vector-object centre with affine tick fit",
                "extraction_uncertainty": {
                    "status": "controlled_digitization_bound",
                    "composition": "conservative_l1",
                    "coordinate_absolute": coordinate_max,
                    "value_absolute": value_max,
                    "experimental_uncertainty_status":
                        "not_reported_separately",
                },
                "data": data,
            })
            raw_points.append({
                "figure": panel["figure"], "page": panel["page"],
                "quantity": panel["quantity"],
                "station_index": station_index,
                "station_x_over_d": station,
                "presentation_offset": offset,
                "display_data": display_data,
                "objects": mapped_objects,
            })
        calibration_panels.append({
            "figure": panel["figure"], "page": panel["page"],
            "quantity": panel["quantity"], "panel_box_pt": panel["box"],
            "presentation_offsets": panel["offsets"],
            "tick_inventory": panel_ticks,
            "calibration": calibration, "profile_bounds": panel_bounds,
        })
    if len(profiles) != 15:
        raise ExtractionError("profile matrix is not 3 x 5")
    common_low = max(profile["data"][0][0] for profile in profiles)
    common_high = min(profile["data"][-1][0] for profile in profiles)
    if not common_low < common_high:
        raise ExtractionError("no common y/D interval")
    return profiles, calibration_panels, raw_points, [common_low, common_high]


def render_panel(source, output, panel, dpi):
    # Preserve printed tick labels around the complete panel for human review.
    margin_pt = 24.0
    x0, x1, y0, y1 = panel["box"]
    scale = dpi / 72.0
    crop_x = max(0, int(math.floor((x0 - margin_pt) * scale)))
    crop_y = max(0, int(math.floor((y0 - margin_pt) * scale)))
    crop_w = int(math.ceil((x1 - x0 + 2.0 * margin_pt) * scale))
    crop_h = int(math.ceil((y1 - y0 + 2.0 * margin_pt) * scale))
    command = [
        "pdftocairo", "-png", "-singlefile", "-r", str(dpi),
        "-f", str(panel["page"]), "-l", str(panel["page"]),
        "-x", str(crop_x), "-y", str(crop_y),
        "-W", str(crop_w), "-H", str(crop_h),
        str(source), str(output.with_suffix("")),
    ]
    run(command, output.with_suffix(".log"))
    if not output.exists():
        raise ExtractionError("missing panel rendering {}".format(output))
    return {
        "figure": panel["figure"], "page": panel["page"], "dpi": dpi,
        "crop_pixels": [crop_x, crop_y, crop_w, crop_h],
        "path": str(output), "sha256": sha256(output),
        "log_path": str(output.with_suffix(".log")),
        "log_sha256": sha256(output.with_suffix(".log")),
    }


def _target_marker_tags(text):
    return [
        match for match in PATH_TAG_RE.finditer(text)
        if dict(ATTRIBUTE_RE.findall(match.group(1))).get("style") ==
        BLACK_CIRCLE_STYLE
    ]


def _mutate_marker_text(text, mutation):
    targets = _target_marker_tags(text)
    if len(targets) < 2:
        raise ExtractionError("self-test source has too few marker tags")
    first = targets[0]
    if mutation == "duplicate-marker":
        return text[:first.end()] + first.group(0) + text[first.end():]
    if mutation == "style-change":
        replacement = first.group(0).replace(
            "fill:rgb(0%,0%,0%);", "fill:rgb(1%,1%,1%);", 1)
        return text[:first.start()] + replacement + text[first.end():]
    if mutation == "delete-marker":
        return text[:first.start()] + text[first.end():]
    if mutation == "sequence-boundary-shift":
        attributes = dict(ATTRIBUTE_RE.findall(first.group(1)))
        path_d = attributes["d"]
        number_index = 0

        def shift_x(match):
            nonlocal number_index
            value = float(match.group(0))
            if number_index % 2 == 0:
                value -= 50.0
            number_index += 1
            return "{:.12g}".format(value)

        shifted_d = NUMBER_RE.sub(shift_x, path_d)
        replacement = first.group(0).replace(path_d, shifted_d, 1)
        return text[:first.start()] + replacement + text[first.end():]
    if mutation == "ordering-swap":
        second = targets[1]
        return (text[:first.start()] + second.group(0) +
                text[first.end():second.start()] + first.group(0) +
                text[second.end():])
    raise ExtractionError("unknown self-test mutation {}".format(mutation))


def _mutate_tick_text(text):
    """Move one non-endpoint Fig. 11 major tick in a common-mode copy."""
    for match in PATH_TAG_RE.finditer(text):
        element = match.group(0).encode("utf-8")
        if (hashlib.sha256(element).hexdigest() !=
                EXPECTED_AXIS_ELEMENT_SHA256[11]):
            continue
        attributes = dict(ATTRIBUTE_RE.findall(match.group(1)))
        path_d = attributes["d"]
        old_coordinate = "121.453125"
        new_coordinate = "121.453625"
        if old_coordinate not in path_d:
            raise ExtractionError("self-test axis tick coordinate not found")
        shifted_d = path_d.replace(old_coordinate, new_coordinate)
        replacement = match.group(0).replace(path_d, shifted_d, 1)
        return text[:match.start()] + replacement + text[match.end():]
    raise ExtractionError("self-test axis/frame path not found")


def self_test(source_pdf):
    source = Path(source_pdf).resolve()
    if (not source.is_file() or source.stat().st_size != SOURCE_BYTES or
            sha256(source) != SOURCE_SHA256):
        raise ExtractionError("source PDF identity mismatch")
    with tempfile.TemporaryDirectory(
            prefix="v04-parnaudeau-self-test-", dir="/tmp") as temporary:
        root = Path(temporary)
        clean = {}
        for page in (11, 12):
            svg = root / "page-{}.svg".format(page)
            run(["pdftocairo", "-f", str(page), "-l", str(page), "-svg",
                 str(source), str(svg)], root / "page-{}.log".format(page))
            clean[page] = svg
        clean_inventory = canonical_inventory(clean)
        clean_ticks = axis_tick_inventory(clean)
        map_inventory(clean_inventory, clean_ticks)
        residuals = []
        for panel_ticks in clean_ticks:
            panel = next(item for item in PANELS
                         if item["figure"] == panel_ticks["figure"])
            calibration = make_calibration(panel, panel_ticks)
            residuals.extend((
                calibration["x_axis"]["max_tick_residual"],
                calibration["y_axis"]["max_tick_residual"],))
        if not any(residual > 1.0e-12 for residual in residuals):
            raise ExtractionError("self-test ticks still look synthetic")
        print("self-test clean: PASS canonical/tick/profile checks")
        source_text = clean[11].read_text(encoding="utf-8")
        for mutation in (
                "duplicate-marker", "style-change", "delete-marker",
                "sequence-boundary-shift", "ordering-swap"):
            mutated = root / "mutated-{}.svg".format(mutation)
            mutated.write_text(
                _mutate_marker_text(source_text, mutation), encoding="utf-8")
            try:
                canonical_inventory({11: mutated, 12: clean[12]})
            except ExtractionError:
                print("self-test {}: RED (rejected)".format(mutation))
            else:
                raise ExtractionError(
                    "self-test {} mutation was accepted".format(mutation))
        tick_mutated = root / "mutated-tick-shift.svg"
        tick_mutated.write_text(
            _mutate_tick_text(source_text), encoding="utf-8")
        tick_inventory = canonical_inventory({
            11: tick_mutated, 12: clean[12]})
        if canonical_bytes(tick_inventory) != canonical_bytes(clean_inventory):
            raise ExtractionError(
                "self-test tick-shift changed marker canonical inventory")
        try:
            axis_tick_inventory({11: tick_mutated, 12: clean[12]})
        except ExtractionError:
            print("self-test tick-shift: RED (rejected)")
        else:
            raise ExtractionError("self-test tick-shift mutation was accepted")


def extract(args):
    source = Path(args.source_pdf).resolve()
    output = Path(args.output_dir).resolve()
    if (not source.is_file() or source.stat().st_size != SOURCE_BYTES or
            sha256(source) != SOURCE_SHA256):
        raise ExtractionError("source PDF identity mismatch")
    if output.exists():
        raise ExtractionError("refusing to overwrite existing output directory")
    output.mkdir(parents=True)
    copied_source = output / "source-primary.pdf"
    shutil.copyfile(str(source), str(copied_source))
    passes = {}
    tick_passes = {}
    conversion_records = []
    for pass_name in ("pass-a", "pass-b"):
        pass_dir = output / pass_name
        pass_dir.mkdir()
        svg_by_page = {}
        for page in (11, 12):
            svg = pass_dir / "page-{}.svg".format(page)
            log = pass_dir / "page-{}.conversion.log".format(page)
            run(["pdftocairo", "-f", str(page), "-l", str(page), "-svg",
                 str(copied_source), str(svg)], log)
            svg_by_page[page] = svg
            conversion_records.append({
                "pass": pass_name, "page": page, "path": str(svg),
                "sha256": sha256(svg), "log_path": str(log),
                "log_sha256": sha256(log),
            })
        passes[pass_name] = canonical_inventory(svg_by_page)
        tick_passes[pass_name] = axis_tick_inventory(svg_by_page)
    if canonical_bytes(passes["pass-a"]) != canonical_bytes(passes["pass-b"]):
        raise ExtractionError("cleared vector inventories disagree")
    if canonical_bytes(tick_passes["pass-a"]) != canonical_bytes(
            tick_passes["pass-b"]):
        raise ExtractionError("cleared vector tick inventories disagree")
    profiles, calibrations, mapped_points, interval = map_inventory(
        passes["pass-a"], tick_passes["pass-a"])
    render_records = []
    render_dir = output / "panel-renders"
    render_dir.mkdir()
    for dpi in (1200, 2400):
        for panel in PANELS:
            path = render_dir / "fig{}-{}dpi.png".format(panel["figure"], dpi)
            render_records.append(render_panel(copied_source, path, panel, dpi))
    raw_trace = {
        "schema": "HUNDUN_V04_PARNAUDEAU_VECTOR_RAW_TRACE_V2",
        "source_sha256": SOURCE_SHA256,
        "canonical_inventory_sha256": CANONICAL_INVENTORY_SHA256,
        "selection": {
            "style": BLACK_CIRCLE_STYLE,
            "diameter_pt_interval": [1.0, 1.15],
            "aspect_difference_max_pt": 0.02,
            "sequence_reset_fraction_of_panel_width": 0.25,
        },
        "pass_a_inventory": passes["pass-a"],
        "pass_b_inventory": passes["pass-b"],
        "pass_a_tick_inventory": tick_passes["pass-a"],
        "pass_b_tick_inventory": tick_passes["pass-b"],
        "mapped_points": mapped_points,
    }
    raw_path = output / "raw-trace.json"
    write_json(raw_path, raw_trace)
    profiles_path = output / "profiles.json"
    write_json(profiles_path, {
        "schema": "HUNDUN_V04_PARNAUDEAU_PROFILES_V1",
        "common_y_over_d_interval": interval,
        "profiles": profiles,
    })
    calibration_path = output / "calibration.json"
    write_json(calibration_path, {
        "schema": "HUNDUN_V04_PARNAUDEAU_VECTOR_CALIBRATION_V2",
        "source": {
            "path": str(copied_source), "bytes": copied_source.stat().st_size,
            "sha256": sha256(copied_source),
        },
        "coordinate_system": "SVG page points; origin top-left; y downward",
        "pdf_page_height_pt": SVG_HEIGHT_PT,
        "decimal_half_quantization_pt": 0.5e-6,
        "panels": calibrations,
        "canonical_inventory_sha256": CANONICAL_INVENTORY_SHA256,
        "tick_inventory_pass_a": tick_passes["pass-a"],
        "tick_inventory_pass_b": tick_passes["pass-b"],
        "common_y_over_d_interval": interval,
        "conversion_artifacts": conversion_records,
        "render_artifacts": render_records,
        "raw_trace_path": str(raw_path), "raw_trace_sha256": sha256(raw_path),
        "profiles_path": str(profiles_path),
        "profiles_sha256": sha256(profiles_path),
    })
    write_json(output / "extraction-summary.json", {
        "source_primary_sha256": sha256(copied_source),
        "raw_trace_sha256": sha256(raw_path),
        "calibration_sha256": sha256(calibration_path),
        "profiles_sha256": sha256(profiles_path),
        "common_y_over_d_interval": interval,
        "profile_point_counts": [
            {"figure": profile["figure"],
             "station_x_over_d": profile["station_x_over_d"],
             "count": len(profile["data"])} for profile in profiles],
    })
    print(str(output / "extraction-summary.json"))


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf")
    parser.add_argument("--output-dir")
    parser.add_argument(
        "--self-test", action="store_true",
        help="run clean and six mutation fail-closed checks in /tmp")
    args = parser.parse_args(argv)
    if args.self_test:
        if not args.source_pdf or args.output_dir:
            parser.error("--self-test requires --source-pdf and no --output-dir")
    elif not args.source_pdf or not args.output_dir:
        parser.error("--source-pdf and --output-dir are required")
    return args


def main(argv=None):
    try:
        args = parse_args(argv)
        if args.self_test:
            self_test(args.source_pdf)
        else:
            extract(args)
        return 0
    except (ExtractionError, OSError, ValueError) as error:
        print("v04_parnaudeau_vector_extract: {}".format(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
