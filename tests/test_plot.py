# SPDX-License-Identifier: Apache-2.0
"""Small format/mask/cache contract tests; no CFD solver is launched."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pyvista as pv

try:
    from app import plot
except ImportError:
    import plot


class FieldPlotTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.run = self.root / "run"
        self.visit = self.run / "Visit"
        self.visit.mkdir(parents=True)
        self.output = self.root / "plots"

    def tearDown(self):
        self.tmp.cleanup()

    def xml(self):
        grid = pv.ImageData(dimensions=(5, 5, 5), spacing=(0.1, 0.1, 0.1))
        centres = grid.cell_centers().points
        grid.cell_data["velocity"] = np.tile([3.0, 4.0, 0.0], (grid.n_cells, 1))
        grid.cell_data["temperature"] = np.where(centres[:, 0] < 0.2, 9999.0, 300.0)
        grid.cell_data["pressure_perturbation"] = np.full(grid.n_cells, 17.0)
        grid.cell_data["cell_activity"] = (centres[:, 0] >= 0.2).astype(float)
        grid.field_data["TimeValue"] = np.array([0.125])
        path = self.visit / "step-00000000000000000002-rank-00000000.vti"
        grid.save(path)
        (self.visit / "step-00000000000000000002.visit").write_text(
            "!NBLOCKS 1\n" + path.name + "\n"
        )
        return path

    def test_short_native_partition_names(self):
        self.xml()
        original = next(self.visit.glob("*.vti"))
        grid = pv.read(original).cast_to_unstructured_grid()
        grid.field_data.clear()
        # Native short names use a VTK leaf per rank and one publication index.
        for rank in range(2):
            grid.save(self.visit / f"s7r{rank}.vtk")
        grid.save(self.visit / "s8r0.vtk")
        for p in list(self.visit.glob("*.vti")) + list(self.visit.glob("step*.visit")):
            p.unlink()
        (self.visit / "solution.visit").write_text(
            "!NBLOCKS 2\n!TIME 0.7\ns7r0.vtk\ns7r1.vtk\n"
        )
        frames = plot.discover_fields(self.run)["frames"]
        self.assertEqual(frames[0]["time"], 0.7)
        self.assertEqual(
            [(f["step"], f["partitions"], f["status"]) for f in frames],
            [(7, 2, "ready")],
        )

    def test_duplicate_interface_points_support_cut(self):
        grid = pv.ImageData(dimensions=(4, 4, 4), spacing=(0.1, 0.1, 0.1))
        grid.point_data["Temperature"] = np.full(grid.n_points, 300.0)
        grid.point_data["vtkGhostType"] = (grid.points[:, 2] > 0.25).astype(np.uint8)
        leaf = self.visit / "step-1-rank-0.vti"
        grid.save(leaf)
        (self.visit / "step-1.visit").write_text("!NBLOCKS 1\n" + leaf.name + "\n")
        rendered = plot.render_slice(
            self.run, {"variable": "temperature", "coordinate": 0.25}, self.output
        )
        self.assertEqual(rendered["minimum"], 300.0)
        grid.point_data["vtkGhostType"] = (grid.points[:, 2] > 0.25).astype(
            np.uint8
        ) * 2
        grid.save(leaf)
        with self.assertRaises(plot.PlotError):
            plot.render_slice(
                self.run, {"variable": "temperature", "coordinate": 0.25}, self.output
            )

    def test_real_xml_mask_units_time_and_cache(self):
        path = self.xml()
        fields = plot.discover_fields(self.run)
        self.assertEqual(
            {v["name"] for v in fields["variables"]},
            {"speed", "temperature", "pressure"},
        )
        rendered = plot.render_slice(self.run, {"variable": "temperature"}, self.output)
        self.assertEqual(rendered["minimum"], 300.0)
        self.assertEqual(rendered["maximum"], 300.0)
        self.assertEqual(rendered["time"], 0.125)
        self.assertEqual(rendered["step"], 2)
        self.assertEqual(rendered["units"], "K")
        self.assertEqual(rendered["provenance"]["mask_fields"], ["cell_activity"])
        self.assertTrue(Path(rendered["path"]).read_bytes().startswith(b"\x89PNG"))
        self.assertTrue(
            plot.render_slice(self.run, {"variable": "temperature"}, self.output)[
                "cache_hit"
            ]
        )
        speed = plot.render_slice(self.run, {"variable": "speed"}, self.output)
        self.assertEqual(speed["minimum"], 5.0)
        before = rendered["path"]
        grid = pv.read(path)
        grid.cell_data["temperature"][:] = 400.0
        grid.save(path)
        changed = plot.render_slice(self.run, {"variable": "temperature"}, self.output)
        self.assertNotEqual(before, changed["path"])
        self.assertEqual(changed["minimum"], 400.0)

    def test_missing_incomplete_and_output_safety(self):
        self.xml()
        with self.assertRaises(plot.PlotError) as caught:
            plot.render_slice(self.run, {}, self.run / "generated")
        self.assertEqual(caught.exception.code, "invalid_output")
        with self.assertRaises(plot.PlotError) as caught:
            plot.render_slice(self.run, {"coordinate": 500.0}, self.output)
        self.assertEqual(caught.exception.code, "invalid_slice")
        index = next(self.visit.glob("*.visit"))
        index.write_text(index.read_text().replace("!NBLOCKS 1", "!NBLOCKS 2"))
        self.assertEqual(plot.discover_fields(self.run)["status"], "no_fields")
        with self.assertRaises(plot.PlotError):
            plot.render_slice(self.run, {}, self.output)

    def test_visit_escape_and_single_task(self):
        (self.visit / "solution.visit").write_text("!NBLOCKS 1\n../../escape.vtk\n")
        with self.assertRaises(plot.PlotError) as caught:
            plot.discover_fields(self.run)
        self.assertEqual(caught.exception.code, "outside_run")
        plot.LOCK.acquire()
        try:
            with self.assertRaises(plot.PlotError) as caught:
                plot.render_slice(self.run, {}, self.output)
            self.assertEqual(caught.exception.code, "busy")
        finally:
            plot.LOCK.release()

    def test_binary_structured_and_xml_multiblock(self):
        points = np.array(
            [
                [x, y, z]
                for z in (0.0, 1.0, 2.0)
                for y in (0.0, 1.0, 2.0)
                for x in (0.0, 1.0, 2.0)
            ]
        )
        file = self.visit / "solution.00000011.domain.000.vtk"
        with file.open("wb") as stream:
            stream.write(
                b"# vtk DataFile Version 2.0\nfixture\nBINARY\nDATASET STRUCTURED_GRID\nDIMENSIONS 3 3 3\nPOINTS 27 float\n"
            )
            stream.write(points.astype(">f4").tobytes())
            stream.write(b"\nPOINT_DATA 27\nVECTORS 01_Velocity float\n")
            stream.write(np.tile([3.0, 4.0, 0.0], (27, 1)).astype(">f4").tobytes())
            stream.write(b"\nSCALARS 07_IBM_cell_type float 1\nLOOKUP_TABLE default\n")
            stream.write(np.ones(27, dtype=">f4").tobytes())
            stream.write(b"\n")
        (self.visit / "solution.visit").write_text("!NBLOCKS 1\n" + file.name + "\n")
        fields = plot.discover_fields(self.run)
        self.assertEqual(fields["variables"][0]["name"], "speed")
        result = plot.render_slice(self.run, {"frame": "11"}, self.output)
        self.assertEqual(result["maximum"], 5.0)
        self.assertIsNone(result["time"])
        (self.run / "diagnostics.jsonl").write_text(
            json.dumps({"step": 11, "time": 0.75}) + "\n"
        )
        timed = plot.render_slice(self.run, {"frame": "11"}, self.output)
        self.assertEqual(timed["time"], 0.75)
        self.assertNotEqual(result["path"], timed["path"])

    def test_multiblock_index(self):
        grid = pv.ImageData(dimensions=(3, 3, 3))
        grid.cell_data["temperature"] = np.full(grid.n_cells, 310.0)
        pv.MultiBlock([grid]).save(self.visit / "step-5.vtm")
        result = plot.render_slice(self.run, {"variable": "temperature"}, self.output)
        self.assertEqual(result["step"], 5)
        self.assertEqual(result["minimum"], 310.0)
        self.assertEqual(len(result["provenance"]["sources"]), 1)

    def test_budget_rejects_before_read(self):
        self.xml()
        with patch.object(plot, "MAX_FRAME_BYTES", 1):
            self.assertEqual(plot.discover_fields(self.run)["status"], "no_fields")
            with self.assertRaises(plot.PlotError):
                plot.render_slice(self.run, {}, self.output)

    def test_native_additional_scalars_and_rejected_markers(self):
        path = self.xml()
        grid = pv.read(path)
        for name, value in {
            "rho": 1.2,
            "h": 12000.0,
            "TCR_custom": 0.4,
            "GlobalCellId": 7.0,
            "solid_mask": 1.0,
        }.items():
            grid.cell_data[name] = np.full(grid.n_cells, value)
        grid.cell_data["custom_vector"] = np.ones((grid.n_cells, 2))
        grid.save(path)
        variables = plot.discover_fields(self.run)["variables"]
        self.assertEqual(
            [v["name"] for v in variables[:3]], ["speed", "pressure", "temperature"]
        )
        extras = {v["field"]: v for v in variables[3:]}
        self.assertEqual(set(extras), {"rho", "h", "TCR_custom"})
        for field, units, value in [
            ("rho", "kg/m³", 1.2),
            ("h", "J/kg", 12000.0),
            ("TCR_custom", "unspecified", 0.4),
        ]:
            result = plot.render_slice(
                self.run, {"variable": extras[field]["name"]}, self.output
            )
            self.assertEqual(result["field"], field)
            self.assertEqual(result["units"], units)
            self.assertAlmostEqual(result["minimum"], value)
            self.assertAlmostEqual(result["maximum"], value)
        for field in ("custom_vector", "GlobalCellId", "solid_mask", "absent"):
            with self.assertRaises(plot.PlotError) as caught:
                plot.render_slice(self.run, {"variable": "field:" + field}, self.output)
            self.assertEqual(caught.exception.code, "invalid_variable")


if __name__ == "__main__":
    unittest.main()
