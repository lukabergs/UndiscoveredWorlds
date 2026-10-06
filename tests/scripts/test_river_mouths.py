"""The river_mouths validation kind of scripts/evaluation/validation.py on a tiny synthetic bundle (8 x 4)."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

from PIL import Image  # noqa: E402

from uwtools import map_rendering  # noqa: E402
from uwtools.bundle import SCHEMA, Bundle  # noqa: E402
from maps import river_mouth_markers  # noqa: E402
from validation import compare_river_mouths, run_validation  # noqa: E402

COLUMNS, ROWS = 8, 4
RADIUS_M = 6371000.0


def write_bundle(directory, fields):
    """A minimal uw-run-bundle/v1 directory holding the given fields ({id: array})."""
    directory = Path(directory)
    entries = []
    for identifier, values in fields.items():
        values = np.asarray(values)
        component, name = identifier.split(".")
        relative = f"fields/{component}/{name}.npy"
        (directory / "fields" / component).mkdir(parents=True, exist_ok=True)
        np.save(directory / relative, values)
        entries.append({"id": identifier, "component": component, "layout": "scalar" if values.ndim == 0 else "grid",
                        "dtype": values.dtype.str, "shape": list(values.shape), "units": "", "file": relative,
                        "hash": "0", "written_by": "test", "description": ""})
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "earth-hydrology",
                "grid": {"columns": COLUMNS, "rows": ROWS, "registration": "cell-centred",
                         "layout": "row-major, north row first"},
                "fields": entries, "tables": []}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return Bundle.open(directory)


def world():
    """Column 0 is sea, columns 1..7 land. Cell centres: longitude -157.5 + 45 x, latitude 67.5 - 45 y.

    Outlet A (x 1, y 1; -112.5, 22.5) flows west, discharge 100; outlet B (x 7, y 2; 157.5, -22.5)
    flows east across the date line, discharge 50. The polar cell (x 1, y 0) points north-west,
    beyond the pole, with discharge 1000: not an outlet, although a wrapped roll would reach sea.
    """
    elevation = np.full((ROWS, COLUMNS), 200.0, dtype=np.float32)
    elevation[:, 0] = -100.0
    direction = np.zeros((ROWS, COLUMNS), dtype=np.uint8)
    discharge = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    direction[1, 1], discharge[1, 1] = 7, 100.0
    direction[2, 7], discharge[2, 7] = 3, 50.0
    direction[0, 1], discharge[0, 1] = 8, 1000.0
    direction[3, 2], discharge[3, 2] = 3, 20.0  # flows to land
    return {"geosphere.base_elevation": elevation, "geosphere.sea_level": np.array(0, dtype=np.int32),
            "hydrosphere.flow_direction": direction, "hydrosphere.discharge": discharge,
            "planetbody.radius": np.array(RADIUS_M, dtype=np.float64)}


RIVERS = [{"name": "A", "lon": -112.5, "lat": 22.5, "discharge_m3_s": 10},
          {"name": "B", "lon": 160.0, "lat": -20.0, "discharge_m3_s": 5}]


def expected_km(lon1, lat1, lon2, lat2):
    """Great-circle distance from the angle between unit vectors (atan2 form, independent of the haversine)."""
    def unit(lon, lat):
        lon, lat = np.deg2rad(lon), np.deg2rad(lat)
        return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    a, b = unit(lon1, lat1), unit(lon2, lat2)
    return RADIUS_M * np.arctan2(np.linalg.norm(np.cross(a, b)), np.dot(a, b)) / 1000


class RiverMouthTests(unittest.TestCase):
    def compare(self, match_km):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / "table.json"
            table.write_text(json.dumps({"description": "test", "rivers": RIVERS}), encoding="utf-8")
            bundle = write_bundle(Path(directory) / "bundle", world())
            spec = {"id": "river_mouths", "kind": "river_mouths", "table": str(table), "top": 2, "match_km": match_km}
            return compare_river_mouths(bundle, spec)["run"]

    def assertDistance(self, actual, expected):
        self.assertLessEqual(abs(actual - expected), 1e-6 * max(abs(expected), 1.0))

    def test_outlets_and_distances(self):
        run = self.compare(500)
        metrics, rivers = run["metrics"], run["rivers"]
        self.assertEqual(metrics["outlets"], 2)
        self.assertEqual(metrics["matched"], 2)
        self.assertEqual(metrics["largest_outlet_discharge_m3_s"], 100.0)
        distance_b = expected_km(160.0, -20.0, 157.5, -22.5)
        self.assertDistance(rivers[0]["distance_km"], 0.0)
        self.assertDistance(rivers[1]["distance_km"], distance_b)
        self.assertDistance(metrics["median_distance_km"], distance_b / 2)
        self.assertEqual([river["outlet_discharge_m3_s"] for river in rivers], [100.0, 50.0])
        self.assertEqual(run["earth_comparison"]["B: matched outlet discharge_m3_s"], 50.0)

    def test_unmatched_beyond_match_distance(self):
        run = self.compare(300)
        self.assertEqual(run["metrics"]["matched"], 1)
        self.assertIsNone(run["rivers"][1]["outlet_discharge_m3_s"])
        self.assertDistance(run["rivers"][1]["distance_km"], expected_km(160.0, -20.0, 157.5, -22.5))

    def test_skipped_without_routing(self):
        fields = world()
        del fields["hydrosphere.discharge"]
        with tempfile.TemporaryDirectory() as directory:
            bundle = write_bundle(directory, fields)
            spec = {"id": "river_mouths", "kind": "river_mouths", "table": "configs/river-mouths.json"}
            result = run_validation(bundle, [spec], None)
        self.assertEqual(result["distributions"], {})
        self.assertEqual(result["skipped"][0]["id"], "river_mouths")

    def test_overlay_markers(self):
        """Rings at the table mouths, squares at the top outlets, on an 8 x 4 raster of 8-pixel cells."""
        markers, _ = river_mouth_markers(self.compare(500))
        scale, background = 8, (40, 60, 80)
        image = Image.new("RGB", (COLUMNS * scale, ROWS * scale), background)
        map_rendering.draw_markers(image, COLUMNS, ROWS, scale, origin=(0, 0), **markers)
        radius = 10  # min(max(2 * scale, 6), 10)
        # Mouth B (160.0, -20.0) has no outlet at its centre: the ring is white at its radius.
        x, y = (160.0 + 180) / 360 * COLUMNS * scale, (90 + 20.0) / 180 * ROWS * scale
        self.assertEqual(image.getpixel((round(x) - radius + 1, round(y))), map_rendering.MOUTH)
        # Outlet A, cell (1, 1): a magenta square at the cell centre.
        self.assertEqual(image.getpixel((round(1.5 * scale), round(1.5 * scale))), map_rendering.OUTLET)
        # A cell far from every marker keeps the raster colour.
        self.assertEqual(image.getpixel((round(4.5 * scale), round(0.5 * scale))), background)


if __name__ == "__main__":
    unittest.main()
