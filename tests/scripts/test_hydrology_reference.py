"""The hydrology_reference validation kinds (erosion plan C1) on tiny synthetic bundles (8 x 4) and
synthetic reference rasters written to a temporary directory."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_river_mouths import COLUMNS, ROWS, write_bundle  # noqa: E402

from PIL import Image  # noqa: E402

from uwtools import map_rendering  # noqa: E402
import hydrology_reference  # noqa: E402
from maps import overlay  # noqa: E402
from validation import run_validation  # noqa: E402

NODATA = -9999.9
C0, C1 = np.cos(np.deg2rad(67.5)), np.cos(np.deg2rad(22.5))  # cosine of the polar and the middle rows
S = 2 * (C0 + C1)  # the cosine weight of one column


def write_reference(directory, name, values, transform=None):
    """A float32 GeoTIFF as the prepared fields: NaN becomes the nodata value."""
    path = Path(directory) / name / f'{COLUMNS}.tif'
    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.where(np.isnan(values), NODATA, values).astype(np.float32)
    with rasterio.open(path, 'w', driver='GTiff', width=COLUMNS, height=ROWS, count=1, dtype='float32',
                       crs='EPSG:4326', transform=transform or from_origin(-180, 90, 45, 45), nodata=NODATA) as target:
        target.write(values, 1)


def rivers():
    """Model (drainage area and discharge alike) and reference rivers.

    Reference cells: (1, 1) 1000 with the same model value (ratio 1); (1, 3) 2e5, a major river with
    no model value at its cell but 5e5 at (1, 4) and 1.5e5 at (2, 2) (best: 0.75); (2, 7) 1e5, major,
    whose only model neighbour is 4e5 at (2, 0) across the date line (factor 4); (0, 6) 3000 on the
    polar row, unmatched: the model's 3000 at (3, 6) would match only if rows wrapped over the pole.
    """
    model = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    model[1, 1], model[1, 4], model[2, 2], model[2, 0], model[3, 6] = 1000, 5e5, 1.5e5, 4e5, 3000
    reference = np.full((ROWS, COLUMNS), np.nan)
    reference[1, 1], reference[1, 3], reference[2, 7], reference[0, 6] = 1000, 2e5, 1e5, 3000
    return model, reference


def routing():
    """Column 0 sea, columns 1..7 land at 200 m (sea level 0).

    Columns 1-2 drain to the outlet (1, 1) (basin A, 8 cells); columns 3-5 drain into the closed
    depression 7 (cells (1..2, 4..5), filled 210, a lake at 205: basin D, 12 cells, endorheic);
    columns 6-7 drain east across the date line to the sea in rows 1 and 2, stop at a pit (0, 7)
    (direction 0) in row 0 and point beyond the south pole from (3, 7): four basins of two cells, the
    first and last endorheic. Depression 9 at (1, 2) overflows (lake at its filled level) and
    depression 11 at (2, 2) has no lake: neither is closed.
    """
    elevation = np.full((ROWS, COLUMNS), 200.0, dtype=np.float32)
    elevation[:, 0] = -100.0
    direction = np.zeros((ROWS, COLUMNS), dtype=np.uint8)
    direction[:, 2] = 7
    direction[0, 1], direction[1, 1], direction[2, 1], direction[3, 1] = 5, 7, 1, 1
    direction[:, 3] = 3
    direction[0, 4], direction[3, 4], direction[0, 5], direction[3, 5] = 5, 1, 5, 1
    direction[1:3, 4:6] = 5  # inside the closed depression: ignored
    direction[:, 6] = 3
    direction[1, 7], direction[2, 7], direction[0, 7], direction[3, 7] = 3, 3, 0, 4
    depression = np.zeros((ROWS, COLUMNS), dtype=np.int32)
    filled = elevation.copy()
    lake = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    depression[1:3, 4:6], filled[1:3, 4:6], lake[2, 4] = 7, 210.0, 205.0
    depression[1, 2], filled[1, 2], lake[1, 2] = 9, 205.0, 205.0
    depression[2, 2], filled[2, 2] = 11, 205.0
    return {'geosphere.base_elevation': elevation, 'geosphere.sea_level': np.array(0, dtype=np.int32),
            'hydrosphere.flow_direction': direction, 'hydrosphere.depression_id': depression,
            'hydrosphere.filled_elevation': filled, 'hydrosphere.lake_surface': lake}


class HydrologyReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.references = self.directory / 'references'
        model, reference = rivers()
        discharge = reference.copy()
        discharge[3, 3] = 0.0  # zero discharge: left out
        write_reference(self.references, 'rivers_upstream_area', reference)
        write_reference(self.references, 'rivers_discharge', discharge)
        basins = np.full((ROWS, COLUMNS), 2.0)
        basins[:, 1:4] = 1.0
        basins[1:, 0] = np.nan  # basin 2 covers the sea cell (0, 0): reference area the model cannot share
        write_reference(self.references, 'basin_level3', basins)
        internal = np.zeros((ROWS, COLUMNS))
        internal[:, 0] = np.nan
        internal[:, 4:6] = 1.0
        internal[0, 7] = 1.0
        write_reference(self.references, 'internal_drainage', internal)
        fields = routing()
        fields['hydrosphere.drainage_area'] = model
        fields['hydrosphere.discharge'] = model
        self.bundle = write_bundle(self.directory / 'bundle', fields)

    def tearDown(self):
        self.temporary.cleanup()

    def compare(self, kind):
        spec = {'id': kind, 'kind': kind, 'reference_directory': str(self.references)}
        return hydrology_reference.compare(kind, self.bundle, spec)

    def assertRiverStatistics(self, run):
        comparison, counts = run['earth_comparison'], run['metrics']
        self.assertAlmostEqual(comparison['median_abs_log10_ratio'], abs(np.log10(0.75)), places=12)
        self.assertEqual(comparison['share_within_factor_2'], 0.5)
        self.assertEqual(comparison['major_river_recall_within_factor_2'], 0.5)
        self.assertEqual(counts, {'reference_cells': 4, 'matched_cells': 3, 'unmatched_cells': 1, 'major_river_cells': 2})

    def test_upstream_area(self):
        result = self.compare('upstream_area')
        self.assertRiverStatistics(result['run'])
        self.assertEqual(result['major_river_cells'], [[-22.5, 22.5], [157.5, -22.5]])

    def test_discharge(self):
        result = self.compare('discharge')
        self.assertIn('not observations', result['reference'])
        self.assertRiverStatistics(result['run'])

    def test_basins(self):
        run = self.compare('basins')['run']
        self.assertAlmostEqual(run['earth_comparison']['best_match_overlap'], 4 * S / (7 * S + C0), places=12)
        self.assertEqual(run['metrics'], {'reference_basins': 2, 'model_basins_on_reference_land': 6, 'model_basins': 6})

    def test_internal_drainage(self):
        run = self.compare('internal_drainage')['run']
        comparison = run['earth_comparison']
        self.assertAlmostEqual(comparison['intersection_over_union'], (2 * S + C0) / (3 * S + 4 * C0), places=12)
        self.assertAlmostEqual(comparison['model_endorheic_land_share'], (3 * S + 4 * C0) / (7 * S), places=12)
        self.assertAlmostEqual(comparison['reference_endorheic_land_share'], (2 * S + C0) / (7 * S), places=12)
        self.assertEqual(run['metrics'], {'model_endorheic_cells': 16, 'reference_endorheic_cells': 9})

    def endorheic(self, fields, name):
        return hydrology_reference.endorheic_land(write_bundle(self.directory / name, fields))

    def test_two_cell_loop_is_exorheic(self):
        """(1, 6) -> (1, 7) -> (1, 6): pointer jumping ends on the loop's own cells, which are no ends."""
        fields = routing()
        fields['hydrosphere.flow_direction'][1, 7] = 7
        endorheic = self.endorheic(fields, 'loop2')
        self.assertFalse(endorheic[1, 6] or endorheic[1, 7])
        self.assertEqual(int(endorheic.sum()), 16)

    def test_three_cell_loop_is_exorheic(self):
        """(1, 6) -> (1, 7) -> (2, 7) -> (1, 6), with (2, 6) draining into the loop."""
        fields = routing()
        fields['hydrosphere.flow_direction'][1, 7] = 5
        fields['hydrosphere.flow_direction'][2, 7] = 8
        endorheic = self.endorheic(fields, 'loop3')
        self.assertFalse(endorheic[1:3, 6:8].any())
        self.assertEqual(int(endorheic.sum()), 16)

    def test_closed_lake_margin(self):
        """Depression 7 (filled 210) is closed when its lake lies more than 0.01 m below the filled level."""
        for level, closed in ((209.995, False), (209.98, True)):
            fields = routing()
            fields['hydrosphere.lake_surface'][2, 4] = level
            bundle = write_bundle(self.directory / f'margin-{level}', fields)
            cells = hydrology_reference.closed_depression_cells(bundle)
            self.assertEqual(bool(cells[1:3, 4:6].all()), closed)
            self.assertEqual(int(cells.sum()), 4 if closed else 0)  # depressions 9 and 11 are never closed

    def test_skipped_without_reference_width(self):
        spec = {'id': 'basins', 'kind': 'basins', 'reference_directory': str(self.directory / 'none')}
        result = run_validation(self.bundle, [spec], None)
        self.assertEqual(result['distributions'], {})
        self.assertIn('no prepared reference at 8 columns', result['skipped'][0]['reason'])

    def test_misaligned_reference(self):
        write_reference(self.references, 'basin_level3', np.ones((ROWS, COLUMNS)), from_origin(-180 + 22.5, 90, 45, 45))
        with self.assertRaisesRegex(ValueError, 'cell centres differ'):
            self.compare('basins')

    def test_skipped_without_field(self):
        fields = routing()
        bundle = write_bundle(self.directory / 'partial', fields)
        spec = {'id': 'discharge_reference', 'kind': 'discharge', 'reference_directory': str(self.references)}
        result = run_validation(bundle, [spec], None)
        self.assertEqual(result['distributions'], {})
        self.assertEqual(result['skipped'][0]['id'], 'discharge_reference')

    def test_overlay_dots(self):
        """A red dot at the centre of each major river cell, on an 8 x 4 raster of 8-pixel cells."""
        result = self.compare('upstream_area')
        specs = [{'id': 'upstream_area', 'kind': 'upstream_area'}]
        markers, note = overlay({'overlay': 'upstream_area'}, specs, {'upstream_area': result}, 'run')
        self.assertIn('100000 km2', note)
        scale, background = 8, (40, 60, 80)
        image = Image.new('RGB', (COLUMNS * scale, ROWS * scale), background)
        map_rendering.draw_markers(image, COLUMNS, ROWS, scale, origin=(0, 0), **markers)
        self.assertEqual(image.getpixel((round(3.5 * scale), round(1.5 * scale))), map_rendering.DOT)
        self.assertEqual(image.getpixel((round(7.5 * scale), round(2.5 * scale))), map_rendering.DOT)
        self.assertEqual(image.getpixel((round(1.5 * scale), round(1.5 * scale))), background)
        self.assertIsNone(overlay({'overlay': 'upstream_area'}, specs, {'upstream_area': result}, 'baseline'))


if __name__ == '__main__':
    unittest.main()
