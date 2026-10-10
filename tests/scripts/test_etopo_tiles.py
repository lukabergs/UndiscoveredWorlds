"""ETOPO tile names, boxes across the date line, the Earth-window catalogue and the mosaic.

No network and no refs/ data: read_box runs on tiny synthetic GeoTIFFs in a temporary
directory (15 degree tiles of 2 x 2 pixels, a globe of 8 x 4 pixels).
"""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from uwtools import etopo

SPEC = importlib.util.spec_from_file_location('fetch_etopo_tiles',
    Path(__file__).resolve().parents[2] / 'scripts/refs/fetch_etopo_tiles.py')
FETCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FETCH)
TIFF = b'II*\x00' + bytes(60)


def write_tiff(path, data, west, north, pixel):
    import rasterio
    from rasterio.transform import Affine
    data = np.asarray(data, dtype=np.float32)
    with rasterio.open(path, 'w', driver='GTiff', height=data.shape[0], width=data.shape[1],
                       count=1, dtype='float32', crs='EPSG:4326', nodata=etopo.NODATA,
                       transform=Affine(pixel, 0, west, 0, -pixel, north)) as dataset:
        dataset.write(data, 1)


class TileNameTests(unittest.TestCase):
    def test_points_inside_and_on_edges(self):
        self.assertEqual(etopo.tile_name(7.5, 37.5), 'N45E000')
        self.assertEqual(etopo.tile_name(0.0, 0.0), 'N00E000')
        # A west edge belongs to its tile, a north edge to the tile below it.
        self.assertEqual(etopo.tile_name(15.0, 45.0), 'N45E015')
        self.assertEqual(etopo.tile_name(14.999, 45.001), 'N60E000')
        self.assertEqual(etopo.tile_name(0.0, 30.0), 'N30E000')

    def test_negative_coordinates(self):
        self.assertEqual(etopo.tile_name(-0.001, -0.001), 'N00W015')
        self.assertEqual(etopo.tile_name(-15.0, -15.0), 'S15W015')
        self.assertEqual(etopo.tile_name(-65.5, -19.0), 'S15W075')
        self.assertEqual(etopo.tile_name(-74.0, -46.5), 'S45W075')

    def test_poles_and_antimeridian(self):
        self.assertEqual(etopo.tile_name(0.0, 90.0), 'N90E000')
        self.assertEqual(etopo.tile_name(0.0, -90.0), 'S75E000')
        self.assertEqual(etopo.tile_name(180.0, 10.0), 'N15W180')
        self.assertEqual(etopo.tile_name(-180.0, 10.0), 'N15W180')
        self.assertEqual(etopo.tile_name(179.999, 10.0), 'N15E165')
        with self.assertRaises(ValueError):
            etopo.tile_name(0.0, 90.5)

    def test_bounds_url_and_path(self):
        self.assertEqual(etopo.tile_bounds('S15W075'), (-75, -30, -60, -15))
        self.assertEqual(etopo.tile_bounds('N90W180'), (-180, 75, -165, 90))
        self.assertTrue(etopo.tile_url('N45E000').endswith(
            '/15s_surface_elev_gtif/ETOPO_2022_v1_15s_N45E000_surface.tif'))
        self.assertEqual(etopo.tile_path('N45E000', 'x').name, 'ETOPO_2022_v1_15s_N45E000_surface.tif')
        for name in ('N45E180', 'S90E000', 'N44E000', 'X45E000', 'N45E00'):
            with self.assertRaises(ValueError):
                etopo.tile_bounds(name)

    def test_every_tile_named_once(self):
        names = {etopo.tile_name(lon + 7.5, lat - 7.5)
                 for lon in range(-180, 180, 15) for lat in range(90, -90, -15)}
        self.assertEqual(len(names), 288)
        self.assertEqual(len(etopo.tiles_for_box(-180, -90, 180, 90)), 288)


class TilesForBoxTests(unittest.TestCase):
    def test_box_on_tile_edges_takes_no_neighbours(self):
        self.assertEqual(etopo.tiles_for_box(0, 30, 15, 45), ['N45E000'])
        self.assertEqual(etopo.tiles_for_box(-1, 29, 1, 31),
                         ['N45W015', 'N45E000', 'N30W015', 'N30E000'])

    def test_across_the_date_line(self):
        self.assertEqual(etopo.tiles_for_box(175, 49, -170, 55), ['N60E165', 'N60W180'])
        self.assertEqual(etopo.tiles_for_box(170, -20, -170, -10),
                         ['N00E165', 'N00W180', 'S15E165', 'S15W180'])
        self.assertEqual(etopo.tiles_for_box(180, 0, -165, 15), ['N15W180'])

    def test_invalid_boxes(self):
        for box in ((0, 10, 10, 10), (0, 20, 10, 10), (5, 0, 5, 10), (0, -91, 10, 0), (0, 0, 181, 10)):
            with self.assertRaises(ValueError):
                etopo.tiles_for_box(*box)


class CatalogueTests(unittest.TestCase):
    KEYS = {'name', 'title', 'centre', 'box', 'metres_per_pixel', 'features'}

    def test_schema(self):
        windows = etopo.load_catalogue()
        self.assertGreaterEqual(len(windows), 13)
        names = [window['name'] for window in windows]
        self.assertEqual(len(names), len(set(names)))
        for window in windows:
            with self.subTest(window=window.get('name')):
                self.assertLessEqual(self.KEYS, set(window))
                self.assertRegex(window['name'], r'^[a-z][a-z0-9]*(_[a-z0-9]+)*$')
                self.assertTrue(window['title'].strip())
                west, south, east, north = window['box']
                etopo.check_box(west, south, east, north)
                width = east - west if west < east else east + 360 - west
                self.assertLessEqual(width, 15)
                self.assertLessEqual(north - south, 15)
                lon, lat = window['centre']
                self.assertTrue(etopo.box_contains(window['box'], lon, lat))
                low, high = window['metres_per_pixel']
                self.assertTrue(0 < low <= high)
                self.assertTrue(window['features'])
                self.assertTrue(all(isinstance(f, str) and f.strip() for f in window['features']))
                self.assertLessEqual(len(etopo.tiles_for_box(*window['box'])), 4)

    def test_the_plan_windows_are_present(self):
        names = {window['name'] for window in etopo.load_catalogue()}
        for expected in ('zagros', 'eastern_tibet_three_rivers', 'bolivian_orocline',
                         'norwegian_fjords', 'patagonian_fjords', 'aleutian_arc', 'mariana_arc',
                         'indonesian_arcs', 'scottish_highlands', 'appalachian_fold_belt',
                         'canadian_shield_lakes', 'dalmatian_coast', 'hawaii'):
            self.assertIn(expected, names)


class ReadBoxTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def tile(self, name, data):
        west, _, _, north = etopo.tile_bounds(name)
        write_tiff(etopo.tile_path(name, self.root), data, west, north, 7.5)

    def test_tiles_mosaic_across_the_date_line(self):
        self.tile('N15E165', [[1, 2], [3, etopo.NODATA]])
        self.tile('N15W180', [[5, 6], [7, 8]])
        with patch.object(etopo, 'TILE_ROOT', self.root):
            elevation, extent = etopo.read_box(166, 1, -166, 14, 'etopo15')
        self.assertEqual(elevation.dtype, np.float32)
        expected = np.array([[1, 2, 5, 6], [3, np.nan, 7, 8]], dtype=np.float32)
        np.testing.assert_array_equal(elevation, expected)
        self.assertEqual(extent, dict(west=165.0, south=0.0, east=-165.0, north=15.0, pixel_deg=7.5))

    def test_tiles_mosaic_north_south_and_snap_outwards(self):
        self.tile('N30E000', [[1, 2], [3, 4]])
        self.tile('N15E000', [[5, 6], [7, 8]])
        with patch.object(etopo, 'TILE_ROOT', self.root):
            elevation, extent = etopo.read_box(8, 6, 14, 20, 'etopo15')
        np.testing.assert_array_equal(elevation, np.array([[4], [6], [8]], dtype=np.float32))
        self.assertEqual(extent, dict(west=7.5, south=0.0, east=15.0, north=22.5, pixel_deg=7.5))

    def test_missing_tile_is_reported(self):
        self.tile('N15E165', [[1, 2], [3, 4]])
        with patch.object(etopo, 'TILE_ROOT', self.root):
            with self.assertRaisesRegex(FileNotFoundError, 'N15W180'):
                etopo.read_box(166, 1, -166, 14, 'etopo15')

    def test_globe_across_the_date_line_and_at_the_poles(self):
        globe = np.arange(32, dtype=np.float32).reshape(4, 8)
        path = self.root / 'globe.tif'
        write_tiff(path, globe, -180, 90, 45)
        with patch.object(etopo, 'GLOBE_60', path):
            elevation, extent = etopo.read_box(100, 10, -100, 90, 'etopo60')
            self.assertEqual(extent, dict(west=90.0, south=0.0, east=-90.0, north=90.0, pixel_deg=45.0))
            np.testing.assert_array_equal(elevation, globe[0:2, [6, 7, 0, 1]])
            whole, extent = etopo.read_box(-180, -90, 180, 90, 'etopo60')
            np.testing.assert_array_equal(whole, globe)
            self.assertEqual(extent, dict(west=-180.0, south=-90.0, east=180.0, north=90.0, pixel_deg=45.0))
        with self.assertRaises(ValueError):
            etopo.read_box(0, 0, 10, 10, 'etopo30')


class DownloaderTests(unittest.TestCase):
    NAME = 'N30W165'

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = etopo.tile_path(self.NAME, self.root)

    def client(self, body, headers=None):
        response = MagicMock(status_code=200, headers=headers or {})
        response.__enter__.return_value = response
        response.url = etopo.tile_url(self.NAME) + '?token=temporary'
        response.iter_content.return_value = iter([body])
        client = MagicMock()
        client.get.return_value = response
        return client

    def fetch(self, client=None, verify_only=False, budget=1024**2):
        return FETCH.fetch(self.NAME, self.root, budget, 0, verify_only, client or MagicMock())

    def test_download_receipt_and_tampered_tile(self):
        self.assertEqual(self.fetch(self.client(TIFF)), len(TIFF))
        receipt = json.loads(self.path.with_name(self.path.name + '.receipt.json').read_text())
        self.assertEqual(receipt['request'], FETCH.request_for(self.NAME))
        self.assertEqual(receipt['bytes'], len(TIFF))
        self.assertEqual(receipt['final_url'], etopo.tile_url(self.NAME))
        client = MagicMock()
        self.assertEqual(self.fetch(client, verify_only=True), len(TIFF))
        client.get.assert_not_called()
        data = bytearray(self.path.read_bytes())
        data[-1] ^= 1
        self.path.write_bytes(bytes(data))
        with self.assertRaisesRegex(ValueError, 'changed'):
            self.fetch(verify_only=True)

    def test_missing_tile_in_verify_only(self):
        client = MagicMock()
        with self.assertRaises(FileNotFoundError):
            self.fetch(client, verify_only=True)
        client.get.assert_not_called()
        self.assertEqual(FETCH.main(['--window', 'hawaii', '--root', str(self.root), '--verify-only']), 1)

    def test_existing_tile_without_receipt_is_preserved(self):
        self.path.write_bytes(b'original')
        client = MagicMock()
        with self.assertRaisesRegex(ValueError, 'no receipt'):
            self.fetch(client)
        client.get.assert_not_called()
        self.assertEqual(self.path.read_bytes(), b'original')

    def test_budget_refusals(self):
        client = MagicMock()
        with patch.object(FETCH, 'stored_bytes', return_value=1024**2):
            with self.assertRaisesRegex(ValueError, 'budget'):
                self.fetch(client)
        client.get.assert_not_called()
        with self.assertRaisesRegex(ValueError, 'remaining budget'):
            self.fetch(self.client(TIFF, {'Content-Length': str(10**9)}))
        with self.assertRaisesRegex(ValueError, 'size budget'):
            self.fetch(self.client(TIFF), budget=len(TIFF) - 1)
        self.assertFalse(self.path.exists())
        self.assertFalse(list(self.root.glob('*.part')))

    def test_html_is_not_saved_as_a_tile(self):
        with self.assertRaisesRegex(ValueError, 'HTML'):
            self.fetch(self.client(b'<html>Not found</html>'))
        self.assertFalse(self.path.exists())

    def test_catalogue_tiles_and_unknown_window(self):
        windows = etopo.load_catalogue()
        self.assertEqual(FETCH.catalogue_tiles(windows, ['aleutian_arc']), ['N60E165', 'N60W180'])
        with self.assertRaisesRegex(ValueError, 'Unknown window'):
            FETCH.catalogue_tiles(windows, ['atlantis'])


if __name__ == '__main__':
    unittest.main()
