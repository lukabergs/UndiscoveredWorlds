"""scripts/evaluation/earth_windows.py and the relief map of maps.py (multiscale plan S1.2): the window's
geometry from a regional manifest, the source at the 1,850 m boundary, the box average, the hillshade's
direction, the rendered window through an injected reader (no refs data, no network) and the entry choice."""

import math
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

import numpy as np  # noqa: E402

from uwtools import palettes  # noqa: E402
from uwtools.bundle import Bundle, FieldInfo  # noqa: E402
import earth_windows  # noqa: E402
from earth_windows import (box_average, choose_entry, global_ocean, map_request, read_window, render_maps,  # noqa: E402
                           render_window, source_for, window_grid, window_sea)
from maps import hillshade, regional_relief, regional_terrain, shade_factor  # noqa: E402


def manifest(columns, rows, west, south, pixel_deg, radius_km=6371.0):
    """A regional bundle manifest: grid, window block in degrees and the planet radius."""
    return {"grid": {"columns": columns, "rows": rows},
            "window": {"west": west, "south": south, "east": west + columns * pixel_deg, "north": south + rows * pixel_deg},
            "parameters": {"world.body.radius_km": radius_km}}


def entry(name="test", centre=(52.0, 32.0), box=(20.0, 0.0, 90.0, 60.0), metres=(1000.0, 100000.0), title="Test window"):
    return {"name": name, "title": title, "centre": list(centre), "box": list(box), "metres_per_pixel": list(metres),
            "features": ["a ridge"]}


class GeometryTests(unittest.TestCase):
    def test_window_has_the_maps_pixel_count_and_extent(self):
        pixel = 22.5 / 256                               # regional run 13: 320 x 256 pixels over 28.125 x 22.5 degrees
        request = map_request(manifest(320, 256, 50.625, 0.0, pixel))
        self.assertAlmostEqual(request["metres_per_pixel"], math.pi * 6371000.0 / 2048, places=6)
        self.assertAlmostEqual(request["latitude"], 11.25)
        grid = window_grid(entry(), request)
        self.assertEqual(grid["pixels"], [320, 256])
        self.assertAlmostEqual(grid["pixel_deg"], pixel, places=12)
        extent = grid["extent"]
        self.assertAlmostEqual(extent["west"], 52.0 - 160 * pixel)
        self.assertAlmostEqual(extent["east"], 52.0 + 160 * pixel)
        self.assertAlmostEqual(extent["north"], 32.0 + 128 * pixel)
        self.assertAlmostEqual(extent["south"], 32.0 - 128 * pixel)

    def test_metres_not_degrees_carry_over_from_another_planet(self):
        request = map_request(manifest(40, 32, 0.0, 0.0, 0.5, radius_km=3185.5))  # half the Earth's radius
        self.assertAlmostEqual(window_grid(entry(), request)["pixel_deg"], 0.25, places=9)

    def test_an_etopo60_window_is_not_cut_to_the_box(self):
        pixel = 22.5 / 256                               # 9,773 m: the 60-arc-second globe
        request = map_request(manifest(320, 256, 50.625, 0.0, pixel))
        grid = window_grid(entry(centre=(52.5, 29.5), box=(47.0, 26.5, 58.0, 34.5)), request)  # Zagros' box
        self.assertEqual(grid["pixels"], [320, 256])
        self.assertAlmostEqual(grid["extent"]["west"], 52.5 - 160 * pixel)

    def test_an_etopo15_window_is_cut_to_the_box(self):
        pixel = 1.0 / 128                                # about 869 m: the 15-arc-second tiles
        request = map_request(manifest(40, 32, 0.0, 0.0, pixel))
        self.assertEqual(source_for(request["metres_per_pixel"]), "etopo15")
        grid = window_grid(entry(centre=(52.0, 32.0), box=(51.875, 31.9375, 52.5, 33.0)), request)
        self.assertEqual(grid["pixels"], [36, 24])       # 16 columns west of the centre, 20 east; 16 rows north, 8 south
        self.assertEqual(grid["extent"], {"west": 51.875, "south": 31.9375, "east": 52.15625, "north": 32.125})

    def test_a_resampled_map_takes_its_pixel_size_and_extent_from_its_pixels(self):
        """A kept box of 8 x 8 degrees in 20 km pixels (0.17986 degrees): 44 x 44 pixels, the partial edge
        pixel dropped, so the pixels end 0.0861 degrees short of the box's south and east edges."""
        resampled = manifest(44, 44, 0.0, 10.0, 8.0 / 44)
        resampled["window"].update({"east": 8.0, "north": 18.0, "pixel_columns": 44, "pixel_rows": 44,
                                    "pixel_cells": 1.11, "resampled": True, "metres_per_pixel": 20000.0})
        request = map_request(resampled)
        pixel = math.degrees(20000.0 / 6371000.0)
        self.assertAlmostEqual(request["metres_per_pixel"], 20000.0, places=6)
        self.assertAlmostEqual(request["pixel_deg"], pixel, places=12)
        self.assertAlmostEqual(request["extent"]["south"], 18.0 - 44 * pixel)
        self.assertAlmostEqual(request["extent"]["east"], 44 * pixel)
        self.assertAlmostEqual(request["latitude"], 18.0 - 22 * pixel)

    def test_across_the_date_line(self):
        request = map_request(manifest(40, 32, 0.0, 0.0, 0.5))
        grid = window_grid(entry(centre=(178.0, 52.0), box=(160.0, 46.0, -160.0, 64.0)), request)
        self.assertEqual(grid["pixels"], [40, 32])
        grid = window_grid(entry(centre=(178.0, 52.0), box=(160.0, 46.0, -160.0, 64.0)), request, clip=True)
        self.assertEqual(grid["pixels"], [40, 28])       # 16 rows north of the centre, 12 south to the box
        self.assertAlmostEqual(grid["extent"]["west"], 168.0)
        self.assertAlmostEqual(grid["extent"]["east"], -172.0)


class SourceTests(unittest.TestCase):
    def test_source_at_the_boundary(self):
        self.assertEqual(source_for(1850.0), "etopo60")
        self.assertEqual(source_for(1849.999), "etopo15")
        self.assertEqual(source_for(9773.0), "etopo60")
        self.assertEqual(source_for(450.0), "etopo15")


class BoxAverageTests(unittest.TestCase):
    def test_aligned_blocks_keep_their_means(self):
        rng = np.random.default_rng(1)
        source = rng.normal(0.0, 1000.0, (8, 12))
        extent = {"west": 10.0, "south": 48.0, "east": 13.0, "north": 50.0, "pixel_deg": 0.25}
        longitudes = 10.0 + 0.5 * np.arange(7)
        latitudes = 50.0 - 0.5 * np.arange(5)
        result = box_average(source, extent, longitudes, latitudes)
        blocks = source.reshape(4, 2, 6, 2).mean(axis=(1, 3))
        np.testing.assert_allclose(result, blocks, rtol=0, atol=1e-9)
        self.assertAlmostEqual(result.mean(), source.mean(), places=9)

    def test_nan_pixels_are_left_out(self):
        source = np.array([[1.0, np.nan], [3.0, 5.0]])
        extent = {"west": 0.0, "north": 1.0, "pixel_deg": 0.5}
        result = box_average(source, extent, [0.0, 1.0], [1.0, 0.0])
        self.assertAlmostEqual(result[0, 0], 3.0)

    def test_a_source_across_the_date_line(self):
        source = np.arange(8.0).reshape(1, 8)            # columns 176..180 and -180..-176 in 1-degree pixels
        extent = {"west": 176.0, "north": 1.0, "pixel_deg": 1.0}
        result = box_average(source, extent, [178.0, 180.0, 182.0], [1.0, 0.0])
        np.testing.assert_allclose(result, [[2.5, 4.5]])

    def test_a_wrapped_source_west_against_edges_below_minus_180(self):
        source = np.arange(8.0).reshape(1, 8)            # columns 175..180 and -180..-177 in 1-degree pixels
        extent = {"west": 175.0, "south": 0.0, "east": -177.0, "north": 1.0, "pixel_deg": 1.0}
        result = box_average(source, extent, [-184.0, -182.0, -180.0, -178.0], [1.0, 0.0])
        np.testing.assert_allclose(result, [[1.5, 3.5, 5.5]])  # -184 is 176 E: the source's second column


class HillshadeTests(unittest.TestCase):
    def test_slopes_facing_the_north_west_light_are_brighter(self):
        rows, columns = 5, 5
        y, x = np.indices((rows, columns))              # x east, y south
        rising_south_east = 100.0 * (x + y)             # faces the north-west
        rising_north_west = -rising_south_east          # faces the south-east
        lit = hillshade(rising_south_east, 1000.0, 1000.0, 1.0)
        shaded = hillshade(rising_north_west, 1000.0, 1000.0, 1.0)
        flat = hillshade(np.zeros((rows, columns)), 1000.0, 1000.0, 1.0)
        self.assertTrue(np.all(lit > flat) and np.all(shaded < flat))
        np.testing.assert_allclose(shade_factor(flat), 1.0)

    def test_exaggeration_strengthens_the_shade(self):
        surface = 30.0 * np.indices((4, 4))[1]          # 300 m over 10 km pixels, rising east
        weak = hillshade(surface, 10000.0, 10000.0, 1.0)
        strong = hillshade(surface, 10000.0, 10000.0, 10.0)
        self.assertGreater(strong[1, 1], weak[1, 1])    # rising east faces west, towards the north-west light


def stripe_reader(calls):
    """A reader of synthetic ETOPO: sea (-100 m) west of 52 E, land rising east of it, pixels of 0.25 degrees
    aligned with the request's edges."""
    def read(west, south, east, north, source):
        calls.append((west, south, east, north, source))
        step = 0.25
        columns, rows = round((east - west) / step), round((north - south) / step)
        longitude = west + step * (np.arange(columns) + 0.5)
        values = np.where(longitude < 52.0, -100.0, 20.0 * (longitude - 52.0) * 100.0)
        return np.tile(values, (rows, 1)).astype(np.float32), {"west": west, "south": south, "east": east,
                                                              "north": north, "pixel_deg": step}
    return read


# An ocean mask of a 1-degree globe that holds every cell.
OCEAN_EVERYWHERE = (np.ones((180, 360), dtype=np.uint8), {"west": -180.0, "north": 90.0, "pixel_deg": 1.0})


def no_files(source, extent):
    """The source files of a synthetic reader: none (the cache records an empty signature)."""
    return []


PROFILE = {"earth": {"maps": ["terrain", "relief"]},
           "maps": [{"name": "terrain", "composite": "regional_terrain", "palette": "tectonics:elevation"},
                    {"name": "relief", "composite": "regional_relief", "options": {"vertical_exaggeration": 10},
                     "palette": "tectonics:elevation"},
                    {"name": "elevation", "derived": "elevation_above_sea", "palette": "tectonics:elevation"}]}


class RenderTests(unittest.TestCase):
    def test_window_renders_at_the_maps_pixels_with_the_palette_and_is_cached(self):
        calls = []
        catalogue = {"schema": "uw-earth-windows/v1", "windows": [entry()]}
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            bundle = Bundle(directory / "bundle", manifest(40, 32, 0.0, 0.0, 0.5))
            written, info = render_maps(bundle, PROFILE, directory / "evaluation", reader=stripe_reader(calls),
                                        catalogue=catalogue, cache_root=directory / "cache", files=no_files, ocean=OCEAN_EVERYWHERE)
            self.assertEqual([item["map"] for item in written], ["terrain", "relief"])
            self.assertTrue(all(item["kind"] == "earth" for item in written))
            self.assertTrue(all((directory / "evaluation" / item["path"]).is_file() for item in written))
            self.assertTrue(all(item["display_fraction"] == 1.0 for item in written))  # as wide as the map
            self.assertEqual(info["source"], "etopo60")
            self.assertEqual(info["pixels"], [40, 32])
            self.assertEqual(written[0]["title"],
                             "Earth: Test window (ETOPO 2022 60 arc-second, 55,597 m/px)")
            cache = directory / "cache" / "test" / "55597m"
            self.assertTrue((cache / "elevation.npy").is_file() and (cache / "manifest.json").is_file())
            self.assertEqual(len(calls), 1)

            request = map_request(bundle.manifest)
            elevation, record = read_window(entry(), request, reader=stripe_reader(calls), cache_root=directory / "cache",
                                            files=no_files)
            self.assertTrue(record["reused"])
            self.assertEqual(len(calls), 1)              # read from the cache
        self.assertEqual(elevation.shape, (32, 40))
        terrain, _ = render_window(elevation, record, "terrain")
        relief, _ = render_window(elevation, record, "relief")
        self.assertEqual(terrain.shape, (32, 40, 3))
        self.assertEqual(relief.shape, (32, 40, 3))
        sea = elevation < 0
        self.assertTrue(sea[:, :20].all() and not sea[:, 20:].any())
        expected = palettes.palette_rgb(np.minimum(elevation, -1.0), "tectonics:elevation")
        np.testing.assert_array_equal(terrain[sea], expected[sea])
        np.testing.assert_array_equal(relief[:, :19], expected[:, :19])  # a flat sea floor keeps its colour
        self.assertTrue(np.all(terrain[:, 20] == 0))     # the landmass outline
        self.assertFalse(np.array_equal(terrain[:, 25:], relief[:, 25:]))  # the slope is shaded
        opposite, _ = render_window(elevation, record, "relief", azimuth_deg=135.0)
        self.assertFalse(np.array_equal(relief[:, 25:], opposite[:, 25:]))  # the light's azimuth is the map's option

    def test_no_entry_is_said(self):
        catalogue = {"schema": "uw-earth-windows/v1", "windows": [entry(metres=(100.0, 1000.0))]}
        with tempfile.TemporaryDirectory() as directory:
            bundle = Bundle(Path(directory), manifest(40, 32, 0.0, 0.0, 0.5))
            written, info = render_maps(bundle, PROFILE, Path(directory), reader=stripe_reader([]), catalogue=catalogue,
                                        cache_root=Path(directory) / "cache", files=no_files, ocean=OCEAN_EVERYWHERE)
        self.assertEqual(written, [])
        self.assertIsNone(info["entry"])
        self.assertIn("55,597 m per pixel", info["note"])

    def test_a_named_entry_outside_its_range_is_drawn_and_said(self):
        catalogue = {"schema": "uw-earth-windows/v1", "windows": [entry(metres=(450.0, 1500.0))]}
        with tempfile.TemporaryDirectory() as directory:
            bundle = Bundle(Path(directory), manifest(40, 32, 0.0, 0.0, 0.5))
            written, info = render_maps(bundle, PROFILE, Path(directory), name="test", reader=stripe_reader([]),
                                        catalogue=catalogue, cache_root=Path(directory) / "cache", files=no_files, ocean=OCEAN_EVERYWHERE)
        self.assertEqual(len(written), 2)
        self.assertTrue(written[0]["title"].endswith("(outside the range where its features show: 450–1,500 m/px)"))

    def test_an_unreadable_window_is_said(self):
        def failing(west, south, east, north, source):
            raise ValueError("pixel grid differs")
        catalogue = {"schema": "uw-earth-windows/v1", "windows": [entry()]}
        with tempfile.TemporaryDirectory() as directory:
            bundle = Bundle(Path(directory), manifest(40, 32, 0.0, 0.0, 0.5))
            written, info = render_maps(bundle, PROFILE, Path(directory), reader=failing, catalogue=catalogue,
                                        cache_root=Path(directory) / "cache", files=no_files, ocean=OCEAN_EVERYWHERE)
        self.assertEqual(written, [])
        self.assertIn("pixel grid differs", info["note"])


class SeaTests(unittest.TestCase):
    """The sea of a window is the world ocean of a global mask: below 0 m and 4-connected, across the date
    line, to the seed at 0 N 160 W, on a synthetic 1-degree globe."""

    def setUp(self):
        below = np.zeros((180, 360), dtype=bool)
        below[70:110, 0:30] = True                       # ocean 180 W to 150 W, 20 S to 20 N, holding the seed
        below[70:110, 350:360] = True                    # and 170 E to 180, joined only across the date line
        below[45:55, 225:235] = True                     # a basin below 0 m at 45 E to 55 E, 35 N to 45 N
        self.extent = {"west": -180.0, "north": 90.0, "pixel_deg": 1.0}
        self.directory = tempfile.TemporaryDirectory()
        self.ocean = global_ocean(self.directory.name, globe=lambda: (below, self.extent), files=lambda: [])

    def tearDown(self):
        self.directory.cleanup()

    def window(self, centre):
        request = map_request(manifest(40, 32, 0.0, 0.0, 0.5))
        return window_grid(entry(centre=centre, box=(centre[0] - 5, centre[1] - 5, centre[0] + 5, centre[1] + 5)),
                           request)

    def test_the_mask(self):
        mask, extent = self.ocean
        self.assertEqual(extent, self.extent)
        self.assertTrue(mask[90, 20] and mask[90, 355])  # the seed's part and the part across the date line
        self.assertFalse(mask[50, 230])                  # the basin
        def failing():
            raise AssertionError("the cached mask is reused")
        cached, _ = global_ocean(self.directory.name, globe=failing, files=lambda: [])
        np.testing.assert_array_equal(cached, mask)

    def test_a_basin_touching_the_windows_edge_stays_land(self):
        record = self.window((50.0, 40.0))
        elevation = np.full((32, 40), -28.0)             # below 0 m everywhere, at every edge of the window
        self.assertFalse(window_sea(elevation, record, self.ocean).any())
        rgb, note = render_window(elevation, record, "terrain", ocean=self.ocean)
        land = palettes.palette_rgb(np.zeros((1, 1)), "tectonics:elevation")[0, 0]
        np.testing.assert_array_equal(rgb[16, 20], land)  # the basin's floor on the land colours (0 m)
        self.assertIn("world ocean", note)

    def test_the_ocean_across_the_date_line_stays_sea(self):
        record = self.window((178.0, 0.0))              # 168 E to 172 W
        sea = window_sea(np.full((32, 40), -100.0), record, self.ocean)
        self.assertFalse(sea[:, :4].any())              # 168 E to 170 E: not ocean on the globe
        self.assertTrue(sea[:, 4:].all())               # 170 E to 180 and 180 to 172 W
        self.assertFalse(window_sea(np.full((32, 40), 10.0), record, self.ocean).any())  # land above 0 m


class CacheTests(unittest.TestCase):
    """A cached window is reused only on the same grid, covering the request, at the same CACHE_VERSION and
    from source files of the same size and modification time."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.source = self.root / "globe.tif"
        self.source.write_bytes(b"0123")
        self.calls = []
        self.request = map_request(manifest(40, 32, 0.0, 0.0, 0.5))
        read_window(entry(), self.request, stripe_reader(self.calls), self.root / "cache", self.files)

    def tearDown(self):
        self.directory.cleanup()

    def files(self, source, extent):
        return [self.source]

    def reused(self, item=None, request=None):
        _, record = read_window(item or entry(), request or self.request, stripe_reader(self.calls), self.root / "cache",
                                self.files)
        return record["reused"]

    def test_the_same_request_is_reused(self):
        self.assertTrue(self.reused())
        self.assertTrue(self.reused(request=map_request(manifest(20, 16, 0.0, 0.0, 0.5))))  # a nested, smaller window

    def test_a_different_centre_is_not_reused(self):
        self.assertFalse(self.reused(item=entry(centre=(52.25, 32.0))))

    def test_a_wider_request_is_not_reused(self):
        self.assertFalse(self.reused(request=map_request(manifest(48, 32, 0.0, 0.0, 0.5))))

    def test_a_changed_version_is_not_reused(self):
        version = earth_windows.CACHE_VERSION
        earth_windows.CACHE_VERSION = version + 1
        try:
            self.assertFalse(self.reused())
        finally:
            earth_windows.CACHE_VERSION = version

    def test_a_changed_source_file_is_not_reused(self):
        self.source.write_bytes(b"012345")
        self.assertFalse(self.reused())


class EntryChoiceTests(unittest.TestCase):
    CATALOGUE = {"schema": "uw-earth-windows/v1", "windows": [
        entry("north", centre=(10.0, 60.0), metres=(1000.0, 20000.0)),
        entry("south", centre=(10.0, 10.0), metres=(1000.0, 20000.0)),
        entry("coarse", centre=(10.0, 30.0), metres=(50000.0, 100000.0))]}

    def test_nearest_latitude_among_the_entries_holding_the_pixel_size(self):
        self.assertEqual(choose_entry(self.CATALOGUE, 5000.0, 25.0)["name"], "south")
        self.assertEqual(choose_entry(self.CATALOGUE, 5000.0, 50.0)["name"], "north")
        self.assertEqual(choose_entry(self.CATALOGUE, 20000.0, 35.0)["name"], "north")   # ranges are closed; tie: first
        self.assertEqual(choose_entry(self.CATALOGUE, 60000.0, -80.0)["name"], "coarse")
        self.assertIsNone(choose_entry(self.CATALOGUE, 30000.0, 30.0))

    def test_a_name_overrides_the_rule(self):
        self.assertEqual(choose_entry(self.CATALOGUE, 5000.0, 60.0, "coarse")["name"], "coarse")
        with self.assertRaises(ValueError):
            choose_entry(self.CATALOGUE, 5000.0, 60.0, "missing")


def regional_bundle(directory, elevation):
    rows, columns = elevation.shape
    values = {"geosphere.base_elevation": elevation.astype(np.float32),
              "geosphere.sea_level": np.array(0.0, dtype=np.float32),
              "geosphere.sea": (elevation < 0).astype(np.uint8),
              "cryosphere.ice_thickness": np.zeros((rows, columns), dtype=np.float32),
              "hydrosphere.lake_surface": np.zeros((rows, columns), dtype=np.float32),
              "hydrosphere.channel_direction": np.zeros((rows, columns), dtype=np.uint8),
              "hydrosphere.discharge": np.zeros((rows, columns), dtype=np.float32),
              "hydrosphere.river_flow_january": np.zeros((rows, columns), dtype=np.int32),
              "hydrosphere.river_flow_july": np.zeros((rows, columns), dtype=np.int32)}
    result = Bundle(Path(directory), manifest(columns, rows, 0.0, 0.0, 0.1))
    for identifier, array in values.items():
        np.save(Path(directory) / f"{identifier}.npy", array)
        result.fields[identifier] = FieldInfo(identifier, identifier.split(".")[0], "grid", str(array.dtype),
                                              array.shape, "", f"{identifier}.npy", "", "", "")
    return result


class RegionalReliefTests(unittest.TestCase):
    def test_flat_ground_keeps_the_terrain_colours_and_slopes_are_shaded(self):
        elevation = np.full((6, 8), 200.0)
        elevation[:, :3] = -50.0
        elevation[3:, 5:] = 200.0 + 300.0 * np.indices((3, 3))[1]
        with tempfile.TemporaryDirectory() as directory:
            item = regional_bundle(directory, elevation)
            terrain, _ = regional_terrain(item)
            relief, note = regional_relief(item, vertical_exaggeration=10)
        np.testing.assert_array_equal(relief[:, 0], terrain[:, 0])        # flat sea
        np.testing.assert_array_equal(relief[0, 4], terrain[0, 4])        # flat land
        self.assertFalse(np.array_equal(relief[4, 6], terrain[4, 6]))     # the slope
        self.assertNotIn("hillshade", note)                               # the title states it; the legend stays short

    def test_the_light_azimuth_option(self):
        """S3.7: a slope rising east faces west: brighter lit from the south-west (225) than from the north-east
        (45); flat land keeps its colour under every light; the default is the north-west light (315)."""
        elevation = np.full((6, 8), 200.0)
        elevation[:, :3] = -50.0
        elevation[3:, 5:] = 200.0 + 300.0 * np.indices((3, 3))[1]
        with tempfile.TemporaryDirectory() as directory:
            item = regional_bundle(directory, elevation)
            default, _ = regional_relief(item, vertical_exaggeration=10)
            north_west, _ = regional_relief(item, vertical_exaggeration=10, azimuth_deg=315)
            north_east, _ = regional_relief(item, vertical_exaggeration=10, azimuth_deg=45)
            south_west, _ = regional_relief(item, vertical_exaggeration=10, azimuth_deg=225)
        np.testing.assert_array_equal(default, north_west)
        self.assertGreater(int(south_west[4, 6].sum()), int(north_east[4, 6].sum()))  # the west-facing slope
        np.testing.assert_array_equal(north_east[0, 4], default[0, 4])                 # flat land keeps its colour
        np.testing.assert_array_equal(south_west[0, 4], default[0, 4])


if __name__ == "__main__":
    unittest.main()
