"""Native-grid ingestion and reference-source download contracts."""

import importlib.util
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

import numpy as np

try:
    import xarray as xr
except ImportError:
    xr = None

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts/refs"
sys.path.insert(0, str(SCRIPTS))


def module(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@unittest.skipIf(xr is None, "Run with --with xarray --with netcdf4 for source ingestion tests")
class AdditionalReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ingest = module("prepare-additional-references")

    def data(self, variable="sst", units="degC"):
        # Native pole-inclusive source, reversed latitude and shifted longitude.
        dates = np.arange("2001-01", "2003-01", dtype="datetime64[M]")
        values = np.broadcast_to(np.arange(24, dtype=np.float32)[:, None, None], (24, 3, 4)).copy()
        return xr.DataArray(values, dims=("valid_time", "latitude", "longitude"),
                            coords=dict(valid_time=dates, latitude=[-90., 0., 90.], longitude=[0., 90., 180., 270.]),
                            name=variable, attrs=dict(units=units))

    def test_native_poles_are_remapped_not_cropped(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, output = root / "source", root / "output"
            source.mkdir()
            self.data().to_dataset().to_netcdf(source / "sst.nc")
            self.ingest.prepare(source, output, "oisst", 8, 2001, 2002)
            bundle = output / "oisst_sst_monthly.uwclim"
            raw = bundle.read_bytes()
            self.assertEqual(struct.unpack("<IIII", raw[8:24]), (1, 8, 4, 12))
            values = np.frombuffer(raw[40:], dtype="<f4").reshape(12, 4, 8)
            for month in range(12):
                np.testing.assert_array_equal(values[month], month+6)
            receipt = json.loads((output / "oisst-additional-preparation.json").read_text())
            self.assertEqual(receipt["height"], 4)
            self.assertEqual(receipt["products"][0]["units"], "degrees C")
            import reference_climate_maps
            maps = module("prepare-reduced-earth-benchmark")
            emitted = []
            land = np.zeros((2, 4), dtype=bool)
            land[:, 0] = True
            names = reference_climate_maps.export_additional(
                output, 4, land, maps.load_uwclim,
                lambda *args, **kwargs: emitted.append((args, kwargs)))
            self.assertEqual(names, {"oisst_sst"})
            self.assertEqual(len(emitted), 4)
            for args, metadata in emitted:
                self.assertEqual(args[1].shape, (2, 4))
                self.assertTrue(np.isnan(args[1][:, 0]).all())
                self.assertTrue(np.isfinite(args[1][:, 1:]).all())
                self.assertIn("day-weighted", metadata["period"])

    def test_temporal_gaps_duplicates_and_units_are_rejected(self):
        data = self.data()
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            list(self.ingest.monthly_climatology([data.isel(valid_time=slice(1, None))], 2001, 2002))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            list(self.ingest.monthly_climatology([data, data], 2001, 2002))
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.data(units="K").to_dataset().to_netcdf(root / "sst.nc")
            with self.assertRaisesRegex(ValueError, "Unexpected units"):
                self.ingest.prepare(root, root / "output", "oisst", 8, 2001, 2002)
        with self.assertRaises(ValueError):
            self.ingest.grid_dimensions(2048, 1025)

    def test_downloaded_group_cannot_silently_omit_a_requested_variable(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.data(variable="viwve", units="kg m**-1 s**-1").to_dataset().to_netcdf(root / "moisture.nc")
            (root / "moisture-download.json").write_text(json.dumps(dict(
                status="downloaded", output=dict(file="moisture.nc"))))
            with self.assertRaisesRegex(ValueError, "lacks expected variables"):
                self.ingest.prepare(root, root / "output", "era5", 8, 2001, 2002)

    def test_divergence_conversion_and_missing_cells(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = self.data(variable="vimdf", units="kg m**-2 s**-1")
            data.values[:] = 1e-5
            data.values[:, 1, 1] = np.nan
            data.to_dataset().to_netcdf(root / "flux.nc")
            self.ingest.prepare(root, root / "output", "era5", 8, 2001, 2002)
            raw = (root / "output/era5_vimd_monthly.uwclim").read_bytes()
            values = np.frombuffer(raw[40:], dtype="<f4")
            np.testing.assert_allclose(values[np.isfinite(values)], -.864, rtol=1e-6)
            self.assertTrue(np.isnan(values).any())

    def test_snow_depth_is_not_converted_to_water_equivalent(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = self.data(variable="sde", units="m")
            data.values[:] = 2
            data.to_dataset().to_netcdf(root / "snow.nc")
            self.ingest.prepare(root, root / "output", "land", 8, 2001, 2002)
            raw = (root / "output/era5land_snow_depth_monthly.uwclim").read_bytes()
            np.testing.assert_array_equal(np.frombuffer(raw[40:], dtype="<f4"), 2)
            receipt = json.loads((root / "output/land-additional-preparation.json").read_text())
            self.assertEqual(receipt["products"][0]["units"], "m snow depth")
            self.assertTrue(receipt["products"][0]["land_only"])
            import reference_climate_maps
            maps = module("prepare-reduced-earth-benchmark")
            land = np.zeros((2, 4), dtype=bool)
            land[:, 0] = True
            emitted = []
            reference_climate_maps.export_additional(root / "output", 4, land, maps.load_uwclim,
                lambda name, values, *args, **kwargs: emitted.append(values))
            self.assertEqual(len(emitted), 4)
            for values in emitted:
                np.testing.assert_array_equal(values[:, 0], 2)
                self.assertTrue(np.isnan(values[:, 1:]).all())

    def test_modern_era5_names_and_ocean_surface_depth(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = self.data(variable="avg_ie", units="kg m**-2 s**-1")
            data.values[:] = -1e-5
            data.to_dataset().to_netcdf(root / "energy.nc")
            self.ingest.prepare(root, root / "output", "era5", 8, 2001, 2002)
            raw = (root / "output/era5_mer_monthly.uwclim").read_bytes()
            np.testing.assert_allclose(np.frombuffer(raw[40:], dtype="<f4"), .864, rtol=1e-6)
            ocean = self.data(variable="uo", units="m s-1").expand_dims(depth=[.494])
            self.assertEqual(self.ingest.normalize(ocean).shape, (24, 3, 4))
            with self.assertRaisesRegex(ValueError, "Ambiguous dimension depth"):
                self.ingest.normalize(xr.concat([ocean, ocean.assign_coords(depth=[1.0])], "depth"))

    def test_ocean_vector_exports_physical_speed_and_masks_land(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            u = self.data(variable="uo", units="m s-1").expand_dims(depth=[.494])
            v = self.data(variable="vo", units="m s-1").expand_dims(depth=[.494])
            u = xr.full_like(u, 3)
            v = xr.full_like(v, 4)
            xr.Dataset(dict(uo=u, vo=v)).to_netcdf(root / "ocean.nc")
            self.ingest.prepare(root, root / "output", "marine", 8, 2001, 2002)
            maps = module("prepare-reduced-earth-benchmark")
            import reference_climate_maps
            emitted = {}
            land = np.zeros((2, 4), dtype=bool)
            land[:, 0] = True
            available = reference_climate_maps.export_additional(root / "output", 4, land, maps.load_uwclim,
                        lambda name, values, *args, **kwargs: emitted.update({name: values}))
            self.assertIn("ocean_currents", available)
            self.assertEqual(len(emitted), 12)
            self.assertTrue(np.isnan(emitted["jan_ocean_current_speed"][:, 0]).all())
            np.testing.assert_array_equal(emitted["jan_ocean_current_speed"][:, 1:], 5)


if __name__ == "__main__":
    unittest.main()
