import base64
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/benchmarks"))
import map_rendering
import tectonics_report
from uwtools import runs
from tectonics_report_data import DTYPES, compare_fields, diagnostics, elevation_m, load_bundle, hypsography_shape, largest_connected_patch, near_convergence
from tectonics_report_data import plate_topology, local_relief, terrain_detail, plate_events, coast_distance, coastal_terrain, island_sizes


def fixture(folder):
    manifest = {"schema": "tectonic-snapshot-bundle/v1", "contract_schema_version": 6,
                "width": 4, "height": 2, "sea_level_m": 31043, "plates": [{}, {}],
                "time_myr": 0, "update_count": 0, "files": {}}
    for field, dtype in DTYPES.items():
        data = np.zeros((2, 4), dtype=dtype)
        if field == "heightmap":
            data[:] = [[.05, .1, .5, .99], [1, 1.1, 1.2, 2.5]]
        filename = field + ".bin"
        data.tofile(folder / filename)
        manifest["files"][field] = filename
    for name in ("boundary_segments", "deforming_regions", "tectonic_junctions"):
        manifest["files"][name] = name + ".json"
        (folder / (name + ".json")).write_text("[]")
    (folder / "manifest.json").write_text(json.dumps(manifest))
    return manifest


class TectonicsReportTests(unittest.TestCase):
    def test_coastline_outline_preserves_grid_and_joins_periodic_seams(self):
        land = np.zeros((6, 8), dtype=bool)
        land[:, :2] = True
        land[:, -2:] = True
        rgba = map_rendering.outline(land)
        self.assertEqual(rgba.shape, (6, 8, 4))
        np.testing.assert_array_equal(rgba[..., 3] > 0, np.tile([False, True, False, False, False, False, True, False], (6, 1)))
        self.assertFalse(map_rendering.outline(np.ones((6, 8), dtype=bool))[..., 3].any())

    def test_external_report_assets_preserve_bytes_deduplicate_and_verify(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            source = folder / 'image.png'
            Image.new('RGB', (2, 2), (20, 80, 130)).save(source)
            pixels = source.read_bytes()
            uri = 'data:image/png;base64,' + base64.b64encode(pixels).decode('ascii')
            report = folder / '1-2.html'
            report.write_text('<script>const DATA=' + json.dumps({'images': [uri, uri]}) + ';</script>', encoding='utf-8')
            report.with_suffix('.json').write_text(json.dumps({'runs': [1, 2]}))
            tectonics_report.externalize_report_assets(report)
            html = report.read_text(encoding='utf-8')
            payload, _ = json.JSONDecoder().raw_decode(html, html.index('const DATA=') + len('const DATA='))
            name = hashlib.sha256(pixels).hexdigest() + '.png'
            self.assertEqual(payload['images'], ['assets/' + name] * 2)
            self.assertEqual(list((folder / 'assets').iterdir()), [folder / 'assets' / name])
            self.assertEqual((folder / 'assets' / name).read_bytes(), pixels)
            metadata = json.loads(report.with_suffix('.json').read_text())
            self.assertEqual(metadata['runs'], [1, 2])
            self.assertEqual(metadata['asset_conversion']['images'], [name])
            self.assertEqual(metadata['asset_conversion']['html_sha256'], hashlib.sha256(report.read_bytes()).hexdigest())
            tectonics_report.externalize_report_assets(report)
            self.assertEqual(report.read_text(encoding='utf-8'), html)
            (folder / 'assets' / name).write_bytes(b'corrupted')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                tectonics_report.externalize_report_assets(report)

    def test_islands_join_both_seams_and_ocean_detail_excludes_coastal_steps(self):
        z = np.full((16,32),-4000.)
        z[0,0] = z[-1,0] = z[0,-1] = 100
        z[8,16] = 100
        islands = island_sizes(z)
        self.assertEqual(islands['landmass_count'],2)
        self.assertEqual(islands['islands_2_to_64_cells'],1)
        self.assertEqual(islands['single_cell_islands'],1)
        detail = terrain_detail({'elevation':z, 'plate_id':np.zeros(z.shape)})
        self.assertEqual(detail['interior_ocean_relief_rms_m'],0)
        self.assertLess(detail['interior_ocean_cells'],z.size)

    def test_coast_distance_and_profiles_join_periodic_seams(self):
        land = np.zeros((4, 16), dtype=bool)
        land[:, :4] = True
        land[:, 12:] = True
        distance = coast_distance(land, 100)
        np.testing.assert_array_equal(distance[0], [350,250,150,50,-1,-1,-1,-1,-1,-1,-1,-1,50,150,250,350])
        z = np.where(land, 1000 + distance*10, -4000)
        m = coastal_terrain({'run_scenario': {'world_width_km': 1600}}, {'elevation': z})
        self.assertEqual(m['coastal_highland_cells'], 0)
        self.assertEqual(m['coastal_face_jump_m']['max'], 5500)
        self.assertEqual(m['coastal_land_m']['max'], 1500)
        self.assertGreater(m['profile'][1]['mean_m'], m['profile'][0]['mean_m'])
        self.assertTrue((coast_distance(np.ones((4,8),dtype=bool), 100) == -1).all())

    def test_zero_elevation_coast_matches_the_land_palette_and_active_sample(self):
        z=np.full((16,32),-4000.)
        z[8,8]=0
        z[8,9]=3000
        boundary=np.zeros(z.shape,dtype=np.uint8);boundary[8,7]=1
        m=coastal_terrain({},dict(elevation=z,boundary_type=boundary))
        self.assertEqual(m['coastal_land_m']['min'],0)
        self.assertEqual(m['active_coastal_cells'],2)
        self.assertEqual(m['active_coastal_land_m']['max'],3000)
        self.assertEqual(island_sizes(z)['islands_2_to_64_cells'],1)

    def test_ancestry_validates_split_fusion_and_rejects_reused_identities(self):
        manifest = {"time_myr": 100, "run_scenario": {"plate_count": 2}, "topology_diagnostics": {
            "stable_ids": [3, 2], "events": [
                {"kind": "split", "time_myr": 20, "parent_a": 0, "parent_b": 0, "result": 2, "cells": 100},
                {"kind": "fusion", "time_myr": 50, "parent_a": 0, "parent_b": 1, "result": 3, "cells": 200}]}}
        result = plate_events(manifest)
        self.assertEqual((result['splits'], result['fusions'], result['invalid']), (1, 1, 0))
        manifest['topology_diagnostics']['events'][1]['result'] = 0
        self.assertGreater(plate_events(manifest)['invalid'], 0)

    def test_identity_map_keeps_colour_key_after_runtime_index_compaction(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            manifest = fixture(folder)
            manifest['topology_diagnostics'] = {'stable_ids': [10, 5]}
            (folder/'manifest.json').write_text(json.dumps(manifest))
            _, fields, _ = load_bundle(folder)
            np.testing.assert_array_equal(fields['plate_identity'], np.full((2,4),10))
            manifest['topology_diagnostics']['stable_ids'] = [10, 10]
            (folder/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'persistent plate identities'):
                load_bundle(folder)

    def test_plate_topology_counts_islands_but_joins_periodic_seams(self):
        owners = np.zeros((8, 16), dtype=np.uint32)
        owners[:, 4:12] = 1
        owners[3, 8] = 0
        metric = plate_topology(owners)
        self.assertEqual(metric['components_per_plate'], {'0': 2, '1': 1})
        self.assertEqual(metric['detached_cells'], 1)
        self.assertEqual(metric['detached_components'], 1)
        self.assertEqual(metric['boundary_faces'], 20)
        self.assertEqual(metric['checkerboard_vertices'], 0)
        self.assertEqual(plate_topology(np.array([[0, 1], [1, 0]]))['checkerboard_vertices'], 4)

    def test_detail_diagnostic_is_translation_invariant_and_not_added_to_height(self):
        z = np.tile(np.array([100, 200, 400, 200, 100, 0, 0, 0.]), (8, 1)) + 100
        detail = local_relief(z)
        np.testing.assert_allclose(local_relief(np.roll(z, (3, 2), (0, 1))), np.roll(detail, (3, 2), (0, 1)))
        self.assertAlmostEqual(float(detail.sum()), 0)
        m = terrain_detail({'elevation': z, 'plate_id': np.zeros(z.shape)})
        self.assertEqual(m['interior_land_cells'], z.size)
        self.assertEqual(m['local_relief_rms_m'], m['interior_land_relief_rms_m'])

    def test_disconnected_owners_fail_topology_check(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fixture(folder)
            manifest, fields, _ = load_bundle(folder)
            fields['plate_id'][:] = [[0, 1, 0, 1], [0, 1, 0, 1]]
            result = diagnostics(manifest, fields)
            self.assertEqual(result['checks']['disconnected_plate_regions'], 2)
            self.assertFalse(result['invariants_pass'])

    def test_age_cannot_exceed_seed_prehistory_plus_elapsed_time(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fixture(folder)
            manifest, fields, _ = load_bundle(folder)
            manifest.update(iteration_count=10, delta_time_myr=1.5)
            fields['crust_age_myr'][:] = 15
            self.assertEqual(diagnostics(manifest, fields)['checks']['crust_age_exceeds_elapsed_and_prehistory'], 0)
            fields['crust_age_myr'][0, 0] = 100
            self.assertEqual(diagnostics(manifest, fields)['checks']['crust_age_exceeds_elapsed_and_prehistory'], 1)

    def test_mass_budget_checks_exported_columns_and_missing_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fixture(folder)
            manifest, fields, _ = load_bundle(folder)
            old = dict(fields)
            fields.update(continental_crust_m=np.full((2, 4), 35000, dtype='<f4'),
                          oceanic_crust_m=np.zeros((2, 4), dtype='<f4'),
                          sediment_m=np.zeros((2, 4), dtype='<f4'))
            mass = 35000 * 2800 * 8 * 1e6
            manifest['run_scenario'] = {'world_width_km': 4}
            manifest['column_budget'] = dict(enabled=True, initial_mass_kg=mass, mass_kg=mass,
                generated_ocean_kg=0, subducted_ocean_kg=0, eroded_kg=0, deposited_kg=0,
                relative_mass_residual=0)
            self.assertEqual(diagnostics(manifest, fields)['checks']['column_mass_budget_invalid'], 0)
            for transfer in ('redistributed_continent_kg', 'redistributed_sediment_kg', 'noise_redistributed_kg'):
                for invalid in (-1, float('nan'), float('inf')):
                    manifest['column_budget'][transfer] = invalid
                    self.assertEqual(diagnostics(manifest, fields)['checks']['column_mass_budget_invalid'], 1)
                manifest['column_budget'][transfer] = mass * 2  # Repeated internal transfers can exceed total mass.
                self.assertEqual(diagnostics(manifest, fields)['checks']['column_mass_budget_invalid'], 0)
            fields['continental_crust_m'][0, 0] *= 2
            self.assertEqual(diagnostics(manifest, fields)['checks']['column_mass_budget_invalid'], 1)
            self.assertNotIn('continental_crust_m', compare_fields(old, fields))
            self.assertIn('continental_crust_m', compare_fields(fields, fields))

    def test_force_budget_is_checked_independently_of_elevation(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fixture(folder)
            manifest, fields, _ = load_bundle(folder)
            manifest["motion_diagnostics"] = {"enabled": True, "relative_force_residual": 1e-15,
                "driving_power_w": 1e12, "dissipation_w": 1e12,
                "mean_speed_cm_yr": 2, "max_speed_cm_yr": 3}
            self.assertTrue(diagnostics(manifest, fields)["invariants_pass"])
            manifest["motion_diagnostics"]["dissipation_w"] = 2e12
            self.assertEqual(diagnostics(manifest, fields)["checks"]["force_balance_invalid"], 1)
            manifest["motion_diagnostics"]["dissipation_w"] = float('nan')
            self.assertFalse(diagnostics(manifest, fields)["invariants_pass"])

    def test_highland_connectivity_and_convergence_include_periodic_seams(self):
        mask = np.zeros((10, 20), dtype=bool)
        mask[0, 0] = mask[-1, 0] = mask[0, -1] = True
        mask[5, 10] = True
        self.assertEqual(largest_connected_patch(mask), 3)
        boundary = np.zeros(mask.shape, dtype=np.uint8)
        boundary[0, 0] = 1
        near = near_convergence(boundary)
        self.assertTrue(near[-1, -1])
        self.assertFalse(near[3, 3])
        self.assertFalse(near[5, 10])

    def test_hypsographic_jump_is_located_without_extreme_tail_bias(self):
        q = np.linspace(0, 1, 201)
        z = np.where(q < .375, -10000., -4000.)
        z[0], z[-1] = -30000, 8000
        result = hypsography_shape({"fraction_below": q, "elevation_m": z})
        self.assertEqual(result["max_central_quantile_step_m"], 6000)
        np.testing.assert_allclose(result["max_step_fraction_below"], [.37, .375])
        self.assertEqual(result["elevation_span_35_40_percent_m"], 6000)

    def test_animation_retains_finalization_at_same_model_time(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            states = [{"run": 2, "label": label, "time": 1, "updates": update}
                      for label, update in (("update000001", 1), ("final", 1))]
            for state in states:
                path = root / f"maps/elevation/{state['label']}.png"
                map_rendering.save_rgb(path, np.zeros((2, 4, 3), dtype=np.uint8))
            with patch.object(tectonics_report, "outputs", lambda run: root):
                metadata = tectonics_report.run_animations(2, states,
                    {"elevation": {"title": "Elevation", "units": "m"}})["elevation"]
            with Image.open(root / metadata["path"]) as animation:
                self.assertEqual(animation.n_frames, 2)
                self.assertEqual(metadata["durations_ms"], [120, 1200])

    def test_animation_keeps_order_timing_and_last_frame(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            sources = []
            for index, colour in enumerate(([30, 60, 90], [180, 130, 60])):
                path = folder / f"{index}.png"
                map_rendering.save_rgb(path, np.broadcast_to(colour, (2, 4, 3)))
                sources.append(path)
            output = folder / "animation.png"
            metadata = map_rendering.save_animation(output, sources, ["0 Myr", "5 Myr"], [600, 1200], width=64)
            self.assertEqual(metadata["frames"], 2)
            with Image.open(output) as animation:
                self.assertEqual(animation.info["loop"], 0)
                self.assertEqual(animation.info["duration"], 600)
                animation.seek(1)
                self.assertEqual(animation.info["duration"], 1200)
                self.assertEqual(animation.convert("RGB").getpixel((32, 55)), (180, 130, 60))

    def test_metric_conversion_preserves_above_codec_heights(self):
        values = elevation_m(np.array([.05, .1, 1, 1.5, 2, 2.5], dtype=np.float32), 31043)
        np.testing.assert_allclose(values, [-31043, -31043, 0, 17246, 34492, 51738])
        self.assertEqual(float(elevation_m(np.array([1]), 65535)[0]), 0)

    def test_codec_loss_is_distinct_from_native_invariants(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fixture(folder)
            manifest, fields, _ = load_bundle(folder)
            result = diagnostics(manifest, fields)
            self.assertTrue(result["invariants_pass"])
            self.assertFalse(result["codec_range_preserved"])
            self.assertEqual(result["checks"]["above_uint16_codec"], 1)
            self.assertEqual(result["checks"]["below_native_physical_floor"], 1)
            self.assertEqual(result["land_elevation_m"]["max"], 51738)
            fields["plate_id"][0, 0] = 2
            self.assertFalse(diagnostics(manifest, fields)["invariants_pass"])

    def test_bundle_checks_hash_size_finiteness_and_containment(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            manifest = fixture(folder)
            _, _, hashes = load_bundle(folder)
            with (folder / "heightmap.bin").open("ab") as stream:
                stream.write(b"x")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                load_bundle(folder, hashes)
            with self.assertRaisesRegex(ValueError, "Wrong field size"):
                load_bundle(folder)
            fixture(folder)
            a = np.ones((2, 4), dtype="<f4");a[0, 0] = np.nan
            a.tofile(folder / "heightmap.bin")
            with self.assertRaisesRegex(ValueError, "Nonfinite"):
                load_bundle(folder)
            manifest["files"]["heightmap"] = "../outside.bin"
            (folder / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "escapes"):
                load_bundle(folder)

    def test_differences_use_numeric_values_even_when_preview_saturates(self):
        a = np.array([[10., 20.]])
        b = np.array([[10., 30.]])
        stops = [(0, [0, 0, 0]), (1, [255, 255, 255])]
        np.testing.assert_array_equal(map_rendering.colourize(a, stops), map_rendering.colourize(b, stops))
        diff = compare_fields({"elevation": a}, {"elevation": b})["elevation"]
        self.assertEqual(diff["mae"], 5)
        self.assertEqual(diff["changed_cells"], 1)
        labels = compare_fields({"plate_id": a}, {"plate_id": b})["plate_id"]
        self.assertNotIn("mae", labels)
        with self.assertRaisesRegex(ValueError, "matching"):
            compare_fields({"elevation": a}, {"elevation": b.T})

    def test_categorical_rendering_and_missing_colour(self):
        rgb = map_rendering.categories(np.array([[0, 1, 9]]), {0: [1, 2, 3], 1: [4, 5, 6]})
        np.testing.assert_array_equal(rgb[0], [[1, 2, 3], [4, 5, 6], [0, 0, 0]])
        with self.assertRaises(ValueError):
            map_rendering.colourize(np.array([1]), [(1, [0]*3), (1, [255]*3)])

    def test_two_run_report_has_real_comparison_and_rejects_tampered_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # Two sealed runs: registry records without a run directory, manifests and bundles in the archive.
            registry = root / "runs/registry/tectonics.json"
            registry.parent.mkdir(parents=True)
            registry.write_text(json.dumps({"runs": [{"id": run, "status": "complete"} for run in (1, 2)]}))
            for run in (1, 2):
                bundle = root / f"runs/archive/checkpoints/tectonics/{run}/final"
                bundle.mkdir(parents=True)
                fixture(bundle)
                heights = np.ones((2, 4), dtype="<f4") * (1 + run * .01)
                heights.tofile(bundle / "heightmap.bin")
                _, _, hashes = load_bundle(bundle)
                record = {"id": run, "status": "complete", "recipe": {"seed": 12345}, "replay": {"exact_match": True},
                          "snapshots": [{"label": "final", "directory": "final"}], "field_sha256": {"final": hashes}}
                path = root / f"runs/archive/manifests/tectonics/{run}.json"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(record))
            outputs = root / "runs/reports/tectonics"
            with patch.object(tectonics_report, "ROOT", root), patch.object(runs, "ROOT", root), \
                    patch.object(runs, "ARCHIVE", root / "runs/archive"):
                output = root / "comparison.html"
                states = tectonics_report.build_report([1, 2], output)
                report = json.loads(output.with_suffix(".json").read_text())
                self.assertEqual(len(states), 2)
                self.assertGreater(report["comparisons"]["0:1"]["elevation"]["mae"], 300)
                with patch.object(map_rendering, 'save_animation', side_effect=AssertionError('Unchanged animation rerendered')):
                    tectonics_report.run_animations(1, states, json.loads(tectonics_report.PALETTE_PATH.read_text())['fields'])
                changed_frame=outputs/'1/analysis/maps/elevation/final.png'
                Image.new('RGB',(4,2),(0,0,0)).save(changed_frame)
                with patch.object(map_rendering, 'save_animation', wraps=map_rendering.save_animation) as render:
                    tectonics_report.run_animations(1, states, json.loads(tectonics_report.PALETTE_PATH.read_text())['fields'])
                    self.assertEqual(render.call_count,1)
                saved = np.load(outputs / "2/analysis/fields/elevation/final.npy")
                self.assertGreater(saved.min(), 600)
                self.assertNotIn("__REPORT_DATA__", output.read_text(encoding="utf-8"))
                maps_only = root / "must-not-exist.html"
                with patch.object(tectonics_report, "compare_fields", side_effect=AssertionError("HTML comparisons requested")):
                    tectonics_report.build_report([1, 2], maps_only, maps_only=True)
                self.assertFalse(maps_only.exists())
                self.assertFalse(maps_only.with_suffix(".json").exists())
                self.assertTrue((outputs / "2/analysis/maps/elevation/animation.png").exists())
                self.assertFalse(any((root / "runs/archive").rglob("*.png")))
                with patch.object(tectonics_report, "compare_fields", side_effect=AssertionError("Quadratic comparisons requested")):
                    lazy = tectonics_report.build_report([1, 2], output, external_assets=True)
                metadata = json.loads(output.with_suffix('.json').read_text())
                self.assertEqual(metadata['comparisons'], {})
                asset = lazy[1]['numeric']['elevation']
                self.assertFalse(asset['categorical'])
                self.assertTrue(lazy[1]['numeric']['plate_id']['categorical'])
                np.testing.assert_array_equal(np.fromfile(root / asset['url'], dtype='<f8'), saved.ravel())
                self.assertNotIn('data:image/png', output.read_text(encoding='utf-8'))
                self.assertNotIn('__MAP_COMPARISON_JS__', output.read_text(encoding='utf-8'))
                self.assertNotIn('__MAP_SELECTION_JS__', output.read_text(encoding='utf-8'))
                (bundle / "heightmap.bin").write_bytes(b"tampered")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    tectonics_report.build_report([1, 2], output)


if __name__ == "__main__":
    unittest.main()
