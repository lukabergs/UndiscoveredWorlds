"""Command-line checks of the tectonics recipe of world_pipeline; all worlds use seed 1.

    uv run --offline python tests/scripts/check_tectonics_cli.py out/build/headless/bin/Release/world_pipeline.exe

Rejected options and the snapshot record are checked here; grid, value-range and
repeat checks of a full run belong to the profile configs/evaluation/tectonics.json.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

SNAPSHOT = "tables/geosphere/tectonic_snapshot.json"


def main(executable):
    executable = str(Path(executable).absolute())
    root_directory = Path(__file__).resolve().parents[2]
    base = [executable, "run", "--recipe", "tectonics", "--seed", "1", "--resolution", "32"]

    def run(extra, **options):
        return subprocess.run(base + extra, cwd=root_directory, capture_output=True, text=True, **options)

    def manifest(bundle):
        return json.loads((bundle / "manifest.json").read_text())

    def record(bundle):
        return json.loads((bundle / SNAPSHOT).read_text())

    def field_hashes(bundle):
        return {entry["id"]: entry["hash"] for entry in manifest(bundle)["fields"]}

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)

        # A checkpoint on a half-Myr step, whatever the order of the options.
        checkpoint = "01-geology.plate_tectonics-t0000.500"
        orders = (["--checkpoint-at", ".5", "--delta-time-myr", ".5"], ["--delta-time-myr", ".5", "--checkpoint-at", ".5"])
        for index, order in enumerate(orders):
            run(order + ["--evolution-duration-myr", "2", "--checkpoint-dir", str(root / f"order{index}")], check=True)
        first, second = root / "order0" / checkpoint, root / "order1" / checkpoint
        assert record(first)["time_myr"] == .5 and record(first)["update_count"] == 1
        assert field_hashes(first) == field_hashes(second), "Argument order changed checkpoint"

        for extra in (
            ["--checkpoint-at", ".25", "--delta-time-myr", ".5"], ["--checkpoint-at", "-1"],
            ["--delta-time-myr", "0"],
            ["--evolution-duration-myr", "0"], ["--evolution-duration-myr", "nan"],
            ["--evolution-duration-myr", "-1"], ["--evolution-duration-myr", "1e20"],
            ["--unknown-option", "1"],
            ["--world-width-km", "0"], ["--mantle-flow-speed-cm-yr", "-1"],
            ["--coastal-gpe-strength", "1.01"], ["--boundary-force-ramp-cm-yr", "-1"],
            ["--initial-margin-width-m", "4001"], ["--sediment-load-feedback", "maybe"],
            ["--noise-pass-wavelength-km", "0"], ["--transport-subdivision", "3"],
            ["--deformation-variation", ".6"], ["--subduction-coupling", "1.1"],
            ["--submarine-erosion-ratio", "1.1"], ["--thermal-subsidence-strength", "-1"],
            ["--initial-ocean-depth-m", "65535"], ["--initial-max-height-m", "30000"],
            ["--plume-count", "65"], ["--plume-tail-radius-km", "0"], ["--plume-weakening", "1.5"],
            ["--plume-count", "2", "--plume-weakening", ".5"],  # weakening needs progressive rifting
            ["--mineral-history", "maybe"],
        ):
            options = extra if "--evolution-duration-myr" in extra else extra + ["--evolution-duration-myr", "2"]
            result = run(options + ["--checkpoint-dir", str(root / "rejected")])
            assert result.returncode != 0, f"Invalid option accepted: {extra}"

        recipe = json.loads((root_directory / "configs/tectonics-baseline.json").read_text())
        explicit = [part for key, value in recipe["scenario"].items() for part in ("--" + key, str(value))]
        records = []
        for name, options in (("defaults", []), ("explicit", explicit)):
            output = root / name
            run(options + ["--evolution-duration-myr", "6", "--checkpoint-at", "3",
                           "--checkpoint-dir", str(root / f"{name}-checkpoints"), "--bundle", str(output)], check=True)
            snapshot = record(output)
            records.append(snapshot)
            scenario = snapshot["run_scenario"]
            assert scenario["motion_model"] == "force-balance"
            assert scenario["surface_model"] == "isostatic-columns"
            assert scenario["transport_subdivision"] == recipe["scenario"]["transport-subdivision"]
            assert scenario["margin_processes"] and scenario["boundary_conditioned_terrain"]
            assert snapshot["time_myr"] == 6 and snapshot["update_count"] == 6 and snapshot["label"] == "final"
            assert snapshot["topology_diagnostics"]["enabled"]
            assert snapshot["column_budget"]["enabled"]
            assert abs(snapshot["column_budget"]["relative_mass_residual"]) < 1e-9
            fields = {entry["id"]: entry for entry in manifest(output)["fields"]}
            for field in ("geosphere.continental_crust", "geosphere.oceanic_crust", "geosphere.sediment",
                          "geosphere.tectonic_height", "geosphere.plate_id"):
                assert fields[field]["shape"] == [16, 32], field
            middle = record(root / f"{name}-checkpoints" / "01-geology.plate_tectonics-t0003.000")
            assert middle["time_myr"] == 3 and middle["update_count"] == 3 and middle["label"] == "update000003"
        assert records[0]["run_scenario"] == records[1]["run_scenario"]
        assert field_hashes(root / "defaults") == field_hashes(root / "explicit"), "Explicit recipe values changed fields"

        # Mineral history is recorded with the run and can be switched off.
        budget = json.loads((root / "defaults/tables/geosphere/mineral_budget.json").read_text())
        rows = (root / "defaults/tables/geosphere/mineral_deposits.csv").read_text().splitlines()
        assert budget["schema"] == "mineral-budget/v1"
        assert rows[0].startswith("id,system,origin,fate,")
        listed = budget["inherited"]["count"] + sum(item["count"] for item in budget["formed"].values())
        assert listed > 0 and len(rows) - 1 == listed, "Deposit list and budget disagree"
        run(["--evolution-duration-myr", "6", "--mineral-history", "off", "--bundle", str(root / "no-minerals")], check=True)
        assert len((root / "no-minerals/tables/geosphere/mineral_deposits.csv").read_text().splitlines()) == 1
        assert field_hashes(root / "no-minerals") == field_hashes(root / "defaults"), "Mineral history changed tectonics"
    print("Tectonics CLI passed: shared recipe, native fields, budgets, continuous time, invalid-option rejection, "
          "mineral history")


if __name__ == "__main__":
    main(sys.argv[1])
