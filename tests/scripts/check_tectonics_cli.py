"""Command-line checks of world_pipeline's tectonics recipe that no bundle or unit test can make; seed 1.

    uv run --offline python tests/scripts/check_tectonics_cli.py out/build/headless/bin/Release/world_pipeline.exe

Here: options the command line refuses before the engine starts, checkpoint times that
do not depend on the order of the options, checkpoint labels and update counts, and
--mineral-history off leaving the tectonic fields unchanged. Elsewhere: the values the
scenario options or the engine refuse (tectonics_tests, TectonicOptions.*); recipe
equality, the column and mineral budgets and the snapshot record (checks of
configs/evaluation/tectonics.json, `uw.py tectonics experiment --repeat --explicit-recipe`).
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

    def record(bundle):
        return json.loads((bundle / SNAPSHOT).read_text())

    def field_hashes(bundle):
        return {entry["id"]: entry["hash"] for entry in json.loads((bundle / "manifest.json").read_text())["fields"]}

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

        # Refused by the command line (one engine refusal stands for the rest).
        for extra in (
            ["--checkpoint-at", ".25", "--delta-time-myr", ".5"], ["--checkpoint-at", "-1"],
            ["--unknown-option", "1"], ["--initial-max-height-m", "30000"], ["--mineral-history", "maybe"],
            ["--transport-subdivision", "3"],
        ):
            result = run(extra + ["--evolution-duration-myr", "2", "--checkpoint-dir", str(root / "rejected")])
            assert result.returncode != 0, f"Invalid option accepted: {extra}"

        defaults = root / "defaults"
        run(["--evolution-duration-myr", "6", "--checkpoint-at", "3", "--checkpoint-dir", str(root / "checkpoints"),
             "--bundle", str(defaults)], check=True)
        final = record(defaults)
        assert final["time_myr"] == 6 and final["update_count"] == 6 and final["label"] == "final"
        middle = record(root / "checkpoints" / "01-geology.plate_tectonics-t0003.000")
        assert middle["time_myr"] == 3 and middle["update_count"] == 3 and middle["label"] == "update000003"

        # Mineral history is recorded with the run and can be switched off without changing tectonics.
        run(["--evolution-duration-myr", "6", "--mineral-history", "off", "--bundle", str(root / "no-minerals")], check=True)
        assert len((root / "no-minerals/tables/geosphere/mineral_deposits.csv").read_text().splitlines()) == 1
        assert len((defaults / "tables/geosphere/mineral_deposits.csv").read_text().splitlines()) > 1
        assert field_hashes(root / "no-minerals") == field_hashes(defaults), "Mineral history changed tectonics"
    print("Tectonics CLI passed: refused options, checkpoint times and labels, mineral history switch")


if __name__ == "__main__":
    main(sys.argv[1])
