"""Public CLI checks for the sole current tectonic engine; all worlds use seed 1."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main(executable):
    executable = str(Path(executable).absolute())
    base = [executable, "inspect", "--seed", "1", "--width", "16", "--height", "8"]

    def invoke(extra):
        return subprocess.run(base + extra, capture_output=True, text=True, check=True)

    first = invoke(["--time-myr", ".5", "--delta-time-myr", ".5"])
    second = invoke(["--delta-time-myr", ".5", "--time-myr", ".5"])
    assert json.loads(first.stdout) == json.loads(second.stdout), "Argument order changed checkpoint"
    for extra in (
        ["--time-myr", ".25", "--delta-time-myr", ".5"],
        ["--delta-time-myr", "0"], ["--time-myr", "1e20"], ["--time-myr", "-1"],
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
    ):
        result = subprocess.run(base + extra, capture_output=True, text=True)
        assert result.returncode != 0, f"Invalid option accepted: {extra}"

    recipe = json.loads((Path(__file__).resolve().parents[2] / "configs/tectonics-baseline.json").read_text())
    explicit = [part for key, value in recipe["scenario"].items() for part in ("--" + key, str(value))]
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        manifests = []
        for name, options in (("defaults", []), ("explicit", explicit)):
            output = root / name
            invoke(options + ["--evolution-duration-myr", "6", "--final", "--bundle-output", str(output)])
            manifest = json.loads((output / "manifest.json").read_text())
            manifests.append(manifest)
            scenario = manifest["run_scenario"]
            assert scenario["motion_model"] == "force-balance"
            assert scenario["surface_model"] == "isostatic-columns"
            assert scenario["transport_subdivision"] == recipe["scenario"]["transport-subdivision"]
            assert scenario["margin_processes"] and scenario["boundary_conditioned_terrain"]
            assert manifest["time_myr"] == 6
            assert manifest["topology_diagnostics"]["enabled"]
            assert manifest["column_budget"]["enabled"]
            assert abs(manifest["column_budget"]["relative_mass_residual"]) < 1e-9
            for field in ("continental_crust_m", "oceanic_crust_m", "sediment_m"):
                assert (output / manifest["files"][field]).stat().st_size == 16 * 8 * 4
        assert manifests[0]["run_scenario"] == manifests[1]["run_scenario"]
        for field, filename in manifests[0]["files"].items():
            assert (root / "defaults" / filename).read_bytes() == (
                root / "explicit" / manifests[1]["files"][field]).read_bytes(), field

        output = root / "batch"
        subprocess.run([executable, "batch-export", "--seed", "1", "--width", "16", "--height", "8",
                        "--evolution-duration-myr", "6", "--output", str(output),
                        "--checkpoint-at-time-myr", "3", "--checkpoint-at-update", "3"],
                       capture_output=True, text=True, check=True)
        batch = json.loads((output / "batch_manifest.json").read_text())
        for label in ("time3.000Myr", "update000003"):
            item = next(s for s in batch["snapshots"] if s["label"] == label)
            assert item["time_myr"] == 3 and item["update_count"] == 3
        assert not (output / "mineral_budget.json").exists(), "Mineral history must be opt-in"

        output = root / "minerals"
        small = ["batch-export", "--seed", "1", "--width", "16", "--height", "8", "--evolution-duration-myr", "6"]
        subprocess.run([executable, *small, "--output", str(output), "--mineral-history", "on"],
                       capture_output=True, text=True, check=True)
        budget = json.loads((output / "mineral_budget.json").read_text())
        rows = (output / "mineral_deposits.csv").read_text().splitlines()
        assert budget["schema"] == "mineral-budget/v1"
        assert rows[0].startswith("id,system,origin,fate,")
        listed = budget["inherited"]["count"] + sum(item["count"] for item in budget["formed"].values())
        assert listed > 0 and len(rows) - 1 == listed, "Deposit list and budget disagree"
        rejected = subprocess.run([executable, *small, "--output", str(root / "rejected"), "--mineral-history", "maybe"],
                                  capture_output=True, text=True)
        assert rejected.returncode != 0, "Invalid --mineral-history value accepted"
    print("Tectonics CLI passed: shared recipe, native fields, budgets, continuous time, invalid-option rejection, "
          "mineral history")


if __name__ == "__main__":
    main(sys.argv[1])
