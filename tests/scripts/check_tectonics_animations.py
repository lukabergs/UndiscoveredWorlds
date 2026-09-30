"""Validate a numbered run's APNG timing and final raster against its static maps."""
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/benchmarks"))
from tectonics_experiment import digest


def main(run_id):
    path = ROOT / f"runs/diagnostics/tectonics/{run_id}/animations.json"
    manifest = json.loads(path.read_text())
    for field, record in manifest.items():
        source = ROOT / record["path"]
        assert digest(source) == record["sha256"], f"Animation hash mismatch: {field}"
        with Image.open(source) as animation:
            assert animation.n_frames == record["frames"] == len(record["times_myr"]), field
            assert animation.info["loop"] == 0, field
            elapsed = 0
            for index, expected_duration in enumerate(record["durations_ms"]):
                animation.seek(index)
                assert abs(animation.info["duration"] - expected_duration) < 1, field
                elapsed += animation.info["duration"]
            assert abs(elapsed - record["duration_ms"]) < 1, field
            with Image.open(ROOT / f"runs/maps/tectonics/{field}/final/{run_id}.png") as final:
                expected = final.convert("RGB").resize((record["width"], record["height"] - 40), Image.Resampling.NEAREST)
            np.testing.assert_array_equal(np.array(animation.convert("RGB"))[40:], np.array(expected), err_msg=field)
    print(f"Run {run_id}: verified {len(manifest)} APNGs, frame order/timing/hashes and exact final-map pixels")


if __name__ == "__main__":
    main(int(sys.argv[1]))
