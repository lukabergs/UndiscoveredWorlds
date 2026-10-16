"""Striping evidence (not a fix): on the global domain's land, the median |dz/dx| over the median |dz/dy| (true km,
dx = dy cos latitude, central differences) by latitude band; Earth's ETOPO 2048 about 0.9 to 1.2 (roadmap v17).

  uv run python scripts/sandbox/striping_ratio.py LABEL1,LABEL2
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path('C:/dev/undiscovered-worlds')
sys.path.insert(0, str(ROOT / 'scripts' / 'sandbox'))
import terrain_metrics as tm  # noqa: E402

for label in sys.argv[1].split(','):
    d = tm.SANDBOX / label
    config = json.loads((d / 'config.json').read_text())
    z = np.load(d / 'global-elevation.npy').astype(float)
    phys = np.load(d / 'global-elevation-physics.npy').astype(float)
    pre = np.load(d / 'global-elevation-procedural.npy').astype(float)
    area, lat, dy = tm.geometry(config, 'global', z.shape)
    land = ~tm.sea_from(phys, seed=tm.sea_from(pre)) & (z > 0)
    gy = np.abs(np.gradient(z, axis=0)) / dy
    gx = np.abs((np.roll(z, -1, 1) - np.roll(z, 1, 1)) / 2.0) / (dy * np.cos(np.radians(lat))[:, None])
    row = []
    for lo, hi in ((0, 30), (30, 50), (50, 60), (60, 70), (70, 80)):
        band = land & (np.abs(lat)[:, None] >= lo) & (np.abs(lat)[:, None] < hi)
        if band.sum() > 100:
            row.append(f'{lo}-{hi}: {np.median(gx[band]) / max(np.median(gy[band]), 1e-9):.2f}')
    print(label, '; '.join(row))
