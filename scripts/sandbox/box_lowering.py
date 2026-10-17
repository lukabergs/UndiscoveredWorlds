"""Mean lowering of the global domain's physics over a window's box: procedural minus physics surface.

The question of 2026-10-17 (terrain roadmap T1 notes, nesting): the global domain removes about 1.7 km from
G's mountains, which a window would inherit from it as a parent. For each label: the mean procedural surface,
the mean physics surface and their difference over the land cells (procedural above 0 m) of each domain's
box (`pixels` in configs/sandbox/terrain.json are world cells; the global preview has the same aspect).

  uv run python scripts/sandbox/box_lowering.py OUT.json LABEL1,LABEL2 G,G-coast
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'


def main():
    out, labels, domains = sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(',')
    report = {}
    for label in labels:
        config = json.loads((SANDBOX / label / 'config.json').read_text(encoding='utf-8'))
        procedural = np.load(SANDBOX / label / 'global-elevation-procedural.npy')
        physics = np.load(SANDBOX / label / 'global-elevation-physics.npy')
        rows, columns = procedural.shape
        world_columns = config.get('world_columns', 256)
        cell = columns / world_columns
        report[label] = {}
        for name in domains:
            x0, y0, x1, y1 = config['domains'][name]['pixels']
            box = (slice(int(y0 * cell), int((y1 + 1) * cell)), slice(int(x0 * cell), int((x1 + 1) * cell)))
            land = procedural[box] > 0
            p, q = procedural[box][land], physics[box][land]
            report[label][name] = {'land_cells': int(land.sum()), 'procedural_m': round(float(p.mean()), 1),
                                   'physics_m': round(float(q.mean()), 1), 'lowering_m': round(float((p - q).mean()), 1),
                                   'lowering_p90_m': round(float(np.percentile(p - q, 90)), 1)}
        print(label, report[label])
    Path(out).write_text(json.dumps(report, indent=1) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
