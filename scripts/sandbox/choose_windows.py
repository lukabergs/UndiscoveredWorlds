"""Candidate window boxes on a sandbox label's global render, by the roles the default windows play.

The windows of configs/sandbox/terrain.json are pixel boxes in world cells, picked on one base world's geography
(base-v2); a new base needs new boxes. Roles (each a 20 x 16 world-cell box at k 64, nested 4 x 4 boxes at k 256):
  B  the highest orogen (mean land height, land >= 70 %), compared with the Himalaya;   A: B's box 16 cells west
  G  a glaciated coastal range at 39-62 N (land 40-75 %, iced land, high land beside the sea), against Norway
  F  the same at 42-65 S, against Patagonia;   O: high-latitude land and coast at 51-75 N
  B-core: B's highest 4 x 4; B-foreland: the 4 x 4 beside it with the steepest fall; G-coast, F-troughs: the 4 x 4
  inside G and F with the most ice-cut coast (iced land share x land height, land 60-95 %, world ocean 10-40 %).
Scores are printed with each box so the choice can be checked; the global wraps east-west.

  uv run python scripts/sandbox/choose_windows.py LABEL [--json OUT]
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'
WORLD_COLUMNS, WORLD_ROWS = 256, 128


def latitude(row):
    return 90.0 - (row + 0.5) * 180.0 / WORLD_ROWS


def block_means(field, cell):
    rows, columns = field.shape[0] // cell, field.shape[1] // cell
    return field[:rows * cell, :columns * cell].reshape(rows, cell, columns, cell).mean(axis=(1, 3))


def box_stat(grid, x0, y0, w, h):
    columns = [(x0 + i) % grid.shape[1] for i in range(w)]
    return grid[y0:y0 + h][:, columns]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('label')
    parser.add_argument('--json')
    args = parser.parse_args()
    run = SANDBOX / args.label
    elevation = np.load(run / 'global-elevation.npy')
    ice = np.load(run / 'global-ice_time_myr.npy')
    cell = elevation.shape[1] // WORLD_COLUMNS
    land = block_means((elevation > 0).astype(float), cell)
    height = block_means(np.where(elevation > 0, elevation, 0.0), cell) / np.maximum(land, 1e-6)
    iced = block_means(((ice > 0) & (elevation > 0)).astype(float), cell) / np.maximum(land, 1e-6)
    # The world ocean: the largest connected water body (lakes and closed basins below 0 m are not sea).
    labels, _ = ndimage.label(elevation <= 0)
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    ocean = block_means((labels == sizes.argmax()).astype(float), cell)
    result = {}

    def best(w, h, score, rows=None, columns=None):
        top = None
        for y in (rows if rows is not None else range(0, WORLD_ROWS - h + 1)):
            for x in (columns if columns is not None else range(WORLD_COLUMNS)):
                value = score(x, y)
                if value is not None and (top is None or value > top[0]):
                    top = (value, x, y)
        return top

    def big(x, y, need_land, lat_range, kind):
        lf = box_stat(land, x, y, 20, 16).mean()
        if not (need_land[0] <= lf <= need_land[1]):
            return None
        lat = latitude(y + 8)
        if lat_range and not (lat_range[0] <= lat <= lat_range[1]):
            return None
        lw = box_stat(land, x, y, 20, 16)
        mean_height = (box_stat(height, x, y, 20, 16) * lw).sum() / max(lw.sum(), 1e-6)
        if kind == 'orogen':
            return mean_height
        icy = (box_stat(iced, x, y, 20, 16) * lw).sum() / max(lw.sum(), 1e-6)
        return mean_height * (0.2 + icy) if kind == 'glacial' else lf * (1 - lf) * (0.2 + icy)

    picks = {
        'B': best(20, 16, lambda x, y: big(x, y, (0.7, 1.0), None, 'orogen')),
        'G': best(20, 16, lambda x, y: big(x, y, (0.4, 0.75), (39, 62), 'glacial')),
        'F': best(20, 16, lambda x, y: big(x, y, (0.4, 0.75), (-65, -42), 'glacial')),
    }
    gx, gy = picks['G'][1], picks['G'][2]
    apart = lambda x, y: min((x - gx) % WORLD_COLUMNS, (gx - x) % WORLD_COLUMNS) >= 20 or abs(y - gy) >= 16
    picks['O'] = best(20, 16, lambda x, y: big(x, y, (0.3, 0.7), (51, 75), 'coast') if apart(x, y) else None)
    for name, top in picks.items():
        value, x, y = top
        result[name] = {'pixels': [x, y, (x + 19) % WORLD_COLUMNS, y + 15], 'score': round(float(value), 1),
                        'lat': [round(latitude(y + 15), 1), round(latitude(y), 1)]}
    x, y = picks['B'][1], picks['B'][2]
    result['A'] = {'pixels': [(x - 16) % WORLD_COLUMNS, y, (x + 3) % WORLD_COLUMNS, y + 15], 'note': "B's box 16 cells west"}

    def inner(parent, score):
        px, py = picks[parent][1], picks[parent][2]
        top = best(4, 4, score, rows=range(py, py + 13), columns=[(px + i) % WORLD_COLUMNS for i in range(17)])
        value, x, y = top
        return {'pixels': [x, y, (x + 3) % WORLD_COLUMNS, y + 3], 'score': round(float(value), 2), 'parent': parent}

    def small(x, y, kind):
        lw = box_stat(land, x, y, 4, 4)
        lf = lw.mean()
        hh = (box_stat(height, x, y, 4, 4) * lw).sum() / max(lw.sum(), 1e-6)
        if kind == 'core':
            return hh if lf >= 0.9 else None
        if not (0.6 <= lf <= 0.95) or not (0.1 <= box_stat(ocean, x, y, 4, 4).mean() <= 0.4):
            return None
        return hh * (box_stat(iced, x, y, 4, 4) * lw).sum() / max(lw.sum(), 1e-6)

    result['B-core'] = inner('B', lambda x, y: small(x, y, 'core'))
    core = result['B-core']['pixels']

    def foreland(x, y):
        if abs(x - core[0]) > 6 or abs(y - core[1]) > 6 or (abs(x - core[0]) < 4 and abs(y - core[1]) < 4):
            return None
        lw = box_stat(land, x, y, 4, 4)
        if lw.mean() < 0.9:
            return None
        h = box_stat(height, x, y, 4, 4)
        return float(h.max() - h.min())

    result['B-foreland'] = inner('B', foreland)
    result['G-coast'] = inner('G', lambda x, y: small(x, y, 'glacial'))
    result['F-troughs'] = inner('F', lambda x, y: small(x, y, 'glacial'))
    for name, item in result.items():
        print(name, item)
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=1) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
