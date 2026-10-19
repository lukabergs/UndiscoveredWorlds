"""Visible stretch of the terrain sandbox's oriented noise banks (W2-SB1 review): the bank alone (alignment 1) at a
constant strike on small patches, its stretch the rms of the derivative across the strike over the rms along it
(a noise stretched by S along the strike gives about S). The sphere's bank (terrain_sandbox.oriented_fbm at unit
vectors: patches on the tangent plane at random points below 70 degrees, random strikes) against the 2-D bank of
another copy of the sandbox (--old-scripts DIR: its scripts/sandbox, e.g. a frozen 2046df8; planar patches at random
offsets, random strikes). Banks: relief (structure: S 4), hardness bands (S 6), sea-floor hills (S 8), with their
configured orientations and sharpness, two octaves (the stretch is per octave), wavelength / 12 per sample.

  uv run python scripts/sandbox/bank_stretch.py [--old-scripts DIR] [--patches 40] [--size 96] [--json OUT]
"""
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import terrain_sandbox as ts  # noqa: E402


class Patch:
    """What oriented_fbm reads of a domain: points, frames, the radius and the sample spacing."""

    def __init__(self, points, fine_km, radius_km):
        self._points, self.fine_km, self.radius_km = points, fine_km, radius_km
        r = np.sqrt((points * points).sum(0))
        latitude, longitude = np.arcsin(points[2] / r), np.arctan2(points[1], points[0])
        so, co, sl, cl = np.sin(longitude), np.cos(longitude), np.sin(latitude), np.cos(latitude)
        self._frame = (np.stack([-so, co, np.zeros_like(so)]), np.stack([-sl * co, -sl * so, cl]))

    def points(self):
        return self._points

    def frame(self):
        return self._frame


def specs(config):
    structure = config['structure']
    options = ts.bank_options(structure)
    bands = config['hardness']['bands']
    return {'relief S 4': {**config['relief'], **structure},
            'bands S 6': {**bands, 'orientations': structure['orientations'], 'sharpness': structure['sharpness'],
                          **options},
            'hills S 8': {**config['seafloor']['hills'], **options}}


def stretch(field):
    """rms of the derivative across (axis 1) over along (axis 0)."""
    along, across = np.diff(field, axis=0), np.diff(field, axis=1)
    return float(np.sqrt((across ** 2).mean()) / np.sqrt((along ** 2).mean()))


def sphere_bank(spec, rng, patches, size, radius_km):
    out = []
    spacing = spec['wavelength_km'] / 12.0
    offsets = (np.arange(size) - size / 2.0) * spacing
    for _ in range(patches):
        latitude = math.asin(rng.uniform(-math.sin(math.radians(70)), math.sin(math.radians(70))))
        longitude = rng.uniform(-math.pi, math.pi)
        theta = rng.uniform(0.0, math.pi)
        n = np.array([math.cos(latitude) * math.cos(longitude), math.cos(latitude) * math.sin(longitude),
                      math.sin(latitude)])
        east = np.array([-math.sin(longitude), math.cos(longitude), 0.0])
        north = np.cross(n, east)
        strike = -math.sin(theta) * east - math.cos(theta) * north
        across = math.cos(theta) * east - math.sin(theta) * north
        u, v = np.meshgrid(offsets, offsets, indexing='ij')  # u along the strike (axis 0), v across (axis 1)
        tangent = n[:, None, None] * radius_km + u * strike[:, None, None] + v * across[:, None, None]
        points = radius_km * tangent / np.sqrt((tangent * tangent).sum(0))
        patch = Patch(points, spacing, radius_km)
        # theta per point in each point's own frame: the strike tangent projected there
        east_p, north_p = patch.frame()
        s = strike[:, None, None]
        theta_p = np.arctan2(-(s * east_p).sum(0), -(s * north_p).sum(0))
        field = ts.oriented_fbm(points, theta_p, np.ones(u.shape), spec, 1, patch)
        out.append(stretch(field))
    return out


def planar_bank(old, spec, rng, patches, size):
    out = []
    spacing = spec['wavelength_km'] / 12.0
    offsets = (np.arange(size) - size / 2.0) * spacing

    class Flat:
        fine_km = spacing
    for _ in range(patches):
        theta = rng.uniform(0.0, math.pi)
        origin = rng.uniform(0.0, 20000.0, size=2)
        strike = np.array([-math.sin(theta), math.cos(theta)])   # (x east, y south)
        across = np.array([math.cos(theta), math.sin(theta)])
        u, v = np.meshgrid(offsets, offsets, indexing='ij')
        x = origin[0] + u * strike[0] + v * across[0]
        y = origin[1] + u * strike[1] + v * across[1]
        field = old.oriented_fbm(x, y, np.full(u.shape, theta), np.ones(u.shape), spec, 1, Flat)
        out.append(stretch(field))
    return out


def percentiles(values):
    return [round(float(np.percentile(values, q)), 2) for q in (10, 50, 90)]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--old-scripts', default=None, help="another sandbox's scripts/sandbox (the 2-D bank)")
    parser.add_argument('--config', default=str(ts.DEFAULT_CONFIG))
    parser.add_argument('--set', action='append', default=[])
    parser.add_argument('--patches', type=int, default=40)
    parser.add_argument('--size', type=int, default=96)
    parser.add_argument('--json', default=None)
    arguments = parser.parse_args()
    config = ts.load_config(arguments.config, arguments.set)
    ts.use_noise_divisor(config)
    old = None
    if arguments.old_scripts:
        sys.path.insert(0, str(Path(arguments.old_scripts).parent))
        loader = importlib.util.spec_from_file_location('terrain_sandbox_old', Path(arguments.old_scripts) / 'terrain_sandbox.py')
        old = importlib.util.module_from_spec(loader)
        loader.loader.exec_module(old)
    record = {'patches': arguments.patches, 'size': arguments.size, 'percentiles': '10, 50, 90', 'banks': {}}
    for name, spec in specs(config).items():
        spec = {**spec, 'octaves': 2}
        entry = {'stretch': spec['stretch'], 'orientations': spec['orientations'], 'sharpness': spec['sharpness'],
                 'sphere': percentiles(sphere_bank(spec, np.random.default_rng(1), arguments.patches, arguments.size,
                                                   6371.0))}
        if old is not None:
            entry['planar_2d'] = percentiles(planar_bank(old, spec, np.random.default_rng(1), arguments.patches,
                                                         arguments.size))
            entry['sphere_over_2d'] = [round(a / b, 3) for a, b in zip(entry['sphere'], entry['planar_2d'])]
        record['banks'][name] = entry
        print(name, json.dumps(entry), flush=True)
    if arguments.json:
        Path(arguments.json).write_text(json.dumps(record, indent=1) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
