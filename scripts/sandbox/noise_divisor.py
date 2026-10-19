"""The terrain sandbox's 3-D noise divisor (configs/sandbox/terrain.json noise.divisor_3d): the standard deviation of
one octave of the 2-D gradient noise (terrain_sandbox.gradient_noise) and of the 3-D one (gradient_noise3) on random
lattice positions, and the divisor 0.22 x std3 / std2 that keeps every fbm's spread what it was with the 2-D noise
and its divisor NOISE_STD 0.22.

  uv run python scripts/sandbox/noise_divisor.py [--seeds 1,7,131] [--samples 2000000] [--rng 1]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import terrain_sandbox as ts  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--seeds', default='1,7,131')
    parser.add_argument('--samples', type=int, default=2_000_000)
    parser.add_argument('--rng', type=int, default=1, help='seed of the random positions')
    arguments = parser.parse_args()
    positions = np.random.default_rng(arguments.rng).uniform(-500.0, 500.0, size=(3, arguments.samples))
    seeds = [int(seed) for seed in arguments.seeds.split(',')]
    two = [float(np.std(ts.gradient_noise(positions[0], positions[1], seed))) for seed in seeds]
    three = [float(np.std(ts.gradient_noise3(positions, seed))) for seed in seeds]
    record = {'samples': arguments.samples, 'rng': arguments.rng, 'seeds': seeds, 'std_2d': two, 'std_3d': three,
              'divisor_2d': ts.NOISE_STD, 'divisor_3d': round(ts.NOISE_STD * float(np.mean(three)) / float(np.mean(two)), 4)}
    print(json.dumps(record))


if __name__ == '__main__':
    main()
