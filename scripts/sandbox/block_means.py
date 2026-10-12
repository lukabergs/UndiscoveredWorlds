"""How much a sandbox run's after-physics filter moves the averages over world cells: for each window, the
mean over every k x k block (one 256-column world cell) of the elevation after the filter minus before it
(<domain>-elevation.npy and <domain>-elevation-physics.npy, both written by the run), with the blocks'
standard deviation of the filter's change for scale.

  uv run python scripts/sandbox/block_means.py --labels fine-filter-3,fine-filter-3m --domains B-core,G-coast --k 256
"""
import argparse
from pathlib import Path

import numpy as np

SANDBOX = Path(__file__).resolve().parents[2] / 'out' / 'scratch' / 'sandbox'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--labels', required=True)
    parser.add_argument('--domains', required=True)
    parser.add_argument('--k', type=int, required=True)
    arguments = parser.parse_args()
    k = arguments.k
    print(f'{"run":<18} {"window":<12} {"blocks":>6} {"mean shift m":>13} {"largest |shift| m":>18} '
          f'{"filter change sd m":>19}')
    for label in arguments.labels.split(','):
        for domain in arguments.domains.split(','):
            after = SANDBOX / label / f'{domain}-elevation.npy'
            before = SANDBOX / label / f'{domain}-elevation-physics.npy'
            if not (after.is_file() and before.is_file()):
                print(f'{label:<18} {domain:<12} missing')
                continue
            change = np.load(after).astype(float) - np.load(before).astype(float)
            rows, columns = change.shape[0] // k, change.shape[1] // k
            blocks = change[:rows * k, :columns * k].reshape(rows, k, columns, k).mean(axis=(1, 3))
            print(f'{label:<18} {domain:<12} {blocks.size:>6} {blocks.mean():>13.1f} {np.abs(blocks).max():>18.1f} '
                  f'{change.std():>19.1f}')


if __name__ == '__main__':
    main()
