"""Statistics of one tectonic snapshot, computed from a run bundle.

    uv run --offline python scripts/evaluation/tectonics_statistics.py BUNDLE [--json OUT]

The cell counts come from the bundle's fields (geologic regime, crust class, boundary
type, nearest boundary, deforming region type, convergence, divergence and shear
scores); the segment and region statistics from the tables
geosphere.tectonic_boundary_segments and geosphere.tectonic_deforming_regions; time,
grid, plates and junctions from the snapshot record (geosphere.tectonic_snapshot).
Junctions are not in the bundle: their count and their stability counts are engine
diagnostics of the record. Field hashes are in the bundle manifest.

The counts are taken from the fields as the bundle holds them. With the tectonics recipe
(and its checkpoints) these are the engine's last update on the world grid; in
procedural-world later terrain steps (tectonic trenches and volcanoes) rewrite the
geologic regime, so its regime counts describe the final world, not the engine. Codes are those of
src/simulations/geology/plate_tectonics/tectonic_contract.hpp.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import runs
from uwtools.bundle import Bundle

SNAPSHOT_TABLE = 'geosphere.tectonic_snapshot'
SEGMENT_TABLE = 'geosphere.tectonic_boundary_segments'
REGION_TABLE = 'geosphere.tectonic_deforming_regions'
REGIMES = ('stable', 'convergent_arc', 'continent_collision', 'divergent_rift', 'transform', 'passive_margin',
           'mid_ocean_ridge', 'trench_adjacent')
CRUST_CLASSES = ('none', 'oceanic', 'transitional', 'continental')
BOUNDARY_TYPES = ('none', 'convergent', 'divergent', 'transform', 'passive_margin')
REGION_TYPES = ('none', 'continental_rift', 'diffuse_collision')
NO_BOUNDARY = 0  # tectonics::contract::kNoBoundaryId


def snapshot_record(bundle):
    return json.loads(bundle.table_text(SNAPSHOT_TABLE))


def mean(values):
    """Mean as the engine's code took it: a running double sum in list order."""
    total = 0.0
    for value in values:
        total += value
    return total / len(values) if values else 0.0


def code_counts(values, names):
    counts = np.bincount(values.ravel().astype(np.int64), minlength=len(names))
    return {name: int(counts[code]) for code, name in enumerate(names)}


def boundary_statistics(bundle):
    types = bundle.load('geosphere.boundary_type')
    cells = code_counts(types, BOUNDARY_TYPES)
    segments = bundle.table_rows(SEGMENT_TABLE)
    segment_types = [BOUNDARY_TYPES[int(row['boundary_type'])] for row in segments]
    ages = [float(row['age_myr']) for row in segments]
    return {
        'segment_count': len(segments),
        **{f'{name}_segment_count': segment_types.count(name) for name in BOUNDARY_TYPES[1:]},
        'active_cell_count': int(np.count_nonzero(types)),
        **{f'{name}_cell_count': cells[name] for name in BOUNDARY_TYPES[1:]},
        'nearest_mapped_cell_count': int(np.count_nonzero(bundle.load('geosphere.nearest_boundary') != NO_BOUNDARY)),
        'max_segment_age_myr': max(ages, default=0.0),
        'mean_segment_age_myr': mean(ages),
        'mean_segment_persistence_steps': mean([float(row['persistence_steps']) for row in segments]),
    }


def deforming_statistics(bundle):
    types = bundle.load('geosphere.deforming_region_type')
    cells = code_counts(types, REGION_TYPES)
    regions = bundle.table_rows(REGION_TABLE)
    region_types = [REGION_TYPES[int(row['type'])] for row in regions]
    ages = [float(row['age_myr']) for row in regions]
    return {
        'region_count': len(regions),
        **{f'{name}_region_count': region_types.count(name) for name in REGION_TYPES[1:]},
        'active_cell_count': int(np.count_nonzero(types)),
        **{f'{name}_cell_count': cells[name] for name in REGION_TYPES[1:]},
        'max_region_age_myr': max(ages, default=0.0),
        'mean_region_age_myr': mean(ages),
        'mean_region_persistence_steps': mean([float(row['persistence_steps']) for row in regions]),
        # The table prints the engine's float32 rates with nine digits: float32 restores them exactly.
        'mean_deformation_rate': mean([float(np.float32(row['deformation_rate'])) for row in regions]),
    }


def statistics(bundle):
    """Counts per regime, crust class and boundary type, active cells, segment, region and junction statistics."""
    bundle = bundle if isinstance(bundle, Bundle) else Bundle.open(bundle)
    record = snapshot_record(bundle)
    return {
        'time_myr': record['time_myr'],
        'iteration_count': record['iteration_count'],
        'width': record['width'],
        'height': record['height'],
        'plate_count': len(record['plates']),
        'tectonic_junction_count': record['tectonic_junction_count'],
        'tectonic_junction_stability_counts': record.get('tectonic_junction_stability_counts'),
        'active_convergence_cells': int(np.count_nonzero(bundle.load('geosphere.tectonic_convergence'))),
        'active_divergence_cells': int(np.count_nonzero(bundle.load('geosphere.tectonic_divergence'))),
        'active_shear_cells': int(np.count_nonzero(bundle.load('geosphere.tectonic_shear'))),
        'regime_counts': code_counts(bundle.load('geosphere.geologic_regime'), REGIMES),
        'crust_class_counts': code_counts(bundle.load('geosphere.crust_class'), CRUST_CLASSES),
        'boundary': boundary_statistics(bundle),
        'deforming': deforming_statistics(bundle),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bundle', type=Path)
    parser.add_argument('--json', type=Path, help='Write the statistics to this file')
    args = parser.parse_args()
    text = json.dumps(statistics(args.bundle), indent=2) + '\n'
    if args.json:
        runs.refuse_archive(args.json).write_text(text, encoding='utf-8')
    print(text, end='')


if __name__ == '__main__':
    main()
