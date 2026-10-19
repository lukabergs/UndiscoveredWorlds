"""Earth hydrology validation (erosion plan, Phase C, step C1): the routing of the earth-hydrology
recipe against the prepared HydroSHEDS fields of refs/processed/physical/fields/hydrology/ at the
run's width. The references validate; they never enter a simulation (world-model rule 4).

compare(kind, bundle, spec, baseline=None) returns {'reference': text, 'run': {'earth_comparison':
{...}, 'metrics': {...}}, 'baseline': {...}} with flat numbers, as landscape_statistics.py does.
The references are float32 W x W/2 GeoTIFFs, north row first, cell centres at the longitudes of
uwtools.grid.longitudes (checked against each file's transform), nodata -9999.9: values below -9000
are missing. A spec key `reference_directory` replaces the default directory (the unit tests use it).

Kinds:

- upstream_area: hydrosphere.drainage_area against rivers_upstream_area on the reference's river
  cells. For each reference cell the best model cell within one cell (3 x 3, longitude wraps, no
  neighbour beyond a pole) is the one with the smallest |log10(model / reference)| among those with a
  positive model value. Reported: the median |log10 ratio| over the matched reference cells, the
  share of all reference cells whose best match is within a factor of 2, the recall of the major
  rivers (reference cells with at least 1e5 km2 upstream: the share within a factor of 2) and the
  counts (a reference cell without a positive model value within one cell is unmatched, a miss).
- discharge: hydrosphere.discharge against rivers_discharge, the same statistics (major rivers by
  rivers_upstream_area); rivers_discharge is a WaterGAP model estimate, not observations. Reference
  cells with zero discharge are left out.
- basins: the routing's drainage basins against basin_level3. Every land cell is labelled by the end
  of its flow_direction path (terminal_cells): its outlet (the last land cell before the sea), the
  closed depression it enters, or the cell without a receiver where it stops. For each reference
  basin the model basin with the largest shared area; reported: the sum of those shared areas over
  the total reference basin area (area-weighted best-match overlap) and the basin counts. A cell on a
  flow_direction loop keeps the cell its pointer jumping stops on as its label; the Priority-Flood
  routing makes no loops, so this does not arise in a run.
- internal_drainage: endorheic land of the routing against internal_drainage (1 internally drained):
  intersection over union (area-weighted) and the endorheic shares of the land of each.

Endorheic land follows erosion plan B3, src/simulations/hydrology/drainage_basins.cpp
(closeddepressioncells and classifyerodedbasins, used by generate_drainage_basins on an eroded
world): a depression (depression_id) is closed when it has a lake cell (lake_surface != 0) and its
highest lake surface lies more than 0.01 m below its highest filled_elevation; a land cell is
endorheic when its path along flow_direction reaches a cell of a closed depression or a land cell
without a receiver (no direction, or one beyond a pole) before the sea; a path that never ends (a
cycle) counts as exorheic. earth-hydrology has no layers.finalize, so basin_class is not written:
the rule is applied here to the routing fields.
"""
import sys
from pathlib import Path

import numpy as np

from checks import DIRECTION_OFFSETS, whole_metres

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from uwtools import grid  # noqa: E402
from uwtools.paths import DATA_ROOT  # noqa: E402

KINDS = ('upstream_area', 'basins', 'internal_drainage', 'discharge')
REFERENCE_DIRECTORY = DATA_ROOT / 'refs' / 'processed' / 'physical' / 'fields' / 'hydrology'
MISSING_BELOW = -9000.0
MAJOR_RIVER_KM2 = 1e5
CLOSED_LAKE_MARGIN_M = np.float32(0.01)  # closedlakemarginm of drainage_basins.cpp
ROUTING = ('hydrosphere.flow_direction', 'geosphere.base_elevation', 'geosphere.sea_level')
DEPRESSIONS = ('hydrosphere.depression_id', 'hydrosphere.filled_elevation', 'hydrosphere.lake_surface')
FIELDS = {'upstream_area': ('hydrosphere.drainage_area',),
          'discharge': ('hydrosphere.discharge',),
          'basins': ROUTING + DEPRESSIONS,
          'internal_drainage': ROUTING + DEPRESSIONS}
REFERENCES = {'upstream_area': ('rivers_upstream_area',),
              'discharge': ('rivers_discharge', 'rivers_upstream_area'),
              'basins': ('basin_level3',),
              'internal_drainage': ('internal_drainage',)}
DESCRIPTION = {
    'upstream_area': 'HydroRIVERS v1.0 upstream area (rivers_upstream_area, reaches of at least 1000 km2; a cell holds '
                     'the largest touching reach), best model drainage area within one cell',
    'discharge': 'HydroRIVERS v1.0 discharge (rivers_discharge): a WaterGAP model estimate, not observations; '
                 'best model discharge within one cell',
    'basins': 'HydroBASINS level 3 (basin_level3) against the routing\'s basins, area-weighted best-match overlap',
    'internal_drainage': 'HydroBASINS internally drained catchments (internal_drainage) against the routing\'s '
                         'endorheic land (erosion plan B3 rule), area-weighted intersection over union'}


def skip_reason(kind, bundle, spec, baseline=None):
    """Why the kind cannot be compared (None when it can): a field the run (or the baseline) lacks,
    or no prepared reference at the run's width (they exist at 128, 256, 512, 1024 and 2048)."""
    missing = [field for field in FIELDS[kind] if field not in bundle.fields or
               (baseline is not None and field not in baseline.fields)]
    if missing:
        return 'model fields absent: ' + ', '.join(missing)
    directory = reference_directory(spec)
    absent = [name for name in REFERENCES[kind] if not (directory / name / f'{bundle.columns}.tif').is_file()]
    if absent:
        return f'no prepared reference at {bundle.columns} columns: ' + ', '.join(absent)
    return None


def reference_directory(spec):
    return Path(spec['reference_directory']) if 'reference_directory' in spec else REFERENCE_DIRECTORY


def load_reference(directory, name, columns):
    """A reference field at the run's width (rows x columns, NaN where missing); its grid must be the
    model's: columns x columns/2, cell centres at grid.longitudes and grid.latitudes."""
    import rasterio
    path = Path(directory) / name / f'{columns}.tif'
    with rasterio.open(path) as source:
        values = source.read(1).astype(float)
        transform = source.transform
    rows = columns // 2
    if values.shape != (rows, columns):
        raise ValueError(f'{path}: {values.shape[1]} x {values.shape[0]}, the run is {columns} x {rows}')
    a, b, c, d, e, f = tuple(transform)[:6]
    centre_lon = c + a * (np.arange(columns) + 0.5)
    centre_lat = f + e * (np.arange(rows) + 0.5)
    if not (b == 0 and d == 0 and np.allclose(centre_lon, grid.longitudes(columns), atol=1e-6) and
            np.allclose(centre_lat, grid.latitudes(rows), atol=1e-6)):
        raise ValueError(f'{path}: cell centres differ from the model grid (transform {tuple(transform)[:6]})')
    return np.where(values < MISSING_BELOW, np.nan, values)


def land_mask(bundle):
    """Land: base elevation in whole metres above sea level (planet::basebelowsea, the routing's sea)."""
    return whole_metres(bundle.load('geosphere.base_elevation').astype(float)) > float(bundle.scalar('geosphere.sea_level'))


def shifted(values, dx, dy, fill):
    """values at (x + dx, y + dy) of every cell; longitude wraps, `fill` beyond a pole."""
    result = np.roll(values, shift=(-dy, -dx), axis=(0, 1))
    if dy < 0:
        result[0, :] = fill
    elif dy > 0:
        result[-1, :] = fill
    return result


def best_log_ratio(model, reference):
    """Per cell: the smallest |log10(model / reference)| over the 3 x 3 neighbourhood with a positive
    model value (inf where there is none or the reference is missing)."""
    best = np.full(reference.shape, np.inf)
    valid = np.isfinite(reference) & (reference > 0)
    safe_reference = np.where(valid, reference, 1.0)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            candidate = shifted(model, dx, dy, 0.0)
            positive = candidate > 0
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.abs(np.log10(np.where(positive, candidate, 1.0) / safe_reference))
            best = np.where(valid & positive, np.minimum(best, ratio), best)
    return np.where(valid, best, np.inf)


def river_statistics(model, reference, upstream):
    """The statistics of upstream_area and discharge (flat numbers)."""
    cells = np.isfinite(reference) & (reference > 0)
    best = best_log_ratio(model, reference)[cells]
    matched = np.isfinite(best)
    within = best <= np.log10(2.0)
    major = (np.isfinite(upstream) & (upstream >= MAJOR_RIVER_KM2))[cells]
    count, major_count = int(cells.sum()), int(major.sum())
    comparison = {'median_abs_log10_ratio': float(np.median(best[matched])) if matched.any() else None,
                  'share_within_factor_2': float(within.sum() / count) if count else None,
                  'major_river_recall_within_factor_2': float(within[major].sum() / major_count) if major_count else None}
    counts = {'reference_cells': count, 'matched_cells': int(matched.sum()), 'unmatched_cells': int(count - matched.sum()),
              'major_river_cells': major_count}
    return comparison, counts


def terminal_cells(bundle):
    """(terminal index per cell, land, closed-depression cells, outlets, converged) along flow_direction.

    A land cell's path ends at an outlet (a land cell whose receiver is sea), at a land cell without a
    receiver (no direction, or one beyond a pole), or at a cell of a closed depression
    (closed_depression_cells); sea cells are their own terminal. Pointer jumping; a cell is converged
    when its terminal is a fixed point of the one-step receiver (a real end). A cell on or leading into
    a loop is not, whatever the loop's length (a loop of a power-of-two length ends the jumping on its
    own cells, which are no ends)."""
    direction = bundle.load('hydrosphere.flow_direction').astype(int)
    land = land_mask(bundle)
    rows, columns = direction.shape
    closed = closed_depression_cells(bundle)
    y, x = np.indices((rows, columns))
    receiver = y * columns + x
    outlet = np.zeros(land.shape, dtype=bool)
    for code, (dx, dy) in DIRECTION_OFFSETS.items():
        ny = y + dy
        cells = land & ~closed & (direction == code) & (ny >= 0) & (ny < rows)
        target = np.clip(ny, 0, rows - 1) * columns + (x + dx) % columns
        to_land = land.ravel()[target]
        receiver = np.where(cells & to_land, target, receiver)
        outlet |= cells & ~to_land
    receiver = receiver.ravel()
    step = receiver
    for _ in range(64):
        following = receiver[receiver]
        if np.array_equal(following, receiver):
            break
        receiver = following
    converged = step[receiver] == receiver
    return receiver.reshape(rows, columns), land, closed, outlet, converged.reshape(rows, columns)


def closed_depression_cells(bundle):
    """Cells of closed depressions (closeddepressioncells of drainage_basins.cpp, in float32)."""
    depression = bundle.load('hydrosphere.depression_id').astype(np.int64)
    filled = bundle.load('hydrosphere.filled_elevation').astype(np.float32)
    surface = bundle.load('hydrosphere.lake_surface').astype(np.float32)
    closed = np.zeros(depression.shape, dtype=bool)
    inside = depression != 0
    if not inside.any():
        return closed
    ids, inverse = np.unique(depression[inside], return_inverse=True)
    highest_filled = np.full(len(ids), -np.inf, dtype=np.float32)
    np.maximum.at(highest_filled, inverse, filled[inside])
    wet = surface[inside] != 0
    highest_lake = np.full(len(ids), -np.inf, dtype=np.float32)
    np.maximum.at(highest_lake, inverse[wet], surface[inside][wet])
    has_lake = np.zeros(len(ids), dtype=bool)
    has_lake[inverse[wet]] = True
    closed_ids = has_lake & (highest_lake < highest_filled - CLOSED_LAKE_MARGIN_M)
    closed[inside] = closed_ids[inverse]
    return closed


def endorheic_land(bundle):
    """Land cells whose path ends in a closed depression or on land without a receiver (erosion plan B3)."""
    terminal, land, closed, outlet, converged = terminal_cells(bundle)
    end = terminal.ravel()
    endorheic = land.ravel() & converged.ravel() & (closed.ravel()[end] | ~outlet.ravel()[end])
    return endorheic.reshape(land.shape)


def basin_labels(bundle):
    """Per land cell its basin: the terminal cell index, or (for a closed depression) size + its id;
    -1 on sea."""
    terminal, land, closed, _, _ = terminal_cells(bundle)
    depression = bundle.load('hydrosphere.depression_id').astype(np.int64)
    end = terminal.ravel()
    labels = np.where(closed.ravel()[end], terminal.size + depression.ravel()[end], end)
    return np.where(land, labels.reshape(land.shape), -1)


def measure_basins(bundle, reference, weights):
    labels = basin_labels(bundle)
    basins = np.isfinite(reference)
    total = float(weights[basins].sum())
    on_land = basins & (labels >= 0)
    reference_ids, reference_index = np.unique(reference[on_land].astype(np.int64), return_inverse=True)
    model_ids, model_index = np.unique(labels[on_land], return_inverse=True)
    shared = np.zeros((len(reference_ids), len(model_ids)))
    np.add.at(shared, (reference_index, model_index), weights[on_land])
    best = shared.max(axis=1).sum() if shared.size else 0.0
    return ({'best_match_overlap': float(best / total) if total else None},
            {'reference_basins': int(len(np.unique(reference[basins]))),
             'model_basins_on_reference_land': int(len(model_ids)),
             'model_basins': int(len(np.unique(labels[labels >= 0])))})


def measure_internal_drainage(bundle, reference, weights):
    model = endorheic_land(bundle)
    land = land_mask(bundle)
    valid = np.isfinite(reference)
    observed = valid & (reference >= 0.5)
    union = float(weights[model | observed].sum())
    return ({'intersection_over_union': float(weights[model & observed].sum() / union) if union else None,
             'model_endorheic_land_share': float(weights[model].sum() / weights[land].sum()) if land.any() else None,
             'reference_endorheic_land_share': float(weights[observed].sum() / weights[valid].sum()) if valid.any() else None},
            {'model_endorheic_cells': int(model.sum()), 'reference_endorheic_cells': int(observed.sum())})


def measure(kind, bundle, directory):
    """(comparison with the reference, counts) of one bundle: flat numbers."""
    columns = bundle.columns
    if kind in ('upstream_area', 'discharge'):
        upstream = load_reference(directory, 'rivers_upstream_area', columns)
        if kind == 'upstream_area':
            return river_statistics(bundle.load('hydrosphere.drainage_area').astype(float), upstream, upstream)
        reference = load_reference(directory, 'rivers_discharge', columns)
        return river_statistics(bundle.load('hydrosphere.discharge').astype(float), reference, upstream)
    weights = bundle.area_weights()
    if kind == 'basins':
        return measure_basins(bundle, load_reference(directory, 'basin_level3', columns), weights)
    if kind == 'internal_drainage':
        return measure_internal_drainage(bundle, load_reference(directory, 'internal_drainage', columns), weights)
    raise ValueError(f'Unknown hydrology validation kind {kind!r}; known: {", ".join(KINDS)}')


def major_river_cells(directory, columns):
    """[lon, lat] of the reference cells with at least 1e5 km2 upstream (the overlay's markers)."""
    upstream = load_reference(directory, 'rivers_upstream_area', columns)
    rows, cols = np.nonzero(np.isfinite(upstream) & (upstream >= MAJOR_RIVER_KM2))
    longitudes, latitudes = grid.longitudes(columns), grid.latitudes(columns // 2)
    return [[float(longitudes[x]), float(latitudes[y])] for y, x in zip(rows, cols)]


def compare(kind, bundle, spec, baseline=None):
    directory = reference_directory(spec)
    result = {'reference': f'{DESCRIPTION[kind]} at {bundle.columns} x {bundle.rows}'}
    for key, source in (('run', bundle), ('baseline', baseline)):
        if source is not None:
            if source.columns != bundle.columns:
                raise ValueError(f"{spec['id']}: run and baseline grids differ")
            comparison, counts = measure(kind, source, directory)
            result[key] = {'earth_comparison': comparison, 'metrics': counts}
    if kind == 'upstream_area':
        result['major_river_cells'] = major_river_cells(directory, bundle.columns)
    return result
