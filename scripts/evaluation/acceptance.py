"""The T8 acceptance page (terrain roadmap section 3, T8 "Acceptance"; step W3.6), rendered from JSON.

    uv run --offline python scripts/evaluation/acceptance.py [--scratch DIR] [--json OUT] [--output DIR]
        [--base-v10 FILE] [--scores STEP=FILE ...] [--no-maps]

The collector reads the evidence the T8 steps left under the main tree's out/scratch/sphere (DIR) and
writes one document (OUT, default DIR/acceptance/t8-acceptance.json) with a section per topic: the engine
(E6/acceptance.json, the gates as amended after E6), striping by band and A4 on the lat-lon view and
natively (E6, W2-FE), cube seam ranks (E6, W3-seam, W2-FE; cube_seams.py's rule: a seam only above the
null's 95th percentile), Earth scores lat-lon against cube (W3.3a, W3.3b, W3.3c, then each --scores
file of a later step in the same layout, recorded and never gated), timings and the D8
face-axis measure (W3.4), the cube windows' bit-identical overlap pairs (parsed from logs/W2-W3a..c.txt),
the sandbox global (W2-SB1) with a placeholder for base-v10 (W2-SB2), and the cube procedural world
(W2-FE) with two maps: the relief of its lat-lon export and its native fields on the cube's cross net,
both on the project palette (configs/tectonics-palettes.json, tectonics:elevation).

Every row is a gate: id, label, value, limit, result and the source file. Results: pass, fail,
accepted (over the original limit but accepted by an amendment of section 3), recorded (no gate: the
value is shown beside its reference) and missing (the input is absent). Rows with gated false keep the
result their limit would give but do not count in the summary (the lat-lon view of the engine's striping,
judged natively since E6; the sandbox global of base-v9, whose gates apply to base-v10). A missing input
gives a placeholder section or row, never an error. The page is scripts/evaluation/templates/acceptance.html
with the document embedded (report.write_page); images are linked by relative path.
"""
import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from uwtools import paths  # noqa: E402
import report  # noqa: E402

SCHEMA = 'uw-t8-acceptance/v1'
ROOT = Path(__file__).resolve().parents[2]
SEAM_RANK_LIMIT = 0.95                     # cube_seams.py: a seam only above the null's 95th percentile
TORUS_RMS_CM_YR = 0.937                    # T8.0: the torus engine at 400 Myr
EARTH_RMS_CM_YR = 4.39                     # configs/evaluation/plate-earth-reference.json
TIME_LIMIT_S = 1.5 * 477.0                 # section 3 "Time": 400 Myr at N 64 at most 1.5 x 477 s
TIME_ACCEPTED_S = 761.5                    # the amendment after E6: 761 s accepted (6 % over)
CELL_AREA_RATIO_LIMIT = 1.41               # section 3 "Cells"
POLAR_SHARE = 1.0 - math.sin(math.radians(60.0))   # the area share poleward of 60 degrees (13.4 %)
WINDOW_TIME_RATIO_LIMIT = 2.5              # W2-W3b decision: the cube-window speed budget relaxed to 2.5 x
STRIPING_BANDS = ('0-30', '30-50', '50-60', '60-70', '70-80')
RELIEF_EXAGGERATION = 100.0                # the maps' hillshade: about 300 km pixels at 128 columns


def default_scratch():
    try:
        return paths.main_root(ROOT) / 'out/scratch/sphere'
    except (OSError, subprocess.CalledProcessError):
        return ROOT / 'out/scratch/sphere'


def code_revision():
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


class Section:
    """One section of the page: its rows (gates), the files it read and the ones it missed."""

    def __init__(self, identifier, title, note=''):
        self.identifier, self.title, self.note = identifier, title, note
        self.sources, self.missing, self.gates, self.maps, self.data = [], [], [], [], {}

    def read(self, path):
        """The JSON document at path, or None (recorded as missing)."""
        path = Path(path)
        document = report.read_json(path)
        if document is None:
            self.missing.append(str(path))
        else:
            self.sources.append(str(path))
        return document

    def text(self, path):
        path = Path(path)
        if not path.is_file():
            self.missing.append(str(path))
            return None
        self.sources.append(str(path))
        return path.read_text(encoding='utf-8', errors='replace')

    def gate(self, identifier, label, value, limit, result, source, gated=True, group=None, note=None):
        if value is None and result not in ('recorded',):
            result = 'missing'
        self.gates.append({'id': identifier, 'label': label, 'value': value, 'limit': limit, 'result': result,
                           'gated': bool(gated) and result in ('pass', 'fail', 'accepted'), 'group': group,
                           'source': str(source), 'note': note})

    def document(self):
        if not self.gates and not self.maps:
            status = 'placeholder'
        elif self.missing:
            status = 'partial'
        else:
            status = 'complete'
        counts = count(self.gates)
        return {'id': self.identifier, 'title': self.title, 'note': self.note, 'status': status,
                'sources': self.sources, 'missing': self.missing, 'counts': counts, 'gates': self.gates,
                'maps': self.maps, 'data': self.data}


def count(gates):
    counts = {'pass': 0, 'fail': 0, 'accepted': 0, 'recorded': 0, 'missing': 0, 'not_gated': 0}
    for row in gates:
        if row['gated']:
            counts[row['result']] += 1
        elif row['result'] in ('pass', 'fail'):
            counts['not_gated'] += 1
        else:
            counts[row['result']] += 1
    return counts


def within(value, low=None, high=None):
    if value is None:
        return 'missing'
    return 'pass' if (low is None or value >= low) and (high is None or value <= high) else 'fail'


def dig(document, *keys):
    """document[k0][k1]..., or None when a key is absent."""
    for key in keys:
        if not isinstance(document, dict) or key not in document:
            return None
        document = document[key]
    return document


def guarded(build):
    """Run a section builder; an unreadable input becomes a placeholder note, not an error."""
    def run(section, *arguments):
        try:
            build(section, *arguments)
        except (KeyError, TypeError, ValueError, IndexError, AttributeError, OSError) as error:
            section.missing.append(f'unreadable input ({type(error).__name__}: {error})')
        return section.document()
    return run


# --- the engine -------------------------------------------------------------------------------------------

@guarded
def engine(section, scratch):
    path = scratch / 'E6/acceptance.json'
    document = section.read(path)
    if document is None:
        return
    run = document['n64_after']
    table = run['table']
    section.data = {'bundle': run.get('bundle'), 'scenario': document.get('scenario'),
                    'retuning': dig(document, 'retuning', 'key_changed')}
    for key, label in (('count', 'plate count'), ('largest_share', 'largest plate, share of the surface'),
                       ('plates_above_5_percent', 'plates above 5 % of the surface'),
                       ('spin_share_mean', 'spin share (mean over plates)')):
        entry = table['plates'][key]
        section.gate(f'plates.{key}', label, entry['value'], entry['limit'], 'recorded', path, group='plates')
    mean = dig(table, 'plates', 'rms_speed_mean_over_run_cm_yr', 'value')
    low = round(0.7 * TORUS_RMS_CM_YR, 4)
    section.gate('plates.rms_speed_mean_over_run', 'run-mean area-weighted RMS speed on the sphere (cm/yr)', mean,
                 f'{low} (0.7 x the torus engine) to {EARTH_RMS_CM_YR} (Earth), as amended after E6',
                 within(mean, low, EARTH_RMS_CM_YR), path, group='plates')
    end = table['plates']['rms_cm_yr_sphere']
    section.gate('plates.rms_end_of_run', 'end-of-run RMS speed (cm/yr)', end['value'],
                 f"{end['limit']}: replaced by the run-mean gate (varies about +-20 % over a run)", 'recorded',
                 path, group='plates')
    entry = table['max_displacement_over_edge']
    section.gate('displacement', 'largest |v| dt over the smallest edge in any update', entry['value'],
                 entry['limit'], entry['result'], path, group='plates')
    for key, entry in table['t412'].items():
        result = entry['result'] if entry['result'] in ('pass', 'fail') else 'recorded'
        section.gate(f't412.{key}', f'T4.12 {key}', entry['value'], entry['limit'], result, path,
                     group='T4.12 gates (true areas)')
    seconds = dig(table, 'time', 'value')
    result = 'missing' if seconds is None else 'pass' if seconds <= TIME_LIMIT_S else \
        'accepted' if seconds <= TIME_ACCEPTED_S else 'fail'
    section.gate('time', '400 Myr at N 64, wall seconds', seconds,
                 f'at most {TIME_LIMIT_S:.1f} s (1.5 x 477 s); 761 s accepted after E6 (transport substeps)',
                 result, path, group='time')
    cells(section, run.get('bundle'))


def cells(section, bundle_path):
    """Section 3 "Cells" from the engine bundle's cell table: largest over smallest area, polar share."""
    group = 'cells'
    if not bundle_path or not Path(bundle_path).is_dir():
        section.missing.append(str(bundle_path))
        section.gate('cells.area_ratio', 'largest over smallest cell area', None, f'at most {CELL_AREA_RATIO_LIMIT}',
                     'missing', bundle_path, group=group)
        return
    from uwtools.bundle import Bundle
    bundle = Bundle.open(bundle_path)
    table = bundle.cells()
    area = table['area_m2'].astype(float)
    latitude = table['latitude_deg'].astype(float)
    ratio = float(area.max() / area.min())
    share = float((np.abs(latitude) > 60.0).mean())
    source = Path(bundle_path) / 'grid'
    section.sources.append(str(source))
    section.gate('cells.area_ratio', 'largest over smallest cell area', round(ratio, 4),
                 f'at most {CELL_AREA_RATIO_LIMIT}', within(ratio, None, CELL_AREA_RATIO_LIMIT), source, group=group)
    section.gate('cells.polar_share', 'share of engine cells poleward of 60 degrees', round(share, 4),
                 f'within 0.02 of the band area share {POLAR_SHARE:.4f} (torus 0.33)',
                 within(share, POLAR_SHARE - 0.02, POLAR_SHARE + 0.02), source, group=group)


# --- striping ---------------------------------------------------------------------------------------------

@guarded
def striping(section, scratch):
    path = scratch / 'E6/acceptance.json'
    document = section.read(path)
    if document is not None:
        table = document['n64_after']['table']
        group = 'engine surface, lat-lon view (N 64, 400 Myr)'
        for key, entry in table['striping_engine_surface'].items():
            result = entry['result'] if entry['result'] in ('pass', 'fail') else 'recorded'
            section.gate(f'engine.latlon.{key}', key.replace('_', ' '), entry['value'], entry['limit'], result, path,
                         gated=False, group=group,
                         note='judged natively since E6: the lat-lon view at 70-80 degrees oversamples cube cells')
        group = 'engine surface, native (face coordinates)'
        for key, entry in table['cube_imprint'].items():
            if not key.startswith('a4_'):
                continue
            result = entry['result'] if entry['result'] in ('pass', 'fail') else 'recorded'
            section.gate(f'engine.native.{key}', key.replace('_', ' '), entry['value'], entry['limit'], result, path,
                         group=group)
    path = scratch / 'E6/striping-r2-n64-400.json'
    document = section.read(path)
    if document is not None:
        field = document['bundle']['fields']['geosphere.tectonic_height']
        for view, bands in (('lat-lon view', field['bands']), ('native', field['native']['bands'])):
            for band, entry in bands.items():
                section.gate(f'engine.a8.{view}.{band}', f'A8 {band}', entry.get('a8'), 'recorded (item 3)',
                             'recorded', path, group=f'engine surface, A8 ({view})')
    path = scratch / 'W2-FE/striping.json'
    document = section.read(path)
    if document is not None:
        for spec, field in document['bundle']['fields'].items():
            name = 'base elevation' if spec == 'geosphere.base_elevation' else 'surface (base + components)'
            for view, bands in (('lat-lon view', field['bands']), ('native', field['native']['bands'])):
                for band, entry in bands.items():
                    if 'ratio' in entry:
                        section.gate(f'fe.{name}.{view}.{band}.ratio', f'ratio {band}', entry['ratio'],
                                     'recorded (minimal world; the 0.8-1.4 gate applies to base-v10)', 'recorded',
                                     path, group=f'cube procedural world 128, {name}, {view}')
                    section.gate(f'fe.{name}.{view}.{band}.a4', f'A4 {band}', entry.get('a4'), 'recorded',
                                 'recorded', path, group=f'cube procedural world 128, {name}, {view}')


# --- seams ------------------------------------------------------------------------------------------------

@guarded
def seams(section, scratch):
    path = scratch / 'E6/acceptance.json'
    document = section.read(path)
    if document is not None:
        for key, entry in document['n64_after']['table']['cube_imprint'].items():
            if key.startswith('a4_'):
                continue
            result = entry['result'] if entry['result'] in ('pass', 'fail') else 'recorded'
            section.gate(f'engine.{key}', key.replace('_', ' '), entry['value'], entry['limit'], result, path,
                         group='engine surface (N 64, 400 Myr)')
    path = scratch / 'W3-seam/seams.json'
    document = section.read(path)
    if document is not None:
        for name, entries in document['nulls'].items():
            for bundle, statistics in entries.items():
                final = '-checkpoints' not in bundle and Path(bundle).name.startswith('before-')
                label = Path(bundle).name if final else '/'.join(Path(bundle).parts[-2:])
                group = f'cube procedural world, {label}' + ('' if final else ' (diagnostic)')
                for statistic, record in statistics.items():
                    rank = record.get('rank_median')
                    limit = (f"rank at most {SEAM_RANK_LIMIT} (ratio {record['ratio_median']:.3f}, "
                             f"null p5/p50/p95 {record['null_median_p5_p50_p95']})")
                    section.gate(f'w3seam.{label}.{statistic}.median', f'{statistic} ratio_median rank', rank, limit,
                                 within(rank, None, SEAM_RANK_LIMIT) if final else 'recorded', path, group=group)
                    section.gate(f'w3seam.{label}.{statistic}.mean', f'{statistic} ratio_mean rank',
                                 record.get('rank_mean'), 'recorded (the rule ranks ratio_median)', 'recorded', path,
                                 group=group)
    for name in ('seams.json', 'seams-land.json'):
        path = scratch / 'W2-FE' / name
        document = section.read(path)
        if document is None:
            continue
        for spec, record in document['fields'].items():
            field = 'base elevation' if spec == 'geosphere.base_elevation' else 'surface'
            part = 'land' if document.get('above') is not None else 'all cells'
            section.gate(f'fe.{name}.{field}', f'{field}, {part}: edge over interior median |grad|',
                         record.get('ratio_median'), 'recorded (before the rank rule; no null)', 'recorded', path,
                         group='cube procedural world 128 (W2-FE)')


# --- Earth scores -----------------------------------------------------------------------------------------

@guarded
def earth_scores(section, scratch, later=()):
    latest = 'W3.3c'
    steps = [(step, scratch / step / 'scores.json') for step in ('W3.3a', 'W3.3b', 'W3.3c')] + list(later)
    for step, path in steps:
        document = section.read(path)
        if document is None:
            continue
        for module, record in document['modules'].items():
            group = f'{module} ({step})'
            verdict = record.get('verdict', {})
            value = ' / '.join(f"{grid} {dig(verdict, grid, 'checks')}" for grid in ('latlon', 'cube'))
            failed = {grid: [entry['check'] for entry in record.get('failed_checks', {}).get(grid) or []]
                      for grid in ('latlon', 'cube')}
            section.gate(f'{step}.{module}.checks', 'evaluation checks (lat-lon / cube)', value,
                         f"failed: lat-lon {failed['latlon'] or 'none'}, cube {failed['cube'] or 'none'}", 'recorded',
                         path, group=group)
            scores = record.get('scores') or {}
            cube_score = dig(scores, 'cube', 'score')
            latlon_score = dig(scores, 'latlon', 'score')
            if cube_score is not None and latlon_score is not None:
                gated = step == latest
                section.gate(f'{step}.{module}.score', 'global-balanced-v1 on the cube (lower is better)',
                             round(cube_score, 3), f'at most the lat-lon run, {latlon_score:.3f}',
                             within(cube_score, None, latlon_score) if gated else 'recorded', path, group=group,
                             note=None if gated else f'superseded by {latest}')
            elif scores.get('common_fields'):
                common = scores['common_fields']
                section.gate(f'{step}.{module}.score_common', 'global-balanced-v1 on the shared fields (cube)',
                             round(common['cube'], 3), f"lat-lon {common['latlon']:.3f} on {len(common['fields'])} "
                             f"fields (weight {common['weight']})", 'recorded', path, group=group,
                             note=f'the cube score was partial; superseded by {latest}')
            for row in record.get('table') or []:
                if row.get('earth') is not None:
                    section.gate(f"{step}.{module}.{row['metric']}", f"{row['validation']}: {row['metric']}",
                                 row['cube'], f"lat-lon {report.fmt(row['latlon'])}, Earth {report.fmt(row['earth'])}",
                                 'recorded', path, group=group)
            for row in record.get('worse_on_cube') or []:
                metric = row['metric'] if 'validation' not in row else f"{row['validation']}: {row['metric']}"
                cube = row.get('cube', row.get('cube_gap'))
                latlon = row.get('latlon', row.get('latlon_gap'))
                section.gate(f'{step}.{module}.worse.{metric}', f'worse on the cube: {metric}', cube,
                             f"lat-lon {report.fmt(latlon)}; gap to Earth +{100 * row['relative_increase']:.0f} % "
                             f"(tolerance {100 * document.get('tolerance_relative', 0.05):.0f} %)",
                             'recorded', path, group=group)


# --- timings and D8 ---------------------------------------------------------------------------------------

def tectonics_seconds(run):
    for step in run.get('steps') or []:
        if step.get('process') == 'geology.plate_tectonics':
            return round(step['seconds'], 1)
    return None


@guarded
def timings(section, scratch):
    path = scratch / 'W3.4/timings.json'
    document = section.read(path)
    if document is None:
        return
    runs = document['runs']
    for name, run in runs.items():
        group = f"{name}: {run['grid']} {run['resolution']}" + (f" (N {run['face_cells']})" if run.get('face_cells') else '')
        section.gate(f'{name}.wall', 'wall seconds, full procedural recipe', run.get('wall_seconds'),
                     run.get('note') or 'recorded', 'recorded', path, group=group)
        section.gate(f'{name}.engine', 'geology.plate_tectonics seconds', tectonics_seconds(run), 'recorded',
                     'recorded', path, group=group)
        if run.get('peak_working_set_gib') is not None:
            section.gate(f'{name}.memory', 'peak working set (GiB)', run['peak_working_set_gib'], 'recorded',
                         'recorded', path, group=group)
    for cube, latlon in (('cube-512', 'latlon-512'), ('cube-256', 'latlon-256 (W3.3a)'),
                         ('cube-256 (W3.3a)', 'latlon-256 (W3.3a)')):
        if cube in runs and latlon in runs:
            ratio = runs[cube]['wall_seconds'] / runs[latlon]['wall_seconds']
            section.gate(f'ratio.{cube}', f'{cube} over {latlon}, wall time', round(ratio, 3),
                         'recorded (climate coupling dominates; machine shared)', 'recorded', path, group='cube over lat-lon')
    section.data['notes'] = document.get('notes')


@guarded
def d8(section, scratch):
    path = scratch / 'W3.4/d8.json'
    document = section.read(path)
    if document is None:
        return
    by_kind = {record['grid']['kind']: record for record in document['bundles'].values()}
    cube, latlon = by_kind.get('cube'), by_kind.get('latlon')
    if cube is None:
        section.missing.append(f'{path}: no cube bundle')
        return
    for group, entry in cube['groups'].items():
        if 'mean_cos4' not in entry:
            continue
        other = dig(latlon, 'groups', group) or {}
        limit = (f"lat-lon {report.fmt(other.get('mean_cos4'))} (random {report.fmt(other.get('random_mean'))}); "
                 f"cube random {report.fmt(entry['random_mean'])}, p95 {report.fmt(entry['random_p5_p50_p95'][2])}")
        section.gate(f'd8.{group}', f'mean cos 4 phi of river directions, {group}', round(entry['mean_cos4'], 4),
                     limit, 'recorded', path, group='face-axis D8 (grid axes)' if not group.startswith('face')
                     else 'by cube face')


# --- windows ----------------------------------------------------------------------------------------------

LEVEL1 = re.compile(r'^(\w+) shared cells (\d+) fields (\d+) identical (\d+)\s*$')
PAIR = re.compile(r'^(\w+) scale (\d+): overlap exit (\d+) \((\d+) shared cells, (\d+)/(\d+) fields identical\)'
                  r'(?:.*?ratio a ([\d.]+) x time ([\d.]+) x memory, b ([\d.]+) x time ([\d.]+) x memory)?')


def parse_windows(step, text):
    """The overlap pairs of a W2-W3 log: pair, scale, shared cells, fields, identical fields, exit code and
    the cube-over-lat-lon time and memory ratios of windows a and b where the log has them. A block under a
    '== Baseline' heading (an earlier build run for comparison) is skipped."""
    rows = []
    baseline = False
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('=='):
            baseline = line[2:].strip().lower().startswith('baseline')
            continue
        if baseline:
            continue
        match = LEVEL1.match(line)
        if match:
            rows.append({'step': step, 'pair': match[1], 'scale': 4, 'shared_cells': int(match[2]),
                         'fields': int(match[3]), 'identical': int(match[4]), 'exit': 0})
            continue
        match = PAIR.match(line)
        if match:
            row = {'step': step, 'pair': match[1], 'scale': int(match[2]), 'exit': int(match[3]),
                   'shared_cells': int(match[4]), 'identical': int(match[5]), 'fields': int(match[6])}
            if match[7]:
                row.update({'time_ratio': max(float(match[7]), float(match[9])),
                            'memory_ratio': max(float(match[8]), float(match[10]))})
            rows.append(row)
    return rows


@guarded
def windows(section, scratch):
    rows = []
    for step in ('W2-W3a', 'W2-W3b', 'W2-W3c'):
        path = scratch / 'logs' / f'{step}.txt'
        text = section.text(path)
        if text is None:
            continue
        found = parse_windows(step, text)
        if not found:
            section.missing.append(f'{path}: no overlap pairs found')
        for row in found:
            identical = row['identical'] == row['fields'] and row['exit'] == 0
            section.gate(f"{step}.{row['pair']}.{row['scale']}", f"{row['pair']} at scale {row['scale']}: "
                         f"{row['shared_cells']} shared cells", f"{row['identical']}/{row['fields']}",
                         'every field bit-identical on the shared cells', 'pass' if identical else 'fail', path,
                         group=f'{step}: overlap pairs')
            if 'time_ratio' in row:
                section.gate(f"{step}.{row['pair']}.{row['scale']}.time", f"{row['pair']} at scale {row['scale']}: "
                             'time over the lat-lon window (worse of a, b)', row['time_ratio'],
                             f'at most {WINDOW_TIME_RATIO_LIMIT} x (relaxed from 1.3 x, W2-W3b)',
                             within(row['time_ratio'], None, WINDOW_TIME_RATIO_LIMIT), path, group=f'{step}: speed')
                section.gate(f"{step}.{row['pair']}.{row['scale']}.memory", f"{row['pair']} at scale {row['scale']}: "
                             'memory over the lat-lon window', row['memory_ratio'], 'recorded', 'recorded', path,
                             group=f'{step}: speed')
        rows.extend(found)
    section.data['pairs'] = rows


# --- the sandbox global -----------------------------------------------------------------------------------

def sandbox_rows(section, run, path, prefix, gated, group):
    note = None if gated else 'base-v9 (torus engine): the section 3 gates apply to base-v10'
    section.gate(f'{prefix}.land_share', 'land share', run.get('land_share'), 'recorded', 'recorded', path,
                 group=group)
    section.gate(f'{prefix}.land_median', 'land median (m)', dig(run, 'land_quantiles_m', '0.5'), 'recorded',
                 'recorded', path, group=group)
    for layer, ratio_high, a4_limit in (('procedural', 1.4, 0.08), ('final', 1.5, None)):
        bands = dig(run, 'striping', layer) or {}
        for band in STRIPING_BANDS:
            entry = bands.get(band)
            if entry is None:
                continue
            low = 0.8 if layer == 'procedural' else None
            limit = f'{low}-{ratio_high}' if low else f'at most {ratio_high}'
            section.gate(f'{prefix}.{layer}.{band}.ratio', f'{layer} ratio {band}', entry.get('ratio'), limit,
                         within(entry.get('ratio'), low, ratio_high), path, gated=gated, group=group, note=note)
            if a4_limit is not None:
                section.gate(f'{prefix}.{layer}.{band}.a4', f'{layer} A4 {band}', entry.get('a4'),
                             f'at most {a4_limit}', within(entry.get('a4'), None, a4_limit), path, gated=gated,
                             group=group, note=note)
            else:
                section.gate(f'{prefix}.{layer}.{band}.a4', f'{layer} A4 {band}', entry.get('a4'), 'recorded',
                             'recorded', path, group=group)


@guarded
def sandbox(section, scratch):
    path = scratch / 'W2-SB1/totals-followup.json'
    document = section.read(path)
    if document is None:
        return
    for name, run in document['runs'].items():
        grid = run.get('grid')
        sandbox_rows(section, run, path, name, False, f"base-v9, {name} ({grid[1]} x {grid[0]})" if grid else name)


@guarded
def base_v10(section, path):
    document = section.read(path)
    if document is None:
        return
    runs = document['runs'] if isinstance(document.get('runs'), dict) else {'base-v10': document}
    for name, run in runs.items():
        sandbox_rows(section, run, path, name, True, f'base-v10, {name}')


# --- the cube procedural world and its maps ---------------------------------------------------------------

def surface_native(bundle):
    """Surface elevation relative to sea level on the native (6, N, N) cells, as validation.surface_elevation."""
    def load(identifier):
        return bundle.load_native(identifier).astype(float) if bundle.has(identifier) else 0.0
    base = load('geosphere.base_elevation')
    sea_level = float(bundle.scalar('geosphere.sea_level'))
    volcano = np.abs(load('geosphere.volcano'))
    land = base + np.maximum(volcano, load('geosphere.mountain_height')) + load('geosphere.extra_elevation') + \
        load('geosphere.crater_rims') + load('geosphere.crater_centres')
    sea = base + load('geosphere.ocean_ridge_height') + volcano
    lakes = load('hydrosphere.lake_surface')
    return np.where(base <= sea_level, sea, land) - sea_level, (base <= sea_level) & (np.asarray(lakes) == 0)


def shaded(elevation, sea, dy, dx):
    import maps
    shade = maps.shade_factor(maps.hillshade(elevation, dy, dx, RELIEF_EXAGGERATION))
    rgb = maps.terrain_base(elevation, sea) * shade[..., None]
    return np.clip(np.floor(rgb + 0.5), 0, 255).astype(np.uint8)


def export_relief(bundle, directory):
    import maps
    from uwtools import map_rendering
    from validation import surface_elevation
    elevation = surface_elevation(bundle)
    rows = elevation.shape[0]
    radius_m = bundle.radius_km * 1000.0
    dy, dx = maps.pixel_sizes_m(90.0, 180.0 / rows, rows, radius_m)
    rgb = shaded(elevation, bundle.sea_mask(), dy, dx)
    title = (f'Relief of the lat-lon export ({bundle.columns} x {rows}): surface elevation on tectonics:elevation '
             f'times a hillshade (light NW, 45 degrees up, exaggeration {RELIEF_EXAGGERATION:g})')
    image = map_rendering.frame(rgb, title, maps.legend('tectonics:elevation', 'm relative to sea level'))
    path = directory / 'relief-export.png'
    image.save(path)
    return path, title


def native_net(bundle, directory):
    """The native fields on the cross net 4 / 0 1 2 3 / 5 (uw-cube-1: faces 0-3 north up, face 4's bottom edge is
    face 0's top edge, face 5's top edge face 0's bottom edge), each face shaded on its own cells."""
    from PIL import Image, ImageDraw, ImageFont
    from uwtools import palettes
    elevation, sea = surface_native(bundle)
    n = elevation.shape[-1]
    cell_m = math.pi / 2 * bundle.radius_km * 1000.0 / n
    missing = palettes.tectonics()['missing']
    net = np.empty((3 * n, 4 * n, 3), dtype=np.uint8)
    net[:] = missing
    places = {4: (0, 0), 0: (1, 0), 1: (1, 1), 2: (1, 2), 3: (1, 3), 5: (2, 0)}
    for face, (row, column) in places.items():
        net[row * n:(row + 1) * n, column * n:(column + 1) * n] = shaded(elevation[face], sea[face], cell_m, cell_m)
    scale = max(1, -(-1024 // (4 * n)))
    raster = Image.fromarray(net).resize((4 * n * scale, 3 * n * scale), Image.Resampling.NEAREST)
    canvas = Image.new('RGB', (raster.width + 20, raster.height + 40), (250, 249, 245))
    canvas.paste(raster, (10, 30))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default(size=13)
    title = (f'Native cube cells (N {n}) on the cross net 4 / 0 1 2 3 / 5: surface elevation on tectonics:elevation '
             f'times a hillshade per face (exaggeration {RELIEF_EXAGGERATION:g})')
    draw.text((10, 8), title, font=font, fill=(31, 35, 40))
    for face, (row, column) in places.items():
        draw.text((10 + column * n * scale + 4, 30 + row * n * scale + 2), f'face {face}', font=font,
                  fill=(255, 255, 255))
    path = directory / 'native-net.png'
    canvas.save(path)
    return path, title


@guarded
def procedural_world(section, scratch, output, render):
    path = scratch / 'W2-FE/measures.json'
    document = section.read(path)
    if document is not None:
        cube, latlon = document['cube-procedural'], document['latlon-procedural']
        for key, label in (('land_share_area', 'land share by area'), ('river_cells', 'river cells'),
                           ('lake_count', 'lakes'), ('lake_area_km2', 'lake area (km2)')):
            section.gate(f'fe.{key}', label, cube.get(key), f'lat-lon {report.fmt(latlon.get(key))}', 'recorded',
                         path, group='cube N 32 against lat-lon 128 x 64 (minimal configuration)')
        for keys, label in ((('hypsometry_land_m', 'q50'), 'land median (m)'),
                            (('plates', 'tectonic_plate_count'), 'plates'),
                            (('temperature', 'mean'), 'mean temperature (C)'),
                            (('precipitation', 'mean'), 'mean precipitation (mm/month)'),
                            (('stage_seconds_total',), 'stage seconds')):
            section.gate(f"fe.{'.'.join(keys)}", label, dig(cube, *keys), f'lat-lon {report.fmt(dig(latlon, *keys))}',
                         'recorded', path, group='cube N 32 against lat-lon 128 x 64 (minimal configuration)')
    bundle_path = scratch / 'W2-FE/cube-procedural'
    if not render:
        return
    if not (bundle_path / 'manifest.json').is_file():
        section.missing.append(str(bundle_path))
        return
    from uwtools import palettes
    from uwtools.bundle import Bundle
    bundle = Bundle.open(bundle_path)
    directory = output / 'maps'
    directory.mkdir(parents=True, exist_ok=True)
    section.sources.append(str(bundle_path))
    for draw in (export_relief, native_net):
        image, title = draw(bundle, directory)
        section.maps.append({'path': report.link(image, output), 'title': title})
    palettes.write_sidecar(directory)


# --- the document and the page ----------------------------------------------------------------------------

SECTIONS = (
    ('engine', 'The engine (E6, baseline v6, N 64, 400 Myr)',
     'Section 3 T8 "Acceptance" as amended after E6: the speed gate is the run-mean RMS between 0.7 x the torus '
     "engine's and Earth's; 761 s accepted for 400 Myr at N 64. Conservation (tests) and identity (step_check) "
     'are not in these inputs.'),
    ('striping', 'Striping by band and A4: the lat-lon view and natively',
     "The engine's striping and A4 are judged natively since E6 (face coordinates; the seam rank below); the "
     'lat-lon view is shown with the result its limit would give, not counted.'),
    ('seams', 'Cube seam ranks',
     'cube_seams.py: the edge band statistics ranked against 200 randomly rotated virtual cubes; a seam only when '
     'the true ratio_median lies above the null\'s 95th percentile (W3-seam decision).'),
    ('earth', 'Earth scores, lat-lon against cube (256, seed 1)',
     'global-balanced-v1 is gated on the latest climate scores (cube at most lat-lon); the rest is recorded beside '
     'its lat-lon value and Earth. W3.3b\'s and W3.3c\'s hydrology scores predate the export lake rule (7769a80, '
     'f0c5209, roadmap W3.3c); a later step\'s rows (--scores) re-measure them after it.'),
    ('timings', 'Timings (W3.4, full procedural recipe, default duration)', ''),
    ('d8', 'D8 face-axis measure (W3.4)',
     'Mean cos 4 phi of river directions against the grid axes; random draws give the null.'),
    ('windows', 'Cube windows: bit-identical overlap pairs (W2-W3a to W2-W3c)',
     'Parsed from the step logs: interior, edge 0|1, edge 0|4 and corner (0, 1, 4) pairs at scales 4, 16 and 64 '
     '(levels 1 to 3) on the cube procedural world.'),
    ('sandbox', 'The sandbox global (W2-SB1, base-v9 at 2048 x 1024)',
     'Unit-vector noise on base-v9 (the torus engine). The section 3 striping gates are shown but not counted: '
     'they apply to base-v10.'),
    ('base-v10', 'The sandbox global of base-v10 (W2-SB2)',
     'Placeholder until W2-SB2 writes its totals (--base-v10 FILE, the W2-SB1 totals layout).'),
    ('world', 'The cube procedural world (W2-FE, N 32, minimal configuration)', ''),
)


def collect(scratch, output, base_v10_path=None, render=True, scores=()):
    """`scores`: (step, path) of later Earth-score files in the W3.3c layout, after W3.3c's."""
    scratch, output = Path(scratch), Path(output)
    base_v10_path = Path(base_v10_path) if base_v10_path else scratch / 'W2-SB2/totals.json'
    later = [(step, Path(path)) for step, path in scores]
    builders = {'engine': lambda s: engine(s, scratch), 'striping': lambda s: striping(s, scratch),
                'seams': lambda s: seams(s, scratch), 'earth': lambda s: earth_scores(s, scratch, later),
                'timings': lambda s: timings(s, scratch), 'd8': lambda s: d8(s, scratch),
                'windows': lambda s: windows(s, scratch), 'sandbox': lambda s: sandbox(s, scratch),
                'base-v10': lambda s: base_v10(s, base_v10_path),
                'world': lambda s: procedural_world(s, scratch, output, render)}
    sections = [builders[identifier](Section(identifier, title, note)) for identifier, title, note in SECTIONS]
    totals = count([row for section in sections for row in section['gates']])
    failing = [f"{section['id']}: {row['id']}" for section in sections for row in section['gates']
               if row['gated'] and row['result'] == 'fail']
    return {'schema': SCHEMA, 'generated': report.now(), 'code': code_revision(), 'scratch': str(scratch),
            'plan': 'docs/architecture/terrain-roadmap.md section 3, T8 "Acceptance"',
            'summary': {**totals, 'failing': failing,
                        'placeholders': [section['id'] for section in sections if section['status'] == 'placeholder']},
            'sections': sections}


def page_data(document, output):
    """The page's copy of the document: source paths as links relative to the page."""
    data = json.loads(json.dumps(report.clean(document)))
    for section in data['sections']:
        section['sources'] = [{'path': source, 'href': report.link(source, output)} for source in section['sources']]
        for row in section['gates']:
            row['href'] = report.link(row['source'], output) if row['source'] not in ('None', '') else None
    return data


def write(document, json_path, output):
    json_path, output = Path(json_path), Path(output)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(report.clean(document), indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    data = page_data(document, output)
    data['json'] = {'path': str(json_path), 'href': report.link(json_path, output)}
    return report.write_page('acceptance.html', data, output / 'index.html')


def main(argv=None):
    scratch = default_scratch()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--scratch', type=Path, default=scratch, help='the T8 scratch directory (out/scratch/sphere)')
    parser.add_argument('--json', type=Path, help='the document (default SCRATCH/acceptance/t8-acceptance.json)')
    parser.add_argument('--output', type=Path, default=scratch.parents[1] / 'evaluation/t8-acceptance',
                        help='the page directory (default out/evaluation/t8-acceptance of the main tree)')
    parser.add_argument('--base-v10', type=Path, help='base-v10 totals (default SCRATCH/W2-SB2/totals.json)')
    parser.add_argument('--scores', action='append', default=[], metavar='STEP=FILE',
                        help='a later step\'s Earth scores (the W3.3c scores.json layout), shown after W3.3c\'s and '
                             'recorded, never gated; repeatable')
    parser.add_argument('--no-maps', action='store_true', help='skip the two map renders')
    args = parser.parse_args(argv)
    if any('=' not in item for item in args.scores):
        parser.error('--scores takes STEP=FILE')
    json_path = args.json or args.scratch / 'acceptance/t8-acceptance.json'
    args.output.mkdir(parents=True, exist_ok=True)
    document = collect(args.scratch, args.output, args.base_v10, not args.no_maps,
                       [tuple(item.split('=', 1)) for item in args.scores])
    page = write(document, json_path, args.output)
    summary = document['summary']
    print(f"{json_path}\n{page}\ngated: {summary['pass']} pass, {summary['fail']} fail, {summary['accepted']} accepted; "
          f"{summary['not_gated']} shown but not gated, {summary['recorded']} recorded, {summary['missing']} missing; "
          f"placeholders: {', '.join(summary['placeholders']) or 'none'}")
    for item in summary['failing']:
        print(f'  fail {item}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
