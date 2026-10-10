"""Report pages rendered from JSON: run page, comparison page and module dashboard.

    uv run --offline python scripts/evaluation/report.py run <evaluation-dir>
    uv run --offline python scripts/evaluation/report.py compare <evaluation-dir> <evaluation-dir>... --output DIR [--earth NAME]
    uv run --offline python scripts/evaluation/report.py dashboard <module> --output DIR

Each page is a template of scripts/evaluation/templates/ (run.html, comparison.html,
dashboard.html) with the page's data embedded as JSON; the template's script builds the
page. Data come from a run's evaluation/evaluation.json, evaluation/review.json and
run.json, the module registry (runs/registry/<module>.json) and the ledger
(runs/ledger/<module>.jsonl). Images are linked by relative path. The run page writes
report.md from the same data, the comparison page a Markdown metric table.
Map entries of kind "earth" (Earth reference windows beside regional maps, earth_windows.py) get
their own column on both pages; `compare --earth NAME` draws that catalogue entry for the first
run's map into the comparison's directory and shows it instead of the runs' own windows.
"""
import argparse
import datetime
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import ledger, runs
from uwtools.bundle import SEASONS
import review

ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = Path(__file__).with_name('templates')
PLACEHOLDER = '__PAGE_DATA__'
PROFILES = ROOT / 'configs/evaluation'
METRIC_NAMES = {'rmse': ('RMSE', 'Vector RMSE'), 'bias': ('Bias', 'Speed bias'),
                'correlation': ('Spatial correlation', 'Speed spatial correlation'),
                'variance_ratio': ('Variance ratio', 'Speed variance ratio')}
KINDS = ('model', 'reference', 'difference', 'earth')


def fmt(value, digits=4):
    if value is None:
        return '–'
    if isinstance(value, float):
        return f'{value:.{digits}g}'
    return str(value)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def read_json(path):
    path = Path(path)
    return json.loads(path.read_text(encoding='utf-8-sig')) if path.is_file() else None


def clean(value):
    """JSON-safe copy: non-finite floats become null."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    return value


def link(path, start):
    """Relative link from directory start to path (a file URI across drives)."""
    try:
        return Path(os.path.relpath(Path(path).absolute(), Path(start).absolute())).as_posix()
    except ValueError:
        return Path(path).absolute().as_uri()


def write_page(template, data, path):
    text = (TEMPLATES / template).read_text(encoding='utf-8')
    if text.count(PLACEHOLDER) != 1:
        raise ValueError(f'{template}: expected one {PLACEHOLDER}')
    payload = json.dumps(clean(data), ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    path = runs.refuse_archive(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.replace(PLACEHOLDER, payload), encoding='utf-8')
    return path


def score_value(entry):
    return entry.get('score', entry.get('partial_score')) if entry else None


def season_label(seasons, season):
    return 'annual' if len(seasons) == 1 else SEASONS[int(season)]


def validation_rows(validation):
    """Per field and season: run loss, baseline loss, delta and the key metrics of the run."""
    rows = []
    for field, seasons in validation['maps'].items():
        for season, entries in seasons.items():
            metrics = {name: value for name, value, _ in entries['run']['rows']}
            loss, base = entries['run']['loss'], entries.get('baseline', {}).get('loss')
            row = {'field': field, 'season': season_label(seasons, season), 'loss': loss, 'baseline': base,
                   'delta': None if base is None else loss - base}
            for key, names in METRIC_NAMES.items():
                row[key] = next((metrics[name] for name in names if name in metrics), None)
            rows.append(row)
    return rows


def field_loss(validation, field, key='run'):
    seasons = validation['maps'].get(field)
    if not seasons or any(key not in entry for entry in seasons.values()):
        return None
    return sum(entry[key]['loss'] for entry in seasons.values()) / len(seasons)


def score_lines(validation):
    lines = []
    for key, entry in validation['scores'].items():
        lines.append({'key': key, 'version': entry['version'], 'inputs': entry['inputs'], 'score': entry.get('score'),
                      'partial_score': entry.get('partial_score'), 'covered_weight': entry['covered_weight'],
                      'missing_fields': entry.get('missing_fields', []), 'note': entry['note']})
    return lines


def distribution_tables(validation):
    tables = []
    for name, entry in validation['distributions'].items():
        keys = [key for key in ('run', 'baseline') if key in entry]
        rows = []
        for group in ('earth_comparison', 'metrics'):
            names = [metric for metric, value in entry[keys[0]][group].items() if not isinstance(value, dict)] if keys else []
            rows += [{'metric': metric, **{key: entry[key][group].get(metric) for key in keys}} for metric in names]
        tables.append({'name': name, 'reference': entry.get('reference'), 'columns': keys, 'rows': rows})
    return tables


def regional_tables(validation):
    """Regional errors (validation.regions) per field and season: rows [region, cells, run bias, run RMSE,
    baseline bias, baseline RMSE] for each region set (values to four significant digits, as displayed),
    or None without a regional table."""
    entry = validation.get('regions')
    if not entry:
        return None
    short = lambda value: None if value is None else float(f'{value:.4g}')
    options = []
    for field, seasons in entry['fields'].items():
        for season, sets in seasons.items():
            tables = {name: [[region, row['cells'], *[short(row.get(key, {}).get(metric)) for key in ('run', 'baseline')
                                                      for metric in ('bias', 'rmse')]] for region, row in rows.items()]
                      for name, rows in sets.items()}
            options.append({'field': field, 'season': season_label(seasons, season), 'sets': tables})
    return {'note': entry['note'], 'method': entry['method'], 'set_names': {'project': 'Project boxes', 'ar6': 'IPCC AR6'},
            'names': entry['sets'], 'options': options}


def map_groups(maps):
    """Static images grouped by map and season ({model, reference, difference}), and the animations."""
    groups, index, animations = [], {}, []
    for item in maps:
        kind = item.get('kind', 'model')
        record = {'title': item['title'], 'path': item['path'], 'half_range': item.get('half_range')}
        if 'display_fraction' in item:
            record['display_fraction'] = item['display_fraction']
        if kind == 'animation':
            animations.append(dict(record, map=item.get('map'), season=item.get('season'), frames=item.get('frames', [])))
            continue
        key = (item.get('map', item['title']), item.get('season'))
        if key not in index:
            index[key] = {'map': key[0], 'season': key[1], 'title': item['title'], 'validation': item.get('validation'),
                          'images': {}}
            groups.append(index[key])
        if kind == 'model':
            index[key]['title'] = item['title']
        index[key]['images'][kind] = record
    return groups, animations


def run_identity(directory, run):
    if run:
        return str(run.get('id') or run.get('name'))
    return directory.parent.name if directory.name == 'evaluation' else directory.name


def run_facts(document, run):
    bundle = document['run']
    facts = [['Module', document['profile']['module']]]
    if run and run.get('label'):
        facts.append(['Label', run['label']])
    facts += [['Seed', bundle['seed']], ['Recipe', bundle['recipe']], ['Grid', f"{bundle['grid'][0]} × {bundle['grid'][1]}"],
              ['Bundle', bundle['bundle']]]
    if document.get('baseline'):
        facts.append(['Baseline bundle', document['baseline']['bundle']])
    if document.get('coarser'):
        facts.append(['Coarser bundle', document['coarser']['bundle']])
    code = (run or {}).get('code') or document['code']
    revision = (code.get('revision') or 'unknown')[:12]
    facts.append(['Code', revision + (f" + diff {code['diff_sha256'][:12]}" if code.get('dirty') else '')])
    steps = (run or {}).get('runs') or []
    command = steps[0]['command'] if steps else bundle.get('command', [])
    facts.append(['Command', ' '.join(map(str, command))])
    if run:
        facts += [[label, run[key]] for label, key in (('Started', 'started'), ('Finished', 'finished')) if run.get(key)]
        facts += [[f"{step.get('label') or Path(step.get('log', 'run')).stem} time", f"{step['seconds']:.1f} s, exit {step['exit_code']}"]
                  for step in steps if 'seconds' in step]
    facts.append(['Evaluated', document['evaluated']])
    return facts


# The per-level diagnostics of a refinement hierarchy bundle (multiscale plan S2.3): column title and
# manifest key, in the order of the run page's table.
LEVEL_COLUMNS = (('Level', 'level'), ('Cell size (km)', 'cellsize_km'), ('Parent-mean error (m)', 'maxparentmeanerror_m'),
                 ('Bound (m)', 'meanbound_m'), ('Global-mean error (m)', 'globalmeanerror_m'),
                 ('Block-mean error before projection (m)', 'blockmeanerror_m'),
                 ('Border step max (m)', 'maxborderstep_m'), ('Border step p99 (m)', 'borderstep99_m'),
                 ('Smoothing passes', 'smoothpasses'), ('λ max before smoothing (m)', 'maxlambdabefore_m'),
                 ('λ step max before smoothing (m)', 'maxlambdastepbefore_m'), ('Smoothing move max (m)', 'maxsmoothmove_m'),
                 ('λ max (m)', 'maxlambda_m'), ('λ step max (m)', 'maxlambdastep_m'), ('Beds moved', 'bedsmovedparents'),
                 ('Forced', 'forcedparents'), ('Forced move max (m)', 'maxforcedmove_m'), ('Dried', 'driedcells'), ('Shores below lakes', 'shoresbelowlake'), ('Land below sea', 'landbelowsea'), ('Land below sea max depth (m)', 'maxlandbelowsea_m'),
                 ('Deepest land below sea (column, row)', 'maxlandbelowsea_position'),
                 ('Rim raised', 'rimraised'), ('Rim raise max (m)', 'maxrimraise_m'),
                 ('Border step max at parent (column, row)', 'maxborderstep_parent'),
                 ('Forced move max at parent (column, row)', 'maxforced_parent'),
                 ('Sea in land parents', 'seainlandparents'), ('Sea beyond one parent cell', 'seabeyondparent'),
                 ('Dry below a lake', 'drybelowlake'),
                 ('Prolongation clamped', 'prolongationclamped'), ('Prolongation clamp max (m)', 'maxprolongationclamp_m'),
                 ('Smoothing move p99 (m)', 'smoothmove99_m'), ('Channel rise above bed max (m)', 'maxchannelrise_m'),
                 ('Land cells', 'landcells'), ('Channel cells', 'channelcells'),
                 ('Drainage density (1/km)', 'drainagedensity_per_km'), ('Tributaries', 'tributaries'),
                 ('Tributary cells', 'tributarycells'), ('Lake edges', 'lakeedges'),
                 ('Lake edges on global cell boundaries (share)', 'lakeedgealigned'),
                 ('Land roughness, RMS Laplacian (m)', 'landroughness_m'), ('Valley cells', 'valleycells'),
                 ('Trough cells', 'troughcells'), ('Cirque cells', 'cirquecells'),
                 ('Seconds', 'seconds'))


def level_value(entry, key):
    """A level diagnostic; the positions (S2.8) as 'column, row', None when absent or -1."""
    if key.endswith('_position') or key.endswith('_parent'):
        stem = key[:-len('_position')] if key.endswith('_position') else key
        column, row = entry.get(stem + '_column'), entry.get(stem + '_row')
        return None if column is None or column < 0 else f'{column}, {row}'
    return entry.get(key)


def level_rows(levels):
    """The run page's per-level table {'columns', 'rows'} from evaluation.json's levels, or None (a single pass)."""
    if not levels:
        return None
    return {'columns': [title for title, _ in LEVEL_COLUMNS],
            'rows': [[level_value(entry, key) for _, key in LEVEL_COLUMNS] for entry in levels],
            'note': 'Per level of the refinement hierarchy, from the bundle\'s diagnostics.levels: the largest parent-mean '
                    'error against its bound (half a float step of the kept box\'s largest |elevation|), the global-mean '
                    'error, the block-mean error the presmoothing left before the projection, the border steps, the '
                    'smoothing before the projection (S2.7: its passes, λ and the largest λ step before it, its largest '
                    'move), λ and the largest λ step between neighbouring parents at the projection (from S2.7 on after the '
                    'smoothing; bundles without smoothing, regional runs 16 to 25, have no "before" columns and there λ '
                    'before = after), the parents whose beds moved or that were forced (with the '
                    'largest forced move), and the dried cells, shores below lakes and land below sea after it; S2.8\'s '
                    'mask checks: sea cells in land parents (level 1: global cells), those beyond one parent cell of the '
                    'parent\'s sea, and dry cells below a 4-neighbouring lake surface (empty before S2.8). Phase S3\'s '
                    '(empty before S3.0): the cells the limited prolongation clamped and its largest move (S3.1), the '
                    'smoothing\'s 99th-percentile move and the largest rise of a dry channel cell above its bed (S3.2), '
                    'the land and channel cells, the drainage density (channel cells over land cells times the cell size) '
                    'and the tributaries (S3.3), the lake edges and their share on global cell boundaries (about 1 / k for '
                    'an edge that ignores the global grid) and the land\'s RMS Laplacian at the level\'s cell size (S3.4), '
                    'the cells the valley floors, troughs and cirques changed (S3.5).'}


def run_data(evaluation_dir):
    directory = Path(evaluation_dir).absolute()
    document = read_json(directory / 'evaluation.json')
    if document is None:
        raise SystemExit(f'{directory} has no evaluation.json')
    run = read_json(directory.parent / 'run.json') if directory.name == 'evaluation' else None
    reviews = review.load(directory)['entries']
    verdict = dict(document['verdict'])
    if reviews:
        verdict['review'] = review.latest(directory)['verdict']  # the user's latest over an agent's
    groups, animations = map_groups(document['maps'])
    files = [{'title': 'evaluation.json', 'path': 'evaluation.json'}, {'title': 'report.md', 'path': 'report.md'}]
    files += [{'title': name, 'path': name} for name in ('review.json',) if (directory / name).is_file()]
    if directory.name == 'evaluation':
        files += [{'title': name, 'path': f'../{name}'} for name in ('run.json', 'log.txt', 'log-repeat.txt', 'log-explicit.txt', 'run.log', 'repeat.log')
                  if (directory.parent / name).is_file()]
    validation = document['validation']
    return {'schema': 'uw-run-page/v1', 'generated': now(), 'module': document['profile']['module'],
            'name': run_identity(directory, run), 'facts': run_facts(document, run), 'verdict': verdict,
            'evidence': document['evidence'], 'checks': document['checks'], 'scores': score_lines(validation),
            'rows': validation_rows(validation), 'has_baseline': document.get('baseline') is not None,
            'distributions': distribution_tables(validation), 'regions': regional_tables(validation),
            'features': validation.get('features'), 'footprint': validation.get('footprint'),
            'skipped': validation['skipped'], 'reviews': reviews,
            'groups': groups, 'animations': animations, 'files': files, 'earth': document.get('earth'),
            'earth_text': earth_text(document['earth']) if document.get('earth') else None,
            'levels': level_rows(document.get('levels')),
            'provenance': {'profile': document['profile'], 'code': document['code'], 'references': document['references']}}


def markdown(data):
    v = data['verdict']
    lines = [f"# Evaluation · {data['module']} · {data['name']}", '', *[f'- {label}: {value}' for label, value in data['facts']],
             '', f"- Checks: {v['checks']}" + (f" ({', '.join(v['failed_checks'])})" if v['failed_checks'] else ''),
             f"- Validation: {v['validation']}", f"- Review: {v['review']}", '',
             f"Evidence weighting: {json.dumps(data['evidence'])}", '', '## Checks', '',
             '| Check | Status | Summary |', '|---|---|---|']
    lines += [f"| {item['check']} | {item['status']} | {item['summary']} |" for item in data['checks']]
    for item in data['checks']:
        if item['details'] and item['status'] in ('fail', 'reported'):
            lines += ['', f"{item['check']}:", *[f'- {detail}' for detail in item['details'][:20]]]
    if data.get('levels'):
        table = data['levels']
        lines += ['', '## Refinement levels', '', table['note'], '', '| ' + ' | '.join(table['columns']) + ' |',
                  '|---|' + '---|' * (len(table['columns']) - 1)]
        lines += ['| ' + ' | '.join(fmt(value) for value in row) + ' |' for row in table['rows']]
    if data['rows'] or data['distributions']:
        lines += ['', '## Validation', '']
        for entry in data['scores']:
            if entry['score'] is not None:
                lines.append(f"- {entry['key']}: {entry['version']} on {entry['inputs']} = {entry['score']:.4f} (lower is better). {entry['note']}")
            else:
                lines.append(f"- {entry['key']}: partial score {fmt(entry['partial_score'])} over weight "
                             f"{entry['covered_weight']:.2f}; {entry['note']}")
        if data['rows']:
            lines += ['', '| Field | Season | Loss | Baseline loss | Delta | RMSE | Bias | Correlation |', '|---|---|---|---|---|---|---|---|']
            lines += [f"| {row['field']} | {row['season']} | {fmt(row['loss'])} | {fmt(row['baseline'])} | {fmt(row['delta'])} | "
                      f"{fmt(row['rmse'])} | {fmt(row['bias'])} | {fmt(row['correlation'])} |" for row in data['rows']]
        for table in data['distributions']:
            lines += ['', f"### {table['name']} ({table['reference']})", '', '| Metric | ' + ' | '.join(table['columns']) + ' |',
                      '|---|' + '---|' * len(table['columns'])]
            lines += [f"| {row['metric']} | " + ' | '.join(fmt(row.get(key)) for key in table['columns']) + ' |' for row in table['rows']]
        lines += [f"- skipped {item['id']}: {item['reason']}" for item in data['skipped']]
    if data.get('regions'):
        lines += ['', '### Regional errors (project boxes)', '',
                  f"{data['regions']['note']} Run bias / RMSE of model minus reference; AR6 regions and the baseline "
                  'are in evaluation.json and on the run page.']
        fields = {}
        for option in data['regions']['options']:
            fields.setdefault(option['field'], []).append(option)
        for field, options in fields.items():
            regions = list(dict.fromkeys(row[0] for option in options for row in option['sets'].get('project', [])))
            values = [{row[0]: f'{fmt(row[2], 3)} / {fmt(row[3], 3)}' for row in option['sets'].get('project', [])} for option in options]
            lines += ['', f'| {field} | ' + ' | '.join(option['season'] for option in options) + ' |',
                      '|---|' + '---|' * len(options)]
            lines += [f"| {data['regions']['names']['project'].get(region, region)} | " +
                      ' | '.join(value.get(region, '–') for value in values) + ' |' for region in regions]
    if data['reviews']:
        lines += ['', '## Review', '']
        lines += [f"- {entry['recorded']} {entry['reviewer']}: **{entry['verdict']}**" +
                  (f" against {entry['against']}" if entry.get('against') else '') + (f". {entry['notes']}" if entry.get('notes') else '')
                  for entry in data['reviews']]
    if data['groups'] or data['animations']:
        lines += ['', '## Maps', '']
        for group in data['groups']:
            lines.append(f"- {group['title']}: " + ' · '.join(f"[{kind}]({group['images'][kind]['path']})" for kind in KINDS
                                                              if kind in group['images']))
        lines += [f"- {item['title']}: [APNG]({item['path']})" for item in data['animations']]
    if data.get('earth'):
        lines += ['', earth_text(data['earth'])]
    provenance = data['provenance']
    lines += ['', '## Provenance', '', f"- Profile: `{provenance['profile']['path']}` ({provenance['profile']['sha256'][:16]})",
              f"- Code: {provenance['code']['revision']} (dirty: {provenance['code']['dirty']})",
              *[f'- Reference `{name}`: {digest[:16]}' for name, digest in provenance['references'].items()]]
    return '\n'.join(lines) + '\n'


def earth_text(earth):
    """One line on a run's Earth reference window: what it is and the features to look for, or why there is none."""
    if not earth.get('entry') or 'caption' not in earth:
        return earth['note']
    features = '; '.join(earth.get('features') or [])
    return (f"{earth['caption']}. {earth['note']} {earth['pixels'][0]} × {earth['pixels'][1]} pixels "
            f"(the map: {earth['map_pixels'][0]} × {earth['map_pixels'][1]})." + (f' Features to look for: {features}.' if features else ''))


def run_page(evaluation_dir):
    """Write index.html and report.md of an evaluation directory from its JSON files."""
    directory = runs.refuse_archive(evaluation_dir)
    data = run_data(directory)
    (directory / 'report.md').write_text(markdown(data), encoding='utf-8')
    return write_page('run.html', data, directory / 'index.html')


def load_profile_maps(module):
    profile = read_json(PROFILES / f'{module}.json') or {}
    return {spec['name']: spec for spec in profile.get('maps', [])}


def season_order(season):
    return -1 if season is None else SEASONS.index(season) if season in SEASONS else len(SEASONS)


def comparison_earth(runs, module, output, name):
    """(map entries, earth block) of catalogue entry `name` drawn for the first run's regional map into output."""
    import earth_windows
    from uwtools.bundle import Bundle
    profile = read_json(PROFILES / f'{module}.json') or {}
    if not profile.get('earth'):
        raise SystemExit(f'--earth: the {module} profile draws no Earth windows')
    bundles = []
    for _, document, _ in runs:
        path = Path(document['run']['bundle'])
        bundles.append(Bundle.open(path if path.is_absolute() else ROOT / path))
    if any(bundle.manifest.get('window') is None for bundle in bundles):
        raise SystemExit('--earth: every run of the comparison must be a regional map (a bundle with a window block)')
    metres = [earth_windows.map_request(bundle.manifest)['metres_per_pixel'] for bundle in bundles]
    if max(metres) - min(metres) > 1e-6 * max(metres):
        raise SystemExit('--earth: the runs have different metres per pixel (' + ', '.join(f'{value:,.0f}' for value in metres) +
                         '); each keeps the Earth window of its evaluation (evaluate --earth NAME)')
    try:
        return earth_windows.render_maps(bundles[0], profile, output, name, tag=f'-{name}')
    except ValueError as error:
        raise SystemExit(f'--earth: {error}')


def comparison_data(evaluation_dirs, output_dir, earth=None):
    if len(evaluation_dirs) < 2:
        raise SystemExit('A comparison needs at least two evaluation directories (the first is the baseline)')
    output = Path(output_dir).absolute()
    runs = []
    for item in evaluation_dirs:
        directory = Path(item).absolute()
        document = read_json(directory / 'evaluation.json')
        if document is None:
            raise SystemExit(f'{directory} has no evaluation.json')
        run = read_json(directory.parent / 'run.json') if directory.name == 'evaluation' else None
        runs.append((directory, document, run))
    modules = {document['profile']['module'] for _, document, _ in runs}
    if len(modules) != 1:
        raise SystemExit(f"Evaluations of different modules: {', '.join(sorted(modules))}")
    module = modules.pop()
    names, header = [], []
    for directory, document, run in runs:
        name = run_identity(directory, run)
        name = name if name not in names else f'{name} ({len(names) + 1})'
        names.append(name)
        counted = review.latest(directory)
        verdict = dict(document['verdict'], review=counted['verdict'] if counted else document['verdict']['review'])
        page = directory / 'index.html'
        header.append({'name': name, 'bundle': document['run']['bundle'], 'seed': document['run']['seed'],
                       'grid': document['run']['grid'], 'revision': (document['code'].get('revision') or '')[:12],
                       'dirty': document['code'].get('dirty'), 'verdict': verdict,
                       'page': link(page, output) if page.is_file() else None})

    def row(metric, values, lower_is_better=True):
        base = values[0]
        deltas = [None if value is None or base is None or not isinstance(value, (int, float)) else value - base for value in values]
        return {'metric': metric, 'values': values, 'deltas': [None, *deltas[1:]], 'lower_is_better': lower_is_better}

    metrics = [row('failed checks', [len(document['verdict']['failed_checks']) for _, document, _ in runs])]
    if any(document['validation']['scores'] for _, document, _ in runs):
        versions = {entry['version'] for _, document, _ in runs for entry in document['validation']['scores'].values()}
        metrics.append(row(f"score ({', '.join(sorted(versions))}; partial where coverage < 1)",
                           [score_value(document['validation']['scores'].get('run')) for _, document, _ in runs]))
        metrics.append(row('covered weight', [document['validation']['scores'].get('run', {}).get('covered_weight')
                                              for _, document, _ in runs], None))
    fields = list(dict.fromkeys(field for _, document, _ in runs for field in document['validation']['maps']))
    metrics += [row(f'{field} loss', [field_loss(document['validation'], field) for _, document, _ in runs]) for field in fields]
    for name in dict.fromkeys(name for _, document, _ in runs for name in document['validation']['distributions']):
        entries = [document['validation']['distributions'].get(name, {}).get('run', {}) for _, document, _ in runs]
        for group in ('earth_comparison', 'metrics'):
            keys = list(dict.fromkeys(key for entry in entries for key, value in entry.get(group, {}).items()
                                      if not isinstance(value, dict)))
            metrics += [row(f'{name}: {key}', [entry.get(group, {}).get(key) for entry in entries], None) for key in keys]

    specs = load_profile_maps(module)
    order = list(specs)
    options, index = [], {}
    # Earth windows: the one named by --earth (drawn here for the first run), else each run's own, one column
    # per distinct window (by caption, which names the entry and the metres per pixel, and pixel count).
    earth_note = None
    sources = [(directory, document['maps']) for directory, document, _ in runs]
    if earth is not None:
        drawn, block = comparison_earth(runs, module, output, earth)
        earth_note = block['note'] if block else None
        sources = [(directory, [item for item in maps if item.get('kind') != 'earth']) for directory, maps in sources]
        sources.append((output, drawn))
    else:
        notes = [f"{name}: {document['earth']['note']}" for name, (_, document, _) in zip(names, runs) if document.get('earth')]
        earth_note = ' '.join(notes) or None
    for position, (directory, items) in enumerate(sources):
        for item in items:
            key = (item.get('map', item['title']), item.get('season'))
            if key not in index:
                spec = specs.get(key[0], {})
                index[key] = {'map': key[0], 'season': key[1], 'season_label': item.get('season_label', key[1]),
                              'title': spec.get('title', item['title']),
                              'runs': [{} for _ in runs], 'reference': None, 'earth': []}
                options.append(index[key])
            kind = item.get('kind', 'model')
            path = link(directory / item['path'], output)
            run_name = names[position] if position < len(names) else 'all runs (--earth)'
            if kind == 'earth':
                shown = next((image for image in index[key]['earth']
                              if (image['title'], image['pixels']) == (item['title'], item.get('pixels'))), None)
                if shown is None:
                    index[key]['earth'].append({'path': path, 'title': item['title'], 'pixels': item.get('pixels'),
                                                'display_fraction': item.get('display_fraction'), 'runs': [run_name]})
                else:
                    shown['runs'].append(run_name)
            elif kind == 'reference':
                if index[key]['reference'] is None:
                    index[key]['reference'] = {'path': path, 'title': item['title'], 'run': names[position]}
            else:
                index[key]['runs'][position][kind] = {'path': path, 'title': item['title'], 'half_range': item.get('half_range'),
                                                      'frames': len(item.get('frames', []))}
    options.sort(key=lambda item: (order.index(item['map']) if item['map'] in order else len(order), season_order(item['season'])))
    return {'schema': 'uw-comparison-page/v1', 'generated': now(), 'module': module, 'runs': header, 'metrics': metrics,
            'options': options, 'note': 'The first run is the baseline; deltas are against it (lower loss is better).',
            'earth_note': earth_note}


def comparison_markdown(data):
    names = [run['name'] for run in data['runs']]
    lines = [f"# Comparison · {data['module']}", '', data['note'], '',
             '| Metric | ' + ' | '.join(names) + ' |', '|---|' + '---|' * len(names)]
    for row in data['metrics']:
        cells = [fmt(row['values'][0])] + [f"{fmt(value)} ({'+' if (delta or 0) > 0 else ''}{fmt(delta)})" if delta is not None else fmt(value)
                                           for value, delta in zip(row['values'][1:], row['deltas'][1:])]
        lines.append(f"| {row['metric']} | " + ' | '.join(cells) + ' |')
    lines += ['', '## Runs', '', *[f"- {run['name']}: `{run['bundle']}` (seed {run['seed']}, {run['grid'][0]} x {run['grid'][1]}, "
                                   f"code {run['revision']}{' dirty' if run['dirty'] else ''}); checks {run['verdict']['checks']}, "
                                   f"review {run['verdict']['review']}" for run in data['runs']]]
    return '\n'.join(lines) + '\n'


def comparison_page(evaluation_dirs, output_dir, earth=None):
    """Write index.html and report.md comparing evaluations (the first is the baseline) into output_dir;
    earth: a catalogue entry drawn beside every run's regional maps instead of their own Earth windows."""
    output = runs.refuse_archive(output_dir)
    data = comparison_data(evaluation_dirs, output, earth)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'report.md').write_text(comparison_markdown(data), encoding='utf-8')
    return write_page('comparison.html', data, output / 'index.html')


def climate_rankings():
    """Run id -> (global-balanced-v1 score, report directory) from the registered ranking batches, latest last."""
    scores = {}
    for path in sorted(ROOT.glob('runs/reports/*/ranking.json'), key=lambda item: item.stat().st_mtime):
        ranking = read_json(path)
        for identifier, value in ranking.get('scores', {}).items():
            scores[int(identifier)] = (value, path.parent)
    return scores


def registered_run(record):
    """A record of the current layout (it names its run directory): run.json and the registry entry."""
    run = read_json(ROOT / record['run'] / 'run.json') or {}
    grid = run.get('grid')
    return {'id': record['id'], 'date': (record.get('date') or '')[:16].replace('T', ' '), 'status': record.get('status'),
            'label': record.get('label', ''), 'grid': f'{grid[0]} × {grid[1]}' if grid else None,
            'score': record.get('score'), 'metrics': {'elapsed_s': run.get('seconds')},
            'links': [('run.json', ROOT / record['run'] / 'run.json'), ('log', ROOT / record['run'] / 'log.txt')]}


def climate_run(record, rankings):
    if record.get('run'):
        return registered_run(record)
    identifier, metrics = record['id'], record.get('metrics') or {}
    stamp = str(record.get('datetime', ''))
    date = f'{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]} {stamp[8:10]}:{stamp[10:12]}' if len(stamp) >= 12 else stamp
    score, source = rankings.get(identifier, (None, None))
    links = [('run page', ROOT / f'runs/reports/{identifier}.html'), ('record', ROOT / f'runs/reports/{identifier}.md'),
             ('gallery', ROOT / f'runs/reports/{identifier}-gallery.html')]
    if source is not None:
        links += [('ranking', next(iter(sorted(source.parent.glob(f'{source.name}*.html'))), source / 'ranking.json'))]
    grid = record.get('horizontal_resolution')
    return {'id': identifier, 'date': date, 'status': 'registered' if metrics else 'no metrics', 'label': record.get('information', ''),
            'grid': f'{grid} × {grid // 2}' if grid else None, 'score': score,
            'metrics': {'weighted_relative_error': metrics.get('weighted_relative_error'),
                        'land_precipitation_correlation': metrics.get('area_weighted_land_precipitation_correlation')},
            'links': links}


def tectonics_run(record):
    identifier = record['id']
    metrics_path = runs.derived('tectonics', identifier, 'metrics.json')
    final = (read_json(metrics_path) or {}).get('final', {}) if metrics_path else {}
    comparison = final.get('earth_comparison', {})
    if record.get('run'):
        entry = registered_run(record)
        entry['metrics'] = {'land_fraction': final.get('land_fraction'),
                            'hypsometry_gap_m': comparison.get('hypsography_mean_absolute_quantile_gap_m'), **entry['metrics']}
        return entry
    manifest = read_json(runs.archived(record['manifest'])) if record.get('manifest') else None
    manifest = manifest or {}
    recipe = manifest.get('recipe', {})
    links = [('run page', ROOT / f'runs/reports/tectonics/{identifier}.html'), ('record', ROOT / f'runs/reports/tectonics/{identifier}.md'),
             ('manifest', runs.archived(record['manifest']) if record.get('manifest') else ROOT / 'missing')]
    label = recipe.get('name', '') + (f" ({record['kind']})" if record.get('kind') else '')
    return {'id': identifier, 'date': (manifest.get('started_utc') or '')[:16].replace('T', ' '), 'status': record.get('status'),
            'label': label, 'grid': f"{recipe['width']} × {recipe['height']}" if 'width' in recipe else None, 'score': None,
            'metrics': {'land_fraction': final.get('land_fraction'),
                        'hypsometry_gap_m': comparison.get('hypsography_mean_absolute_quantile_gap_m'),
                        'elapsed_s': manifest.get('elapsed_seconds')},
            'links': links}


def dashboard_data(module, output_dir):
    output = Path(output_dir).absolute()
    registry = read_json(ROOT / f'runs/registry/{module}.json') or {'runs': []}
    rankings = climate_rankings() if module == 'climate' else {}
    runs = []
    for record in registry['runs']:
        entry = climate_run(record, rankings) if module == 'climate' else tectonics_run(record) if module == 'tectonics' else \
            {'id': record['id'], 'date': None, 'status': record.get('status'), 'label': '', 'grid': None, 'score': None,
             'metrics': {}, 'links': []}
        evaluation = ROOT / f"runs/{module}/{record['id']}/evaluation"
        document = read_json(evaluation / 'evaluation.json')
        entry['evaluation_score'] = score_value(document['validation']['scores'].get('run')) if document else None
        links = [('evaluation', evaluation / 'index.html')] + entry['links']
        entry['links'] = [{'title': title, 'href': link(path, output)} for title, path in links if Path(path).exists()]
        runs.append(entry)
    runs.sort(key=lambda item: item['id'], reverse=True)
    columns = list(dict.fromkeys(key for entry in runs for key in entry['metrics']))
    entries = [{'id': entry['id'], 'date': entry['date'], 'verdict': entry['verdict'], 'decision': entry['decision'],
                'what': entry.get('phase') or entry['hypothesis'], 'hypothesis': entry['hypothesis'],
                'do_not_retry': entry['do_not_retry'], 'lessons': entry['lessons']} for entry in ledger.entries(module)]
    score_note = {'climate': 'Score: global-balanced-v1 of the latest registered ranking batch that contains the run '
                             '(runs/reports/*/ranking.json). Evaluation score: bundle evaluation; the two are not interchangeable.',
                  'tectonics': 'Tectonics runs are judged by visual review; Earth hypsometry is context, not a target.'}
    return {'schema': 'uw-dashboard-page/v1', 'generated': now(), 'module': module, 'baseline': registry.get('baseline'),
            'registry': f'runs/registry/{module}.json', 'columns': columns, 'runs': runs, 'ledger': entries,
            'note': score_note.get(module, '')}


def dashboard(module, output_dir):
    """Write the module dashboard (registry runs and ledger entries) to output_dir/index.html."""
    return write_page('dashboard.html', dashboard_data(module, output_dir), Path(output_dir) / 'index.html')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    run = commands.add_parser('run', help='Rewrite index.html and report.md of an evaluation directory')
    run.add_argument('evaluation', type=Path)
    compare = commands.add_parser('compare', help='Comparison page; the first evaluation is the baseline')
    compare.add_argument('evaluations', type=Path, nargs='+')
    compare.add_argument('--output', type=Path, required=True)
    compare.add_argument('--earth', metavar='NAME', help='Earth reference window (catalogue entry) drawn for the first run')
    board = commands.add_parser('dashboard', help='Module dashboard over the registry and the ledger')
    board.add_argument('module')
    board.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'run':
        path = run_page(args.evaluation)
    elif args.command == 'compare':
        path = comparison_page(args.evaluations, args.output, args.earth)
    else:
        path = dashboard(args.module, args.output)
    print(f'Page: {path}')


if __name__ == '__main__':
    main()
