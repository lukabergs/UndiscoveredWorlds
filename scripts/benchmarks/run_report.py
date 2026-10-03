# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Standard per-run page for every simulation module.

Each module supplies facts, links, exported maps and optional animations; this
module writes projected map views and one HTML page from run_report_template.html.

    uv run --offline scripts/benchmarks/run_report.py climate --run 382 [--compare 381]
    uv run --offline scripts/benchmarks/run_report.py tectonics --run 75 [--compare 74]
"""
import argparse
import html
import json
import os
from pathlib import Path
import re

from map_projection import PROJECTIONS, save_projections

ROOT = Path(__file__).resolve().parents[2]
MAPS = ROOT / 'runs/maps'
TEMPLATE = Path(__file__).with_name('run_report_template.html')
SEASONS = {'1': 'January', '2': 'April', '3': 'July', '4': 'October', 'annual': 'annual'}
LAYERS = {'s': 'surface', 'u': 'upper'}


def relative(path, parent):
    return Path(os.path.relpath(path, parent)).as_posix()


def write_page(module, run_id, title, facts, links, maps, output, compare_run=None, animations=(),
               map_note=''):
    """maps: [(title, path)] for this run; the comparison uses the same path with the other run ID."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    entries = []
    for map_title, path in maps:
        views = {}
        own = save_projections(path, MAPS)
        other = path.with_name(f'{compare_run}{path.suffix}') if compare_run is not None else None
        others = save_projections(other, MAPS) if other is not None and other.exists() else {}
        for name in PROJECTIONS:
            views[name] = {'run': relative(own[name], output.parent),
                           'compare': relative(others[name], output.parent) if name in others else None}
        entries.append({'title': map_title, 'views': views})
    data = {'module': module, 'run': run_id, 'compare_run': compare_run, 'facts': facts,
            'links': [[label, relative(target, output.parent) if isinstance(target, Path) else target]
                      for label, target in links],
            'maps': entries, 'projections': PROJECTIONS, 'map_note': map_note,
            'animations': [{'title': a_title, 'href': relative(path, output.parent)} for a_title, path in animations]}
    encoded = json.dumps(data, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')
    output.write_text(TEMPLATE.read_text(encoding='utf-8').replace('@@TITLE@@', html.escape(title))
                      .replace('@@DATA@@', encoded), encoding='utf-8')
    print(f'Run page: {output}')
    return output


def climate_map_title(path):
    parts = list(path.relative_to(MAPS / 'climate').parts[:-1])
    words = []
    for part in parts:
        words.append(SEASONS.get(part, LAYERS.get(part, part.replace('_', ' '))))
    return ' · '.join(words)


def climate_page(run_id, compare_run=None):
    registry = json.loads((ROOT / 'runs/registry/climate.json').read_text(encoding='utf-8-sig'))['runs']
    record = next(r for r in registry if r['id'] == run_id)
    width = record.get('horizontal_resolution')
    if compare_run is None:
        earlier = [r['id'] for r in registry if r['id'] < run_id and r.get('horizontal_resolution') == width]
        compare_run = max(earlier) if earlier else None
    log = ROOT / f'runs/logs/climate/{run_id}.log'
    text = log.read_text(encoding='utf-8', errors='replace') if log.exists() else ''
    maps = []
    for line in text.splitlines():
        if line.startswith('Climate benchmark map: ') and line.endswith('.png'):
            path = ROOT / line.removeprefix('Climate benchmark map: ').replace('\\', '/')
            maps.append((climate_map_title(path), path))
    receipt = ROOT / f'runs/reports/{run_id}-record.json'
    elapsed = json.loads(receipt.read_text()).get('elapsed_seconds') if receipt.exists() else None
    drag = re.findall(r'Prescribed surface drag input=(\S+)', text)
    facts = [['Module', 'climate (Earth benchmark)'], ['Run', run_id], ['Compared with', compare_run or 'none'],
             ['World grid', f'{width} × {width // 2}' if width else 'unrecorded'],
             ['Configuration', record.get('information', '')],
             ['Surface drag input', Path(drag[-1]).name if drag else 'uniform land drag'],
             ['Status', 'Not converged' if 'Not converged' in text else 'see log']]
    if elapsed:
        facts.append(['Elapsed', f'{float(elapsed):.1f} s'])
    reports = ROOT / 'runs/reports'
    links = [['Run record', reports / f'{run_id}.md'], ['Raw log', log],
             ['Diagnostics', ROOT / f'runs/diagnostics/climate/{run_id}'],
             ['Latest batch comparison', reports / 'climate-latest-comparison.html']]
    if (reports / f'{run_id}-gallery.html').exists():
        links.insert(1, ['Reference gallery', reports / f'{run_id}-gallery.html'])
    return write_page('climate', run_id, f'Climate run {run_id}', facts, links, maps, reports / f'{run_id}.html',
                      compare_run, map_note='Maps have the world grid size. Coarser atmosphere, ocean and moisture '
                      'grids are remapped onto it. Black is masked or missing data.')


def tectonics_page(run_id, compare_run=None):
    manifest = json.loads((ROOT / f'runs/manifests/tectonics/{run_id}.json').read_text())
    if compare_run is None:
        registry = json.loads((ROOT / 'runs/registry/tectonics.json').read_text())['runs']
        earlier = [r['id'] for r in registry if r['id'] < run_id and r.get('status') == 'complete']
        compare_run = max(earlier) if earlier else None
    fields = json.loads((ROOT / 'configs/tectonics-palettes.json').read_text())['fields']
    maps, animations = [], []
    for label in ('final', 'initial'):
        for field, spec in fields.items():
            path = MAPS / f'tectonics/{field}/{label}/{run_id}.png'
            if path.exists():
                maps.append((f"{spec['title']} · {label}", path))
    for field, spec in fields.items():
        path = MAPS / f'tectonics/{field}/animation/{run_id}.png'
        if path.exists():
            animations.append((spec['title'], path))
    recipe = manifest['recipe']
    scenario = recipe.get('scenario', {})
    final = next((s for s in manifest.get('snapshots', []) if s['label'] == 'final'), {})
    facts = [['Module', 'tectonics'], ['Run', run_id], ['Compared with', compare_run or 'none'],
             ['Recipe', recipe.get('name', '')], ['Seed', recipe['seed']],
             ['Native grid', f"{recipe['width']} × {recipe['height']} (planar torus)"],
             ['Plates', scenario.get('plates', '')], ['Duration', f"{scenario.get('evolution-duration-myr', '')} Myr"],
             ['Final time', f"{final.get('time_myr', '')} Myr"], ['Status', manifest['status']]]
    if manifest.get('elapsed_seconds') is not None:
        facts.append(['Elapsed', f"{manifest['elapsed_seconds']:.1f} s"])
    reports = ROOT / 'runs/reports/tectonics'
    links = [['Run record', reports / f'{run_id}.md'], ['Raw log', ROOT / f'runs/logs/tectonics/{run_id}.log'],
             ['Manifest', ROOT / f'runs/manifests/tectonics/{run_id}.json'],
             ['Metrics', ROOT / f'runs/diagnostics/tectonics/{run_id}/metrics.json']]
    if compare_run is not None and (reports / f'{compare_run}-{run_id}.html').exists():
        links.insert(1, ['Detailed comparison', reports / f'{compare_run}-{run_id}.html'])
    return write_page('tectonics', run_id, f'Tectonics run {run_id}', facts, links, maps, reports / f'{run_id}.html',
                      compare_run, animations, map_note='The tectonics domain is a planar torus. Projected views '
                      'wrap it onto the globe the way the application ingests it; poles are not physical poles '
                      'of the native simulation.')


PAGES = {'climate': climate_page, 'tectonics': tectonics_page}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('module', choices=PAGES)
    parser.add_argument('--run', type=int, required=True)
    parser.add_argument('--compare', type=int)
    args = parser.parse_args()
    PAGES[args.module](args.run, args.compare)


if __name__ == '__main__':
    main()
