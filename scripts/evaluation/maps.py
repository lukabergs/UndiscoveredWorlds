"""Review maps of bundle fields in the project palettes (climate and tectonics)."""
import colorsys
import json
import sys
from pathlib import Path

import numpy as np

from bundle import SEASONS
from validation import load_selected, surface_elevation

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/benchmarks'))
import map_rendering  # noqa: E402
from climate_palettes import field_rgb  # noqa: E402

TECTONICS = json.loads((ROOT / 'configs/tectonics-palettes.json').read_text(encoding='utf-8'))
DERIVED = {'surface_elevation': surface_elevation}


def tectonics_rgb(values, name):
    """The rule of tectonics_report.py: identities, labelled classes or interpolated stops."""
    spec = TECTONICS['fields'][name]
    if spec.get('ids'):
        colours = {int(identity): [round(c * 255) for c in colorsys.hsv_to_rgb((int(identity) * .61803398875) % 1, .52, .92)]
                   for identity in np.unique(values[np.isfinite(values)])}
        return map_rendering.categories(values, colours, TECTONICS['missing'])
    if 'colours' in spec:
        return map_rendering.categories(values, dict(enumerate(spec['colours'])), TECTONICS['missing'])
    if 'stops' in spec:
        stops = spec['stops']
    else:
        low, high = spec['range']
        stops = [(low + t * (high - low), rgb) for t, rgb in TECTONICS['sequential']]
    return map_rendering.colourize(values, stops, TECTONICS['missing'])


def palette_rgb(values, palette):
    family, _, name = palette.partition(':')
    if family == 'climate':
        return field_rgb(values, name)
    if family == 'tectonics':
        return tectonics_rgb(values, name)
    raise ValueError(f'Unknown palette {palette!r}; use climate:<scale> or tectonics:<field>')


def map_values(bundle, spec):
    if 'derived' in spec:
        values = DERIVED[spec['derived']](bundle)
    elif 'vector' in spec:
        east, north = (load_selected(bundle, selector) for selector in spec['vector'])
        values = np.hypot(east, north)
    else:
        values = load_selected(bundle, spec['field'])
    if spec.get('reduce') == 'season-mean':
        values = values.mean(axis=0)
    return values * spec.get('scale', 1)


def available(bundle, spec):
    selectors = spec.get('vector', [spec['field']] if 'field' in spec else [])
    if spec.get('derived') == 'surface_elevation':
        selectors = ['geosphere.base_elevation', 'geosphere.mountain_height', 'geosphere.volcano']
    return all(selector.split('[')[0] in bundle.fields for selector in selectors)


def render(bundle, specs, directory):
    """Write PNGs under directory/maps; returns [{'title', 'path', 'season'}] relative to directory."""
    written = []
    for spec in specs:
        if not available(bundle, spec):
            continue
        values = map_values(bundle, spec).astype(float)
        if spec.get('mask') in ('sea', 'land'):
            sea = bundle.sea_fraction(*values.shape[-2:]) > 0.9
            values = np.where(sea if spec['mask'] == 'sea' else ~sea, values, np.nan)
        layers = [(None, values)] if values.ndim == 2 else [(season, values[season]) for season in spec.get('seasons', range(values.shape[0]))]
        for season, layer in layers:
            suffix = '' if season is None else f'-{SEASONS[season].lower()}'
            relative = Path('maps') / f"{spec['name']}{suffix}.png"
            map_rendering.save_rgb(directory / relative, palette_rgb(layer, spec['palette']))
            title = spec['title'] + ('' if season is None else f' · {SEASONS[season]}')
            written.append({'title': title, 'path': relative.as_posix(), 'palette': spec['palette']})
    return written
