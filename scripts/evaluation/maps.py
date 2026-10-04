"""Review maps of bundle fields in the project palettes (climate and tectonics).

A profile map with "validation": "<id>" is drawn from that validation entry when it was
compared: model and reference reduced by validation.py's season rule, conversion and
remap, and model minus reference on the score's mask, so the pictures show what the loss
measures. Vector fields are drawn as speed; their difference is the speed difference.
Differences use the diverging scale of configs/climate-palettes.json ("difference").
Runs with checkpoint bundles get one APNG per map and season.
"""
import colorsys
import io
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import palettes
from uwtools.bundle import SEASONS
from validation import field_mask, load_selected, model_field, reference_field, seasons_for, surface_elevation

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/benchmarks'))
import map_rendering  # noqa: E402
from climate_palettes import DEFINITION as CLIMATE, field_rgb  # noqa: E402

TECTONICS = palettes.tectonics()
DERIVED = {'surface_elevation': surface_elevation}
QUARTERS = ('Jan–Mar', 'Apr–Jun', 'Jul–Sep', 'Oct–Dec')
SNAPSHOT_TABLE = 'geosphere.tectonic_snapshot'


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
    if family == 'climate-difference':
        return difference_rgb(values, name)
    if family == 'tectonics':
        return tectonics_rgb(values, name)
    raise ValueError(f'Unknown palette {palette!r}; use climate:<scale> or tectonics:<field>')


def difference_rgb(values, scale):
    """Diverging colours symmetric about zero over +-half_range of the scale."""
    spec = CLIMATE['difference']
    half = float(spec['half_ranges'][scale])
    stops = [(-half + 2 * half * t, colour) for t, colour in spec['stops']]
    return map_rendering.colourize(values, stops, CLIMATE['missing'])


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


def masked(bundle, spec, values):
    if spec.get('mask') in ('sea', 'land'):
        sea = bundle.sea_fraction(*values.shape[-2:]) > 0.9
        values = np.where(sea if spec['mask'] == 'sea' else ~sea, values, np.nan)
    return values


def model_layers(bundle, spec):
    """[(season or None, values)] of a profile map drawn from its own field."""
    values = masked(bundle, spec, map_values(bundle, spec).astype(float))
    if values.ndim == 2:
        return [(None, values)]
    return [(season, values[season]) for season in spec.get('seasons', range(values.shape[0]))]


def speed(values):
    return np.hypot(*values) if values.ndim == 3 else values


def season_title(rule, season):
    if season is None or rule == 'annual':
        return None
    return QUARTERS[season] if rule == 'quarter' else SEASONS[season]


def save(directory, name, rgb):
    relative = Path('maps') / f'{name}.png'
    map_rendering.save_rgb(directory / relative, rgb)
    return relative.as_posix()


def entry(spec, kind, season, title, path, palette, **extra):
    """season: the model season index (None when annual); season_label (in extra) names the months shown."""
    name = None if season is None else SEASONS[season]
    return {'map': spec['name'], 'kind': kind, 'season': name, 'season_label': name, 'title': title, 'path': path,
            'palette': palette, **extra}


def validated(bundle, spec, validation, references, directory):
    """Model, reference and difference images of a map linked to a compared validation entry."""
    scale = validation.get('metric_scale', validation['id'])
    written = []
    seasons = seasons_for(validation)
    if 'seasons' in spec and len(seasons) > 1:
        seasons = [season for season in seasons if season in spec['seasons']]
    for season in seasons:
        model = speed(model_field(bundle, validation, season))
        reference = speed(reference_field(validation, references, season, model.shape[-1]))
        mask = np.isfinite(reference) & np.isfinite(model) & field_mask(bundle, validation, season, model.shape)
        annual = validation['season'] == 'annual' or validation.get('reduce') == 'season-mean'
        layer = None if annual else season
        label = season_title(validation['season'], layer)
        title = spec['title'] + ('' if label is None else f' · {label}')
        suffix = '' if layer is None else f'-{SEASONS[layer].lower()}'
        source = ', '.join(validation['reference'])
        common = {'validation': validation['id'], 'season_label': label}
        written.append(entry(spec, 'model', layer, title,
                             save(directory, spec['name'] + suffix, palette_rgb(masked(bundle, spec, model), spec['palette'])),
                             spec['palette'], **common))
        written.append(entry(spec, 'reference', layer, f'{title} · reference ({source})',
                             save(directory, f"{spec['name']}{suffix}-reference", palette_rgb(reference, spec['palette'])),
                             spec['palette'], **common))
        difference = np.where(mask, model - reference, np.nan)
        palette = f'climate-difference:{scale}'
        written.append(entry(spec, 'difference', layer, f'{title} · model − reference',
                             save(directory, f"{spec['name']}{suffix}-difference", palette_rgb(difference, palette)),
                             palette, **common, half_range=CLIMATE['difference']['half_ranges'][scale]))
    return written


def checkpoint_time(bundle):
    if SNAPSHOT_TABLE in bundle.tables:
        return json.loads(bundle.table_text(SNAPSHOT_TABLE)).get('time_myr')
    return None


def animations(bundle, specs, directory, checkpoints):
    """One APNG per map and season over the checkpoint bundles and the run bundle."""
    frames = [item for item in checkpoints if item.path.resolve() != bundle.path.resolve()]
    times = [checkpoint_time(item) for item in frames]
    final = checkpoint_time(bundle)
    # The end of a run can be both the last periodic checkpoint and the run bundle: keep the run bundle.
    if times and final is not None and times[-1] == final:
        frames.pop()
        times.pop()
    frames.append(bundle)
    times.append(final)
    labels = [f'{time:g} Myr' if time is not None else item.path.name for item, time in zip(frames, times)]
    if all(time is not None for time in times):
        durations = [max(120, round((b - a) * 120)) for a, b in zip(times, times[1:])] + [1200]
    else:
        durations = [400] * (len(frames) - 1) + [1200]
    width = max(512, bundle.columns)
    written = []
    for spec in specs:
        if not all(available(item, spec) for item in frames):
            continue
        layers = {}
        for item in frames:
            for season, values in model_layers(item, spec):
                buffer = io.BytesIO()
                Image.fromarray(palette_rgb(values, spec['palette'])).save(buffer, format='PNG')
                buffer.seek(0)
                layers.setdefault(season, []).append(buffer)
        for season, buffers in layers.items():
            suffix = '' if season is None else f'-{SEASONS[season].lower()}'
            relative = Path('maps') / f"{spec['name']}{suffix}-animation.png"
            title = spec['title'] + ('' if season is None else f' · {SEASONS[season]}')
            metadata = map_rendering.save_animation(directory / relative, buffers, [f'{title} · {label}' for label in labels],
                                                    durations, width)
            written.append(entry(spec, 'animation', season, f'{title} · animation', relative.as_posix(), spec['palette'],
                                 frames=labels, frame_count=metadata['frames'], duration_ms=metadata['duration_ms']))
    return written


def render(bundle, specs, directory, validation=(), compared=(), references=None, checkpoints=()):
    """Write PNGs under directory/maps; returns the maps entries of evaluation.json (paths relative to directory).

    validation: the profile's validation specs; compared: ids that validation compared
    for this bundle (only those get reference and difference maps, so the maps read no
    reference the score did not use); checkpoints: bundles in time order.
    """
    linked = {item['id']: item for item in validation}
    written = []
    for spec in specs:
        target = linked.get(spec.get('validation'))
        if target is not None and target['id'] in compared and references is not None:
            written += validated(bundle, spec, target, references, directory)
            continue
        if not available(bundle, spec):
            continue
        for season, layer in model_layers(bundle, spec):
            suffix = '' if season is None else f'-{SEASONS[season].lower()}'
            title = spec['title'] + ('' if season is None else f' · {SEASONS[season]}')
            written.append(entry(spec, 'model', season, title,
                                 save(directory, spec['name'] + suffix, palette_rgb(layer, spec['palette'])), spec['palette']))
    if checkpoints:
        written += animations(bundle, specs, directory, checkpoints)
    return written
