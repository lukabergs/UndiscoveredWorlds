"""Review maps of bundle fields in the project palettes (climate and tectonics).

Every image is framed for reading (uwtools.map_rendering.frame): the raster enlarged
without smoothing, a graticule every 30 degrees, the title and the legend of its palette.
A profile map with "validation": "<id>" is drawn from that validation entry when it was
compared: model and reference reduced by validation.py's season rule, conversion and
remap, and model minus reference on the score's mask, so the pictures show what the loss
measures. Vector fields are coloured by speed with arrows for direction; their difference
is the speed difference. Differences use the diverging scale of
configs/climate-palettes.json ("difference"). With a parent bundle the model side of
every map is drawn for the parent too (kind "parent"). Runs with checkpoint bundles get
one APNG per map and season; plate colours in them follow identities tracked from frame
to frame (track_identities), because the engine's plate index is not persistent.
"""
import colorsys
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import map_rendering, palettes
from uwtools.bundle import SEASONS
from checks import whole_metres
from validation import field_mask, load_selected, model_field, reference_field, seasons_for, surface_elevation

CLIMATE = palettes.climate()
TECTONICS = palettes.tectonics()


def lake_elevation(bundle):
    """Lake surfaces relative to sea level (m); NaN where there is no lake."""
    surface = bundle.load('hydrosphere.lake_surface').astype(float)
    return np.where(surface > 0, surface - float(bundle.scalar('geosphere.sea_level')), np.nan)


def lake_depth(bundle):
    """Depth of the routing's filled surface above the terrain (m); zero outside depressions."""
    return bundle.load('hydrosphere.filled_elevation').astype(float) - bundle.load('geosphere.base_elevation').astype(float)


def depression_class(bundle):
    """Depression id folded onto 12 colour classes (1..12); 0 outside depressions."""
    identity = bundle.load('hydrosphere.depression_id').astype(np.int64)
    return np.where(identity != 0, identity % 12 + 1, 0).astype(float)


# Flow direction codes as river_direction: (dx, dy) with east = +x, south = +y (north row first).
FLOW_OFFSETS = {1: (0, -1), 2: (1, -1), 3: (1, 0), 4: (1, 1), 5: (0, 1), 6: (-1, 1), 7: (-1, 0), 8: (-1, -1)}


def drainage_basin(bundle):
    """The outlet of every land cell along flow_direction, folded onto 12 colour classes (1..12); 0 on sea."""
    direction = bundle.load('hydrosphere.flow_direction')
    land = whole_metres(bundle.load('geosphere.base_elevation').astype(float)) > float(bundle.scalar('geosphere.sea_level'))
    rows, columns = direction.shape
    y, x = np.indices((rows, columns))
    receiver = (y * columns + x).ravel()
    for code, (dx, dy) in FLOW_OFFSETS.items():
        ny = y + dy
        cells = (direction == code) & (ny >= 0) & (ny < rows)
        receiver[cells.ravel()] = (ny * columns + (x + dx) % columns)[cells]
    for _ in range(64):
        following = receiver[receiver]
        if np.array_equal(following, receiver):
            break
        receiver = following
    return np.where(land, receiver.reshape(rows, columns) % 12 + 1, 0).astype(float)


# Derived maps: the function and the fields it reads.
DERIVED = {'surface_elevation': surface_elevation, 'lake_elevation': lake_elevation, 'lake_depth': lake_depth,
           'depression_class': depression_class, 'drainage_basin': drainage_basin}
DERIVED_FIELDS = {'surface_elevation': ['geosphere.base_elevation', 'geosphere.mountain_height', 'geosphere.volcano'],
                  'lake_elevation': ['hydrosphere.lake_surface', 'geosphere.sea_level'],
                  'lake_depth': ['hydrosphere.filled_elevation', 'geosphere.base_elevation'],
                  'depression_class': ['hydrosphere.depression_id'],
                  'drainage_basin': ['hydrosphere.flow_direction', 'geosphere.base_elevation', 'geosphere.sea_level']}
QUARTERS = ('Jan–Mar', 'Apr–Jun', 'Jul–Sep', 'Oct–Dec')
SNAPSHOT_TABLE = 'geosphere.tectonic_snapshot'
TRACKED = 'colours follow plate identities tracked by cell overlap between consecutive frames, not the engine index'


def tectonics_stops(spec):
    if 'stops' in spec:
        return spec['stops']
    low, high = spec['range']
    return [(low + t * (high - low), rgb) for t, rgb in TECTONICS['sequential']]


def tectonics_rgb(values, name):
    """The rule of tectonics_report.py: identities, labelled classes or interpolated stops."""
    spec = TECTONICS['fields'][name]
    if spec.get('ids'):
        colours = {int(identity): [round(c * 255) for c in colorsys.hsv_to_rgb((int(identity) * .61803398875) % 1, .52, .92)]
                   for identity in np.unique(values[np.isfinite(values)])}
        return map_rendering.categories(values, colours, TECTONICS['missing'])
    if 'colours' in spec:
        return map_rendering.categories(values, dict(enumerate(spec['colours'])), TECTONICS['missing'])
    return map_rendering.colourize(values, tectonics_stops(spec), TECTONICS['missing'])


def palette_rgb(values, palette):
    family, _, name = palette.partition(':')
    if family == 'climate':
        return palettes.climate_field_rgb(values, name)
    if family == 'climate-difference':
        return map_rendering.colourize(values, difference_stops(name), CLIMATE['missing'])
    if family == 'tectonics':
        return tectonics_rgb(values, name)
    raise ValueError(f'Unknown palette {palette!r}; use climate:<scale> or tectonics:<field>')


def difference_stops(scale):
    """Diverging colours symmetric about zero over +-half_range of the scale."""
    spec = CLIMATE['difference']
    half = float(spec['half_ranges'][scale])
    return [(-half + 2 * half * t, colour) for t, colour in spec['stops']]


def legend(palette, note=None):
    """The legend of a palette for map_rendering.frame (the title of a map carries its units)."""
    family, _, name = palette.partition(':')
    if family == 'climate':
        low, high = CLIMATE['scales'][name]
        return {'kind': 'scale', 'stops': palettes.climate_anchors(low, high, low < 0, CLIMATE), 'note': note}
    if family == 'climate-difference':
        return {'kind': 'scale', 'stops': difference_stops(name), 'units': 'model − reference', 'note': note}
    spec = TECTONICS['fields'][name]
    if spec.get('ids'):
        return {'kind': 'classes', 'entries': [], 'note': note or 'one colour per plate (engine index)'}
    if 'colours' in spec:
        return {'kind': 'classes', 'entries': list(zip(spec['labels'], spec['colours'])), 'note': note}
    stops = tectonics_stops(spec)
    ticks = [value for value, _ in stops]
    ticks = ticks if 'stops' in spec else None
    return {'kind': 'scale', 'stops': stops, 'ticks': ticks, 'note': note}


def available(bundle, spec):
    selectors = spec.get('vector', [spec['field']] if 'field' in spec else [])
    if 'derived' in spec:
        selectors = DERIVED_FIELDS[spec['derived']]
    return all(selector.split('[')[0] in bundle.fields for selector in selectors)


def masked(bundle, spec, values):
    if spec.get('mask') in ('sea', 'land'):
        sea = bundle.sea_fraction(*values.shape[-2:]) > 0.9
        values = np.where(sea if spec['mask'] == 'sea' else ~sea, values, np.nan)
    return values


def transformed(spec, values):
    """Map spec options: "hide_zero" (zeros become NaN) and "transform": "log10" (values <= 0 become NaN)."""
    if spec.get('hide_zero'):
        values = np.where(values == 0, np.nan, values)
    if spec.get('transform') == 'log10':
        with np.errstate(divide='ignore', invalid='ignore'):
            values = np.where(values > 0, np.log10(np.where(values > 0, values, 1.0)), np.nan)
    elif spec.get('transform') is not None:
        raise ValueError(f"Unknown map transform {spec['transform']!r}; known: log10")
    return values


def model_layers(bundle, spec):
    """[(season or None, values, arrows or None)] of a profile map drawn from its own field."""
    components = None
    if 'derived' in spec:
        values = DERIVED[spec['derived']](bundle)
    elif 'vector' in spec:
        components = [load_selected(bundle, selector).astype(float) for selector in spec['vector']]
        if spec.get('reduce') == 'season-mean':
            components = [component.mean(axis=0) for component in components]
        values = np.hypot(*components)
    else:
        values = load_selected(bundle, spec['field'])
    if spec.get('reduce') == 'season-mean' and components is None:
        values = values.mean(axis=0)
    values = transformed(spec, masked(bundle, spec, values.astype(float) * spec.get('scale', 1)))
    if values.ndim == 2:
        return [(None, values, arrows(components, values))]
    return [(season, values[season], arrows(components and [component[season] for component in components], values[season]))
            for season in spec.get('seasons', range(values.shape[0]))]


def arrows(components, shown):
    """(east, north) where the map shows a value, or None for a scalar map."""
    if components is None:
        return None
    return tuple(np.where(np.isfinite(shown), component, np.nan) for component in components)


def speed(values):
    return np.hypot(*values) if values.ndim == 3 else values


def season_title(rule, season):
    if season is None or rule == 'annual':
        return None
    return QUARTERS[season] if rule == 'quarter' else SEASONS[season]


def validated_seasons(spec, validation):
    seasons = seasons_for(validation)
    if 'seasons' in spec and len(seasons) > 1:
        seasons = [season for season in seasons if season in spec['seasons']]
    annual = validation['season'] == 'annual' or validation.get('reduce') == 'season-mean'
    return [(season, None if annual else season) for season in seasons]


def model_images(bundle, spec, validation):
    """(season or None, the label of the months shown, values, arrows) of the model side of a map:
    by the validation entry's reduction when the map is linked to one, else from its own field."""
    if validation is None:
        for season, values, vectors in model_layers(bundle, spec):
            yield season, None if season is None else SEASONS[season], values, vectors
        return
    for season, layer in validated_seasons(spec, validation):
        model = model_field(bundle, validation, season)
        values = masked(bundle, spec, speed(model))
        yield layer, season_title(validation['season'], layer), values, arrows(model if model.ndim == 3 else None, values)


def save(directory, name, rgb, title, palette, vectors=None, note=None):
    relative = Path('maps') / f'{name}.png'
    (directory / relative).parent.mkdir(parents=True, exist_ok=True)
    map_rendering.frame(rgb, title, legend(palette, note), vectors).save(directory / relative)
    return relative.as_posix()


def entry(spec, kind, season, title, path, palette, **extra):
    """season: the model season index (None when annual); season_label (in extra) names the months shown."""
    name = None if season is None else SEASONS[season]
    return {'map': spec['name'], 'kind': kind, 'season': name, 'season_label': name, 'title': title, 'path': path,
            'palette': palette, **extra}


def suffix(season):
    return '' if season is None else f'-{SEASONS[season].lower()}'


def validated(bundle, spec, validation, references, directory):
    """Model, reference and difference images of a map linked to a compared validation entry."""
    scale = validation.get('metric_scale', validation['id'])
    written = []
    for season, layer in validated_seasons(spec, validation):
        model = model_field(bundle, validation, season)
        reference = reference_field(validation, references, season, model.shape[-1])
        mask = np.isfinite(speed(reference)) & np.isfinite(speed(model)) & field_mask(bundle, validation, season, speed(model).shape)
        label = season_title(validation['season'], layer)
        title = spec['title'] + ('' if label is None else f' · {label}')
        source = ', '.join(validation['reference'])
        common = {'validation': validation['id'], 'season_label': label}
        shown = masked(bundle, spec, speed(model))
        written.append(entry(spec, 'model', layer, title,
                             save(directory, spec['name'] + suffix(layer), palette_rgb(shown, spec['palette']), title, spec['palette'],
                                  arrows(model if model.ndim == 3 else None, shown)),
                             spec['palette'], **common))
        written.append(entry(spec, 'reference', layer, f'{title} · reference ({source})',
                             save(directory, f"{spec['name']}{suffix(layer)}-reference", palette_rgb(speed(reference), spec['palette']),
                                  f'{title} · reference ({source})', spec['palette'],
                                  arrows(reference if reference.ndim == 3 else None, speed(reference))),
                             spec['palette'], **common))
        difference = np.where(mask, speed(model) - speed(reference), np.nan)
        palette = f'climate-difference:{scale}'
        written.append(entry(spec, 'difference', layer, f'{title} · model − reference',
                             save(directory, f"{spec['name']}{suffix(layer)}-difference", palette_rgb(difference, palette),
                                  f'{title} · model − reference', palette),
                             palette, **common, half_range=CLIMATE['difference']['half_ranges'][scale]))
    return written


def checkpoint_time(bundle):
    if SNAPSHOT_TABLE in bundle.tables:
        return json.loads(bundle.table_text(SNAPSHOT_TABLE)).get('time_myr')
    return None


def is_identity(spec):
    family, _, name = spec['palette'].partition(':')
    return family == 'tectonics' and bool(TECTONICS['fields'][name].get('ids'))


def track_identities(frames):
    """Plate-index rasters in time order -> (rasters of tracked identities, events per frame).

    A plate keeps the identity of the plate of the previous frame it shares most cells
    with; when two plates claim the same one, the larger overlap keeps it and the other
    gets a new identity. Events name the identities that appear and vanish at each frame,
    so a split, a merger and a plate that disappears are visible as such.
    """
    tracked, events, previous, following = [], [], None, 0
    for values in frames:
        values = np.asarray(values, dtype=float)
        indices = [int(index) for index in np.unique(values[np.isfinite(values)])]
        appeared = []
        if previous is None:
            mapping = {index: index for index in indices}
            following = max(indices, default=-1) + 1
        else:
            claims = []
            for index in indices:
                identities, counts = np.unique(previous[values == index], return_counts=True)
                known = identities >= 0
                best = counts[known].argmax() if known.any() else None
                claims.append((int(counts[known][best]) if best is not None else 0, index,
                               int(identities[known][best]) if best is not None else -1))
            mapping, taken = {}, set()
            for _, index, identity in sorted(claims, reverse=True):
                if identity < 0 or identity in taken:
                    identity, following = following, following + 1
                    appeared.append(identity)
                taken.add(identity)
                mapping[index] = identity
        current = np.full(values.shape, -1, dtype=int)
        for index, identity in mapping.items():
            current[values == index] = identity
        before = set() if previous is None else {int(identity) for identity in np.unique(previous) if identity >= 0}
        events.append({'appeared': sorted(appeared), 'vanished': sorted(before - set(mapping.values()))})
        tracked.append(np.where(current >= 0, current, np.nan))
        previous = current
    return tracked, events


def timeline(bundle, checkpoints):
    """The checkpoint bundles and the run bundle in time order, with frame labels and durations."""
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
    return frames, labels, durations


def animations(bundle, specs, directory, checkpoints):
    """One APNG per map and season over the checkpoint bundles and the run bundle.

    Returns the entries and, per identity map, the tracked raster of the run bundle,
    so its static map has the colours of the last frame.
    """
    frames, labels, durations = timeline(bundle, checkpoints)
    width = max(512, bundle.columns)
    written, final = [], {}
    for spec in specs:
        if not all(available(item, spec) for item in frames):
            continue
        layers = {}
        for item in frames:
            for season, values, vectors in model_layers(item, spec):
                layers.setdefault(season, []).append((values, vectors))
        for season, series in layers.items():
            extra, note = {}, None
            if is_identity(spec):
                tracked, events = track_identities([values for values, _ in series])
                series = [(values, None) for values in tracked]
                final[spec['name']] = tracked[-1]
                note = TRACKED
                extra = {'identity': note, 'identity_events': [dict(event, frame=label) for event, label in zip(events, labels)
                                                              if event['appeared'] or event['vanished']]}
            relative = Path('maps') / f"{spec['name']}{suffix(season)}-animation.png"
            title = spec['title'] + ('' if season is None else f' · {SEASONS[season]}')
            images = [map_rendering.frame(palette_rgb(values, spec['palette']), f'{title} · {label}', legend(spec['palette'], note),
                                          vectors, width)
                      for (values, vectors), label in zip(series, labels)]
            metadata = map_rendering.save_frames(directory / relative, images, durations)
            written.append(entry(spec, 'animation', season, f'{title} · animation', relative.as_posix(), spec['palette'],
                                 frames=labels, frame_count=metadata['frames'], duration_ms=metadata['duration_ms'], **extra))
    return written, final


def render(bundle, specs, directory, validation=(), compared=(), references=None, checkpoints=(), parent=None):
    """Write PNGs under directory/maps; returns the maps entries of evaluation.json (paths relative to directory).

    validation: the profile's validation specs; compared: ids that validation compared
    for this bundle (only those get reference and difference maps, so the maps read no
    reference the score did not use); checkpoints: bundles in time order; parent: the
    bundle the run is compared with, drawn as kind "parent" with the reduction of the run's maps.
    """
    linked = {item['id']: item for item in validation}
    animated, tracked = animations(bundle, specs, directory, checkpoints) if checkpoints else ([], {})
    written = []
    for spec in specs:
        target = linked.get(spec.get('validation'))
        if target is not None and not (target['id'] in compared and references is not None):
            target = None
        if target is not None:
            written += validated(bundle, spec, target, references, directory)
        elif available(bundle, spec):
            for season, label, values, vectors in model_images(bundle, spec, None):
                note = None
                if spec['name'] in tracked and season is None:
                    values, note = tracked[spec['name']], TRACKED
                title = spec['title'] + ('' if label is None else f' · {label}')
                written.append(entry(spec, 'model', season, title,
                                     save(directory, spec['name'] + suffix(season), palette_rgb(values, spec['palette']), title,
                                          spec['palette'], vectors, note), spec['palette']))
        else:
            continue
        if parent is not None and available(parent, spec):
            for season, label, values, vectors in model_images(parent, spec, target):
                title = spec['title'] + ('' if label is None else f' · {label}') + ' · parent run'
                written.append(entry(spec, 'parent', season, title,
                                     save(directory, f"parent-{spec['name']}{suffix(season)}", palette_rgb(values, spec['palette']),
                                          title, spec['palette'], vectors), spec['palette'], season_label=label))
    written += animated
    if written:
        palettes.write_sidecar(directory / 'maps')
    return written
