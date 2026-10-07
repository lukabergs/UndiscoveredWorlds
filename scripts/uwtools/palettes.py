"""Map palettes: configs/climate-palettes.json and configs/tectonics-palettes.json.

Runs and references are rendered with the same palettes. Each call reads the file afresh.
A maps folder carries palette.json, a sidecar naming the palette files it was rendered
with (version and SHA-256), so a later palette change stays visible beside the images.
"""
import json

from .hashing import sha256
from .paths import ROOT

CLIMATE_PATH = ROOT / 'configs/climate-palettes.json'
TECTONICS_PATH = ROOT / 'configs/tectonics-palettes.json'


def climate():
    return json.loads(CLIMATE_PATH.read_text(encoding='utf-8'))


def tectonics():
    return json.loads(TECTONICS_PATH.read_text(encoding='utf-8'))


def climate_anchors(low, high, signed=False, definition=None):
    """(value, RGB) anchors of a climate scale over low..high."""
    if not high > low:
        raise ValueError('Display range must increase')
    entries = (definition or climate())['signed' if signed else 'positive']
    # Every signed scale gives zero the same cyan colour, including asymmetric T ranges.
    return [(low + 2*t*(0-low) if t <= .5 else 2*(t-.5)*high, colour)
            if signed and low < 0 < high else (low+t*(high-low), colour)
            for t, colour in entries]


def climate_rgb(values, low=0, high=1, signed=False):
    from .map_rendering import colourize
    definition = climate()
    return colourize(values, climate_anchors(low, high, signed, definition), definition['missing'])


def climate_field_rgb(values, field):
    low, high = climate()['scales'][field]
    return climate_rgb(values, low, high, low < 0)


def climate_scale(field):
    """(low, high) of a climate scale."""
    low, high = climate()['scales'][field]
    return low, high


def tectonics_stops(spec, definition=None):
    """(value, RGB) anchors of a continuous tectonics field: its stops, or the sequential ramp over its range."""
    if 'stops' in spec:
        return spec['stops']
    low, high = spec['range']
    return [(low + t * (high - low), rgb) for t, rgb in (definition or tectonics())['sequential']]


def identity_colour(identity):
    """Colour of an identity (plate index or id) on an 'ids' tectonics field."""
    import colorsys
    return [round(c * 255) for c in colorsys.hsv_to_rgb((int(identity) * .61803398875) % 1, .52, .92)]


def tectonics_rgb(values, name):
    """Identities, labelled classes or interpolated stops of tectonics field `name`."""
    import numpy as np
    from .map_rendering import categories, colourize
    definition = tectonics()
    spec = definition['fields'][name]
    if spec.get('ids'):
        colours = {int(identity): identity_colour(identity) for identity in np.unique(values[np.isfinite(values)])}
        return categories(values, colours, definition['missing'])
    if 'colours' in spec:
        return categories(values, dict(enumerate(spec['colours'])), definition['missing'])
    return colourize(values, tectonics_stops(spec, definition), definition['missing'])


def palette_stops(palette):
    """(value, RGB) anchors of a continuous palette 'climate:<scale>' or 'tectonics:<field>', in its own units."""
    family, _, name = palette.partition(':')
    if family == 'climate':
        definition = climate()
        low, high = definition['scales'][name]
        return climate_anchors(low, high, low < 0, definition)
    if family == 'tectonics':
        definition = tectonics()
        spec = definition['fields'][name]
        if spec.get('ids') or 'colours' in spec:
            raise ValueError(f'{palette} is categorical')
        return tectonics_stops(spec, definition)
    raise ValueError(f'Unknown palette {palette!r}; use climate:<scale> or tectonics:<field>')


def palette_rgb(values, palette):
    """RGB of values in palette 'climate:<scale>' or 'tectonics:<field>' (the generated maps' rule)."""
    family, _, name = palette.partition(':')
    if family == 'climate':
        return climate_field_rgb(values, name)
    if family == 'tectonics':
        return tectonics_rgb(values, name)
    raise ValueError(f'Unknown palette {palette!r}; use climate:<scale> or tectonics:<field>')


SIDECAR = 'palette.json'
SIDECAR_SCHEMA = 'uw-palette/v1'


def contract():
    """The palette files in force: {family: {file, version, sha256}}."""
    return {family: {'file': path.relative_to(ROOT).as_posix(), 'version': load()['version'], 'sha256': sha256(path)}
            for family, path, load in (('climate', CLIMATE_PATH, climate), ('tectonics', TECTONICS_PATH, tectonics))}


def write_sidecar(directory):
    """Record the palettes in force in directory/palette.json; called by every renderer of a maps folder."""
    document = {'schema': SIDECAR_SCHEMA, 'palettes': contract()}
    (directory / SIDECAR).write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')


def sidecar_issue(directory):
    """None when directory/palette.json is sound, else what is wrong with it.

    A sidecar of an earlier palette version is sound: the images are evidence of that
    version. One that names the current version with another hash is not: the palette
    file changed without a new version.
    """
    path = directory / SIDECAR
    if not path.is_file():
        return 'no palette.json'
    document = json.loads(path.read_text(encoding='utf-8'))
    if document.get('schema') != SIDECAR_SCHEMA:
        return f"palette.json schema is {document.get('schema')!r}"
    current = contract()
    for family, recorded in document.get('palettes', {}).items():
        if family in current and recorded['version'] == current[family]['version'] and recorded['sha256'] != current[family]['sha256']:
            return f"palette.json names {recorded['version']} with another hash than {current[family]['file']}"
    return None
