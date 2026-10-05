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
