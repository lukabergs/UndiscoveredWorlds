"""Map palettes: configs/climate-palettes.json and configs/tectonics-palettes.json.

Runs and references are rendered with the same palettes. Each call reads the file afresh.
"""
import json

from .paths import ROOT

CLIMATE_PATH = ROOT / 'configs/climate-palettes.json'
TECTONICS_PATH = ROOT / 'configs/tectonics-palettes.json'


def climate():
    return json.loads(CLIMATE_PATH.read_text(encoding='utf-8'))


def tectonics():
    return json.loads(TECTONICS_PATH.read_text(encoding='utf-8'))
