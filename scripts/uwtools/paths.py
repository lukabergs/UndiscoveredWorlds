"""Build output locations: the one place scripts take executable paths from.

Presets (CMakePresets.json) configure into out/build/<preset>; binaries land in
bin/<configuration>. Scripts run the headless preset's Release binaries.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRESET = 'headless'
CONFIGURATION = 'Release'


def build_directory(preset=PRESET):
    return ROOT / 'out/build' / preset


def binary_directory(preset=PRESET, configuration=CONFIGURATION):
    return build_directory(preset) / 'bin' / configuration


def executable(name, preset=PRESET, configuration=CONFIGURATION):
    return binary_directory(preset, configuration) / f'{name}.exe'
