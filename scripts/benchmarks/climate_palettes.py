"""The climate palette functions are in scripts/uwtools/palettes.py. These names stay for the scripts that import them from here."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import palettes
from uwtools.palettes import climate_anchors as anchors, climate_field_rgb as field_rgb, climate_rgb as rgb

DEFINITION = palettes.climate()
VERSION = DEFINITION['version']
