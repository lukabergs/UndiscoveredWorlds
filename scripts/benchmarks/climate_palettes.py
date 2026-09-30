"""One palette contract for numerical exports, reference maps and comparisons."""
import json
from pathlib import Path
from map_rendering import colourize

ROOT = Path(__file__).resolve().parents[2]
DEFINITION = json.loads((ROOT / 'configs/climate-palettes.json').read_text())
VERSION = DEFINITION['version']

def anchors(low, high, signed=False):
    if not high > low:
        raise ValueError('Display range must increase')
    entries = DEFINITION['signed' if signed else 'positive']
    # Every signed scale gives zero the same cyan colour, including asymmetric T ranges.
    return [(low + 2*t*(0-low) if t <= .5 else 2*(t-.5)*high, colour)
            if signed and low < 0 < high else (low+t*(high-low), colour)
            for t, colour in entries]

def rgb(values, low=0, high=1, signed=False):
    return colourize(values, anchors(low, high, signed), DEFINITION['missing'])

def field_rgb(values, field):
    low, high = DEFINITION['scales'][field]
    return rgb(values, low, high, low < 0)
