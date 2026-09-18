"""One palette contract for numerical exports, reference maps and comparisons."""
import json
from pathlib import Path
import numpy as np

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
    values = np.asarray(values)
    stops = anchors(low, high, signed)
    safe = np.where(np.isfinite(values), values, 0)
    result = np.stack([np.interp(safe, [v for v, _ in stops], [c[k] for _, c in stops])
                       for k in range(3)], axis=-1)
    result[~np.isfinite(values)] = DEFINITION['missing']
    return np.floor(result+.5).astype(np.uint8)

def field_rgb(values, field):
    low, high = DEFINITION['scales'][field]
    return rgb(values, low, high, low < 0)
