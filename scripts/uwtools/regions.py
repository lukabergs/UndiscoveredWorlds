"""Regions of configs/regions.json and their masks on cell-centred W x W/2 grids.

Two sets: the project's latitude/longitude boxes ('project') and the IPCC AR6 WGI
reference regions ('ar6'). A box selects the cells whose centre has
south < latitude <= north and west <= longitude < east (west > east crosses the date
line), the cells the climate report's map crop shows. AR6 regions come from the
rasterised mask (pixel value = id + 1, 0 outside every region) reduced to the target
grid by area-mean overlap fractions, as the validation summary reduces them.
"""
import json

import numpy as np

from .grid import latitudes, longitudes, remap_area_mean
from .paths import DATA_ROOT, ROOT

PATH = ROOT / 'configs/regions.json'
SCHEMA = 'uw-regions/v1'


def load(path=PATH):
    config = json.loads(path.read_text(encoding='utf-8'))
    if config.get('schema') != SCHEMA:
        raise ValueError(f'{path}: expected schema {SCHEMA}')
    return config


def boxes(config=None):
    """Project boxes: dicts with id, name, south, north, west, east (degrees)."""
    return (config or load())['sets']['project']['regions']


def box_mask(box, rows, columns=None):
    """Boolean rows x columns mask of the cells whose centre lies in the box."""
    columns = 2 * rows if columns is None else columns
    lat, lon = latitudes(rows), longitudes(columns)
    in_lat = (lat > box['south']) & (lat <= box['north'])
    west, east = box['west'], box['east']
    in_lon = (lon >= west) & (lon < east) if west < east else (lon >= west) | (lon < east)
    return in_lat[:, None] & in_lon[None, :]


def box_masks(rows, columns=None, config=None):
    """{box id: boolean mask} in config order."""
    return {box['id']: box_mask(box, rows, columns) for box in boxes(config)}


def ar6(config=None):
    """AR6 regions: dicts with id (int), acronym, name, continent, type."""
    return (config or load())['sets']['ar6']['regions']


def ar6_mask_path(config=None):
    return DATA_ROOT / (config or load())['sets']['ar6']['mask']


def ar6_labels(config=None):
    """Label grid of the rasterised AR6 mask: region id per cell, -1 outside every region."""
    from PIL import Image
    path = ar6_mask_path(config)
    if not path.is_file():
        raise FileNotFoundError(f'Missing AR6 region mask {path.relative_to(ROOT)}; see refs/README.md')
    with Image.open(path) as image:
        return np.asarray(image.convert('RGB'))[:, :, 0].astype(int) - 1


def ar6_fractions(rows, columns=None, config=None, labels=None):
    """{region id: fraction of each cell of a rows x columns grid inside the region}, clipped to [0, 1]."""
    columns = 2 * rows if columns is None else columns
    config = config or load()
    labels = ar6_labels(config) if labels is None else labels
    return {int(item['id']): np.clip(remap_area_mean((labels == int(item['id'])).astype(float), rows, columns), 0, 1)
            for item in ar6(config)}
