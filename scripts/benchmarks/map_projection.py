"""Reproject cell-centred equirectangular map images for viewing.

Every simulator map is a global W x W/2 grid: longitude -180..180 left to right,
latitude 90..-90 top to bottom, cell centred. These functions only resample the
coloured image (nearest cell), so no value is interpolated or invented.
"""
from pathlib import Path

import numpy as np
from PIL import Image

PROJECTIONS = {
    'equirectangular': 'Equirectangular (native grid)',
    'mollweide': 'Mollweide (equal area)',
    'north-polar': 'North polar (azimuthal equal area, 90°N–0°)',
    'south-polar': 'South polar (azimuthal equal area, 90°S–0°)',
}


def _sample(rgb, longitude, latitude, valid):
    height, width = rgb.shape[:2]
    column = np.floor((longitude + 180.0) / 360.0 * width).astype(int) % width
    row = np.clip(np.floor((90.0 - latitude) / 180.0 * height).astype(int), 0, height - 1)
    result = np.zeros((*longitude.shape, 4), dtype=np.uint8)
    result[..., :3] = rgb[row, column, :3]
    result[..., 3] = np.where(valid, 255, 0)
    result[~valid, :3] = 0
    return result


def _mollweide(rgb, width):
    height = width // 2
    y, x = np.indices((height, width), dtype=float)
    # Ellipse with semi-axes 2*sqrt(2) by sqrt(2).
    px = ((x + 0.5) / width * 2 - 1) * 2 * np.sqrt(2)
    py = (1 - (y + 0.5) / height * 2) * np.sqrt(2)
    valid = (px / (2 * np.sqrt(2))) ** 2 + (py / np.sqrt(2)) ** 2 <= 1
    theta = np.arcsin(np.clip(py / np.sqrt(2), -1, 1))
    latitude = np.degrees(np.arcsin(np.clip((2 * theta + np.sin(2 * theta)) / np.pi, -1, 1)))
    cosine = np.cos(theta)
    longitude = np.degrees(np.pi * px / (2 * np.sqrt(2) * np.where(np.abs(cosine) < 1e-12, 1e-12, cosine)))
    valid &= np.abs(longitude) <= 180
    return _sample(rgb, np.clip(longitude, -180, 180), latitude, valid)


def _polar(rgb, size, north):
    y, x = np.indices((size, size), dtype=float)
    # Unit disc reaches the equator: Lambert azimuthal radius r = 2 sin(c/2), c = 90° -> sqrt(2).
    px = ((x + 0.5) / size * 2 - 1) * np.sqrt(2)
    py = (1 - (y + 0.5) / size * 2) * np.sqrt(2)
    radius = np.hypot(px, py)
    valid = radius <= np.sqrt(2)
    colatitude = 2 * np.degrees(np.arcsin(np.clip(radius / 2, 0, 1)))
    if north:
        latitude = 90 - colatitude
        longitude = np.degrees(np.arctan2(px, -py))
    else:
        latitude = colatitude - 90
        longitude = np.degrees(np.arctan2(px, py))
    return _sample(rgb, longitude, latitude, valid)


def project(rgb, name, width=None):
    """Return an RGBA view of an equirectangular RGB(A) image; transparent outside the globe."""
    rgb = np.asarray(rgb)
    if rgb.ndim != 3 or rgb.shape[1] != 2 * rgb.shape[0]:
        raise ValueError('Expected a global W x W/2 equirectangular image')
    width = width or rgb.shape[1]
    if name == 'equirectangular':
        result = np.zeros((rgb.shape[0], rgb.shape[1], 4), dtype=np.uint8)
        result[..., :3] = rgb[..., :3]
        result[..., 3] = 255
        return result
    if name == 'mollweide':
        return _mollweide(rgb, width)
    if name in ('north-polar', 'south-polar'):
        return _polar(rgb, width // 2, name == 'north-polar')
    raise ValueError(f'Unknown projection {name}; choose from {", ".join(PROJECTIONS)}')


def projected_path(map_path, maps_root, name):
    """runs/maps/<module>/<rest> -> runs/maps/<module>/projections/<name>/<rest>."""
    relative = Path(map_path).relative_to(maps_root)
    return Path(maps_root) / relative.parts[0] / 'projections' / name / Path(*relative.parts[1:])


def save_projections(map_path, maps_root, names=tuple(PROJECTIONS)):
    """Write projected copies of one exported map; reuse them while newer than the source."""
    map_path = Path(map_path)
    with Image.open(map_path) as image:
        rgb = np.asarray(image.convert('RGB'))
    written = {}
    for name in names:
        if name == 'equirectangular':
            written[name] = map_path
            continue
        target = projected_path(map_path, maps_root, name)
        if not target.exists() or target.stat().st_mtime < map_path.stat().st_mtime:
            target.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(project(rgb, name)).save(target)
        written[name] = target
    return written
