"""Numeric reference rasters and deterministic circulation previews.

Palettes and LIC settings follow climate_circulation_diagnostics.cpp. This is a
headless NumPy renderer, not a pixel-identical SFML renderer.
"""

from pathlib import Path
import sys

import numpy as np
from PIL import Image, ImageDraw, TiffImagePlugin
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from climate_palettes import rgb as shared_rgb

from cell_grid import latitude_centres, latitude_edges

EARTH_RADIUS_M = 6371000.0
NODATA = -9999.9


def write_float32_geotiff(path: Path, values: np.ndarray) -> None:
    """Single-band, unclamped IEEE float32; same grid/nodata as benchmark TIFFs."""
    height, width = values.shape
    tags = TiffImagePlugin.ImageFileDirectory_v2()
    tags[277] = 1
    tags[284] = 1
    tags[33550] = (360.0 / width, 180.0 / height, 0.0)
    tags[33922] = (0.0, 0.0, 0.0, -180.0, 90.0, 0.0)
    tags[34735] = (1, 1, 0, 3, 1024, 0, 1, 2, 1025, 0, 1, 1, 2048, 0, 1, 4326)
    tags[42113] = str(NODATA)
    path.parent.mkdir(parents=True, exist_ok=True)
    stored = np.where(np.isfinite(values), values, NODATA).astype(np.float32)
    Image.fromarray(stored).save(path, compression="raw", tiffinfo=tags)


def scalar_rgb(values, limit, palette="speed"):
    return shared_rgb(values, 0 if palette == 'speed' else -limit, limit, palette != 'speed')


def divergence(east, north):
    """Spherical finite-volume divergence, per second; east/north components.

    Centre values are averaged to faces. Polar faces have zero length and
    longitude is periodic, so the global area integral telescopes to zero.
    """
    height, width = east.shape
    lat = np.deg2rad(latitude_edges(height))
    dlambda = 2.0 * np.pi / width
    east_face = (east + np.roll(east, -1, axis=1)) * 0.5
    north_face = np.zeros((height + 1, width), dtype=np.float64)
    north_face[1:-1] = (north[:-1] + north[1:]) * 0.5
    north_face *= np.cos(lat)[:, None]
    net = ((east_face - np.roll(east_face, 1, axis=1)) * (np.pi / height)
           + (north_face[:-1] - north_face[1:]) * dlambda)
    area = dlambda * (np.sin(lat[:-1]) - np.sin(lat[1:]))[:, None]
    return (net / (EARTH_RADIUS_M * area)).astype(np.float32)


def sample(field, x, y):
    height, width = field.shape
    x = x % width
    y = np.clip(y, 0, height - 1)
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    x1, y1 = (x0 + 1) % width, np.minimum(y0 + 1, height - 1)
    fx, fy = x - x0, y - y0
    # Zero-weight neighbours must not contaminate exact samples with NaN.
    result = np.zeros_like(x, dtype=np.float64)
    for xx, yy, weight in ((x0, y0, (1-fx)*(1-fy)), (x1, y0, fx*(1-fy)),
                           (x0, y1, (1-fx)*fy), (x1, y1, fx*fy)):
        result += np.where(weight > 0, field[yy, xx], 0.0) * weight
    return result


def grid_velocity(east, north, x, y):
    height, width = east.shape
    latitude = np.deg2rad(90.0 - (y + 0.5) * 180.0 / height)
    dx = 2 * np.pi * EARTH_RADIUS_M / width * np.maximum(np.cos(latitude), 1e-6)
    dy = np.pi * EARTH_RADIUS_M / height
    return sample(east, x, y) / dx, -sample(north, x, y) / dy


def direction(east, north, x, y):
    dx, dy = grid_velocity(east, north, x, y)
    norm = np.hypot(dx, dy)
    valid = np.isfinite(norm) & (norm > 1e-12)
    return (np.divide(dx, norm, out=np.zeros_like(dx), where=valid),
            np.divide(dy, norm, out=np.zeros_like(dy), where=valid), valid)


def lic_luminance(east, north, seed):
    """12 steps each way, 0.65-cell RK2, cosine weights, fixed MT19937 seed."""
    height, width = east.shape
    noise = np.random.Generator(np.random.MT19937(seed)).random(east.shape)
    yy, xx = np.indices(east.shape, dtype=np.float64)
    total, weights = noise.copy(), np.ones(east.shape)
    for sign in (-1, 1):
        x, y = xx.copy(), yy.copy()
        active = np.isfinite(east) & np.isfinite(north)
        for step in range(1, 13):
            dx, dy, valid = direction(east, north, x, y)
            mx, my = (x + sign * dx * 0.325) % width, y + sign * dy * 0.325
            active &= valid & (my >= 0) & (my <= height - 1)
            dx, dy, valid = direction(east, north, mx, my)
            x, y = (x + sign * dx * 0.65) % width, y + sign * dy * 0.65
            active &= valid & (y >= 0) & (y <= height - 1)
            weight = (0.5 + 0.5 * np.cos(np.pi * step / 13)) * active
            total += weight * sample(noise, x, y)
            weights += weight
    result = np.clip((total / weights - 0.5) * 2.2 + 0.5, 0, 1).astype(np.float32)
    result[~np.isfinite(east) | ~np.isfinite(north)] = np.nan
    return result


def flow_background(east, north, land, limit=25):
    rgb = scalar_rgb(np.hypot(east, north), limit).astype(float)
    tint = np.where(land[..., None], (32, 27, 23), (4, 13, 25))
    rgb = rgb * 0.45 + tint * 0.55
    coast = (land != np.roll(land, 1, axis=1)) | (land != np.roll(land, -1, axis=1))
    coast[1:] |= land[1:] != land[:-1]
    coast[:-1] |= land[:-1] != land[1:]
    rgb[coast] = rgb[coast] * 0.55 + np.array((155, 155, 145)) * 0.45
    return rgb


def vector_preview(rgb, east, north):
    image = Image.fromarray(np.clip(rgb, 0, 255).astype(np.uint8))
    draw = ImageDraw.Draw(image)
    height, width = east.shape
    spacing = max(1, width // 16)
    length = min(7.0, spacing * 0.8)
    for y in range(spacing // 2, height, spacing):
        for x in range(spacing // 2, width, spacing):
            dx, dy, valid = direction(east, north, np.array(float(x)), np.array(float(y)))
            if not valid:
                continue
            dx, dy = float(dx), float(dy)
            tip = (x + dx * length / 2, y + dy * length / 2)
            draw.line((x - dx * length / 2, y - dy * length / 2, *tip), fill=(245, 225, 95))
            for sign in (-1, 1):
                draw.line((*tip, tip[0] - dx * length * .35 + sign * dy * length * .2,
                           tip[1] - dy * length * .35 - sign * dx * length * .2), fill=(245, 225, 95))
    return image


def lic_preview(east, north, land, luminance, limit=25):
    rgb = flow_background(east, north, land, limit)
    light = np.nan_to_num(luminance)[..., None]
    rgb = np.where(light < .5, rgb * (1 - (.5-light)*1.4),
                   rgb + (255-rgb)*(light-.5)*1.1)
    rgb[~np.isfinite(luminance)] = 0
    return vector_preview(rgb, east, north)


def particle_intensity(east, north, seed):
    """Area-uniform seeds, 36 three-hour intervals; adaptive half-cell RK2."""
    height, width = east.shape
    rng = np.random.Generator(np.random.MT19937(seed))
    count = min(20000, max(5000, width * height // 10))
    x = rng.uniform(0, width, count)
    y = np.clip((90 - np.rad2deg(np.arcsin(rng.uniform(-.995, .995, count))))
                * height / 180 - .5, 0, height - 1)
    remaining = np.full(count, 36 * 10800.0)
    hits = np.zeros((height, width), dtype=np.float64)
    while x.size:
        dx, dy = grid_velocity(east, north, x, y)
        rate = np.hypot(dx, dy)
        valid = np.isfinite(rate) & (rate > 1e-12) & (remaining > 0)
        x, y, dx, dy, rate, remaining = [a[valid] for a in (x, y, dx, dy, rate, remaining)]
        if not x.size:
            break
        dt = np.minimum(remaining, np.minimum(10800.0, 0.5 / rate))
        mx, my = (x + dx * dt / 2) % width, y + dy * dt / 2
        dx, dy = grid_velocity(east, north, mx, my)
        nx, ny = (x + dx * dt) % width, y + dy * dt
        valid = (np.isfinite(nx) & np.isfinite(ny) & (my >= 0) & (my <= height-1)
                 & (ny >= 0) & (ny <= height-1))
        x, y, remaining = nx[valid], ny[valid], (remaining-dt)[valid]
        np.add.at(hits, (np.rint(y).astype(int), np.rint(x).astype(int) % width), 1)
    intensity = (1 - np.power(.84, hits)).astype(np.float32)
    intensity[~np.isfinite(east) | ~np.isfinite(north)] = np.nan
    return intensity


def particle_preview(east, north, land, intensity, limit=25):
    rgb = flow_background(east, north, land, limit)
    alpha = np.nan_to_num(intensity)[..., None]
    rgb += (np.array((245, 247, 245)) - rgb) * alpha
    rgb[~np.isfinite(intensity)] = 0
    return vector_preview(rgb, east, north)
