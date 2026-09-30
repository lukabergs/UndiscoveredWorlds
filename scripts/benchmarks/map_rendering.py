"""Origin-independent numeric raster colouring. No run, module or file-layout knowledge."""
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

VERSION = "numeric-raster-v1"


def colourize(values, stops, missing=(0, 0, 0)):
    """Interpolate explicit (value, RGB) anchors; saturation affects previews only."""
    values = np.asarray(values)
    positions = np.array([v for v, _ in stops], dtype=float)
    if len(positions) < 2 or not np.all(np.diff(positions) > 0):
        raise ValueError("Colour anchors must increase strictly")
    safe = np.where(np.isfinite(values), values, 0)
    result = np.stack([np.interp(safe, positions, [c[k] for _, c in stops])
                       for k in range(3)], axis=-1)
    result[~np.isfinite(values)] = missing
    return np.floor(result + .5).astype(np.uint8)


def categories(values, colours, missing=(0, 0, 0)):
    """Explicit labels preserve class identity, without interpolating between classes."""
    values = np.asarray(values)
    result = np.empty((*values.shape, 3), dtype=np.uint8)
    result[:] = missing
    for key, colour in colours.items():
        result[values == int(key)] = colour
    return result


def save_rgb(path, rgb):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(np.asarray(rgb, dtype=np.uint8)).save(path)


def outline(mask, colour=(245, 244, 218), alpha=230):
    """RGBA outline on the true side of a Boolean mask; both grid axes wrap."""
    mask = np.asarray(mask, dtype=bool)
    edge = mask & (~np.roll(mask, 1, 0) | ~np.roll(mask, -1, 0)
                   | ~np.roll(mask, 1, 1) | ~np.roll(mask, -1, 1))
    result = np.zeros((*mask.shape, 4), dtype=np.uint8)
    result[..., :3] = colour
    result[..., 3] = np.where(edge, alpha, 0)
    return result


def save_scalar(path, values, stops, missing=(0, 0, 0)):
    save_rgb(path, colourize(values, stops, missing))


def save_animation(path, map_paths, labels, durations_ms, width=1024):
    """Write a looping APNG from ordered previews; callers supply labels and timing."""
    if not map_paths or len(map_paths) != len(labels) or len(labels) != len(durations_ms):
        raise ValueError("Animation requires one label and duration per frame")
    if any(duration <= 0 for duration in durations_ms):
        raise ValueError("Animation durations must be positive")
    frames = []
    font = ImageFont.load_default(size=18)
    for source, label in zip(map_paths, labels):
        with Image.open(source) as image:
            height = round(width * image.height / image.width)
            raster = image.convert("RGB").resize((width, height), Image.Resampling.NEAREST)
        frame = Image.new("RGB", (width, height + 40), (16, 23, 32))
        frame.paste(raster, (0, 40))
        ImageDraw.Draw(frame).text((14, 10), label, font=font, fill=(225, 232, 237))
        frames.append(frame)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(path, format="PNG", save_all=True, append_images=frames[1:],
                   duration=durations_ms, loop=0, disposal=0, blend=0)
    with Image.open(path) as animation:
        return {"format": "APNG", "frames": animation.n_frames,
                "width": animation.width, "height": animation.height,
                "duration_ms": sum(durations_ms)}
