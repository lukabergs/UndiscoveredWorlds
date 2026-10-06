"""Origin-independent numeric raster colouring. No run, module or file-layout knowledge."""
from pathlib import Path
import unicodedata

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


# Framed maps: what a reader needs to judge a raster without the palette file.
FRAME_VERSION = "framed-map-v1"
INK, PAPER, GRID = (225, 232, 237), (16, 23, 32), (255, 255, 255, 70)
LEFT, TOP, RIGHT, AXIS = 46, 30, 14, 22


def _plain(text):
    """The built-in font has no minus sign or en dash."""
    return text.replace("−", "-").replace("–", "-")


def _legend_height(legend, width, font):
    if legend is None:
        return 0
    if legend["kind"] == "classes":
        x, rows = LEFT, 1
        for label, _ in legend["entries"]:
            step = 22 + round(font.getlength(str(label))) + 14
            if x + step > width - RIGHT and x > LEFT:
                x, rows = LEFT, rows + 1
            x += step
        return 12 + 20 * rows + 18
    return 48 + 18


def _draw_legend(draw, canvas, legend, top, width, font):
    if legend["kind"] == "classes":
        x, y = LEFT, top + 8
        for label, colour in legend["entries"]:
            step = 22 + round(font.getlength(str(label))) + 14
            if x + step > width - RIGHT and x > LEFT:
                x, y = LEFT, y + 20
            draw.rectangle([x, y, x + 16, y + 13], fill=tuple(colour), outline=INK)
            draw.text((x + 22, y - 1), _plain(str(label)), font=font, fill=INK)
            x += step
        y += 20
    else:
        stops = legend["stops"]
        low, high = stops[0][0], stops[-1][0]
        length = min(width - LEFT - RIGHT, 640)
        bar = colourize(np.linspace(low, high, length)[None, :].repeat(14, 0), stops)
        canvas.paste(Image.fromarray(bar), (LEFT, top + 8))
        draw.rectangle([LEFT, top + 8, LEFT + length - 1, top + 21], outline=INK)
        edge = -1e9
        for value in legend.get("ticks") or np.linspace(low, high, 5):
            x = LEFT + (value - low) / (high - low) * (length - 1)
            draw.line([x, top + 21, x, top + 25], fill=INK)
            text = f"{value:g}"
            start = x - font.getlength(text) / 2
            if start > edge + 6:  # a label that would run into the previous one is left out; its tick stays
                draw.text((start, top + 27), text, font=font, fill=INK)
                edge = start + font.getlength(text)
        if legend.get("units"):
            draw.text((LEFT + length + 12, top + 7), _plain(legend["units"]), font=font, fill=INK)
        y = top + 46
    if legend.get("note"):
        draw.text((LEFT, y), _plain(legend["note"]), font=font, fill=INK)


def _draw_arrows(draw, east, north, scale, columns=36):
    """Direction arrows on a regular subset of cells; length follows speed up to its 90th percentile."""
    rows, width = east.shape
    stride = max(1, round(width / columns))
    speed = np.hypot(east, north)
    finite = speed[np.isfinite(speed)]
    if not finite.size or not finite.max() > 0:
        return None
    full = float(np.percentile(finite, 90)) or float(finite.max())
    reach = stride * scale * 0.9
    for row in range(stride // 2, rows, stride):
        for column in range(stride // 2, width, stride):
            u, v, s = east[row, column], north[row, column], speed[row, column]
            if not np.isfinite(s) or s <= 0:
                continue
            length = reach * min(1.0, s / full)
            if length < 2:
                continue
            x, y = LEFT + (column + .5) * scale, TOP + (row + .5) * scale
            dx, dy = u / s * length, -v / s * length  # rows run north to south
            tail, tip = (x - dx / 2, y - dy / 2), (x + dx / 2, y + dy / 2)
            head = max(3.0, length * .35)
            angle = np.arctan2(dy, dx)
            wings = [(tip[0] - head * np.cos(angle + turn), tip[1] - head * np.sin(angle + turn)) for turn in (.45, -.45)]
            for colour, line in (((0, 0, 0), 3), ((255, 255, 255), 1)):
                draw.line([tail, tip], fill=colour, width=line)
                draw.line([wings[0], tip, wings[1]], fill=colour, width=line)
    return full


MOUTH, OUTLET, DARK = (255, 255, 255), (255, 0, 255), (16, 16, 16)


def _position(lon, lat, columns, rows, scale):
    return (lon + 180) / 360 * columns * scale, (90 - lat) / 180 * rows * scale


def _ascii(text):
    """The built-in font has Latin letters without accents: fold them (Paraná -> Parana)."""
    return unicodedata.normalize("NFKD", _plain(text)).encode("ascii", "ignore").decode("ascii")


def _label_box(x, y, length, radius, height, place):
    """(left, top, right, bottom) of a label beside a ring: right, left, above or below it."""
    if place == "right":
        left, top = x + radius + 3, y - 7
    elif place == "left":
        left, top = x - radius - 3 - length, y - 7
    else:
        left, top = x - length / 2, y - radius - 15 if place == "above" else y + radius + 2
    return left, top, left + length, top + 13


def draw_markers(canvas, columns, rows, scale, rings=(), squares=(), links=(), origin=(LEFT, TOP)):
    """Point markers on the raster of a framed map (cell size `scale` pixels, raster at `origin`).

    rings: (lon, lat, label), a white ring with a dark outline and its label; squares: (lon, lat),
    a filled magenta square; links: (lon, lat, lon, lat), a thin line drawn the short way in
    longitude. Positions in degrees, columns east from 180 W. Markers keep a few pixels at
    small cell sizes and are drawn again one raster width east and west, so the left and
    right edges join; everything is clipped to the raster. A label goes right of its ring,
    left, above or below it: the first of those places inside the raster that covers the
    least area of the markers and of the labels placed before it.
    """
    width, height = columns * scale, rows * scale
    layer = Image.new("RGBA", (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(layer)
    radius, half = min(max(2 * scale, 6), 10), min(max(scale, 3), 5)
    shifts = (0, -width, width)
    for lon, lat, end_lon, end_lat in links:
        x, y = _position(lon, lat, columns, rows, scale)
        reach = ((end_lon - lon + 180) % 360 - 180) / 360 * width
        end_y = _position(end_lon, end_lat, columns, rows, scale)[1]
        for shift in shifts:
            draw.line([x + shift, y, x + shift + reach, end_y], fill=DARK, width=3)
            draw.line([x + shift, y, x + shift + reach, end_y], fill=MOUTH, width=1)
    for lon, lat, _ in rings:
        x, y = _position(lon, lat, columns, rows, scale)
        for shift in shifts:
            draw.ellipse([x + shift - radius - 1, y - radius - 1, x + shift + radius + 1, y + radius + 1], outline=DARK, width=4)
            draw.ellipse([x + shift - radius, y - radius, x + shift + radius, y + radius], outline=MOUTH, width=2)
    for lon, lat in squares:
        x, y = _position(lon, lat, columns, rows, scale)
        for shift in shifts:
            draw.rectangle([x + shift - half, y - half, x + shift + half, y + half], fill=OUTLET, outline=DARK)
    font = ImageFont.load_default(size=11)
    placed = []
    for lon, lat, _ in rings:
        x, y = _position(lon, lat, columns, rows, scale)
        placed.append((x - radius, y - radius, x + radius, y + radius))
    for lon, lat in squares:
        x, y = _position(lon, lat, columns, rows, scale)
        placed.append((x - half, y - half, x + half, y + half))
    covered = lambda box: sum(max(0, min(box[2], other[2]) - max(box[0], other[0])) *
                              max(0, min(box[3], other[3]) - max(box[1], other[1])) for other in placed)
    for lon, lat, label in rings:
        if not label:
            continue
        x, y = _position(lon, lat, columns, rows, scale)
        text = _ascii(str(label))
        length = font.getlength(text)
        boxes = [_label_box(x, y, length, radius, height, place) for place in ("right", "left", "above", "below")]
        inside = [box for box in boxes if box[0] >= 0 and box[2] <= width and box[1] >= 0 and box[3] <= height] or boxes[:2]
        box = min(inside, key=covered)
        placed.append(box)
        draw.text((box[0], box[1]), text, font=font, fill=MOUTH, stroke_width=1, stroke_fill=DARK)
    canvas.paste(layer, tuple(origin), layer)
    return canvas


def frame(rgb, title=None, legend=None, arrows=None, width=1024, markers=None):
    """A map for reading: the raster enlarged without smoothing to at least `width` pixels,
    a graticule every 30 degrees with its labels, the title, and the legend.

    The grid is cell-centred with rows from north to south and columns east from 180 W.
    legend: {"kind": "scale", "stops": [(value, rgb), ...], "units", "ticks", "note"} or
    {"kind": "classes", "entries": [(label, rgb), ...], "note"}. arrows: (east, north)
    components on the grid of the raster. markers: the keyword arguments of draw_markers
    (rings, squares, links), drawn over the raster and the graticule.
    """
    rgb = np.asarray(rgb, dtype=np.uint8)
    rows, columns = rgb.shape[:2]
    scale = max(1, -(-width // columns))
    raster = Image.fromarray(rgb).resize((columns * scale, rows * scale), Image.Resampling.NEAREST)
    font = ImageFont.load_default(size=13)
    legend = dict(legend) if legend else None
    total = LEFT + raster.width + RIGHT
    canvas = Image.new("RGB", (total, TOP + raster.height + AXIS + _legend_height(legend, total, font)), PAPER)
    canvas.paste(raster, (LEFT, TOP))
    draw = ImageDraw.Draw(canvas, "RGBA")
    for longitude in range(-180, 181, 30):
        x = LEFT + (longitude + 180) / 360 * raster.width
        if -180 < longitude < 180:
            draw.line([x, TOP, x, TOP + raster.height - 1], fill=GRID)
        text = "0" if longitude == 0 else "180" if abs(longitude) == 180 else f"{abs(longitude)}{'E' if longitude > 0 else 'W'}"
        draw.text((min(x - font.getlength(text) / 2, total - 2 - font.getlength(text)), TOP + raster.height + 4), text, font=font, fill=INK)
    for latitude in range(-90, 91, 30):
        y = TOP + (90 - latitude) / 180 * raster.height
        if -90 < latitude < 90:
            draw.line([LEFT, y, LEFT + raster.width - 1, y], fill=GRID)
        text = "0" if latitude == 0 else f"{abs(latitude)}{'N' if latitude > 0 else 'S'}"
        draw.text((LEFT - 6 - font.getlength(text), min(max(y - 8, TOP), TOP + raster.height - 14)), text, font=font, fill=INK)
    if arrows is not None:
        full = _draw_arrows(draw, np.asarray(arrows[0], float), np.asarray(arrows[1], float), scale)
        if legend is not None and full is not None:
            note = f"arrows: direction; full length at {full:.3g} or more"
            legend["note"] = f"{legend['note']}; {note}" if legend.get("note") else note
    if markers:
        draw_markers(canvas, columns, rows, scale, **markers)
    if title:
        draw.text((LEFT, 7), _plain(title), font=ImageFont.load_default(size=15), fill=INK)
    if legend is not None:
        _draw_legend(draw, canvas, legend, TOP + raster.height + AXIS, total, font)
    return canvas


def save_frames(path, images, durations_ms):
    """Write a looping APNG from framed images of one size."""
    if not images or len(images) != len(durations_ms) or any(duration <= 0 for duration in durations_ms):
        raise ValueError("Animation requires one positive duration per frame")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(path, format="PNG", save_all=True, append_images=images[1:],
                   duration=durations_ms, loop=0, disposal=0, blend=0)
    return {"format": "APNG", "frames": len(images), "width": images[0].width, "height": images[0].height,
            "duration_ms": sum(durations_ms)}
