"""Scale references drawn on the terrain sandbox's and the Earth atlas's maps: a faint grid of the world
cells a window was refined from (one line every k pixels) and a scale bar whose length is the longest of
SCALE_LENGTHS_KM that fits in a quarter of the image width."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

SCALE_LENGTHS_KM = (1, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000)


def grid(rgb, step_px, phase_x=0.0, phase_y=0.0, alpha=0.18, colour=(255, 255, 255)):
    """Lines every step_px pixels, at x = phase_x + n step_px and y = phase_y + n step_px (fractional phases
    are rounded), blended at `alpha` over a float or uint8 RGB array; returns float RGB."""
    out = np.asarray(rgb, dtype=float).copy()
    rows, columns = out.shape[:2]
    xs = np.round(np.arange(phase_x % step_px, columns, step_px)).astype(int)
    ys = np.round(np.arange(phase_y % step_px, rows, step_px)).astype(int)
    colour = np.array(colour, dtype=float)
    for x in xs[(xs >= 0) & (xs < columns)]:
        out[:, x] = out[:, x] * (1 - alpha) + colour * alpha
    for y in ys[(ys >= 0) & (ys < rows)]:
        out[y, :] = out[y, :] * (1 - alpha) + colour * alpha
    return out


def scale_length_km(km_per_px, width_px, fraction=0.25):
    """The longest SCALE_LENGTHS_KM whose bar fits in fraction x the width (the shortest if none fits)."""
    fitting = [length for length in SCALE_LENGTHS_KM if length / km_per_px <= fraction * width_px]
    return fitting[-1] if fitting else SCALE_LENGTHS_KM[0]


def scale_bar(rgb, km_per_px, note=''):
    """A scale bar at the bottom left: a white bar with a dark edge, end ticks and its length in km (`note`
    after it, e.g. the latitude the east-west pixel size holds at). Returns uint8 RGB."""
    image = Image.fromarray(np.clip(np.asarray(rgb), 0, 255).astype(np.uint8))
    rows, columns = image.height, image.width
    length = scale_length_km(km_per_px, columns)
    bar = max(1, int(round(length / km_per_px)))
    size = max(12, columns // 70)
    try:
        font = ImageFont.load_default(size=size)
    except TypeError:  # Pillow without FreeType sizes
        font = ImageFont.load_default()
    draw = ImageDraw.Draw(image)
    margin, thickness = max(8, columns // 80), max(3, size // 4)
    x0, y1 = margin, rows - margin
    y0 = y1 - thickness
    label = f'{length:g} km' + (f'  {note}' if note else '')
    draw.rectangle((x0 - 1, y0 - 1, x0 + bar + 1, y1 + 1), fill=(20, 20, 20))
    draw.rectangle((x0, y0, x0 + bar, y1), fill=(250, 250, 250))
    for x in (x0, x0 + bar):
        draw.rectangle((x - 1, y0 - thickness, x + 1, y1), fill=(250, 250, 250), outline=(20, 20, 20))
    text_y = y0 - thickness - size - 4
    for dx, dy in ((-1, -1), (1, -1), (-1, 1), (1, 1)):
        draw.text((x0 + dx, text_y + dy), label, font=font, fill=(20, 20, 20))
    draw.text((x0, text_y), label, font=font, fill=(250, 250, 250))
    return np.asarray(image)
