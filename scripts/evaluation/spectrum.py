"""The directional anisotropy of a terrain's 2D power spectrum (multiscale plan S3.7): a measure of the
striping of the regional relief maps that does not depend on the hillshade's light.

The field is the land's elevation above the sea level, zero at sea (land meets the sea near zero, so the
coast is a kink, not a step), detrended by the least-squares plane over the whole grid and tapered by a
periodic Hann window in both axes. Its spectrum is weighted by |k|^2: the slope spectrum, the power of
the gradient that a hillshade draws, so every scale counts as the eye sees it rather than the largest.

Directions are those of the crests (perpendicular to the wavevector), in four sectors of 45 degrees on
the ground: crests E-W, SW-NE, N-S, NW-SE. Light from the north-west (azimuth 315) brightens and darkens
the slopes of crests running SW-NE most, so a SW-NE excess that the spectrum also shows is the
terrain's, one it does not show is the light's.

Pixels are square in degrees (as maps.pixel_sizes_m): the north-south size dy and the east-west size
dy * cos(latitude). The spectrum takes one east-west size for the whole grid, at the window's central
latitude; the rows' true sizes differ from it by `ew_spread` (reported), which tilts the sector edges by
at most atan(spread) and scales east-west wavenumbers by as much.
"""
import math

import numpy as np

# Crest orientation -> the bearing of the wavevector (degrees clockwise from north, modulo 180).
SECTORS = {'E-W': 0.0, 'NW-SE': 45.0, 'N-S': 90.0, 'SW-NE': 135.0}
SECTOR_HALF_WIDTH = 22.5
# The first-order points of a parent lattice in (east, north) multiples of its frequency, by crest orientation.
LATTICE = {'N-S': (1, 0), 'E-W': (0, 1), 'NW-SE': (1, 1), 'SW-NE': (1, -1)}
# The parent band in |k| times the parent period: it holds the lattice's axis points (1) and diagonal points (sqrt 2).
PARENT_BAND = (0.6, 1.6)
# The ring around a lattice point whose mean is the peak's background, in |k| relative to the point's.
PEAK_RING = 1.15
# The overall band starts at this many periods over the grid's shorter side (the largest scales are the window's shape).
OVERALL_PERIODS = 4.0


def periodic_hann(count):
    """The periodic Hann window of `count` points (its DFT has three nonzero terms: an exact-bin sinusoid
    leaks only into the bins beside it)."""
    return 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(count) / count) if count > 1 else np.ones(count)


def prepared(elevation, land):
    """The land's elevation (zero at sea), minus the least-squares plane of the whole grid, times the
    two-dimensional periodic Hann window."""
    values = np.where(land, np.asarray(elevation, dtype=np.float64), 0.0)
    rows, columns = values.shape
    y, x = np.indices(values.shape, dtype=np.float64)
    basis = np.stack([np.ones(values.size), x.ravel(), y.ravel()], axis=1)
    coefficients, *_ = np.linalg.lstsq(basis, values.ravel(), rcond=None)
    values = values - (basis @ coefficients).reshape(values.shape)
    return values * periodic_hann(rows)[:, None] * periodic_hann(columns)[None, :]


def slope_spectrum(elevation, land, dy_km, dx_km):
    """(slope power |k|^2 |F|^2, |k| in cycles per km, bearing of the wavevector in degrees modulo 180,
    weight) on rfft2's half plane (columns run east, rows south). The weight is the number of bins of the
    full plane a bin stands for: 2 (itself and its conjugate), 1 in column 0 and the Nyquist column of an
    even column count, whose conjugates are stored in the same column."""
    values = prepared(elevation, land)
    rows, columns = values.shape
    power = np.abs(np.fft.rfft2(values)) ** 2
    k_east = np.broadcast_to(np.fft.rfftfreq(columns)[None, :] / dx_km, power.shape)
    k_north = np.broadcast_to(-np.fft.fftfreq(rows)[:, None] / dy_km, power.shape)
    magnitude = np.hypot(k_east, k_north)
    bearing = np.degrees(np.arctan2(k_east, k_north)) % 180.0
    weight = np.full(power.shape, 2.0)
    weight[:, 0] = 1.0
    if columns % 2 == 0:
        weight[:, -1] = 1.0
    return power * magnitude ** 2, magnitude, bearing, weight


def sector_of(bearing):
    """The index into SECTORS of each bearing's crest orientation (degrees modulo 180)."""
    labels = np.full(bearing.shape, -1, dtype=np.int8)
    for index, centre in enumerate(SECTORS.values()):
        offset = (bearing - centre + 90.0) % 180.0 - 90.0
        labels[(offset >= -SECTOR_HALF_WIDTH) & (offset < SECTOR_HALF_WIDTH)] = index
    return labels


def weighted_mean(slope, weight, cells):
    total = float(weight[cells].sum())
    return float((slope[cells] * weight[cells]).sum()) / total if total > 0 else None


def sector_means(slope, magnitude, labels, weight, low, high):
    """Mean slope power per crest orientation over |k| in [low, high) on the full plane, and the full
    plane's bin count of each."""
    inside = (magnitude >= low) & (magnitude < high)
    means, counts = {}, {}
    for index, name in enumerate(SECTORS):
        cells = inside & (labels == index)
        counts[name] = int(weight[cells].sum())
        means[name] = weighted_mean(slope, weight, cells)
    return means, counts


def ratios(means):
    """(SW-NE over NW-SE, E-W over N-S) of the sector means; None where a denominator is absent or zero."""
    divide = lambda a, b: None if means[a] is None or not means[b] else means[a] / means[b]
    return divide('SW-NE', 'NW-SE'), divide('E-W', 'N-S')


def lattice_peaks(slope, magnitude, weight, shape, dy_km, period_km):
    """The parent lattice's first-order points (frequency dy / period cycles per pixel along both pixel
    axes: parent cells and pixels are both square in degrees) on a grid of `shape` (rows, columns): for
    each crest orientation the mean slope power of the 3 x 3 bins around the point over the mean of the
    ring of |k| within PEAK_RING of the point's (the 3 x 3 bins of every lattice point left out; means on
    the full plane), and the point's period in km (1 / |k|). None where the point lies beyond the grid's
    frequencies."""
    rows, columns = shape
    half = slope.shape[1]
    frequency = dy_km / period_km
    points, excluded = {}, np.zeros(slope.shape, dtype=bool)
    for name, (east, north) in LATTICE.items():
        column, row = int(round(east * frequency * columns)), int(round(-north * frequency * rows)) % rows
        if column > half - 2 or abs(int(round(north * frequency * rows))) >= rows // 2 or (column == 0 and row == 0):
            points[name] = None
            continue
        neighbourhood = (np.arange(row - 1, row + 2) % rows)[:, None], np.arange(max(column - 1, 0), column + 2)[None, :]
        points[name] = (row, column, neighbourhood)
        excluded[neighbourhood] = True
        if column == 0:  # crests E-W: the conjugate point's neighbourhood is stored in the same columns
            excluded[(-np.arange(row - 1, row + 2)) % rows, 0:2] = True
    peaks = {}
    for name, point in points.items():
        if point is None:
            peaks[name] = None
            continue
        row, column, neighbourhood = point
        radius = magnitude[row, column]
        ring = (magnitude >= radius / PEAK_RING) & (magnitude <= radius * PEAK_RING) & ~excluded
        background = weighted_mean(slope, weight, ring) or 0.0
        cells = np.zeros(slope.shape, dtype=bool)
        cells[neighbourhood] = True
        strength = weighted_mean(slope, weight, cells) / background if background > 0 else None
        peaks[name] = {'strength': strength, 'period_km': float(1.0 / radius)}
    return peaks


def spectral_anisotropy(elevation, land, dy_km, dx_km, parents=(), band=PARENT_BAND):
    """The anisotropy of the land's slope spectrum (module docstring): 'overall' over |k| from
    OVERALL_PERIODS over the grid's shorter side to the Nyquist circle (half a cycle per north-south or
    east-west pixel, whichever is longer); and for each (label, parent period in km) of `parents` the band
    |k| * period in `band` (clipped at the Nyquist circle) with the lattice peaks of lattice_peaks. Each
    entry gives the sector means, the diagonal ratio (crests SW-NE over NW-SE) and the axis ratio (crests
    E-W over N-S). A parent period under 2 pixels or over a quarter of the grid's shorter side is not
    resolved (entry 'resolved' False)."""
    land = np.asarray(land, dtype=bool)
    slope, magnitude, bearing, weight = slope_spectrum(elevation, land, dy_km, dx_km)
    labels = sector_of(bearing)
    rows, columns = land.shape
    shorter = min(rows * dy_km, columns * dx_km)
    nyquist = 0.5 / max(dy_km, dx_km)
    low = OVERALL_PERIODS / shorter
    means, counts = sector_means(slope, magnitude, labels, weight, low, nyquist)
    diagonal, axis = ratios(means)
    result = {'land_share': float(land.mean()) if land.size else 0.0,
              'overall': {'band_km': [1.0 / nyquist, 1.0 / low], 'diagonal_ratio': diagonal, 'axis_ratio': axis,
                          'sector_means': means, 'sector_bins': counts},
              'levels': []}
    for label, period in parents:
        entry = {'level': label, 'parent_period_km': float(period)}
        if period < 2.0 * max(dy_km, dx_km) or period > shorter / OVERALL_PERIODS:
            entry['resolved'] = False
            result['levels'].append(entry)
            continue
        low_k, high_k = band[0] / period, min(band[1] / period, nyquist)
        means, counts = sector_means(slope, magnitude, labels, weight, low_k, high_k)
        diagonal, axis = ratios(means)
        peaks = lattice_peaks(slope, magnitude, weight, land.shape, dy_km, period)
        strongest = max((name for name in peaks if peaks[name] and peaks[name]['strength'] is not None),
                        key=lambda name: peaks[name]['strength'], default=None)
        entry.update({'resolved': True, 'band_km': [1.0 / high_k, 1.0 / low_k], 'diagonal_ratio': diagonal,
                      'axis_ratio': axis, 'sector_means': means, 'sector_bins': counts, 'peaks': peaks,
                      'peak': None if strongest is None else dict(peaks[strongest], crests=strongest)})
        result['levels'].append(entry)
    return result


def ew_spread(north, south, mid_cos):
    """The largest relative difference of a row's east-west pixel size from the one at the central latitude."""
    cosines = [math.cos(math.radians(latitude)) for latitude in (north, south)]
    edge = max(abs(value - mid_cos) for value in cosines)
    if south < 0.0 < north:
        edge = max(edge, abs(1.0 - mid_cos))
    return edge / mid_cos if mid_cos > 0 else float('inf')
