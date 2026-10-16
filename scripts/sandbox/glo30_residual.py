"""The Earth's detail below a parent cell (terrain roadmap T1.5): how much relief a 100 m surface has that its 611 m
box means do not carry, measured on GLO-30 windows, so nested_window.py's detail amplitude and spectrum are
calibrated against Earth data (a calibration reference, recorded in configs/sandbox/nested.json, not a forcing).

Per window (glo30_windows.py's STEM.npy, STEM-sea.npy, STEM.json, refs/processed/earth-atlas):
  1. The coarse surface: box means over factor x factor pixels (6: the windows have the nested child's pixels,
     101.8 x 63.6 m, so a block is the parent's 611 x 382 m cell in km, as the parent's cells relate to the child's).
     Sea pixels keep their ETOPO depth, as the parent's surface keeps its sea floor.
  2. Prolonged back with nested_window.prolong (Catmull-Rom, clamped to the 2 x 2 bracketing coarse cells).
  3. The residual: the surface minus the prolonged coarse surface, on land; its "detail" part is the residual minus
     its mean over each block's land pixels (what nested_window.py's detail is: zero land mean per parent cell; the
     block-mean part is the prolongation's error, which the child's projection removes). Statistics leave out the
     outer `edge` coarse cells (the prolongation's clamped indices).
  4. The detail's rms against the coarse surface's local relief (nested_window.local_relief: max - min over 3 x 3
     coarse cells, sea at 0 m, prolonged as nested_window.py prolongs it) and against its slope (central differences
     of the coarse surface, sea at 0 m, prolonged) times the coarse cell's north-south size: land pixels in `bins`
     quantile bins of the predictor, the rms per bin, a least-squares line rms = a x predictor + b through the bins
     (equal counts), and its R^2; per window and both windows pooled.
  5. Spectra: tiles of tile_rows x tile_columns pixels wholly on land (inside the edge), a plane removed, a Hann
     window, the 2D power spectrum radially averaged in true wavenumber (the pixels' km sizes), log-spaced bins;
     the slope beta of log power against log wavelength between the two wavelengths (default 1.2 and 0.25 km: the
     north-south Nyquist wavelength is 0.204 km and the box-averaged pixels attenuate the last octave), for the
     surface and for the detail. A self-affine surface with Hurst exponent H has beta = 2 + 2H, an octave gain of
     2^-H in fBm.

Writes the JSON (default refs/processed/earth-atlas/glo30-residual.json, beside its inputs).

  uv run python scripts/sandbox/glo30_residual.py --items glo30-sognefjord-98km,glo30-patagonia-98km
"""
import argparse
import datetime
import json
import math
import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS / 'sandbox'))
import nested_window as nw  # noqa: E402

ATLAS = nw.ROOT / 'refs' / 'processed' / 'earth-atlas'


def load(stem):
    elevation = np.load(ATLAS / f'{stem}.npy').astype(np.float64)
    sea = np.load(ATLAS / f'{stem}-sea.npy').astype(bool)
    receipt = json.loads((ATLAS / f'{stem}.json').read_text(encoding='utf-8'))
    return elevation, sea, receipt


def fields(elevation, sea, factor, ns_km, ew_km, window):
    """Per pixel: the detail (zero land mean per block), the full residual, the prolonged relief and slope x cell."""
    rows, columns = (n // factor * factor for n in elevation.shape)
    elevation, sea = elevation[:rows, :columns], sea[:rows, :columns]
    land = ~sea
    coarse = nw.block_mean(elevation, factor)
    cr, cc = coarse.shape
    prolonged = nw.prolong(coarse, factor, 0, 0, cr, cc)
    residual = elevation - prolonged
    share = nw.block_mean(land.astype(float), factor)
    with np.errstate(invalid='ignore', divide='ignore'):
        per_land = np.where(share > 0, nw.block_mean(np.where(land, residual, 0.0), factor) / share, 0.0)
    detail = residual - nw.nearest(per_land, factor, 0, 0, cr, cc)
    relief = np.maximum(nw.prolong(nw.local_relief(coarse, window), factor, 0, 0, cr, cc), 0.0)
    gy, gx = np.gradient(np.maximum(coarse, 0.0))
    slope = np.hypot(gx / (ew_km * factor * 1000.0), gy / (ns_km * factor * 1000.0))
    slope_cell = np.maximum(nw.prolong(slope, factor, 0, 0, cr, cc), 0.0) * ns_km * factor * 1000.0
    return {'land': land, 'elevation': elevation, 'detail': detail, 'residual': residual, 'relief': relief,
            'slope_cell': slope_cell, 'coarse_shape': [cr, cc], 'mean_error': per_land}


def binned_fit(predictor, values, bins):
    """rms of `values` in quantile bins of `predictor`; the line rms = a x + b through the bins and its R^2."""
    edges = np.unique(np.percentile(predictor, np.linspace(0, 100, bins + 1)))
    index = np.clip(np.searchsorted(edges, predictor, side='right') - 1, 0, edges.size - 2)
    x = np.array([predictor[index == i].mean() for i in range(edges.size - 1)])
    rms = np.array([math.sqrt((values[index == i] ** 2).mean()) for i in range(edges.size - 1)])
    a, b = np.polyfit(x, rms, 1)
    fitted = a * x + b
    r2 = 1.0 - ((rms - fitted) ** 2).sum() / ((rms - rms.mean()) ** 2).sum()
    through_origin = float((x * rms).sum() / (x * x).sum())
    return {'a': round(float(a), 4), 'b_m': round(float(b), 2), 'r2': round(float(r2), 4),
            'a_through_origin': round(through_origin, 4),
            'bins': [[round(float(u), 2), round(float(v), 2)] for u, v in zip(x, rms)]}


def tile_spectrum(values, ns_km, ew_km):
    """(wavelength km, power) of one tile: plane removed, Hann window, 2D power by |k|."""
    rows, columns = values.shape
    yy, xx = np.mgrid[0:rows, 0:columns]
    design = np.column_stack([np.ones(values.size), xx.ravel(), yy.ravel()])
    coefficients, *_ = np.linalg.lstsq(design, values.ravel(), rcond=None)
    flat = values - (design @ coefficients).reshape(values.shape)
    window = np.outer(np.hanning(rows), np.hanning(columns))
    power = np.abs(np.fft.fft2(flat * window)) ** 2
    ky = np.fft.fftfreq(rows, d=ns_km)
    kx = np.fft.fftfreq(columns, d=ew_km)
    k = np.hypot(*np.meshgrid(kx, ky))
    keep = k > 0
    return 1.0 / k[keep], power[keep]


def spectrum_slope(values, land, ns_km, ew_km, tile, edge_px, longest_km, shortest_km, bins=16):
    """beta: the slope of log10 power against log10 wavelength between the wavelengths, over the land tiles."""
    rows, columns = values.shape
    edges = np.logspace(math.log10(shortest_km), math.log10(longest_km), bins + 1)
    sums, counts, tiles = np.zeros(bins), np.zeros(bins), 0
    for r in range(edge_px, rows - edge_px - tile[0] + 1, tile[0]):
        for c in range(edge_px, columns - edge_px - tile[1] + 1, tile[1]):
            if not land[r:r + tile[0], c:c + tile[1]].all():
                continue
            wavelength, power = tile_spectrum(values[r:r + tile[0], c:c + tile[1]], ns_km, ew_km)
            index = np.searchsorted(edges, wavelength) - 1
            inside = (index >= 0) & (index < bins)
            np.add.at(sums, index[inside], power[inside])
            np.add.at(counts, index[inside], 1)
            tiles += 1
    if tiles == 0:
        return {'tiles': 0}
    used = counts > 0
    centres = np.sqrt(edges[:-1] * edges[1:])[used]
    mean_power = sums[used] / counts[used]
    beta, _ = np.polyfit(np.log10(centres), np.log10(mean_power), 1)
    hurst = beta / 2.0 - 1.0
    return {'tiles': tiles, 'beta': round(float(beta), 3), 'hurst': round(float(hurst), 3),
            'octave_gain': round(2.0 ** -hurst, 3)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--items', default='glo30-sognefjord-98km,glo30-patagonia-98km')
    parser.add_argument('--factor', type=int, default=6)
    parser.add_argument('--relief-window', type=int, default=3, help='coarse cells')
    parser.add_argument('--edge', type=int, default=2, help='coarse cells left out at each side')
    parser.add_argument('--bins', type=int, default=20)
    parser.add_argument('--tile', default='96,160', help='spectrum tile rows,columns (pixels)')
    parser.add_argument('--band-km', default='1.2,0.25', help='longest,shortest wavelength of the spectrum fit')
    parser.add_argument('--out', default=str(ATLAS / 'glo30-residual.json'))
    arguments = parser.parse_args()
    tile = tuple(int(value) for value in arguments.tile.split(','))
    longest, shortest = (float(value) for value in arguments.band_km.split(','))
    f = arguments.factor
    windows, pooled = {}, {key: [] for key in ('detail', 'residual', 'relief', 'slope_cell')}
    for stem in arguments.items.split(','):
        elevation, sea, receipt = load(stem)
        ns_km, ew_km = float(receipt['pixel_ns_km']), float(receipt['pixel_ew_km'])
        values = fields(elevation, sea, f, ns_km, ew_km, arguments.relief_window)
        inner = np.zeros(values['land'].shape, bool)
        e = arguments.edge * f
        inner[e:-e, e:-e] = True
        use = values['land'] & inner
        for key in pooled:
            pooled[key].append(values[key][use])
        detail = values['detail'][use]
        windows[stem] = {
            'coarse_cell_km': [round(ns_km * f, 4), round(ew_km * f, 4)], 'coarse_cells': values['coarse_shape'],
            'land_pixels': int(use.sum()),
            'residual_rms_m': round(float(np.sqrt((values['residual'][use] ** 2).mean())), 2),
            'detail_rms_m': round(float(np.sqrt((detail ** 2).mean())), 2),
            'detail_abs_p50_p90_p99_m': [round(float(v), 1) for v in np.percentile(np.abs(detail), [50, 90, 99])],
            'prolongation_block_error_rms_m': round(float(np.sqrt((values['mean_error'][
                arguments.edge:-arguments.edge, arguments.edge:-arguments.edge] ** 2).mean())), 2),
            'coarse_relief_p50_p90_m': [round(float(v), 1) for v in np.percentile(values['relief'][use], [50, 90])],
            'detail_vs_relief': binned_fit(values['relief'][use], detail, arguments.bins),
            'detail_vs_slope_x_cell': binned_fit(values['slope_cell'][use], detail, arguments.bins),
            'spectrum_surface': spectrum_slope(values['elevation'], values['land'], ns_km, ew_km, tile,
                                               e, longest, shortest),
            'spectrum_detail': spectrum_slope(np.where(values['land'], values['detail'], 0.0), values['land'], ns_km,
                                              ew_km, tile, e, longest, shortest)}
        print(stem, json.dumps(windows[stem], indent=None)[:1500], flush=True)
    joined = {key: np.concatenate(parts) for key, parts in pooled.items()}
    document = {
        'schema': 'uw-glo30-residual/v1', 'module': 'scripts/sandbox/glo30_residual.py',
        'method': __doc__.split('\n\n')[1].replace('\n', ' '),
        'factor': f, 'relief_window_coarse_cells': arguments.relief_window, 'edge_coarse_cells': arguments.edge,
        'spectrum_band_km': [longest, shortest], 'spectrum_tile_pixels': list(tile),
        'windows': windows,
        'pooled': {'detail_rms_m': round(float(np.sqrt((joined['detail'] ** 2).mean())), 2),
                   'detail_vs_relief': binned_fit(joined['relief'], joined['detail'], arguments.bins),
                   'detail_vs_slope_x_cell': binned_fit(joined['slope_cell'], joined['detail'], arguments.bins)},
        'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    out = Path(arguments.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(document, indent=1) + '\n', encoding='utf-8')
    print(json.dumps(document['pooled'], indent=None))
    print(out)


if __name__ == '__main__':
    main()
