"""Hypsographic curves of terrain sandbox worlds against Earth's (ETOPO 2022 at 2048 columns, the cached
reference of scripts/uwtools/tectonics_earth_reference.py), with the summary numbers of compare_to_earth.

  uv run python scripts/sandbox/hypsometry.py --name sea-v7 --labels sea-v3b,sea-v7
      [--array run62=runs/archive/fields/tectonics/elevation/final/62.npy]

A label reads out/scratch/sandbox/<label>/global-elevation.npy (the world after the physics, or the procedural
surface of a --no-physics run); --array adds any 2:1 elevation array (metres from the sea level). Every curve
is weighted by spherical cell area (measure_earth); land is elevation above 0 (Earth's mask also counts land
below the sea and lakes as land). Writes out/scratch/sandbox/hypsometry-<name>/index.html, hypsometry.png and
summary.json.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from uwtools.tectonics_earth_reference import compare_to_earth, load_earth_reference, measure_earth  # noqa: E402

SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'
SURFACE, TEXT, MUTED, GRID = '#1a1a19', '#f0efec', '#c3c2b7', '#3a3a38'
SERIES = ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300']  # categorical slots, dark mode

PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Hypsometry {name}</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}img{{max-width:100%}}
table{{border-collapse:collapse;margin:12px 0}}td,th{{border-bottom:1px solid #333;padding:4px 10px;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body><h1>Hypsometry: {name}</h1>
<p>Fraction of the surface below each elevation, weighted by spherical area. Earth: {earth}. Land is
elevation above 0 m in the worlds (Earth's mask also counts land below the sea and lakes). The quantile gap
is the mean absolute difference from Earth's curve over all quantiles.</p>
<img src="hypsometry.png"><table><tr><th>world</th><th>land %</th><th>land below 50 m, % of land</th>
<th>land median m</th><th>ocean median depth m</th><th>mean elevation m</th><th>quantile gap m</th></tr>{rows}</table>
</body></html>"""


def measure(elevation):
    elevation = np.asarray(elevation, dtype=float)
    metrics = measure_earth(elevation, elevation > 0.0)
    land = elevation[elevation > 0.0]
    metrics['land_below_50'] = float((land < 50.0).mean()) if land.size else 0.0
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--name', required=True)
    parser.add_argument('--labels', default='')
    parser.add_argument('--array', action='append', default=[], metavar='NAME=PATH')
    arguments = parser.parse_args()

    reference = load_earth_reference(ROOT, 2048)
    earth = reference['metrics']
    worlds = []
    for label in [item for item in arguments.labels.split(',') if item]:
        worlds.append((label, np.load(SANDBOX / label / 'global-elevation.npy')))
    for item in arguments.array:
        name, _, path = item.partition('=')
        worlds.append((name, np.load(ROOT / path)))
    if len(worlds) > len(SERIES):
        raise SystemExit(f'at most {len(SERIES)} worlds per chart')

    figure, axes = plt.subplots(1, 2, figsize=(13, 5.2), facecolor=SURFACE)
    panels = [(axes[0], (-11000, 9000), 'Whole surface'), (axes[1], (-6500, 1500), 'Sea floor and lowlands')]
    q = np.asarray(earth['hypsography']['fraction_below'])
    summary, rows = {}, []
    curves = [('Earth (ETOPO 2048)', np.asarray(earth['hypsography']['elevation_m']), TEXT, 2.5)]
    for (name, elevation), colour in zip(worlds, SERIES):
        metrics = measure(elevation)
        comparison = compare_to_earth(metrics, reference)
        summary[name] = {'shape': list(np.shape(elevation)), 'land_fraction': metrics['land_fraction'],
                         'land_below_50_m_share': metrics['land_below_50'],
                         'land_median_m': metrics['land_elevation_m']['median'] if metrics['land_elevation_m'] else None,
                         'ocean_median_depth_m': metrics['ocean_depth_m']['median'] if metrics['ocean_depth_m'] else None,
                         **comparison}
        curves.append((name, np.asarray(metrics['hypsography']['elevation_m']), colour, 2.0))
    for axis, limits, title in panels:
        axis.set_facecolor(SURFACE)
        for name, values, colour, width in curves:
            axis.plot(q * 100.0, values, color=colour, linewidth=width, label=name)
        axis.axhline(0.0, color=MUTED, linewidth=0.8, linestyle=':')
        axis.set_ylim(*limits)
        axis.set_xlim(0, 100)
        axis.set_title(title, color=TEXT, fontsize=12, loc='left')
        axis.set_xlabel('surface below this elevation (%)', color=MUTED)
        axis.set_ylabel('elevation (m)', color=MUTED)
        axis.grid(True, color=GRID, linewidth=0.6)
        axis.tick_params(colors=MUTED)
        for spine in axis.spines.values():
            spine.set_color(GRID)
    legend = axes[0].legend(loc='upper left', frameon=False)
    for text in legend.get_texts():
        text.set_color(TEXT)
    figure.tight_layout()
    directory = SANDBOX / f'hypsometry-{arguments.name}'
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(directory / 'hypsometry.png', dpi=110, facecolor=SURFACE)

    land_earth = 100 * earth['land_fraction']
    rows.append(f'<tr><td>Earth</td><td>{land_earth:.1f}</td><td>7.1</td><td>{earth["land_elevation_m"]["median"]:.0f}</td>'
                f'<td>{earth["ocean_depth_m"]["median"]:.0f}</td><td>{earth["elevation_m"]["mean"]:.0f}</td><td>0</td></tr>')
    for name, values in summary.items():
        rows.append(f'<tr><td>{name}</td><td>{100 * values["land_fraction"]:.1f}</td>'
                    f'<td>{100 * values["land_below_50_m_share"]:.1f}</td><td>{values["land_median_m"]:.0f}</td>'
                    f'<td>{values["ocean_median_depth_m"]:.0f}</td>'
                    f'<td>{values["mean_elevation_delta_m"] + earth["elevation_m"]["mean"]:.0f}</td>'
                    f'<td>{values["hypsography_mean_absolute_quantile_gap_m"]:.0f}</td></tr>')
    (directory / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n', encoding='utf-8')
    (directory / 'index.html').write_text(PAGE.format(name=arguments.name, earth=reference.get('title', 'ETOPO 2022, 2048 x 1024'),
                                                      rows=''.join(rows)), encoding='utf-8')
    print(directory / 'index.html')


if __name__ == '__main__':
    main()
