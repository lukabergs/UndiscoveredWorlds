"""Side-by-side page of terrain sandbox runs: one row per window, one column per run label, the chosen map
(relief by default), with Earth atlas images (refs/processed/earth-atlas) beside them for reference.

  uv run python scripts/sandbox/compare_runs.py --name fine-filter --labels fine-nofilter-2,fine-filter-2
      --domains B-core,B-foreland,F-troughs,G-coast [--map relief] [--earth close-everest,close-sognefjord]

Writes out/scratch/sandbox/compare-<name>/index.html; the images are referenced, not copied.
"""
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'
ATLAS = ROOT / 'refs' / 'processed' / 'earth-atlas'

PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Sandbox comparison {name}</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}h2{{margin:24px 0 8px}}
.row{{display:flex;gap:8px;flex-wrap:wrap}}figure{{margin:0;flex:1 1 30%;min-width:280px}}
img{{width:100%;image-rendering:pixelated;cursor:zoom-in}}figcaption{{color:#999;margin:4px 0}}</style></head><body>
<h1>Sandbox comparison: {name}</h1>{sections}</body></html>"""


def figure(path, caption):
    relative = Path('..') / path.relative_to(SANDBOX) if path.is_relative_to(SANDBOX) else path.as_uri()
    source = relative.as_posix() if isinstance(relative, Path) else relative
    return f'<figure><a href="{source}"><img src="{source}"></a><figcaption>{caption}</figcaption></figure>'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--name', required=True)
    parser.add_argument('--labels', required=True)
    parser.add_argument('--domains', required=True)
    parser.add_argument('--map', default='relief')
    parser.add_argument('--earth-map', default=None, help='atlas map name when it differs from --map')
    parser.add_argument('--earth', default='', help='atlas image stems, e.g. close-everest,window-norway')
    arguments = parser.parse_args()

    sections = []
    for domain in arguments.domains.split(','):
        figures = [figure(SANDBOX / label / f'{domain}-{arguments.map}.png', f'{label}: {domain} {arguments.map}')
                   for label in arguments.labels.split(',')]
        sections.append(f'<h2>{domain}</h2><div class="row">{"".join(figures)}</div>')
    earth_map = arguments.earth_map or arguments.map
    if arguments.earth:
        figures = [figure(ATLAS / f'{stem}-{earth_map}.png', f'Earth: {stem} {earth_map}')
                   for stem in arguments.earth.split(',')]
        sections.append(f'<h2>Earth references (refs/processed/earth-atlas)</h2><div class="row">{"".join(figures)}</div>')
    directory = SANDBOX / f'compare-{arguments.name}'
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'index.html').write_text(PAGE.format(name=arguments.name, sections=''.join(sections)), encoding='utf-8')
    print(directory / 'index.html')


if __name__ == '__main__':
    main()
