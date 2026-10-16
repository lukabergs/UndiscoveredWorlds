"""Side-by-side page of terrain sandbox runs: one row per window, one column per run label, the chosen map
(relief by default), with Earth atlas images (refs/processed/earth-atlas) beside them for reference.

  uv run python scripts/sandbox/compare_runs.py --name fine-filter --labels fine-nofilter-2,fine-filter-2
      --domains B-core,B-foreland,F-troughs,G-coast [--map relief] [--earth close-everest,close-sognefjord]
      [--extra LABEL/DOMAIN,...] [--table metrics.json]

--map takes a comma-separated list (one block of rows per map). --extra adds a row of single images
(SANDBOX/LABEL/DOMAIN-<map>.png), e.g. a nested window's parent box once instead of once per label.
--table renders a window_metrics.py JSON (its items as rows) under the images and copies it beside the page.
Writes out/scratch/sandbox/compare-<name>/index.html; the images are referenced, not copied.
"""
import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'
ATLAS = ROOT / 'refs' / 'processed' / 'earth-atlas'

PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Sandbox comparison {name}</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}h2{{margin:24px 0 8px}}
.row{{display:flex;gap:8px;flex-wrap:wrap}}figure{{margin:0;flex:1 1 30%;min-width:280px}}
img{{width:100%;image-rendering:pixelated;cursor:zoom-in}}figcaption{{color:#999;margin:4px 0}}
table{{border-collapse:collapse;font-size:12px}}td,th{{border:1px solid #333;padding:3px 6px;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body>
<h1>Sandbox comparison: {name}</h1>{sections}</body></html>"""


def figure(path, caption):
    relative = Path('..') / path.relative_to(SANDBOX) if path.is_relative_to(SANDBOX) else path.as_uri()
    source = relative.as_posix() if isinstance(relative, Path) else relative
    return f'<figure><a href="{source}"><img src="{source}"></a><figcaption>{caption}</figcaption></figure>'


def table(path):
    """An HTML table of a window_metrics.py JSON: one row per item, one column per metric."""
    document = json.loads(Path(path).read_text(encoding='utf-8'))
    items = document['items']
    keys = [key for key in next(iter(items.values())) if key not in ('pixels', 'pixel_km')]
    head = ''.join(f'<th>{key}</th>' for key in ['item'] + keys)
    rows = ''.join('<tr>' + f'<td>{item}</td>' + ''.join(f'<td>{values.get(key, "")}</td>' for key in keys) + '</tr>'
                   for item, values in items.items())
    return (f'<h2>Metrics ({Path(path).name})</h2><p>{document.get("method", "")}</p>'
            f'<table><tr>{head}</tr>{rows}</table>')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--name', required=True)
    parser.add_argument('--labels', required=True)
    parser.add_argument('--domains', required=True)
    parser.add_argument('--map', default='relief', help='one map name or a comma-separated list')
    parser.add_argument('--earth-map', default=None, help='atlas map name when it differs from --map')
    parser.add_argument('--earth', default='', help='atlas image stems, e.g. close-everest,window-norway')
    parser.add_argument('--extra', default='', help='LABEL/DOMAIN images in a row of their own')
    parser.add_argument('--table', default=None, help='a window_metrics.py JSON to show as a table')
    arguments = parser.parse_args()

    directory = SANDBOX / f'compare-{arguments.name}'
    directory.mkdir(parents=True, exist_ok=True)
    sections = []
    for map_name in arguments.map.split(','):
        for domain in arguments.domains.split(','):
            figures = [figure(SANDBOX / label / f'{domain}-{map_name}.png', f'{label}: {domain} {map_name}')
                       for label in arguments.labels.split(',')]
            sections.append(f'<h2>{domain} ({map_name})</h2><div class="row">{"".join(figures)}</div>')
        if arguments.extra:
            figures = [figure(SANDBOX / f'{item}-{map_name}.png', f'{item} {map_name}')
                       for item in arguments.extra.split(',')]
            sections.append(f'<h2>Also ({map_name})</h2><div class="row">{"".join(figures)}</div>')
        earth_map = arguments.earth_map or map_name
        if arguments.earth:
            figures = [figure(ATLAS / f'{stem}-{earth_map}.png', f'Earth: {stem} {earth_map}')
                       for stem in arguments.earth.split(',')]
            sections.append(f'<h2>Earth references (refs/processed/earth-atlas; {earth_map})</h2>'
                            f'<div class="row">{"".join(figures)}</div>')
    if arguments.table:
        sections.append(table(arguments.table))
        if Path(arguments.table).resolve() != (directory / Path(arguments.table).name).resolve():
            shutil.copyfile(arguments.table, directory / Path(arguments.table).name)
    (directory / 'index.html').write_text(PAGE.format(name=arguments.name, sections=''.join(sections)), encoding='utf-8')
    print(directory / 'index.html')


if __name__ == '__main__':
    main()
