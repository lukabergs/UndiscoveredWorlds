"""Side-by-side page of terrain sandbox runs: one row per window, one column per run label, the chosen map
(relief by default), with Earth atlas images (refs/processed/earth-atlas) beside them for reference.

  uv run python scripts/sandbox/compare_runs.py --name fine-filter --labels fine-nofilter-2,fine-filter-2
      --domains B-core,B-foreland,F-troughs,G-coast [--map relief] [--earth close-everest,close-sognefjord]
      [--extra LABEL/DOMAIN,...] [--table metrics.json]

A label is a sandbox run of the main tree (naming.find: an id, <id>-<tag>, an earlier label or an absolute path).
--map takes a comma-separated list (one block of rows per map). --extra adds a row of single images
(<run LABEL>/DOMAIN-<map>.png), e.g. a nested window's parent box once instead of once per label.
--table renders a window_metrics.py JSON (its items as rows) under the images and copies it beside the page.
Writes the main tree's out/evaluation/<id>-compare-<name>/index.html (scripts/uwtools/naming.py, kind page: the
tag compare-<name> follows the tag rule; earlier pages of that tag are then marked superseded, naming.finish); the
images are linked by relative path, not copied.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'scripts' / 'evaluation'))
from uwtools import naming  # noqa: E402
from uwtools.paths import main_root  # noqa: E402
import report  # noqa: E402  (link)

# The main tree's atlas (also from a git worktree, which holds none); runs are found through naming.find.
ATLAS = main_root() / 'refs' / 'processed' / 'earth-atlas'
TAG_PREFIX = 'compare-'

PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Sandbox comparison {name}</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}h2{{margin:24px 0 8px}}
.row{{display:flex;gap:8px;flex-wrap:wrap}}figure{{margin:0;flex:1 1 30%;min-width:280px}}
img{{width:100%;image-rendering:pixelated;cursor:zoom-in}}figcaption{{color:#999;margin:4px 0}}
table{{border-collapse:collapse;font-size:12px}}td,th{{border:1px solid #333;padding:3px 6px;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body>
<h1>Sandbox comparison: {name}</h1>{sections}</body></html>"""


def figure(path, caption, page):
    source = report.link(path, page)
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
    parser.add_argument('--name', required=True,
                        help=f'the page tag is {TAG_PREFIX}<name> (the naming.py tag rule: lower-case words joined by '
                             f'"-", no version suffix, at most {naming.TAG_LENGTH - len(TAG_PREFIX)} characters)')
    parser.add_argument('--labels', required=True, help='sandbox runs: ids, <id>-<tag> names or earlier labels')
    parser.add_argument('--domains', required=True)
    parser.add_argument('--map', default='relief', help='one map name or a comma-separated list')
    parser.add_argument('--earth-map', default=None, help='atlas map name when it differs from --map')
    parser.add_argument('--earth', default='', help='atlas image stems, e.g. close-everest,window-norway')
    parser.add_argument('--extra', default='', help='LABEL/DOMAIN images in a row of their own')
    parser.add_argument('--table', default=None, help='a window_metrics.py JSON to show as a table')
    arguments = parser.parse_args()
    try:
        tag = naming.check_tag(f'{TAG_PREFIX}{arguments.name}')
        # Every run through naming.find, before the page folder exists: a wrong label leaves nothing behind.
        labels = {label: naming.find('sandbox', label) for label in arguments.labels.split(',')}
        extras = [(item, naming.find('sandbox', item.partition('/')[0]) / item.partition('/')[2])
                  for item in arguments.extra.split(',') if item]
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))

    named = naming.new_directory('page', tag)
    directory = named.path
    sections = []
    for map_name in arguments.map.split(','):
        for domain in arguments.domains.split(','):
            figures = [figure(run / f'{domain}-{map_name}.png', f'{label}: {domain} {map_name}', directory)
                       for label, run in labels.items()]
            sections.append(f'<h2>{domain} ({map_name})</h2><div class="row">{"".join(figures)}</div>')
        if extras:
            figures = [figure(Path(f'{stem}-{map_name}.png'), f'{item} {map_name}', directory)
                       for item, stem in extras]
            sections.append(f'<h2>Also ({map_name})</h2><div class="row">{"".join(figures)}</div>')
        earth_map = arguments.earth_map or map_name
        if arguments.earth:
            figures = [figure(ATLAS / f'{stem}-{earth_map}.png', f'Earth: {stem} {earth_map}', directory)
                       for stem in arguments.earth.split(',')]
            sections.append(f'<h2>Earth references (refs/processed/earth-atlas; {earth_map})</h2>'
                            f'<div class="row">{"".join(figures)}</div>')
    if arguments.table:
        sections.append(table(arguments.table))
        if Path(arguments.table).resolve() != (directory / Path(arguments.table).name).resolve():
            shutil.copyfile(arguments.table, directory / Path(arguments.table).name)
    (directory / 'index.html').write_text(PAGE.format(name=arguments.name, sections=''.join(sections)), encoding='utf-8')
    print(directory / 'index.html')
    naming.finish(named, True)  # earlier pages of the tag marked superseded


if __name__ == '__main__':
    main()
