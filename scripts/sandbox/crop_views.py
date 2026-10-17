"""1:1 crops of a sandbox domain's maps, saved as views that compare_runs.py can show as domains.

For each LABEL and each VIEW `name=x0,y0,x1,y1` (pixels of the domain's PNGs, which are the cropped domain
arrays at 1:1), the maps named by --maps are cropped from `out/scratch/sandbox/LABEL/DOMAIN-<map>.png`,
scaled up by --scale (nearest neighbour, so pixels stay pixels) and written as
`out/scratch/sandbox/LABEL/DOMAIN-<name>-<map>.png`; `compare_runs.py --domains DOMAIN-<name> --map <map>`
then shows them side by side. Used for the straight-river question on the global domain (2048 px wide,
one-pixel rivers), where the whole-domain page shows nothing.

Usage: uv run python scripts/sandbox/crop_views.py --labels A,B --domain global --views basinA=1200,700,1456,956
       [--maps relief,discharge] [--scale 2]
"""
from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--labels", required=True)
    parser.add_argument("--domain", default="global")
    parser.add_argument("--views", required=True, help="name=x0,y0,x1,y1;name2=...")
    parser.add_argument("--maps", default="relief,discharge")
    parser.add_argument("--scale", type=int, default=2)
    args = parser.parse_args()
    views = {}
    for item in args.views.split(";"):
        name, _, box = item.partition("=")
        views[name] = tuple(int(v) for v in box.split(","))
    for label in args.labels.split(","):
        directory = ROOT / "out/scratch/sandbox" / label
        for name, box in views.items():
            for map_name in args.maps.split(","):
                source = directory / f"{args.domain}-{map_name}.png"
                if not source.exists():
                    print("missing", source)
                    continue
                image = Image.open(source).crop(box)
                image = image.resize((image.width * args.scale, image.height * args.scale), Image.NEAREST)
                target = directory / f"{args.domain}-{name}-{map_name}.png"
                image.save(target)
                print(target)


if __name__ == "__main__":
    main()
