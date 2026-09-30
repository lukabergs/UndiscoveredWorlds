"""Inspect a fixed native-grid coastal window without rerunning the simulation."""
import argparse
import json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from map_rendering import colourize
from tectonics_report_data import load_bundle

ROOT = Path(__file__).resolve().parents[2]


def review(run):
    record = json.loads((ROOT/f'runs/manifests/tectonics/{run}.json').read_text())
    palette = json.loads((ROOT/'configs/tectonics-palettes.json').read_text())['fields']['elevation']['stops']
    canvas = Image.new('RGB', (1280, 930), '#101b28')
    draw = ImageDraw.Draw(canvas)
    draw.text((15, 10), f'Run {run}: fixed location x=150..190, y=50..90; yellow box is x=160..180, y=60..80. Native cells = 156.25 km.', fill='white')
    results = []
    for panel, time in enumerate([0, 100, 300, 600]):
        snap = next(s for s in record['snapshots'] if s['time_myr'] == time)
        m, f, _ = load_bundle(ROOT/f'runs/checkpoints/tectonics/{run}'/snap['directory'])
        z = f['elevation']
        x0, y0 = (panel % 2)*640, 40+(panel//2)*440
        im = Image.fromarray(colourize(z[50:91, 150:191], palette)).resize((410, 410), Image.Resampling.NEAREST)
        canvas.paste(im, (x0, y0))
        draw.rectangle((x0+100, y0+100, x0+309, y0+309), outline='yellow', width=2)
        draw.line((x0, y0+205, x0+410, y0+205), fill='white')
        draw.text((x0+420, y0), f'{time} Myr; row y=70', fill='white')
        for level in [-6000, 0, 6000]:
            yy=y0+200-level/60
            draw.line((x0+420, yy, x0+635, yy), fill='#566779')
            draw.text((x0+420, yy+2), str(level), fill='#b9c9da')
        points=[(x0+420+k*5.25,y0+200-float(v)/60) for k,v in enumerate(z[70,150:191])]
        draw.line(points, fill='#ffc96b', width=2)
        roi=z[60:81,160:181]
        coast=(roi>=0)&(f['coast_distance_km'][60:81,160:181]<156.25)
        results.append(dict(time_myr=time, mean_m=float(roi.mean()), min_m=float(roi.min()), max_m=float(roi.max()),
            land_fraction=float((roi>=0).mean()), coast_mean_m=float(roi[coast].mean()) if coast.any() else None,
            lowland_cells=int(((roi>=0)&(roi<300)).sum()), shallow_cells=int(((roi<0)&(roi>=-200)).sum()),
            row70={name:np.round(f[name][70,150:191],2).tolist() for name in
                ['elevation','continental_crust_m','oceanic_crust_m','sediment_m','plate_identity']}))
    out=ROOT/f'runs/diagnostics/tectonics/{run}'
    out.mkdir(parents=True,exist_ok=True)
    canvas.save(out/'coastal-window.png')
    (out/'coastal-window.json').write_text(json.dumps(dict(run=run, coordinate_convention='zero-based, top-left; fixed location, not material tracking',states=results),indent=2))
    print(json.dumps([dict((k,v) for k,v in s.items() if k!='row70') for s in results]))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',type=int,nargs='+')
    for run in parser.parse_args().runs:
        review(run)
