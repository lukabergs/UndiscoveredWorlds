# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.2", "pillow==12.3.0", "rasterio==1.5.1", "netcdf4==1.7.4", "cmocean==4.0.3", "matplotlib==3.11.1", "pyhdf==0.11.7", "openpyxl==3.1.5", "pyshp==3.1.6", "geopandas==1.1.4", "pyshtools==4.14.1", "pygplates==1.0.0"]
# ///
"""Prepare physical reference maps at 128–2048 px from retained, verified sources."""
import argparse
import json
import sys
import traceback

from physical_map_core import Atlas, REPO, atlas_lock
import physical_map_rasters as rasters
import physical_map_points as points
import physical_map_vectors as vectors
import physical_map_fields as fields

FAMILIES={name:getattr(rasters,name) for name in ('glim','crust','seafloor','emit_minerals','soils','ocean','regional_rasters','surface_water','moon_topography','productivity')}
FAMILIES.update({name:getattr(points,name) for name in ('heatflow','volcanoes','deposits','geochemistry','discharge')})
FAMILIES.update({name:getattr(vectors,name) for name in ('geology_uk','groundwater','basins','rivers','lakes','glaciers')})
FAMILIES.update({name:getattr(fields,name) for name in ('gravity','magnetism','plates')})


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',default=REPO/'refs/source/physical')
    parser.add_argument('--output',default=REPO/'refs/processed/physical')
    parser.add_argument('--width',type=int,nargs='+',default=[128,256,512,1024,2048])
    parser.add_argument('--family',nargs='+',choices=list(FAMILIES))
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--verify-only',action='store_true')
    parser.add_argument('--index-only',action='store_true')
    parser.add_argument('--render-only',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        import unittest
        suite=unittest.defaultTestLoader.discover(str(REPO/'tests/scripts'),pattern='test_physical_reference_maps.py')
        return int(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    if any(w not in (128,256,512,1024,2048) for w in args.width):
        parser.error('Use the standard widths: 128 256 512 1024 2048')
    with atlas_lock(args.output):
        return run(args)


def run(args):
    from physical_map_index import index, verify, render_existing
    if args.verify_only:return verify(args.output,args.width)
    if args.render_only:
        render_existing(args.output)
        index(args.output)
        return 0
    if args.index_only:
        index(args.output)
        return 0
    atlas=Atlas(args.source,args.output,args.width)
    failed=[]
    for family in args.family or FAMILIES:
        print('FAMILY '+family,flush=True)
        status=dict(family=family,status='complete')
        try:
            FAMILIES[family](atlas)
        except Exception as exc:
            status.update(status='failed',error=str(exc))
            failed.append(family)
            traceback.print_exc()
        atlas.status=[x for x in atlas.status if x['family']!=family]+[status]
        atlas.save()
    print(json.dumps(dict(products=len(atlas.products),failed=failed)),flush=True)
    index(args.output)
    return int(bool(failed))


if __name__=='__main__':
    raise SystemExit(main())
