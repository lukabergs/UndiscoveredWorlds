"""Download the Copernicus DEM GLO-30 tiles (30 m) of the fjord-coast areas.

Tiles come from the public AWS bucket (copernicus-dem-30m). Each tile gets
<file>.receipt.json beside it (request, bytes, sha256, retrieved_utc, final_url,
content_type, etag); a tile without a receipt is never overwritten. A tile the bucket does
not have (all sea: HTTP 403 or 404) is recorded as absent in <file>.absent.json, not as an
error. Data stay local-only under refs/ (never tracked).

    uv run scripts/refs/fetch_glo30.py                     # every area
    uv run scripts/refs/fetch_glo30.py --area sognefjord   # one area
    uv run scripts/refs/fetch_glo30.py --verify-only       # hashes and receipts only

Copernicus DEM GLO-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018
provided under COPERNICUS by the European Union and ESA; all rights reserved.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
from urllib.parse import urlsplit, urlunsplit

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import glo30
from uwtools.hashing import sha256
from uwtools.protected_io import staged

TILE_MAX_MIB = 120
TIFF_SIGNATURES = (b'II*\x00', b'MM\x00*', b'II+\x00', b'MM\x00+')
ABSENT_STATUS = (403, 404)


def request_for(name):
    return dict(source=glo30.SOURCE, tile=name, url=glo30.tile_url(name), kind='tiff')


def validate(path):
    with path.open('rb') as stream:
        head = stream.read(512)
    if not head:
        raise ValueError('Empty download')
    if b'<html' in head.lower() or b'<!doctype html' in head.lower():
        raise ValueError('Provider returned an HTML page instead of data')
    if not head.startswith(TIFF_SIGNATURES):
        raise ValueError('Invalid TIFF signature')


def _clean(url):
    parts = urlsplit(url)
    return urlunsplit((parts.scheme, parts.netloc, parts.path, '', ''))


def fetch(lat, lon, root, minimum_free, verify_only, client=None):
    """Verify or download one tile; returns ('fetched' | 'verified' | 'absent', bytes)."""
    path = glo30.tile_path(lat, lon, root)
    name = path.stem
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Tile path escapes the output directory')
    receipt = path.with_name(path.name + '.receipt.json')
    absent = glo30.absent_path(path)
    request = request_for(name)
    if path.exists():
        if not receipt.exists():
            raise ValueError('Existing file has no receipt; refusing to overwrite')
        saved = json.loads(receipt.read_text(encoding='utf-8'))
        if saved['request'] != request or saved['sha256'] != sha256(path):
            raise ValueError('Cached tile or request changed')
        validate(path)
        print(f'Verified {path.name}', flush=True)
        return 'verified', path.stat().st_size
    if absent.exists():
        if json.loads(absent.read_text(encoding='utf-8'))['request'] != request:
            raise ValueError('Absent record and request differ')
        print(f'Absent {name} (recorded)', flush=True)
        return 'absent', 0
    if verify_only:
        raise FileNotFoundError('Tile has not been downloaded')
    root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(root).free < minimum_free:
        raise ValueError('Free-space floor reached')
    limit = TILE_MAX_MIB * 1024**2
    size = 0
    client = client or requests.Session()
    with client.get(request['url'], stream=True, timeout=(20, 120)) as response:
        if response.status_code in ABSENT_STATUS:
            record = dict(request=request, http_status=response.status_code,
                          checked_utc=datetime.now(timezone.utc).isoformat(),
                          final_url=_clean(response.url))
            absent.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
            print(f'Absent {name}: HTTP {response.status_code}', flush=True)
            return 'absent', 0
        if response.status_code != 200:
            raise ValueError(f'HTTP {response.status_code}')
        expected = int(response.headers.get('Content-Length', 0))
        if expected > limit:
            raise ValueError(f'Tile exceeds {TILE_MAX_MIB} MiB: {expected} bytes')
        with staged(path) as partial:
            with partial.open('wb') as stream:
                for block in response.iter_content(1024**2):
                    size += len(block)
                    if size > limit or shutil.disk_usage(root).free - len(block) < minimum_free:
                        raise ValueError('Transfer exceeds the tile limit or the free-space floor')
                    stream.write(block)
            if expected and not response.headers.get('Content-Encoding') and size != expected:
                raise ValueError('Incomplete response body')
            validate(partial)
            checksum = sha256(partial)
        result = dict(request=request, bytes=size, sha256=checksum,
                      retrieved_utc=datetime.now(timezone.utc).isoformat(),
                      final_url=_clean(response.url),
                      content_type=response.headers.get('Content-Type'),
                      etag=response.headers.get('ETag'))
    receipt.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(f'Downloaded {path.name}: {size / 1024**2:.2f} MiB', flush=True)
    return 'fetched', size


def area_tiles(names=None):
    """(lat, lon) corners of the tiles meeting the chosen areas (all when names is empty)."""
    unknown = sorted(set(names or ()) - set(glo30.AREAS))
    if unknown:
        raise ValueError(f'Unknown area(s): {", ".join(unknown)}')
    tiles = []
    for name, box in glo30.AREAS.items():
        if names and name not in names:
            continue
        for tile in glo30.tiles_for_box(*box):
            if tile not in tiles:
                tiles.append(tile)
    return tiles


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path, default=glo30.TILE_ROOT)
    parser.add_argument('--area', action='append', help='Area name (repeatable)')
    parser.add_argument('--minimum-free-gib', type=float, default=5)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args(argv)
    tiles = area_tiles(args.area)
    client = requests.Session()
    counts = dict(fetched=0, verified=0, absent=0)
    failed = []
    total = 0
    for lat, lon in tiles:
        try:
            outcome, size = fetch(lat, lon, args.root, int(args.minimum_free_gib * 1024**3),
                                  args.verify_only, client)
            counts[outcome] += 1
            total += size
        except Exception as exc:
            message = str(exc) if isinstance(exc, (ValueError, FileNotFoundError)) else type(exc).__name__
            failed.append(glo30.tile_name(lat, lon))
            print(f'Unavailable {glo30.tile_name(lat, lon)}: {message}', flush=True)
    print(f'{len(tiles)} tiles: {counts["fetched"]} fetched, {counts["verified"]} verified, '
          f'{counts["absent"]} absent, {len(failed)} failed; {total / 1024**2:.1f} MiB under {args.root}')
    return int(bool(failed))


if __name__ == '__main__':
    raise SystemExit(main())
