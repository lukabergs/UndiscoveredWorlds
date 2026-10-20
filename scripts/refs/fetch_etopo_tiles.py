"""Download the ETOPO 2022 15-arc-second surface tiles the Earth-window catalogue needs.

Only tiles whose 15 degree square meets a catalogue box are fetched (not the globe). Each
tile gets <file>.receipt.json beside it (request, bytes, sha256, retrieved_utc, final_url,
content_type, etag); a tile without a receipt is never overwritten. Data stay local-only
under refs/ (never tracked).

    uv run scripts/refs/fetch_etopo_tiles.py                  # every catalogue window
    uv run scripts/refs/fetch_etopo_tiles.py --window zagros  # one window
    uv run scripts/refs/fetch_etopo_tiles.py --verify-only    # hashes and receipts only
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
from uwtools import etopo
from uwtools.hashing import sha256
from uwtools.protected_io import staged

SOURCE = 'ETOPO 2022 v1 15 arc-second surface elevation (DOI 10.25921/fd45-gt74)'
TILE_MAX_MIB = 64
TIFF_SIGNATURES = (b'II*\x00', b'MM\x00*', b'II+\x00', b'MM\x00+')


def request_for(name):
    return dict(source=SOURCE, tile=name, url=etopo.tile_url(name), kind='tiff')


def validate(path):
    with path.open('rb') as stream:
        head = stream.read(512)
    if not head:
        raise ValueError('Empty download')
    if b'<html' in head.lower() or b'<!doctype html' in head.lower():
        raise ValueError('Provider returned an HTML page instead of data')
    if not head.startswith(TIFF_SIGNATURES):
        raise ValueError('Invalid TIFF signature')


def stored_bytes(root):
    return sum(p.stat().st_size for p in root.glob('*.tif') if p.is_file())


def fetch(name, root, budget, minimum_free, verify_only, client=None):
    """Verify or download one tile; returns its size in bytes."""
    path = etopo.tile_path(name, root)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Tile path escapes the output directory')
    receipt = path.with_name(path.name + '.receipt.json')
    request = request_for(name)
    if path.exists():
        if not receipt.exists():
            raise ValueError('Existing file has no receipt; refusing to overwrite')
        saved = json.loads(receipt.read_text(encoding='utf-8'))
        if saved['request'] != request or saved['sha256'] != sha256(path):
            raise ValueError('Cached tile or request changed')
        validate(path)
        print(f'Verified {path.name}', flush=True)
        return path.stat().st_size
    if verify_only:
        raise FileNotFoundError('Tile has not been downloaded')
    root.mkdir(parents=True, exist_ok=True)
    limit = min(TILE_MAX_MIB * 1024**2, budget - stored_bytes(root))
    if limit <= 0 or shutil.disk_usage(root).free < minimum_free:
        raise ValueError('Size budget or free-space floor reached')
    size = 0
    client = client or requests.Session()
    with client.get(request['url'], stream=True, timeout=(20, 120)) as response:
        if response.status_code != 200:
            raise ValueError(f'HTTP {response.status_code}')
        expected = int(response.headers.get('Content-Length', 0))
        if expected > limit:
            raise ValueError(f'Tile exceeds the remaining budget: {expected} > {limit}')
        with staged(path) as partial:
            with partial.open('wb') as stream:
                for block in response.iter_content(1024**2):
                    size += len(block)
                    if size > limit or shutil.disk_usage(root).free - len(block) < minimum_free:
                        raise ValueError('Transfer exceeds the size budget')
                    stream.write(block)
            if expected and not response.headers.get('Content-Encoding') and size != expected:
                raise ValueError('Incomplete response body')
            validate(partial)
            checksum = sha256(partial)
        url = urlsplit(response.url)
        result = dict(request=request, bytes=size, sha256=checksum,
                      retrieved_utc=datetime.now(timezone.utc).isoformat(),
                      final_url=urlunsplit((url.scheme, url.netloc, url.path, '', '')),
                      content_type=response.headers.get('Content-Type'),
                      etag=response.headers.get('ETag'))
    receipt.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(f'Downloaded {path.name}: {size / 1024**2:.2f} MiB', flush=True)
    return size


def catalogue_tiles(windows, names=None):
    """Tiles meeting the boxes of the chosen windows (all when names is empty), in order."""
    known = {window['name'] for window in windows}
    unknown = sorted(set(names or ()) - known)
    if unknown:
        raise ValueError(f'Unknown window(s): {", ".join(unknown)}')
    tiles = []
    for window in windows:
        if names and window['name'] not in names:
            continue
        for tile in etopo.tiles_for_box(*window['box']):
            if tile not in tiles:
                tiles.append(tile)
    return tiles


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--catalogue', type=Path, default=etopo.CATALOGUE)
    parser.add_argument('--root', type=Path, default=etopo.TILE_ROOT)
    parser.add_argument('--window', action='append', help='Catalogue window name (repeatable)')
    parser.add_argument('--max-mib', type=float, default=1500,
                        help='Size budget of all tiles under --root')
    parser.add_argument('--minimum-free-gib', type=float, default=5)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args(argv)
    tiles = catalogue_tiles(etopo.load_catalogue(args.catalogue), args.window)
    client = requests.Session()
    failed = []
    total = 0
    for name in tiles:
        try:
            total += fetch(name, args.root, int(args.max_mib * 1024**2),
                           int(args.minimum_free_gib * 1024**3), args.verify_only, client)
        except Exception as exc:
            message = str(exc) if isinstance(exc, (ValueError, FileNotFoundError)) else type(exc).__name__
            failed.append(name)
            print(f'Unavailable {name}: {message}', flush=True)
    verified = len(tiles) - len(failed)
    print(f'{verified} of {len(tiles)} tiles verified, {total / 1024**2:.1f} MiB under {args.root}')
    return int(bool(failed))


if __name__ == '__main__':
    raise SystemExit(main())
