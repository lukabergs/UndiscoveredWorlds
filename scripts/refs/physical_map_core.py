"""Shared grid, provenance, palettes and exports for the physical reference atlas."""
from __future__ import annotations

from dataclasses import dataclass, field
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import textwrap
import zlib
import zipfile

import cmocean
import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import rasterio
from rasterio.transform import from_bounds
from rasterio.warp import reproject, Resampling

from cell_grid import conservative_remap, grid_dimensions

REPO = Path(__file__).resolve().parents[2]
NODATA = -9999.9
GLOBAL = (-180., -90., 180., 90.)
MOON_CRS = '+proj=longlat +a=1737400 +b=1737400 +no_defs'
CMO_SOURCE = 'https://matplotlib.org/cmocean/'
WORLD_COVER = {
    10: ('Tree cover', '#006400'), 20: ('Shrubland', '#ffbb22'),
    30: ('Grassland', '#ffff4c'), 40: ('Cropland', '#f096ff'),
    50: ('Built-up', '#fa0000'), 60: ('Bare / sparse vegetation', '#b4b4b4'),
    70: ('Snow and ice', '#f0f0f0'), 80: ('Permanent water', '#0064c8'),
    90: ('Herbaceous wetland', '#0096a0'), 95: ('Mangroves', '#00cf75'),
    100: ('Moss and lichen', '#fae6a0'),
}


def slug(value):
    return re.sub(r'[^a-z0-9]+', '_', str(value).lower()).strip('_')


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


@contextmanager
def atlas_lock(root):
    """OS lock prevents concurrent manifest writers; process exit releases it."""
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    with (root/'.build.lock').open('a+b') as stream:
        if stream.tell()==0:stream.write(b'0');stream.flush()
        stream.seek(0)
        try:
            if os.name=='nt':
                import msvcrt
                msvcrt.locking(stream.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError('Another process is using this atlas output directory') from exc
        yield


def as_float(value):
    return np.asarray(np.ma.filled(np.ma.asarray(value).astype('float32'), np.nan), dtype='float32')


def coverage_bounds(values):
    valid=np.isfinite(values)
    rows=np.flatnonzero(valid.any(axis=1));cols=np.flatnonzero(valid.any(axis=0))
    return [int(cols[0]),int(rows[0]),int(cols[-1]+1),int(rows[-1]+1)] if rows.size else None


def warp(values, width, bounds=GLOBAL, src_crs='EPSG:4326', dst_crs='EPSG:4326', method='average'):
    """GDAL resampling with explicit bounds; categorical values use mode/nearest."""
    height = grid_dimensions(width)[1]
    source = as_float(values)
    result = np.full((height, width), np.nan, 'float32')
    reproject(source, result, src_transform=from_bounds(*bounds, source.shape[1], source.shape[0]),
              src_crs=src_crs, src_nodata=np.nan, dst_transform=from_bounds(*GLOBAL, width, height),
              dst_crs=dst_crs, dst_nodata=np.nan, resampling=getattr(Resampling, method),
              num_threads=2, warp_mem_limit=64)
    return result


def regular_grid(values, lat, lon, width, categorical=False):
    """Preserve partial coverage; drop only a verified duplicate periodic seam."""
    data, lat, lon = as_float(values), np.asarray(lat, float), np.asarray(lon, float)
    if lat[0] < lat[-1]:
        lat, data = lat[::-1], data[::-1]
    if np.isclose(abs(lon[-1] - lon[0]), 360, atol=1e-4):
        if not np.allclose(data[:, 0], data[:, -1], equal_nan=True):
            raise ValueError('Conflicting duplicated longitude seam')
        lon, data = lon[:-1], data[:, :-1]
    dx = np.median(np.diff(lon))
    if not np.allclose(np.diff(lon), dx, atol=1e-4):
        raise ValueError('Longitude coordinates are not regular')
    if np.isclose(dx * len(lon), 360, atol=1e-3) and not categorical:
        return conservative_remap(data, width, source_latitudes=lat, source_longitudes=lon)
    dy = abs(float(np.median(np.diff(lat))))
    bounds = (float(lon[0]-dx/2), max(-90., float(lat[-1]-dy/2)),
              float(lon[-1]+dx/2), min(90., float(lat[0]+dy/2)))
    return warp(data, width, bounds, method='mode' if categorical else 'average')


def aggregate_points(lon, lat, width, values=None, reducer='mean'):
    """One contribution per valid input record; no spatial interpolation."""
    height = grid_dimensions(width)[1]
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    valid = np.isfinite(lon) & np.isfinite(lat) & (abs(lon) <= 180) & (abs(lat) <= 90)
    if values is not None:
        values = np.asarray(values, float)
        valid &= np.isfinite(values)
    lon, lat = lon[valid], lat[valid]
    x = np.floor(((lon + 180) % 360) * width / 360).astype(int)
    y = np.minimum(np.floor((90-lat)*height/180).astype(int), height-1)
    index = y*width+x
    count = np.bincount(index, minlength=height*width).reshape(height,width)
    if values is None or reducer == 'count':
        return count.astype('float32')
    if reducer == 'max':
        output = np.full(height*width, -np.inf)
        np.maximum.at(output, index, values[valid])
        output[~np.isfinite(output)] = np.nan
        return output.reshape(height,width).astype('float32')
    total = np.bincount(index, weights=values[valid], minlength=height*width).reshape(height,width)
    if reducer == 'sum':
        return total.astype('float32')
    return np.divide(total, count, out=np.full_like(total, np.nan), where=count>0).astype('float32')


@dataclass
class Layer:
    key: str
    title: str
    units: str
    dataset: str
    limits: tuple[float,float]
    palette: str = 'viridis'
    scale: str = 'linear'
    classes: dict = field(default_factory=dict)
    notes: str = ''
    method: str = 'Spherical-area weighted mean of valid source cells'
    period: str = 'Provider reference epoch'
    body: str = 'Earth'
    coordinate_frame: str = 'WGS84 geographic'
    palette_status: str = 'Author-selected published scientific palette; no universal official colour standard claimed'
    palette_source: str = CMO_SOURCE
    colours: list | None = None
    bands: list | None = None
    zero_colour: str | None = None

    def __post_init__(self):
        if self.palette_source==CMO_SOURCE and not self.palette.startswith('cmo.'):
            self.palette_source=('manifest.json' if self.classes else
                                 'https://matplotlib.org/stable/users/explain/colors/colormaps.html')


def provider_cpt(name):
    bands=[]
    for line in (REPO/'assets/palettes/physical'/name).read_text().splitlines():
        parts=line.split(';')[0].replace('/',' ').split()
        if len(parts)!=8 or parts[0].startswith(('#','B','F','N')):continue
        numbers=list(map(float,parts))
        bands.append([numbers[0],numbers[4],numbers[1:4],numbers[5:8]])
    if not bands: raise ValueError('Empty provider CPT')
    return bands


def band_rgb(values,bands):
    output=np.zeros((*values.shape,3),dtype=float)
    output[:]=bands[0][2]
    for lo,hi,a,b in bands:
        selected=values>=lo
        ratio=np.clip((values-lo)/(hi-lo),0,1)[...,None]
        colours=np.asarray(a)+(np.asarray(b)-a)*ratio
        output[selected]=colours[selected]
    return np.rint(output).astype('uint8')


def colour_table(layer):
    if layer.bands:
        return band_rgb(np.linspace(*layer.limits,256),layer.bands)
    if layer.colours is not None:
        return np.asarray(layer.colours, dtype='uint8')
    name = layer.palette
    cmap = getattr(cmocean.cm, name[4:]) if name.startswith('cmo.') else matplotlib.colormaps[name]
    return np.rint(cmap(np.linspace(0,1,256))[:,:3]*255).astype('uint8')


def rgba(values, layer):
    valid = np.isfinite(values)
    output = np.zeros((*values.shape,4), 'uint8')
    if layer.classes:
        for code, (_, colour) in layer.classes.items():
            mask = valid & (values == int(code))
            output[mask,:3] = tuple(bytes.fromhex(colour.lstrip('#')))
        valid &= np.isin(values, list(map(int,layer.classes)))
    elif layer.bands:
        output[:,:,:3]=band_rgb(np.where(valid,values,layer.limits[0]),layer.bands)
    else:
        low, high = layer.limits
        safe = np.where(valid, values, low)
        if layer.scale == 'log':
            ratio = (np.log10(np.maximum(safe,low))-np.log10(low))/(np.log10(high)-np.log10(low))
        else:
            ratio = (safe-low)/(high-low)
        index = np.rint(np.clip(ratio,0,1)*255).astype('uint8')
        output[:,:,:3] = colour_table(layer)[index]
    output[:,:,3] = valid*255
    output[~valid,:3] = 0
    if layer.zero_colour:
        output[valid & (values==0),:3]=tuple(bytes.fromhex(layer.zero_colour.lstrip('#')))
    return output


def font(size=15):
    try:
        return ImageFont.truetype('C:/Windows/Fonts/segoeui.ttf', size)
    except OSError:
        return ImageFont.load_default(size=size)


def legend(path, layer):
    width = 780
    class_lines={key:textwrap.wrap(f'{key}: {name}',width=88) or [''] for key,(name,_) in layer.classes.items()}
    height = 100 + (sum(24*len(lines)+4 for lines in class_lines.values()) if layer.classes else 80)
    canvas = Image.new('RGB', (width,height), '#111820')
    draw = ImageDraw.Draw(canvas)
    draw.text((18,12), layer.title, font=font(20), fill='white')
    draw.text((18,43), layer.units + ' | ' + layer.period, font=font(14), fill='#c7d3df')
    if layer.classes:
        y=75
        for key,(name,colour) in layer.classes.items():
            draw.rectangle((18,y,42,y+18),fill=colour)
            for line in class_lines[key]:
                draw.text((54,y),line,font=font(),fill='white')
                y+=24
            y+=4
    else:
        strip = np.repeat(colour_table(layer)[None,:,:], 24, axis=0)
        canvas.paste(Image.fromarray(strip).resize((720,24)),(24,78))
        for i in range(5):
            a,b = layer.limits
            value = 10**(np.log10(a)+(np.log10(b)-np.log10(a))*i/4) if layer.scale=='log' else a+(b-a)*i/4
            draw.text((24+i*180,108),f'{value:g}',font=font(13),fill='white',anchor='mt')
        extra=(' Zero = '+('white' if layer.zero_colour=='#ffffff' else 'charcoal')+'.') if layer.zero_colour else ''
        draw.text((24,140),f'{layer.scale} scale; GeoTIFF values are unclamped.'+extra,font=font(12),fill='#c7d3df')
    path.parent.mkdir(parents=True,exist_ok=True)
    canvas.save(path)


class Atlas:
    def __init__(self, source, output, widths):
        self.source, self.output, self.widths = Path(source), Path(output), sorted(set(widths))
        self.products, self.used, self.status = [], {}, []
        self.output.mkdir(parents=True,exist_ok=True)
        self.cache = self.output / 'cache'
        self.cache.mkdir(exist_ok=True)
        self.catalog = json.loads((REPO/'configs/physical-reference-sources.json').read_text())
        self.entries = {(x['dataset'],x['file']):x for x in self.catalog['files']}
        previous=self.output/'manifest.json'
        if previous.exists():
            old=json.loads(previous.read_text(encoding='utf-8'))
            self.products=old.get('products',[])
            self.used=old.get('sources',{})
            self.status=old.get('datasets',[])
        self.checked=set()
        self.checked_cache=set()

    def paths(self, dataset, pattern='*'):
        paths = [p for p in sorted((self.source/dataset).glob(pattern))
                 if p.is_file() and not p.name.endswith(('.receipt.json','.part'))]
        for path in paths:
            key = f'{dataset}/{path.name}'
            if key not in self.checked:
                receipt = json.loads(path.with_name(path.name+'.receipt.json').read_text())
                if receipt['sha256'] != sha256(path) or receipt['request']!=self.entries[(dataset,path.name)]:
                    raise ValueError('Source checksum mismatch: '+key)
                self.used[key] = dict(sha256=receipt['sha256'],bytes=path.stat().st_size,
                                     provenance=self.entries[(dataset,path.name)])
                self.checked.add(key)
        return paths

    def one(self, dataset, pattern):
        paths = self.paths(dataset,pattern)
        if len(paths)!=1:
            raise ValueError(f'Expected one {dataset}/{pattern}, found {len(paths)}')
        return paths[0]

    def unzip(self, path, member):
        destination = self.cache / path.stem / member
        if not destination.resolve().is_relative_to(self.cache.resolve()):
            raise ValueError('Archive path escapes cache')
        with zipfile.ZipFile(path) as archive:
            info = archive.getinfo(member)
            if destination.exists() and destination not in self.checked_cache:
                crc=0
                with destination.open('rb') as stream:
                    while block:=stream.read(1024**2):crc=zlib.crc32(block,crc)
                if destination.stat().st_size!=info.file_size or crc!=info.CRC:
                    raise ValueError('Extracted source cache CRC mismatch: '+str(destination))
            if not destination.exists() or destination.stat().st_size != info.file_size:
                if info.file_size > 3*1024**3 or shutil.disk_usage(self.output).free-info.file_size < 20*1024**3:
                    raise ValueError('Extraction exceeds storage safety limits')
                destination.parent.mkdir(parents=True,exist_ok=True)
                with archive.open(member) as source, destination.open('wb') as output:
                    shutil.copyfileobj(source,output,1024**2)
            self.checked_cache.add(destination)
        return destination

    def emit(self, layer, make):
        print('MAP '+layer.key,flush=True)
        records = []
        legend_path = self.output/'legends'/f'{layer.key}.png'
        legend(legend_path,layer)
        for width in self.widths:
            if shutil.disk_usage(self.output).free < 20*1024**3:
                raise ValueError('Destination free-space floor reached')
            values = as_float(make(width))
            if values.shape != (width//2,width):
                raise ValueError(f'Unexpected output shape for {layer.key}')
            png = self.output/'maps'/layer.key/f'{width}.png'
            tif = self.output/'fields'/layer.key/f'{width}.tif'
            png.parent.mkdir(parents=True,exist_ok=True)
            tif.parent.mkdir(parents=True,exist_ok=True)
            Image.fromarray(rgba(values,layer)).save(png)
            crs = 'EPSG:4326' if layer.body=='Earth' else MOON_CRS
            with rasterio.open(tif,'w',driver='GTiff',width=width,height=width//2,count=1,
                               dtype='float32',crs=crs,transform=from_bounds(*GLOBAL,width,width//2),
                               nodata=NODATA,compress='deflate',predictor=3,tiled=True) as out:
                out.write(np.where(np.isfinite(values),values,NODATA),1)
                out.set_band_description(1,layer.title)
                out.update_tags(AREA_OR_POINT='Area',units=layer.units,body=layer.body,
                                coordinate_frame=layer.coordinate_frame,
                                source_dataset=layer.dataset,period=layer.period,method=layer.method)
            # PNG coordinates belong to the same body-specific grid, not automatically WGS84.
            png.with_suffix('.pgw').write_text(f'{360/width}\n0\n0\n{-360/width}\n{-180+180/width}\n{90-180/width}\n')
            png.with_suffix('.prj').write_text(rasterio.crs.CRS.from_user_input(crs).to_wkt())
            valid = values[np.isfinite(values)]
            records.append(dict(width=width,height=width//2,valid_cells=int(valid.size),
                                coverage_pixel_bounds=coverage_bounds(values),
                                min=float(valid.min()) if valid.size else None,
                                max=float(valid.max()) if valid.size else None,
                                png=str(png.relative_to(self.output)).replace('\\','/'),
                                tiff=str(tif.relative_to(self.output)).replace('\\','/'),
                                png_sha256=sha256(png),tiff_sha256=sha256(tif)))
        entry = dict(vars(layer),legend=str(legend_path.relative_to(self.output)).replace('\\','/'),outputs=records)
        self.products=[p for p in self.products if p['key']!=layer.key]+[entry]
        self.save()

    def save(self):
        manifest = dict(version=1,grid='Cell centres, 2:1, longitude [-180,180), north-up',
                        widths=self.widths,nodata=NODATA,png_nodata='Transparent black',
                        products=self.products,sources=self.used,datasets=self.status)
        path = self.output/'manifest.json'
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        temporary.replace(path)
