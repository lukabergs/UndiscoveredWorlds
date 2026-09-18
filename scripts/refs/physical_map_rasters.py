"""Raster and crust-model adapters; numeric values remain separate from colour."""
import io
import re
import tarfile
import zipfile

import netCDF4
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.warp import reproject, Resampling
from pyhdf.SD import SD, SDC

from physical_map_core import Layer, WORLD_COVER, as_float, regular_grid, warp, GLOBAL, provider_cpt


def glim(atlas):
    path=atlas.one('glim-2012','*.zip')
    with zipfile.ZipFile(path) as z:
        text=z.read('glim_wgs84_0point5deg.txt.asc').decode()
    lines=text.splitlines()
    header={s.split()[0].lower():float(s.split()[1]) for s in lines[:6]}
    data=np.loadtxt(io.StringIO('\n'.join(lines[6:])),dtype='float32')
    data[data==header['nodata_value']]=np.nan
    # Class codes follow the provider. The colours below are explicitly author-selected.
    labels=['Unconsolidated sediment','Basic volcanic rock','Siliciclastic sedimentary rock',
            'Basic plutonic rock','Mixed sedimentary rock','Carbonate sedimentary rock',
            'Acid volcanic rock','Metamorphic rock','Acid plutonic rock','Intermediate volcanic rock',
            'Water body','Pyroclastic rock','Intermediate plutonic rock','Evaporite','No data','Ice / glacier']
    colours=['eed9a3','774b8b','d7aa72','5a3475','bbb977','76b8c9','ef8262','bc82bc',
             'ed5650','b16699','397ac4','dbaa9b','e39698','ddb5dc','555555','eaf4ff']
    data[data==15]=np.nan
    classes={i+1:(label,'#'+colour) for i,(label,colour) in enumerate(zip(labels,colours)) if i!=14}
    bounds=(header['xllcorner'],header['yllcorner'],header['xllcorner']+header['cellsize']*data.shape[1],
            header['yllcorner']+header['cellsize']*data.shape[0])
    layer=Layer('geology/lithology','Dominant surface lithology','GLiM category','glim-2012',(1,16),
                classes=classes,method='GDAL mode for categories; no averaging of class IDs',
                notes='Provider class codes; selected rock-family colours are not an official GLiM palette.',
                palette_status='Author-selected categorical geology colours',
                palette_source='https://doi.org/10.1029/2012GC004370')
    atlas.emit(layer,lambda w:warp(data,w,bounds,method='mode'))


def crust(atlas):
    path=atlas.one('crust1-2013','crust1.0.tar.gz')
    with tarfile.open(path) as archive:
        arrays={name:np.loadtxt(archive.extractfile('crust1.'+name)).reshape(180,360,9) for name in ('bnds','vp','vs','rho')}
    lat=89.5-np.arange(180);lon=-179.5+np.arange(360)
    b=arrays['bnds']
    specs=[('moho_depth',-b[:,:,8],'Moho depth below sea level','km',(0,80),'cmo.deep'),
           ('solid_crust_thickness',b[:,:,2]-b[:,:,8],'Solid crust thickness (sediments + crystalline crust)','km',(0,80),'cmo.deep'),
           ('sediment_thickness',b[:,:,2]-b[:,:,5],'Crustal sediment thickness','km',(0,15),'cmo.turbid')]
    for i,name in enumerate(('upper','middle','lower'),5):
        specs.append((name+'_crystalline_thickness',b[:,:,i]-b[:,:,i+1],name.title()+' crystalline crust thickness','km',(0,35),'cmo.deep'))
    for param,units,limits,palette in [('rho','kg m-3',(2400,3400),'cmo.dense'),('vp','km s-1',(3,9),'viridis'),('vs','km s-1',(1,5),'viridis')]:
        for i,name in [(5,'upper_crust'),(6,'middle_crust'),(7,'lower_crust'),(8,'uppermost_mantle')]:
            values=arrays[param][:,:,i]*(1000 if param=='rho' else 1)
            specs.append((name+'_'+param,values,name.replace('_',' ').title()+' '+{'rho':'density','vp':'P-wave speed','vs':'S-wave speed'}[param],units,limits,palette))
    for key,data,title,units,limits,palette in specs:
        layer=Layer('crust/'+key,title,units,'crust1-2013',limits,palette,
                    notes='CRUST1.0 1-degree layered model. Solid thickness excludes water and ice; Moho depth uses sea level.',period='2013 model')
        atlas.emit(layer,lambda w,a=data:regular_grid(a,lat,lon,w))


def seafloor(atlas):
    names={'age':('Ocean-floor age','Ma',(0,200),'cmo.thermal_r'),
           'full_rate':('Full seafloor spreading rate','mm yr-1',(0,180),'cmo.speed'),
           'asym':('Spreading asymmetry','percent',(0,100),'cmo.balance'),
           'conf':('Ocean-floor age confidence rating','provider rating',(0,5),'viridis'),
           'dir':('Spreading direction','degrees',(0,360),'cmo.phase'),
           'obliq':('Spreading obliquity','degrees',(0,90),'viridis')}
    for path in atlas.paths('earthbyte-age-2020','*.nc'):
        key=path.name.split('.')[0]
        with netCDF4.Dataset(path) as ds:
            data=as_float(ds['z'][:]);lat=ds['lat'][:];lon=ds['lon'][:]
            attrs={k:str(v) for k,v in ds['z'].__dict__.items()}
        title,units,limits,palette=names[key]
        if key=='dir':
            # Axial orientations wrap at 180 degrees; average doubled angles.
            def remap(w):
                theta=2*np.deg2rad(data)
                return (np.rad2deg(np.arctan2(regular_grid(np.sin(theta),lat,lon,w),regular_grid(np.cos(theta),lat,lon,w)))/2)%180
        elif key=='conf':
            remap=lambda w:regular_grid(data,lat,lon,w,categorical=True)
        else:
            remap=lambda w:regular_grid(data,lat,lon,w)
        cpt={'age':'age_2020.cpt','asym':'asymgrid.cpt','conf':'confidence_rating.cpt',
             'dir':'spreaddirgrid.cpt','full_rate':'fullrategrid.cpt','obliq':'spreadobliqgrid.cpt'}[key]
        bands=provider_cpt(cpt)
        limits=(bands[0][0],bands[-1][1])
        classes={0:('Higher confidence: code 0','#ffffff'),1:('Higher confidence: code 1','#ffffff'),2:('Lower confidence: code 2','#e41a1c')} if key=='conf' else {}
        atlas.emit(Layer('tectonics/seafloor_'+key,title,units,'earthbyte-age-2020',limits,palette,bands=bands,classes=classes,
                        period='2020 model',notes='Source variable metadata: '+str(attrs),
                        palette_status='Official EarthByte published CPT colour table',
                        palette_source='https://www.earthbyte.org/webdav/ftp/earthbyte/agegrid/2020/Grids/cpt/'+cpt,
                        method='Axial circular mean (180-degree periodic orientation)' if key=='dir' else 'Mode of provider confidence codes' if key=='conf' else 'Spherical-area weighted mean'),remap)
    path=atlas.one('globsed-v3','*.zip')
    with zipfile.ZipFile(path) as z: member=next(n for n in z.namelist() if n.endswith('.nc'))
    p=atlas.unzip(path,member)
    with netCDF4.Dataset(p) as ds:
        lon=ds['lon'][:];lat=ds['lat'][:];data=as_float(ds['z'][:])
        units=getattr(ds['z'],'units','m')
    atlas.emit(Layer('geology/marine_sediment_thickness','Marine sediment thickness',units,'globsed-v3',(0,10000),'cmo.turbid',period='GlobSed v3'),
               lambda w:regular_grid(data,lat,lon,w))


def emit_minerals(atlas):
    path=atlas.one('emit-l3-v002','*.nc')
    with netCDF4.Dataset(path) as ds:
        lat=ds['lat'][:];lon=ds['lon'][:]
        display={}
        for base,var in ds.variables.items():
            if var.ndim==2 and '_UQ_' not in base:
                q=float(np.nanquantile(as_float(var[:]),.99))
                step=10**np.floor(np.log10(max(q,1e-6)))
                display[base]=min(1.,float(np.ceil(q/step)*step))
        for name,var in ds.variables.items():
            if var.ndim!=2: continue
            data=as_float(var[:]);key=name.lower().replace('+','_').replace(' ','_')
            base=name.split('_UQ_')[0]
            layer=Layer('mineralogy/'+key,name.replace('_',' ')+' spectral abundance','fraction','emit-l3-v002',(0,display[base]),'cmo.matter',
                        period='2022–2024 observations / v002',
                        notes='Exposed-surface spectral abundance with provider vegetation/water/cloud masks. UQ fields reflect grain-size assumptions; not underground ore fraction. Display maximum is the rounded 99th percentile of valid baseline values, shared with its UQ maps; numeric values are never clipped.',
                        method='GDAL area resampling over native 0.5-degree observed footprint; no extrapolation')
            atlas.emit(layer,lambda w,a=data:regular_grid(a,lat,lon,w))


def soils(atlas):
    specs={'clay':('Clay content','percent',.1,(0,70),'cmo.turbid'),
           'silt':('Silt content','percent',.1,(0,70),'cmo.turbid'),
           'sand':('Sand content','percent',.1,(0,100),'cmo.turbid'),
           'phh2o':('Soil pH in water','pH',.1,(3,10),'cmo.balance'),
           'soc':('Soil organic carbon','g kg-1',.1,(0,200),'cmo.matter'),
           'nitrogen':('Total soil nitrogen','g kg-1',.01,(0,15),'cmo.algae'),
           'bdod':('Soil bulk density','kg m-3',10,(0,1800),'cmo.dense'),
           'cec':('Cation exchange capacity','cmol(c) kg-1',.1,(0,60),'viridis')}
    for p in atlas.paths('soilgrids-2-wcs','*.tif'):
        key=p.name.split('_')[0];title,units,factor,limits,palette=specs[key]
        with rasterio.open(p) as ds:
            data=as_float(ds.read(1,masked=True))*factor
            lat=ds.transform.f+(np.arange(ds.height)+.5)*ds.transform.e
            lon=ds.transform.c+(np.arange(ds.width)+.5)*ds.transform.a
        atlas.emit(Layer('soils/'+key,title+' (0–5 cm)',units,'soilgrids-2-wcs',limits,palette,
                        notes=f'Native encoded values multiplied by {factor:g}. Input is provider-resampled 0.25 degrees, 60S–85N; deeper soil is not represented.',
                        palette_source='https://docs.isric.org/globaldata/soilgrids/SoilGrids_faqs_01.html'),
                   lambda w,a=data:regular_grid(a,lat,lon,w))


def ocean(atlas):
    specs={'o':('Dissolved oxygen',(0,450),'cmo.ice_r'), 'n':('Nitrate',(0,50),'cmo.haline'),
           'p':('Phosphate',(0,4),'cmo.haline'),'i':('Silicate',(0,180),'cmo.haline')}
    for p in atlas.paths('woa23-biogeochemistry','*.nc'):
        key=p.name.split('_')[2][0]; title,limits,palette=specs[key]
        with netCDF4.Dataset(p) as ds:
            lat=ds['lat'][:];lon=ds['lon'][:];depths=ds['depth'][:];v=ds[key+'_an']
            for depth in (0,100,1000,3000):
                index=int(np.argmin(abs(depths-depth)))
                if not np.isclose(depths[index],depth): raise ValueError('Requested WOA depth not present')
                data=as_float(v[0,index,:,:])
                atlas.emit(Layer(f'ocean/{key}/{depth}m',f'{title} at {depth} m',v.units,'woa23-biogeochemistry',limits,palette,
                                period='WOA23 annual climatology; 1965–2022',notes='Objectively analysed mean; land, seafloor and unsupported cells retain nodata.'),
                           lambda w,a=data:regular_grid(a,lat,lon,w))
    p=atlas.one('modis-chlorophyll-2020','*.nc')
    with netCDF4.Dataset(p) as ds:
        v=ds['chlor_a'];data=as_float(v[:]);lat=ds['lat'][:];lon=ds['lon'][:]
        palette=np.asarray(ds['palette'][:])
        if palette.shape==(3,256):palette=palette.T
        layer=Layer('ocean/chlorophyll','Ocean chlorophyll-a',v.units,'modis-chlorophyll-2020',
                    (float(v.display_min),float(v.display_max)),scale=str(v.display_scale),
                    colours=palette.tolist(),period='2020 annual composite',
                    palette_status='Official provider palette and display limits embedded in source NetCDF',
                    palette_source='https://oceancolor.gsfc.nasa.gov/data/aqua/',notes='Annual composite, not a multiyear climatology.')
    atlas.emit(layer,lambda w:regular_grid(data,lat,lon,w))


def regional_rasters(atlas):
    for p in atlas.paths('worldcover-2021-samples','*.tif'):
        tile=p.name.split('_')[-2]
        def make(w):
            with rasterio.open(p) as ds:
                result=np.full((w//2,w),np.nan,'float32')
                reproject(rasterio.band(ds,1),result,src_nodata=0,dst_nodata=np.nan,
                          dst_transform=from_bounds(*GLOBAL,w,w//2),dst_crs='EPSG:4326',
                          resampling=Resampling.mode,warp_mem_limit=64,num_threads=2)
            return result
        atlas.emit(Layer('land_cover/'+tile,'ESA WorldCover '+tile,'land-cover class','worldcover-2021-samples',(10,100),
                        classes=WORLD_COVER,period='2021',method='Mode of source classes intersecting target cells',
                        palette_status='Official ESA WorldCover class colours',palette_source='https://esa-worldcover.org/en/data-access',
                        notes='Regional 3x3-degree tile only. Outside this tile is nodata.'),make)
    p=atlas.one('gbif-coverage-2026-07','*.tiff')
    with rasterio.open(p) as ds:
        data=as_float(ds.read(1,masked=True));bounds=tuple(ds.bounds)
    atlas.emit(Layer('biosphere/gbif_record_density','GBIF observation density','records per source 0.1-degree cell','gbif-coverage-2026-07',(1,10000),'cmo.amp',scale='log',
                    period='2026-07-01 snapshot',method='Mean source-cell observation counts; not total records per output cell',
                    notes='Sampling-effort reference only; neither species richness nor abundance.'),lambda w:warp(data,w,bounds))


def surface_water(atlas):
    p=atlas.one('jrc-gsw-2024-sample','*.tif')
    def make_water(w):
        # Provider 255 encodes no-data even though this TIFF omits a nodata tag.
        with rasterio.open(p) as ds:
            result=np.full((w//2,w),np.nan,'float32')
            reproject(rasterio.band(ds,1),result,src_nodata=255,dst_nodata=np.nan,
                      dst_transform=from_bounds(*GLOBAL,w,w//2),dst_crs='EPSG:4326',
                      resampling=Resampling.average,warp_mem_limit=64,num_threads=2)
        return result
    # JRC defines red-to-blue with opacity proportional to occurrence. Composite
    # on white so valid low occurrence is distinguishable from transparent nodata.
    value=np.linspace(0,100,256);ratio=np.clip((value-1)/99,0,1)
    base=np.stack([255*(1-ratio),np.zeros(256),255*ratio],axis=1)
    colours=np.rint(255+(base-255)*(value/100)[:,None]).astype('uint8').tolist()
    atlas.emit(Layer('hydrology/surface_water_occurrence','Surface-water occurrence, western Europe','percent','jrc-gsw-2024-sample',(0,100),
                    colours=colours,period='1984–2024',method='GDAL average of valid occurrence percentages',
                    palette_status='Official JRC colour/opacity convention composited on white; nodata transparent',
                    palette_source='https://storage.googleapis.com/water-world/downloads_ancillary/DataUsersGuidev2024_v.5.pdf',
                    zero_colour='#ffffff',notes='One regional tile, not global coverage. 0 means observed absence of water, 255 means no data. Provider red-to-blue and 1–100% opacity convention flattened on white for valid pixels; missing data stays transparent.'),make_water)


def moon_topography(atlas):
    path=atlas.one('lola-ldem4','*.img')
    atlas.one('lola-ldem4','*.lbl')
    data=np.fromfile(path,dtype='<i2').reshape(720,1440).astype('float32')*.5
    data=np.roll(data,720,axis=1)
    atlas.emit(Layer('moon/elevation','Lunar elevation relative to 1737.4 km sphere','m','lola-ldem4',(-9000,11000),'viridis',
                    body='Moon',coordinate_frame='DE421 mean-Earth / polar-axis',period='LOLA v3 / 2009–2016',notes='DN multiplied by 0.5 m. Longitudes rotated from 0–360 to -180–180; mean-Earth/polar-axis DE421 coordinates. GRAIL uses a different principal-axes frame; no cross-frame registration is implied.',
                    method='Spherical-area mean; native 0.25-degree resolution'),
               lambda w:regular_grid(data,89.875-np.arange(720)*.25,-179.875+np.arange(1440)*.25,w))


def sinusoidal_pixel_cells(h,v,shape,width):
    """Assign equal-area MODIS pixel centres to output cells, excluding the map lobes."""
    rows,cols=shape
    latitude=np.pi/2-(v+(np.arange(rows)+.5)/rows)*np.pi/18
    x=-np.pi+(h+(np.arange(cols)+.5)/cols)*np.pi/18
    longitude=x[None,:]/np.cos(latitude[:,None])
    valid=(abs(longitude)<=np.pi)&(abs(latitude[:,None])<np.pi/2)
    col=np.floor((longitude+np.pi)*width/(2*np.pi)).astype('int32')
    col=np.clip(col,0,width-1)
    row=np.floor((np.pi/2-latitude)*width/(2*np.pi)).astype('int32')
    return np.where(valid,row[:,None]*width+col,-1).astype('int32')


def accumulate_equal_area(total,weights,values,cells):
    valid=np.isfinite(values)&(cells>=0)
    index=cells[valid]
    total+=np.bincount(index,weights=values[valid],minlength=total.size).reshape(total.shape)
    weights+=np.bincount(index,minlength=total.size).reshape(total.shape)


def productivity(atlas):
    """Stream equal-area pixel-centre sums; no full-resolution global mosaic."""
    width=max(atlas.widths)
    shape=(width//2,width)
    sums={name:np.zeros(shape,'float64') for name in ('Gpp_500m','Npp_500m','Npp_QC_500m')}
    weights={name:np.zeros(shape,'float64') for name in sums}
    files=atlas.paths('mod17-npp-2020','*.hdf')
    if len(files)!=290: raise ValueError('Expected all 290 MOD17 tiles')
    for k,path in enumerate(files):
        match=re.search(r'\.h(\d+)v(\d+)\.',path.name);h,v=map(int,match.groups())
        cells=sinusoidal_pixel_cells(h,v,(2400,2400),width)
        ds=SD(str(path),SDC.READ)
        try:
            for name in sums:
                s=ds.select(name);a=as_float(s[:]);attrs=s.attributes();s.endaccess()
                valid=np.isfinite(a)&(a!=attrs['_FillValue'])
                low,high=attrs['valid_range'];valid&=(a>=low)&(a<=high)
                if not valid.any(): continue
                a*=float(attrs.get('scale_factor',1));a[~valid]=np.nan
                accumulate_equal_area(sums[name],weights[name],a,cells)
        finally: ds.end()
        if k%10==0: print(f'MOD17 tiles {k+1}/{len(files)}',flush=True)
    for name,total in sums.items():
        def make(w):
            factor=width//w
            if width%w: raise ValueError('MOD17 output widths must divide the largest width')
            num=total.reshape(w//2,factor,w,factor).sum(axis=(1,3))
            den=weights[name].reshape(w//2,factor,w,factor).sum(axis=(1,3))
            return np.divide(num,den,out=np.full_like(num,np.nan),where=den>0)
        key={'Gpp_500m':'gpp','Npp_500m':'npp','Npp_QC_500m':'npp_qc'}[name]
        qc=name=='Npp_QC_500m'
        atlas.emit(Layer('biosphere/'+key,{'gpp':'Gross primary productivity','npp':'Net primary productivity','npp_qc':'NPP gap-filled input percentage'}[key],
                        'percent' if qc else 'kg C m-2 yr-1','mod17-npp-2020',(0,100) if qc else (0,3 if key=='gpp' else 1.5),
                        'cmo.amp' if qc else 'cmo.algae',period='2020 annual',method='Mean of equal-area source pixels assigned by their geographic centres; coarser maps combine sums and counts',
                        notes='290 native 500m MOD17 tiles. Reserved codes removed; provider scales applied. Pixel-centre binning approximates cell-edge overlap at 500m precision. QC is percentage of growing-season days with gap-filled LAI/FPAR, not generic confidence.'),make)
