"""Vector sources rendered on the atlas grid with explicit categorical/line semantics."""
import colorsys
import json
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd
from rasterio.features import rasterize
from rasterio.transform import from_bounds

from physical_map_core import Layer, REPO, GLOBAL, slug
from physical_map_points import point_layer


def vector_file(atlas,archive,member):
    stem=member[:-4]
    with zipfile.ZipFile(archive) as z:
        names=set(z.namelist())
    for ext in ('.shp','.shx','.dbf','.prj','.cpg'):
        if stem+ext in names: atlas.unzip(archive,stem+ext)
    return atlas.unzip(archive,member)


def read_vector(atlas,archive,member,**kwargs):
    path=vector_file(atlas,archive,member)
    data=gpd.read_file(path,engine='pyogrio',**kwargs)
    if data.crs is None: raise ValueError('Vector source has no CRS: '+member)
    return data.to_crs(4326)


def draw_geometry(geometry,values,width,all_touched=False,fill=np.nan,maximum=False):
    values=np.asarray(values,float)
    valid=np.isfinite(values)&np.asarray([g is not None and not g.is_empty for g in geometry])
    order=np.flatnonzero(valid)
    if maximum: order=order[np.argsort(values[order],kind='stable')]
    if not len(order): return np.full((width//2,width),fill,'float32')
    return rasterize(((geometry[i],values[i]) for i in order),out_shape=(width//2,width),
                     transform=from_bounds(*GLOBAL,width,width//2),fill=fill,dtype='float32',all_touched=all_touched)


def class_colours(labels,identities=None):
    identities=range(len(labels)) if identities is None else identities
    return {i+1:(label,'#'+''.join(f'{round(c*255):02x}' for c in colorsys.hsv_to_rgb((identity*.61803398875)%1,.5,.9)))
            for i,(label,identity) in enumerate(zip(labels,identities))}


def geology_uk(atlas):
    archive=atlas.one('bgs-geology-625k','*.zip')
    with zipfile.ZipFile(archive) as z: names=z.namelist()
    for part,pattern in [('bedrock','BEDROCK_Geology_Polygons.shp'),('superficial','SUPERFICIAL_Geology_Polygons.shp')]:
        member=next(n for n in names if n.endswith(pattern))
        data=read_vector(atlas,archive,member)
        labels=sorted(data.RCS_D.fillna('Unspecified').unique())
        codes={name:i+1 for i,name in enumerate(labels)}
        values=data.RCS_D.fillna('Unspecified').map(codes).to_numpy()
        geoms=data.geometry.to_numpy()
        atlas.emit(Layer('geology/uk_'+part+'_lithology','UK '+part+' lithology','BGS rock-class index','bgs-geology-625k',(1,len(codes)),
                        classes=class_colours(labels),palette_status='Author-selected distinct categorical colours; original BGS rock labels preserved',
                        palette_source='https://www.bgs.ac.uk/datasets/bgs-geology-625k/',
                        method='Source polygon at target cell centre',notes='UK-only 1:625000 reference. Category lookup is in this manifest; categories are not geological ages.'),
                   lambda w:draw_geometry(geoms,values,w))
        if part=='bedrock':
            official=json.loads((REPO/'assets/palettes/physical/ics-periods.json').read_text())
            periods=data.MAX_PERIOD.fillna('UNSPECIFIED').replace({'PALAEOGENE':'PALEOGENE'})
            period_names=sorted(periods.unique())
            period_codes={name:i+1 for i,name in enumerate(period_names)}
            classes={period_codes[name]:(name.title(),official['colours'].get(name,'#777777')) for name in period_names}
            period_values=periods.map(period_codes).to_numpy()
            atlas.emit(Layer('geology/uk_oldest_period','Oldest mapped geological period, UK bedrock','indexed period','bgs-geology-625k',(1,len(period_codes)),
                            classes=classes,method='BGS MAX_PERIOD attribute at cell centre',
                            palette_status='Official ICS period colours; unspecified/unmatched intervals use neutral grey',
                            palette_source=official['source'],notes='Oldest end of each BGS unit age range; a unit may extend across more than one period. Grey intervals have no matched ICS period colour.'),
                       lambda w:draw_geometry(geoms,period_values,w))
            for field,key,title in [('MAX_TIME_Y','oldest_age','Oldest mapped bedrock age'),('MIN_TIME_Y','youngest_age','Youngest mapped bedrock age')]:
                age=pd.to_numeric(data[field],errors='coerce').to_numpy()/1e6
                atlas.emit(Layer('geology/uk_'+key,title,'Ma','bgs-geology-625k',(0,3000),'cmo.thermal_r',
                                method='Source polygon attribute at cell centre',notes='Limits of mapped unit age from BGS, not a single crystallisation age.'),
                           lambda w,a=age:draw_geometry(geoms,a,w))
    member=next(n for n in names if n.endswith('FAULT_Geology_Lines.shp'))
    data=read_vector(atlas,archive,member);geoms=data.geometry.to_numpy()
    atlas.emit(Layer('geology/uk_faults','Mapped UK geological faults','line-cell presence','bgs-geology-625k',(0,1),
                    classes={1:('Mapped fault','#ec8477')},method='Any source fault line touching the cell',
                    notes='Regional geological fault map; not earthquake hazard or active-fault classification.',palette_status='Author-selected line colour'),
               lambda w:draw_geometry(geoms,np.ones(len(geoms)),w,all_touched=True))


def groundwater(atlas):
    archive=atlas.one('whymap-gwr-v1','*.zip')
    with zipfile.ZipFile(archive) as z: member=next(n for n in z.namelist() if n.endswith('GW_aquifers_v1_poly.shp'))
    data=read_vector(atlas,archive,member);geoms=data.geometry.to_numpy();values=data.HYGEO2.to_numpy()
    metadata=json.loads((REPO/'assets/palettes/physical/whymap_aquifers.json').read_text())
    classes={}
    family={1:'Major groundwater basin',2:'Complex hydrogeology',3:'Local / shallow aquifers'}
    for item in metadata['renderer']['uniqueValueInfos']:
        colour=item['symbol']['color'];code=int(item['value'])
        if colour[3]==0:continue
        classes[code]=(family[code//10]+': recharge '+item['label']+' mm/yr','#'+''.join(f'{v:02x}' for v in colour[:3]))
    values=np.where(np.isin(values,list(classes)),values,np.nan)
    atlas.emit(Layer('hydrology/groundwater_resources','Groundwater resources and recharge classes','WHYMAP category','whymap-gwr-v1',(11,34),classes=classes,
                    method='Provider hydrogeological polygons at cell centres',
                    palette_status='Official BGR WHYMAP ArcGIS renderer colours and class labels',palette_source=metadata['source'],
                    notes='Regional aquifer/recharge classes, not local water-table depth or available extractable storage.'),lambda w:draw_geometry(geoms,values,w))


def basins(atlas):
    frames=[]
    for archive in atlas.paths('hydrobasins-v1c','*.zip'):
        with zipfile.ZipFile(archive) as z:member=next(n for n in z.namelist() if n.endswith('lev03_v1c.shp'))
        frames.append(read_vector(atlas,archive,member,columns=['HYBAS_ID','ENDO','UP_AREA']))
    data=pd.concat(frames,ignore_index=True).sort_values('HYBAS_ID').drop_duplicates('HYBAS_ID')
    geoms=data.geometry.to_numpy();ids=data.HYBAS_ID.astype('int64').to_list();codes=np.arange(1,len(data)+1)
    classes=class_colours([str(v) for v in ids])
    lookup=atlas.output/'basin-id-lookup.json'
    lookup.write_text(json.dumps({int(code):int(original) for code,original in zip(codes,ids)},indent=2)+'\n')
    atlas.emit(Layer('hydrology/basin_level3','Drainage basins, HydroBASINS level 3','indexed basin ID','hydrobasins-v1c',(1,len(ids)),classes=classes,
                    palette_status='Author-selected distinct basin colours; colour has no physical magnitude',
                    method='Polygon at cell centre; reversible sequential IDs avoid float32 rounding of 11-digit identifiers',
                    notes='Original HYBAS_ID values retained in basin-id-lookup.json. Smaller basins can disappear at coarse resolutions.'),
               lambda w:draw_geometry(geoms,codes,w))
    internal=(data.ENDO.to_numpy()>0).astype('float32')
    atlas.emit(Layer('hydrology/internal_drainage','Internally drained catchments','drainage class','hydrobasins-v1c',(0,1),
                    classes={0:('Externally drained','#78bde8'),1:('Internally drained','#e8bc70')},
                    method='Provider ENDO > 0, source polygon at cell centre',palette_status='Author-selected categorical colours'),
               lambda w:draw_geometry(geoms,internal,w))


def rivers(atlas):
    archive=atlas.one('hydrorivers-v10','*.zip')
    with zipfile.ZipFile(archive) as z:member=next(n for n in z.namelist() if n.endswith('.shp'))
    data=read_vector(atlas,archive,member,columns=['UPLAND_SKM','DIS_AV_CMS','ORD_STRA'],where='UPLAND_SKM >= 1000')
    geoms=data.geometry.to_numpy()
    note='HydroRIVERS reaches with upstream area >=1000 km2, selected for global 128–2048 px display. Numeric cells contain the maximum touching-reach value, not an area-mean field. Discharge is a provider model estimate, not gauge observations.'
    for field,key,title,units,limits,palette,scale in [('DIS_AV_CMS','discharge','Major-river estimated discharge','m3 s-1',(1,100000),'cmo.deep','log'),
                                                  ('UPLAND_SKM','upstream_area','Major-river upstream drainage area','km2',(1000,7000000),'cmo.deep','log'),
                                                  ('ORD_STRA','strahler_order','Major-river Strahler order','order',(1,12),'viridis','linear')]:
        values=data[field].to_numpy()
        atlas.emit(Layer('hydrology/rivers_'+key,title,units,'hydrorivers-v10',limits,palette,scale,
                        method='Maximum attribute among retained source reaches touching each target cell',notes=note),
                   lambda w,a=values:draw_geometry(geoms,a,w,all_touched=True,maximum=True))


def lakes(atlas):
    archive=atlas.one('hydrolakes-v10','*.zip')
    with zipfile.ZipFile(archive) as z: member=next(n for n in z.namelist() if n.endswith('.shp'))
    path=vector_file(atlas,archive,member)
    data=gpd.read_file(path,engine='pyogrio',columns=['Pour_long','Pour_lat','Lake_area','Depth_avg','Vol_total','Dis_avg'],ignore_geometry=True)
    lon=data.Pour_long.to_numpy();lat=data.Pour_lat.to_numpy()
    notes='HydroLAKES pour-point locations; point attributes do not represent a raster of lake footprints. Depth, volume and discharge include provider estimates.'
    point_layer(atlas,'hydrology/lake_count','Catalogued lake pour-point count','hydrolakes-v10',lon,lat,notes=notes)
    for field,key,title,units,limits in [('Lake_area','area','Lake surface area at pour points','km2',(1,400000)),
                                       ('Depth_avg','depth','Estimated mean lake depth at pour points','m',(1,1000)),
                                       ('Vol_total','volume','Estimated lake volume at pour points','million m3',(1,100000000))]:
        values=data[field].to_numpy();values=np.where(values>=0,values,np.nan)
        point_layer(atlas,'hydrology/lake_'+key,title,'hydrolakes-v10',lon,lat,values,units,limits,'cmo.deep','log',notes,reducer='max')


def glaciers(atlas):
    archive=atlas.one('rgi7','*.zip')
    masks={w:np.zeros((w//2,w),'float32') for w in atlas.widths}
    with zipfile.ZipFile(archive) as z: members=[n for n in z.namelist() if n.endswith('.zip')]
    for i,member in enumerate(members):
        regional=atlas.unzip(archive,member)
        with zipfile.ZipFile(regional) as z: shp=next(n for n in z.namelist() if n.endswith('.shp'))
        data=read_vector(atlas,regional,shp,columns=[])
        geoms=data.geometry.to_numpy()
        for w in atlas.widths:
            masks[w]=np.maximum(masks[w],draw_geometry(geoms,np.ones(len(geoms)),w,all_touched=True,fill=0))
        print(f'RGI regions {i+1}/{len(members)}',flush=True)
    atlas.emit(Layer('cryosphere/glacier_presence','Glacier polygons intersecting each cell','polygon-cell presence','rgi7',(0,1),
                    classes={0:('No mapped glacier intersection','#101820'),1:('Mapped glacier intersection','#bcecff')},
                    palette_status='Author-selected ice-blue convention',period='RGI7 nominal year 2000',
                    method='Any RGI glacier polygon touching the output cell',
                    notes='Presence map, not glacier area fraction or ice thickness. All-touched rendering retains narrow glaciers and enlarges their apparent footprint at coarse resolutions. Ice sheets excluded.'),
               lambda w:masks[w])
