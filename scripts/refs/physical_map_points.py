"""Observation maps: counts and within-cell summaries, with no interpolation."""
import calendar
import csv
import io
import json
import re
import zipfile

import numpy as np
import openpyxl
import pandas as pd

from physical_map_core import Layer, aggregate_points


def number(value):
    try: return float(value)
    except (ValueError,TypeError): return np.nan


def point_layer(atlas,key,title,dataset,lon,lat,values=None,units='records per output cell',limits=(1,100),palette='cmo.amp',scale='log',notes='',period='Provider reference epoch',reducer='mean'):
    layer=Layer(key,title,units,dataset,limits,palette,scale,notes=notes,period=period,
                zero_colour='#111820' if values is None else None,
                method='Record count in each cell; zero = no catalogued observation' if values is None else f'Within-cell {reducer} of valid observations; unsampled cells are nodata')
    atlas.emit(layer,lambda w:aggregate_points(lon,lat,w,values,reducer))


def heatflow(atlas):
    path=atlas.one('ihfc-2024','*.xlsx')
    workbook=openpyxl.load_workbook(path,read_only=True,data_only=True)
    sheet=workbook['data release2024']
    rows=sheet.iter_rows(min_row=6,max_col=74,values_only=True)
    names=next(rows);positions={n:i for i,n in enumerate(names) if n}
    selected=['long_EW','lat_NS','q','q_uncertainty','T_grad_mean','tc_mean']
    records=[]
    for row in rows:
        record=[number(row[positions[n]]) for n in selected]
        if np.isfinite(record[:3]).all(): records.append(record)
    workbook.close()
    a=np.asarray(records);lon,lat=a[:,:2].T
    notes='IHFC measurements, uneven spatial coverage and review status. Repeated measurements retained as records; maps do not estimate unsampled heat flow.'
    point_layer(atlas,'geothermal/observation_count','Heat-flow observation count','ihfc-2024',lon,lat,notes=notes)
    specs=[(2,'heat_flow','Measured surface heat flow','mW m-2',(0,200),'cmo.thermal'),
           (3,'heat_flow_uncertainty','Reported heat-flow uncertainty','mW m-2',(0,50),'cmo.amp'),
           (4,'temperature_gradient','Measured geothermal temperature gradient','K km-1',(0,100),'cmo.thermal'),
           (5,'thermal_conductivity','Measured thermal conductivity','W m-1 K-1',(0,6),'cmo.dense')]
    for i,key,title,units,limits,palette in specs:
        point_layer(atlas,'geothermal/'+key,title,'ihfc-2024',lon,lat,a[:,i],units,limits,palette,'linear',notes)


def volcanoes(atlas):
    p=atlas.one('smithsonian-votw-2026-09','Holocene_Volcanoes.geojson')
    fs=json.loads(p.read_text(encoding='utf-8'))['features']
    lon=[f['properties']['Longitude'] for f in fs];lat=[f['properties']['Latitude'] for f in fs]
    base='Recorded Holocene volcanoes; not an eruption probability or hazard assessment.'
    point_layer(atlas,'volcanism/holocene_count','Holocene volcano count','smithsonian-votw-2026-09',lon,lat,notes=base)
    for field,key,title,units,limits in [('Elevation','elevation','Volcano summit elevation','m',(0,6500)),
                                        ('Last_Eruption_Year','last_eruption','Latest recorded eruption year','calendar year',(1500,2026))]:
        point_layer(atlas,'volcanism/'+key,title,'smithsonian-votw-2026-09',lon,lat,[number(f['properties'][field]) for f in fs],
                    units,limits,'cmo.thermal','linear',base,reducer='max')
    coords={f['properties']['Volcano_Number']:(f['properties']['Longitude'],f['properties']['Latitude']) for f in fs}
    p=atlas.one('smithsonian-votw-2026-09','Holocene_Eruptions.geojson')
    rows=[f['properties'] for f in json.loads(p.read_text(encoding='utf-8'))['features']]
    rows=[r for r in rows if r['Volcano_Number'] in coords]
    xy=np.array([coords[r['Volcano_Number']] for r in rows])
    point_layer(atlas,'volcanism/eruption_record_count','Recorded Holocene eruptions','smithsonian-votw-2026-09',xy[:,0],xy[:,1],notes='Catalogue event counts at volcano locations; recording completeness varies strongly with time and place.')
    point_layer(atlas,'volcanism/max_vei','Maximum recorded volcanic explosivity index','smithsonian-votw-2026-09',xy[:,0],xy[:,1],
                [number(r['ExplosivityIndexMax']) for r in rows],'VEI',(0,8),'cmo.thermal','linear',base,reducer='max')
    p=atlas.one('smithsonian-votw-2026-09','Pleistocene_Volcanoes.geojson')
    fs=json.loads(p.read_text(encoding='utf-8'))['features']
    point_layer(atlas,'volcanism/pleistocene_count','Pleistocene volcano count','smithsonian-votw-2026-09',
                [f['properties']['Longitude'] for f in fs],[f['properties']['Latitude'] for f in fs],notes='Provider Pleistocene volcano inventory.')


def deposits(atlas):
    p=atlas.one('mrds-legacy','*.zip')
    with zipfile.ZipFile(p) as z, z.open('mrds.csv') as stream:
        data=pd.read_csv(stream,usecols=['latitude','longitude','commod1','commod2','commod3','score'],low_memory=False)
    lon=pd.to_numeric(data.longitude,errors='coerce').to_numpy();lat=pd.to_numeric(data.latitude,errors='coerce').to_numpy()
    note='Legacy MRDS catalogued occurrences; observation counts do not measure ore reserves, undiscovered deposits or resource absence.'
    point_layer(atlas,'resources/mrds_count','Known mineral occurrences','mrds-legacy',lon,lat,notes=note)
    commodities=data[['commod1','commod2','commod3']].fillna('').agg(','.join,axis=1)
    for commodity in ('Copper','Iron','Gold','Uranium','Lithium'):
        selected=commodities.str.contains(r'\b'+commodity+r'\b',case=False,regex=True).to_numpy()
        point_layer(atlas,'resources/mrds_'+commodity.lower(),'Known '+commodity.lower()+' occurrences','mrds-legacy',lon[selected],lat[selected],notes=note)
    # The CSV embeds the USGS record-quality score; preserve the separate explanation as provenance.
    atlas.paths('mrds-quality-2019','mrds-grade.*')
    for score in ('A','B','C','D','E'):
        selected=(data.score==score).to_numpy()
        point_layer(atlas,'resources/mrds_quality_'+score.lower(),'MRDS records: quality '+score,'mrds-legacy',lon[selected],lat[selected],
                    notes=note+' Record-quality group '+score+' retained without inventing a numeric quality average.')
    p=atlas.one('usgs-porphyry-copper-2025','Porphyry_datasheet.csv')
    data=pd.read_csv(p,encoding='utf-8',encoding_errors='replace')
    lon=pd.to_numeric(data.LONGITUDE,errors='coerce');lat=pd.to_numeric(data.LATITUDE,errors='coerce')
    point_layer(atlas,'resources/porphyry_copper_count','Porphyry copper deposits','usgs-porphyry-copper-2025',lon,lat,notes='Complete downloaded release; points represent known deposits.')
    for field,key,title,units,limits in [('CU_PERCENT','copper_grade','Reported porphyry copper grade','percent',(0,2)),
                                       ('ORE_TONNAGE_MT','ore_tonnage','Reported porphyry ore tonnage','million tonnes',(1,10000)),
                                       ('ASSIGNED_AGE_MA','formation_age','Assigned porphyry formation age','Ma',(0,300))]:
        point_layer(atlas,'resources/porphyry_'+key,title,'usgs-porphyry-copper-2025',lon,lat,pd.to_numeric(data[field],errors='coerce'),
                    units,limits,'cmo.matter','log' if key=='ore_tonnage' else 'linear',
                    'Within-cell mean of reported deposit attributes; reporting dates and resource classifications differ.')
    p=atlas.one('usmin-lithium-2020','*_CSV.zip')
    with zipfile.ZipFile(p) as z:
        member=next(n for n in z.namelist() if n.endswith('/Loc_Pt.csv'))
        with z.open(member) as f: data=pd.read_csv(f)
    latcol=next(c for c in data if c.startswith('Lat_WGS84'));loncol=next(c for c in data if c.startswith('Long_WGS84'))
    point_layer(atlas,'resources/usmin_lithium_features','Mapped features in USMIN lithium release','usmin-lithium-2020',
                pd.to_numeric(data[loncol],errors='coerce'),pd.to_numeric(data[latcol],errors='coerce'),
                notes='Features associated with the selected US lithium release; multiple features can belong to one site and some have associated commodities. Not a count of separate lithium deposits.')


GEOROC_COLUMNS=['LATITUDE MIN','LATITUDE MAX','LONGITUDE MIN','LONGITUDE MAX','MATERIAL','SIO2(WT%)','MGO(WT%)','CAO(WT%)']


def whole_rock_rows(data):
    # CSV usecols preserves file order, which can place CaO before MgO.
    numeric=data[[c for c in GEOROC_COLUMNS if c!='MATERIAL']].apply(pd.to_numeric,errors='coerce')
    keep=data.MATERIAL.fillna('').str.match(r'^WR\b')
    keep&=(numeric['LATITUDE MAX']-numeric['LATITUDE MIN']).abs()<=.5
    keep&=(numeric['LONGITUDE MAX']-numeric['LONGITUDE MIN']).abs()<=.5
    return numeric.loc[keep].to_numpy()


def geochemistry(atlas):
    for rock in ('BASALT','GRANITE','GABBRO','PERIDOTITE','RHYOLITE'):
        batches=[]
        for path in atlas.paths('georoc-rocks-2026-09','*'+rock+'*.csv'):
            for data in pd.read_csv(path,usecols=GEOROC_COLUMNS,chunksize=20000,encoding='utf-8',encoding_errors='replace',low_memory=False):
                batches.append(whole_rock_rows(data))
        a=np.concatenate(batches)
        lat=(a[:,0]+a[:,1])*.5;lon=(a[:,2]+a[:,3])*.5
        notes='GEOROC selected '+rock.lower()+' compilation, whole-rock (WR) records only; coordinate bounds no wider than 0.5 degrees; midpoint locations. Individual records are not guaranteed independent specimens.'
        point_layer(atlas,'geochemistry/'+rock.lower()+'/count',rock.title()+' whole-rock observations','georoc-rocks-2026-09',lon,lat,notes=notes)
        for i,oxide,limits in [(4,'silica',(0,80)),(5,'magnesium_oxide',(0,50)),(6,'calcium_oxide',(0,20))]:
            values=a[:,i].copy();values[(values<0)|(values>100)]=np.nan
            point_layer(atlas,'geochemistry/'+rock.lower()+'/'+oxide,rock.title()+' '+oxide.replace('_',' '),'georoc-rocks-2026-09',lon,lat,values,
                        'wt%',limits,'cmo.matter','linear',notes)


def discharge(atlas):
    path=atlas.one('grdc-wmo-2024','*.zip')
    records=[]
    with zipfile.ZipFile(path) as z:
        members=[n for n in z.namelist() if 'Month' in n and n.endswith('.txt')]
        if not members: raise ValueError('Monthly GRDC station files not found')
        for name in members:
            text=z.read(name).decode('utf-8',errors='replace')
            def header(label):
                match=re.search(r'^#\s*'+label+r'[^:]*:\s*([-\d.]+)',text,re.M)
                return float(match.group(1)) if match else np.nan
            lat,lon=header('Latitude'),header('Longitude')
            lines=[line for line in text.splitlines() if line and not line.startswith('#')]
            if not lines:continue
            total=days=count=0
            for row in csv.DictReader(lines,delimiter=';'):
                row={k.strip():v.strip() for k,v in row.items() if k}
                date=row.get('YYYY-MM-DD','')
                if not re.match(r'20\d\d-\d\d-\d\d',date):continue
                year,month=map(int,date[:7].split('-'))
                value=number(row.get('Calculated',row.get('Original')))
                if 2001<=year<=2020 and np.isfinite(value) and value>=0:
                    duration=calendar.monthrange(year,month)[1]
                    total+=value*duration;days+=duration;count+=1
            if count>=120: records.append((lon,lat,total/days,count))
    if not records: raise ValueError('No GRDC stations meet the 120-month reference threshold')
    a=np.asarray(records)
    notes='GRDC stations with at least 120 valid months in 2001–2020; calendar-day weighted monthly mean. Unequal station coverage. CC BY-NC 4.0 and contributing-provider terms apply.'
    point_layer(atlas,'hydrology/gauged_discharge','Observed mean river discharge at gauges','grdc-wmo-2024',a[:,0],a[:,1],a[:,2],
                'm3 s-1',(1,100000),'cmo.deep','log',notes,period='2001–2020 available months')
    point_layer(atlas,'hydrology/gauge_valid_months','Gauge reference coverage','grdc-wmo-2024',a[:,0],a[:,1],a[:,3],
                'valid months',(120,240),'viridis','linear',notes,period='2001–2020')
