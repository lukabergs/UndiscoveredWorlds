"""Aggregate one benchmark for the validation workbook; never reads old run results.

Usage: python scripts/summarize-climate-validation.py --run-id 159 --output DIR
Numeric references are area-averaged onto each diagnostic grid. Excel receives
weighted moments, quantiles and paired reference values, not rendered map colours
(except the source categorical Koppen/AR6 masks).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import struct
from pathlib import Path
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "extra/reference/climate/processed"
RADIUS = 6371000.0
DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31.], float)
QDAYS = DAYS.reshape(4, 3).sum(1)
QUARTERS = ["Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec"]
MONTHS = ["Jan", "Apr", "Jul", "Oct"]
RECORDS, VECTORS, FEATURES = [], [], []
REF_CACHE, REGION_CACHE = {}, {}

def read_csv(p):
    with Path(p).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def bounds(n, pole=False):
    if not pole:
        return np.linspace(90, -90, n + 1)
    lat = np.linspace(90, -90, n)
    return np.r_[90, (lat[:-1] + lat[1:]) / 2, -90]

def remap(a, ny, nx, pole=True, source_pole=True):
    sy, sx = a.shape[-2:]
    sb, tb = np.sin(np.deg2rad(bounds(sy, source_pole))), np.sin(np.deg2rad(bounds(ny, pole)))
    lat = np.maximum(0, np.minimum(tb[:-1, None], sb[None, :-1]) - np.maximum(tb[1:, None], sb[None, 1:]))
    tx, ss = np.linspace(0, 1, nx + 1), np.linspace(0, 1, sx + 1)
    lon = np.maximum(0, np.minimum(tx[1:, None], ss[None, 1:]) - np.maximum(tx[:-1, None], ss[None, :-1]))
    def one(v):
        finite = np.isfinite(v)
        if sx % nx == 0:
            k = sx // nx
            num = np.where(finite, v, 0).reshape(sy, nx, k).sum(2)
            den = finite.reshape(sy, nx, k).sum(2)
            num, den = lat @ num, lat @ den
        else:
            num, den = lat @ np.where(finite, v, 0) @ lon.T, lat @ finite.astype(float) @ lon.T
        return np.divide(num, den, out=np.full((ny, nx), np.nan), where=den > 0)
    if a.ndim == 2:
        return one(a)
    return np.stack([one(v) for v in a])

def reference(name, ny, nx, pole=True):
    key = (name, ny, nx, pole)
    if key not in REF_CACHE:
        with (REF / name).open('rb') as f:
            assert f.read(8) == b'UWCLIM1\0'
            version, w, h, months = struct.unpack('<IIII', f.read(16))
            assert version == 1
            f.read(16)
            a = np.fromfile(f, '<f4').reshape(months, h, w)
        REF_CACHE[key] = remap(a, ny, nx, pole)
    return REF_CACHE[key]

def grid(ny, nx, pole):
    b = bounds(ny, pole)
    lat = np.linspace(90, -90, ny) if pole else 90 - (np.arange(ny) + .5) * 180 / ny
    lon = -180 + (np.arange(nx) + .5) * 360 / nx
    area = np.broadcast_to((np.sin(np.deg2rad(b[:-1])) - np.sin(np.deg2rad(b[1:])))[:, None] * 2 * np.pi * RADIUS**2 / nx / 1e6, (ny, nx)).copy()
    assert abs(area.sum() / (4*np.pi*RADIUS**2/1e6) - 1) < 1e-12
    return lat, lon, area

def regions(ny, nx, pole, land, mask_key='output'):
    key = ny, nx, pole, mask_key
    if key in REGION_CACHE:
        return REGION_CACHE[key]
    lat, lon, area = grid(ny, nx, pole)
    ocean = 1-land
    b = bounds(ny, pole)
    band = np.sin(np.deg2rad(b[:-1]))-np.sin(np.deg2rad(b[1:]))
    def latitude_fraction(south, north):
        overlap = np.maximum(0,np.sin(np.deg2rad(np.minimum(b[:-1],north)))-np.sin(np.deg2rad(np.maximum(b[1:],south))))
        return np.broadcast_to((overlap/band)[:,None],land.shape)
    def longitude_fraction(west,east):
        edges=np.linspace(-180,180,nx+1)
        def part(a,z):return np.maximum(0,np.minimum(edges[1:],z)-np.maximum(edges[:-1],a))/(360/nx)
        return part(west,east) if west<east else part(west,180)+part(-180,east)
    masks = {'GLOBAL': np.ones_like(land), 'LAND': land, 'OCEAN': ocean,
             'NH': latitude_fraction(0,90), 'SH': latitude_fraction(-90,0)}
    for lo, hi, label in [(-10,10,'EQUATOR'),(-35,-10,'S_SUBTROPICS'),(10,35,'N_SUBTROPICS'),(-60,-35,'S_MIDLAT'),(35,60,'N_MIDLAT'),(-90,-60,'S_POLAR'),(60,90,'N_POLAR')]:
        m = latitude_fraction(lo,hi)
        masks[label] = m
        masks[label+'_LAND'] = m*land
    for name, west, east, south, north in [('N_ATLANTIC',-80,20,10,65),('S_ATLANTIC',-70,20,-60,0),('N_PACIFIC',120,-70,10,65),('S_PACIFIC',150,-70,-60,0),('INDIAN',20,120,-60,30),('SOUTHERN_OCEAN',-180,180,-90,-50),('ARCTIC_OCEAN',-180,180,65,90),('TROPICAL_ATLANTIC',-70,20,-20,20),('TROPICAL_PACIFIC',120,-70,-20,20),('TROPICAL_INDIAN',20,120,-20,20)]:
        masks[name] = ocean*latitude_fraction(south,north)*longitude_fraction(west,east)[None,:]
    source = np.asarray(Image.open(REF/'ipcc_ar6_regions.png').convert('RGB'))[:,:,0].astype(int)-1
    for item in AR6:
        fraction=np.clip(remap((source==int(item['id'])).astype(float),ny,nx,pole),0,1)
        if item['type']=='Land':fraction*=land
        elif item['type']=='Ocean':fraction*=ocean
        masks['AR6_'+item['acronym']] = fraction
    REGION_CACHE[key] = (masks,area,lat,lon)
    return REGION_CACHE[key]

def quantiles(a,w):
    order = np.argsort(a)
    cumulative = np.cumsum(w[order]) / w.sum()
    return [float(np.interp(p,cumulative,a[order])) for p in (.1,.5,.9,.95)]

def moments(s,r,w):
    good = np.isfinite(s) & (w>0)
    if r is not None:
        good &= np.isfinite(r)
    if not good.any(): return None
    ss,ww = s[good],w[good]
    total = ww.sum(); norm = ww/total
    sm = norm@ss
    out = dict(sim_mean=float(sm),sim_variance=float(norm@((ss-sm)**2)),sim_quantiles=quantiles(ss,ww),
               ref_mean=None,ref_variance=None,covariance=None,mae=None,mse=None,ref_quantiles=[None]*4,
               area_km2=float(total),cells=int(good.sum()))
    if r is not None:
        rr = r[good]; rm=norm@rr
        out.update(ref_mean=float(rm),ref_variance=float(norm@((rr-rm)**2)),covariance=float(norm@((ss-sm)*(rr-rm))),
                   mae=float(norm@abs(ss-rr)),mse=float(norm@((ss-rr)**2)),ref_quantiles=quantiles(rr,ww))
    return out

def add_field(domain,name,unit,s,r,period,refid,region_data,level='',notes='',positive=False,selected=None):
    masks,area,_,_=region_data
    for region,fraction in masks.items():
        if selected is not None and region not in selected: continue
        m = moments(s,r,area*fraction)
        if m is None: continue
        effective = float((fraction*(np.isfinite(s) & (np.isfinite(r) if r is not None else True))).sum())
        status = 'Internal check' if refid == 'INTERNAL_CHECK' else ('Comparable' if r is not None else 'Reference needed')
        if effective<6: status += '; coarse coverage'
        RECORDS.append(dict(domain=domain,variable=name,units=unit,region=region,period=period,level=level,
            reference_id=refid,status=status,notes=notes,positive=positive,effective_cells=effective,**m))

def add_summary(domain,name,unit,s,r,period,refid,region,area_km2,cells,notes,level=''):
    RECORDS.append(dict(domain=domain,variable=name,units=unit,region=region,period=period,level=level,
        reference_id=refid,status='Regional feature summary',notes=notes,positive=False,effective_cells=float(cells),
        sim_mean=s,ref_mean=r,sim_variance=None,ref_variance=None,covariance=None,mae=None,mse=None,
        sim_quantiles=[None]*4,ref_quantiles=[None]*4,area_km2=float(area_km2),cells=int(cells)))

def wettest_share(a,w):
    order=np.argsort(a)[::-1];a=a[order];w=w[order]
    allowed=np.clip(.1*w.sum()-np.r_[0,np.cumsum(w[:-1])],0,w)
    total=np.sum(a*w)
    return float(np.sum(a*allowed)/total) if total>0 else None

def periods(a, annual=True):
    result = list(a)
    if annual: result.append(np.average(a,axis=0,weights=QDAYS))
    return result

def quarter_means(a):
    return np.stack([np.average(a[s*3:s*3+3],axis=0,weights=DAYS[s*3:s*3+3]) for s in range(4)])

def gradients(a,lat):
    # Spherical metric derivatives on the common comparison grid. Polar caps excluded.
    ny,nx=a.shape
    phi=np.deg2rad(lat); dl=2*np.pi/nx
    dx=(np.roll(a,-1,axis=1)-np.roll(a,1,axis=1))/(2*dl*RADIUS*np.maximum(np.cos(phi[:,None]),.03))
    dy=np.gradient(a,phi,axis=0)/RADIUS
    dx[abs(lat)>80]=np.nan; dy[abs(lat)>80]=np.nan
    return dx,dy

def rotation(u,v,lat):
    cos=np.cos(np.deg2rad(lat[:,None])); dxv,_=gradients(v,lat); dxu,_=gradients(u,lat)
    _,dyu=gradients(u*cos,lat); _,dyv=gradients(v*cos,lat)
    return dxu+dyv/np.maximum(cos,.03), dxv-dyu/np.maximum(cos,.03)

def corr(a,b,w):
    a=a-w@a;b=b-w@b
    d=math.sqrt(max(0,float(w@(a*a)))*max(0,float(w@(b*b))))
    return float(w@(a*b))/d if d>1e-15 else None

def wind_vectors(u,v,ru,rv,period,level,region_data,mode):
    masks,area,_,_=region_data
    for name,m in masks.items():
        good=(m>0)&np.isfinite(u)&np.isfinite(v)&np.isfinite(ru)&np.isfinite(rv)
        if not good.any():continue
        w=(area*m)[good];total=w.sum();w=w/total
        uu,vv,rr,ss=(a[good] for a in (u,v,ru,rv))
        speed,rspeed=np.hypot(uu,vv),np.hypot(rr,ss)
        calm=(speed>=.5)&(rspeed>=.5)
        angle=np.rad2deg(np.arccos(np.clip((uu*rr+vv*ss)/np.maximum(speed*rspeed,1e-30),-1,1)))
        angular=float(w[calm]@angle[calm]/w[calm].sum()) if calm.any() else None
        d=math.sqrt(float(w@(uu*uu+vv*vv))*float(w@(rr*rr+ss*ss)))
        VECTORS.append(dict(region=name,period=period,level=level,mode=mode,sim_speed=float(w@speed),ref_speed=float(w@rspeed),
            vector_mse=float(w@((uu-rr)**2+(vv-ss)**2)),vector_alignment=float(w@(uu*rr+vv*ss))/d if d else None,
            angle_error=angular,u_correlation=corr(uu,rr,w),v_correlation=corr(vv,ss,w),sim_u=float(w@uu),ref_u=float(w@rr),sim_v=float(w@vv),ref_v=float(w@ss),
            sim_rms2=float(w@(uu*uu+vv*vv)),ref_rms2=float(w@(rr*rr+ss*ss)),area_km2=float(total),cells=int(good.sum()),
            effective_cells=float(m[good].sum()),reference_id='ERA5_2001_2020',status='Comparable' if m[good].sum()>=6 else 'Coarse coverage'))

def load_field_table(p):
    rows=read_csv(p);ny=len({r['latitude'] for r in rows});nx=len({r['longitude'] for r in rows})
    keys=[k for k in rows[0] if k not in ('season','layer','latitude','longitude')]
    layers=sorted({int(r.get('layer',0)) for r in rows})
    out={}
    for layer in layers:
        selected=[r for r in rows if int(r.get('layer',0))==layer]
        out[layer]={k:np.array([float(r[k]) for r in selected]).reshape(4,ny,nx) for k in keys}
    return ny,nx,out

def identified_ocean_mask(fields):
    # The solver leaves these fields zero on land. A completely stationary,
    # ice-free wet cell is ambiguous until the native wet mask is exported.
    return np.any((fields['east_current_mps'] != 0) |
                  (fields['south_current_mps'] != 0) |
                  (fields['ice_thickness_m'] > 0), axis=0)

def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.generic):v=v.item()
    if isinstance(v,float) and not math.isfinite(v):return None
    return v

def main():
    global AR6
    with (REF/'ipcc_ar6_regions.tsv').open(encoding='utf-8-sig',newline='') as f:
        AR6=list(csv.DictReader(f,delimiter='\t'))
    parser=argparse.ArgumentParser();parser.add_argument('--run-id',required=True,type=int);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args(); run=ROOT/'extra/climate/benchmarks/runs'/str(args.run_id)
    data=read_csv(run/'climate_comparison_cells.csv')
    ny=max(int(r['y']) for r in data)+1;nx=max(int(r['x']) for r in data)+1
    cell={k:np.array([float(r[k]) for r in data]).reshape(4,ny,nx) for k in data[0]}
    land=cell['land'][0]; reg=regions(ny,nx,True,land); masks,area,lat,lon=reg
    wc_t=reference('worldclim_tavg_monthly.uwclim',ny,nx)
    wc_p=reference('worldclim_prec_monthly.uwclim',ny,nx)
    im_p=reference('imerg_prec_monthly.uwclim',ny,nx)
    era_p=reference('era5_pr_monthly.uwclim',ny,nx)
    pressure=reference('era5_slp_anom_monthly.uwclim',ny,nx)[[0,3,6,9]]
    simtemp=cell['temperature_c'];rain=cell['precipitation_mm_month']
    # Land-only WorldClim temperature and precipitation; native stored temperatures labelled.
    for s,name in enumerate(MONTHS):
        add_field('Temperature','Air temperature','degC',np.where(land>0,simtemp[s],np.nan),wc_t[s*3],name,'WORLDCLIM_1970_2000',reg,notes='Stored near-surface temperature; representative month; quantized to 1 degC.')
    for name,s,r in [('Four-month mean',simtemp.mean(0),wc_t[[0,3,6,9]].mean(0)),('Warmest sampled month',simtemp.max(0),wc_t[[0,3,6,9]].max(0)),('Coldest sampled month',simtemp.min(0),wc_t[[0,3,6,9]].min(0)),('Sampled annual range',np.ptp(simtemp,axis=0),np.ptp(wc_t[[0,3,6,9]],axis=0))]:
        add_field('Temperature',name,'degC',np.where(land>0,s,np.nan),r,'Annual proxy','WORLDCLIM_1970_2000',reg,notes='Both sides use Jan/Apr/Jul/Oct; not twelve independent monthly states.')
    for refid,refp in [('WORLDCLIM_1970_2000',wc_p),('IMERG_2001_2022',im_p),('ERA5_2001_2020',era_p)]:
        for season,name in enumerate(QUARTERS):
            sim=rain[season]*3; rr=refp[season*3:season*3+3].sum(0)
            if refid.startswith('WORLDCLIM'):sim=np.where(land>0,sim,np.nan)
            add_field('Rainfall','Precipitation total','mm',sim,rr,name,refid,reg,positive=True)
        sim=rain.sum(0)*3
        if refid.startswith('WORLDCLIM'):sim=np.where(land>0,sim,np.nan)
        add_field('Rainfall','Annual precipitation','mm/year',sim,refp.sum(0),'Annual',refid,reg,positive=True)
    impq=im_p.reshape(4,3,ny,nx).sum(1)
    for region,fraction in masks.items():
        s=rain.sum(0)*3;r=im_p.sum(0);good=(fraction>0)&np.isfinite(s)&np.isfinite(r)
        if not good.any():continue
        w=(area*fraction)[good]
        add_summary('Rainfall','Wettest 10% area rainfall share','fraction',wettest_share(s[good],w),wettest_share(r[good],w),'Annual','IMERG_2001_2022',region,w.sum(),fraction[good].sum(),'Rank cells by rainfall on each side separately, using identical valid area; split the cutoff cell to exactly 10% area.')
    sr=rain*3/QDAYS[:,None,None];rr=impq/QDAYS[:,None,None]
    sm=np.average(sr,axis=0,weights=QDAYS);rm=np.average(rr,axis=0,weights=QDAYS)
    scv=np.divide(np.sqrt(np.average((sr-sm)**2,axis=0,weights=QDAYS)),sm,out=np.full_like(sm,np.nan),where=sm>0)
    rcv=np.divide(np.sqrt(np.average((rr-rm)**2,axis=0,weights=QDAYS)),rm,out=np.full_like(rm,np.nan),where=rm>0)
    add_field('Rainfall','Seasonality coefficient of variation','ratio',scv,rcv,'Annual','IMERG_2001_2022',reg,notes='Duration-weighted variability of four quarterly precipitation rates / annual mean rate.')
    for season,period in enumerate(QUARTERS):
        for region in ['TROPICAL_ATLANTIC','TROPICAL_PACIFIC','TROPICAL_INDIAN']:
            fraction=masks[region];a=rain[season]*3/QDAYS[season];b=impq[season]/QDAYS[season]
            good=np.isfinite(a)&np.isfinite(b)&(fraction>0);w=area*fraction*good
            if not good.any():continue
            def belt(v):
                mass=np.where(good,v,0)*w;total=mass.sum()
                if total<=0:return None,None
                centre=float((mass*lat[:,None]).sum()/total)
                width=float(np.sqrt((mass*(lat[:,None]-centre)**2).sum()/total))
                return centre,width
            sa,sw=belt(a);ra,rw=belt(b)
            for label,sv,rv in [('Rain-belt centroid latitude',sa,ra),('Rain-belt RMS latitude width',sw,rw)]:
                add_summary('Rainfall',label,'degrees',sv,rv,period,'IMERG_2001_2022',region,w.sum(),fraction[good].sum(),'Precipitation-area-weighted centroid/standard deviation within fixed +/-20 degree ocean basin; multiple rain bands may share a centroid.')
            add_summary('Rainfall','Rain-belt mean precipitation','mm/day',float(np.sum(np.where(good,a,0)*w)/w.sum()),float(np.sum(np.where(good,b,0)*w)/w.sum()),period,'IMERG_2001_2022',region,w.sum(),fraction[good].sum(),'Mean precipitation rate within fixed +/-20 degree ocean basin.')
    for name,s,r in [('Wettest quarter',rain.max(0)*3,impq.max(0)),('Driest quarter',rain.min(0)*3,impq.min(0)),('Seasonal rainfall range',np.ptp(rain*3,axis=0),np.ptp(impq,axis=0))]:
        add_field('Rainfall',name,'mm',s,r,'Annual','IMERG_2001_2022',reg,positive=True)
    for threshold in [100,250,1000,2000]:
        add_field('Rainfall',f'Area with annual P < {threshold} mm','fraction',(rain.sum(0)*3<threshold).astype(float),np.where(np.isfinite(im_p.sum(0)),(im_p.sum(0)<threshold).astype(float),np.nan),'Annual','IMERG_2001_2022',reg,notes='Area fraction; spatial error scores compare binary threshold masks.')
    for threshold,seasonal,refseasonal in [(0,simtemp.min(0),wc_t[[0,3,6,9]].min(0)),(10,simtemp.max(0),wc_t[[0,3,6,9]].max(0)),(18,simtemp.min(0),wc_t[[0,3,6,9]].min(0))]:
        add_field('Temperature',f'Area above {threshold} C sampled threshold','fraction',np.where(land>0,(seasonal>threshold).astype(float),np.nan),np.where(np.isfinite(refseasonal),(refseasonal>threshold).astype(float),np.nan),'Annual proxy','WORLDCLIM_1970_2000',reg,notes='0/18 C use sampled coldest month; 10 C uses sampled warmest month. No daily frost claim.')
    # Pressure decomposition preserves both latitude belts and basin departures.
    for p,name in enumerate(MONTHS+['Four-month mean']):
        sp=cell['pressure_anomaly_hpa'][p] if p<4 else cell['pressure_anomaly_hpa'].mean(0)
        rp=pressure[p] if p<4 else pressure.mean(0)
        for mode in ('Full','Zonal mean','Stationary'):
            sm,rm=sp.copy(),rp.copy()
            if mode=='Zonal mean':sm=np.broadcast_to(sp.mean(1)[:,None],sp.shape);rm=np.broadcast_to(rp.mean(1)[:,None],rp.shape)
            if mode=='Stationary':sm=sp-sp.mean(1)[:,None];rm=rp-rp.mean(1)[:,None]
            add_field('Pressure',mode+' pressure anomaly','hPa',sm,rm,name,'ERA5_2001_2020',reg,notes='Stored pressure quantized to 1 hPa; stationarity means remove full latitude-row mean.')
        for component,ss,rr in zip(['Eastward gradient','Northward gradient'],gradients(sp*100,lat),gradients(rp*100,lat)):
            add_field('Pressure',component,'Pa/m',ss,rr,name,'ERA5_2001_2020',reg,selected={'GLOBAL','N_ATLANTIC','N_PACIFIC','INDIAN','SOUTHERN_OCEAN'},notes='Spherical centred differences after common-grid area averaging; exclude abs(latitude)>80.')
        for basin in ['N_ATLANTIC','N_PACIFIC','INDIAN','SOUTHERN_OCEAN']:
            for sign,feature in [(1,'Pressure high'),(-1,'Pressure low')]:
                m=masks[basin]>0
                if not m.any():continue
                si=np.unravel_index(np.nanargmax(np.where(m,sp*sign,np.nan)),sp.shape)
                ri=np.unravel_index(np.nanargmax(np.where(m,rp*sign,np.nan)),rp.shape)
                a,b,c,d=np.deg2rad([lat[si[0]],lon[si[1]],lat[ri[0]],lon[ri[1]]])
                dist=2*RADIUS/1000*np.arcsin(np.sqrt(np.clip(np.sin((c-a)/2)**2+np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2,0,1)))
                FEATURES.append([args.run_id,basin,name,feature,float(lat[si[0]]),float(lon[si[1]]),float(lat[ri[0]]),float(lon[ri[1]]),float(sp[si]),float(rp[ri]),float(dist),'hPa','Grid-cell extremum; broad fixed basin box, not a tracked pressure system.'])
    # Quarter-mean wind transport: use corresponding quarter means from ERA5.
    wy,wx,weather=load_field_table(run/'maps/weather_statistics.csv')
    wland=np.clip(remap(land,wy,wx,False),0,1);wreg=regions(wy,wx,False,wland);wlat=wreg[2]
    windpairs=[]
    for layer,reflevel,label in [(0,'10m','Surface / 10 m proxy'),(1,'500','Upper / 500 hPa proxy')]:
        u,v=weather[layer]['mean_u_mps'],-weather[layer]['mean_v_south_mps']
        ru,rv=(quarter_means(reference(f'era5_{c}{reflevel}_monthly.uwclim',wy,wx,False)) for c in ['u','v'])
        windpairs.append((u,v,ru,rv))
        for i,period in enumerate(QUARTERS+['Annual']):
            base=[a[i] if i<4 else np.average(a,axis=0,weights=QDAYS) for a in [u,v,ru,rv]]
            for mode in ['Full','Zonal mean','Stationary']:
                fields=[a.copy() for a in base]
                if mode=='Zonal mean':fields=[np.broadcast_to(a.mean(1)[:,None],a.shape) for a in fields]
                if mode=='Stationary':fields=[a-a.mean(1)[:,None] for a in fields]
                wind_vectors(*fields,period,label,wreg,mode)
            if layer==0:
                for region,lo,hi,sign,beltname in [('N_SUBTROPICS',10,35,-1,'Trade-wind'),('S_SUBTROPICS',-35,-10,-1,'Trade-wind'),('N_MIDLAT',35,60,1,'Westerly'),('S_MIDLAT',-60,-35,1,'Westerly')]:
                    mask=wreg[0][region];validrow=mask.sum(1)>0
                    def belt(a):
                        speed=np.mean(a,axis=1)*sign
                        at=np.argmax(np.where(validrow,speed,-np.inf));strength=float(speed[at])
                        if strength<=0:return None,0.,None
                        b=bounds(wy,False);widths=np.maximum(0,np.minimum(b[:-1],hi)-np.maximum(b[1:],lo))
                        return float(wlat[at]),strength,float(widths[(speed>=strength/2)&validrow].sum())
                    sa,ss,sw=belt(base[0]);ra,rs,rw=belt(base[2])
                    for var,unit,sim_summary,ref_summary in [('peak latitude','degrees',sa,ra),('peak strength','m/s',ss,rs),('half-maximum width','degrees',sw,rw)]:
                        add_summary('Circulation',beltname+' '+var,unit,sim_summary,ref_summary,period,'ERA5_2001_2020',region,(mask*wreg[1]).sum(),mask.sum(),'Zonal u profile; fixed hemisphere search band. Trade strength is easterly magnitude; width is within the search band.',label)
            sd,sz=rotation(*base[:2],wlat);rd,rz=rotation(*base[2:],wlat)
            for variable,ss,rr in [('Horizontal divergence',sd,rd),('Relative vorticity',sz,rz)]:
                add_field('Circulation',variable,'1/s',ss,rr,period,'ERA5_2001_2020',wreg,level=label,notes='Spherical common-grid derivative; exclude abs(latitude)>80.')
        scalar_ref=np.average(reference('worldclim_wind_monthly.uwclim',wy,wx,False),axis=0,weights=DAYS) if layer==0 else None
        scalar_sim=np.average(weather[layer]['mean_speed_mps'],axis=0,weights=QDAYS)
        if layer==0:scalar_sim=np.where(wland>=.5,scalar_sim,np.nan)
        add_field('Circulation','Mean scalar wind speed','m/s',scalar_sim,scalar_ref,'Annual','WORLDCLIM_1970_2000' if layer==0 else 'NOT_AVAILABLE',wreg,level=label,notes='Time-mean scalar speed; WorldClim land reference for surface layer, with effective model-height caveat.',selected={'GLOBAL','LAND','OCEAN'})
        for key,name,unit in [('directional_consistency','Directional consistency','fraction'),('speed_stddev_mps','Within-run speed variability','m/s'),('correlated_stderr_mps','Correlated standard error','m/s')]:
            add_field('Circulation',name,unit,np.average(weather[layer][key],axis=0,weights=QDAYS),None,'Annual','NOT_AVAILABLE',wreg,level=label,notes='Internal weather-sampling statistic; monthly mean references cannot validate it.',selected={'GLOBAL','LAND','OCEAN'})
    for i,period in enumerate(QUARTERS+['Annual']):
        fields=[a[i] if i<4 else np.average(a,axis=0,weights=QDAYS) for pair in windpairs for a in pair]
        u,v,ru,rv,U,V,RU,RV=fields
        add_field('Circulation','Upper-minus-surface vector shear','m/s',np.hypot(U-u,V-v),np.hypot(RU-ru,RV-rv),period,'ERA5_2001_2020',wreg,notes='Effective two-layer shear, not an independently resolved vertical profile.')
    # Actual hydrology transport and heating. Mixing applied separately is not in advective fluxes.
    py,px,process=load_field_table(run/'maps/climate_process_fields.csv');process=process[0]
    preg=regions(py,px,False,np.clip(remap(land,py,px,False),0,1))
    tcwv=quarter_means(reference('era5_tcwv_monthly.uwclim',py,px,False))
    ascent=quarter_means(reference('era5_w500_ascent_monthly.uwclim',py,px,False))
    small={'GLOBAL','LAND','OCEAN','N_ATLANTIC','N_PACIFIC','INDIAN','SOUTHERN_OCEAN','EQUATOR_LAND','N_MIDLAT_LAND','S_MIDLAT_LAND'}
    for key,name,unit,domain,rr,refid in [
        ('mean_column_water_mm','Total column water vapour','mm','Moisture',tcwv,'ERA5_2001_2020'),
        ('ascent_hpa_day','Midlevel ascent','hPa/day','Ascent',ascent,'ERA5_2001_2020'),
        ('transport_convergence_mm_day','Advective moisture convergence','mm/day','Transport',None,'NOT_AVAILABLE'),
        ('lower_radiative_wm2','Lower radiative heating','W/m2','Energy',None,'NOT_AVAILABLE'),
        ('upper_radiative_wm2','Upper radiative heating','W/m2','Energy',None,'NOT_AVAILABLE'),
        ('lower_latent_wm2','Lower latent heating','W/m2','Energy',None,'NOT_AVAILABLE'),
        ('upper_latent_wm2','Upper latent heating','W/m2','Energy',None,'NOT_AVAILABLE'),
        ('sensible_wm2','Surface sensible heat','W/m2','Energy',None,'NOT_AVAILABLE'),
        ('surface_net_heating_wm2','Net surface heat exchange','W/m2','Energy',None,'NOT_AVAILABLE'),
        ('column_energy_residual_wm2','Column energy residual','W/m2','Energy',None,'INTERNAL_CHECK')]:
        for i,period in enumerate(QUARTERS+['Annual']):
            ss=process[key][i] if i<4 else np.average(process[key],axis=0,weights=QDAYS)
            rf=(rr[i] if i<4 else np.average(rr,axis=0,weights=QDAYS)) if rr is not None else None
            add_field(domain,name,unit,ss,rf,period,refid,preg,selected=small if domain=='Energy' else None,
                notes='Quarter-mean actual transport diagnostic. Convergence excludes separately applied mixing.' if domain=='Transport' else '')
            if domain=='Ascent':
                for sign,label in [(1,'Ascending area fraction'),(-1,'Descending area fraction')]:
                    add_field('Ascent',label,'fraction',(ss*sign>0).astype(float),(rf*sign>0).astype(float),period,refid,preg,notes='Area fraction with strictly positive ascent / descent; zero threshold in hPa/day.')
    for layer in ['lower','upper','column']:
        if layer=='column':fu=process['lower_flux_east_kg_m_s']+process['upper_flux_east_kg_m_s'];fv=-process['lower_flux_south_kg_m_s']-process['upper_flux_south_kg_m_s']
        else:fu=process[layer+'_flux_east_kg_m_s'];fv=-process[layer+'_flux_south_kg_m_s']
        for direction,a in [('Eastward',fu),('Northward',fv),('Magnitude',np.hypot(fu,fv))]:
            for i,period in enumerate(QUARTERS+['Annual']):
                ss=a[i] if i<4 else np.average(a,axis=0,weights=QDAYS)
                add_field('Transport',direction+' moisture flux','kg/m/s',ss,None,period,'NOT_AVAILABLE',preg,level=layer,
                    selected=None if layer=='column' else small,notes='Actual donor plus corrective face transport; regional boundary inflow/outflow still requires face-boundary integration.')
    # Final terrain can change after climate (fjords/rivers), so it cannot recover
    # the ocean solver's earlier mask. Retain only demonstrably wet solver cells.
    oy,ox,ocean=load_field_table(run/'maps/ocean_fields.csv');ocean=ocean[0]
    wet=identified_ocean_mask(ocean);oreg=regions(oy,ox,False,1-wet.astype(float),'identified-ocean')
    for key,name,unit,domain,scale in [('sst_c','Liquid SST','degC','Ocean',1),('east_current_mps','Eastward current','m/s','Ocean',1),('south_current_mps','Northward current','m/s','Ocean',-1),('ekman_upwelling_mps','Ekman upwelling','m/day','Ocean',86400),('ice_thickness_m','Sea-ice thickness','m','Snow and ice',1),('surface_skin_temperature_c','Ocean or ice surface temperature','degC','Snow and ice',1)]:
        for i,period in enumerate(QUARTERS+['Annual']):
            ss=(ocean[key][i] if i<4 else np.average(ocean[key],axis=0,weights=QDAYS))*scale
            ss=np.where(wet,ss,np.nan)
            start=len(RECORDS)
            add_field(domain,name,unit,ss,None,period,'NOT_AVAILABLE',oreg,selected={x for x in small if 'LAND' not in x}|{'ARCTIC_OCEAN','S_POLAR','N_POLAR'},notes='Partial coverage: cells identified as wet by nonzero solver currents or ice. Native wet mask missing; completely stationary ice-free cells may be omitted. Final terrain is not the earlier solver mask. Liquid SST differs from ice skin.')
            for record in RECORDS[start:]:record['status']+='; partial wet-cell mask'
    # The stored ocean evaporation field is an earlier bulk-aerodynamic diagnostic
    # in mm/day, NOT the actual later hydrology evaporation. Never use it in P-E.
    evaporation_key='ocean_evaporation_proxy_mm_day' if 'ocean_evaporation_proxy_mm_day' in cell else 'evaporation_stored_mm_month'
    for i,period in enumerate(QUARTERS+['Annual']):
        ss=cell[evaporation_key][i] if i<4 else np.average(cell[evaporation_key],axis=0,weights=QDAYS)
        add_field('Surface water','Ocean bulk evaporation proxy','mm/day',np.where(land<.5,ss,np.nan),None,period,'NOT_AVAILABLE',reg,selected={'GLOBAL','OCEAN','N_ATLANTIC','N_PACIFIC','INDIAN'},notes='Rounded earlier ocean bulk-aerodynamic diagnostic, not actual hydrology evaporation. Actual totals are in Water budget.')
    # Global classification only, on a common valid land mask. Aw/As merged for scoring.
    sim=cell['koppen_id'][0].astype(int); ref=cell['reference_koppen_id'][0].astype(int)
    valid=(land>0)&(sim>=1)&(sim<=31)&(ref>=1)&(ref<=31)
    sim=np.where(sim==4,3,sim);ref=np.where(ref==4,3,ref)
    codes=['Af','Am','Aw/As','BWh','BWk','BSh','BSk','Csa','Csb','Csc','Cwa','Cwb','Cwc','Cfa','Cfb','Cfc','Dsa','Dsb','Dsc','Dsd','Dwa','Dwb','Dwc','Dwd','Dfa','Dfb','Dfc','Dfd','ET','EF']
    ids=[1,2,3]+list(range(5,32))
    confusion=np.array([[float(area[valid&(sim==i)&(ref==j)].sum()) for j in ids] for i in ids])
    assert abs(confusion.sum()-area[valid].sum())<1e-6
    koppen=dict(codes=codes,ids=ids,confusion_km2=confusion.tolist(),common_area_km2=float(area[valid].sum()),
        model_land_area_km2=float(area[land>0].sum()),valid_cells=int(valid.sum()),
        raw_counts=[int(((cell['koppen_id'][0]==i)&(land>0)).sum()) for i in range(1,32)],
        raw_areas=[float(area[(cell['koppen_id'][0]==i)&(land>0)].sum()) for i in range(1,32)])
    diagnostics={}
    for name in ['climate_coupling','climate_hydrology_spinup','climate_energy_budget','climate_water_budget_area_weighted','climate_precipitation_processes','climate_precipitation_distribution','climate_condensation_activity','climate_atmosphere_budget','climate_circulation_precision']:
        p=run/(name+'.csv')
        if p.exists():diagnostics[name]=read_csv(p)
    registry=json.loads((ROOT/'data/climate/benchmark_runs.json').read_text())
    runs=registry['runs']; metadata=next(r for r in runs if int(r['id'])==args.run_id)
    refs=[]
    for p in sorted(REF.glob('*.uwclim')):
        refs.append(dict(file=p.name,path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    logpath=args.output/'benchmark.log'
    detailed_coupling=[]
    if logpath.exists():
        for line in logpath.read_text(encoding='utf-8',errors='replace').splitlines():
            if line.startswith('Climate coupling iteration='):
                detailed_coupling.append({k:float(v) for k,v in re.findall(r'(\w+)=([-+\deE.]+)',line)})
    payload=clean(dict(run_id=args.run_id,resolution=[nx,ny],wind_grid=[wx,wy],moisture_grid=[px,py],metadata=metadata,
        metrics=RECORDS,vectors=VECTORS,features=FEATURES,koppen=koppen,diagnostics=diagnostics,
        detailed_coupling=detailed_coupling,ar6=AR6,regions=list(masks),reference_files=refs,
        hydrology_weight_sum=float(2*py*px/np.pi),
        method='Spherical cell-area weights; conservative area-average references; fractional AR6 region overlap; common finite-value mask.'))
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'validation_data.json').write_text(json.dumps(payload,allow_nan=False),encoding='utf-8')
    for name,items in [('validation_metrics',RECORDS),('validation_vectors',VECTORS)]:
        with (run/(name+'.csv')).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(clean(items))
    # Flat, worksheet-ordered inputs for pasting subsequent runs into Excel.
    paired_headers=['Input ID','Run ID','Domain','Region','Period','Variable','Layer','Units','Reference ID','Sim mean','Ref mean','Sim variance','Ref variance','Covariance','MAE','MSE','Sim P10','Ref P10','Sim P50','Ref P50','Sim P90','Ref P90','Sim P95','Ref P95','Area km2','Valid cells','Effective cells','Positive quantity','Status','Definition / source']
    paired_rows=[]
    for index,r in enumerate(RECORDS,1):
        quantile_pairs=[v for pair in zip(r['sim_quantiles'],r['ref_quantiles']) for v in pair]
        paired_rows.append([index,args.run_id,*[r[k] for k in ['domain','region','period','variable','level','units','reference_id','sim_mean','ref_mean','sim_variance','ref_variance','covariance','mae','mse']],*quantile_pairs,r['area_km2'],r['cells'],r['effective_cells'],int(r['positive']),r['status'],r['notes']])
    vector_headers=['Run ID','Region','Period','Level','Component','Sim speed','Ref speed','Vector MSE','Alignment','Angle error deg','u correlation','v correlation','Sim u','Ref u','Sim v','Ref v','Sim RMS squared','Ref RMS squared','Area km2','Valid cells','Eff. cells','Reference ID','Status']
    vector_keys=['region','period','level','mode','sim_speed','ref_speed','vector_mse','vector_alignment','angle_error','u_correlation','v_correlation','sim_u','ref_u','sim_v','ref_v','sim_rms2','ref_rms2','area_km2','cells','effective_cells','reference_id','status']
    for name,headers,rows in [('excel_paired_data',paired_headers,paired_rows),
                             ('excel_vector_data',vector_headers,[[args.run_id,*[r[k] for k in vector_keys]] for r in VECTORS]),
                             ('excel_koppen_matrix',['Simulation / reference',*codes],[[code,*row] for code,row in zip(codes,confusion.tolist())])]:
        with (run/(name+'.csv')).open('w',newline='',encoding='utf-8-sig') as f:
            writer=csv.writer(f);writer.writerow(headers);writer.writerows(clean(rows))
    print(json.dumps(dict(run_id=args.run_id,scalar_rows=len(RECORDS),vector_rows=len(VECTORS),regions=len(masks),common_koppen_area_km2=koppen['common_area_km2'])))

AR6=[]
if __name__=='__main__':main()
