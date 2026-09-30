"""Deterministic, matched ocean LIC and physical-time particle comparisons.

Inputs are cell-centred east/north velocities in m/s on W x W/2 grids.
Rendering interpolates the available fields; it does not resolve new dynamics.
Both panels use the same mask, seed, speed scale and integration duration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'refs'))
from reference_map_rendering import grid_velocity, lic_luminance, sample, scalar_rgb


def resample(east, north, wet, width):
    h, w = east.shape
    if north.shape != (h,w) or wet.shape != (h,w) or w != h*2:
        raise ValueError('Expected matching cell-centred W x W/2 fields')
    if width < 8 or width % 4:
        raise ValueError('Render width must be a positive multiple of four')
    yy,xx=np.indices((width//2,width),dtype=float)
    x=(xx+.5)*w/width-.5; y=(yy+.5)*h/(width//2)-.5
    # Conservative coast stopping: do not interpolate across any contributing
    # dry/absent corner. Periodic longitude; nearest edge row at the poles.
    e=sample(np.where(wet,east,np.nan),x,y)
    n=sample(np.where(wet,north,np.nan),x,y)
    valid=np.isfinite(e)&np.isfinite(n)
    return e,n,valid


def advance(east,north,x,y,remaining):
    """One midpoint step, capped at 0.45 pixel and six hours; coast/pole stop."""
    h,w=east.shape
    dx,dy=grid_velocity(east,north,x,y)
    rate=np.hypot(dx,dy)
    valid=np.isfinite(rate)&(rate>1e-15)&(remaining>0)
    dt=np.minimum(remaining,np.minimum(21600.,.45/np.maximum(np.nan_to_num(rate),1e-15)))
    dx=np.nan_to_num(dx);dy=np.nan_to_num(dy)
    mx=(x+dx*dt*.5)%w;my=y+dy*dt*.5
    mdx,mdy=grid_velocity(east,north,mx,my)
    nx=(x+np.nan_to_num(mdx)*dt)%w;ny=y+np.nan_to_num(mdy)*dt
    # Check start, midpoint and endpoint so no trail can enter masked water.
    ex,ey=grid_velocity(east,north,nx,ny)
    valid &= np.isfinite(mdx)&np.isfinite(mdy)&np.isfinite(ex)&np.isfinite(ey)
    valid &= (my>=0)&(my<=h-1)&(ny>=0)&(ny<=h-1)
    return nx,ny,np.maximum(0.,remaining-dt),valid


def particles(east,north,seed,days=45,count=20000):
    if days<=0 or days>180: raise ValueError('Particle duration must be in (0,180] days')
    h,w=east.shape
    rng=np.random.Generator(np.random.MT19937(seed))
    x=rng.uniform(0,w,count)
    y=(90-np.degrees(np.arcsin(rng.uniform(-.995,.995,count))))*h/180-.5
    y=np.clip(y,0,h-1)
    remaining=np.full(count,days*86400.)
    hits=np.zeros((h,w),float);heads=[];iterations=0
    while x.size:
        nx,ny,nr,valid=advance(east,north,x,y,remaining)
        # Deposits brighten with time; terminal points identify flow direction.
        if valid.any():
            alpha=.2+.8*(1-nr[valid]/(days*86400.))
            np.add.at(hits,(np.rint(ny[valid]).astype(int),np.rint(nx[valid]).astype(int)%w),alpha)
        finished=valid&(nr<=0)
        if finished.any():heads.extend(zip(nx[finished].tolist(),ny[finished].tolist()))
        keep=valid&(nr>0);x,y,remaining=nx[keep],ny[keep],nr[keep]
        iterations+=1
        if iterations>20000: raise RuntimeError('Particle integration exceeded safety limit')
    intensity=1-np.exp(-.5*hits)
    intensity[~np.isfinite(east)|~np.isfinite(north)]=np.nan
    return intensity,heads,iterations


def background(east,north,valid,speed_limit):
    rgb=scalar_rgb(np.hypot(east,north),speed_limit).astype(float)
    rgb[~valid]=(38,39,40)
    return rgb


def arrows(image,east,north,valid):
    draw=ImageDraw.Draw(image);h,w=east.shape;spacing=max(12,w//48)
    for y in range(spacing//2,h,spacing):
        for x in range(spacing//2,w,spacing):
            if not valid[y,x]:continue
            dx,dy=grid_velocity(east,north,np.array(float(x)),np.array(float(y)))
            length=float(np.hypot(dx,dy))
            if not np.isfinite(length) or length<=1e-15:continue
            dx=float(dx)/length;dy=float(dy)/length
            size=min(9.,spacing*.35);tip=(x+dx*size/2,y+dy*size/2)
            draw.line((x-dx*size/2,y-dy*size/2,*tip),fill=(245,232,170))
            for side in (-1,1):
                draw.line((*tip,tip[0]-dx*size*.4+side*dy*size*.22,
                           tip[1]-dy*size*.4-side*dx*size*.22),fill=(245,232,170))


def render_pair(east,north,reference_east,reference_north,wet,output_prefix,
                label='Simulation',width=2048,seed=20260906,speed_limit=1.,days=45):
    if speed_limit<=0:raise ValueError('Speed scale must be positive')
    common=wet & np.isfinite(east)&np.isfinite(north)&np.isfinite(reference_east)&np.isfinite(reference_north)
    panels={};stats=[]
    for title,e,n in ((label,east,north),('GLORYS reference',reference_east,reference_north)):
        e,n,valid=resample(e,n,common,width)
        rgb=background(e,n,valid,speed_limit)
        light=lic_luminance(e,n,seed)
        lic=rgb*(.65+.7*np.nan_to_num(light)[...,None]);lic[~valid]=(38,39,40)
        intensity,heads,iterations=particles(e,n,seed,days)
        particle=rgb*.7+(255-rgb*.7)*np.nan_to_num(intensity)[...,None]
        particle[~valid]=(38,39,40)
        for kind,array in (('lic',lic),('particles',particle)):
            img=Image.fromarray(np.clip(array,0,255).astype(np.uint8));arrows(img,e,n,valid)
            if kind=='particles':
                draw=ImageDraw.Draw(img)
                for x,y in heads:draw.ellipse((x-1,y-1,x+1,y+1),fill=(255,246,184))
            panels[(title,kind)]=img
        stats.append(dict(label=title,valid_pixels=int(valid.sum()),particle_iterations=iterations,
                          terminal_particles=len(heads),speed_above_scale_pixels=int((np.hypot(e,n)>speed_limit).sum())))
    prefix=Path(output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    files=[]
    for kind in ('lic','particles'):
        canvas=Image.new('RGB',(width*2,width//2+88),(17,25,35));draw=ImageDraw.Draw(canvas)
        for i,title in enumerate((label,'GLORYS reference')):
            draw.text((i*width+12,12),f'{title} | ocean {kind} | native {east.shape[1]} x {east.shape[0]}',fill='white')
            canvas.paste(panels[(title,kind)],(i*width,36))
        y=width//2+44
        draw.text((12,y),f'Identical scale 0 to {speed_limit:g} m/s; seed {seed}; gray = common land/missing/coast mask. Arrows show direction.',fill='white')
        draw.text((12,y+17),f'Particles: {days:g} physical days, fixed seasonal flow, coast stopping. LIC: normalized streamlines; texture length is not speed.',fill='white')
        ramp=np.linspace(0,speed_limit,256)[None,:].repeat(10,axis=0)
        canvas.paste(Image.fromarray(scalar_rgb(ramp,speed_limit)),(width,y+6))
        draw.text((width+265,y+4),f'0 -> {speed_limit:g} m/s (higher values saturate)',fill='white')
        path=Path(str(prefix)+f'-{kind}.png');canvas.save(path);files.append(path)
    receipt=dict(native_grid=[east.shape[1],east.shape[0]],render_panel=[width,width//2],seed=seed,
                 particle_days=days,scale_mps=[0,speed_limit],common_wet_native_cells=int(common.sum()),panels=stats,
                 caveat='GLORYS near-surface current versus model effective mixed-layer transport velocity; interpolation adds no resolved structure.',
                 output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    Path(str(prefix)+'-render.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path,help='NPZ: east,north,reference_east,reference_north,wet')
    p.add_argument('output_prefix',type=Path);p.add_argument('--label',default='Simulation')
    p.add_argument('--width',type=int,default=2048);p.add_argument('--days',type=float,default=45)
    args=p.parse_args();data=np.load(args.input)
    result=render_pair(*(data[k] for k in ('east','north','reference_east','reference_north','wet')),
                       args.output_prefix,label=args.label,width=args.width,days=args.days)
    print(json.dumps(result))


if __name__=='__main__':main()
