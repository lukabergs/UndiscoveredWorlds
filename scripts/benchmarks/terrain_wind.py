"""Experimental diagnostic terrain-wind downscaling; SI units, east/north vectors.

No moisture/heat feedback or mass-consistency solve is implied. Terrain response
at coarse resolution is subtracted as an approximate scale separation; this is
not an exact inversion of the global circulation model's orographic forcing.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class TerrainWindConfig:
    shelter: float = .35
    exposure: float = .20
    turning: float = .30
    channel: float = .10
    stability_per_second: float = .012
    wake_length_m: float = 100000.
    terrain_scale_m: float = 25000.
    maximum_relative_correction: float = .5
    distances_m: tuple = (12500.,25000.,50000.,100000.,200000.)


def coordinates(shape):
    rows,columns = shape
    return (np.deg2rad(90-(np.arange(rows,dtype=np.float32)+.5)*180/rows)[:,None],
            np.deg2rad(-180+(np.arange(columns,dtype=np.float32)+.5)*360/columns)[None,:])


def destination(lat,lon,east,north,distance,radius):
    """Great-circle offset in the local tangent basis, including polar crossings."""
    delta = distance/radius
    sine = np.sin(lat)*np.cos(delta)+np.cos(lat)*np.sin(delta)*north
    target_lat = np.arcsin(np.clip(sine,-1,1))
    target_lon = lon+np.arctan2(east*np.sin(delta)*np.cos(lat),
                              np.cos(delta)-np.sin(lat)*np.sin(target_lat))
    return target_lat,target_lon


def sample(field,lat,lon):
    """Periodic longitude; cell-centre interpolation. Polar caps use last row."""
    rows,columns = field.shape
    x = np.remainder((lon+np.pi)*(columns/(2*np.pi))-.5,columns)
    y = np.clip((np.pi/2-lat)*(rows/np.pi)-.5,0,rows-1)
    x0,y0 = np.floor(x).astype(int),np.floor(y).astype(int)
    fx,fy = x-x0,y-y0
    x1,y1 = (x0+1)%columns,np.minimum(y0+1,rows-1)
    return ((1-fy)*((1-fx)*field[y0,x0]+fx*field[y0,x1])+
            fy*((1-fx)*field[y1,x0]+fx*field[y1,x1])).astype(np.float32)


class TerrainWind:
    def __init__(self,fine_height,coarse_height,config=TerrainWindConfig(),radius=6371000.):
        if radius<=0 or config.terrain_scale_m<=0 or config.wake_length_m<=0:
            raise ValueError("Physical lengths must be positive")
        self.fine = np.maximum(np.asarray(fine_height,dtype=np.float32),0)
        self.coarse = np.maximum(np.asarray(coarse_height,dtype=np.float32),0)
        if not np.isfinite(self.fine).all() or not np.isfinite(self.coarse).all():
            raise ValueError("Terrain must be finite")
        self.config,self.radius = config,radius
        self.lat,self.lon = coordinates(self.fine.shape)
        self.heights = [self.fine,self.coarse if self.coarse.shape==self.fine.shape else sample(self.coarse,self.lat,self.lon)]
        offsets = [destination(self.lat,self.lon,e,n,config.terrain_scale_m,radius)
                   for e,n in ((1,0),(-1,0),(0,1),(0,-1))]
        self.gradients = []
        for field in (self.fine,self.coarse):
            e,w,n,s = [sample(field,*point) for point in offsets]
            self.gradients.append(((e-w)/(2*config.terrain_scale_m),(n-s)/(2*config.terrain_scale_m)))

    def apply(self,east,north):
        u,v = np.asarray(east,dtype=np.float32),np.asarray(north,dtype=np.float32)
        if u.shape!=self.fine.shape or v.shape!=u.shape or not np.isfinite(u).all() or not np.isfinite(v).all():
            raise ValueError("Winds must be finite and match terrain shape")
        c = self.config
        speed = np.hypot(u,v)
        e = np.divide(u,speed,out=np.zeros_like(u),where=speed>0)
        n = np.divide(v,speed,out=np.zeros_like(v),where=speed>0)
        wake = [np.zeros_like(u),np.zeros_like(u)]
        ahead = [np.zeros_like(u),np.zeros_like(u)]
        for distance in c.distances_m:
            back = destination(self.lat,self.lon,-e,-n,distance,self.radius)
            front = destination(self.lat,self.lon,e,n,distance,self.radius)
            for j,field in enumerate((self.fine,self.coarse)):
                wake[j] = np.maximum(wake[j],np.maximum(sample(field,*back)-self.heights[j],0)*np.exp(-distance/c.wake_length_m))
                if distance<=2*c.terrain_scale_m:
                    ahead[j] = np.maximum(ahead[j],np.maximum(sample(field,*front)-self.heights[j],0))
        points = [destination(self.lat,self.lon,ee,nn,c.terrain_scale_m,self.radius)
                  for ee,nn in ((e,n),(-e,-n),(-n,e),(n,-e))]
        corrections = []
        for j,field in enumerate((self.fine,self.coarse)):
            forward,back,left,right = [sample(field,*point) for point in points]
            h = self.heights[j]
            # Convex terrain along flow exposes a ridge; higher ground on both
            # sides indicates a channel. These are bounded diagnostic closures.
            ridge = np.tanh(np.maximum(h-.5*(forward+back),0)/(c.terrain_scale_m*.04))
            gap = np.tanh(np.maximum(np.minimum(left,right)-h,0)/(c.terrain_scale_m*.04))
            denom = np.maximum(speed,2)
            sheltered = -np.expm1(-np.square(c.stability_per_second*wake[j]/denom))
            blocked = -np.expm1(-np.square(c.stability_per_second*ahead[j]/denom))
            gx,gy = self.gradients[j]
            norm = np.hypot(gx,gy)
            gx = np.divide(gx,norm,out=np.zeros_like(gx),where=norm>1e-7)
            gy = np.divide(gy,norm,out=np.zeros_like(gy),where=norm>1e-7)
            into = np.maximum(e*gx+n*gy,0)
            scaling = c.exposure*ridge+c.channel*gap-c.shelter*sheltered
            turning = c.turning*blocked*into
            corrections.append((scaling*e-turning*gx,scaling*n-turning*gy))
        du,dv = [corrections[0][i]-corrections[1][i] for i in (0,1)]
        length = np.hypot(du,dv)
        bound = np.minimum(1,c.maximum_relative_correction/np.maximum(length,1e-20))
        du,dv = speed*du*bound,speed*dv*bound
        return (u+du).astype(np.float32),(v+dv).astype(np.float32)
