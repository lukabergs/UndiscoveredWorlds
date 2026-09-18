"""Advected, laterally mixed terrain shelter; optional diagnostic wind closure.

Height memory (m) is carried downwind and relaxed over the legacy wake length.
It is an obstruction proxy, not conserved momentum. Positive interpolation and
mixing keep it bounded by the highest source terrain. Old slope/turning and
look-ahead diversion remain unchanged; no climate feedback is added.
"""
from dataclasses import dataclass
import numpy as np
from terrain_wind import TerrainWind, destination, sample
from terrain_wind_obstacles import ObstacleWind, ObstacleConfig


@dataclass(frozen=True)
class WakeConfig:
    step_m: float = 20000.
    reach_m: float = 240000.
    crosswind_mixing_m: float = 1500.


class Stencil:
    def __init__(self,shape,lat,lon):
        rows,columns=shape
        x=np.remainder((lon+np.pi)*(columns/(2*np.pi))-.5,columns)
        y=np.clip((np.pi/2-lat)*(rows/np.pi)-.5,0,rows-1)
        x,y=np.broadcast_arrays(x,y)
        x0=np.floor(x).astype(np.int32);y0=np.floor(y).astype(np.int32)
        fx,fy=(x-x0).astype(np.float32),(y-y0).astype(np.float32)
        x1=(x0+1)%columns;y1=np.minimum(y0+1,rows-1)
        self.indices=[y0*columns+x0,y0*columns+x1,y1*columns+x0,y1*columns+x1]
        self.weights=[(1-fx)*(1-fy),fx*(1-fy),(1-fx)*fy,fx*fy]

    def apply(self,field):
        flat=field.ravel();out=np.zeros_like(self.weights[0])
        for indices,weight in zip(self.indices,self.weights):out+=flat[indices]*weight
        return out


class SpreadTerrainWind(TerrainWind):
    def __init__(self,fine,coarse,config,wake=WakeConfig(),radius=6371000.):
        if (not np.isfinite([wake.step_m,wake.reach_m,wake.crosswind_mixing_m]).all()
                or wake.step_m<=0 or wake.reach_m<wake.step_m or wake.crosswind_mixing_m<0):
            raise ValueError('Invalid wake lengths or mixing coefficient')
        super().__init__(fine,coarse,config,radius)
        self.wake_config=wake

    def shelter(self,e,n):
        c=self.wake_config;steps=int(np.ceil(c.reach_m/c.step_m));step=c.reach_m/steps
        lateral=2*np.sqrt(c.crosswind_mixing_m*step)
        stencils=[]
        for side,weight in ((0,.5),(-1,.25),(1,.25)) if lateral else ((0,1.),):
            offset=side*lateral;distance=np.hypot(step,offset)
            point=destination(self.lat,self.lon,(-step*e-offset*n)/distance,
                (-step*n+offset*e)/distance,distance,self.radius)
            stencils.append((Stencil(self.fine.shape,*point),weight))
        decay=np.exp(-step/self.config.wake_length_m)
        memory=[h.copy() for h in self.heights]
        for _ in range(steps):
            for j,h in enumerate(self.heights):
                carried=np.zeros_like(h)
                for stencil,weight in stencils:carried+=weight*stencil.apply(memory[j])
                memory[j]=np.maximum(h,decay*carried)
        return [np.maximum(q-h,0) for q,h in zip(memory,self.heights)]

    def apply(self,east,north):
        # Frozen legacy response algebra, with only the shelter calculation
        # replaced. Keeping the old module intact preserves previous receipts.
        u,v=np.asarray(east,dtype=np.float32),np.asarray(north,dtype=np.float32)
        if u.shape!=self.fine.shape or v.shape!=u.shape or not np.isfinite(u).all() or not np.isfinite(v).all():
            raise ValueError('Winds must be finite and match terrain shape')
        c=self.config;speed=np.hypot(u,v)
        e=np.divide(u,speed,out=np.zeros_like(u),where=speed>0)
        n=np.divide(v,speed,out=np.zeros_like(v),where=speed>0)
        wake=self.shelter(e,n);ahead=[np.zeros_like(u),np.zeros_like(u)]
        for distance in c.distances_m:
            if distance>2*c.terrain_scale_m:continue
            front=destination(self.lat,self.lon,e,n,distance,self.radius)
            for j,field in enumerate((self.fine,self.coarse)):
                ahead[j]=np.maximum(ahead[j],np.maximum(sample(field,*front)-self.heights[j],0))
        points=[destination(self.lat,self.lon,ee,nn,c.terrain_scale_m,self.radius)
            for ee,nn in ((e,n),(-e,-n),(-n,e),(n,-e))]
        corrections=[]
        for j,field in enumerate((self.fine,self.coarse)):
            forward,back,left,right=[sample(field,*point) for point in points]
            h=self.heights[j]
            ridge=np.tanh(np.maximum(h-.5*(forward+back),0)/(c.terrain_scale_m*.04))
            gap=np.tanh(np.maximum(np.minimum(left,right)-h,0)/(c.terrain_scale_m*.04))
            sheltered=-np.expm1(-np.square(c.stability_per_second*wake[j]/np.maximum(speed,2)))
            blocked=-np.expm1(-np.square(c.stability_per_second*ahead[j]/np.maximum(speed,2)))
            gx,gy=self.gradients[j];norm=np.hypot(gx,gy)
            gx=np.divide(gx,norm,out=np.zeros_like(gx),where=norm>1e-7)
            gy=np.divide(gy,norm,out=np.zeros_like(gy),where=norm>1e-7)
            into=np.maximum(e*gx+n*gy,0)
            scaling=c.exposure*ridge+c.channel*gap-c.shelter*sheltered
            turning=c.turning*blocked*into
            corrections.append((scaling*e-turning*gx,scaling*n-turning*gy))
        du,dv=[corrections[0][i]-corrections[1][i] for i in (0,1)]
        bound=np.minimum(1,c.maximum_relative_correction/np.maximum(np.hypot(du,dv),1e-20))
        return (u+speed*du*bound).astype(np.float32),(v+speed*dv*bound).astype(np.float32)


class SpreadObstacleWind(ObstacleWind):
    def __init__(self,fine,coarse,config=ObstacleConfig(),wake=WakeConfig(),radius=6371000.):
        # Parent validates diversion bounds. Its final turning algorithm is reused.
        super().__init__(fine,coarse,config,radius)
        self.base=SpreadTerrainWind(fine,coarse,config.terrain,wake,radius)
