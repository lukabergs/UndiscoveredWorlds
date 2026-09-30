"""Diagnostic obstacle diversion on a fixed background wind, using SI units.

The baseline terrain closure is unchanged. Look-ahead steering responds to a
height difference between the two flanks of an upcoming barrier. It is a
bounded steady heuristic, not a mass/momentum solver or a turbulence model.
"""
from dataclasses import dataclass
import numpy as np
from terrain_wind import TerrainWind, TerrainWindConfig, destination, sample


@dataclass(frozen=True)
class ObstacleConfig:
    terrain: TerrainWindConfig = TerrainWindConfig(
        shelter=.85, exposure=.35, turning=.9, channel=.25, maximum_relative_correction=.9)
    diversion_degrees: float = 55.
    lookahead_m: float = 75000.
    flank_distance_m: float = 40000.
    lateral_height_scale_m: float = 500.
    coarse_response_fraction: float = 1.
    maximum_turn_degrees: float = 75.
    minimum_speed_ratio: float = .1
    maximum_speed_ratio: float = 1.6


class ObstacleWind:
    def __init__(self, fine_height, coarse_height, config=ObstacleConfig(), radius=6371000.):
        c=config
        if not (0 <= c.diversion_degrees <= 90 and 0 < c.maximum_turn_degrees < 90 and
                0 < c.minimum_speed_ratio <= 1 <= c.maximum_speed_ratio and
                min(c.lookahead_m,c.flank_distance_m,c.lateral_height_scale_m) > 0 and
                0 <= c.coarse_response_fraction <= 1):
            raise ValueError('Invalid obstacle geometry or response bounds')
        self.config=c
        self.base=TerrainWind(fine_height,coarse_height,c.terrain,radius)

    def apply(self,east,north):
        u,v=self.base.apply(east,north)
        b,c=self.base,self.config
        original=np.hypot(east,north)
        # Probe along the imposed background, not successively along each
        # corrected result: repeated calls must not compound artificial drag.
        e=np.divide(east,original,out=np.zeros_like(u),where=original>0)
        n=np.divide(north,original,out=np.zeros_like(v),where=original>0)
        middle=destination(b.lat,b.lon,e,n,c.lookahead_m,b.radius)
        distance=np.hypot(c.lookahead_m,c.flank_distance_m)
        along,across=c.lookahead_m/distance,c.flank_distance_m/distance
        left=destination(b.lat,b.lon,along*e-across*n,along*n+across*e,distance,b.radius)
        right=destination(b.lat,b.lon,along*e+across*n,along*n-across*e,distance,b.radius)
        responses=[]
        for field,height in zip((b.fine,b.coarse),b.heights):
            front=np.maximum(sample(field,*middle)-height,0)
            lh=np.maximum(sample(field,*left)-height,0)
            rh=np.maximum(sample(field,*right)-height,0)
            blocked=-np.expm1(-np.square(c.terrain.stability_per_second*front/np.maximum(original,2)))
            # Positive rotation turns left, toward the lower flank. A centred
            # symmetric barrier leaves the centreline direction unchanged.
            responses.append(blocked*np.tanh((rh-lh)/c.lateral_height_scale_m))
        steer=np.deg2rad(c.diversion_degrees)*(responses[0]-c.coarse_response_fraction*responses[1])
        speed=np.clip(np.hypot(u,v),c.minimum_speed_ratio*original,c.maximum_speed_ratio*original)
        angle=np.arctan2(e*v-n*u,e*u+n*v)+steer
        angle=np.clip(angle,-np.deg2rad(c.maximum_turn_degrees),np.deg2rad(c.maximum_turn_degrees))
        cosine,sine=np.cos(angle),np.sin(angle)
        return (speed*(e*cosine-n*sine)).astype(np.float32),(speed*(n*cosine+e*sine)).astype(np.float32)
