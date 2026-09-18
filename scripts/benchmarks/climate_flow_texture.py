"""Optional compiled equivalent of the reference NumPy flow texture kernels."""
import ctypes
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
_library=None
for directory in ('x64-Debug','x64-Release'):
    path=ROOT/'out/build'/directory/'Release/climate_flow_texture.dll'
    if path.exists():
        _library=ctypes.CDLL(str(path));break
if _library is not None:
    D=ctypes.POINTER(ctypes.c_double);F=ctypes.POINTER(ctypes.c_float);I=ctypes.c_int
    _library.uw_wind_lic.argtypes=[D,D,D,I,I,F];_library.uw_wind_lic.restype=None
    _library.uw_wind_particles.argtypes=[D,D,I,I,D,D,I,ctypes.c_double,D];_library.uw_wind_particles.restype=None

def pointer(a,kind=ctypes.c_double):return a.ctypes.data_as(ctypes.POINTER(kind))
def lic_luminance(east,north,seed):
    if _library is None:
        from reference_map_rendering import lic_luminance as fallback
        return fallback(east,north,seed)
    east=np.ascontiguousarray(east,dtype=float);north=np.ascontiguousarray(north,dtype=float)
    h,w=east.shape;assert north.shape==east.shape
    noise=np.random.Generator(np.random.MT19937(seed)).random(east.shape)
    result=np.empty(east.shape,dtype=np.float32)
    _library.uw_wind_lic(pointer(east),pointer(north),pointer(noise),w,h,pointer(result,ctypes.c_float))
    return result

def particle_intensity(east,north,seed):
    if _library is None:
        from reference_map_rendering import particle_intensity as fallback
        return fallback(east,north,seed)
    east=np.ascontiguousarray(east,dtype=float);north=np.ascontiguousarray(north,dtype=float)
    h,w=east.shape;assert north.shape==east.shape
    rng=np.random.Generator(np.random.MT19937(seed));count=min(20000,max(5000,w*h//10))
    x=rng.uniform(0,w,count)
    y=np.clip((90-np.rad2deg(np.arcsin(rng.uniform(-.995,.995,count))))*h/180-.5,0,h-1)
    hits=np.empty(east.shape,dtype=float)
    _library.uw_wind_particles(pointer(east),pointer(north),w,h,pointer(x),pointer(y),count,36*10800.,pointer(hits))
    result=(1-np.power(.84,hits)).astype(np.float32)
    result[~np.isfinite(east)|~np.isfinite(north)]=np.nan
    return result
