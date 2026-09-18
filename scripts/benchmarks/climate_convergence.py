"""Wind-only tropical convergence indicator; not a complete convective ITCZ diagnosis."""
import numpy as np

def surface_convergence(east, north, radius_m=6371000.0):
    """Closed spherical finite-volume -div(u), s^-1, cell-centred 2:1 grid."""
    east=np.asarray(east,dtype=float);north=np.asarray(north,dtype=float)
    if east.ndim!=2 or north.shape!=east.shape or east.shape[1]!=2*east.shape[0] or east.shape[0]<2:
        raise ValueError('Expected matching cell-centred W × W/2 vector grids')
    if not np.isfinite(east).all() or not np.isfinite(north).all() or not np.isfinite(radius_m) or radius_m<=0:
        raise ValueError('Vectors and positive radius must be finite')
    ny,nx=east.shape;dphi=np.pi/ny;dlambda=2*np.pi/nx
    edges=np.linspace(np.pi/2,-np.pi/2,ny+1)
    area=radius_m**2*dlambda*(np.sin(edges[:-1])-np.sin(edges[1:]))
    zonal=(np.roll(east,-1,1)-np.roll(east,1,1))*.5*radius_m*dphi
    faces=np.zeros((ny+1,nx));faces[1:-1]=(north[:-1]+north[1:])*.5*radius_m*dlambda*np.cos(edges[1:-1,None])
    return -(zonal+faces[:-1]-faces[1:])/area[:,None]

def smooth_vectors(east,north,passes=2):
    """Fixed binomial smoothing before diagnosis, cyclic longitude, reflected latitude."""
    result=[]
    for values in (east,north):
        a=np.asarray(values,dtype=float).copy()
        for _ in range(passes):
            a=(np.roll(a,1,1)+2*a+np.roll(a,-1,1))/4
            p=np.pad(a,((1,1),(0,0)),mode='edge');a=(p[:-2]+2*p[1:-1]+p[2:])/4
        result.append(a)
    return result

def tropical_band(convergence_per_second,latitude_limit=30.0,minimum_per_day=.10,peak_fraction=.5):
    """Keep all sufficiently strong tropical branches; weak columns remain blank."""
    c=np.asarray(convergence_per_second,dtype=float)*86400
    if c.ndim!=2 or not np.isfinite(c).all() or not 0<latitude_limit<90 or minimum_per_day<=0 or not 0<peak_fraction<=1:
        raise ValueError('Invalid convergence-band input or threshold')
    lat=90-(np.arange(c.shape[0])+.5)*180/c.shape[0];tropics=np.abs(lat[:,None])<=latitude_limit
    peaks=np.maximum(0,np.max(np.where(tropics,c,0),axis=0))
    threshold=np.maximum(minimum_per_day,peak_fraction*peaks)
    return tropics&(c>=threshold),threshold
