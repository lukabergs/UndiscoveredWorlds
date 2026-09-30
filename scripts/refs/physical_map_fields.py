"""Evaluate published harmonic and tectonic models; retain reference frames."""
import zipfile

import numpy as np
import pygplates
import pyshtools
from shapely.geometry import Polygon, LineString

from physical_map_core import Layer, regular_grid
from physical_map_vectors import draw_geometry, class_colours


def geodetic_latitude(geocentric,flattening):
    phi=np.deg2rad(geocentric)
    return np.rad2deg(np.arctan2(np.sin(phi),(1-flattening)**2*np.cos(phi)))


def magnetic_local_components(radial,theta,east,lat_gc,flattening):
    lat=geodetic_latitude(lat_gc,flattening)
    delta=np.deg2rad(lat-np.asarray(lat_gc))[:,None]
    north=-theta*np.cos(delta)-radial*np.sin(delta)
    down=-radial*np.cos(delta)+theta*np.sin(delta)
    return lat,north,east,down


def gravity(atlas):
    archive=atlas.one('goco06s','*.zip')
    path=atlas.unzip(archive,'GOCO06s.gfc')
    model=pyshtools.SHGravCoeffs.from_file(path,format='icgem')
    ref=pyshtools.constants.Earth.wgs84
    a,f,omega=ref.a.value,ref.f.value,ref.omega.value
    model.set_omega(omega)
    degree=max(511,model.lmax)
    grid=model.expand(a=a,f=f,lmax=degree,lmax_calc=model.lmax,normal_gravity=True,
                      normal_gravity_gm=ref.gm.value,omega=omega,extend=False)
    lat=geodetic_latitude(grid.total.lats(),f);lon=grid.total.lons()
    data=grid.total.data*1e5
    note=f'GOCO06s full degree {model.lmax}; field evaluated on WGS84 ellipsoid, centrifugal force included, WGS84 normal gravity subtracted. Geocentric sample latitudes converted to geodetic latitude.'
    atlas.emit(Layer('gravity/earth_disturbance','Earth gravity disturbance relative to WGS84','mGal','goco06s',(-150,150),'cmo.balance',
                    notes=note,method='SHTOOLS harmonic synthesis followed by spherical-area remapping'),lambda w:regular_grid(data,lat,lon,w))
    # The cubic solver produces nonfinite roots at some GOCO06s sample points.
    # The quadratic approximation is finite globally at this display scale.
    geoid=model.geoid(potref=ref.u0.value,a=a,f=f,r=model.r0,omega=omega,order=2,lmax=degree,lmax_calc=model.lmax,extend=False)
    data=geoid.geoid.data
    if not np.isfinite(data).all():raise ValueError('Geoid harmonic evaluation contains nonfinite values')
    atlas.emit(Layer('gravity/earth_geoid','GOCO06s geoid height relative to WGS84','m','goco06s',(-110,110),'cmo.balance',
                    notes=note+f' Reference potential W0={ref.u0.value} m2/s2; second-order geoid approximation for numerical stability. This does not replace the EGM2008 vertical datum of existing terrain.',
                    method='SHTOOLS second-order geoid solution; spherical-area remapping'),lambda w:regular_grid(data,lat,lon,w))
    path=atlas.one('grail-grgm660pr','*.tab');atlas.one('grail-grgm660pr','*.lbl')
    lunar=pyshtools.SHGravCoeffs.from_file(path,format='shtools',header_units='km',errors=True)
    grid=lunar.expand(a=lunar.r0,f=0,lmax=1023,lmax_calc=660,normal_gravity=False,omega=0,extend=False)
    lat=grid.rad.lats();lon=grid.rad.lons()
    data=(-grid.rad.data-lunar.gm/lunar.r0**2)*1e5
    note=f'GRAIL GRGM660PRIM, degree 660; radius {lunar.r0} m; radial inward gravity minus GM/r². No topographic/Bouguer correction or centrifugal acceleration. DE421 principal-axes frame; LOLA uses mean-Earth/polar-axis coordinates. The two lunar layers have not been rotated into a common frame.'
    atlas.emit(Layer('moon/gravity_radial_anomaly','Lunar radial gravity anomaly relative to spherical monopole','mGal','grail-grgm660pr',(-600,600),'cmo.balance',body='Moon',coordinate_frame='DE421 principal axes',notes=note,
                    method='SHTOOLS degree-660 synthesis; spherical-area mean'),lambda w:regular_grid(data,lat,lon,w))


def magnetism(atlas):
    path=atlas.one('igrf14','*.txt')
    model=pyshtools.SHMagCoeffs.from_file(path,format='igrf',r0=6371200,year=2025,units='nT')
    ref=pyshtools.constants.Earth.wgs84
    grid=model.expand(a=ref.a.value,f=ref.f.value,lmax=511,lmax_calc=13,extend=False)
    lat_gc=grid.rad.lats();lon=grid.rad.lons()
    lat,north,east,down=magnetic_local_components(grid.rad.data,grid.theta.data,grid.phi.data,lat_gc,ref.f.value)
    note='IGRF14 epoch 2025.0, sea level on WGS84; Schmidt-normalised coefficients evaluated with SHTOOLS. Vector components rotated into geodetic north/east/down. Main field only; no crustal anomalies or space-weather variation.'
    for key,title,units,limits,palette in [('intensity','Geomagnetic total intensity','nT',(20000,65000),'cmo.thermal'),
                                         ('north','Magnetic north component','nT',(-45000,45000),'cmo.balance'),
                                         ('east','Magnetic east component','nT',(-20000,20000),'cmo.balance'),
                                         ('down','Magnetic downward component','nT',(-65000,65000),'cmo.balance'),
                                         ('declination','Magnetic declination (east positive)','degrees',(-180,180),'cmo.phase'),
                                         ('inclination','Magnetic inclination (down positive)','degrees',(-90,90),'cmo.balance')]:
        def make(w):
            n,e,d=[regular_grid(v,lat,lon,w) for v in (north,east,down)]
            if key=='intensity':return np.sqrt(n*n+e*e+d*d)
            if key=='declination':return np.where(np.hypot(n,e)>1,np.rad2deg(np.arctan2(e,n)),np.nan)
            if key=='inclination':return np.rad2deg(np.arctan2(d,np.hypot(n,e)))
            return {'north':n,'east':e,'down':d}[key]
        atlas.emit(Layer('magnetism/'+key,title,units,'igrf14',limits,palette,period='2025.0',notes=note,
                        method='Degree-13 synthesis; remap vector components before deriving intensity and angles'),make)


def wrap_geometry(geometry):
    """Split spherical geometries at the dateline before planar rasterisation."""
    result=[]
    for item in pygplates.DateLineWrapper().wrap(geometry,tessellate_degrees=1):
        if isinstance(item,pygplates.DateLineWrapper.LatLonPolygon):
            ring=[(p.get_longitude(),p.get_latitude()) for p in item.get_exterior_points()]
            if len(ring)>=3:result.append(Polygon(ring))
        elif isinstance(item,pygplates.DateLineWrapper.LatLonPolyline):
            ring=[(p.get_longitude(),p.get_latitude()) for p in item.get_points()]
            if len(ring)>=2:result.append(LineString(ring))
    return result


def plates(atlas):
    archive=atlas.one('earthbyte-plates-2019','*.zip')
    with zipfile.ZipFile(archive) as z: names=z.namelist()
    def extract(suffix): return atlas.unzip(archive,next(n for n in names if n.endswith(suffix)))
    rotations=pygplates.RotationModel(str(extract('Muller_etal_2019_CombinedRotations.rot')))
    topology=pygplates.FeatureCollection(str(extract('Muller_etal_2019_PlateBoundaries_DeformingNetworks.gpmlz')))
    coast=pygplates.FeatureCollection(str(extract('Muller_etal_2019_Global_Coastlines.gpmlz')))
    for epoch in (0,100):
        resolved=[];sections=[]
        pygplates.resolve_topologies(topology,rotations,resolved,epoch,sections,
                                     resolve_topology_types=pygplates.ResolveTopologyType.boundary)
        geometries=[];codes=[]
        for boundary in resolved:
            if not isinstance(boundary,pygplates.ResolvedTopologicalBoundary):continue
            pid=boundary.get_feature().get_reconstruction_plate_id()
            pieces=wrap_geometry(boundary.get_resolved_boundary())
            geometries.extend(pieces);codes.extend([pid]*len(pieces))
        if not geometries:raise ValueError('No plate topologies resolved')
        unique=sorted(set(codes));colours=class_colours([str(p) for p in unique],unique)
        classes={p:colours[i+1] for i,p in enumerate(unique)}
        note='Müller et al. 2019 resolved rigid plate topologies; deformation-network interiors can be uncovered. One-degree great-circle tessellation and dateline splitting; anchor plate 0.'
        atlas.emit(Layer(f'tectonics/plates/{epoch}ma',f'Resolved plate IDs at {epoch} Ma','plate ID','earthbyte-plates-2019',(min(unique),max(unique)),classes=classes,
                        period=f'{epoch} Ma',notes=note,method='Resolved plate polygon at target cell centre',palette_status='Author-selected stable distinct plate colours'),
                   lambda w:draw_geometry(geometries,codes,w))
        outlines=[g.boundary for g in geometries]
        atlas.emit(Layer(f'tectonics/plate_outlines/{epoch}ma',f'Plate polygon outlines at {epoch} Ma','line-cell presence','earthbyte-plates-2019',(0,1),classes={1:('Plate outline','#ffb86b')},
                        period=f'{epoch} Ma',notes=note+' Includes clipping edges at the map seam; outlines are not classified fault mechanisms.',
                        method='Resolved polygon boundary touching target cell',palette_status='Author-selected boundary colour'),
                   lambda w:draw_geometry(outlines,np.ones(len(outlines)),w,all_touched=True))
        reconstructed=[]
        pygplates.reconstruct(coast,rotations,reconstructed,epoch)
        coast_geoms=[]
        for feature in reconstructed:coast_geoms.extend(wrap_geometry(feature.get_reconstructed_geometry()))
        atlas.emit(Layer(f'tectonics/reconstructed_land/{epoch}ma',f'Reconstructed coastline geometry at {epoch} Ma','land polygon presence','earthbyte-plates-2019',(0,1),
                        classes={0:('Outside reconstructed land geometry','#193445'),1:('Reconstructed land geometry','#bbc68b')},
                        period=f'{epoch} Ma',notes='Rotated present-day coastline geometry; not a reconstruction of ancient sea level, palaeoshorelines or palaeotopography.',
                        method='Reconstructed source polygon at cell centre',palette_status='Author-selected land/ocean colours'),
                   lambda w:draw_geometry(coast_geoms,np.ones(len(coast_geoms)),w,fill=0))
        if epoch==0:
            # A compact source mesh is sufficient for smoothly varying rigid rotations.
            width=512;mesh=draw_geometry(geometries,codes,width)
            lat=90-(np.arange(width//2)+.5)*360/width;lon=-180+(np.arange(width)+.5)*360/width
            north=np.full(mesh.shape,np.nan);east=north.copy()
            for pid in unique:
                yy,xx=np.where(mesh==pid)
                if not len(xx):continue
                points=list(zip(lat[yy],lon[xx]))
                rotation=rotations.get_rotation(0,pid,1)
                velocities=pygplates.calculate_velocities(points,rotation,1,velocity_units=pygplates.VelocityUnits.cms_per_yr)
                local=pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(points,velocities)
                vectors=np.array([v.to_xyz() for v in local])*10
                north[yy,xx]=vectors[:,0];east[yy,xx]=vectors[:,1]
            for key,data in [('north',north),('east',east),('speed',np.hypot(north,east))]:
                atlas.emit(Layer('tectonics/plate_velocity_'+key,'Rigid plate velocity '+key,'mm yr-1','earthbyte-plates-2019',
                                (0,150) if key=='speed' else (-150,150),'cmo.speed' if key=='speed' else 'cmo.balance',
                                period='1 Ma to present stage rotation',notes=note+' Velocities sampled on a 512x256 mesh; deformation within plates is excluded.',
                                method='GPlates finite-rotation velocity, geocentric north/east, then spherical-area remap'),
                           lambda w,a=data:regular_grid(a,lat,lon,w))
