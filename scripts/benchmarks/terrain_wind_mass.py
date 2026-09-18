"""Regional, diagnostic mass projection on terrain-following finite volumes.

Constant-density kinematics, not a momentum/turbulence solver. Horizontal
faces carry east/south volume flux (m3/s); vertical faces carry contravariant
W = w - u dh/dx - v_south dh/dy. The ground has W=0, lateral fluxes are
fixed, and the top is open. A diagonal metric approximation omits the full
terrain-coordinate cross terms in the least-change norm, not in continuity.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MassConfig:
    layer_depths_m: tuple = (150., 350., 1000.)
    terrain_decay_m: float = 350.
    vertical_mobility: float = 0.0001
    relative_tolerance: float = 1e-8
    maximum_iterations: int = 2000


class TerrainMassGrid:
    def __init__(self, height, latitudes_deg, longitude_step_deg, config=MassConfig(), radius=6371000.):
        self.h = np.asarray(height, dtype=float)
        self.lat = np.deg2rad(np.asarray(latitudes_deg, dtype=float))
        self.config = config
        if self.h.ndim != 2 or min(self.h.shape) < 3 or self.lat.shape != (self.h.shape[0],):
            raise ValueError('Need a regional height grid and descending cell-centre latitudes')
        dphi = self.lat[0] - self.lat[1]
        dlambda = np.deg2rad(longitude_step_deg)
        dz = np.asarray(config.layer_depths_m, dtype=float)
        if (not np.isfinite(self.h).all() or not np.isfinite(self.lat).all()
                or not np.isfinite(dz).all() or not np.isfinite(dlambda)
                or dlambda <= 0 or dphi <= 0 or not np.allclose(np.diff(self.lat), -dphi)
                or (dz <= 0).any() or config.vertical_mobility <= 0
                or config.terrain_decay_m <= 0 or config.relative_tolerance <= 0
                or config.maximum_iterations < 1):
            raise ValueError('Invalid geometry or solver configuration')
        edges = np.r_[self.lat + dphi/2, self.lat[-1]-dphi/2]
        if np.max(np.abs(edges)) >= np.pi/2:
            raise ValueError('Regional solver excludes polar caps')
        ny, nx = self.h.shape
        self.dz = dz
        self.area = np.broadcast_to(radius**2*dlambda*(np.sin(edges[:-1])-np.sin(edges[1:]))[:, None], (ny, nx))
        self.volume = dz[:, None, None]*self.area
        self.dx = radius*np.cos(self.lat)*dlambda
        self.dy = radius*dphi
        self.ax = np.broadcast_to(dz[:, None, None]*self.dy, (len(dz), ny, nx+1))
        self.ay = np.broadcast_to(dz[:, None, None]*radius*np.cos(edges)[None, :, None]*dlambda, (len(dz), ny+1, nx))
        self.hx = np.gradient(self.h, axis=1)/self.dx[:, None]
        self.hy = np.gradient(self.h, axis=0)/self.dy
        self.cx = self.ax[:, :, 1:-1]/self.dx[None, :, None]
        self.cy = self.ay[:, 1:-1, :]/self.dy
        self.cz = self.area[None]*config.vertical_mobility/((dz[:-1]+dz[1:])/2)[:, None, None]
        self.ct = self.area*config.vertical_mobility/(dz[-1]/2)
        self.diagonal = np.zeros_like(self.volume)
        self.diagonal[:, :, :-1] += self.cx; self.diagonal[:, :, 1:] += self.cx
        self.diagonal[:, :-1, :] += self.cy; self.diagonal[:, 1:, :] += self.cy
        self.diagonal[:-1] += self.cz; self.diagonal[1:] += self.cz
        self.diagonal[-1] += self.ct

    @staticmethod
    def divergence(qx, qy, qz):
        return np.diff(qx, axis=2)+np.diff(qy, axis=1)+np.diff(qz, axis=0)

    def correction(self, pressure):
        qx = np.zeros_like(self.ax); qy = np.zeros_like(self.ay)
        qz = np.zeros((len(self.dz)+1, *self.h.shape))
        qx[:, :, 1:-1] = self.cx*np.diff(pressure, axis=2)
        qy[:, 1:-1, :] = self.cy*np.diff(pressure, axis=1)
        qz[1:-1] = self.cz*np.diff(pressure, axis=0)
        qz[-1] = -self.ct*pressure[-1]
        return qx, qy, qz

    def operator(self, pressure):
        return -self.divergence(*self.correction(pressure))

    def initialize(self, east, north, parent_east, parent_north):
        fields = [np.asarray(a, dtype=float) for a in (east, north, parent_east, parent_north)]
        if any(a.shape != self.h.shape or not np.isfinite(a).all() for a in fields):
            raise ValueError('Wind fields must be finite and match the terrain grid')
        u, v, bu, bv = fields
        centres = np.cumsum(self.dz)-self.dz/2
        taper = np.exp(-(centres-centres[0])/self.config.terrain_decay_m)[:, None, None]
        u = bu+taper*(u-bu); south = -(bv+taper*(v-bv))
        ux = np.concatenate((u[:, :, :1], (u[:, :, :-1]+u[:, :, 1:])/2, u[:, :, -1:]), axis=2)
        sy = np.concatenate((south[:, :1], (south[:, :-1]+south[:, 1:])/2, south[:, -1:]), axis=1)
        # Physical w=0 initially above the ground; impermeability replaces the
        # lower boundary. A solve is therefore needed even for uniform ridge flow.
        W = -(u*self.hx+south*self.hy)
        a = self.dz[1:, None, None]; b = self.dz[:-1, None, None]
        faceW = np.concatenate((np.zeros_like(W[:1]), (a*W[:-1]+b*W[1:])/(a+b), W[-1:]), axis=0)
        return ux*self.ax, sy*self.ay, faceW*self.area

    def solve(self, east, north, parent_east, parent_north):
        initial = self.initialize(east, north, parent_east, parent_north)
        rhs = self.divergence(*initial)
        p = np.zeros_like(rhs); residual = rhs.copy()
        norm0 = np.sqrt(np.sum(rhs*rhs))
        z = residual/self.diagonal; direction = z.copy(); rz = np.sum(residual*z)
        iterations = 0
        while np.sqrt(np.sum(residual*residual)) > self.config.relative_tolerance*norm0 and iterations < self.config.maximum_iterations:
            ap = self.operator(direction)
            alpha = rz/np.sum(direction*ap)
            p += alpha*direction; residual -= alpha*ap
            z = residual/self.diagonal; next_rz = np.sum(residual*z)
            direction = z+(next_rz/rz)*direction; rz = next_rz
            iterations += 1
        flux = tuple(q+c for q,c in zip(initial,self.correction(p)))
        final = self.divergence(*flux)
        relative = np.sqrt(np.sum(final*final))/max(norm0,1e-30)
        if relative > self.config.relative_tolerance*1.1 and norm0 > 1e-20:
            raise RuntimeError(f'Mass projection failed: {iterations} iterations, residual {relative}')
        qx,qy,qz = flux
        u = (qx[:, :, :-1]+qx[:, :, 1:])/(2*self.ax[:, :, :-1])
        south = (qy[:, :-1]/self.ay[:, :-1]+qy[:, 1:]/self.ay[:, 1:])/2
        w = (qz[:-1]+qz[1:])/(2*self.area)+u*self.hx+south*self.hy
        lateral = float(np.sum(qx[:, :, -1]-qx[:, :, 0])+np.sum(qy[:, -1]-qy[:, 0]))
        stats = dict(iterations=iterations, relative_flux_residual=float(relative),
            initial_rms_divergence_per_second=float(np.sqrt(np.mean((rhs/self.volume)**2))),
            final_rms_divergence_per_second=float(np.sqrt(np.mean((final/self.volume)**2))),
            maximum_divergence_per_second=float(np.max(np.abs(final/self.volume))),
            lateral_outflow_m3s=lateral, top_outflow_m3s=float(qz[-1].sum()),
            budget_error_m3s=float(lateral+qz[-1].sum()),
            maximum_physical_vertical_speed_mps=float(np.abs(w).max()))
        return dict(east=u[0], north=-south[0], east_layers=u, north_layers=-south,
                    vertical_layers=w, qx=qx, qy=qy, qz=qz, stats=stats)
