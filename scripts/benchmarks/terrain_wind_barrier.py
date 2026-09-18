"""Optional global Q3 diagnostic: steady shallow-layer volume projection.

Periodic longitude, closed polar faces, cell-centred spherical finite volumes.
This redistributes volume flux, not momentum; it does not feed the climate engine.
An open top regularizes connected wet/lowland regions. Mountains close faces.
Requires scipy and pyamg only for this offline experiment.
"""
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg
import pyamg


class GlobalBarrier:
    def __init__(self, height, layer_top_m=500., vertical_mobility=1e-8):
        start = time.perf_counter()
        h = np.asarray(height, dtype=float)
        if (h.ndim != 2 or h.shape[1] != 2*h.shape[0] or h.shape[0] < 4
                or not np.isfinite(h).all() or not np.isfinite([layer_top_m, vertical_mobility]).all()
                or min(layer_top_m, vertical_mobility) <= 0):
            raise ValueError('Finite, global cell-centred W x W/2 terrain and positive parameters required')
        ny, nx = h.shape
        radius = 6371000.
        step = np.pi / ny
        lat = np.pi/2 - (np.arange(ny)+.5)*step
        edges = np.pi/2 - np.arange(ny+1)*step
        self.depth = np.maximum(layer_top_m-np.maximum(h, 0), 0)
        self.dx = radius*np.cos(lat)*step
        self.dy = radius*step
        self.area = np.broadcast_to((radius**2*step*(np.sin(edges[:-1])-np.sin(edges[1:])))[:, None], h.shape)
        self.ax = np.minimum(self.depth, np.roll(self.depth, -1, axis=1))*self.dy
        self.ay = np.zeros((ny+1, nx))
        self.ay[1:-1] = np.minimum(self.depth[:-1], self.depth[1:])*(radius*np.cos(edges[1:-1])*step)[:, None]
        self.cx = self.ax/self.dx[:, None]
        self.cy = self.ay[1:-1]/self.dy
        self.ct = self.area*vertical_mobility/layer_top_m*(self.depth > 0)
        diagonal = self.ct+self.cx+np.roll(self.cx, 1, axis=1)
        diagonal[:-1] += self.cy
        diagonal[1:] += self.cy
        diagonal[self.depth == 0] = 1
        ids = np.arange(ny*nx).reshape(h.shape)
        first = np.r_[ids.ravel(), ids[:-1].ravel()]
        second = np.r_[np.roll(ids, -1, axis=1).ravel(), ids[1:].ravel()]
        coeff = np.r_[self.cx.ravel(), self.cy.ravel()]
        nonzero = coeff > 0
        first, second, coeff = first[nonzero], second[nonzero], coeff[nonzero]
        self.matrix = sparse.coo_matrix((np.r_[diagonal.ravel(), -coeff, -coeff],
            (np.r_[ids.ravel(), first, second], np.r_[ids.ravel(), second, first])),
            shape=(ny*nx, ny*nx)).tocsr()
        self.multigrid = pyamg.ruge_stuben_solver(self.matrix,
            presmoother=('gauss_seidel', {'sweep': 'symmetric'}),
            postsmoother=('gauss_seidel', {'sweep': 'symmetric'}), max_coarse=64)
        self.setup_seconds = time.perf_counter()-start

    @staticmethod
    def divergence(qx, qy, top):
        return qx-np.roll(qx, 1, axis=1)+qy[1:]-qy[:-1]+top

    def solve(self, east, north, tolerance=1e-8, limit=500):
        start = time.perf_counter()
        u, v = [np.asarray(a, dtype=float) for a in (east, north)]
        if (u.shape != self.depth.shape or v.shape != u.shape or not np.isfinite(u).all()
                or not np.isfinite(v).all() or not 0 < tolerance < 1 or limit < 1):
            raise ValueError('Invalid winds or solver settings')
        qx = .5*(u+np.roll(u, -1, axis=1))*self.ax
        qy = np.zeros_like(self.ay)
        qy[1:-1] = -.5*(v[:-1]+v[1:])*self.ay[1:-1]
        rhs = self.divergence(qx, qy, np.zeros_like(u))
        iterations = [0]
        def tick(_): iterations[0] += 1
        p, info = cg(self.matrix, rhs.ravel(), M=self.multigrid.aspreconditioner(),
            rtol=tolerance, atol=0, maxiter=limit, callback=tick)
        p = p.reshape(u.shape)
        qx += self.cx*(np.roll(p, -1, axis=1)-p)
        qy[1:-1] += self.cy*(p[1:]-p[:-1])
        top = -self.ct*p
        norm0 = np.linalg.norm(rhs)
        error = self.divergence(qx, qy, top)
        relative = float(np.linalg.norm(error)/max(norm0, 1e-30))
        if info or (norm0 > 1e-20 and relative > 1.2*tolerance):
            raise RuntimeError(f'Global barrier did not converge: info={info}, relative={relative}')
        east = np.divide(qx+np.roll(qx, 1, axis=1), 2*self.depth*self.dy,
            out=np.zeros_like(u), where=self.depth > 0)
        north = -np.divide(qy[:-1]+qy[1:], 2*self.depth*self.dx[:, None],
            out=np.zeros_like(v), where=self.depth > 0)
        assert np.all(qx[self.ax == 0] == 0) and np.all(qy[self.ay == 0] == 0)
        return dict(east=east, north=north, qx=qx, qy=qy, top=top, depth=self.depth,
            stats=dict(iterations=iterations[0], relative_flux_residual=relative,
                setup_seconds=self.setup_seconds, solve_seconds=time.perf_counter()-start,
                top_outflow_m3s=float(top.sum()), budget_error_m3s=float(error.sum()),
                absolute_flux_scale_m3s=float(np.abs(qx).sum()+np.abs(qy).sum()+np.abs(top).sum()),
                maximum_top_speed_mps=float(np.max(np.abs(top)/self.area))))
