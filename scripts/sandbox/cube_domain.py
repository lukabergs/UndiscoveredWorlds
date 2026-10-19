"""Cube grids and domains of the terrain sandbox (docs/architecture/sphere-mesh.md sections 2, 3, 5, 7, 11, 12).

A CubeGrid holds a cube's cell table as C++ writes it (a bundle v2's grid/*.npy, or `terrain_sandbox grid`
for a domain's finer cube: Python never builds the topology) and the vectorised operations the sandbox's
raster code had on its W x W/2 grid, here as gathers over the eight neighbour slots (-1 none; a missing
slot reads the cell itself):
- smooth(values, sigma): repeated 3 x 3 means (each pass adds a variance of 2/3 cell^2 per axis), the
  Gaussian of the raster code;
- label(mask): 8-connected components (scipy.sparse.csgraph);
- distance_km(sources, cutoff_km): the great-circle distance to the nearest source cell (its id carried through a
  relaxation over the slots; Dijkstra seeds it without a cutoff), the raster code's distance transform;
- minimum and maximum over (2 r + 1)-cell neighbourhoods (r passes over the slots), dilation and erosion;
- gradient(values): (d/d east, d/d north) per km, least squares over the neighbours' tangent offsets;
- locate(points) and bilinear(values, points): the closed-form point location of section 2 and the bilinear
  stencil of section 5 (in face coordinates between the four surrounding centres, cells across a face edge from
  the table, a missing corner diagonal's weight shared by the other three);
- priority_flood and sloped_fill: the raster code's Priority-Floods on the slots.
A CubeDomain is the sandbox's global domain on a cube of N_d = k N cells per face edge over a base cube of N: its
cells' unit vectors (km) carry the noise (section 11); world fields reach it by the bilinear stencil from the base
followed by smooth(0.41 k) (the cubic B-spline's variance, 1/6 base cell^2, the raster code's sampler); the
physics runs on it in terrain_sandbox.exe (--face-cells N_d) and its renders come from the exe's lat-lon export.
"""
import heapq
from collections import deque
import math
import subprocess
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph

# uw-cube-1 (section 2): centre c, right axis u, up axis v of faces 0-5.
FACES = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
                  [[-1, 0, 0], [0, -1, 0], [0, 0, 1]], [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                  [[0, 0, 1], [0, 1, 0], [-1, 0, 0]], [[0, 0, -1], [0, 1, 0], [1, 0, 0]]], dtype=np.float64)
CARDINAL = (0, 2, 4, 6)  # columns of slots 1 N, 3 E, 5 S, 7 W


class CubeGrid:
    """A cube's cell table (centre, latitude_deg, longitude_deg, area_m2, neighbours, reverse) and its operations."""

    def __init__(self, cells, face_cells, radius_km):
        self.N = int(face_cells)
        self.n = 6 * self.N * self.N
        self.radius_km = float(radius_km)
        self.centre = np.asarray(cells['centre'], dtype=np.float64)
        if self.centre.shape != (self.n, 3):
            raise ValueError(f'cell table of {self.centre.shape[0]} cells, expected 6 x {self.N}^2')
        self.latitude = np.radians(np.asarray(cells['latitude_deg'], dtype=np.float64))
        self.longitude = np.radians(np.asarray(cells['longitude_deg'], dtype=np.float64))
        self.area_km2 = np.asarray(cells['area_m2'], dtype=np.float64) / 1e6
        self.neighbours = np.asarray(cells['neighbours'], dtype=np.int64)
        self.valid = self.neighbours >= 0
        self.safe = np.where(self.valid, self.neighbours, np.arange(self.n)[:, None])
        sl, cl = np.sin(self.latitude), np.cos(self.latitude)
        so, co = np.sin(self.longitude), np.cos(self.longitude)
        self.east = np.stack([-so, co, np.zeros_like(so)])
        self.north = np.stack([-sl * co, -sl * so, cl])
        # Great-circle distances (km) and tangent offsets (east, north; km) of the slots; inf and 0 where none.
        other = self.centre[self.safe]                                   # (n, 8, 3)
        cosine = np.clip(np.einsum('nkd,nd->nk', other, self.centre), -1.0, 1.0)
        self.slot_km = np.where(self.valid, self.radius_km * np.arccos(cosine), np.inf)
        oe = np.einsum('nkd,dn->nk', other, self.east)
        on = np.einsum('nkd,dn->nk', other, self.north)
        length = np.maximum(np.hypot(oe, on), 1e-300)
        scale = np.where(self.valid, self.slot_km / length, 0.0)
        self.offset_east, self.offset_north = oe * scale, on * scale
        see = (self.offset_east ** 2).sum(1)
        snn = (self.offset_north ** 2).sum(1)
        sen = (self.offset_east * self.offset_north).sum(1)
        determinant = see * snn - sen * sen
        self.inverse = (snn / determinant, -sen / determinant, see / determinant)
        self.count = 1.0 + self.valid.sum(1)
        self._graph = None

    @classmethod
    def from_bundle(cls, bundle):
        return cls(bundle.cells(), bundle.face_cells, bundle.radius_km or 6371.0)

    @classmethod
    def from_directory(cls, directory, face_cells, radius_km):
        names = ('centre', 'latitude_deg', 'longitude_deg', 'area_m2', 'neighbours', 'reverse')
        return cls({name: np.load(Path(directory) / f'{name}.npy') for name in names}, face_cells, radius_km)

    # -------------------------------------------------------------- neighbourhood operations

    def gather(self, values):
        """(n, 8, ...) neighbour values; a missing slot reads the cell itself."""
        return np.asarray(values)[self.safe]

    def smooth(self, values, sigma_cells):
        """A Gaussian of sigma_cells cells (about): round(1.5 sigma^2) passes of the 3 x 3 mean over the existing slots."""
        out = np.asarray(values, dtype=np.float64)
        passes = int(round(1.5 * float(sigma_cells) ** 2))
        valid = self.valid if out.ndim == 1 else self.valid[..., None]
        count = self.count if out.ndim == 1 else self.count[:, None]
        for _ in range(passes):
            out = (out + np.where(valid, out[self.safe], 0.0).sum(1)) / count
        return out

    def minimum(self, values, radius):
        out = np.asarray(values)
        for _ in range(int(radius)):
            out = np.minimum(out, out[self.safe].min(1))
        return out

    def maximum(self, values, radius):
        out = np.asarray(values)
        for _ in range(int(radius)):
            out = np.maximum(out, out[self.safe].max(1))
        return out

    def dilate(self, mask, iterations=1):
        return self.maximum(np.asarray(mask, dtype=bool), iterations)

    def erode(self, mask, iterations=1):
        return self.minimum(np.asarray(mask, dtype=bool), iterations)

    def gradient(self, values):
        """(d/d east, d/d north) of a cell field per km: least squares over the neighbours' tangent offsets."""
        values = np.asarray(values, dtype=np.float64)
        change = np.where(self.valid, values[self.safe] - values[:, None], 0.0)
        be = (self.offset_east * change).sum(1)
        bn = (self.offset_north * change).sum(1)
        a, b, c = self.inverse
        return a * be + b * bn, b * be + c * bn

    def graph(self):
        """The slots as a sparse matrix of great-circle lengths (km)."""
        if self._graph is None:
            rows = np.repeat(np.arange(self.n), 8).reshape(self.n, 8)[self.valid]
            self._graph = sparse.csr_matrix((self.slot_km[self.valid], (rows, self.neighbours[self.valid])),
                                            shape=(self.n, self.n))
        return self._graph

    def label(self, mask):
        """8-connected components of a mask: (labels, count), 0 off the mask, 1.. in the order of their lowest cell."""
        mask = np.asarray(mask, dtype=bool)
        cells = np.flatnonzero(mask)
        labels = np.zeros(self.n, dtype=np.int64)
        if cells.size == 0:
            return labels, 0
        local = np.full(self.n, -1, dtype=np.int64)
        local[cells] = np.arange(cells.size)
        nbr = self.neighbours[cells]
        keep = (nbr >= 0) & mask[np.maximum(nbr, 0)]
        rows = np.repeat(np.arange(cells.size), 8).reshape(-1, 8)[keep]
        graph = sparse.csr_matrix((np.ones(rows.size, dtype=np.int8), (rows, local[nbr[keep]])),
                                  shape=(cells.size, cells.size))
        count, component = csgraph.connected_components(graph, directed=False)
        # Renumber by the lowest cell of each component (deterministic, ascending ids).
        first = np.full(count, cells.size, dtype=np.int64)
        np.minimum.at(first, component, np.arange(cells.size))
        order = np.empty(count, dtype=np.int64)
        order[np.argsort(first, kind='stable')] = np.arange(1, count + 1)
        labels[cells] = order[component]
        return labels, int(count)

    def great_circle_km(self, cells, others):
        """Great-circle distance (km) between the centres of two arrays of cell ids (2 asin(chord / 2): exact near 0)."""
        chord = np.linalg.norm(self.centre[cells] - self.centre[others], axis=-1)
        return 2.0 * self.radius_km * np.arcsin(np.minimum(chord / 2.0, 1.0))

    def distance_km(self, sources, cutoff_km=None):
        """Great-circle distance (km) from each cell's centre to the centre of its nearest source cell (0 on one, inf
        with no source): the nearest source is carried through a relaxation over the slots (each cell takes the source
        of a neighbour when that source is nearer to it), started from the sources themselves within cutoff_km (cells
        farther than it inf), or without a cutoff from Dijkstra's nearest source over the slot paths, then relaxed
        until no cell changes."""
        sources = np.asarray(sources, dtype=bool)
        if not sources.any():
            return np.full(self.n, np.inf)
        if cutoff_km is None:
            _, _, nearest = csgraph.dijkstra(self.graph(), directed=True, indices=np.flatnonzero(sources),
                                             min_only=True, return_predecessors=True)
            nearest = np.where(nearest >= 0, nearest, -1).astype(np.int64)
            passes = 64
        else:
            nearest = np.where(sources, np.arange(self.n), -1)
            passes = int(math.ceil(cutoff_km / float(np.min(self.slot_km)))) + 2
        cells = np.arange(self.n)
        distance = np.where(nearest >= 0, self.great_circle_km(cells, np.maximum(nearest, 0)), np.inf)
        for _ in range(passes):
            offered = nearest[self.safe]                                         # (n, 8) the neighbours' sources
            reach = np.where(offered >= 0, self.great_circle_km(cells[:, None], np.maximum(offered, 0)), np.inf)
            best = np.argmin(reach, axis=1)
            better = reach[cells, best] < distance
            if not better.any():
                break
            nearest = np.where(better, offered[cells, best], nearest)
            distance = np.where(better, reach[cells, best], distance)
        if cutoff_km is not None:
            distance = np.where(distance <= cutoff_km, distance, np.inf)
        return distance

    # -------------------------------------------------------------- location and sampling

    def locate(self, points):
        """Cell ids of points (3, ...) and their positions (s_a, s_b) in cell units from the cell's centre, in
        [-0.5, 0.5): the closed form of section 2 (the face of the largest |coordinate|, ties to the lower face)."""
        p = np.asarray(points, dtype=np.float64)
        shape = p.shape[1:]
        p = p.reshape(3, -1)
        p = p / np.linalg.norm(p, axis=0)
        face = np.argmax(np.stack([p[0], p[1], -p[0], -p[1], p[2], -p[2]]), axis=0)
        c, u, v = FACES[face, 0].T, FACES[face, 1].T, FACES[face, 2].T
        pc = (p * c).sum(0)
        a = np.arctan((p * u).sum(0) / pc)
        b = np.arctan((p * v).sum(0) / pc)
        step = math.pi / (2 * self.N)
        sa = (a + math.pi / 4) / step
        sb = (math.pi / 4 - b) / step
        i = np.clip(np.floor(sa), 0, self.N - 1).astype(np.int64)
        j = np.clip(np.floor(sb), 0, self.N - 1).astype(np.int64)
        ids = face * self.N * self.N + j * self.N + i
        return ids.reshape(shape), (sa - i - 0.5).reshape(shape), (sb - j - 0.5).reshape(shape)

    def bilinear(self, values, points):
        """values (n,) at points (3, ...): the bilinear stencil of section 5 on the cell table."""
        values = np.asarray(values, dtype=np.float64)
        ids, ta, tb = self.locate(points)
        east = ta >= 0.0
        south = tb >= 0.0
        wa, wb = np.abs(ta), np.abs(tb)
        horizontal = np.where(east, 2, 6)                              # slot 3 E or 7 W (columns)
        vertical = np.where(south, 4, 0)                               # slot 5 S or 1 N
        diagonal = np.where(east, np.where(south, 3, 1), np.where(south, 5, 7))   # 4 SE, 2 NE, 6 SW, 8 NW
        h, vv, d = (self.neighbours[ids, slot] for slot in (horizontal, vertical, diagonal))
        weights = [(1 - wa) * (1 - wb), wa * (1 - wb), (1 - wa) * wb, wa * wb]
        missing = d < 0                                                # a cube corner's diagonal
        share = np.where(missing, weights[3], 0.0)
        rest = weights[0] + weights[1] + weights[2]
        weights = [w * (1.0 + share / np.maximum(rest, 1e-300)) for w in weights[:3]] + [np.where(missing, 0.0, weights[3])]
        cells = [ids, h, vv, np.where(missing, ids, d)]
        first = values[ids]
        total = first.copy()
        for w, cell in zip(weights[1:], cells[1:]):
            total = total + w * (values[cell] - first)
        return total

    def resample(self, values, points_km, east_km, north_km):
        """values at the cells' points moved by (east_km, north_km) on the tangent plane (bilinear)."""
        moved = points_km + east_km * self.east + north_km * self.north
        return self.bilinear(values, moved)

    # -------------------------------------------------------------- Priority-Floods

    def priority_flood(self, elevation, sea):
        """The filled surface (Priority-Flood, Barnes et al. 2014, over the eight slots): every land cell raised to the
        lowest level at which water on it reaches a sea cell; the sea cells the outlets at their own elevation."""
        elevation = np.asarray(elevation, dtype=np.float64)
        sea = np.asarray(sea, dtype=bool)
        neighbours = self.neighbours.tolist()
        height = elevation.tolist()
        filled = list(height)
        closed = sea.tolist()
        shore = sea & ~self.minimum(sea, 1)
        heap = [(height[i], i) for i in np.flatnonzero(shore).tolist()]
        heapq.heapify(heap)
        pit = deque()
        while heap or pit:
            if pit:
                cell = pit.popleft()
                level = filled[cell]
            else:
                level, cell = heapq.heappop(heap)
            for j in neighbours[cell]:
                if j < 0 or closed[j]:
                    continue
                closed[j] = True
                if height[j] <= level:
                    filled[j] = level
                    pit.append(j)
                else:
                    heapq.heappush(heap, (height[j], j))
        return np.where(sea, elevation, np.asarray(filled))

    def sloped_fill(self, elevation, sea, gradient, noise=None):
        """The raster code's sloped_fill on the slots: a cell lower than the cell reaching it is raised to that level
        plus max(0.1 g d, g d + n_j - n_i) (d the great-circle distance, m) and goes back into the heap."""
        elevation = np.asarray(elevation, dtype=np.float64)
        sea = np.asarray(sea, dtype=bool)
        neighbours = self.neighbours.tolist()
        rises = (gradient * np.where(self.valid, self.slot_km, 0.0) * 1000.0).tolist()
        height = elevation.tolist()
        lift = (np.zeros(self.n) if noise is None else np.asarray(noise, dtype=np.float64)).tolist()
        filled = list(height)
        closed = sea.tolist()
        shore = sea & ~self.minimum(sea, 1)
        heap = [(height[i], i) for i in np.flatnonzero(shore).tolist()]
        heapq.heapify(heap)
        while heap:
            level, cell = heapq.heappop(heap)
            run = rises[cell]
            for k, j in enumerate(neighbours[cell]):
                if j < 0 or closed[j]:
                    continue
                closed[j] = True
                if height[j] <= level:
                    filled[j] = level + max(0.1 * run[k], run[k] + lift[j] - lift[cell])
                    heapq.heappush(heap, (filled[j], j))
                else:
                    heapq.heappush(heap, (height[j], j))
        return np.where(sea, elevation, np.asarray(filled))


def cube_table(exe, face_cells, radius_km, directory):
    """The cell table of a cube of face_cells from `terrain_sandbox grid` (cached in directory/cells-<N>)."""
    target = Path(directory) / f'cells-{face_cells}'
    if not (target / 'reverse.npy').is_file():
        target.mkdir(parents=True, exist_ok=True)
        result = subprocess.run([str(exe), 'grid', '--face-cells', str(face_cells), '--radius-km', str(radius_km),
                                 '--out', str(target)], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f'terrain_sandbox grid failed: {result.stderr.strip()}')
    return CubeGrid.from_directory(target, face_cells, radius_km)


class CubeDomain:
    """The global domain on a cube of k N cells per face edge (the attributes the sandbox's code reads of a Domain)."""

    kind = 'cube'

    def __init__(self, name, grid, base_face_cells, radius_km, exaggeration=10.0, physics=None, smoothing=0.41):
        self.name, self.grid = name, grid
        self.k = grid.N // base_face_cells
        self.base_face_cells = base_face_cells
        self.radius_km = radius_km
        self.cell_km = 2.0 * math.pi * radius_km / (4 * base_face_cells)
        self.fine_km = self.cell_km / self.k
        self.world_columns, self.world_rows = 4 * base_face_cells, 2 * base_face_cells
        self.columns, self.rows = 4 * grid.N, 2 * grid.N  # the lat-lon export's
        self.crop = None
        self.whole = True
        self.exaggeration = exaggeration
        self.physics = physics or {}
        self.smoothing = smoothing
        self._points = None

    def points(self):
        """The cells' unit vectors times the radius (3, n), km: the noise's positions (section 11)."""
        if self._points is None:
            self._points = self.radius_km * self.grid.centre.T
        return self._points

    def frame(self):
        return self.grid.east, self.grid.north

    def abs_latitude(self):
        return np.abs(np.degrees(self.grid.latitude))

    def area_km2(self):
        return self.grid.area_km2

    def sampler(self, base):
        """World fields on the base cube to this domain: the bilinear stencil, then smooth(smoothing x k)."""
        points = self.points()
        sigma = self.smoothing * self.k

        def sample(field):
            return self.grid.smooth(base.bilinear(np.asarray(field, dtype=np.float64).reshape(-1), points), sigma)
        return sample


def export(exe, domain, arrays, directory, columns=None):
    """Lat-lon rasters (columns x columns / 2) of cube fields by `terrain_sandbox export`: arrays maps a name to
    (values (n,), rule: bilinear, mean or nearest). The fields go stacked through <domain>-export-in.f32 and
    <domain>-export-out.f32 in directory (overwritten by the next export). Returns name -> (rows, columns) arrays."""
    directory = Path(directory)
    columns = columns or domain.columns
    names = list(arrays)
    source = directory / f'{domain.name}-export-in.f32'
    target = directory / f'{domain.name}-export-out.f32'
    np.concatenate([np.asarray(arrays[name][0], dtype='<f4').reshape(-1) for name in names]).tofile(source)
    result = subprocess.run([str(exe), 'export', '--face-cells', str(domain.grid.N), '--columns', str(columns),
                             '--radius-km', str(domain.radius_km), '--input', str(source), '--output', str(target),
                             '--rules', ','.join(arrays[name][1] for name in names)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'terrain_sandbox export failed: {result.stderr.strip()}')
    stacked = np.fromfile(target, dtype='<f4').reshape(len(names), columns // 2, columns).astype(np.float64)
    return {name: stacked[index] for index, name in enumerate(names)}
