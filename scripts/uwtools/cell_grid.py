"""Shared spherical remapping for processed global reference grids, and measures on the native fields of
a cube bundle through its cell table (docs/architecture/sphere-mesh.md sections 2, 3 and 7)."""

from __future__ import annotations

import numpy as np

from . import grid

# Columns of the cell table's neighbours array for the edge slots 1 N, 3 E, 5 S, 7 W (slot k in column k - 1).
SLOT_NORTH, SLOT_EAST, SLOT_SOUTH, SLOT_WEST = 0, 2, 4, 6


def face_edge_distance(face_cells: int) -> np.ndarray:
    """(6, N, N): each cell's distance in cells to the nearest edge of its face (0 on the edge row)."""
    index = np.arange(face_cells)
    edge = np.minimum(index, face_cells - 1 - index)
    one = np.minimum(edge[:, None], edge[None, :])
    return np.broadcast_to(one, (6, face_cells, face_cells)).copy()


def _arc_km(centre: np.ndarray, first: np.ndarray, second: np.ndarray, radius_km: float) -> np.ndarray:
    a, b = centre[first], centre[second]
    return radius_km * np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.einsum('ij,ij->i', a, b))


def _tangent(centre: np.ndarray, at: np.ndarray, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """Unit tangent at the centres `at` of the chord from `start` to `end`."""
    c = centre[at]
    t = centre[end] - centre[start]
    t -= np.einsum('ij,ij->i', t, c)[:, None] * c
    return t / np.linalg.norm(t, axis=1, keepdims=True)


def face_gradients(cells: dict, values: np.ndarray, radius_km: float):
    """(g_u, g_v, valid) of a native field in face coordinates, per km, in id order (flat, n). The derivatives
    along the face's two grid lines through each cell are central differences over the edge-slot neighbours
    of the cell table (E minus W over the two great-circle distances between centres, N minus S likewise;
    slot directions are the cell's own face frame, also across face edges); the gradient they determine is
    returned in the orthonormal frame of the cell's u line (W to E) and the perpendicular towards v (S to N),
    so hypot(g_u, g_v) is the true |grad| where the grid lines meet at an angle (near face corners).
    valid is False where an edge neighbour is missing (-1)."""
    z = np.asarray(values, dtype=np.float64).reshape(-1)
    neighbours = np.asarray(cells['neighbours'])
    centre = np.asarray(cells['centre'], dtype=np.float64)
    if neighbours.shape != (z.size, 8) or centre.shape != (z.size, 3):
        raise ValueError(f'cell table of {neighbours.shape[0]} cells does not match a field of {z.size}')
    north, east, south, west = (neighbours[:, slot] for slot in (SLOT_NORTH, SLOT_EAST, SLOT_SOUTH, SLOT_WEST))
    valid = (north >= 0) & (east >= 0) & (south >= 0) & (west >= 0)
    ids = np.arange(z.size)
    g_u = np.full(z.size, np.nan)
    g_v = np.full(z.size, np.nan)
    c, n, e, s, w = ids[valid], north[valid], east[valid], south[valid], west[valid]
    d_u = (z[e] - z[w]) / (_arc_km(centre, c, e, radius_km) + _arc_km(centre, c, w, radius_km))
    d_v = (z[n] - z[s]) / (_arc_km(centre, c, n, radius_km) + _arc_km(centre, c, s, radius_km))
    k = np.einsum('ij,ij->i', _tangent(centre, c, w, e), _tangent(centre, c, s, n))  # cosine between the lines
    g_u[valid] = d_u
    g_v[valid] = (d_v - k * d_u) / np.sqrt(1.0 - k * k)
    return g_u, g_v, valid


def aspect_harmonics(theta: np.ndarray) -> tuple:
    """(A8, A4) of aspect angles: |mean e^(8 i theta)| and mean cos(4 theta) (None for no angles)."""
    theta = np.asarray(theta)
    if not theta.size:
        return None, None
    return float(np.abs(np.mean(np.exp(8j * theta)))), float(np.mean(np.cos(4.0 * theta)))


def grid_dimensions(width: int, height: int | None = None) -> tuple[int, int]:
    if width < 2 or width % 2:
        raise ValueError(f"Width must be even and at least 2: {width}")
    expected_height = width // 2
    if height is not None and height != expected_height:
        raise ValueError(
            f"Cell-centred global grids must be {width}x{expected_height}, not "
            f"{width}x{height}"
        )
    return width, expected_height


def latitude_centres(height: int) -> np.ndarray:
    if height < 1:
        raise ValueError("Height must be positive")
    return grid.latitudes(height)


def latitude_edges(height: int) -> np.ndarray:
    if height < 1:
        raise ValueError("Height must be positive")
    return np.linspace(90.0, -90.0, height + 1, dtype=np.float64)


def _latitude_edges_from_centres(centres: np.ndarray) -> np.ndarray:
    centres = np.asarray(centres, dtype=np.float64)
    if centres.ndim != 1 or centres.size < 2 or not np.isfinite(centres).all():
        raise ValueError("Source latitude centres must be a finite one-dimensional array")
    if centres[0] < centres[-1]:
        centres = centres[::-1]
    if np.any(np.diff(centres) >= 0.0):
        raise ValueError("Source latitude centres must be strictly monotonic")
    edges = np.empty(centres.size + 1, dtype=np.float64)
    # Preserve the actual latitude coverage (e.g. an ocean product ending at
    # 80 S); extending a regional edge to the pole invents reference data.
    edges[0] = min(90.0, centres[0] + (centres[0] - centres[1]) * 0.5)
    edges[-1] = max(-90.0, centres[-1] - (centres[-2] - centres[-1]) * 0.5)
    edges[1:-1] = 0.5 * (centres[:-1] + centres[1:])
    if np.any(np.diff(edges) >= 0.0):
        raise ValueError("Source latitude cells do not form a north-to-south grid")
    return edges


def _latitude_plan(source_edges: np.ndarray, destination_height: int):
    source_north = np.sin(np.deg2rad(source_edges[:-1]))
    source_south = np.sin(np.deg2rad(source_edges[1:]))
    destination_edges = latitude_edges(destination_height)
    destination_north = np.sin(np.deg2rad(destination_edges[:-1]))
    destination_south = np.sin(np.deg2rad(destination_edges[1:]))
    plan = []
    for north, south in zip(destination_north, destination_south):
        overlap = np.maximum(
            0.0,
            np.minimum(north, source_north) - np.maximum(south, source_south),
        )
        indices = np.flatnonzero(overlap > 1.0e-15)
        plan.append((indices, overlap[indices] / (north - south)))
    return plan


def _longitude_plan(source_width: int, destination_width: int, first_face: float = 0.0):
    source_step = 360.0 / source_width
    destination_step = 360.0 / destination_width
    segments = []
    for source in range(source_width):
        west = first_face + source * source_step
        east = west + source_step
        while west < 0.0:
            west += 360.0
            east += 360.0
        while west >= 360.0:
            west -= 360.0
            east -= 360.0
        if east <= 360.0:
            segments.append((west, east, source))
        else:
            segments.append((west, 360.0, source))
            segments.append((0.0, east - 360.0, source))
    segments.sort()
    plan = []
    cursor = 0
    for destination in range(destination_width):
        west = destination * destination_step
        east = west + destination_step
        indices = []
        weights = []
        while cursor < len(segments) and segments[cursor][1] <= west + 1.0e-12:
            cursor += 1
        candidate = cursor
        while candidate < len(segments) and segments[candidate][0] < east - 1.0e-12:
            source_west, source_east, source = segments[candidate]
            overlap = min(east, source_east) - max(west, source_west)
            if overlap > 1.0e-12:
                indices.append(source)
                weights.append(overlap / destination_step)
            candidate += 1
        plan.append((np.asarray(indices, dtype=np.int64), np.asarray(weights)))
    return plan


def _normalize_source(
    values: np.ndarray,
    source_latitudes: np.ndarray | None,
    source_longitudes: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, float]:
    source = np.asarray(values, dtype=np.float32)
    if source.ndim != 2:
        raise ValueError("Expected one two-dimensional global field")
    height, width = source.shape

    if source_latitudes is None:
        latitude_source_edges = latitude_edges(height)
    else:
        latitudes = np.asarray(source_latitudes, dtype=np.float64)
        if latitudes.size != height:
            raise ValueError("Latitude coordinate length does not match the source")
        if latitudes[0] < latitudes[-1]:
            source = source[::-1]
            latitudes = latitudes[::-1]
        latitude_source_edges = _latitude_edges_from_centres(latitudes)

    first_face = 0.0
    if source_longitudes is not None:
        longitudes = np.asarray(source_longitudes, dtype=np.float64)
        if longitudes.ndim != 1 or longitudes.size != width or not np.isfinite(longitudes).all():
            raise ValueError("Longitude coordinate length does not match the source")
        normalized = np.mod(longitudes + 180.0, 360.0)
        order = np.argsort(normalized)
        normalized = normalized[order]
        source = source[:, order]
        step = 360.0 / width
        tolerance = max(1.0e-6, 4 * np.finfo(np.float32).eps * 360.0)
        if not np.allclose(np.diff(normalized), step, rtol=0.0, atol=tolerance):
            raise ValueError("Source longitudes must be a regular global grid")
        first_face = float(normalized[0] - 0.5 * step)

    return source, latitude_source_edges, first_face


def conservative_remap(
    values: np.ndarray,
    width: int,
    height: int | None = None,
    *,
    source_latitudes: np.ndarray | None = None,
    source_longitudes: np.ndarray | None = None,
) -> np.ndarray:
    """Area-average a global field onto the canonical cell-centred grid."""
    width, height = grid_dimensions(width, height)
    source, source_edges, first_face = _normalize_source(
        values, source_latitudes, source_longitudes
    )
    source_height, source_width = source.shape
    longitude_plan = _longitude_plan(source_width, width, first_face)
    latitude_plan = _latitude_plan(source_edges, height)
    finite = np.isfinite(source)

    def apply(field: np.ndarray) -> np.ndarray:
        longitude = np.empty((source_height, width), dtype=np.float64)
        for x, (indices, weights) in enumerate(longitude_plan):
            longitude[:, x] = field[:, indices] @ weights
        destination = np.empty((height, width), dtype=np.float64)
        for y, (indices, weights) in enumerate(latitude_plan):
            destination[y] = weights @ longitude[indices]
        return destination

    numerator = apply(np.where(finite, source, 0.0))
    denominator = apply(finite.astype(np.float64))
    result = np.full((height, width), np.nan, dtype=np.float32)
    np.divide(numerator, denominator, out=result, where=denominator > 1.0e-12)
    return result


def spherical_area_weights(width: int, height: int | None = None) -> np.ndarray:
    width, height = grid_dimensions(width, height)
    bands = grid.band_area(height, degree_edges=True)
    return np.broadcast_to((bands * 2.0 * np.pi / width)[:, None], (height, width)).copy()
