"""Seeded synthetic geometry for examples and tests; no vendor or client data.

Every function is deterministic for a given seed. Coordinates sit at UTM scale (easting near 700 000,
northing near 7 500 000, metres) so examples exercise the precision rules of real projects.
"""

from __future__ import annotations

import numpy as np

from ._drusa import Mesh, PointSet, Polylines, dtm, loft

ORIGIN = np.array([700_000.0, 7_500_000.0, 0.0])

_BOX_FACES = np.array(
    [[0, 2, 1], [1, 2, 3], [4, 5, 6], [5, 7, 6], [0, 1, 4], [1, 5, 4],
     [2, 6, 3], [3, 6, 7], [0, 4, 2], [2, 4, 6], [1, 3, 5], [3, 7, 5]]
)  # fmt: skip


def box(lo, hi, *, name: str | None = None, length_unit: str | None = "m") -> Mesh:
    """An axis-aligned box solid with 12 outward triangles.

    Parameters
    ----------
    lo, hi : array_like
        Opposite corners, ``lo < hi`` on every axis.
    name : str, optional
        Value of a face column ``name`` on every triangle.
    length_unit : str, optional
        Length unit label.

    Returns
    -------
    Mesh
        A solid of volume ``prod(hi - lo)``.
    """
    (x0, y0, z0), (x1, y1, z1) = np.asarray(lo, float), np.asarray(hi, float)
    corners = [[x, y, z] for z in (z0, z1) for y in (y0, y1) for x in (x0, x1)]
    faces = {"name": [name] * 12} if name is not None else None
    return Mesh(corners, _BOX_FACES, face_attributes=faces, length_unit=length_unit)


def orebody(seed: int = 0, *, sections: int = 9, points: int = 24, length: float = 400.0) -> Mesh:
    """A lens-shaped orebody lofted from noisy elliptical cross sections.

    The lens strikes east, centred at ``ORIGIN + (0, 0, 200)``. Section ``k`` lies in the plane
    ``x = const`` and its semi-axes taper towards both ends.

    Parameters
    ----------
    seed : int
        Seed of the radius noise.
    sections : int
        Number of cross sections, at least 2.
    points : int
        Vertices per section ring, at least 3.
    length : float
        Strike length in metres.

    Returns
    -------
    Mesh
        A validated solid with face column ``section``.
    """
    rng = np.random.default_rng(seed)
    xs = np.linspace(-length / 2, length / 2, sections)
    angles = np.linspace(0.0, 2 * np.pi, points, endpoint=False)
    rings = []
    for x in xs:
        taper = np.sqrt(max(1.0 - (x / (0.55 * length)) ** 2, 0.05))
        r = 1.0 + 0.08 * rng.standard_normal(points)
        y = 60.0 * taper * r * np.cos(angles)
        z = 25.0 * taper * r * np.sin(angles) + 200.0
        rings.append(ORIGIN + np.column_stack([np.full(points, x), y, z]))
    return loft(Polylines(rings, closed=True, length_unit="m"))


def topography(seed: int = 0, *, n: int = 400, size: float = 800.0) -> Mesh:
    """A rolling topographic surface triangulated from seeded survey points.

    Parameters
    ----------
    seed : int
        Seed of the point positions and the relief.
    n : int
        Number of survey points inside the square, plus its four corners.
    size : float
        Side of the square, centred on ``ORIGIN``.

    Returns
    -------
    Mesh
        A height-field surface with vertex column ``z_noise``.
    """
    rng = np.random.default_rng(seed)
    half = size / 2
    xy = rng.uniform(-half, half, (n, 2))
    xy = np.vstack([xy, [[-half, -half], [half, -half], [half, half], [-half, half]]])
    phase = rng.uniform(0, 2 * np.pi, 2)
    noise = 2.0 * rng.standard_normal(len(xy))
    z = 300.0 + 15.0 * np.sin(xy[:, 0] / 150 + phase[0]) + 10.0 * np.cos(xy[:, 1] / 110 + phase[1]) + noise
    pts = ORIGIN + np.column_stack([xy, z])
    surface = dtm(PointSet(pts, length_unit="m"))
    return Mesh(
        surface.coords,
        surface.triangles,
        vertex_attributes={"z_noise": noise},
        length_unit="m",
    )


def pit(
    ground: Mesh | None = None,
    *,
    floor: float = 180.0,
    radius: float = 60.0,
    slope: float = 45.0,
    size: float = 800.0,
    n: int = 41,
) -> Mesh:
    """A conical pit design surface on a regular plan grid: a flat floor and walls at `slope` degrees.

    Outside the pit the design follows `ground`, so the design never rises above it.

    Parameters
    ----------
    ground : Mesh, optional
        Surface the walls daylight on; None uses ``topography()``.
    floor : float
        Floor elevation.
    radius : float
        Floor radius in plan.
    slope : float
        Wall angle from horizontal in degrees.
    size : float
        Side of the square, centred on ``ORIGIN``.
    n : int
        Grid nodes per side.

    Returns
    -------
    Mesh
        A height-field surface.
    """
    g = np.linspace(-size / 2, size / 2, n)
    gx, gy = np.meshgrid(g, g)
    ground = topography() if ground is None else ground
    plan = ORIGIN[:2] + np.column_stack([gx.ravel(), gy.ravel()])
    cone = floor + np.maximum(np.hypot(gx, gy).ravel() - radius, 0.0) * np.tan(np.radians(slope))
    z = np.fmin(cone, ground.surface_z(plan))
    return dtm(PointSet(np.column_stack([plan, z]), length_unit="m"))


def drillholes(seed: int = 0, *, n: int = 12, depth: float = 350.0) -> Polylines:
    """Straight inclined drillhole traces collared on a line across the orebody.

    Parameters
    ----------
    seed : int
        Seed of the collar spacing, azimuth and dip.
    n : int
        Number of holes.
    depth : float
        Length of each hole.

    Returns
    -------
    Polylines
        One open two-vertex part per hole, with feature columns ``hole`` (str) and ``depth``.
    """
    rng = np.random.default_rng(seed)
    xs = np.linspace(-170.0, 170.0, n) + rng.uniform(-10, 10, n)
    parts = []
    for x in xs:
        dip = np.radians(rng.uniform(55, 75))
        azimuth = np.radians(rng.uniform(175, 185))
        collar = ORIGIN + [x, 60.0 + rng.uniform(-10, 10), 320.0]
        direction = np.array([np.sin(azimuth) * np.cos(dip), np.cos(azimuth) * np.cos(dip), -np.sin(dip)])
        parts.append(np.vstack([collar, collar + depth * direction]))
    names = [f"DH{i + 1:03d}" for i in range(n)]
    return Polylines(parts, attributes={"hole": names, "depth": np.full(n, depth)}, length_unit="m")
