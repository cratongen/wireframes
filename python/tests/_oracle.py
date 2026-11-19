"""Seeded solids and exact tie-breakers shared by the oracle tests (``test_oracle_*.py``).

The oracles (Manifold through ``manifold3d``, ``trimesh`` and ``scipy.spatial``) are dev-only
dependencies; every oracle test skips when its package is missing.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import wireframes as wf

EPS = np.finfo(float).eps


def rotation(rng, angle=None):
    """A random rotation matrix (uniform axis; uniform angle unless given)."""
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    t = rng.uniform(0, 2 * np.pi) if angle is None else angle
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * (k @ k)


def rigid(mesh, r=None, t=(0.0, 0.0, 0.0)):
    """``mesh`` rotated by ``r`` about the origin, then translated by ``t``."""
    m = np.eye(4)
    if r is not None:
        m[:3, :3] = r
    m[:3, 3] = t
    return mesh.transform(m)


def from_arrays(vertices, faces):
    """A Mesh from arrays, keeping only the vertices the faces use."""
    faces = np.asarray(faces, np.int64)
    used, inverse = np.unique(faces, return_inverse=True)
    return wf.Mesh(np.asarray(vertices, float)[used], inverse.reshape(-1, 3))


def outward(vertices, faces):
    """Faces of a convex polytope wound outward (normal away from the vertex mean)."""
    v, f = np.asarray(vertices, float), np.array(faces, np.int64)
    c = v.mean(axis=0)
    a, b, d = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    flip = np.einsum("ij,ij->i", np.cross(b - a, d - a), a - c) < 0
    f[flip] = f[flip][:, ::-1]
    return f


def hull_solid(rng, n=24, radius=1.0, center=(0.0, 0.0, 0.0)):
    """Convex hull of ``n`` random points near a sphere: a convex solid."""
    from scipy.spatial import ConvexHull

    p = rng.normal(size=(n, 3))
    p *= (radius * rng.uniform(0.7, 1.0, n) / np.linalg.norm(p, axis=1))[:, None]
    p += center
    h = ConvexHull(p)
    return from_arrays(p, outward(p, h.simplices))


def icosphere(radius=1.0, subdivisions=2):
    import trimesh

    t = trimesh.creation.icosphere(subdivisions, radius)
    return wf.Mesh(t.vertices, t.faces)


def torus(major=1.0, minor=0.35):
    import trimesh

    t = trimesh.creation.torus(major, minor, major_sections=20, minor_sections=10)
    return wf.Mesh(t.vertices, t.faces)


def star_prism(rng, points=7, height=1.0):
    """A non-convex star polygon extruded along z, base centred on the origin."""
    k = np.arange(2 * points)
    r = np.where(k % 2 == 0, 1.0, rng.uniform(0.35, 0.6))
    t = np.pi * k / points + rng.uniform(0, np.pi)
    ring = np.column_stack([r * np.cos(t), r * np.sin(t), np.full(len(k), -height / 2)])
    return wf.extrude(wf.Polylines([ring], closed=True), [0.0, 0.0, height])


def l_prism(height=1.0):
    ring = [[0, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [1, 2, 0], [0, 2, 0]]
    m = wf.extrude(wf.Polylines([np.array(ring, float) - [1, 1, height / 2]], closed=True), [0, 0, height])
    return m


def lens(seed):
    """The lofted orebody of ``wf.datasets.orebody`` scaled to unit size at the origin."""
    o = wf.datasets.orebody(seed, sections=5, points=10)
    lo, hi = (np.array(b) for b in o.bounds)
    s = 2.0 / (hi - lo).max()
    m = np.diag([s, s, s, 1.0])
    m[:3, 3] = -s * (lo + hi) / 2
    return o.transform(m)


def shape(kind, rng):
    if kind == "box":
        return wf.datasets.box(-rng.uniform(0.5, 1, 3), rng.uniform(0.5, 1, 3), length_unit=None)
    if kind == "hull":
        return hull_solid(rng, n=int(rng.integers(8, 40)))
    if kind == "sphere":
        return icosphere(rng.uniform(0.6, 1.0))
    if kind == "torus":
        return torus()
    if kind == "star":
        return star_prism(rng)
    if kind == "l":
        return l_prism()
    if kind == "lens":
        return lens(int(rng.integers(0, 1000)))
    raise ValueError(kind)


CONVEX = ("box", "hull", "sphere")
NONCONVEX = ("torus", "star", "l", "lens")


def random_pair(rng, kinds_a, kinds_b, arrangement):
    """A seeded pair of solids in one of the arrangements the oracle suite covers.

    ``overlap``: b rotated at random and translated within the reach of a.
    ``nested``: b is a with coordinates halved about the origin: nested in a when a is star-shaped
    about the origin (box, hull, sphere, star, lens), touching or crossing a otherwise (torus, L).
    ``coplanar``: b is a copy of a rotated by 1e-9 rad about a random axis and moved 1e-9.
    ``utm``: an ``overlap`` pair moved to ``wf.datasets.ORIGIN`` (easting 7e5, northing 7.5e6).
    ``far``: b translated well away from a (disjoint).
    ``touching``: a's bounding box and a box b that touches it along part of a face (1 shared
    coordinate), an edge (2) or a corner (3), with exactly equal coordinates.
    """
    a = shape(rng.choice(kinds_a), rng)
    if arrangement == "nested":
        b = a.transform(np.diag([0.5, 0.5, 0.5, 1.0]))
        return a, b
    if arrangement == "coplanar":
        return a, rigid(a, rotation(rng, 1e-9), rng.normal(size=3) * 1e-9)
    if arrangement == "touching":
        lo = np.array(a.bounds[0])
        hi = np.array(a.bounds[1])
        size = rng.uniform(0.5, 1.5, 3)
        contact = rng.integers(1, 4)
        start = lo + rng.uniform(0, 1, 3) * (hi - lo)
        axes = rng.permutation(3)[:contact]
        start[axes] = hi[axes]
        return wf.datasets.box(lo, hi, length_unit=None), wf.datasets.box(
            start, start + size, length_unit=None
        )
    b = shape(rng.choice(kinds_b), rng)
    if arrangement == "overlap":
        return a, rigid(b, rotation(rng), rng.uniform(-0.8, 0.8, 3))
    if arrangement == "utm":
        return rigid(a, t=wf.datasets.ORIGIN), rigid(
            b, rotation(rng), wf.datasets.ORIGIN + rng.uniform(-0.8, 0.8, 3)
        )
    if arrangement == "far":
        return a, rigid(b, rotation(rng), [5.0, 0, 0])
    raise ValueError(arrangement)


def edge_lengths(mesh):
    """Sum of every face's perimeter."""
    v, f = mesh.coords, mesh.triangles
    if len(f) == 0:
        return 0.0
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    return float(sum(np.linalg.norm(p - q, axis=1).sum() for p, q in ((a, b), (b, c), (c, a))))


def vol(m):
    return 0.0 if m.n_faces == 0 else m.volume


def area(m):
    return 0.0 if m.n_faces == 0 else m.area


def _q(p):
    return [Fraction(float(x)) for x in p]


def _sub(a, b):
    return [x - y for x, y in zip(a, b, strict=True)]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def exact_ray_hits(vertices, faces, origin, direction):
    """Exact rational hits of the ray ``origin + t * direction`` (t >= 0) with the triangles.

    Returns the sorted list of distinct exact ``t`` (Fractions) where the ray meets a closed
    triangle that it crosses (not lies in). Used to settle a disagreement, never as the oracle.
    """
    o, d = _q(origin), _q(direction)
    ts = set()
    for f in faces:
        a, b, c = (_q(vertices[i]) for i in f)
        n = _cross(_sub(b, a), _sub(c, a))
        den = _dot(n, d)
        if den == 0:
            continue
        t = _dot(n, _sub(a, o)) / den
        if t < 0:
            continue
        p = [oi + t * di for oi, di in zip(o, d, strict=True)]
        s = [_dot(_cross(_sub(v1, v0), _sub(p, v0)), n) for v0, v1 in ((a, b), (b, c), (c, a))]
        if all(x >= 0 for x in s) or all(x <= 0 for x in s):
            ts.add(t)
    return sorted(ts)


def exact_inside(vertices, faces, point, direction=(0.2718281828, 0.3141592653, 0.9)):
    """Exact parity test along a ray from a point off the surface; ``None`` when the ray meets an edge
    or vertex, or runs inside a face's plane, where parity is not decided."""
    o, d = _q(point), _q(direction)
    count = 0
    for f in faces:
        a, b, c = (_q(vertices[i]) for i in f)
        n = _cross(_sub(b, a), _sub(c, a))
        den = _dot(n, d)
        if den == 0:
            continue
        t = _dot(n, _sub(a, o)) / den
        if t <= 0:
            continue
        p = [oi + t * di for oi, di in zip(o, d, strict=True)]
        s = [_dot(_cross(_sub(v1, v0), _sub(p, v0)), n) for v0, v1 in ((a, b), (b, c), (c, a))]
        if any(x == 0 for x in s) and (all(x >= 0 for x in s) or all(x <= 0 for x in s)):
            return None
        if all(x > 0 for x in s) or all(x < 0 for x in s):
            count += 1
    return count % 2 == 1


def exact_incircle(a, b, c, d):
    """Sign of the in-circle determinant of plan points a, b, c (counter-clockwise) and d."""
    rows = []
    for p in (a, b, c):
        dx, dy = Fraction(float(p[0])) - Fraction(float(d[0])), Fraction(float(p[1])) - Fraction(float(d[1]))
        rows.append((dx, dy, dx * dx + dy * dy))
    (a0, a1, a2), (b0, b1, b2), (c0, c1, c2) = rows
    det = a0 * (b1 * c2 - b2 * c1) - a1 * (b0 * c2 - b2 * c0) + a2 * (b0 * c1 - b1 * c0)
    return (det > 0) - (det < 0)


def _orient2d(a, b, c):
    a, b, c = (tuple(Fraction(float(x)) for x in p[:2]) for p in (a, b, c))
    d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (d > 0) - (d < 0)


def delaunay_violations(xy, triangles, only=None):
    """Interior edges whose opposite vertex lies strictly inside the other triangle's circumcircle.

    Exact on the ``f64`` plan coordinates. ``only``: restrict to edges of these triangle indices.
    """
    xy = np.asarray(xy, float)
    owner = {}
    for f, t in enumerate(np.asarray(triangles)):
        for k in range(3):
            owner.setdefault(frozenset((int(t[k]), int(t[(k + 1) % 3]))), []).append((f, int(t[(k + 2) % 3])))
    bad = []
    for e, fs in owner.items():
        if len(fs) != 2 or (only is not None and fs[0][0] not in only and fs[1][0] not in only):
            continue
        (f, c), (_, d) = fs
        p, q = sorted(e)
        a, b = (p, q) if _orient2d(xy[p], xy[q], xy[c]) > 0 else (q, p)
        if exact_incircle(xy[a], xy[b], xy[c], xy[d]) > 0:
            bad.append(tuple(sorted(e)))
    return bad


def orient3d_sign(a, b, c, d):
    """Exact sign of the orientation of d against the plane through a, b, c."""
    a, b, c, d = (_q(p) for p in (a, b, c, d))
    x = _dot(_cross(_sub(b, a), _sub(c, a)), _sub(d, a))
    return (x > 0) - (x < 0)


def coplanar_only_to_rounding(mesh):
    """Number of edges whose two faces are not exactly coplanar but are within rounding of it.

    Such a pair is what a quad becomes when its corners are rotated or computed in ``f64``: the
    case in which the spec's tolerance policy (rule 5) says rounding can leave self-intersections.
    """
    v, t = mesh.coords, mesh.triangles
    scale = 64 * EPS * max(1.0, float(np.abs(v).max()))
    owner = {}
    for f, tri in enumerate(t):
        for k in range(3):
            owner.setdefault(frozenset((int(tri[k]), int(tri[(k + 1) % 3]))), []).append(f)
    count = 0
    for fs in owner.values():
        if len(fs) != 2:
            continue
        f, g = fs
        d = next(int(i) for i in t[g] if i not in t[f])
        a, b, c = v[t[f]]
        n = np.cross(b - a, c - a)
        if abs(n @ (v[d] - a)) <= scale * np.linalg.norm(n) and orient3d_sign(a, b, c, v[d]) != 0:
            count += 1
    return count


def _edges(mesh):
    t = mesh.triangles
    return np.unique(np.sort(np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1), axis=0)


def _segment_distance(p, s0, s1):
    d = s1 - s0
    t = np.clip(((p - s0) * d).sum(1) / (d * d).sum(1), 0, 1)
    return np.linalg.norm(s0 + t[:, None] * d - p, axis=1)


def crossing_error(vertices, a, b, near=1e-6, edges=4, faces=12):
    """Distance of each intersection vertex of a boolean result to the exact crossing it stands for.

    For every vertex of ``vertices`` that is not an input vertex, the exact rational crossing
    points of the input edges and faces nearest to it (an edge of one input within ``near``, one of
    the ``faces`` nearest faces of the other) are computed with Fractions, and the distance to the
    closest one is returned. Vertices with no crossing nearby are skipped.
    """
    inputs = {tuple(p) for p in np.vstack([a.coords, b.coords]).tolist()}
    sides = [(x, _edges(x), y, y.coords[y.triangles].mean(axis=1)) for x, y in ((a, b), (b, a))]
    out = []
    for v in np.asarray(vertices, float):
        if tuple(v.tolist()) in inputs:
            continue
        best = None
        for x, ex, y, cy in sides:
            dist = _segment_distance(v, x.coords[ex[:, 0]], x.coords[ex[:, 1]])
            near_faces = np.argsort(np.linalg.norm(cy - v, axis=1))[:faces]
            for e in np.argsort(dist)[:edges]:
                if dist[e] > near:
                    continue
                p0, p1 = (_q(x.coords[i]) for i in ex[e])
                d = _sub(p1, p0)
                for f in near_faces:
                    fa, fb, fc = (_q(y.coords[i]) for i in y.triangles[f])
                    n = _cross(_sub(fb, fa), _sub(fc, fa))
                    den = _dot(n, d)
                    if den == 0:
                        continue
                    t = _dot(n, _sub(fa, p0)) / den
                    if not 0 <= t <= 1:
                        continue
                    p = [s + t * w for s, w in zip(p0, d, strict=True)]
                    s = [
                        _dot(_cross(_sub(q1, q0), _sub(p, q0)), n)
                        for q0, q1 in ((fa, fb), (fb, fc), (fc, fa))
                    ]
                    if all(z >= 0 for z in s):
                        err = float(max(abs(Fraction(float(c)) - e_) for c, e_ in zip(v, p, strict=True)))
                        best = err if best is None else min(best, err)
        if best is not None:
            out.append(best)
    return np.array(out)
