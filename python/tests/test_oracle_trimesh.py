"""Measures, containment and ray hits against trimesh, an independent implementation used only as an oracle.

Agreement bounds (``docs/design/spec.md`` section 12, "Oracle suite"), with ``n`` faces, ``R`` the
largest coordinate magnitude and ``eps`` the f64 machine epsilon:

- Volume ``|dV| <= 8 n eps (sqrt(3) R)^3``, area ``|dA| <= 8 n eps (2 sqrt(3) R)^2`` and centroid
  ``|dc| <= 8 n eps (sqrt(3) R)^4 / V``: trimesh sums one f64 term per face from raw coordinates,
  each term bounded by the power of the corner norms shown and rounded with a few eps; the
  constant 8 covers the per-term rounding and numpy's pairwise summation. wireframes subtracts the
  bounding-box minimum first, so its own error is smaller.
- Containment and ray hits are decisions; a third of the points share two coordinates with a vertex and
  half the rays aim at a vertex, where ties are hardest. Decisions must agree exactly. Every disagreement
  is settled by an exact rational test (``_oracle.exact_inside``, ``_oracle.exact_ray_hits``), which
  must side with wireframes.
- Hit distances: ``|dt| |d| <= 1e-9 R`` after both sides agree on the number of hits. trimesh solves
  each hit in f64 barycentrics, whose error grows as the ray grazes a face, so its points carry more
  than rounding error; wireframes rounds exact hits once.
"""

import _oracle as o
import numpy as np
import pytest
import wireframes as wf

trimesh = pytest.importorskip("trimesh")

pytestmark = pytest.mark.oracle

KINDS = o.CONVEX + o.NONCONVEX


def to_trimesh(m):
    return trimesh.Trimesh(np.asarray(m.coords), np.asarray(m.triangles), process=False)


def solids(seed, n):
    """Seeded shapes of every kind, turned at random, boolean results of overlapping pairs, and
    integer boxes whose vertices and edges line up with the coordinate axes."""
    rng = np.random.default_rng(seed)
    out = [
        wf.datasets.box([0, 0, 0], [2, 3, 1], length_unit=None)
        | wf.datasets.box([1, 1, 0], [3, 2, 2], length_unit=None)
    ]
    for k in range(n):
        m = o.rigid(o.shape(KINDS[k % len(KINDS)], rng), o.rotation(rng), rng.uniform(-0.5, 0.5, 3))
        out.append(m)
        a, b = o.random_pair(rng, KINDS, KINDS, "overlap")
        try:
            r = [wf.union, wf.intersection, wf.difference][k % 3](a, b)
        except wf.KernelError as e:
            r = e.mesh
        if r.n_faces:
            out.append(r)
    return out


def test_volume_area_centroid_agree_with_trimesh():
    shapes = solids(11, 28)
    for m in shapes:
        t = to_trimesh(m)
        n, r = m.n_faces, float(np.abs(m.coords).max())
        assert abs(m.volume - t.volume) <= 8 * n * o.EPS * (np.sqrt(3) * r) ** 3
        assert abs(m.area - t.area) <= 8 * n * o.EPS * (2 * np.sqrt(3) * r) ** 2
        assert np.abs(m.centroid - t.center_mass).max() <= 8 * n * o.EPS * (np.sqrt(3) * r) ** 4 / m.volume
    print(
        f"\noracle trimesh measures: {len(shapes)} solids, {3 * len(shapes)} measures compared, 0 outside the bound"
    )


def test_trimesh_centroid_drifts_at_utm_coordinates():
    """Recorded oracle disagreement: a box's centroid is its centre (closed form). At UTM
    coordinates trimesh's centroid is off by 3e-4 m in northing; wireframes' is exact."""
    lo = wf.datasets.ORIGIN + [0.123, 0.456, 300.7]
    m = wf.datasets.box(lo, lo + [10, 20, 5])
    centre = lo + [5, 10, 2.5]
    assert np.abs(m.centroid - centre).max() <= np.spacing(centre).max()
    assert np.abs(to_trimesh(m).center_mass - centre).max() > 1e-5


def settle_inside(m, p):
    for d in [(0.2718281828, 0.3141592653, 0.9), (-0.577, 0.211, -0.3), (0.1, -0.9, 0.05)]:
        inside = o.exact_inside(m.coords, m.triangles, p, d)
        if inside is not None:
            return inside
    raise AssertionError(f"every settling ray grazes an edge from {p}")


def test_containment_agrees_with_trimesh_or_exact_parity():
    rng = np.random.default_rng(12)
    shapes = solids(12, 14)
    compared = trimesh_wrong = 0
    for m in shapes:
        lo, hi = (np.array(b) for b in m.bounds)
        p = rng.uniform(lo - 0.2, hi + 0.2, (1500, 3))
        v = m.coords[rng.integers(0, m.n_vertices, 500)]
        axis = rng.integers(0, 3, 500)
        v[np.arange(500), axis] += rng.uniform(-0.5, 0.5, 500)
        p[:500] = v
        loc = m.locate(p)
        p = p[loc != 0]
        ours = m.contains(p)
        theirs = to_trimesh(m).contains(p)
        compared += len(p)
        for i in np.flatnonzero(ours != theirs):
            assert settle_inside(m, p[i]) == ours[i], f"wireframes wrong at {p[i]!r}"
            trimesh_wrong += 1
    print(
        f"\noracle trimesh containment: {len(shapes)} solids, {compared} points compared, "
        f"{trimesh_wrong} disagreements, all settled exactly in favour of wireframes"
    )


def distinct(ts, scale):
    out = []
    for t in sorted(ts):
        if not out or t - out[-1] > 1e-9 * scale:
            out.append(t)
    return out


def test_ray_hits_agree_with_trimesh_or_exact_hits():
    rng = np.random.default_rng(13)
    shapes = solids(13, 14)
    rays = count_wrong = t_compared = 0
    for m in shapes:
        lo, hi = (np.array(b) for b in m.bounds)
        r = float(np.abs(m.coords).max())
        origins = rng.uniform(lo - 1, hi + 1, (200, 3))
        targets = rng.uniform(lo, hi, (200, 3))
        targets[::2] = m.coords[rng.integers(0, m.n_vertices, 100)]
        dirs = targets - origins
        table = m.ray_hits(origins, dirs).to_dict()
        tm = to_trimesh(m)
        loc, ray, _ = tm.ray.intersects_location(origins, dirs, multiple_hits=True)
        for k in range(len(origins)):
            d = dirs[k]
            ours = sorted(table["t"][table["ray"] == k].tolist())
            theirs_t = ((loc[ray == k] - origins[k]) @ d) / (d @ d)
            theirs = distinct(theirs_t[theirs_t >= 0].tolist(), r / np.linalg.norm(d))
            rays += 1
            if len(ours) != len(theirs):
                exact = o.exact_ray_hits(m.coords, m.triangles, origins[k], d)
                assert len(exact) == len(ours), f"wireframes wrong on ray {k}"
                assert all(
                    abs(float(e) - t) <= np.spacing(max(1.0, abs(float(e)))) for e, t in zip(exact, ours)
                )
                count_wrong += 1
                continue
            t_compared += len(ours)
            assert np.all(np.abs(np.subtract(ours, theirs)) * np.linalg.norm(d) <= 1e-9 * r)
    print(
        f"\noracle trimesh rays: {len(shapes)} solids, {rays} rays, {t_compared} hit distances compared, "
        f"{count_wrong} hit-count disagreements, all settled exactly in favour of wireframes"
    )
