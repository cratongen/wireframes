"""DTM and convex solids against scipy.spatial (Qhull), an independent implementation used only as an oracle.

Agreement bounds (``docs/design/spec.md`` section 12, "Oracle suite"):

- Delaunay triangles are decisions: the triangle sets must be equal. A Delaunay triangulation of
  points with no four cocircular is unique, so any difference means one side is not Delaunay; each
  one is settled by the exact in-circle test on every edge of the differing triangles
  (``_oracle.delaunay_violations``), which must find wireframes clean.
- Plan area of the DTM against the area of the convex hull: ``|dA| <= 8 n eps A`` relative, the
  rounding of ``n`` f64 triangle areas summed against Qhull's own sum.
- Convex intersections: ``vol(a & b)`` against Qhull's halfspace intersection of the face planes
  of both solids, ``|dV| <= 1e-12 (vol(a) + vol(b))``. Qhull solves each vertex from three planes
  in f64, an error of about ``eps R / sin^2(angle)`` where the planes meet at a small angle, so the
  bound leaves room above the largest relative difference seen in this suite (1.4e-16).
- The hull of the vertices of a convex result has the volume of the result:
  ``|dV| <= 8 n eps R^3`` with ``n`` the result's faces and ``R`` its largest coordinate.
"""

import _oracle as o
import numpy as np
import pytest
import wireframes as wf

spatial = pytest.importorskip("scipy.spatial")

pytestmark = pytest.mark.oracle


def plan_triangles(xy, tris):
    return {tuple(sorted(map(tuple, xy[t].tolist()))) for t in tris}


def index_of(xy):
    return {tuple(p): i for i, p in enumerate(np.asarray(xy).tolist())}


@pytest.mark.parametrize("utm", [False, True])
def test_dtm_is_the_delaunay_triangulation_of_scipy(utm):
    compared = scipy_wrong = 0
    for seed in range(30):
        rng = np.random.default_rng(seed)
        n = int(rng.integers(10, 400))
        off = wf.datasets.ORIGIN[:2] if utm else np.zeros(2)
        xy = off + rng.uniform(-500, 500, (n, 2))
        m = wf.dtm(np.column_stack([xy, rng.uniform(0, 100, n)]))
        assert np.array_equal(np.sort(m.coords[:, :2], axis=0), np.sort(xy, axis=0))
        d = spatial.Delaunay(xy)
        ours, theirs = plan_triangles(m.coords[:, :2], m.triangles), plan_triangles(xy, d.simplices)
        compared += len(ours)
        assert o.delaunay_violations(m.coords[:, :2], m.triangles) == []
        if ours != theirs:
            ix = index_of(xy)
            only = [f for f, t in enumerate(d.simplices) if plan_triangles(xy, [t]) <= theirs - ours]
            assert o.delaunay_violations(xy, d.simplices, only=set(only)), (
                f"seed {seed}: two Delaunay answers"
            )
            assert all(all(tuple(p) in ix for p in t) for t in ours - theirs)
            scipy_wrong += 1
        hull = spatial.ConvexHull(xy - off)
        q = m.coords[m.triangles][:, :, :2] - off
        u, w = q[:, 1] - q[:, 0], q[:, 2] - q[:, 0]
        plan = 0.5 * np.abs(u[:, 0] * w[:, 1] - u[:, 1] * w[:, 0]).sum()
        assert abs(plan - hull.volume) <= 8 * m.n_faces * o.EPS * hull.volume
    print(
        f"\noracle scipy delaunay {'utm' if utm else 'origin'}: 30 point sets, {compared} triangles compared, "
        f"{scipy_wrong} sets where Qhull returned a non-Delaunay edge (settled by exact in-circle)"
    )


def test_dtm_of_cocircular_grid_is_delaunay_with_the_hull_area():
    """A grid has four cocircular points in every cell, so the Delaunay triangulation is not unique;
    both sides must be Delaunay (no point strictly inside a circumcircle) and cover the hull."""
    g = np.arange(12.0)
    xy = np.column_stack([a.ravel() for a in np.meshgrid(g, g)]) + wf.datasets.ORIGIN[:2]
    m = wf.dtm(np.column_stack([xy, np.zeros(len(xy))]))
    d = spatial.Delaunay(xy)
    assert o.delaunay_violations(m.coords[:, :2], m.triangles) == []
    assert o.delaunay_violations(xy, d.simplices) == []
    assert m.area == 121.0 and m.n_faces == len(d.simplices) == 242


def halfspaces(m):
    v, t = m.coords, m.triangles
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1)[:, None]
    return np.column_stack([n, -np.einsum("ij,ij->i", n, a)])


def test_convex_intersections_agree_with_qhull_halfspaces():
    rng = np.random.default_rng(21)
    compared = 0
    for _ in range(60):
        a, b = o.random_pair(rng, ("hull", "sphere", "box"), ("hull", "sphere", "box"), "overlap")
        got = a & b
        if got.n_faces == 0:
            continue
        inner = got.coords.mean(axis=0)
        assert got.locate(inner[None])[0] == -1
        h = spatial.HalfspaceIntersection(np.vstack([halfspaces(a), halfspaces(b)]), inner)
        ref = spatial.ConvexHull(h.intersections).volume
        assert abs(got.volume - ref) <= 1e-12 * (a.volume + b.volume)
        hull = spatial.ConvexHull(got.coords).volume
        r = float(np.abs(got.coords).max())
        assert abs(got.volume - hull) <= 8 * got.n_faces * o.EPS * r**3
        compared += 1
    print(
        f"\noracle scipy convex: {compared} intersections compared with Qhull halfspaces and hulls, 0 outside the bound"
    )
