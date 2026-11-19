"""Booleans against Manifold (``manifold3d``), an independent implementation used only as an oracle.

Agreement bounds (``docs/design/spec.md`` section 12, "Oracle suite"):

- ``delta = max(tol_manifold, 8 eps R)``: Manifold may move or collapse geometry by its tolerance
  (``get_tolerance()``, 1e-12 times the bounding-box scale); wireframes moves each output vertex by
  at most half an ulp of ``R``, the largest coordinate magnitude.
- Volume: ``|dV| <= 4 delta (area(a) + area(b))``. Moving a surface by ``delta`` changes the
  enclosed volume by at most ``delta`` times its area, and a result's surface is part of the two
  input surfaces; both sides move, plus summation slack.
- Area: ``|dA| <= 4 delta P`` with ``P`` the sum of face perimeters of both inputs and both results:
  moving a corner by ``delta`` changes a triangle's area by at most ``delta`` times its perimeter.
  For ``coplanar`` pairs (b is a turned by ``theta = 1e-9`` rad) a crossing point slides along the
  surface by up to ``delta / theta`` when a plane moves by ``delta``, so ``delta`` becomes
  ``delta / theta`` there; ``test_coplanar_crossings_are_exact`` checks wireframes' crossing points
  exactly instead.
- Manifold returns a zero-volume sheet where solids touch along a face (the regularized
  intersection is empty); components with ``|volume| <= delta * area`` are dropped from its result.
"""

import _oracle as o
import numpy as np
import pytest
import wireframes as wf

mf = pytest.importorskip("manifold3d")

pytestmark = pytest.mark.oracle

KINDS = o.CONVEX + o.NONCONVEX
PAIRS = 40
THETA = 1e-9
OPS = {
    "union": (wf.union, lambda a, b: a + b),
    "intersection": (wf.intersection, lambda a, b: a ^ b),
    "difference": (wf.difference, lambda a, b: a - b),
    "symmetric_difference": (wf.symmetric_difference, lambda a, b: (a - b) + (b - a)),
}


def to_manifold(m):
    out = mf.Manifold(
        mf.Mesh64(np.array(m.coords, np.float64, order="C"), np.array(m.triangles, np.uint64, order="C"))
    )
    assert out.status() == mf.Error.NoError
    return out


def solid_part(m, delta):
    """Volume, area and perimeter sum of the components of a Manifold result that hold volume."""
    vol = area = perimeter = 0.0
    dropped = 0
    for c in m.decompose():
        if abs(c.volume()) <= delta * c.surface_area():
            dropped += 1
            continue
        vol, area = vol + c.volume(), area + c.surface_area()
        g = c.to_mesh64()
        perimeter += o.edge_lengths(
            wf.Mesh(np.asarray(g.vert_properties)[:, :3], np.asarray(g.tri_verts).astype(np.int64))
        )
    return vol, area, perimeter, dropped


def compare(arrangement, seed):
    rng = np.random.default_rng(seed)
    tally = {
        "pairs": 0,
        "strict": 0,
        "comparisons": 0,
        "sheets": 0,
        "not_solid": 0,
        "errors": [],
        "dv": 0.0,
        "da": 0.0,
    }
    for _ in range(PAIRS):
        a, b = o.random_pair(rng, KINDS, KINDS, arrangement)
        assert a.is_solid and b.is_solid
        ma, mb = to_manifold(a), to_manifold(b)
        r = float(max(np.abs(a.coords).max(), np.abs(b.coords).max()))
        delta = max(ma.get_tolerance(), mb.get_tolerance(), 8 * o.EPS * r)
        slide = delta / THETA if arrangement == "coplanar" else delta
        rounding_case = (
            arrangement == "coplanar" or o.coplanar_only_to_rounding(a) + o.coplanar_only_to_rounding(b) > 0
        )
        tally["pairs"] += 1
        tally["strict"] += not rounding_case
        for name, (ours, theirs) in OPS.items():
            try:
                got = ours(a, b)
            except wf.KernelError as e:
                got = e.mesh
                assert rounding_case, f"{arrangement} {name}: {e}"
            ref = theirs(ma, mb)
            v, ar, per, dropped = solid_part(ref, delta)
            tally["sheets"] += dropped
            tally["comparisons"] += 1
            dv = abs(o.vol(got) - v)
            da = abs(o.area(got) - ar)
            bound_v = 4 * delta * (o.area(a) + o.area(b))
            bound_a = 4 * slide * (o.edge_lengths(a) + o.edge_lengths(b) + o.edge_lengths(got) + per)
            tally["dv"], tally["da"] = max(tally["dv"], dv), max(tally["da"], da)
            if dv > bound_v or da > bound_a:
                tally["errors"].append((name, dv, bound_v, da, bound_a))
            assert got.n_faces == 0 or got.is_closed
            if got.n_faces and not got.is_solid:
                tally["not_solid"] += 1
                assert rounding_case, f"{arrangement} {name}: not a solid without a rounding cause"
    return tally


@pytest.mark.parametrize(
    ("arrangement", "seed"),
    [("overlap", 1), ("nested", 2), ("touching", 3), ("utm", 4), ("far", 5), ("coplanar", 6)],
)
def test_boolean_volumes_and_areas_agree_with_manifold(arrangement, seed):
    t = compare(arrangement, seed)
    print(
        f"\noracle manifold {arrangement}: {t['pairs']} pairs, {t['comparisons']} op results compared, "
        f"{len(t['errors'])} outside the bound (largest |dV| {t['dv']:.1e}, |dA| {t['da']:.1e}), {t['sheets']} Manifold zero-volume sheets dropped, "
        f"{t['strict']} pairs without faces coplanar only to rounding (results must be solids), "
        f"{t['not_solid']} wireframes results refused with KernelError, closed but not solid (tolerance policy rule 5)"
    )
    assert t["errors"] == []


def test_manifold_keeps_a_flat_sheet_where_boxes_touch_along_a_face():
    """Recorded oracle disagreement: the regularized intersection of face-touching boxes is empty.

    Manifold returns a zero-volume double-sided sheet over the contact for one argument order and
    nothing for the other, so its intersection is not even commutative here.
    """
    a = wf.datasets.box([0, 0, 0], [1, 1, 1], length_unit=None)
    b = wf.datasets.box([-1, 0.25, 0.25], [0, 0.75, 0.5], length_unit=None)
    assert (a & b).n_faces == 0 and (b & a).n_faces == 0
    ab, ba = to_manifold(a) ^ to_manifold(b), to_manifold(b) ^ to_manifold(a)
    assert ab.volume() == 0.0 and ba.volume() == 0.0
    assert ab.surface_area() == pytest.approx(2 * 0.5 * 0.25, rel=1e-12)
    assert ba.surface_area() == 0.0


def test_coplanar_crossings_are_exact():
    """Crossings of solids turned 1e-9 rad apart: wireframes rounds the exact crossing points once
    (error at most half an ulp per coordinate); Manifold's points slide by up to ``delta / theta``.
    """
    rng = np.random.default_rng(7)
    worst_ours, worst_theirs = 0.0, 0.0
    for _ in range(3):
        a, b = o.random_pair(rng, ("hull",), ("hull",), "coplanar")
        got = wf.difference(a, b)
        ref = (to_manifold(a) - to_manifold(b)).to_mesh64()
        ours = o.crossing_error(got.coords, a, b)
        theirs = o.crossing_error(np.asarray(ref.vert_properties)[:, :3], a, b)
        assert len(ours) > 0
        half_ulp = np.spacing(np.abs(got.coords).max()) / 2
        assert ours.max() <= half_ulp
        worst_ours, worst_theirs = max(worst_ours, ours.max()), max(worst_theirs, theirs.max(initial=0.0))
    print(
        f"\noracle manifold coplanar crossings: wireframes worst {worst_ours:.1e}, Manifold worst {worst_theirs:.1e}"
    )
