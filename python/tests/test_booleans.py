import drusa as dr
import numpy as np
import pytest
from conftest import close


def pairs(box, rng, n=6):
    for _ in range(n):
        lo1 = rng.integers(-3, 3, 3).astype(float) + rng.uniform(0, 0.5, 3)
        lo2 = lo1 + rng.uniform(-1.5, 1.5, 3)
        yield box(lo1, lo1 + rng.uniform(1, 3, 3), name="a"), box(lo2, lo2 + rng.uniform(1, 3, 3), name="b")


def pieces_touch_only_along_curves(m):
    """``a ^ b`` is ``a - b`` and ``b - a`` in one mesh; the two pieces meet along the intersection curve,
    so validation finds duplicate vertices there and nothing else."""
    s = m.validate().summary
    others = {
        k: v for k, v in s.items() if k not in ("duplicate_vertices", "shells", "is_closed", "is_solid")
    }
    return s["is_closed"] and not any(others.values())


def vol(m):
    """Volume, 0 for an empty result (nothing kept)."""
    return 0.0 if m.n_faces == 0 else m.volume


def test_volume_identity_and_closed_output(box, rng):
    for a, b in pairs(box, rng):
        u, i, d, x = a | b, a & b, a - b, a ^ b
        for m in (u, i, d):
            assert m.n_faces == 0 or m.is_solid
        assert x.n_faces == 0 or pieces_touch_only_along_curves(x)
        assert close(vol(u) + vol(i), vol(a) + b.volume)
        assert close(vol(d), vol(a) - vol(i))
        assert close(vol(x), vol(u) - vol(i))
        assert close(vol(b | a), vol(u)) and close(vol(b & a), vol(i))


def test_methods_and_functions_agree(box):
    a, b = box([0, 0, 0], [2, 2, 2]), box([1, 1, 1], [3, 3, 3])
    for name, op in [("union", "__or__"), ("intersection", "__and__"), ("difference", "__sub__"),
                     ("symmetric_difference", "__xor__")]:  # fmt: skip
        v = getattr(dr, name)(a, b).volume
        assert getattr(a, name)(b).volume == v == getattr(a, op)(b).volume
    assert (a | b).volume == pytest.approx(15.0)
    assert (a & b).volume == pytest.approx(1.0)
    assert (a - b).volume == pytest.approx(7.0)
    assert (a ^ b).volume == pytest.approx(14.0)


def test_self_operations(orebody):
    assert close((orebody | orebody).volume, orebody.volume)
    assert close((orebody & orebody).volume, orebody.volume)
    assert (orebody - orebody).n_faces == 0


def test_orebody_against_a_block_at_utm_scale(orebody, box):
    o = dr.datasets.ORIGIN
    block = box(o + [-50, -80, 150], o + [50, 80, 250], name="block")
    inside, outside, union = orebody & block, orebody - block, orebody | block
    assert close(inside.volume + outside.volume, orebody.volume, rel=1e-8)
    assert close(union.volume + inside.volume, orebody.volume + block.volume, rel=1e-8)
    from_orebody = inside.face_attributes.mask("section")
    assert from_orebody.any() and not from_orebody.all()
    assert set(inside["name"][~from_orebody].tolist()) == {"block"}
    assert set(inside["name"][from_orebody].tolist()) == {None}


def test_locate_agrees_with_boolean_of_locates(box, rng):
    a, b = box([0, 0, 0], [2, 2, 2]), box([1, 0.5, 0.5], [3, 2.5, 2.5])
    p = rng.uniform(-0.5, 3.5, (300, 3))
    la, lb = a.locate(p), b.locate(p)
    off = (la != 0) & (lb != 0)
    ia, ib = la == -1, lb == -1
    for result, expected in [(a | b, ia | ib), (a & b, ia & ib), (a - b, ia & ~ib), (b - a, ib & ~ia)]:
        got = result.locate(p[off]) == -1
        assert np.array_equal(got, expected[off])
    x = a ^ b
    assert pieces_touch_only_along_curves(x)
    assert x.volume == pytest.approx((a - b).volume + (b - a).volume)


def test_coplanar_cases(box):
    a = box([0, 0, 0], [1, 1, 1])
    assert (a | box([1, 0, 0], [2, 1, 1])).volume == pytest.approx(2.0)
    assert (a | box([0.5, 0, 0], [1.5, 1, 1])).volume == pytest.approx(1.5)
    assert (a & box([0.5, 0, 0], [1.5, 1, 1])).volume == pytest.approx(0.5)
    assert (a - box([0.25, 0.25, 0.25], [0.75, 0.75, 0.75])).volume == pytest.approx(1 - 0.125)


def test_invalid_inputs_raise_with_the_report(box):
    a = box([0, 0, 0], [1, 1, 1])
    holed = a.filter_faces(np.arange(12) != 0)
    with pytest.raises(dr.InvalidMesh, match="a") as e:
        holed | a
    assert e.value.report.summary["boundary_edges"] == 3
    with pytest.raises(dr.InvalidMesh):
        dr.union(a, holed)
    welded = dr.Mesh(np.vstack([a.coords, a.coords]), np.vstack([a.triangles[:6], a.triangles[6:] + 8]))
    with pytest.raises(dr.InvalidMesh):
        welded | box([0.5, 0.5, 0.5], [2, 2, 2])
    mesh, (ra, rb) = dr.union(welded, box([0.5, 0.5, 0.5], [2, 2, 2]), repair=True)
    assert mesh.volume == pytest.approx(1 + 1.5**3 - 0.125)
    assert ra.summary["vertex_merged"] == 8 and len(rb) == 0
    with pytest.raises(TypeError):
        a | 3


def test_clip_by_surface_splits_volume(orebody):
    o = dr.datasets.ORIGIN
    g = np.linspace(-300, 300, 7)
    gx, gy = np.meshgrid(g, g)
    z = 200 + 0.05 * gx + 0.03 * gy
    surface = dr.dtm(o + np.column_stack([gx.ravel(), gy.ravel(), z.ravel()]))
    below, above = orebody.clip_by_surface(surface)
    assert below.is_solid and above.is_solid
    assert close(below.volume + above.volume, orebody.volume, rel=1e-8)
    assert below.centroid[2] < above.centroid[2]


def test_split_and_slice_by_plane(orebody, box):
    o = dr.datasets.ORIGIN + [13.0, 0, 200]
    neg, pos = orebody.split_by_plane(o, [1, 0.2, 0.1])
    assert neg.is_solid and pos.is_solid
    assert close(neg.volume + pos.volume, orebody.volume, rel=1e-8)
    b = box([0, 0, 0], [2, 3, 4])
    rings = b.slice_by_plane([0, 0, 0], [0, 0, 1], offsets=[1.0, 2.5])
    assert len(rings) == 2 and rings.closed.all()
    assert rings.area().tolist() == [6.0, 6.0]
    assert rings.attributes["offset"].tolist() == [1.0, 2.5]
    top = box([0, 0, 0], [1, 1, 1]).filter_faces(np.arange(12) < 2)
    lo, hi = top.split_by_plane([0.5, 0, 0], [1, 0, 0])
    assert lo.area + hi.area == pytest.approx(top.area)


def test_cut_fill_identities():
    t = dr.datasets.topography(n=120)
    zero = dr.cut_fill(t, t)
    assert zero["cut"] == 0 and zero["fill"] == 0
    h = 2.5
    up = dr.cut_fill(t, t.translate([0, 0, h]))
    assert close(up["fill"], h * up["area"], rel=1e-8) and up["cut"] == 0
    down = dr.cut_fill(t.translate([0, 0, h]), t)
    assert close(down["cut"], up["fill"], rel=1e-8) and down["fill"] == up["cut"]
    pit = dr.datasets.pit(t, n=21)
    result = dr.cut_fill(t, pit)
    assert all(np.isfinite(v) for v in result.values())
    assert close(result["net"], result["fill"] - result["cut"])
    o = dr.datasets.ORIGIN
    ring = dr.Polylines([o[:2] + [[-100, -100], [100, -100], [100, 100], [-100, 100]]], closed=True)
    inside = dr.cut_fill(t, t.translate([0, 0, 1.0]), boundary=ring)
    assert close(inside["area"], 40_000.0, rel=1e-9) and close(inside["fill"], 40_000.0, rel=1e-8)
