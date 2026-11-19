import numpy as np
import pytest
import wireframes as wf
from conftest import close


def pairs(box, rng, n=6):
    for _ in range(n):
        lo1 = rng.integers(-3, 3, 3).astype(float) + rng.uniform(0, 0.5, 3)
        lo2 = lo1 + rng.uniform(-1.5, 1.5, 3)
        yield box(lo1, lo1 + rng.uniform(1, 3, 3), name="a"), box(lo2, lo2 + rng.uniform(1, 3, 3), name="b")


def pieces_touch_only_along_curves(m):
    """``a ^ b`` holds ``a - b`` and ``b - a`` as separate shells that touch along the
    intersection curve: a valid solid whose report counts the contact."""
    s = m.validate().summary
    return s["is_solid"] and s["shells"] >= 2 and s["touching_shell_edges"] > 0


def vol(m):
    """Volume, 0 for an empty result (nothing kept)."""
    return 0.0 if m.n_faces == 0 else m.volume


def test_volume_identity_and_closed_output(box, rng):
    for a, b in pairs(box, rng):
        u, i, d, x = a | b, a & b, a - b, a ^ b
        for m in (u, i, d):
            assert m.n_faces == 0 or m.is_solid
        assert x.n_faces == 0 or x.is_solid
        assert close(vol(u) + vol(i), vol(a) + b.volume)
        assert close(vol(d), vol(a) - vol(i))
        assert close(vol(x), vol(u) - vol(i))
        assert close(vol(x), vol(a) + vol(b) - 2 * vol(i))
        assert close(vol(b | a), vol(u)) and close(vol(b & a), vol(i))


def test_methods_and_functions_agree(box):
    a, b = box([0, 0, 0], [2, 2, 2]), box([1, 1, 1], [3, 3, 3])
    for name, op in [("union", "__or__"), ("intersection", "__and__"), ("difference", "__sub__"),
                     ("symmetric_difference", "__xor__")]:  # fmt: skip
        v = getattr(wf, name)(a, b).volume
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
    o = wf.datasets.ORIGIN
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


def test_touching_shells_are_valid_and_chain(box):
    a = box([0, 0, 0], [1, 1, 1])
    edge, corner, face = box([1, 1, 0], [2, 2, 1]), box([1, 1, 1], [2, 2, 2]), box([1, 0, 0], [2, 1, 1])
    for b, shells, edges, points in [(edge, 2, 1, 2), (corner, 2, 0, 1), (face, 1, 0, 0)]:
        for m in (a | b, a ^ b):
            s = m.validate().summary
            assert s["is_solid"] and m.volume == pytest.approx(2.0)
            assert (s["shells"], s["touching_shell_edges"], s["touching_shell_vertices"]) == (
                shells,
                edges,
                points,
            )
        assert (a & b).n_faces == 0 and (a - b).volume == pytest.approx(1.0)
    two = a | edge
    joined = two | box([0.5, 0.5, 0], [1.5, 1.5, 1])
    assert joined.is_solid and joined.validate().summary["shells"] == 1
    assert joined.volume == pytest.approx(2.0 + 1.0 - 0.5)
    p = np.array([[0.5, 0.5, 0.5], [1.5, 1.5, 0.5], [1.5, 0.5, 0.5], [1.0, 1.0, 0.5]])
    assert two.locate(p).tolist() == [-1, -1, 1, 0]
    d = two.signed_distance(p)
    assert d[0] == pytest.approx(-0.5) and d[1] == pytest.approx(-0.5)
    assert d[2] == pytest.approx(0.5) and d[3] == pytest.approx(0.0)
    l_shape = box([0, 0, 0], [2, 2, 2]) - box([1, 1, 1], [3, 3, 3])
    pinched = l_shape - box([-1, -1, -1], [1, 1, 3])
    s = pinched.validate().summary
    assert pinched.is_solid and s["shells"] == 1 and s["touching_shell_edges"] == 1
    assert pinched.volume == pytest.approx(7.0 - 2.0)


def test_coplanar_cases(box):
    a = box([0, 0, 0], [1, 1, 1])
    assert (a | box([1, 0, 0], [2, 1, 1])).volume == pytest.approx(2.0)
    assert (a | box([0.5, 0, 0], [1.5, 1, 1])).volume == pytest.approx(1.5)
    assert (a & box([0.5, 0, 0], [1.5, 1, 1])).volume == pytest.approx(0.5)
    assert (a - box([0.25, 0.25, 0.25], [0.75, 0.75, 0.75])).volume == pytest.approx(1 - 0.125)


def test_invalid_inputs_raise_with_the_report(box):
    a = box([0, 0, 0], [1, 1, 1])
    holed = a.filter_faces(np.arange(12) != 0)
    with pytest.raises(wf.InvalidMesh, match="a") as e:
        holed | a
    assert e.value.report.summary["boundary_edges"] == 3
    with pytest.raises(wf.InvalidMesh):
        wf.union(a, holed)
    welded = wf.Mesh(np.vstack([a.coords, a.coords]), np.vstack([a.triangles[:6], a.triangles[6:] + 8]))
    with pytest.raises(wf.InvalidMesh):
        welded | box([0.5, 0.5, 0.5], [2, 2, 2])
    mesh, (ra, rb) = wf.union(welded, box([0.5, 0.5, 0.5], [2, 2, 2]), repair=True)
    assert mesh.volume == pytest.approx(1 + 1.5**3 - 0.125)
    assert ra.summary["vertex_merged"] == 8 and len(rb) == 0
    with pytest.raises(TypeError):
        a | 3


def test_clip_by_surface_splits_volume(orebody):
    o = wf.datasets.ORIGIN
    g = np.linspace(-300, 300, 7)
    gx, gy = np.meshgrid(g, g)
    z = 200 + 0.05 * gx + 0.03 * gy
    surface = wf.dtm(o + np.column_stack([gx.ravel(), gy.ravel(), z.ravel()]))
    below, above = orebody.clip_by_surface(surface)
    assert below.is_solid and above.is_solid
    assert close(below.volume + above.volume, orebody.volume, rel=1e-8)
    assert below.centroid[2] < above.centroid[2]


def test_split_and_slice_by_plane(orebody, box):
    o = wf.datasets.ORIGIN + [13.0, 0, 200]
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
    corner, empty = box([0, 0, 0], [1, 1, 1]).split_by_plane([1, 1, 1], [1, 1, 1])
    assert corner.n_faces == 12 and empty.n_faces == 0


def test_split_of_a_closed_dtm_at_its_local_minima_gives_solids():
    t = wf.datasets.topography()
    s = t.close_to_elevation(float(t.z.min()) - 50)
    v, f = t.coords, np.asarray(t.triangles)
    lowest = np.full(len(v), np.inf)
    for k in range(3):
        for j in (1, 2):
            np.minimum.at(lowest, f[:, k], v[f[:, (k + j) % 3], 2])
    minima = np.flatnonzero(lowest > v[:, 2])
    assert len(minima) >= 5
    for i in minima[:8]:
        lower, upper = s.split_by_plane(v[i], [0, 0, 1])
        assert lower.is_solid and upper.is_solid, i
        assert close(lower.volume + upper.volume, s.volume, rel=1e-9)


def test_cut_fill_identities():
    t = wf.datasets.topography(n=120)
    zero = wf.cut_fill(t, t)
    assert zero["cut"] == 0 and zero["fill"] == 0
    h = 2.5
    up = wf.cut_fill(t, t.translate([0, 0, h]))
    assert close(up["fill"], h * up["area"], rel=1e-8) and up["cut"] == 0
    down = wf.cut_fill(t.translate([0, 0, h]), t)
    assert close(down["cut"], up["fill"], rel=1e-8) and down["fill"] == up["cut"]
    pit = wf.datasets.pit(t, n=21)
    result = wf.cut_fill(t, pit)
    assert all(np.isfinite(v) for v in result.values())
    assert close(result["net"], result["fill"] - result["cut"])
    o = wf.datasets.ORIGIN
    ring = wf.Polylines([o[:2] + [[-100, -100], [100, -100], [100, 100], [-100, 100]]], closed=True)
    inside = wf.cut_fill(t, t.translate([0, 0, 1.0]), boundary=ring)
    assert close(inside["area"], 40_000.0, rel=1e-9) and close(inside["fill"], 40_000.0, rel=1e-8)


def test_rounding_report_comes_last_and_its_bounds_hold():
    topography = wf.datasets.topography()
    sky = topography.translate([0, 0, 200])
    closed1 = wf.datasets.pit(sky, radius=60).close_to_elevation(600)
    closed2 = wf.datasets.pit(sky, radius=90).close_to_elevation(600)
    pit1 = closed1.clip_by_surface(topography)[0]
    pit2 = closed2.clip_by_surface(topography)[0]
    (u, ru), (i, ri) = (
        wf.union(pit1, pit2, rounding_report=True),
        pit1.intersection(pit2, rounding_report=True),
    )
    for m, r in ((u, ru), (i, ri)):
        assert isinstance(r, wf.RoundingReport) and m.is_solid
        assert r.steps <= r.step_bound
        assert all(x <= b for x, b in zip(r.max_move, r.move_bound, strict=True))
    gap = abs(u.volume + i.volume - pit1.volume - pit2.volume)
    assert gap <= ru.volume_bound + ri.volume_bound + 1e-12 * pit2.volume
    m, repairs, r = wf.difference(pit2, pit1, repair=True, rounding_report=True)
    assert m.is_solid and len(repairs) == 2 and isinstance(r, wf.RoundingReport)
    assert r.inserted_points > 0 and r.moved_input_vertices == 0
    below, above, reports = pit2.clip_by_surface(topography, rounding_report=True)
    assert below.is_solid and above.n_faces == 0 and len(reports) == 2
    lo, hi, reports = pit2.split_by_plane(pit2.centroid, [1, 2, 3], rounding_report=True)
    assert lo.is_solid and hi.is_solid and all(isinstance(r, wf.RoundingReport) for r in reports)
    assert "volume_bound" in repr(r)


def test_rounding_report_has_a_runnable_example():
    import doctest

    finder = doctest.DocTestFinder(recurse=False)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    tests = [t for t in finder.find(wf.RoundingReport, "RoundingReport", globs={"wf": wf}) if t.examples]
    assert tests and all(runner.run(t, out=lambda s: None).failed == 0 for t in tests)
