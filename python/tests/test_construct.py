import itertools

import numpy as np
import pytest
import wireframes as wf
from conftest import close


def grid_points(n=6, size=10.0, f=lambda x, y: 0 * x):
    g = np.linspace(0, size, n)
    gx, gy = np.meshgrid(g, g)
    x, y = gx.ravel(), gy.ravel()
    return np.column_stack([x, y, f(x, y)])


def sampled_hausdorff(a, b):
    c = a.coords[a.triangles]
    samples = np.vstack([a.coords, c.mean(axis=1), ((c + np.roll(c, -1, axis=1)) / 2).reshape(-1, 3)])
    return b.distance(samples).max()


def test_dtm_keeps_points_and_winds_up(rng):
    p = rng.uniform(0, 100, (200, 3))
    s = wf.dtm(p)
    assert s.n_vertices == 200 and (s.face_normals[:, 2] > 0).all()
    assert np.array_equal(np.sort(s.coords, axis=0), np.sort(p, axis=0))
    near = np.vstack([p, p[:5] + [1e-3, 0, 0]])
    assert wf.dtm(near, min_spacing=0.01).n_vertices == 200
    pts = wf.PointSet(p, length_unit="m", crs="EPSG:31983")
    labelled = wf.dtm(pts)
    assert labelled.length_unit == "m" and labelled.crs == "EPSG:31983"
    with pytest.raises(wf.InvalidInput):
        wf.dtm([[0, 0, 0], [1, 1, 0], [2, 2, 0]])


def test_dtm_breaklines_and_boundary():
    p = grid_points(11)
    brk = wf.Polylines([[[0.5, 0.3, 4.0], [9.5, 9.1, 4.0]]])
    s = wf.dtm(p, breaklines=brk)
    edges = {tuple(sorted(e)) for t in s.triangles for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
    a, b = brk.coords
    d = b[:2] - a[:2]
    rel = s.coords[:, :2] - a[:2]
    t = rel @ d / (d @ d)
    on = np.flatnonzero(
        np.isclose(rel[:, 0] * d[1] - rel[:, 1] * d[0], 0, atol=1e-9) & (t >= -1e-12) & (t <= 1 + 1e-12)
    )
    chain = on[np.argsort(t[on])]
    assert np.array_equal(s.coords[chain[[0, -1]]], brk.coords)
    assert np.allclose(s.coords[chain, 2], 4.0)
    assert all(tuple(sorted((int(p), int(q)))) in edges for p, q in itertools.pairwise(chain))
    ring = wf.Polylines([[[2, 2], [8, 2], [8, 8], [2, 8]]], closed=True)
    clipped = wf.dtm(p, boundary=ring)
    assert clipped.area == pytest.approx(36.0)


def test_triangulate_area_matches_polygon_area():
    outer = [[0, 0], [4, 0], [4, 3], [0, 3]]
    hole = [[1, 1], [2, 1], [2, 2], [1, 2]]
    p = wf.Polylines([outer, hole, [[10, 0], [12, 0], [11, 2]]], closed=True, features=[0, 0, 1])
    m = wf.triangulate(p)
    assert m.area == pytest.approx(p.area().sum())
    assert sorted(set(m["feature"].tolist())) == [0, 1]


def test_extrude_volume_is_area_times_height():
    ring = wf.Polylines([[[0, 0], [3, 0], [3, 1], [1, 2], [0, 1]]], closed=True)
    m = wf.extrude(ring, [0, 0, 4])
    assert m.is_solid and close(m.volume, ring.area()[0] * 4)
    with pytest.raises(wf.InvalidInput):
        wf.extrude(ring, [1, 0, 0])


def test_loft_of_two_squares_is_a_box():
    sq = np.array([[0, 0, 0], [2, 0, 0], [2, 3, 0], [0, 3, 0]], float)
    m = wf.loft(wf.Polylines([sq, sq + [0, 0, 5]], closed=True))
    assert m.is_solid and close(m.volume, 30.0)
    assert set(m["section"].tolist()) == {0, -1, -2}
    open_ends = wf.loft(wf.Polylines([sq, sq + [0, 0, 5]], closed=True), cap_ends=False)
    assert not open_ends.is_closed


def test_close_to_elevation_matches_cut_fill(topography):
    floor = 200.0
    solid = topography.close_to_elevation(floor)
    assert solid.is_solid
    lo, hi = topography.bounds
    plane = wf.dtm(
        [[lo[0], lo[1], floor], [hi[0], lo[1], floor], [hi[0], hi[1], floor], [lo[0], hi[1], floor]]
    )
    cf = wf.cut_fill(plane, topography)
    assert close(solid.volume, cf["fill"], rel=1e-8)


def test_contours_of_a_plane_are_straight():
    s = wf.dtm(grid_points(6, f=lambda x, y: 0.3 * x + 0.1 * y + 2))
    lines = s.contours([3.0, 4.0])
    assert lines.attributes["elevation"].tolist() == [3.0, 4.0]
    for part in lines.parts:
        assert np.allclose(0.3 * part[:, 0] + 0.1 * part[:, 1] + 2, part[:, 2])
        d = part[1:] - part[:1]
        assert np.allclose(d[:, 0] * d[-1, 1] - d[:, 1] * d[-1, 0], 0, atol=1e-9)


def test_decimate_smooth_offset_remesh(orebody):
    flat = wf.dtm(grid_points(9))
    d = flat.decimate(target_faces=20)
    assert d.n_faces <= 20 and d.area == pytest.approx(flat.area)
    sd = orebody.decimate(target_faces=200)
    summary = sd.validate().summary
    assert sd.is_closed and sd.n_faces <= 200
    assert summary["non_manifold_edges"] == 0 and summary["non_manifold_vertices"] == 0
    smooth = orebody.smooth(iterations=3)
    assert smooth.n_faces == orebody.n_faces and smooth.is_closed
    bumpy = wf.dtm(grid_points(7, f=lambda x, y: np.sin(x) * np.cos(y)))
    sm = bumpy.smooth()
    boundary = (bumpy.x == 0) | (bumpy.x == 10) | (bumpy.y == 0) | (bumpy.y == 10)
    assert np.array_equal(sm.coords[boundary], bumpy.coords[boundary])
    grown = orebody.offset(1.0)
    assert grown.is_solid and grown.volume > orebody.volume
    assert np.allclose(np.linalg.norm(grown.coords - orebody.coords, axis=1), 1.0)
    cube = wf.datasets.box([0, 0, 0], [1, 1, 1]).offset(0.5)
    assert np.isclose(cube.volume, (1 + 1 / np.sqrt(3)) ** 3, rtol=1e-12)
    assert np.allclose(cube.coords.max(axis=0), 1 + 0.5 / np.sqrt(3), rtol=1e-12)
    box = wf.datasets.box([0, 0, 0], [4, 4, 4])
    re = box.remesh(1.0, iterations=2)
    assert re.is_closed and re.n_faces > 12
    assert max(sampled_hausdorff(box, re), sampled_hausdorff(re, box)) < 1.0
    with pytest.raises(wf.InvalidInput):
        orebody.smooth(lambda_=-1.0)
