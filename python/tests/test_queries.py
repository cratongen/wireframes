import drusa as dr
import numpy as np
import pytest
from conftest import close


def test_locate_contains_and_signed_distance(orebody, rng):
    lo, hi = (np.array(b) for b in orebody.bounds)
    p = rng.uniform(lo - 10, hi + 10, (400, 3))
    loc = orebody.locate(p)
    sd = orebody.signed_distance(p)
    assert set(np.unique(loc)) <= {-1, 0, 1} and (loc == -1).any() and (loc == 1).any()
    assert np.array_equal(np.sign(sd).astype(np.int8), loc)
    assert np.array_equal(orebody.contains(p), loc <= 0)
    assert np.allclose(np.abs(sd), orebody.distance(p))
    assert (orebody.locate(orebody.coords) == 0).all()
    assert orebody.contains(orebody.coords).all()
    assert not orebody.contains(orebody.coords, boundary=False).any()
    assert np.array_equal(orebody.locate(dr.PointSet(p)), loc)
    with pytest.raises(dr.InvalidInput):
        orebody.locate([[0, 0, np.nan]])
    with pytest.raises(dr.InvalidMesh):
        dr.datasets.topography(n=20).locate(p)


def test_closest_point_never_beats_a_vertex(orebody, rng):
    p = rng.uniform(-50, 50, (100, 3)) + dr.datasets.ORIGIN + [0, 0, 200]
    q, face, d = orebody.closest_point(p)
    assert q.shape == (100, 3) and face.dtype == np.int64
    assert np.allclose(np.linalg.norm(q - p, axis=1), d)
    nearest_vertex = np.min(np.linalg.norm(p[:, None, :] - orebody.coords[None], axis=2), axis=1)
    assert (d <= nearest_vertex + 1e-9).all()


def test_ray_hits_parity_matches_locate(box, rng):
    b = box([0, 0, 0], [1, 1, 1])
    origins = rng.uniform(-2, -1, (50, 3))
    ends = rng.uniform(-0.5, 1.5, (50, 3))
    hits = b.ray_hits(origins, ends - origins)
    assert hits.column_names == ["ray", "t", "face", "x", "y", "z", "kind", "entering"]
    counts = np.bincount(hits["ray"][hits["t"] <= 1], minlength=50)
    inside = b.locate(ends) == -1
    assert np.array_equal(counts % 2 == 1, inside)
    seg = b.segment_hits(origins, ends)
    assert np.array_equal(np.bincount(seg["segment"], minlength=50), counts)
    open_box = b.filter_faces(np.arange(12) > 0)
    assert not open_box.ray_hits([[0.5, 0.5, -1]], [[0, 0, 1]]).mask("entering").any()


def test_trace_intervals_give_the_chord(box):
    b = box([0, 0, 0], [2, 2, 2])
    trace = np.array([[0.3, 0.7, -1.0], [0.3, 0.7, 1.0], [0.3, 0.7, 5.0]])
    iv = b.trace_intervals(trace, [10.0, 12.0, 16.0])
    assert iv["from"].tolist() == [10.0, 11.0, 13.0]
    assert iv["to"].tolist() == [11.0, 13.0, 16.0]
    assert iv["inside"].tolist() == [False, True, False]
    assert iv.mask("entry_face").tolist() == [False, True, True]


def test_drillholes_through_the_orebody(orebody):
    holes = dr.datasets.drillholes()
    total = 0.0
    for part in holes.parts:
        iv = orebody.trace_intervals(part, [0.0, 350.0])
        inside = iv.filter(iv["inside"])
        hits = orebody.segment_hits(part[:1], part[1:])
        total += (inside["to"] - inside["from"]).sum()
        assert len(hits) == 2 * len(inside)
    assert total > 0


def test_block_fraction_sums_to_solid_volume(orebody):
    lo, hi = (np.array(b) for b in orebody.bounds)
    size = np.array([25.0, 20.0, 10.0])
    axes = [np.arange(lo[k] - size[k] / 3, hi[k] + size[k], size[k]) for k in range(3)]
    centres = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    f = orebody.block_fraction(centres, size)
    assert ((f >= 0) & (f <= 1)).all() and (f == 1).any() and (f == 0).any()
    assert close(f.sum() * size.prod(), orebody.volume, rel=1e-8)
    theta = np.radians(30.0)
    rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    c = orebody.centroid
    span = np.max(hi - lo)
    g = np.arange(-span / 2 - 25, span / 2 + 25, 25.0)
    local = np.stack(np.meshgrid(g, g, np.arange(-40, 41, 10.0), indexing="ij"), -1).reshape(-1, 3)
    rotated = c + local @ rot.T
    fr = orebody.block_fraction(rotated, [25, 25, 10], rotation=rot)
    assert close(fr.sum() * 25 * 25 * 10, orebody.volume, rel=1e-8)
    with pytest.raises(dr.InvalidInput):
        orebody.block_fraction(centres[:2], [1, 1, -1])
    with pytest.raises(dr.InvalidInput):
        orebody.block_fraction(centres[:2], [1, 1, 1], rotation=np.diag([1.0, 1.0, -1.0]))


def test_block_fraction_is_one_inside_zero_outside(box):
    b = box([0, 0, 0], [4, 4, 4])
    f = b.block_fraction([[2, 2, 2], [9, 9, 9], [4, 2, 2]], [1, 1, 1])
    assert f.tolist() == [1.0, 0.0, 0.5]
