import numpy as np
import pytest
import wireframes as wf


def test_clean_solid_reports_nothing(box, orebody):
    for m in (box([0, 0, 0], [1, 1, 1]), orebody):
        r = m.validate()
        assert r.is_solid and m.is_solid and m.is_closed and len(r) == 0
        assert r.summary["shells"] == 1 and r.summary["is_solid"] is True
        assert r.problems.column_names == ["kind", "face", "vertex", "edge", "other"]


def test_removing_one_face_opens_three_edges(box):
    m = box([0, 0, 0], [1, 1, 1]).filter_faces(np.arange(12) != 0)
    r = m.validate()
    assert not r.is_solid and not m.is_closed
    assert r.summary["boundary_edges"] == 3
    problems = r.problems
    assert problems["kind"].tolist() == ["boundary_edge"] * 3
    assert problems.mask("edge").all() and problems.mask("face").all()


def test_duplicate_vertices_within_tolerance(box):
    b = box([0, 0, 0], [1, 1, 1])
    near = wf.Mesh(np.vstack([b.coords, b.coords[:1] + 1e-6]), b.triangles)
    assert near.validate().summary["duplicate_vertices"] == 0
    assert near.validate(tolerance=1e-5).summary["duplicate_vertices"] == 1
    assert near.validate().summary["unused_vertices"] == 1
    with pytest.raises(wf.InvalidInput):
        b.validate(tolerance=-1.0)


def test_overlapping_cubes_self_intersect(box):
    m = wf.concat([box([0, 0, 0], [2, 2, 2]), box([1, 1, 1], [3, 3, 3])])
    r = m.validate()
    assert r.summary["self_intersections"] > 0 and not r.is_solid
    assert r.summary["shells"] == 2


def test_inward_shell_and_flipped_face(box):
    b = box([0, 0, 0], [1, 1, 1])
    assert b.flip().validate().summary["inward_shells"] == 1
    t = b.triangles.copy()
    t[0] = t[0][::-1]
    assert wf.Mesh(b.coords, t).validate().summary["inconsistent_edges"] > 0


def test_counts_survive_face_permutation(rng, orebody):
    m = orebody.filter_faces(np.arange(orebody.n_faces) % 7 != 0)
    order = rng.permutation(m.n_faces)
    shuffled = wf.Mesh(m.coords, m.triangles[order])
    assert shuffled.validate().summary == m.validate().summary
