import importlib.util

import numpy as np
import pytest
import wireframes as wf


def test_table_types_nulls_and_units():
    t = wf.Table(
        {"f": [1.0, np.nan, 3.0], "i": [1, None, 3], "b": [True, None, False], "s": ["a", None, "c"]}
    )
    assert t.num_rows == len(t) == 3
    assert t.column_names == ["f", "i", "b", "s"]
    assert t["f"].dtype == np.float64 and np.isnan(t["f"][1])
    assert t["i"].dtype == np.int64 and t["b"].dtype == np.bool_ and t["s"].dtype == object
    assert t["s"][1] is None
    for name in t.column_names:
        assert t.mask(name).tolist() == [True, False, True]
    assert t.drop_null().num_rows == 2
    assert t.drop_null("f").num_rows == 2
    u = t.with_column("grade", [0.5, 1.0, 1.5], unit="g/t")
    assert u.units == {"grade": "g/t"}
    assert "grade" in u and "grade" not in t
    assert t.filter([True, False, True])["i"].tolist() == [1, 3]
    assert set(t.with_columns({"k": [1, 2, 3]}).column_names) == {"f", "i", "b", "s", "k"}
    assert wf.Table().with_column("a", [1, 2]).num_rows == 2
    assert t.to_dict()["s"].tolist() == ["a", None, "c"]
    with pytest.raises(wf.MissingColumn, match="columns are"):
        t.column("nope")
    with pytest.raises(wf.InvalidInput):
        wf.Table({"a": [1, 2], "b": [1]})
    with pytest.raises(wf.InvalidInput):
        wf.Table({"a": [object(), 1]})


@pytest.mark.skipif(importlib.util.find_spec("pandas") is None, reason="needs pandas")
def test_table_to_pandas_keeps_nulls():
    df = wf.Table({"i": [1, None], "s": ["a", None]}).to_pandas()
    assert df["i"].isna().tolist() == [False, True]


NULLED = {"f": [1.0, np.nan, 3.0], "i": [1, None, 3], "b": [True, None, False], "s": ["a", None, "c"]}


def null_lists(table, frame, kind):
    if kind == "pandas":
        return {n: frame[n].isna().tolist() for n in table.column_names}
    if kind == "polars":
        return {n: frame[n].is_null().to_list() for n in table.column_names}
    return {n: frame.column(n).is_null().to_pylist() for n in table.column_names}


@pytest.mark.parametrize("kind", ["pandas", "polars", "pyarrow"])
def test_table_conversions_keep_values_and_nulls(kind):
    pytest.importorskip(kind)
    t = wf.Table(NULLED)
    frame = getattr(t, f"to_{kind}")()
    assert list(frame.columns if kind != "pyarrow" else frame.column_names) == t.column_names
    assert len(frame) == 3
    assert null_lists(t, frame, kind) == {n: [False, True, False] for n in t.column_names}
    row0 = [frame[n][0] if kind != "pyarrow" else frame[n][0].as_py() for n in t.column_names]
    assert row0 == [1.0, 1, True, "a"]
    full = wf.Table({"a": [1, 2]})
    assert null_lists(full, getattr(full, f"to_{kind}")(), kind) == {"a": [False, False]}


@pytest.mark.parametrize("kind", ["pandas", "polars", "pyarrow"])
def test_container_conversions_match_to_table(kind, box, rng):
    pytest.importorskip(kind)
    mesh = box([0, 0, 0], [1, 1, 1], name="b").with_vertex_column("id", np.arange(8))
    lines = wf.Polylines([[[0, 0, 0], [1, 0, 0]], [[2, 2, 2], [3, 3, 3]]], attributes={"n": ["x", None]})
    points = wf.PointSet(rng.uniform(0, 1, (5, 3)), {"g": [0.1, None, 0.3, 0.4, 0.5]})
    cases = [
        (mesh, {}, mesh.to_table()),
        (mesh, {"on": "faces"}, mesh.to_table(on="faces")),
        (lines, {}, lines.to_table()),
        (points, {}, points.to_table()),
    ]
    for obj, kw, table in cases:
        frame = getattr(obj, f"to_{kind}")(**kw)
        names = list(frame.column_names if kind == "pyarrow" else frame.columns)
        assert names == table.column_names and len(frame) == table.num_rows
        assert null_lists(table, frame, kind) == {n: (~table.mask(n)).tolist() for n in table.column_names}


def test_mesh_construction_and_attributes(box):
    b = box([0, 0, 0], [1, 2, 3], name="b")
    assert b.n_vertices == len(b) == 8 and b.n_faces == 12
    assert b.bounds == ([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])
    assert np.array_equal(b.x, b.coords[:, 0])
    assert b.column_names == ["name"]
    assert b["name"].tolist() == ["b"] * 12
    m = b.with_vertex_column("id", np.arange(8)).with_face_column("w", np.ones(12), unit="m")
    assert m.column_names == ["id", "name", "w"]
    assert m.face_attributes.units == {"w": "m"}
    both = m.with_face_column("id", np.arange(12))
    with pytest.raises(wf.InvalidInput, match="both"):
        both["id"]
    assert m.to_table().column_names == ["x", "y", "z", "id"]
    assert m.to_table(on="faces").column_names == ["v0", "v1", "v2", "name", "w"]
    with pytest.raises(wf.InvalidInput):
        m.to_table(on="edges")
    moved = m.with_vertices(m.coords + 1)
    assert np.array_equal(moved.triangles, m.triangles) and moved.column_names == m.column_names
    half = m.filter_faces(np.arange(12) < 6)
    assert half.n_faces == 6 and half.n_vertices <= 8
    t = m.with_attributes(face={"k": np.arange(12)})
    assert t.column_names == ["id", "k"]
    with pytest.raises(wf.InvalidInput):
        wf.Mesh([[0, 0, 0]], [[0, 1, 2]])
    with pytest.raises(wf.InvalidInput):
        wf.Mesh([[0, 0]], [])
    with pytest.raises(wf.InvalidInput):
        b.with_face_column("x", [1, 2])


def test_mesh_is_immutable(box):
    b = box([0, 0, 0], [1, 1, 1])
    before = b.coords.copy()
    with pytest.raises(ValueError, match="read-only"):
        b.coords[0, 0] = 99.0
    b.translate([5, 5, 5])
    assert np.array_equal(b.coords, before)
    with pytest.raises(AttributeError):
        b.crs = "EPSG:4326"


def test_concat_offsets_indices(box):
    a, b = box([0, 0, 0], [1, 1, 1], name="a"), box([2, 0, 0], [3, 1, 1])
    c = wf.concat([a, b])
    assert c.n_vertices == 16 and c.n_faces == 24
    assert c.triangles[12:].min() == 8
    assert c["name"].tolist() == ["a"] * 12 + [None] * 12
    assert c.volume == pytest.approx(2.0)


def test_polylines(rng):
    square = [[0, 0], [2, 0], [2, 2], [0, 2]]
    hole = [[0.5, 0.5], [1.5, 0.5], [1.5, 1.5], [0.5, 1.5]]
    line = [[5, 5, 1], [6, 5, 2]]
    p = wf.Polylines(
        [square, hole, line],
        closed=[True, True, False],
        features=[0, 0, 1],
        attributes={"name": ["pad", "trace"]},
        vertex_attributes={"k": np.arange(10)},
        length_unit="m",
    )
    assert len(p) == 2 and p.feature.tolist() == [0, 0, 1]
    assert p.closed.tolist() == [True, True, False]
    assert p.area().tolist() == [3.0, 0.0]
    assert p.length()[0] == pytest.approx(12.0)
    assert p.contains([[0.25, 0.25], [1, 1], [9, 9]]).tolist() == [True, False, False]
    assert p.contains([[1, 1]], feature=0).tolist() == [False]
    assert p.locate([[0.25, 0.25], [9, 9]]).tolist() == [0, -1]
    d = p.distance([[1, 1], [-1, 1]], feature=0, signed=True)
    assert d.tolist() == [0.5, 1.0]
    assert p.column_names == ["name", "k"]
    table = p.to_table()
    assert table.column_names == ["feature", "part", "vertex", "x", "y", "z", "closed", "k", "name"]
    assert table["name"].tolist() == ["pad"] * 8 + ["trace"] * 2
    points = p.to_points()
    assert len(points) == 10 and "part" in points.column_names
    assert p.filter([False, True]).coords.tolist() == [[5, 5, 1], [6, 5, 2]]
    assert p.with_column("w", [1.0, 2.0]).attributes["w"].tolist() == [1.0, 2.0]
    assert wf.Polylines([[[0, 0]]]).coords.tolist() == [[0, 0, 0]]
    with pytest.raises(wf.InvalidInput):
        wf.Polylines([[[0, 0]]], closed=True)
    with pytest.raises(wf.InvalidInput):
        p.contains([[0, 0]], feature=7)


def test_polylines_from_table():
    t = {"ID": ["a", "a", "a", "b", "b"], "X": [0, 1, 1, 5, 6], "Y": [0, 0, 1, 0, 0], "Z": [1, 1, 1, 2, 2]}
    p = wf.Polylines.from_table(t, z="Z", crs="EPSG:31983")
    assert len(p) == 2 and p.attributes["ID"].tolist() == ["a", "b"]
    assert p.parts[1].tolist() == [[5, 0, 2], [6, 0, 2]] and p.crs == "EPSG:31983"
    parts = wf.Polylines.from_table({**t, "P": [0, 0, 0, 1, 1]}, z="Z", part="P")
    assert len(parts) == 2


def test_point_set(rng):
    xyz = rng.uniform(0, 10, (20, 3))
    p = wf.PointSet(xyz, {"grade": rng.uniform(0, 1, 20)}, crs="EPSG:31983")
    assert len(p) == 20 and p.column_names == ["grade"] and p.crs == "EPSG:31983"
    assert np.array_equal(p.coords, xyz)
    assert len(p.filter(p["grade"] > 0.5)) == int((p["grade"] > 0.5).sum())
    assert p.to_table().column_names == ["x", "y", "z", "grade"]
    assert p.with_column("k", np.arange(20))["k"].dtype == np.int64


def test_arrays_copied_from_containers_are_read_only(box):
    m = box([0, 0, 0], [1, 1, 1]).with_face_column("d", np.arange(12.0))
    holes = wf.datasets.drillholes(n=2)
    pts = holes.to_points()
    arrays = [m.coords, m.x, m.triangles, m["d"], m.face_attributes["d"], pts.coords]
    arrays += [holes.coords, holes.parts[0], holes.closed, holes.feature]
    for a in arrays:
        assert not a.flags.writeable
        with pytest.raises(ValueError):
            a[0] = 0
    c = m.coords.copy()
    c[:, 2] += 10
    assert m.with_vertices(c).bounds[0][2] == 10
