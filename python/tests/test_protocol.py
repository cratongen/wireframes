"""Protocol version 1, built by hand from the contract and read back.

Cross-library cases are separate and skip when that sibling is not installed.
"""

import numpy as np
import pytest
import wireframes as wf
from conftest import close


def spec(kind, **extra):
    out = {
        "protocol": 1,
        "kind": kind,
        "columns": {},
        "units": {},
        "categories": {},
        "crs": None,
        "length_unit": None,
    }
    out.update(extra)
    return out


def grid(**extra):
    out = spec(
        "grid",
        origin=np.array([0.0, 0.0, 0.0]),
        size=np.array([2.0, 2.0, 2.0]),
        count=(1, 1, 1),
        rotation=np.array([0.0, 0.0, 0.0]),
    )
    out.update(extra)
    return out


def block_model(layout, **extra):
    out = grid()
    out["kind"] = "block_model"
    out["layout"] = layout
    out["index"] = None
    out["parent"] = None
    out["offset"] = None
    out["scale"] = None
    out["subgrid"] = None
    out.update(extra)
    return out


class Carry:
    """A foreign object. ``block_fraction`` must read ``__minetypes__`` from it."""

    def __init__(self, payload):
        self.payload = payload

    def __minetypes__(self):
        return self.payload


def test_mesh_views_and_round_trip():
    src = wf.datasets.box([0, 0, 0], [1, 1, 1])
    mesh = wf.Mesh(src.coords, src.triangles, crs="local", length_unit="m")
    mesh = mesh.with_face_column("zone", ["ore"] * mesh.n_faces)
    d = mesh.__minetypes__()
    assert d["protocol"] == 1 and d["kind"] == "mesh"
    assert d["columns"] == {} and d["units"] == {} and d["categories"] == {}
    assert d["crs"] == "local" and d["length_unit"] == "m"
    assert d["coords"].dtype == np.float64 and d["coords"].shape == (mesh.n_vertices, 3)
    assert d["triangles"].dtype == np.uint32 and d["triangles"].shape == (mesh.n_faces, 3)
    assert d["coords"].flags.c_contiguous and not d["coords"].flags.writeable
    assert d["triangles"].flags.c_contiguous and not d["triangles"].flags.writeable
    assert np.shares_memory(d["coords"], mesh.__minetypes__()["coords"])
    assert np.shares_memory(d["triangles"], mesh.__minetypes__()["triangles"])
    assert d["coords"].base is mesh and d["triangles"].base is mesh
    assert not np.shares_memory(d["triangles"], mesh.triangles)
    assert d["face_columns"]["zone"].dtype.kind in "UT"
    assert list(d["face_columns"]["zone"]) == ["ore"] * mesh.n_faces
    assert "valid" not in d and "mask" not in d
    with pytest.raises(ValueError):
        d["coords"][0, 0] = 9.0
    back = wf.Mesh.from_minetypes(d)
    assert back.volume == 1.0 and back.crs == "local" and back.length_unit == "m"
    assert list(back.face_attributes["zone"]) == ["ore"] * mesh.n_faces
    again = wf.Mesh.from_minetypes(mesh)
    assert again.volume == 1.0


def test_mesh_consumer_names_the_broken_key():
    mesh = wf.datasets.box([0, 0, 0], [1, 1, 1])
    d = mesh.__minetypes__()
    broken = dict(d)
    broken["protocol"] = True
    with pytest.raises(wf.InvalidInput, match="protocol"):
        wf.Mesh.from_minetypes(broken)
    broken = dict(d)
    del broken["triangles"]
    with pytest.raises(wf.InvalidInput, match="triangles"):
        wf.Mesh.from_minetypes(broken)
    broken = dict(d)
    broken["triangles"] = np.array(d["triangles"], dtype=np.int64)
    with pytest.raises(wf.InvalidInput, match="triangles"):
        wf.Mesh.from_minetypes(broken)
    foreign = ">f8" if np.little_endian else "<f8"
    broken = dict(d)
    broken["coords"] = np.array(d["coords"], dtype=foreign)
    with pytest.raises(wf.InvalidInput, match="coords"):
        wf.Mesh.from_minetypes(broken)
    wide = np.zeros((d["coords"].shape[0], 6))
    broken = dict(d)
    broken["coords"] = wide[:, ::2]
    assert not broken["coords"].flags.c_contiguous
    with pytest.raises(wf.InvalidInput, match="coords"):
        wf.Mesh.from_minetypes(broken)
    broken = dict(d)
    broken["kind"] = "point_set"
    with pytest.raises(wf.InvalidInput, match="point_set"):
        wf.Mesh.from_minetypes(broken)
    broken = dict(d)
    broken["columns"] = {"x": np.zeros(mesh.n_vertices)}
    with pytest.raises(wf.InvalidInput, match="columns"):
        wf.Mesh.from_minetypes(broken)
    snapshot = np.array(d["coords"])
    wf.Mesh.from_minetypes(d)
    assert np.array_equal(d["coords"], snapshot)


def test_point_set_and_table_columns():
    points = wf.PointSet([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]], {"au": [0.5, 1.5]}, length_unit="m")
    d = points.__minetypes__()
    assert d["kind"] == "point_set" and d["coords"].shape == (2, 3)
    assert np.shares_memory(d["coords"], points.__minetypes__()["coords"])
    assert d["coords"].base is points
    assert d["columns"]["au"].dtype == np.float64
    assert np.shares_memory(d["columns"]["au"], points.__minetypes__()["columns"]["au"])
    back = wf.PointSet.from_minetypes(Carry(d))
    assert back.length_unit == "m" and back["au"].tolist() == [0.5, 1.5]

    plain = wf.Table({"b": [1.0], "a": [2.0]})
    plain = plain.with_column("au", [4.0], unit="g/t")
    produced = plain.__minetypes__()
    assert produced["kind"] == "table" and produced["rows"] == 1
    assert list(produced["columns"]) == ["b", "a", "au"]
    assert produced["units"] == {"au": "g/t"} and produced["categories"] == {}
    assert produced["crs"] is None and produced["length_unit"] is None
    assert np.shares_memory(produced["columns"]["b"], plain.__minetypes__()["columns"]["b"])
    assert produced["columns"]["b"].base is plain

    masked = wf.Table({"au": [1.0, None], "n": [3, None], "flag": [True, None], "rock": ["fr", None]})
    raw = masked.__minetypes__()
    assert "valid" not in raw and "mask" not in raw
    assert np.isnan(raw["columns"]["au"][1]) and raw["columns"]["au"][0] == 1.0
    assert not np.shares_memory(raw["columns"]["au"], masked.__minetypes__()["columns"]["au"])
    assert raw["columns"]["n"].dtype == np.float64 and np.isnan(raw["columns"]["n"][1])
    assert raw["columns"]["n"][0] == 3.0 and "n" not in raw["categories"]
    assert raw["columns"]["flag"].dtype == np.int32 and raw["columns"]["flag"].tolist() == [1, -1]
    assert raw["categories"]["flag"]["names"] == ["false", "true"]
    assert raw["categories"]["flag"]["colors"] is None and raw["categories"]["flag"]["mapping"] == {}
    assert raw["columns"]["rock"].tolist() == [0, -1]
    assert raw["categories"]["rock"]["names"] == ["fr"]
    restored = wf.Table.from_minetypes(raw)
    assert restored.mask("au").tolist() == [True, False]
    assert restored["n"].dtype == np.float64 and restored.mask("n").tolist() == [True, False]
    assert restored["flag"].dtype == np.bool_ and restored.mask("flag").tolist() == [True, False]
    assert bool(restored["flag"][0])
    assert restored.mask("rock").tolist() == [True, False] and restored["rock"][0] == "fr"

    whole = wf.Table({"n": [-1, 4], "rock": ["ox", "fr"], "flag": [False, True]})
    whole_d = whole.__minetypes__()
    assert whole_d["columns"]["n"].dtype == np.int64 and whole_d["columns"]["n"].tolist() == [-1, 4]
    assert whole_d["columns"]["n"].base is whole
    assert whole_d["columns"]["rock"].dtype.kind in "UT"
    assert whole_d["columns"]["flag"].dtype == np.bool_
    assert np.shares_memory(whole_d["columns"]["flag"], whole.__minetypes__()["columns"]["flag"])
    text = wf.Table.from_minetypes(
        spec(
            "table",
            rows=2,
            columns={"rock": np.array(["ox", "fr"])},
        )
    )
    assert text["rock"].tolist() == ["ox", "fr"]


def test_table_refuses_blank_text_infinity_and_labels():
    assert wf.Table({"rock": ["ok", ""]}).__minetypes__()["columns"]["rock"].tolist() == ["ok", ""]
    with pytest.raises(wf.InvalidInput, match="empty"):
        wf.Table({"rock": ["", None]}).__minetypes__()
    blank = wf.Table({"rock": [None, None]}).__minetypes__()
    assert blank["columns"]["rock"].dtype == np.float64 and np.isnan(blank["columns"]["rock"]).all()
    assert blank["categories"] == {}
    with pytest.raises(wf.InvalidInput, match="finite|infinity"):
        wf.Table({"au": [1.0, np.inf]}).__minetypes__()
    kept = wf.Table({"n": [-1]})
    assert kept.__minetypes__()["columns"]["n"].dtype == np.int64
    empty = wf.Table.from_minetypes(spec("table", rows=3))
    assert empty.num_rows == 3 and empty.column_names == []
    labelled = spec("table", rows=0, crs="EPSG:31982")
    with pytest.raises(wf.InvalidInput, match="crs"):
        wf.Table.from_minetypes(labelled)
    wide = spec(
        "table",
        rows=1,
        columns={"n": np.array([3], dtype=np.int32), "z": np.array([1.5], dtype=np.float32)},
    )
    back = wf.Table.from_minetypes(wide)
    assert back["n"].dtype == np.int64 and back["z"].dtype == np.float64
    named = spec(
        "table",
        rows=1,
        columns={"rock": np.array([0], dtype=np.int32)},
        categories={
            "rock": {"names": ["sulphide"], "colors": ["#abc"], "mapping": {}, "other": None},
        },
    )
    assert wf.Table.from_minetypes(named)["rock"].tolist() == ["sulphide"]


def test_polylines_need_two_vertices_when_open():
    line = wf.Polylines([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]], attributes={"name": ["drive"]})
    line = line.with_vertex_column("k", [1, 2])
    d = line.__minetypes__()
    assert d["kind"] == "polylines"
    assert d["part_offsets"].dtype == np.uint32 and d["part_offsets"].tolist() == [0, 2]
    assert d["feature_offsets"].dtype == np.uint32 and d["feature_offsets"].tolist() == [0, 1]
    assert d["closed"].dtype == np.bool_ and d["closed"].tolist() == [False]
    assert np.shares_memory(d["coords"], line.__minetypes__()["coords"])
    assert np.shares_memory(d["part_offsets"], line.__minetypes__()["part_offsets"])
    assert np.shares_memory(d["closed"], line.__minetypes__()["closed"])
    assert d["coords"].base is line
    assert d["columns"]["name"].tolist() == ["drive"]
    assert d["vertex_columns"]["k"].tolist() == [1, 2]
    back = wf.Polylines.from_minetypes(d)
    assert back.parts[0].shape == (2, 3) and back["name"].tolist() == ["drive"]
    assert back.vertex_attributes["k"].tolist() == [1, 2]

    ring = wf.Polylines([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]], closed=True)
    assert ring.__minetypes__()["part_offsets"].tolist() == [0, 3]
    assert ring.__minetypes__()["closed"].tolist() == [True]

    one = wf.Polylines([[[0.0, 0.0, 0.0]]])
    with pytest.raises(wf.InvalidInput, match="PointSet"):
        one.__minetypes__()
    hand = spec(
        "polylines",
        coords=np.array([[0.0, 0.0, 0.0]]),
        part_offsets=np.array([0, 1], dtype=np.uint32),
        feature_offsets=np.array([0, 1], dtype=np.uint32),
        closed=np.array([False]),
        vertex_columns={},
        vertex_units={},
        vertex_categories={},
    )
    with pytest.raises(wf.InvalidInput, match="PointSet") as info:
        wf.Polylines.from_minetypes(hand)
    assert "part 0" in str(info.value) and "vertex" in str(info.value)


def test_own_types_stay_on_the_array_path():
    solid = wf.datasets.box([0, 0, 0], [2, 2, 2])
    fractions = solid.block_fraction(solid, [1.0, 1.0, 1.0])
    assert len(fractions) == solid.n_vertices
    one = solid.block_fraction(wf.PointSet([[1.0, 1.0, 1.0]]), [2.0, 2.0, 2.0])
    assert one.tolist() == [1.0]
    runs = solid.trace_intervals([[1.0, 1.0, -1.0], [1.0, 1.0, 3.0]], [0.0, 4.0])
    assert "hole_index" not in runs.column_names


def test_grid_and_block_model_fractions():
    solid = wf.datasets.box([0, 0, 0], [2, 2, 2])
    cells = grid(origin=np.array([-1.0, -1.0, -1.0]), count=(2, 2, 2))
    fractions = solid.block_fraction(cells)
    assert fractions.shape == (8,)
    assert close(float(fractions.sum()) * 8.0, solid.volume, rel=1e-8)
    direct = solid.block_fraction([[1.0, 1.0, 1.0]], [2.0, 2.0, 2.0])
    single = solid.block_fraction(grid())
    assert single.tobytes() == direct.tobytes()

    masked = block_model(
        "masked", index=np.arange(8, dtype=np.uint64), origin=cells["origin"], count=(2, 2, 2)
    )
    assert solid.block_fraction(masked).tobytes() == fractions.tobytes()
    with pytest.raises(wf.InvalidInput, match="index"):
        solid.block_fraction(block_model("masked", index=np.array([0, 0], dtype=np.uint64)))
    stepped = np.arange(4, dtype=np.uint64)[::2]
    with pytest.raises(wf.InvalidInput, match="index"):
        solid.block_fraction(block_model("masked", index=stepped))

    offset = np.array([[0.25, 0.25, 0.25]], dtype=np.float32)
    scale = np.array([[0.5, 0.5, 0.5]], dtype=np.float32)
    sub = block_model(
        "subblocked", parent=np.array([0], dtype=np.uint64), offset=offset, scale=scale, subgrid=(2, 2, 2)
    )
    cube = wf.datasets.box([0, 0, 0], [1, 1, 1])
    got = cube.block_fraction(sub)
    at = cube.block_fraction([[0.5, 0.5, 0.5]], [1.0, 1.0, 1.0])
    assert got.tobytes() == at.tobytes() == np.array([1.0]).tobytes()

    centres = [[ox, oy, oz] for oz in (0.25, 0.75) for oy in (0.25, 0.75) for ox in (0.25, 0.75)]
    octants = block_model(
        "subblocked",
        parent=np.zeros(8, dtype=np.uint64),
        offset=np.array(centres, dtype=np.float32),
        scale=np.full((8, 3), 0.5, dtype=np.float32),
        subgrid=(2, 2, 2),
    )
    cover = solid.block_fraction(octants)
    assert close(float(cover.sum()) * 1.0, solid.volume, rel=1e-8)

    with pytest.raises(wf.InvalidInput, match="offset"):
        solid.block_fraction(
            block_model(
                "subblocked",
                parent=np.array([0], dtype=np.uint64),
                offset=np.array([[0.1, 0.5, 0.5]], dtype=np.float32),
                scale=np.array([[0.5, 0.5, 0.5]], dtype=np.float32),
                subgrid=None,
            )
        )
    with pytest.raises(wf.InvalidInput, match="count"):
        solid.block_fraction(grid(count=(2**31, 2**31, 2)))
    with pytest.raises(wf.InvalidInput, match="block_model|size"):
        solid.block_fraction(Carry(grid()), [2.0, 2.0, 2.0])
    with pytest.raises(wf.InvalidInput, match="grid or 'block_model'|kind"):
        solid.block_fraction(
            spec("mesh", coords=np.zeros((1, 3)), triangles=np.zeros((0, 3), dtype=np.uint32))
        )


def test_rotation_matches_the_array_call():
    # Azimuth 90: F.T columns are the block axes. Local centre (1, 2, 3) lands at (12, 19, 33).
    rotation = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    solid = wf.datasets.box([12.0, 18.0, 30.0], [14.0, 20.0, 36.0])
    model = grid(
        origin=np.array([10.0, 20.0, 30.0]),
        size=np.array([2.0, 4.0, 6.0]),
        rotation=np.array([90.0, 0.0, 0.0]),
    )
    via = solid.block_fraction(model)
    direct = solid.block_fraction([[12.0, 19.0, 33.0]], [2.0, 4.0, 6.0], rotation=rotation)
    assert via.tobytes() == direct.tobytes()
    assert close(float(via[0]) * 48.0, solid.volume, rel=1e-8)


def _holes():
    return spec(
        "drillholes",
        holes=["DH1", "DH2"],
        coords=np.array([[1.0, 1.0, -1.0], [9.0, 9.0, -1.0]]),
        survey=None,
        intervals=None,
        tables={},
        path_coords=np.array([[1.0, 1.0, -1.0], [1.0, 1.0, 3.0], [9.0, 9.0, -1.0], [9.0, 9.0, 3.0]]),
        path_depths=np.array([0.0, 4.0, 0.0, 4.0]),
        path_offsets=np.array([0, 2, 4], dtype=np.uint32),
        method="minimum_curvature",
    )


def test_drillholes_add_hole_index_and_keep_masks():
    solid = wf.datasets.box([0.0, 0.0, 0.0], [2.0, 2.0, 2.0])
    runs = solid.trace_intervals(Carry(_holes()))
    assert runs["hole_index"].dtype == np.int64
    assert runs["hole_index"].tolist() == sorted(runs["hole_index"].tolist())
    assert set(runs["hole_index"].tolist()) == {0, 1}
    first = runs.filter(runs["hole_index"] == 0)
    assert first["from"].tolist() == [0.0, 1.0, 3.0]
    assert first["to"].tolist() == [1.0, 3.0, 4.0]
    assert first["inside"].tolist() == [False, True, False]
    assert not runs.mask("entry_face")[0]
    assert not runs.mask("exit_face")[2]
    assert "hole_index" in runs.column_names
    outside = runs.filter(runs["hole_index"] == 1)
    assert not outside["inside"].any()

    bad = _holes()
    bad["survey"] = {
        "hole_index": np.array([0], dtype=np.int32),
        "depth": np.array([0.0]),
        "azimuth": np.array([0.0]),
        "dip": np.array([91.0]),
    }
    with pytest.raises(wf.InvalidInput, match="survey|dip"):
        solid.trace_intervals(bad)
    short = _holes()
    short["path_offsets"] = np.array([0, 1, 4], dtype=np.uint32)
    with pytest.raises(wf.InvalidInput, match="hole 0|vertices"):
        solid.trace_intervals(short)
    with pytest.raises(wf.InvalidInput, match="depths"):
        solid.trace_intervals(Carry(_holes()), [0.0, 4.0])


def test_boitata_block_model_sums_to_the_solid_volume():
    bt = pytest.importorskip("boitata")
    model = bt.BlockModel(origin=(-1.0, -1.0, -1.0), size=(2.0, 2.0, 2.0), count=(2, 2, 2))
    if not callable(getattr(model, "__minetypes__", None)):
        pytest.skip("installed boitata BlockModel has no __minetypes__")
    solid = wf.datasets.box([0.0, 0.0, 0.0], [2.0, 2.0, 2.0])
    fractions = solid.block_fraction(model)
    assert close(float(np.asarray(fractions).sum()) * 8.0, solid.volume, rel=1e-8)
