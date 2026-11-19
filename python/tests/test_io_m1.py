"""OBJ, PLY, OFF and STL through the Python API."""

import numpy as np
import pytest
import wireframes as wf


@pytest.fixture
def solid(box):
    b = box([0.5, 0.25, 0.125], [2.5, 1.25, 4.0])
    return b.with_face_column("rock", np.arange(12) % 3).with_vertex_column("t", np.linspace(0, 1, 8))


@pytest.mark.parametrize("ext", ["obj", "off", "ply", "stl"])
@pytest.mark.parametrize("ascii", [False, True])
def test_roundtrip_keeps_geometry(tmp_path, solid, ext, ascii):
    path = tmp_path / f"m.{ext}"
    plain = wf.Mesh(solid.coords, solid.triangles)
    wf.write(path, plain, ascii=ascii)
    back = wf.read_mesh(path)
    assert np.array_equal(back.coords[back.triangles], plain.coords[plain.triangles])
    if ext != "stl":
        assert np.array_equal(back.coords, plain.coords)
    assert back.volume == plain.volume


def test_ply_roundtrip_keeps_numeric_attributes(tmp_path, solid):
    for ascii in (False, True):
        path = tmp_path / f"a{ascii}.ply"
        wf.write(path, solid, ascii=ascii)
        back = wf.read_mesh(path)
        assert back["rock"].tolist() == solid["rock"].tolist()
        assert back["t"].tolist() == solid["t"].tolist()


def test_obj_roundtrip_polygons_and_lines(tmp_path):
    path = tmp_path / "q.obj"
    path.write_text("v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 5 5 5\nv 6 5 5\ng pad\nf 1 2 3 4\nl 5 6\n")
    data = wf.read(path)
    assert data.mesh.n_faces == 2 and data.mesh["group"].tolist() == ["pad", "pad"]
    assert len(data.polylines) == 1
    wf.write(tmp_path / "again.obj", data)
    again = wf.read(tmp_path / "again.obj")
    assert np.array_equal(again.mesh.triangles, data.mesh.triangles)


def test_read_labels_and_errors(tmp_path, solid):
    path = tmp_path / "m.stl"
    wf.write(path, solid, skip_unsupported_columns=True)
    m = wf.read_mesh(path, crs="EPSG:31983", length_unit="m")
    assert m.crs == "EPSG:31983" and m.length_unit == "m"
    with pytest.raises(wf.InvalidInput, match="no polylines"):
        wf.read_polylines(path)
    with pytest.raises(wf.InvalidInput):
        wf.write(path, solid)
    with pytest.raises(wf.InvalidInput):
        wf.write(tmp_path / "x.stl", wf.datasets.drillholes(n=2))
    with pytest.raises(wf.InvalidInput):
        wf.write(tmp_path / "x.obj", "not data")


def test_binary_stl_warns_only_when_float32_moves_a_vertex(tmp_path):
    import warnings

    exact = wf.datasets.box([0, 0, 0], [1, 2, 3])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        wf.write(tmp_path / "exact.stl", exact)
        wf.write(tmp_path / "utm_ascii.stl", exact.translate([700000.1, 0, 0]), ascii=True)
    with pytest.warns(UserWarning, match=r"8 of 8 vertices moved.*ascii=True"):
        wf.write(tmp_path / "utm.stl", exact.translate([700000.1, 0, 0]))
