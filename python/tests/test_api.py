"""Stub and runtime agree; errors map to the documented classes."""

import ast
import subprocess
import sys
from pathlib import Path

import drusa as dr
import numpy as np
import pytest

STUB = Path(dr.__file__).with_name("_drusa.pyi")


def stub_names():
    tree = ast.parse(STUB.read_text())
    top, members = set(), {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            top.add(node.name)
        if isinstance(node, ast.ClassDef):
            members[node.name] = {n.name for n in node.body if isinstance(n, ast.FunctionDef)}
        is_alias = isinstance(node, ast.AnnAssign) and ast.unparse(node.annotation) == "TypeAlias"
        if isinstance(node, ast.AnnAssign) and not is_alias:
            top.add(node.target.id)
    return top, members


def test_every_stub_name_exists_at_runtime():
    top, members = stub_names()
    for name in top:
        assert hasattr(dr._drusa, name), name
    for cls, names in members.items():
        for name in names:
            assert hasattr(getattr(dr._drusa, cls), name), f"{cls}.{name}"


def test_every_public_runtime_name_is_in_the_stub():
    top, members = stub_names()
    public = {n for n in dir(dr._drusa) if not n.startswith("__") or n == "__version__"}
    assert public - top == set()
    for cls, names in members.items():
        runtime = {n for n in vars(getattr(dr._drusa, cls)) if not n.startswith("_")}
        runtime |= {
            n for n in ("__or__", "__and__", "__sub__", "__xor__") if n in vars(getattr(dr._drusa, cls))
        }
        assert runtime - names == set(), cls


def test_every_exported_name_exists():
    for name in dr.__all__:
        assert getattr(dr, name) is not None, name
    assert dr.__version__.count(".") == 2


def test_unknown_attribute_raises():
    with pytest.raises(AttributeError):
        dr.no_such_name  # noqa: B018


def test_errors_subclass_builtins():
    assert issubclass(dr.InvalidMesh, dr.InvalidInput)
    assert issubclass(dr.InvalidInput, ValueError)
    assert issubclass(dr.FileError, OSError)
    assert issubclass(dr.MissingColumn, KeyError)
    assert all(issubclass(e, dr.DrusaError) for e in (dr.InvalidInput, dr.FileError, dr.MissingColumn))


def test_each_rust_error_maps_to_its_class(tmp_path, box):
    b = box([0, 0, 0], [1, 1, 1])
    with pytest.raises(dr.InvalidInput):
        dr.Mesh([[0, 0, np.inf]], [])
    with pytest.raises(dr.MissingColumn):
        b["nope"]
    with pytest.raises(dr.FileError):
        dr.read(tmp_path / "missing.obj")
    bad = tmp_path / "bad.off"
    bad.write_text("OFF\n3 1 0\n0 0 0\n1 0 0\n")
    with pytest.raises(dr.FileError):
        dr.read(bad)
    with pytest.raises(dr.InvalidInput, match="unknown file extension"):
        dr.read(tmp_path / "x.unknown")
    with pytest.raises(dr.InvalidInput, match="crs mismatch"):
        dr.Mesh(b.coords, b.triangles, crs="EPSG:1") | dr.Mesh(b.coords + 0.5, b.triangles, crs="EPSG:2")
    with pytest.raises(dr.InvalidInput, match="length unit mismatch"):
        dr.concat(
            [
                dr.Mesh(b.coords, b.triangles, length_unit="m"),
                dr.Mesh(b.coords, b.triangles, length_unit="ft"),
            ]
        )
    with pytest.raises(dr.InvalidMesh) as e:
        b.filter_faces(np.arange(12) > 0).volume  # noqa: B018
    assert isinstance(e.value.report, dr.MeshReport)
    assert e.value.report.summary["boundary_edges"] == 3


def test_matplotlib_stays_optional():
    code = "import sys, drusa; drusa.Mesh([[0,0,0]], []); assert 'matplotlib' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)


def test_property_dtypes_and_shapes(box, orebody):
    m = orebody
    for name, dtype, shape in [
        ("coords", np.float64, (m.n_vertices, 3)),
        ("triangles", np.int64, (m.n_faces, 3)),
        ("x", np.float64, (m.n_vertices,)),
        ("face_areas", np.float64, (m.n_faces,)),
        ("face_normals", np.float64, (m.n_faces, 3)),
        ("face_centroids", np.float64, (m.n_faces, 3)),
        ("vertex_normals", np.float64, (m.n_vertices, 3)),
        ("dip", np.float64, (m.n_faces, 2)),
        ("centroid", np.float64, (3,)),
    ]:
        value = getattr(m, name)
        assert value.dtype == dtype and value.shape == shape, name
    assert m.locate(m.coords).dtype == np.int8
    assert m.contains(m.coords).dtype == np.bool_
    assert m.component_labels().dtype == np.int64
    p = dr.datasets.drillholes(n=3)
    assert p.closed.dtype == np.bool_ and p.feature.dtype == np.int64
    assert p.locate(p.coords).dtype == np.int64
