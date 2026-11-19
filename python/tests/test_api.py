"""Stub and runtime agree; errors map to the documented classes."""

import ast
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import wireframes as wf

STUB = Path(wf.__file__).with_name("_wireframes.pyi")


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
        assert hasattr(wf._wireframes, name), name
    for cls, names in members.items():
        for name in names:
            assert hasattr(getattr(wf._wireframes, cls), name), f"{cls}.{name}"


def test_every_public_runtime_name_is_in_the_stub():
    top, members = stub_names()
    public = {n for n in dir(wf._wireframes) if not n.startswith("__") or n == "__version__"}
    assert public - top == set()
    for cls, names in members.items():
        runtime = {n for n in vars(getattr(wf._wireframes, cls)) if not n.startswith("_")}
        runtime |= {
            n for n in ("__or__", "__and__", "__sub__", "__xor__") if n in vars(getattr(wf._wireframes, cls))
        }
        assert runtime - names == set(), cls


def test_every_exported_name_exists():
    for name in wf.__all__:
        assert getattr(wf, name) is not None, name
    assert wf.__version__.count(".") == 2


def test_unknown_attribute_raises():
    with pytest.raises(AttributeError):
        wf.no_such_name  # noqa: B018


def test_errors_subclass_builtins():
    assert issubclass(wf.InvalidMesh, wf.InvalidInput)
    assert issubclass(wf.InvalidInput, ValueError)
    assert issubclass(wf.FileError, OSError)
    assert issubclass(wf.MissingColumn, KeyError)
    assert all(issubclass(e, wf.WireframesError) for e in (wf.InvalidInput, wf.FileError, wf.MissingColumn))


def test_each_rust_error_maps_to_its_class(tmp_path, box):
    b = box([0, 0, 0], [1, 1, 1])
    with pytest.raises(wf.InvalidInput):
        wf.Mesh([[0, 0, np.inf]], [])
    with pytest.raises(wf.MissingColumn):
        b["nope"]
    with pytest.raises(wf.FileError):
        wf.read(tmp_path / "missing.obj")
    bad = tmp_path / "bad.off"
    bad.write_text("OFF\n3 1 0\n0 0 0\n1 0 0\n")
    with pytest.raises(wf.FileError):
        wf.read(bad)
    with pytest.raises(wf.InvalidInput, match="unknown file extension"):
        wf.read(tmp_path / "x.unknown")
    with pytest.raises(wf.InvalidInput, match="crs mismatch"):
        wf.Mesh(b.coords, b.triangles, crs="EPSG:1") | wf.Mesh(b.coords + 0.5, b.triangles, crs="EPSG:2")
    with pytest.raises(wf.InvalidInput, match="length unit mismatch"):
        wf.concat(
            [
                wf.Mesh(b.coords, b.triangles, length_unit="m"),
                wf.Mesh(b.coords, b.triangles, length_unit="ft"),
            ]
        )
    with pytest.raises(wf.InvalidMesh) as e:
        b.filter_faces(np.arange(12) > 0).volume  # noqa: B018
    assert isinstance(e.value.report, wf.MeshReport)
    assert e.value.report.summary["boundary_edges"] == 3


def test_matplotlib_stays_optional():
    code = "import sys, wireframes; wireframes.Mesh([[0,0,0]], []); assert 'matplotlib' not in sys.modules"
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
    p = wf.datasets.drillholes(n=3)
    assert p.closed.dtype == np.bool_ and p.feature.dtype == np.int64
    assert p.locate(p.coords).dtype == np.int64


def test_every_dataset_has_a_runnable_example():
    import doctest

    from wireframes import datasets

    finder = doctest.DocTestFinder(recurse=False)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    for name in ("box", "orebody", "topography", "pit", "drillholes"):
        tests = [
            t for t in finder.find(getattr(datasets, name), name, globs={"wf": wf, "np": np}) if t.examples
        ]
        assert tests, f"datasets.{name} has no example"
        assert all(runner.run(t, out=lambda s: None).failed == 0 for t in tests), name


def test_every_container_class_has_a_runnable_example():
    import doctest
    import inspect

    finder = doctest.DocTestFinder(recurse=False)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    for cls in (wf.Mesh, wf.Polylines, wf.PointSet, wf.Table):
        tests = [t for t in finder.find(cls, cls.__name__, globs={"wf": wf, "np": np}) if t.examples]
        assert tests, f"{cls.__name__} has no example"
        assert all(runner.run(t, out=lambda s: None).failed == 0 for t in tests), cls.__name__
    assert inspect.signature(wf.Polylines).parameters["closed"].default is False
