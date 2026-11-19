"""DXF, VTK, VTU/VTP and Datamine pairs through the Python API."""

from pathlib import Path

import numpy as np
import pytest
import wireframes as wf

FIXTURES = Path(__file__).resolve().parents[2] / "crates" / "io" / "tests" / "fixtures"


@pytest.fixture
def lens():
    return wf.datasets.orebody(points=8, sections=4)


@pytest.mark.parametrize("name", ["m.dxf", "m.vtk", "m.vtu", "m.vtp", "m_tr.csv"])
def test_mesh_round_trip(tmp_path, lens, name):
    wf.write(tmp_path / name, lens, skip_unsupported_columns=True)
    back = wf.read_mesh(tmp_path / name)
    assert np.array_equal(back.coords[back.triangles], lens.coords[lens.triangles])
    if name != "m.dxf":
        assert np.array_equal(back.coords, lens.coords)
    assert back.volume == lens.volume


def test_dxf_polyface_and_vtu_zlib(tmp_path, lens):
    wf.write(tmp_path / "p.dxf", lens, dxf_entity="polyface", skip_unsupported_columns=True)
    assert wf.read_mesh(tmp_path / "p.dxf").n_faces == lens.n_faces
    wf.write(tmp_path / "z.vtu", lens, zlib=True)
    assert wf.read_mesh(tmp_path / "z.vtu")["section"].tolist() == lens["section"].tolist()
    with pytest.raises(wf.InvalidInput):
        wf.write(tmp_path / "q.dxf", lens, dxf_entity="mesh")


def test_dxf_vtk_polylines_and_points_round_trip(tmp_path):
    holes = wf.datasets.drillholes(n=4)
    points = holes.to_points()
    for name in ("h.dxf", "h.vtk", "h.vtp"):
        wf.write(tmp_path / name, holes, skip_unsupported_columns=True)
        back = wf.read_polylines(tmp_path / name)
        assert np.array_equal(back.coords, holes.coords) and len(back) == len(holes)
        wf.write(tmp_path / f"p{name}", points, skip_unsupported_columns=True)
        assert np.array_equal(wf.read_points(tmp_path / f"p{name}").coords, points.coords)


def test_vtk_dataset_round_trip(tmp_path, lens):
    data = wf.Dataset(lens, wf.datasets.drillholes(n=2), None, None, None, [])
    wf.write(tmp_path / "all.vtp", data, skip_unsupported_columns=True)
    back = wf.read(tmp_path / "all.vtp")
    assert back.mesh.n_faces == lens.n_faces and len(back.polylines) == 2 and back.warnings == []


def test_vtu_grid_fixtures():
    data = wf.read(FIXTURES / "two_pieces.vtu")
    grid = data.grid if data.grid is not None else None
    assert grid is not None and grid.n_cells >= 2
    assert grid.cell_types.dtype == np.uint8 and grid.offsets[0] == 0
    assert len(grid.offsets) == grid.n_cells + 1
    assert grid.points.shape[1] == 3 and "k" in grid.cell_data.column_names
    surface = grid.surface()
    assert surface.n_vertices == len(grid.points)
    poly = wf.read(FIXTURES / "polyhedron.vtu").grid
    assert poly is not None and poly.n_cells >= 1


def test_vtu_vtk_grid_round_trips(tmp_path):
    grid = wf.read(FIXTURES / "two_pieces.vtu").grid
    for name in ("g.vtu", "g.vtk"):
        wf.write(tmp_path / name, grid)
        back = wf.read(tmp_path / name).grid
        assert np.array_equal(back.cell_types, grid.cell_types)
        assert np.array_equal(back.connectivity, grid.connectivity)
        assert np.array_equal(back.offsets, grid.offsets)
        assert back.cell_data["k"].tolist() == grid.cell_data["k"].tolist()


def test_dxf_unit_and_cp1252():
    data = wf.read(FIXTURES / "cp1252_layers.dxf")
    assert data.length_unit == "ft"
    objects = [o for o in (data.mesh, data.polylines, data.points) if o is not None]
    assert objects and all(o.length_unit == "ft" for o in objects)
    with pytest.raises(wf.InvalidInput, match="disagrees"):
        wf.read(FIXTURES / "cp1252_layers.dxf", length_unit="m")
    assert wf.read(FIXTURES / "cp1252_layers.dxf", length_unit="foot").length_unit == "ft"


def test_datamine_pair_from_either_name():
    tr = wf.read_mesh(FIXTURES / "gapped_semicolon_tr.csv")
    pt = wf.read_mesh(FIXTURES / "gapped_semicolon_pt.csv")
    assert np.array_equal(tr.coords, pt.coords) and tr.n_faces > 0
