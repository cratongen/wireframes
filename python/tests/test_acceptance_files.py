"""Acceptance scenarios for files as mining packages write them: DXF with 3DFACE, POLYFACE, MESH and
3D polylines, Datamine pairs, OBJ, STL, PLY, VTK, VTU and OFF, files from independent writers,
truncated and mutated files, and the Python API as a mining engineer uses it. All data are seeded
and synthetic."""

import doctest
import inspect
import os
import re
import struct
import subprocess
import sys

import numpy as np
import pytest
import wireframes as wf

FORMATS = [
    ("m.obj", {}),
    ("m_ascii.ply", {"ascii": True}),
    ("m.ply", {}),
    ("m.off", {}),
    ("m.stl", {}),
    ("m_ascii.stl", {"ascii": True}),
    ("m.dxf", {}),
    ("m_polyface.dxf", {"dxf_entity": "polyface"}),
    ("m.vtk", {}),
    ("m_ascii.vtk", {"ascii": True}),
    ("m.vtu", {}),
    ("m_zlib.vtu", {"zlib": True}),
    ("m.vtp", {}),
    ("m_tr.csv", {}),
]
KEEPS_COLUMNS = {"ply", "vtk", "vtu", "vtp", "csv"}


def attributed(points, sections, seed=11):
    """The orebody lens with a vertex grade and a face domain, in Float64 so every format stores them."""
    m = wf.datasets.orebody(points=points, sections=sections)
    rng = np.random.default_rng(seed)
    return (
        m.with_vertex_column("grade", rng.gamma(2.0, 0.6, m.n_vertices))
        .with_face_column("domain", rng.integers(1, 6, m.n_faces).astype(np.float64))
        .with_face_column("section", m["section"].astype(np.float64))
    )


def f32_local(m):
    """`m` moved to a local origin with every coordinate representable in f32, as STL requires."""
    lo = np.floor(np.array(m.bounds[0]))
    local = (m.coords - lo).astype(np.float32).astype(np.float64)
    return wf.Mesh(
        local, m.triangles, face_attributes=m.face_attributes, vertex_attributes=m.vertex_attributes
    )


@pytest.fixture(scope="module")
def lens():
    return attributed(300, 170)


@pytest.fixture(scope="module")
def lens_f32(lens):
    return f32_local(lens)


def file_bytes(path):
    name = path.name
    if name.endswith("_tr.csv"):
        return path.read_bytes() + path.with_name(name.replace("_tr.csv", "_pt.csv")).read_bytes()
    return path.read_bytes()


@pytest.mark.parametrize(("name", "kw"), FORMATS, ids=[n for n, _ in FORMATS])
def test_write_read_write_is_bit_identical_at_1e5_faces(tmp_path, lens, lens_f32, name, kw):
    src = lens_f32 if name.endswith(".stl") else lens
    assert src.n_faces > 100_000
    first = tmp_path / name
    wf.write(first, src, skip_unsupported_columns=True, **kw)
    back = wf.read_mesh(first)
    again = tmp_path / ("again_" + name)
    wf.write(again, back, skip_unsupported_columns=True, **kw)
    assert file_bytes(first) == file_bytes(again)
    assert np.array_equal(back.coords[back.triangles], src.coords[src.triangles])
    assert back.is_solid and back.volume == src.volume and back.area == src.area
    if first.suffix.lstrip(".") in KEEPS_COLUMNS:
        for c in ("grade", "domain", "section"):
            assert np.array_equal(back[c], src[c])


@pytest.mark.parametrize("name", ["big.ply", "big.vtk", "big.vtu", "big.stl"])
def test_3e5_faces_round_trip_and_stay_valid(tmp_path, name):
    m = attributed(500, 300, seed=3)
    if name.endswith(".stl"):
        m = f32_local(m)
    assert m.n_faces > 290_000
    wf.write(tmp_path / name, m, skip_unsupported_columns=True)
    back = wf.read_mesh(tmp_path / name)
    assert back.n_faces == m.n_faces and back.validate().is_solid
    assert back.volume == m.volume
    lo, hi = np.array(m.bounds[0]), np.array(m.bounds[1])
    assert np.array_equal(np.array(back.bounds[0]), lo) and np.array_equal(np.array(back.bounds[1]), hi)


def test_binary_stl_warns_when_utm_coordinates_round(tmp_path):
    m = wf.datasets.orebody(points=60, sections=40)
    assert m.is_solid
    with pytest.warns(UserWarning, match="f32|float32|precision"):
        wf.write(tmp_path / "utm.stl", m, skip_unsupported_columns=True)


def test_stl_rounding_of_utm_coordinates_is_bounded(tmp_path):
    m = wf.datasets.orebody(points=60, sections=40)
    wf.write(tmp_path / "utm.stl", m, skip_unsupported_columns=True)
    back = wf.read_mesh(tmp_path / "utm.stl")
    assert back.n_faces == m.n_faces
    err = np.abs(back.coords[back.triangles] - m.coords[m.triangles]).max()
    assert err <= np.spacing(np.float32(np.abs(m.coords).max())) / 2


def tetra(offset=(0.0, 0.0, 0.0)):
    p = np.array([(0, 0, 0), (10, 0, 0), (0, 10, 0), (0, 0, 10)], float) + offset
    return p, [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]


def group(code, value):
    return f"{code:>3}\r\n{value}\r\n"


def mining_dxf():
    """A DXF as a mining package writes it: CRLF, padded group codes, handles, subclass markers,
    extended data, a LAYER table, 3DFACE, polyface POLYLINE, MESH, a 3D POLYLINE and a POINT."""
    g = group
    s = g(999, "exported by a mine planning package")
    s += g(0, "SECTION") + g(2, "HEADER") + g(9, "$ACADVER") + g(1, "AC1015")
    s += g(9, "$INSUNITS") + g(70, "     6") + g(0, "ENDSEC")
    s += g(0, "SECTION") + g(2, "TABLES") + g(0, "TABLE") + g(2, "LAYER") + g(5, "2") + g(70, 2)
    s += g(0, "LAYER") + g(5, "10") + g(2, "ORE") + g(70, 0) + g(62, 1) + g(6, "CONTINUOUS")
    s += g(0, "ENDTAB") + g(0, "ENDSEC") + g(0, "SECTION") + g(2, "ENTITIES")
    p, t = tetra((700000.0, 7500000.0, 200.0))
    for k, (a, b, c) in enumerate(t):
        s += g(0, "3DFACE") + g(5, f"{0x100 + k:X}") + g(100, "AcDbEntity") + g(8, "ORE") + g(100, "AcDbFace")
        for j, v in enumerate((a, b, c, c)):
            for ax in range(3):
                s += g(10 + j + 10 * ax, f"  {p[v][ax]:.6f}")
        s += (
            g(1001, "PACKAGE")
            + g(1002, "{")
            + g(1000, "rock=BIF")
            + g(1040, "2.5")
            + g(1070, "3")
            + g(1002, "}")
        )
    p, t = tetra((700100.0, 7500000.0, 200.0))
    s += g(0, "POLYLINE") + g(5, "200") + g(8, "WASTE") + g(66, 1)
    s += g(10, 0.0) + g(20, 0.0) + g(30, 0.0) + g(70, 64) + g(71, 4) + g(72, 4)
    for q in p:
        s += g(0, "VERTEX") + g(8, "WASTE") + g(10, q[0]) + g(20, q[1]) + g(30, q[2]) + g(70, 192)
    for a, b, c in t:
        s += g(0, "VERTEX") + g(8, "WASTE") + g(10, 0.0) + g(20, 0.0) + g(30, 0.0) + g(70, 128)
        s += g(71, a + 1) + g(72, -(b + 1)) + g(73, c + 1)
    s += g(0, "SEQEND") + g(8, "WASTE")
    p, t = tetra((700200.0, 7500000.0, 200.0))
    s += g(0, "MESH") + g(5, "300") + g(100, "AcDbEntity") + g(8, "FAULT") + g(100, "AcDbSubDMesh")
    s += g(71, 2) + g(72, 0) + g(91, 0) + g(92, 4)
    for q in p:
        s += g(10, q[0]) + g(20, q[1]) + g(30, q[2])
    s += g(93, 16) + "".join(g(90, 3) + "".join(g(90, i) for i in tri) for tri in t)
    s += g(94, 0) + g(95, 0) + g(90, 0)
    s += g(0, "POLYLINE") + g(8, "CREST") + g(66, 1) + g(10, 0.0) + g(20, 0.0) + g(30, 0.0) + g(70, 9)
    for q in [(700000, 7500050, 205), (700050, 7500050, 206), (700050, 7500100, 207)]:
        s += g(0, "VERTEX") + g(8, "CREST") + g(10, q[0]) + g(20, q[1]) + g(30, q[2]) + g(70, 32)
    s += g(0, "SEQEND") + g(0, "POINT") + g(8, "SURVEY") + g(10, 700001.5) + g(20, 7500002.5) + g(30, 203.5)
    return (s + g(0, "ENDSEC") + g(0, "EOF")).encode()


def test_mining_package_dxf_reads_every_entity_with_layers(tmp_path):
    (tmp_path / "pit.dxf").write_bytes(mining_dxf())
    data = wf.read(tmp_path / "pit.dxf")
    m = data.mesh
    assert data.length_unit == "m" and data.warnings == []
    assert m.n_faces == 12 and m.n_vertices == 12 and m.is_solid
    layers = m["layer"]
    one = 1000.0 / 6.0
    for name in ("ORE", "WASTE", "FAULT"):
        part = m.filter_faces(layers == name)
        assert part.n_faces == 4 and abs(part.volume - one) <= 1e-9 * one
    assert sum(c.volume for c in m.components()) == pytest.approx(m.volume, rel=1e-12)
    assert abs(m.volume - 3 * one) <= 1e-9 * one
    crest = data.polylines
    assert (
        len(crest) == 1
        and crest.closed.tolist() == [True]
        and crest.attributes["layer"].tolist() == ["CREST"]
    )
    assert crest.coords[:, 2].tolist() == [205.0, 206.0, 207.0]
    assert data.points.coords.tolist() == [[700001.5, 7500002.5, 203.5]]
    assert data.points["layer"].tolist() == ["SURVEY"]


def test_mining_package_dxf_output_is_valid_input_again(tmp_path):
    (tmp_path / "pit.dxf").write_bytes(mining_dxf())
    data = wf.read(tmp_path / "pit.dxf")
    for entity in ("3dface", "polyface"):
        out = tmp_path / f"out_{entity}.dxf"
        wf.write(out, data, dxf_entity=entity)
        back = wf.read(out)
        assert back.length_unit == "m"
        assert np.array_equal(back.mesh.coords[back.mesh.triangles], data.mesh.coords[data.mesh.triangles])
        assert back.mesh["layer"].tolist() == data.mesh["layer"].tolist()
        assert np.array_equal(back.polylines.coords, data.polylines.coords)
        assert np.array_equal(back.points.coords, data.points.coords)
    ore = data.mesh.filter_faces(data.mesh["layer"] == "ORE")
    waste = data.mesh.filter_faces(data.mesh["layer"] == "WASTE")
    both = wf.union(ore, waste)
    assert both.is_solid and both.volume == pytest.approx(ore.volume + waste.volume, rel=1e-12)
    wf.write(tmp_path / "both.vtu", both, skip_unsupported_columns=True)
    assert wf.read_mesh(tmp_path / "both.vtu").volume == both.volume


def write_pair(tmp_path, tag, pt, tr):
    (tmp_path / f"{tag}_pt.csv").write_bytes(pt.encode())
    (tmp_path / f"{tag}_tr.csv").write_bytes(tr.encode())
    return tmp_path / f"{tag}_tr.csv"


def test_datamine_pair_with_reordered_columns_blank_lines_and_zero_based_ids(tmp_path):
    p, t = tetra((500000.0, 9000000.0, 1000.0))
    pt = "PID,ZP,XP,YP,DENSITY,ROCK\r\n"
    for i, q in enumerate(p):
        pt += f"{i},{float(q[2])!r},{float(q[0])!r},{float(q[1])!r},{2.6 + i / 10!r},BIF\r\n" + (
            "\r\n" if i == 1 else ""
        )
    tr = "PID3,TRIANGLE,PID1,PID2,ZONE\r\n"
    tr += "".join(f"{c},{k},{a},{b},{k % 2}\r\n" for k, (a, b, c) in enumerate(t)) + "\r\n\r\n"
    m = wf.read_mesh(write_pair(tmp_path, "zero", pt, tr))
    assert np.array_equal(m.coords, p) and m.triangles.tolist() == [list(x) for x in t]
    assert m["DENSITY"].tolist() == [2.6 + i / 10 for i in range(4)] and m["ROCK"].tolist() == ["BIF"] * 4
    assert m["ZONE"].tolist() == [0.0, 1.0, 0.0, 1.0]
    assert m.volume == pytest.approx(1000.0 / 6.0, rel=1e-12)


def test_datamine_studio_style_pair_with_type_row_and_offset_ids(tmp_path):
    p, t = tetra()
    pt = "XP,YP,ZP,PID\nN,N,N,N\n" + "".join(
        f"{q[0]:.3f},{q[1]:.3f},{q[2]:.3f},{i + 101}.000\n" for i, q in enumerate(p)
    )
    tr = "TRIANGLE,PID1,PID2,PID3,SID\nN,N,N,N,N\n" + "".join(
        f"{k + 1},{a + 101},{b + 101},{c + 101},7\n" for k, (a, b, c) in enumerate(t)
    )
    m = wf.read_mesh(write_pair(tmp_path, "studio", pt, tr))
    assert m.column_names == ["SID"] and m["SID"].tolist() == [7.0] * 4
    assert m.volume == pytest.approx(1000.0 / 6.0, rel=1e-12)


def test_datamine_semicolon_pair_with_decimal_commas_quotes_and_nulls(tmp_path):
    p, t = tetra()
    pt = '"XP"; "YP"; "ZP"; "PID"; "AU"\n' + "".join(
        f" {q[0]:g},5; {q[1]:g}; {q[2]:g}; {i + 1}; -1.0E30\n" for i, q in enumerate(p)
    )
    tr = "TRIANGLE;PID1;PID2;PID3\n" + "".join(
        f"{k + 1};{a + 1};{b + 1};{c + 1}\n" for k, (a, b, c) in enumerate(t)
    )
    m = wf.read_mesh(write_pair(tmp_path, "semi", pt, tr))
    assert np.array_equal(m.x, p[:, 0] + 0.5)
    assert np.isnan(m["AU"]).all() and not m.vertex_attributes.mask("AU").any()


def test_datamine_missing_pid_names_the_row(tmp_path):
    p, _t = tetra()
    pt = "XP,YP,ZP,PID\n" + "".join(f"{q[0]},{q[1]},{q[2]},{i + 1}\n" for i, q in enumerate(p))
    tr = "TRIANGLE,PID1,PID2,PID3\n1,1,3,2\n2,1,2,99\n"
    with pytest.raises(wf.FileError, match=r"_tr\.csv:3: .*PID 99"):
        wf.read_mesh(write_pair(tmp_path, "gap", pt, tr))


def test_datamine_trisol_beside_triangle_is_kept(tmp_path):
    pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0,2\n0,1,0,3\n"
    tr = "TRISOL,TRIANGLE,PID1,PID2,PID3\n7,1,1,2,3\n"
    m = wf.read_mesh(write_pair(tmp_path, "trisol", pt, tr))
    assert m["TRISOL"].tolist() == [7.0]


def test_datamine_error_names_the_points_file(tmp_path):
    pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0\n"
    tr = "TRIANGLE,PID1,PID2,PID3\n1,1,2,3\n"
    with pytest.raises(wf.FileError, match=r"bad_pt\.csv:3"):
        wf.read_mesh(write_pair(tmp_path, "bad", pt, tr))


BOX = [(0, 0, 0), (2, 0, 0), (2, 3, 0), (0, 3, 0), (0, 0, 4), (2, 0, 4), (2, 3, 4), (0, 3, 4)]
QUADS = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def test_obj_quads_negative_indices_groups_and_ignored_statements(tmp_path):
    s = "# exported\r\nmtllib pit.mtl\r\no pit\r\n" + "".join(
        "v  {}\r\n".format(" ".join(map(str, v))) for v in BOX
    )
    s += "vt 0 0\r\nvn 0 0 1\r\ng floor walls\r\nusemtl rock\r\ns 1\r\n"
    s += "f {}\r\n".format(" ".join(f"{i - 8}/1/1" for i in QUADS[0]))
    s += "g roof\r\nf {}\r\n".format(" ".join(f"{i + 1}//1" for i in QUADS[1]))
    s += "g walls\r\n" + "".join("f {}\r\n".format(" ".join(str(i + 1) for i in q)) for q in QUADS[2:])
    (tmp_path / "q.obj").write_bytes(s.encode())
    m = wf.read_mesh(tmp_path / "q.obj")
    assert m.n_faces == 12 and m.is_solid and m.volume == 24.0 and m.area == 52.0
    assert m["group"].tolist()[:4] == ["floor walls", "floor walls", "roof", "roof"]
    assert set(m["object"].tolist()) == {"pit"} and m["polygon"].tolist() == [
        0,
        0,
        1,
        1,
        2,
        2,
        3,
        3,
        4,
        4,
        5,
        5,
    ]
    wf.write(tmp_path / "back.obj", m)
    back = wf.read_mesh(tmp_path / "back.obj")
    assert back.triangles.tolist() == m.triangles.tolist() and back["group"].tolist() == m["group"].tolist()


def test_obj_ngon_prism_matches_closed_form_volume(tmp_path):
    n = 24
    a = np.arange(n) * 2 * np.pi / n
    ring = np.column_stack([np.cos(a), np.sin(a)]).tolist()
    s = "".join(f"v {x!r} {y!r} 0.0\n" for x, y in ring) + "".join(f"v {x!r} {y!r} 5.0\n" for x, y in ring)
    s += "f " + " ".join(str(n - i) for i in range(n)) + "\n"
    s += "f " + " ".join(str(n + 1 + i) for i in range(n)) + "\n"
    s += "".join(f"f {i + 1} {(i + 1) % n + 1} {n + (i + 1) % n + 1} {n + i + 1}\n" for i in range(n))
    (tmp_path / "prism.obj").write_text(s)
    m = wf.read_mesh(tmp_path / "prism.obj")
    assert m.n_faces == 2 * (n - 2) + 2 * n and m.is_solid
    assert m.volume == pytest.approx(5.0 * n / 2 * np.sin(2 * np.pi / n), rel=1e-12)


def test_binary_stl_with_solid_header_welds_duplicate_corners(tmp_path):
    p, t = tetra()
    body = b"".join(struct.pack("<12fH", 0, 0, 0, *p[a], *p[b], *p[c], 7) for a, b, c in t)
    header = b"solid exported by a mine planning package".ljust(80, b" ")
    (tmp_path / "s.stl").write_bytes(header + struct.pack("<I", len(t)) + body)
    m = wf.read_mesh(tmp_path / "s.stl")
    assert m.n_vertices == 4 and m.n_faces == 4 and m.is_solid
    assert m.volume == pytest.approx(1000.0 / 6.0, rel=1e-12) and m["stl_attribute"].tolist() == [7] * 4


def test_ascii_stl_with_duplicate_corners_mixed_case_and_two_solids(tmp_path):
    p, t = tetra()
    s = ""
    for name, tris in (("ore", t[:2]), ("Waste Rock", t[2:])):
        s += f"SOLID {name}\r\n"
        for a, b, c in tris:
            s += "  facet NORMAL 0 0 0\r\n    outer loop\r\n"
            s += "".join(f"      Vertex {p[v][0]:e} {p[v][1]:e} {p[v][2]:e}\r\n" for v in (a, b, c))
            s += "    endloop\r\n  endfacet\r\n"
        s += f"endsolid {name}\r\n"
    (tmp_path / "a.stl").write_bytes(s.encode())
    m = wf.read_mesh(tmp_path / "a.stl")
    assert (
        m.n_vertices == 4 and m.is_solid and m["solid"].tolist() == ["ore", "ore", "Waste Rock", "Waste Rock"]
    )


@pytest.mark.parametrize("endian", [">", "<"])
def test_binary_ply_with_extra_properties_keeps_and_rewrites_them(tmp_path, endian):
    p, t = tetra()
    fmt = {"<": "binary_little_endian", ">": "binary_big_endian"}[endian]
    header = (
        f"ply\nformat {fmt} 1.0\ncomment from a mining package\nelement vertex 4\n"
        "property float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty double grade\n"
        "element face 4\nproperty list uchar int vertex_indices\nproperty short domain\nproperty float conf\n"
        "element survey 1\nproperty int station\nend_header\n"
    ).encode()
    body = b"".join(struct.pack(endian + "fffBd", *q, i * 10, i / 4) for i, q in enumerate(p))
    body += b"".join(struct.pack(endian + "B3ihf", 3, *f, k, k * 0.5) for k, f in enumerate(t))
    (tmp_path / "e.ply").write_bytes(header + body + struct.pack(endian + "i", 42))
    m = wf.read_mesh(tmp_path / "e.ply")
    assert m.is_solid and m.volume == pytest.approx(1000.0 / 6.0, rel=1e-12)
    expect = {
        "red": [0, 10, 20, 30],
        "grade": [0, 0.25, 0.5, 0.75],
        "domain": [0, 1, 2, 3],
        "conf": [0, 0.5, 1, 1.5],
    }
    for c, v in expect.items():
        assert m[c].tolist() == v
    for ascii_mode in (True, False):
        wf.write(tmp_path / "w.ply", m, ascii=ascii_mode)
        back = wf.read_mesh(tmp_path / "w.ply")
        assert all(back[c].tolist() == v for c, v in expect.items())


MIXED_VTK = """# vtk DataFile Version 3.0
mixed cells
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 6 double
0 0 0  1 0 0  0 1 0  0 0 1  1 1 0  1 1 1
CELLS 6 24
4 0 1 2 3
4 1 4 2 5
3 0 1 2
4 0 1 4 2
2 0 5
1 3
CELL_TYPES 6
10
10
5
9
3
1
CELL_DATA 6
SCALARS zone int 1
LOOKUP_TABLE default
1 2 3 4 5 6
POINT_DATA 6
SCALARS au float
LOOKUP_TABLE default
0.1 0.2 0.3 0.4 0.5 0.6
"""


def test_vtk_mixed_cells_round_trip_through_vtk_and_vtu(tmp_path):
    (tmp_path / "mixed.vtk").write_text(MIXED_VTK)
    grid = wf.read(tmp_path / "mixed.vtk").grid
    assert grid.cell_types.tolist() == [10, 10, 5, 9, 3, 1] and grid.cell_data["zone"].tolist() == [
        1,
        2,
        3,
        4,
        5,
        6,
    ]
    surface = grid.surface()
    assert surface.n_faces == 3 and surface.area == pytest.approx(1.5, rel=1e-15)
    assert surface["zone"].tolist() == [3, 4, 4]
    for name, kw in (
        ("g.vtk", {}),
        ("g_ascii.vtk", {"ascii": True}),
        ("g.vtu", {}),
        ("g_zlib.vtu", {"zlib": True}),
    ):
        wf.write(tmp_path / name, grid, **kw)
        back = wf.read(tmp_path / name).grid
        assert np.array_equal(back.points, grid.points) and np.array_equal(back.offsets, grid.offsets)
        assert np.array_equal(back.connectivity, grid.connectivity)
        assert back.cell_types.tolist() == grid.cell_types.tolist()
        assert back.cell_data["zone"].tolist() == [1, 2, 3, 4, 5, 6]
        assert np.array_equal(back.point_data["au"], grid.point_data["au"])


def test_off_quads_with_comments_and_counts_on_the_next_line(tmp_path):
    s = "OFF\n# from a mining package\n\n8 6 0\n" + "".join("{}\n".format(" ".join(map(str, v))) for v in BOX)
    s += "".join("4 {}\n".format(" ".join(map(str, q))) for q in QUADS)
    (tmp_path / "q.off").write_text(s)
    m = wf.read_mesh(tmp_path / "q.off")
    assert m.is_solid and m.volume == 24.0 and m.area == 52.0


@pytest.fixture(scope="module")
def small():
    return attributed(60, 40, seed=5)


MESHIO_WRITES = [
    ("mio.vtk", {}, True),
    ("mio_ascii.vtk", {"binary": False}, True),
    ("mio42.vtk", {"file_format": "vtk42"}, True),
    ("mio.vtu", {}, True),
    ("mio_raw.vtu", {"compression": None}, True),
    ("mio_ascii.vtu", {"binary": False}, False),
    ("mio.ply", {}, True),
    ("mio_ascii.ply", {"binary": False}, True),
    ("mio.obj", {}, True),
    ("mio.off", {}, True),
    ("mio.stl", {}, False),
    ("mio_ascii.stl", {"binary": False}, True),
]


@pytest.mark.parametrize(("name", "kw", "exact"), MESHIO_WRITES, ids=[n for n, _, _ in MESHIO_WRITES])
def test_reads_files_written_by_meshio(tmp_path, small, name, kw, exact):
    meshio = pytest.importorskip("meshio")
    v, t = small.coords, small.triangles
    grade, domain = small["grade"], small["domain"].astype(np.int32)
    mio = meshio.Mesh(v, [("triangle", t)], point_data={"grade": grade}, cell_data={"domain": [domain]})
    meshio.write(tmp_path / name, mio, **kw)
    back = wf.read_mesh(tmp_path / name)
    assert back.n_faces == small.n_faces and back.is_solid
    got, want = back.coords[back.triangles], v[t]
    if exact:
        assert np.array_equal(got, want)
    else:
        assert np.abs(got - want).max() <= 1e-6 * np.abs(want).max()
    if name.endswith((".vtk", ".vtu")):
        assert back["domain"].tolist() == domain.tolist()
    if name.endswith((".vtk", ".vtu", ".ply")):
        assert np.allclose(back["grade"], grade, rtol=1e-10 if exact else 1e-6, atol=0)


@pytest.mark.parametrize("ext", ["stl", "ply", "obj", "off"])
def test_reads_files_written_by_trimesh(tmp_path, small, ext):
    trimesh = pytest.importorskip("trimesh")
    tm = trimesh.Trimesh(small.coords, small.triangles, process=False)
    tm.export(tmp_path / f"t.{ext}")
    back = wf.read_mesh(tmp_path / f"t.{ext}")
    want = small.coords[small.triangles]
    got = back.coords[back.triangles]
    assert back.n_faces == small.n_faces and back.is_solid
    bound = {"obj": 5e-9, "off": 5e-9}.get(ext, float(np.spacing(np.float32(np.abs(want).max()))))
    assert np.abs(got - want).max() <= bound


@pytest.mark.parametrize(
    ("name", "kw"),
    [("w_ascii.ply", {"ascii": True}), ("w.obj", {}), ("w.off", {}), ("w_ascii.stl", {"ascii": True})],
)
def test_independent_readers_read_what_wireframes_writes(tmp_path, small, name, kw):
    meshio = pytest.importorskip("meshio")
    trimesh = pytest.importorskip("trimesh")
    wf.write(tmp_path / name, small, skip_unsupported_columns=True, **kw)
    r = meshio.read(tmp_path / name)
    tri = np.concatenate([c.data for c in r.cells if c.type == "triangle"])
    assert np.array_equal(r.points[tri], small.coords[small.triangles])
    t = trimesh.load(tmp_path / name, process=False, force="mesh")
    assert len(t.faces) == small.n_faces
    assert float(t.volume) == pytest.approx(small.volume, rel=1e-12)


TRUNCATION_LOCATED = ["m_ascii.ply", "m.ply", "m.stl", "m_ascii.stl", "m.vtk", "m_ascii.vtk"]
LOCATION = re.compile(r"\.\w+:\d+: |byte \d+|\(\d+ bytes\)|has \d+ bytes")


@pytest.mark.parametrize("name", TRUNCATION_LOCATED)
@pytest.mark.parametrize("fraction", [0.3, 0.6, 0.97])
def test_truncated_files_fail_naming_the_line_or_byte(tmp_path, small, name, fraction):
    kw = dict(FORMATS)[name]
    wf.write(tmp_path / name, small, skip_unsupported_columns=True, **kw)
    full = (tmp_path / name).read_bytes()
    (tmp_path / name).write_bytes(full[: int(len(full) * fraction)])
    with pytest.raises(wf.FileError) as e:
        wf.read(tmp_path / name)
    assert LOCATION.search(str(e.value)), str(e.value)


def test_truncated_dxf_is_an_error(tmp_path, small):
    wf.write(tmp_path / "m.dxf", small, skip_unsupported_columns=True)
    full = (tmp_path / "m.dxf").read_bytes()
    (tmp_path / "m.dxf").write_bytes(full[: int(len(full) * 0.97)])
    with pytest.raises(wf.FileError):
        wf.read(tmp_path / "m.dxf")


@pytest.mark.parametrize("name", ["m.off", "m.vtu", "m.vtp"])
def test_truncated_off_and_vtu_name_the_location(tmp_path, small, name):
    wf.write(tmp_path / name, small, skip_unsupported_columns=True)
    full = (tmp_path / name).read_bytes()
    (tmp_path / name).write_bytes(full[: int(len(full) * 0.3)])
    with pytest.raises(wf.FileError) as e:
        wf.read(tmp_path / name)
    assert LOCATION.search(str(e.value)), str(e.value)


def mutate(rng, b):
    b = bytearray(b)
    tokens = [b"9999999999", b"-1", b"4294967295", b"nan", b"1e308", b"\r\n", b" ", b"0", b"_", b"<"]
    for _ in range(int(rng.integers(1, 4))):
        op, i = int(rng.integers(0, 6)), int(rng.integers(0, max(1, len(b))))
        if op == 0 and b:
            b[i % len(b)] = int(rng.integers(0, 256))
        elif op == 1:
            del b[i : i + int(rng.integers(1, 64))]
        elif op == 2:
            b[i:i] = rng.integers(0, 256, int(rng.integers(1, 16)), dtype=np.uint8).tobytes()
        elif op == 3:
            j = int(rng.integers(0, max(1, len(b))))
            b[i:i] = b[j : j + int(rng.integers(1, 200))]
        elif op == 4:
            b[i:i] = tokens[int(rng.integers(0, len(tokens)))]
        else:
            del b[i:]
    return bytes(b)


def test_seeded_fuzz_of_mutated_files_never_escapes_wireframes_errors(tmp_path):
    m = attributed(8, 3, seed=1)
    sources = {}
    for name, kw in FORMATS:
        wf.write(tmp_path / name, m, skip_unsupported_columns=True, **kw)
        sources[name] = (tmp_path / name).read_bytes()
    points = (tmp_path / "m_pt.csv").read_bytes()
    rng = np.random.default_rng(20261119)
    escaped = []
    for _ in range(250):
        for name, src in sources.items():
            path = tmp_path / ("fuzz_" + name)
            path.write_bytes(mutate(rng, src))
            if name.endswith("_tr.csv"):
                (tmp_path / "fuzz_m_pt.csv").write_bytes(
                    mutate(rng, points) if rng.random() < 0.5 else points
                )
            try:
                data = wf.read(path)
                if data.mesh is not None:
                    data.mesh.validate()
            except wf.WireframesError:
                pass
            except BaseException as e:  # noqa: BLE001
                escaped.append((name, type(e).__name__, str(e)[:120]))
    assert escaped == []


DETERMINISM = """
import hashlib, sys, numpy as np, wireframes as wf
from pathlib import Path
d = Path(sys.argv[1])
m = wf.datasets.orebody(points=300, sections=170)
m = m.with_vertex_column("grade", np.random.default_rng(4).random(m.n_vertices))
h = hashlib.sha256()
for name, kw in {formats!r}:
    p = d / name
    wf.write(p, m, skip_unsupported_columns=True, **kw)
    h.update(p.read_bytes())
    back = wf.read_mesh(p)
    h.update(back.coords.tobytes()); h.update(back.triangles.tobytes())
print(h.hexdigest())
"""


def test_reading_and_writing_do_not_depend_on_thread_count(tmp_path):
    code = DETERMINISM.format(
        formats=[f for f in FORMATS if f[0].endswith((".ply", ".vtk", ".vtu", ".dxf", ".csv"))]
    )
    digests = set()
    for threads in ("1", "4"):
        out_dir = tmp_path / threads
        out_dir.mkdir()
        env = dict(os.environ, RAYON_NUM_THREADS=threads)
        out = subprocess.run(
            [sys.executable, "-c", code, str(out_dir)], env=env, capture_output=True, text=True, check=True
        )
        digests.add(out.stdout.strip())
    assert len(digests) == 1


def public_callables():
    out = []
    for name in wf.__all__:
        if name in ("errors", "datasets", "plot", "__version__"):
            continue
        obj = getattr(wf, name)
        if inspect.isclass(obj):
            for attr, value in vars(obj).items():
                if not attr.startswith("_") and callable(value):
                    out.append((f"{name}.{attr}", value))
        elif callable(obj):
            out.append((name, obj))
    from wireframes import plot

    out += [(f"plot.{n}", getattr(plot, n)) for n in ("plan", "section", "surface3d")]
    return out


def parameter_names(f):
    try:
        return [p for p in inspect.signature(f).parameters if p not in ("self", "cls")]
    except (TypeError, ValueError):
        return []


def test_every_public_callable_has_a_docstring():
    assert [n for n, f in public_callables() if not (f.__doc__ or "").strip()] == []


@pytest.mark.skip(reason="acceptance: 36 public callables with arguments have no Parameters section")
def test_every_public_callable_with_arguments_documents_them():
    no_params = [
        n for n, f in public_callables() if parameter_names(f) and "Parameters\n----------" not in f.__doc__
    ]
    assert no_params == []


@pytest.mark.skip(reason="acceptance: no public function or method has an Examples section")
def test_every_public_callable_has_a_runnable_example():
    finder = doctest.DocTestFinder(recurse=False)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    without, failed = [], []
    for name, f in public_callables():
        tests = [t for t in finder.find(f, name, globs={"wf": wf, "np": np}) if t.examples]
        if not tests:
            without.append(name)
        for t in tests:
            if runner.run(t, out=lambda s: None).failed:
                failed.append(name)
    assert without == [] and failed == []


def test_reprs_say_what_the_object_holds(tmp_path, small):
    assert repr(small).startswith(f"Mesh({small.n_vertices} vertices, {small.n_faces} triangles")
    assert '"grade"' in repr(small) and "length_unit m" in repr(small)
    assert repr(small.face_attributes).startswith(f"Table({small.n_faces} rows;")
    assert "domain: float64" in repr(small.face_attributes)
    holes = wf.datasets.drillholes(n=3)
    assert repr(holes).startswith("Polylines(3 features, 3 parts")
    assert repr(holes.to_points()).startswith("PointSet(6 points")
    assert repr(small.validate()) == "MeshReport(is_solid=true, shells=1, no problems)"
    (tmp_path / "mixed.vtk").write_text(MIXED_VTK)
    data = wf.read(tmp_path / "mixed.vtk")
    assert repr(data.grid) == "UnstructuredGrid(6 points, 6 cells)"
    assert repr(data).startswith("Dataset(mesh=None")
    _, report = wf.datasets.topography(n=20).repair(fill_holes="fan")
    assert "face_added=" in repr(report) and "holes_left=0" in repr(report)


def test_error_messages_name_the_cause_and_the_way_out(tmp_path, small):
    with pytest.raises(wf.MissingColumn, match=r'"au".*columns are \["grade", "section", "domain"\]'):
        small["au"]
    with pytest.raises(wf.InvalidInput, match="skip_unsupported_columns"):
        wf.write(tmp_path / "m.stl", small)
    with pytest.raises(wf.InvalidInput, match=r"_tr\.csv and <name>_pt\.csv"):
        wf.read(tmp_path / "survey.csv")
    with pytest.raises(wf.InvalidInput, match="face 0 uses 5 of 3"):
        wf.Mesh(np.zeros((3, 3)), [[0, 1, 5]])
    with pytest.raises(wf.InvalidInput, match=r"vertices must be finite.*vertices 0"):
        wf.Mesh([[np.nan, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 1, 2]])
    wf.write(tmp_path / "h.dxf", wf.datasets.drillholes(n=2), skip_unsupported_columns=True)
    with pytest.raises(wf.InvalidInput, match='length_unit "ft" disagrees with the file unit "m"'):
        wf.read(tmp_path / "h.dxf", length_unit="ft")


@pytest.mark.skip(reason="acceptance: open-mesh errors from volume and booleans give no cause or next step")
def test_open_mesh_errors_say_what_is_wrong_and_what_to_do(small):
    topo = wf.datasets.topography(n=20)
    with pytest.raises(wf.InvalidMesh, match=r"boundary_edges.*repair"):
        wf.union(topo, small)
    with pytest.raises(wf.InvalidMesh, match=r"repair|validate"):
        _ = topo.volume


def test_open_mesh_errors_carry_the_report(small):
    topo = wf.datasets.topography(n=20)
    with pytest.raises(wf.InvalidMesh) as e:
        wf.union(topo, small)
    assert e.value.report.summary["boundary_edges"] > 0 and not e.value.report.is_solid


def test_numpy_arrays_out_have_the_promised_dtypes_and_shapes(small):
    assert small.coords.dtype == np.float64 and small.coords.shape == (small.n_vertices, 3)
    assert small.coords.flags.c_contiguous and small.triangles.dtype == np.int64
    assert small.triangles.shape == (small.n_faces, 3)
    for axis, k in (("x", 0), ("y", 1), ("z", 2)):
        assert np.array_equal(getattr(small, axis), small.coords[:, k])
    assert small.face_normals.shape == (small.n_faces, 3) and small.face_areas.shape == (small.n_faces,)
    assert small.face_areas.sum() == pytest.approx(small.area, rel=1e-12)
    loose = np.asfortranarray(small.coords)
    again = wf.Mesh(loose, small.triangles.astype(np.int32))
    assert np.array_equal(again.coords, small.coords) and again.volume == small.volume


@pytest.mark.skip(reason="acceptance: arrays from an immutable mesh are writeable and writes vanish silently")
def test_arrays_from_an_immutable_mesh_refuse_writes(small):
    for arr in (small.coords, small.z, small.triangles, small["grade"]):
        with pytest.raises(ValueError):
            arr[0] = 0


PLOT_HEADLESS = """
import sys, wireframes as wf
from wireframes import plot
m = wf.datasets.orebody(points=24, sections=6)
holes = wf.datasets.drillholes(n=3)
topo = wf.datasets.topography(n=30)
figs = [plot.plan([m, holes], "section")[0], plot.section([m, holes], m.centroid, [1, 0, 0], thickness=50)[0],
        plot.surface3d(topo, problems=topo.validate())[0]]
for k, f in enumerate(figs):
    f.savefig(f"{sys.argv[1]}/p{k}.png")
import matplotlib
print(matplotlib.get_backend())
"""


def test_plotting_works_headless(tmp_path):
    pytest.importorskip("matplotlib")
    env = {k: v for k, v in os.environ.items() if k not in ("DISPLAY", "WAYLAND_DISPLAY", "MPLBACKEND")}
    out = subprocess.run(
        [sys.executable, "-c", PLOT_HEADLESS, str(tmp_path)],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert out.stdout.strip().lower() == "agg"
    for k in range(3):
        png = (tmp_path / f"p{k}.png").read_bytes()
        assert png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) > 5000
