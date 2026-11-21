"""Second-round acceptance scenarios for files as mining packages write them: odd but legal DXF,
Datamine pairs, OBJ, STL, OFF, PLY, VTK and VTU at realistic size, files from independent writers,
huge declared counts, attribute fidelity, and the Python API as a mining engineer uses it. All data
are seeded and synthetic."""

import re
import struct
import sys
import warnings
import xml.etree.ElementTree as ET
import zlib
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
import wireframes as wf

try:
    import resource
except ImportError:
    resource = None

TET = np.array([(0, 0, 0), (10, 0, 0), (0, 10, 0), (0, 0, 10)], float)
TET_FACES = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
TET_VOLUME = 1000.0 / 6.0
BOX = [(0, 0, 0), (2, 0, 0), (2, 3, 0), (0, 3, 0), (0, 0, 4), (2, 0, 4), (2, 3, 4), (0, 3, 4)]
QUADS = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
LOCATED = re.compile(r":\d+: |byte \d+")


def is_box(m):
    return m.is_solid and m.volume == 24.0 and m.area == 52.0


def lens(points=300, sections=170, seed=11):
    m = wf.datasets.orebody(points=points, sections=sections)
    rng = np.random.default_rng(seed)
    return wf.Mesh(m.coords, m.triangles).with_vertex_column("grade", rng.gamma(2.0, 0.6, m.n_vertices))


def quad_box(n, step):
    """A closed box of n x n quads per side with integer lattice corners scaled by `step`."""
    ids, coords, quads = {}, [], []

    def vid(p):
        if p not in ids:
            ids[p] = len(coords)
            coords.append(p)
        return ids[p]

    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        for side in (0, n):
            for u in range(n):
                for v in range(n):
                    corners = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        p = [0, 0, 0]
                        p[a], p[b], p[c] = side, u + du, v + dv
                        corners.append(vid(tuple(p)))
                    quads.append(corners if side == n else corners[::-1])
    return np.array(coords, float) * step, np.array(quads)


def g(code, value):
    return f"{code:>3}\r\n{value}\r\n"


def dxf(entities, header="", blocks=""):
    s = g(0, "SECTION") + g(2, "HEADER") + header + g(0, "ENDSEC")
    if blocks:
        s += g(0, "SECTION") + g(2, "BLOCKS") + blocks + g(0, "ENDSEC")
    return s + g(0, "SECTION") + g(2, "ENTITIES") + entities + g(0, "ENDSEC") + g(0, "EOF")


def face3d(p, tri, layer="ORE", extra="", reverse=False):
    corners = [p[tri[0]], p[tri[1]], p[tri[2]], p[tri[2]]]
    groups = [(10 + j + 10 * ax, repr(float(corners[j][ax]))) for j in range(4) for ax in range(3)]
    if reverse:
        groups = groups[::-1]
    return g(0, "3DFACE") + g(8, layer) + extra + "".join(g(c, v) for c, v in groups)


def tet_dxf(**kw):
    return "".join(face3d(TET, t, **kw) for t in TET_FACES)


def polyface(points, faces, layer="BOX", counts=None, vertex_extra=""):
    n_v, n_f = counts if counts is not None else (len(points), len(faces))
    s = g(0, "POLYLINE") + g(8, layer) + g(66, 1) + g(10, 0.0) + g(20, 0.0) + g(30, 0.0)
    s += g(70, 64) + g(71, n_v) + g(72, n_f)
    for q in points:
        s += (
            g(0, "VERTEX") + vertex_extra + g(8, layer) + g(10, q[0]) + g(20, q[1]) + g(30, q[2]) + g(70, 192)
        )
    for f in faces:
        s += g(0, "VERTEX") + vertex_extra + g(8, layer) + g(10, 0.0) + g(20, 0.0) + g(30, 0.0) + g(70, 128)
        s += "".join(g(71 + k, i + 1) for k, i in enumerate(f))
    return s + g(0, "SEQEND")


def subd_mesh(points, faces, level=0, edges=True):
    s = g(0, "MESH") + g(5, "3A") + g(100, "AcDbEntity") + g(8, "FAULT") + g(100, "AcDbSubDMesh")
    s += g(71, 2) + g(72, 0) + g(91, level) + g(92, len(points))
    s += "".join(g(10, q[0]) + g(20, q[1]) + g(30, q[2]) for q in points)
    s += g(93, 4 * len(faces)) + "".join(g(90, 3) + "".join(g(90, i) for i in f) for f in faces)
    pairs = sorted({tuple(sorted((f[k], f[(k + 1) % 3]))) for f in faces for k in range(3)}) if edges else []
    s += g(94, len(pairs)) + "".join(g(90, a) + g(90, b) for a, b in pairs)
    s += g(95, len(pairs)) + "".join(g(140, 0.0) for _ in pairs)
    return s + g(90, 0)


def write_pair(tmp_path, tag, pt, tr, encoding="utf-8"):
    (tmp_path / f"{tag}_pt.csv").write_bytes(pt.encode(encoding))
    (tmp_path / f"{tag}_tr.csv").write_bytes(tr.encode(encoding))
    return tmp_path / f"{tag}_tr.csv"


TET_PT = "XP,YP,ZP,PID\n" + "".join(
    f"{q[0]!r},{q[1]!r},{q[2]!r},{i + 1}\n" for i, q in enumerate(TET.tolist())
)
TET_TR = "TRIANGLE,PID1,PID2,PID3\n" + "".join(
    f"{k + 1},{a + 1},{b + 1},{c + 1}\n" for k, (a, b, c) in enumerate(TET_FACES)
)


# Realistic sizes: one closed solid through every polygon format, a package DXF and a Datamine pair.


@pytest.fixture(scope="module")
def quad_box_files(tmp_path_factory):
    d = tmp_path_factory.mktemp("quadbox")
    p, q = quad_box(130, 0.5)
    obj = ["# exported by a mine planning package", "mtllib box.mtl", "o box"]
    obj += ["v  " + " ".join(repr(float(c)) for c in v) for v in p]
    n_v = len(p)
    for k, f in enumerate(q):
        if k % 16900 == 0:
            obj.append(f"g side{k // 16900}")
        if k % 2:
            obj.append("f " + " ".join(f"{i - n_v}//1" for i in f))
        else:
            obj.append("f " + " ".join(f"{i + 1}/1" for i in f) + " # quad")
    (d / "box.obj").write_bytes(("\r\n".join(obj) + "\r\n").encode())
    off = ["OFF # quads", f"{len(p)} {len(q)} 0", "# vertices"] + [
        " ".join(repr(float(c)) for c in v) for v in p
    ]
    off += ["4 " + " ".join(map(str, f)) for f in q]
    (d / "box.off").write_text("\n".join(off) + "\n")
    head = (
        f"ply\r\nformat binary_big_endian 1.0\r\ncomment package\r\nelement vertex {len(p)}\r\n"
        "property float32 confidence\r\nproperty float64 x\r\nproperty float64 y\r\nproperty float64 z\r\n"
        f"element face {len(q)}\r\nproperty uint8 flag\r\nproperty list uint16 int32 vertex_indices\r\nend_header\r\n"
    )
    vert = np.zeros(len(p), [("c", ">f4"), ("x", ">f8"), ("y", ">f8"), ("z", ">f8")])
    vert["c"], vert["x"], vert["y"], vert["z"] = 0.5, p[:, 0], p[:, 1], p[:, 2]
    fac = np.zeros(len(q), [("flag", "u1"), ("n", ">u2"), ("i", ">i4", 4)])
    fac["flag"], fac["n"], fac["i"] = 1, 4, q
    (d / "box.ply").write_bytes(head.encode() + vert.tobytes() + fac.tobytes())
    vtk = ["# vtk DataFile Version 4.2", "quads", "ASCII", "DATASET POLYDATA", f"POINTS {len(p)} double"]
    vtk += [" ".join(repr(float(c)) for c in v) for v in p]
    vtk += [f"POLYGONS {len(q)} {5 * len(q)}"] + ["4 " + " ".join(map(str, f)) for f in q]
    (d / "box.vtk").write_text("\n".join(vtk) + "\n")
    return d, p, q


@pytest.mark.parametrize("name", ["box.obj", "box.off", "box.ply", "box.vtk"])
def test_quad_box_at_2e5_triangles_matches_the_closed_form(quad_box_files, name):
    d, p, q = quad_box_files
    m = wf.read_mesh(d / name)
    assert m.n_faces == 2 * len(q) > 200_000 and m.n_vertices == len(p)
    assert m.is_solid
    side = 130 * 0.5
    assert m.volume == pytest.approx(side**3, rel=1e-12)
    assert m.area == pytest.approx(6 * side**2, rel=1e-12)
    assert np.array_equal(m["polygon"], np.repeat(np.arange(len(q)), 2))


def test_quad_box_reads_bit_identically_through_every_polygon_format(quad_box_files):
    d, _, _ = quad_box_files
    meshes = [wf.read_mesh(d / n) for n in ("box.obj", "box.off", "box.ply", "box.vtk")]
    first = meshes[0]
    for m in meshes[1:]:
        assert np.array_equal(m.coords, first.coords) and np.array_equal(m.triangles, first.triangles)
        assert m.volume == first.volume
    assert first["group"][0] == "side0" and first["group"][-1] == "side5"
    assert set(first["object"].tolist()) == {"box"}
    assert wf.read_mesh(d / "box.ply")["flag"].tolist()[:2] == [1, 1]


def test_package_dxf_at_1e5_faces_on_two_layers(tmp_path):
    m = lens()
    p, t = m.coords, m.triangles
    half = len(t) // 2
    ents = []
    for k, tri in enumerate(t):
        handle = g(5, f"{0x400 + k:X}") + g(100, "AcDbEntity")
        xdata = g(1001, "PKG") + g(1000, "EOF") + g(1040, 2.7) if k % 997 == 0 else ""
        ents.append(face3d(p, tri, layer="ORE_HG" if k < half else "ORE_LG", extra=handle) + xdata)
    decoy = g(0, "BLOCK") + g(8, "0") + g(2, "DECOY") + g(70, 0) + g(10, 0.0) + g(20, 0.0) + g(30, 0.0)
    decoy += tet_dxf(layer="DECOY") + g(0, "ENDBLK") + g(8, "0")
    text = dxf("".join(ents), header=g(9, "$INSUNITS") + g(70, "     6"), blocks=decoy)
    (tmp_path / "lens.dxf").write_bytes(text.encode())
    data = wf.read(tmp_path / "lens.dxf")
    got = data.mesh
    assert got.n_faces == m.n_faces > 100_000 and data.length_unit == "m"
    assert np.array_equal(got.coords[got.triangles], p[t])
    assert got.is_solid and got.volume == m.volume
    hg = got.filter_faces(got["layer"] == "ORE_HG")
    lg = got.filter_faces(got["layer"] == "ORE_LG")
    assert hg.n_faces + lg.n_faces == got.n_faces
    assert hg.area + lg.area == pytest.approx(got.area, rel=1e-12)


def test_datamine_pair_at_1e5_faces_shuffled_with_gapped_ids(tmp_path):
    m = lens()
    rng = np.random.default_rng(29)
    pid = 1000 + 3 * np.arange(m.n_vertices)
    density = rng.normal(2.7, 0.1, m.n_vertices)
    order = rng.permutation(m.n_vertices)
    rows = ["ZP,DENSITY,PID,XP,YP,ROCK"]
    for k, i in enumerate(order):
        x, y, z = m.coords[i].tolist()
        rows.append(f"{z!r},{float(density[i])!r},{pid[i]},{x!r},{y!r},BIF")
        if k % 5000 == 0:
            rows.append("")
    tri = ["PID2,SID,PID1,TRIANGLE,PID3"] + [
        f"{pid[b]},{k % 9},{pid[a]},{k},{pid[c]}" for k, (a, b, c) in enumerate(m.triangles)
    ]
    path = write_pair(tmp_path, "lens", "\r\n".join(rows) + "\r\n\r\n", "\n".join(tri) + "\n")
    got = wf.read_mesh(path)
    assert got.n_faces == m.n_faces > 100_000
    assert np.array_equal(got.coords[got.triangles], m.coords[m.triangles])
    assert got.is_solid and got.volume == m.volume
    back = {tuple(c): v for c, v in zip(got.coords.tolist(), got["DENSITY"].tolist(), strict=True)}
    assert all(back[tuple(c)] == density[i] for i, c in enumerate(m.coords.tolist()))
    assert got["SID"].tolist() == [float(k % 9) for k in range(m.n_faces)]


CHAIN = [
    ("a.ply", {}),
    ("b.vtu", {"zlib": True}),
    ("c.vtk", {}),
    ("d.dxf", {"dxf_entity": "polyface"}),
    ("e_tr.csv", {}),
    ("f.obj", {}),
    ("g.off", {}),
    ("h.vtp", {}),
]


def test_a_lens_survives_a_chain_of_every_format_with_its_volume(tmp_path):
    m = lens()
    expected = m.coords[m.triangles]
    current = m
    for name, kw in CHAIN:
        wf.write(tmp_path / name, current, skip_unsupported_columns=True, **kw)
        current = wf.read_mesh(tmp_path / name)
        assert np.array_equal(current.coords[current.triangles], expected), name
        assert current.volume == m.volume and current.is_solid, name
    assert current.n_faces == m.n_faces


def test_reading_in_threads_matches_reading_in_sequence(tmp_path):
    m = lens(120, 60)
    names = []
    for name, kw in CHAIN:
        wf.write(tmp_path / name, m, skip_unsupported_columns=True, **kw)
        names.append(tmp_path / name)
    seq = [wf.read_mesh(p) for p in names * 3]
    with ThreadPoolExecutor(max_workers=6) as pool:
        par = list(pool.map(wf.read_mesh, names * 3))
    for a, b in zip(seq, par, strict=True):
        assert np.array_equal(a.coords, b.coords) and np.array_equal(a.triangles, b.triangles)


# DXF written with odd but legal formatting.

PF = polyface(BOX, QUADS)
M, N = 2, 24
CYL = g(0, "POLYLINE") + g(8, "CYL") + g(66, 1) + g(70, 16 + 32) + g(71, M) + g(72, N)
for _i in range(M):
    for _j in range(N):
        _a = 2 * np.pi * _j / N
        CYL += g(0, "VERTEX") + g(8, "CYL") + g(10, repr(float(np.cos(_a)))) + g(20, repr(float(np.sin(_a))))
        CYL += g(30, float(_i)) + g(70, 64)
CYL += g(0, "SEQEND")

DXF_VARIANTS = {
    "reversed_corner_groups": (dxf(tet_dxf(reverse=True)), "tet"),
    "xdata_holding_keywords": (
        dxf(
            "".join(face3d(TET, t) + g(1001, "ACAD") + g(1000, "EOF") + g(1000, "SECTION") for t in TET_FACES)
        ),
        "tet",
    ),
    "comments_between_entities": (dxf("".join(face3d(TET, t) + g(999, "note") for t in TET_FACES)), "tet"),
    "lf_line_ends": (dxf(tet_dxf()).replace("\r\n", "\n"), "tet"),
    "codes_padded_with_trailing_spaces": (dxf(tet_dxf()).replace(" 10\r\n", " 10  \r\n"), "tet"),
    "codes_padded_with_tabs": (dxf(tet_dxf()).replace(" 10\r\n", "\t10\r\n"), "tet"),
    "entities_section_only": (
        g(0, "SECTION") + g(2, "ENTITIES") + tet_dxf() + g(0, "ENDSEC") + g(0, "EOF"),
        "tet",
    ),
    "block_definitions_ignored": (
        dxf(
            tet_dxf(),
            blocks=g(0, "BLOCK") + g(2, "B") + g(70, 0) + tet_dxf(layer="BLK") + g(0, "ENDBLK"),
        ),
        "tet",
    ),
    "text_holding_eof_skipped": (dxf(tet_dxf() + g(0, "TEXT") + g(8, "T") + g(1, "EOF")), "tet"),
    "polyface_of_quads": (dxf(PF), "box"),
    "polyface_with_zero_counts": (dxf(polyface(BOX, QUADS, counts=(0, 0))), "box"),
    "polyface_vertices_with_handles_and_subclasses": (
        dxf(polyface(BOX, QUADS, vertex_extra=g(5, "AB") + g(100, "AcDbEntity") + g(100, "AcDbVertex"))),
        "box",
    ),
    "subd_mesh_with_edges_and_creases": (dxf(subd_mesh(TET, TET_FACES)), "tet"),
}


@pytest.mark.parametrize("name", list(DXF_VARIANTS))
def test_dxf_odd_but_legal_formatting(tmp_path, name):
    text, shape = DXF_VARIANTS[name]
    (tmp_path / f"{name}.dxf").write_bytes(text.encode())
    data = wf.read(tmp_path / f"{name}.dxf")
    m = data.mesh
    if shape == "tet":
        assert m.n_faces == 4 and m.is_solid and m.volume == pytest.approx(TET_VOLUME, rel=1e-15)
    else:
        assert is_box(m)


def test_dxf_polygon_mesh_closed_in_n_matches_the_prism_area(tmp_path):
    (tmp_path / "cyl.dxf").write_bytes(dxf(CYL).encode())
    m = wf.read_mesh(tmp_path / "cyl.dxf")
    assert m.n_faces == 2 * N and m.validate().summary["boundary_edges"] == 2 * N
    assert m.area == pytest.approx(2 * N * np.sin(np.pi / N), rel=1e-14)


def test_dxf_strings_from_lwpolyline_ocs_and_line(tmp_path):
    lw = g(0, "LWPOLYLINE") + g(8, "TOE") + g(90, 3) + g(70, 1) + g(38, 50.0)
    lw += (
        g(10, 1)
        + g(20, 2)
        + g(40, 0)
        + g(41, 0)
        + g(42, 0)
        + g(10, 3)
        + g(20, 2)
        + g(91, 7)
        + g(10, 3)
        + g(20, 5)
    )
    lw += g(210, 0) + g(220, 0) + g(230, -1)
    ln = g(0, "LINE") + g(8, "FAULT") + g(10, 1) + g(20, 2) + g(30, 3) + g(11, 4) + g(21, 5) + g(31, 6)
    (tmp_path / "s.dxf").write_bytes(dxf(lw + ln).encode())
    pl = wf.read_polylines(tmp_path / "s.dxf")
    assert pl.closed.tolist() == [True, False]
    assert pl.coords.tolist() == [[-1, 2, -50], [-3, 2, -50], [-3, 5, -50], [1, 2, 3], [4, 5, 6]]
    assert pl.attributes["layer"].tolist() == ["TOE", "FAULT"]


def test_dxf_layers_and_colours_round_trip_at_utm_scale(tmp_path):
    m = lens(60, 40)
    layers = np.where(np.arange(m.n_faces) % 3 == 0, "ORE", "ORE HG").tolist()
    colour = (np.arange(m.n_faces) % 7 + 1).astype(np.int64)
    m = wf.Mesh(m.coords, m.triangles).with_face_column("layer", layers).with_face_column("color", colour)
    for entity in ("3dface", "polyface"):
        wf.write(tmp_path / f"{entity}.dxf", m, dxf_entity=entity)
        back = wf.read_mesh(tmp_path / f"{entity}.dxf")
        assert back["layer"].tolist() == layers and np.array_equal(back["color"], colour)
        assert back.volume == m.volume


def test_dxf_cp1252_layer_names_before_2007_and_utf8_after(tmp_path):
    for version, enc in (("AC1009", "cp1252"), ("AC1027", "utf-8")):
        text = dxf(tet_dxf(layer="MINÉRIO"), header=g(9, "$ACADVER") + g(1, version))
        (tmp_path / f"{version}.dxf").write_bytes(text.encode(enc))
        assert wf.read_mesh(tmp_path / f"{version}.dxf")["layer"].tolist() == ["MINÉRIO"] * 4


def test_dxf_paper_space_entities_stay_out_of_the_model(tmp_path):
    paper = "".join(face3d(TET + 100, t, layer="TITLE", extra=g(67, 1)) for t in TET_FACES)
    (tmp_path / "p.dxf").write_bytes(dxf(tet_dxf() + paper).encode())
    m = wf.read_mesh(tmp_path / "p.dxf")
    assert m.n_faces == 4 and m["layer"].tolist() == ["ORE"] * 4


def test_dxf_with_a_byte_order_mark_reads(tmp_path):
    (tmp_path / "bom.dxf").write_bytes(("﻿" + dxf(tet_dxf())).encode("utf-8"))
    assert wf.read_mesh(tmp_path / "bom.dxf").volume == pytest.approx(TET_VOLUME, rel=1e-15)


def test_dxf_smoothed_mesh_warns_that_only_the_cage_is_read(tmp_path):
    (tmp_path / "s.dxf").write_bytes(dxf(subd_mesh(TET, TET_FACES, level=2)).encode())
    data = wf.read(tmp_path / "s.dxf")
    assert any("MESH" in w or "smooth" in w or "subdivision" in w for w in data.warnings)


@pytest.mark.parametrize(
    "text",
    [
        dxf(subd_mesh(TET, TET_FACES).replace(g(92, 4), g(92, 2_000_000_000))),
        dxf(CYL.replace(g(71, M) + g(72, N), g(71, 100_000) + g(72, 100_000))),
        dxf(PF.replace(g(71, 1) + g(72, 4), g(71, 1) + g(72, 99), 1)),
    ],
    ids=["mesh_vertex_count", "polygon_mesh_size", "polyface_index"],
)
def test_dxf_entity_errors_name_the_line(tmp_path, text):
    (tmp_path / "bad.dxf").write_bytes(text.encode())
    with pytest.raises(wf.FileError, match=r"bad\.dxf:\d+: "):
        wf.read(tmp_path / "bad.dxf")


# Datamine pairs as packages and spreadsheets write them.

DM_VARIANTS = {
    "byte_order_marks": ("﻿" + TET_PT, "﻿" + TET_TR),
    "spaces_around_fields": (TET_PT.replace(",", " , "), TET_TR.replace(",", " , ")),
    "lowercase_headers": (TET_PT.lower(), TET_TR.lower()),
    "crlf_and_whitespace_lines": (
        TET_PT.replace("\n", "\r\n") + "   \r\n",
        TET_TR.replace("\n", "\r\n") + " \r\n\r\n",
    ),
    "exponent_notation": (TET_PT.replace("10.0", "1.0E+01"), TET_TR),
    "triang_column": (TET_PT, TET_TR.replace("TRIANGLE", "TRIANG")),
    "no_triangle_column": (
        TET_PT,
        "PID1,PID2,PID3\n" + "".join(f"{a + 1},{b + 1},{c + 1}\n" for a, b, c in TET_FACES),
    ),
    "points_out_of_order": (
        "XP,YP,ZP,PID\n" + "".join(f"{TET[i][0]},{TET[i][1]},{TET[i][2]},{i + 1}\n" for i in (3, 1, 0, 2)),
        TET_TR,
    ),
    "large_gapped_ids": (
        "XP,YP,ZP,PID\n" + "".join(f"{q[0]},{q[1]},{q[2]},{(i + 1) * 100000}\n" for i, q in enumerate(TET)),
        "TRIANGLE,PID1,PID2,PID3\n"
        + "".join(
            f"{k + 1},{(a + 1) * 10**5},{(b + 1) * 10**5},{(c + 1) * 10**5}\n"
            for k, (a, b, c) in enumerate(TET_FACES)
        ),
    ),
    "negative_and_zero_ids": (
        "XP,YP,ZP,PID\n" + "".join(f"{q[0]},{q[1]},{q[2]},{i - 1}\n" for i, q in enumerate(TET)),
        "TRIANGLE,PID1,PID2,PID3\n"
        + "".join(f"{k},{a - 1},{b - 1},{c - 1}\n" for k, (a, b, c) in enumerate(TET_FACES)),
    ),
}


@pytest.mark.parametrize("name", list(DM_VARIANTS))
def test_datamine_pair_variants_read_the_same_tetrahedron(tmp_path, name):
    pt, tr = DM_VARIANTS[name]
    m = wf.read_mesh(write_pair(tmp_path, name, pt, tr))
    assert m.is_solid and m.volume == TET_VOLUME


def test_datamine_pair_reads_from_either_name_and_cp1252(tmp_path):
    pt = "XP,YP,ZP,PID,ROCK\n" + "".join(f"{q[0]},{q[1]},{q[2]},{i + 1},MINÉRIO\n" for i, q in enumerate(TET))
    write_pair(tmp_path, "pit", pt, TET_TR, encoding="cp1252")
    a = wf.read_mesh(tmp_path / "pit_tr.csv")
    b = wf.read_mesh(tmp_path / "pit_pt.csv")
    assert a.volume == b.volume == TET_VOLUME and a["ROCK"].tolist() == ["MINÉRIO"] * 4


def test_datamine_extra_triangle_columns_in_any_order(tmp_path):
    tr = "SID,PID2,ZONE,PID1,TRIANGLE,PID3,NAME\n" + "".join(
        f"9,{b + 1},{k},{a + 1},{k + 1},{c + 1},Z{k}\n" for k, (a, b, c) in enumerate(TET_FACES)
    )
    m = wf.read_mesh(write_pair(tmp_path, "cols", TET_PT, tr))
    assert m.volume == TET_VOLUME
    assert m["NAME"].tolist() == ["Z0", "Z1", "Z2", "Z3"] and m["ZONE"].tolist() == [0.0, 1.0, 2.0, 3.0]


def test_datamine_writer_quotes_text_and_keeps_nulls(tmp_path):
    m = lens(60, 40)
    names = [f'rock "{i % 3}", ox' if i % 5 else None for i in range(m.n_faces)]
    au = np.random.default_rng(3).gamma(2, 1, m.n_vertices)
    au[::7] = np.nan
    m = wf.Mesh(m.coords, m.triangles).with_face_column("NAME", names).with_vertex_column("AU", au)
    wf.write(tmp_path / "w_tr.csv", m)
    back = wf.read_mesh(tmp_path / "w_tr.csv")
    assert back["NAME"].tolist() == names and np.array_equal(back["AU"], au, equal_nan=True)
    assert back.volume == m.volume


@pytest.mark.parametrize(("column", "where"), [("PID", "v"), ("XP", "v"), ("PID1", "f"), ("TRIANGLE", "f")])
def test_datamine_writer_refuses_columns_named_like_its_fields(tmp_path, column, where):
    m = wf.datasets.box([0, 0, 0], [1, 1, 1])
    m = m.with_vertex_column(column, np.ones(8)) if where == "v" else m.with_face_column(column, np.ones(12))
    with pytest.raises(wf.InvalidInput, match=f'"{column}" clashes'):
        wf.write(tmp_path / "r_tr.csv", m)


def test_datamine_trailing_delimiter_adds_no_columns(tmp_path):
    m = wf.read_mesh(write_pair(tmp_path, "trail", TET_PT.replace("\n", ",\n"), TET_TR.replace("\n", ",\n")))
    assert m.volume == TET_VOLUME and m.column_names == []


def test_datamine_text_codes_round_trip_as_text(tmp_path):
    m = wf.datasets.box([0, 0, 0], [1, 1, 1]).with_face_column("CODE", ["007", "1010"] * 6)
    wf.write(tmp_path / "code_tr.csv", m)
    assert wf.read_mesh(tmp_path / "code_tr.csv")["CODE"].tolist() == ["007", "1010"] * 6


def test_datamine_large_integers_round_trip_or_are_refused(tmp_path):
    ids = np.array([2**53 + 1 + 2 * k for k in range(12)], np.int64)
    m = wf.datasets.box([0, 0, 0], [1, 1, 1]).with_face_column("BLOCKID", ids)
    try:
        wf.write(tmp_path / "ids_tr.csv", m)
    except wf.InvalidInput:
        return
    back = wf.read_mesh(tmp_path / "ids_tr.csv")["BLOCKID"]
    assert [int(v) for v in back] == ids.tolist()


# OBJ, STL, OFF and PLY as other packages write them.

OBJ_V = "".join("v {} {} {}\n".format(*v) for v in BOX)
OBJ_F = "".join("f {}\n".format(" ".join(str(i + 1) for i in q)) for q in QUADS)
OBJ_VARIANTS = {
    "faces_before_vertices": OBJ_F + OBJ_V,
    "w_and_vertex_colours": "".join("v {} {} {} 1.0 0.5 0.5 0.5\n".format(*v) for v in BOX) + OBJ_F,
    "tabs_and_trailing_blanks": OBJ_V.replace(" ", "\t") + OBJ_F.replace("\n", "  \t\n"),
    "comments_after_faces": OBJ_V + OBJ_F.replace("\n", " # face\n"),
    "line_continuation": OBJ_V + OBJ_F.replace("f 1 4", "f 1 \\\n4"),
}


@pytest.mark.parametrize("name", list(OBJ_VARIANTS))
def test_obj_variants_read_the_box(tmp_path, name):
    (tmp_path / f"{name}.obj").write_text(OBJ_VARIANTS[name])
    assert is_box(wf.read_mesh(tmp_path / f"{name}.obj"))


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (OBJ_V + "f 0 1 2\n", r"bad\.obj:9: index 0"),
        (OBJ_V + "f -1 -2 -99\n", r"bad\.obj:9: index -99"),
        ("v a 0 0\n" + OBJ_V, r"bad\.obj:1: "),
        ("v nan 0 0\n" + OBJ_V, r"bad\.obj:1: .*not finite"),
        (OBJ_V + "f 1 2 99999999999999999999\n", r"bad\.obj:9: "),
        (OBJ_V + "f 1 2 99\n", r"bad\.obj: index 99 .* 8 vertices"),
    ],
)
def test_obj_errors_name_the_line_or_the_count(tmp_path, text, message):
    (tmp_path / "bad.obj").write_text(text)
    with pytest.raises(wf.FileError, match=message):
        wf.read_mesh(tmp_path / "bad.obj")


def box_triangles():
    out = []
    for q in QUADS:
        out += [(BOX[q[0]], BOX[q[1]], BOX[q[2]]), (BOX[q[0]], BOX[q[2]], BOX[q[3]])]
    return out


def binary_stl(header=b"binary", count=None, tail=b""):
    tris = box_triangles()
    body = b"".join(struct.pack("<12fH", 0, 0, 0, *a, *b, *c, 0) for a, b, c in tris)
    return header.ljust(80, b" ") + struct.pack("<I", len(tris) if count is None else count) + body + tail


ASCII_STL = (
    "solid box\n"
    + "".join(
        "  facet normal 0 0 0\n    outer loop\n"
        + "".join(f"      vertex {x:.1e} {y} {z}\n" for x, y, z in t)
        + "    endloop\n  endfacet\n"
        for t in box_triangles()
    )
    + "endsolid box\n"
)


@pytest.mark.parametrize(
    "data",
    [
        binary_stl(),
        binary_stl(header=b"solid box\n  facet normal 0 0 0\n"),
        ASCII_STL.encode(),
        ASCII_STL.replace("\n", "\r\n").replace("  ", "\t").encode(),
        ASCII_STL.replace("endsolid box", "endsolid").encode(),
        ASCII_STL.replace("normal 0 0 0", "normal nan nan nan").encode(),
    ],
    ids=[
        "binary",
        "binary_solid_header",
        "ascii",
        "ascii_crlf_tabs",
        "ascii_bare_endsolid",
        "ascii_nan_normals",
    ],
)
def test_stl_variants_read_the_box(tmp_path, data):
    (tmp_path / "s.stl").write_bytes(data)
    assert is_box(wf.read_mesh(tmp_path / "s.stl"))


def test_binary_stl_with_padding_reads_and_warns(tmp_path):
    (tmp_path / "pad.stl").write_bytes(binary_stl(tail=b"\0" * 7))
    data = wf.read(tmp_path / "pad.stl")
    assert is_box(data.mesh) and data.warnings == [
        "binary stl declares 12 triangles; ignored 7 trailing bytes"
    ]


def test_binary_stl_with_solid_header_and_padding_is_not_read_as_empty(tmp_path):
    (tmp_path / "pad.stl").write_bytes(binary_stl(header=b"solid exported by package", tail=b"\0" * 7))
    try:
        m = wf.read_mesh(tmp_path / "pad.stl")
    except wf.WireframesError:
        return
    assert is_box(m)


def test_ascii_stl_cut_after_a_facet_is_an_error(tmp_path):
    cut = ASCII_STL[: ASCII_STL.index("endfacet", len(ASCII_STL) // 2) + len("endfacet\n")]
    (tmp_path / "cut.stl").write_text(cut)
    with pytest.raises(wf.FileError, match=r"cut\.stl"):
        wf.read_mesh(tmp_path / "cut.stl")


OFF_BOX = (
    "OFF\n8 6 0\n"
    + "".join("{} {} {}\n".format(*v) for v in BOX)
    + "".join("4 {}\n".format(" ".join(map(str, q))) for q in QUADS)
)


@pytest.mark.parametrize(
    "text",
    [
        OFF_BOX.replace("\n", "\r\n"),
        OFF_BOX.replace("OFF\n", "OFF # made by a package\n# comment\n\n"),
        "COFF\n8 6 12\n"
        + "".join("{} {} {} 255 0 0 255\n".format(*v) for v in BOX)
        + "".join("4 {} 1 0 0\n".format(" ".join(map(str, q))) for q in QUADS),
        "NOFF\n8 6 0\n"
        + "".join("{} {} {} 0 0 1\n".format(*v) for v in BOX)
        + "".join("4 {}\n".format(" ".join(map(str, q))) for q in QUADS),
    ],
    ids=["crlf", "comments", "coff", "noff"],
)
def test_off_variants_read_the_box(tmp_path, text):
    (tmp_path / "b.off").write_text(text)
    assert is_box(wf.read_mesh(tmp_path / "b.off"))


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        ("8 6 0", "4000000000 4000000000 0", r"b\.off:16: .*4000000000 vertices"),
        ("4 0 3 2 1", "4 0 3 2 99", r"b\.off:11: face index 99"),
        ("4 0 3 2 1", "-4 0 3 2 1", r"b\.off:11: corner count -4"),
        ("4 0 3 2 1", "4000000000 0 3 2 1", r"b\.off:11: "),
    ],
)
def test_off_errors_name_the_line(tmp_path, old, new, message):
    (tmp_path / "b.off").write_text(OFF_BOX.replace(old, new))
    with pytest.raises(wf.FileError, match=message):
        wf.read_mesh(tmp_path / "b.off")


def ascii_ply(eol="\n", faces_first=False):
    head = ["ply", "format ascii 1.0", "comment package", "obj_info x"]
    verts = ["element vertex 8"] + [
        f"property {t} {n}"
        for t, n in (("float", "nx"), ("double", "x"), ("double", "z"), ("double", "y"), ("uchar", "red"))
    ]
    faces = [
        "element face 6",
        "property uchar flag",
        "property list uint int vertex_indices",
        "property float qual",
    ]
    head += (faces + verts) if faces_first else (verts + faces)
    vb = "".join(f"0.5 {x} {z} {y} 7\n" for x, y, z in BOX)
    fb = "".join(f"1 4 {' '.join(map(str, q))} 0.25\n" for q in QUADS)
    return (eol.join(head + ["end_header"]) + eol + ((fb + vb) if faces_first else (vb + fb))).encode()


def binary_ply(endian, eol="\n", extra_element=False):
    fmt = "binary_little_endian" if endian == "<" else "binary_big_endian"
    head = [
        "ply",
        f"format {fmt} 1.0",
        "element vertex 8",
        "property float32 nx",
        "property float64 x",
        "property float64 y",
        "property float64 z",
        "property uint8 red",
        "property int16 s",
        "element face 6",
        "property uint8 flag",
        "property list uint16 int32 vertex_indices",
        "property float32 qual",
    ]
    if extra_element:
        head += ["element edge 2", "property int vertex1", "property list uchar int extra"]
    out = (eol.join(head + ["end_header"]) + eol).encode()
    out += b"".join(struct.pack(endian + "fdddBh", 0.5, *v, 7, -3) for v in BOX)
    out += b"".join(struct.pack(endian + "BH4if", 1, 4, *q, 0.25) for q in QUADS)
    if extra_element:
        out += struct.pack(endian + "iBii", 0, 2, 5, 6) + struct.pack(endian + "iB", 1, 0)
    return out


@pytest.mark.parametrize(
    "data",
    [
        ascii_ply(),
        ascii_ply(eol="\r\n"),
        ascii_ply(faces_first=True),
        binary_ply("<"),
        binary_ply(">"),
        binary_ply("<", eol="\r\n"),
        binary_ply(">", extra_element=True),
    ],
    ids=[
        "ascii",
        "ascii_crlf",
        "faces_first",
        "binary_le",
        "binary_be",
        "binary_crlf_header",
        "binary_edge_element",
    ],
)
def test_ply_variants_read_the_box_and_keep_properties(tmp_path, data):
    (tmp_path / "b.ply").write_bytes(data)
    m = wf.read_mesh(tmp_path / "b.ply")
    assert is_box(m)
    assert m["red"].tolist() == [7] * 8 and m["flag"].tolist() == [1] * 12
    assert np.array_equal(m["qual"], np.full(12, 0.25)) and np.array_equal(m["nx"], np.full(8, 0.5))


def test_binary_ply_errors_give_the_file_byte_offset(tmp_path):
    data = binary_ply("<")
    (tmp_path / "cut.ply").write_bytes(data[:-3])
    with pytest.raises(wf.FileError) as e:
        wf.read_mesh(tmp_path / "cut.ply")
    offset = int(re.search(r"byte (\d+)", str(e.value)).group(1))
    face_start = len(data) - 6 * struct.calcsize("<BH4if")
    assert face_start <= offset <= len(data)


# VTK legacy and VTU.

VTK_BOX = (
    "# vtk DataFile Version 2.0\nbox\nASCII\nDATASET POLYDATA\nPOINTS 8 float\n"
    + " ".join(f"{x} {y} {z}" for x, y, z in BOX)
    + "\nPOLYGONS 6 30\n"
    + "".join("4 " + " ".join(map(str, q)) + "\n" for q in QUADS)
)


@pytest.mark.parametrize(
    ("text", "columns"),
    [
        (VTK_BOX.replace("\n", "\r\n"), []),
        (VTK_BOX.replace("POLYGONS", "polygons").replace("POINTS", "points"), []),
        (
            VTK_BOX
            + "CELL_DATA 6\nFIELD FieldData 2\nzone 1 6 int\n1 2 3 4 5 6\nv 3 6 double\n"
            + "1 0 0 " * 6
            + "\n",
            ["zone", "v_0", "v_1", "v_2"],
        ),
        (
            VTK_BOX
            + "CELL_DATA 6\nSCALARS zone int 1\nLOOKUP_TABLE t\n1 2 3 4 5 6\nLOOKUP_TABLE t 2\n0 0 0 1\n1 1 1 1\n",
            ["zone"],
        ),
        (
            VTK_BOX + "POINT_DATA 8\nSCALARS my%20grade double\nLOOKUP_TABLE default\n" + "0.5 " * 8 + "\n",
            ["my grade"],
        ),
    ],
    ids=["crlf", "lowercase_keywords", "field_data", "custom_lookup_table", "escaped_name"],
)
def test_vtk_variants_read_the_box_and_data(tmp_path, text, columns):
    (tmp_path / "b.vtk").write_text(text)
    m = wf.read_mesh(tmp_path / "b.vtk")
    assert is_box(m) and m.column_names == columns + ["polygon"]
    if "zone" in columns:
        assert m["zone"].tolist() == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]


MESHIO_TYPES = {
    "tetra": 10,
    "triangle": 5,
    "hexahedron": 12,
    "wedge": 13,
    "pyramid": 14,
    "quad": 9,
    "line": 3,
}


@pytest.mark.parametrize(
    ("name", "kw"),
    [
        ("mio.vtu", {}),
        ("mio_raw.vtu", {"compression": None}),
        ("mio.vtk", {}),
        ("mio42.vtk", {"file_format": "vtk42"}),
    ],
)
def test_reads_mixed_volume_cells_written_by_meshio(tmp_path, name, kw):
    meshio = pytest.importorskip("meshio")
    rng = np.random.default_rng(7)
    points = rng.random((2000, 3)) * [500.0, 500.0, 100.0] + [700000.0, 7500000.0, 200.0]
    sizes = {"tetra": 4, "triangle": 3, "hexahedron": 8, "wedge": 6, "pyramid": 5, "quad": 4, "line": 2}
    cells = [
        (t, rng.integers(0, 2000, (n, sizes[t])))
        for t, n in (
            ("tetra", 5000),
            ("triangle", 800),
            ("hexahedron", 1000),
            ("wedge", 300),
            ("pyramid", 200),
            ("quad", 400),
            ("line", 50),
        )
    ]
    rock = [np.full(len(c), k, np.int32) for k, (_, c) in enumerate(cells)]
    grade = [rng.random(len(c)) for _, c in cells]
    meshio.write(
        tmp_path / name,
        meshio.Mesh(
            points, cells, point_data={"au": rng.random(2000)}, cell_data={"rock": rock, "grade": grade}
        ),
        **kw,
    )
    data = wf.read(tmp_path / name)
    grid = data.grid
    oracle = meshio.read(tmp_path / name)
    assert np.array_equal(grid.points, oracle.points)
    assert grid.n_cells == sum(len(c) for _, c in cells)
    assert grid.cell_types.tolist() == [MESHIO_TYPES[t] for t, c in cells for _ in range(len(c))]
    sizes_back = np.diff(grid.offsets)
    assert sizes_back.tolist() == [sizes[t] for t, c in cells for _ in range(len(c))]
    # meshio permutes wedge corners between its order and the VTK order on write and on read
    vtk_order = {"wedge": [0, 2, 1, 3, 5, 4]}
    expected = np.concatenate([c[:, vtk_order[t]].ravel() if t in vtk_order else c.ravel() for t, c in cells])
    assert np.array_equal(grid.connectivity, expected)
    assert np.array_equal(grid.cell_data["rock"], np.concatenate(rock))
    assert np.array_equal(grid.cell_data["grade"], np.concatenate(grade))
    assert np.array_equal(grid.point_data["au"], oracle.point_data["au"])
    assert grid.surface().n_faces == 800 + 2 * 400
    out = tmp_path / ("again_" + name)
    wf.write(out, data)
    again = wf.read(out).grid
    assert np.array_equal(again.connectivity, grid.connectivity) and np.array_equal(
        again.offsets, grid.offsets
    )
    assert again.cell_data.column_names == ["rock", "grade"]


def parse_vtu(path):
    """A minimal reader of appended raw VTU (UInt64 headers, optional zlib) after the VTK XML file format."""
    raw = path.read_bytes()
    tag = re.search(rb"<AppendedData[^>]*>", raw)
    start = raw.index(b"_", tag.end()) + 1
    root = ET.fromstring(raw[: tag.end()] + b"</AppendedData></VTKFile>")
    compressed = root.get("compressor") == "vtkZLibDataCompressor"
    assert root.get("header_type") == "UInt64" and root.get("byte_order") == "LittleEndian"
    kinds = {"Float64": "<f8", "Int64": "<i8", "UInt8": "u1", "Int32": "<i4", "UInt64": "<u8"}

    def array(da):
        at = start + int(da.get("offset"))
        if not compressed:
            (n,) = struct.unpack_from("<Q", raw, at)
            body = raw[at + 8 : at + 8 + n]
        else:
            (blocks,) = struct.unpack_from("<Q", raw, at)
            sizes = struct.unpack_from(f"<{blocks}Q", raw, at + 24)
            pos, body = at + 24 + 8 * blocks, b""
            for s in sizes:
                body += zlib.decompress(raw[pos : pos + s])
                pos += s
        return np.frombuffer(body, kinds[da.get("type")])

    return {da.get("Name"): array(da) for da in root.iter("DataArray")}


@pytest.mark.parametrize("zlib_on", [False, True])
def test_vtu_output_parses_with_an_independent_reader(tmp_path, zlib_on):
    m = lens(120, 60)
    wf.write(tmp_path / "m.vtu", m, zlib=zlib_on)
    arrays = parse_vtu(tmp_path / "m.vtu")
    assert np.array_equal(arrays["Points"].reshape(-1, 3), m.coords)
    assert np.array_equal(arrays["connectivity"].reshape(-1, 3), m.triangles)
    assert np.array_equal(arrays["offsets"], 3 * np.arange(1, m.n_faces + 1))
    assert set(arrays["types"].tolist()) == {5}
    assert np.array_equal(arrays["grade"], m["grade"])


@pytest.mark.parametrize("ext", ["ply", "stl"])
def test_independent_readers_read_wireframes_binary_output(tmp_path, ext):
    trimesh = pytest.importorskip("trimesh")
    m = lens(60, 40)
    local = wf.Mesh((m.coords - np.floor(m.coords.min(0))).astype(np.float32).astype(np.float64), m.triangles)
    wf.write(tmp_path / f"m.{ext}", local)
    t = trimesh.load(tmp_path / f"m.{ext}", process=False)
    assert np.array_equal(np.asarray(t.vertices)[np.asarray(t.faces)], local.coords[local.triangles])


# Non-finite coordinates, huge declared counts and attribute fidelity.

NONFINITE_LOCATED = [
    ("a.obj", {}),
    ("a.ply", {"ascii": True}),
    ("b.ply", {}),
    ("a.off", {}),
    ("a.stl", {"ascii": True}),
]
NONFINITE_UNLOCATED = [
    ("a.dxf", {}),
    ("p.dxf", {"dxf_entity": "polyface"}),
    ("a.vtk", {"ascii": True}),
    ("b.vtk", {}),
    ("a.vtu", {}),
    ("a_tr.csv", {}),
]


def write_with_nan(tmp_path, name, kw):
    m = wf.Mesh(TET + [0.5, 0.25, 0.125], TET_FACES)
    path = tmp_path / name
    wf.write(path, m, **kw)
    if name.endswith("_tr.csv"):
        path = path.with_name(name.replace("_tr", "_pt"))
    data = path.read_bytes()
    if name.startswith("b.") or name.endswith(".vtu"):
        for fmt in ("<d", ">d"):
            k = data.find(struct.pack(fmt, 0.5))
            if k >= 0:
                data = data[:k] + struct.pack(fmt, float("nan")) + data[k + 8 :]
                break
    else:
        data = data.replace(b"0.5", b"nan", 1)
    path.write_bytes(data)
    return tmp_path / name


@pytest.mark.parametrize(("name", "kw"), NONFINITE_LOCATED, ids=[n for n, _ in NONFINITE_LOCATED])
def test_nonfinite_coordinates_fail_naming_the_file_and_place(tmp_path, name, kw):
    path = write_with_nan(tmp_path, name, kw)
    with pytest.raises(wf.FileError, match=re.escape(name) + r"(:\d+: |: .*byte \d+)"):
        wf.read(path)


@pytest.mark.parametrize(("name", "kw"), NONFINITE_UNLOCATED, ids=[n for n, _ in NONFINITE_UNLOCATED])
def test_nonfinite_coordinates_fail_naming_the_file_and_place_everywhere(tmp_path, name, kw):
    path = write_with_nan(tmp_path, name, kw)
    with pytest.raises(wf.FileError) as e:
        wf.read(path)
    assert name.replace("_tr", "_pt")[:-4] in str(e.value) and LOCATED.search(str(e.value)), str(e.value)


def peak_mb():
    if resource is None:
        pytest.skip("the resource module is POSIX only")
    # ru_maxrss is in bytes on macOS and in kibibytes on Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (
        1024**2 if sys.platform == "darwin" else 1024
    )


def test_huge_binary_counts_fail_cleanly_without_allocating(tmp_path):
    m = wf.datasets.orebody(points=20, sections=6)
    m = wf.Mesh(m.coords, m.triangles).with_vertex_column("au", np.arange(m.n_vertices, dtype=float))
    before = peak_mb()
    escaped = []
    for name, kw in [
        ("z.vtu", {"zlib": True}),
        ("r.vtu", {}),
        ("r.vtp", {}),
        ("b.vtk", {}),
        ("b.ply", {}),
        ("b.stl", {}),
    ]:
        path = tmp_path / name
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wf.write(path, m, skip_unsupported_columns=True, **kw)
        src = path.read_bytes()
        start = src.index(b"_", src.index(b'encoding="raw"')) + 1 if name.endswith((".vtu", ".vtp")) else 0
        for k in range(start, min(len(src) - 8, start + 1500)):
            for value, fmt in ((2**31 - 1, "<I"), (2**32 - 1, ">I"), (2**62, "<Q"), (2**64 - 1, ">Q")):
                size = struct.calcsize(fmt)
                (tmp_path / ("f_" + name)).write_bytes(src[:k] + struct.pack(fmt, value) + src[k + size :])
                try:
                    wf.read(tmp_path / ("f_" + name))
                except wf.WireframesError:
                    pass
                except BaseException as e:  # noqa: BLE001
                    escaped.append((name, k, value, type(e).__name__))
    assert escaped == []
    assert peak_mb() - before < 1024


def test_huge_text_counts_fail_cleanly_without_allocating(tmp_path):
    m = wf.datasets.orebody(points=12, sections=4)
    m = wf.Mesh(m.coords, m.triangles).with_vertex_column("au", np.arange(m.n_vertices, dtype=float))
    integer = re.compile(rb"(?<![\w.])-?\d+(?![\w.])")
    before = peak_mb()
    escaped = []
    for name, kw in [
        ("a.obj", {}),
        ("a.ply", {"ascii": True}),
        ("a.off", {}),
        ("a.vtk", {"ascii": True}),
        ("a.dxf", {}),
        ("p.dxf", {"dxf_entity": "polyface"}),
        ("a_tr.csv", {}),
    ]:
        path = tmp_path / name
        wf.write(path, m, skip_unsupported_columns=True, **kw)
        src = path.read_bytes()
        for hit in list(integer.finditer(src))[:800]:
            for value in (
                b"2147483647",
                b"4294967295",
                b"9223372036854775807",
                b"18446744073709551615",
                b"-2147483648",
            ):
                (tmp_path / ("f_" + name)).write_bytes(src[: hit.start()] + value + src[hit.end() :])
                try:
                    wf.read(tmp_path / ("f_" + name))
                except wf.WireframesError:
                    pass
                except BaseException as e:  # noqa: BLE001
                    escaped.append((name, hit.start(), value, type(e).__name__))
    assert escaped == []
    assert peak_mb() - before < 1024


ATTRIBUTE_FORMATS = [
    ("a.ply", {"ascii": True}),
    ("b.ply", {}),
    ("a.vtk", {"ascii": True}),
    ("b.vtk", {}),
    ("a.vtu", {}),
    ("a.vtp", {}),
    ("a_tr.csv", {}),
]


@pytest.mark.parametrize(("name", "kw"), ATTRIBUTE_FORMATS, ids=[n for n, _ in ATTRIBUTE_FORMATS])
def test_extreme_attribute_values_keep_their_bits_and_nulls(tmp_path, name, kw):
    m = wf.datasets.box([0, 0, 0], [2, 3, 4])
    m = wf.Mesh(m.coords, m.triangles)
    columns = {
        "extreme": ("v", np.array([5e-324, -0.0, 1e-310, np.inf, -np.inf, -1e300, 1e300, 1.0])),
        "with_null": ("v", [1.5, None] * 4),
        "nan": ("v", np.array([np.nan, 1.0] * 4)),
        "i32_edge": ("f", np.array([2**31 - 1 - k for k in range(12)], np.int64)),
        "big_int": ("f", np.array([2**40 + k for k in range(12)], np.int64)),
    }
    for column, (where, values) in columns.items():
        m = m.with_vertex_column(column, values) if where == "v" else m.with_face_column(column, values)
    wf.write(tmp_path / name, m, **kw)
    back = wf.read_mesh(tmp_path / name)
    for column, (where, _) in columns.items():
        table, got = (
            (back.vertex_attributes, m.vertex_attributes)
            if where == "v"
            else (back.face_attributes, m.face_attributes)
        )
        assert np.array_equal(table.mask(column), got.mask(column)), column
        a, b = np.asarray(m[column], float), np.asarray(back[column], float)
        valid = got.mask(column)
        assert np.array_equal(a[valid].view(np.int64), b[valid].view(np.int64)), column


def test_strings_with_multipart_features_round_trip(tmp_path):
    ring = np.array([[0, 0, 100], [10, 0, 101], [10, 10, 102], [0, 10, 103]], float)
    line = np.array([[0, 0, 0], [5, 5, 5], [9, 9, 9]], float)
    pl = wf.Polylines(
        [ring, ring + [50, 0, 0], line],
        closed=[True, True, False],
        features=[0, 0, 1],
        attributes={"bench": [3.0, 4.0]},
        vertex_attributes={"z0": np.arange(11.0)},
    )
    for name in ("a.vtp", "a.vtk", "a_ascii.vtk"):
        wf.write(tmp_path / name, pl, ascii="ascii" in name)
        back = wf.read_polylines(tmp_path / name)
        assert np.array_equal(back.coords, pl.coords) and back.closed.tolist() == [True, True, False]
        assert back.feature.tolist() == [0, 0, 1] and back.attributes["bench"].tolist() == [3.0, 4.0]
        assert np.array_equal(back.vertex_attributes["z0"], np.arange(11.0))
    wf.write(tmp_path / "a.dxf", pl, skip_unsupported_columns=True)
    back = wf.read_polylines(tmp_path / "a.dxf")
    assert np.array_equal(back.coords, pl.coords) and back.feature.tolist() == [0, 0, 1]


# The Python API as a mining engineer uses it.


def test_plot_functions_draw_into_the_axes_given(tmp_path):
    plt = pytest.importorskip("matplotlib.pyplot")
    from wireframes import plot

    m = lens(300, 170)
    fig, ax = plt.subplots()
    assert plot.plan(m, "grade", ax=ax) == (fig, ax)
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(projection="3d")
    assert plot.surface3d(m, ax=ax3) == (fig3, ax3)
    f, _ = plot.section(m, m.centroid, [1, 0, 0], thickness=5.0)
    f.savefig(tmp_path / "s.png")
    assert (tmp_path / "s.png").stat().st_size > 0
    plt.close("all")


def test_plot_by_a_text_column_colours_by_category_or_says_why_not():
    pytest.importorskip("matplotlib")
    from wireframes import plot

    m = wf.datasets.box([0, 0, 0], [1, 1, 1]).with_face_column("layer", ["ORE", "WASTE"] * 6)
    points = wf.datasets.drillholes(n=3).to_points()
    for call in (
        lambda: plot.plan(m, "layer"),
        lambda: plot.surface3d(m, "layer"),
        lambda: plot.plan(points, "hole"),
    ):
        try:
            call()
        except wf.InvalidInput as e:
            assert "layer" in str(e) or "hole" in str(e)


def test_plot_missing_column_lists_every_column():
    pytest.importorskip("matplotlib")
    from wireframes import plot

    m = (
        wf.datasets.box([0, 0, 0], [1, 1, 1])
        .with_face_column("rock", [1.0, 2.0] * 6)
        .with_vertex_column("au", [0.5] * 8)
    )
    with pytest.raises(wf.MissingColumn, match=r'\["au", "rock"\]'):
        plot.plan(m, "gold")


@pytest.mark.parametrize(
    "triangles", [[[0.0, 1.5, 2.0]], np.array([[True, False, True]])], ids=["fraction", "bool"]
)
def test_mesh_refuses_triangle_indices_that_are_not_integers(triangles):
    with pytest.raises(wf.InvalidInput):
        wf.Mesh(np.eye(3), triangles)


def test_mesh_accepts_any_integer_dtype_and_layout():
    for tri in (
        np.array([[0, 1, 2]], np.uint64),
        np.asfortranarray(np.array([[0, 1, 2]], np.int32)),
        [(0, 1, 2)],
    ):
        m = wf.Mesh(np.asfortranarray(np.eye(3)), tri)
        assert m.triangles.dtype == np.int64 and m.triangles.tolist() == [[0, 1, 2]]
        assert m.coords.dtype == np.float64 and m.area == pytest.approx(np.sqrt(3) / 2, rel=1e-15)


def test_unknown_extension_error_lists_the_supported_formats(tmp_path):
    (tmp_path / "pit.xyz").write_text("1 2 3\n")
    with pytest.raises(wf.InvalidInput) as e:
        wf.read(tmp_path / "pit.xyz")
    for ext in ("obj", "ply", "off", "stl", "dxf", "vtk", "vtu", "_tr.csv"):
        assert ext in str(e.value)


def test_container_constructors_have_a_runnable_example():
    for cls in (wf.Mesh, wf.Polylines, wf.PointSet, wf.Table):
        assert ">>>" in (cls.__doc__ or ""), cls.__name__


def test_polylines_columns_are_reachable_like_mesh_columns():
    holes = wf.datasets.drillholes(n=3)
    assert holes["hole"].tolist() == holes.attributes["hole"].tolist()
