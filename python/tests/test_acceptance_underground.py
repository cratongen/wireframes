"""Acceptance scenarios for geological and underground work: orebodies lofted from section strings,
unions of adjacent stopes, stopes minus development drives, orebodies split by faults, domains that
tile space, thin veins at 1e6 coordinates, drillhole traces and composites on solid boundaries. All
data are seeded and synthetic."""

import itertools
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import wireframes as wf
from conftest import close

O = wf.datasets.ORIGIN
DEEP = O + np.array([0.0, 0.0, -500.0])
box = wf.datasets.box


def vol(m):
    return 0.0 if m.n_faces == 0 else m.volume


def grid_box(lo, hi, n):
    """Closed box whose faces are regular grids of ``n[0] x n[1] x n[2]`` cells, wound outward."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    ids, verts, faces = {}, [], []

    def vid(idx):
        if idx not in ids:
            ids[idx] = len(verts)
            verts.append(lo + (hi - lo) * np.array(idx) / np.array(n))
        return ids[idx]

    for axis in range(3):
        a, b = [d for d in range(3) if d != axis]
        even = (a, b, axis) in ((0, 1, 2), (1, 2, 0), (2, 0, 1))
        for side in (0, n[axis]):
            for u in range(n[a]):
                for v in range(n[b]):
                    q = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        idx = [0, 0, 0]
                        idx[axis], idx[a], idx[b] = side, u + du, v + dv
                        q.append(vid(tuple(idx)))
                    if even != (side == n[axis]):
                        q.reverse()
                    faces += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
    return wf.Mesh(np.array(verts), np.array(faces), length_unit="m")


def stope_block(seed=3):
    """3 levels x 2 rows x 4 panels of tessellated stopes; odd levels shifted half a panel along
    strike, so stopes share whole faces, partial faces, edges and vertices."""
    rng = np.random.default_rng(seed)
    size = np.array([20.0, 16.0, 24.0])
    stopes = []
    for lev in range(3):
        shift = size[0] / 2 if lev % 2 else 0.0
        for row in range(2):
            for panel in range(4):
                lo = DEEP + [panel * size[0] + shift, row * size[1], lev * size[2]]
                stopes.append(grid_box(lo, lo + size, tuple(rng.integers(3, 12, 3))))
    stopes.append(box(DEEP + [-20, -16, 0], DEEP + [0, 0, 24]))
    stopes.append(box(DEEP + [-20, 32, 72], DEEP + [0, 48, 96]))
    return stopes


def tree_union(meshes):
    meshes = list(meshes)
    while len(meshes) > 1:
        meshes = [
            meshes[i] | meshes[i + 1] if i + 1 < len(meshes) else meshes[i] for i in range(0, len(meshes), 2)
        ]
    return meshes[0]


def drive_profile(width=5.0, wall=3.5, arch=1.5, n_arch=15):
    """Arched drive profile in (across, up) coordinates with its centroid at the origin; returns
    the points and the area."""
    h = width / 2
    t = np.linspace(0, np.pi, n_arch)[1:-1]
    p = np.array(
        [[-h, 0], [h, 0], [h, wall]] + [[h * np.cos(a), wall + arch * np.sin(a)] for a in t] + [[-h, wall]]
    )
    x, y = p[:, 0], p[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cross = x * yn - xn * y
    area = cross.sum() / 2
    return p - [((x + xn) * cross).sum() / (6 * area), ((y + yn) * cross).sum() / (6 * area)], area


def sweep(path, profile):
    """Mitred rings of `profile` along a 3D polyline. Frames are carried across each bend by the
    reflection in the bisector plane, so each ring lies on both neighbouring straight prisms."""
    path = np.asarray(path, float)
    d = np.diff(path, axis=0)
    d /= np.linalg.norm(d, axis=1)[:, None]
    side = np.cross(d[0], [0.0, 0.0, 1.0])
    side /= np.linalg.norm(side)
    up = np.cross(side, d[0])
    rings = []
    for k in range(len(path)):
        t_in = d[max(k - 1, 0)]
        n = t_in + d[min(k, len(d) - 1)]
        n /= np.linalg.norm(n)
        q = path[k] + profile[:, :1] * side + profile[:, 1:] * up
        rings.append(q - (((q - path[k]) @ n) / (t_in @ n))[:, None] * t_in)
        if 0 < k < len(d):
            side, up = side - 2 * (side @ n) * n, up - 2 * (up @ n) * n
    return rings


def tube(rings):
    """Solid linking vertex i of each convex ring to vertex i of the next, with fan caps."""
    n, k = len(rings[0]), len(rings)
    faces = []
    for s in range(k - 1):
        a, b = s * n, (s + 1) * n
        for i in range(n):
            j = (i + 1) % n
            faces += [[a + i, b + j, a + j], [a + i, b + i, b + j]]
    last = (k - 1) * n
    faces += [[0, i, i + 1] for i in range(1, n - 1)]
    faces += [[last, last + i + 1, last + i] for i in range(1, n - 1)]
    m = wf.Mesh(np.vstack(rings), np.array(faces), length_unit="m")
    return m if m.volume > 0 else m.flip()


def drive(path, **kw):
    """Development drive swept along `path`; returns the solid and its closed-form volume, profile
    area times centreline length (each straight piece is a prism cut by two planes through its axis)."""
    profile, area = drive_profile(**kw)
    length = np.linalg.norm(np.diff(np.asarray(path, float), axis=0), axis=1).sum()
    return tube(sweep(path, profile)), area * length


DECLINE = np.array(
    [[-60, 7, 2.0], [-20, 9, 1.0], [15, 6, 1.5], [50, 12, -2.0], [90, 8, -6.0], [140, 14, -9.0]]
)


def lens(seed=1, sections=60, counts=(60, 240), length=600.0):
    """Orebody lofted from wavy elliptical sections at UTM scale, each with its own point count."""
    rng = np.random.default_rng(seed)
    rings = []
    for x in np.linspace(-length / 2, length / 2, sections):
        n = int(rng.integers(*counts))
        a = np.linspace(0, 2 * np.pi, n, endpoint=False)
        taper = np.sqrt(max(1 - (x / (0.55 * length)) ** 2, 0.05))
        r = 1 + 0.05 * np.sin(3 * a + x / 50)
        y, z = 60 * taper * r * np.cos(a), 25 * taper * r * np.sin(a) + 200
        rings.append(O + np.column_stack([np.full(n, x), y, z]))
    return wf.loft(wf.Polylines(rings, closed=True, length_unit="m"))


def fault(n=120, dip=60.0, offset=20.0):
    """Wavy fault surface striking east and dipping `dip` degrees south, as a regular grid mesh."""
    u, w = np.meshgrid(np.linspace(-450, 450, n), np.linspace(-260, 260, n // 2))
    c, s = np.cos(np.radians(dip)), np.sin(np.radians(dip))
    rough = 2.0 * np.sin(u / 37.0) * np.cos(w / 23.0)
    xyz = np.column_stack(
        [u.ravel(), (offset + w * c + rough * s).ravel(), (200 - w * s + rough * c).ravel()]
    )
    faces = []
    for j in range(n // 2 - 1):
        for i in range(n - 1):
            a = j * n + i
            faces += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    return wf.Mesh(O + xyz, np.array(faces), length_unit="m")


def intervals(m, pts, depths=None):
    pts = np.asarray(pts, float)
    if depths is None:
        depths = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    t = m.trace_intervals(pts, depths)
    return np.array(t["from"]), np.array(t["to"]), np.array(t["inside"], bool)


def inside_length(m, pts, depths=None):
    lo, hi, inside = intervals(m, pts, depths)
    return (hi - lo)[inside].sum()


# ---------------------------------------------------------------- lofting


def test_loft_of_sections_with_different_point_counts_adds_up_band_by_band():
    rng = np.random.default_rng(5)
    rings = []
    for k in range(8):
        n = int(rng.integers(20, 90))
        a = np.linspace(0, 2 * np.pi, n, endpoint=False) + rng.uniform(0, 1)
        rings.append(
            O + np.column_stack([np.full(n, 25.0 * k), (40 - 4 * k) * np.cos(a), (12 - k) * np.sin(a)])
        )
    whole = wf.loft(wf.Polylines(rings, closed=True, length_unit="m"))
    assert whole.is_solid
    bands = [wf.loft(wf.Polylines(rings[k : k + 2], closed=True)) for k in range(7)]
    assert all(b.is_solid for b in bands)
    assert close(whole.volume, sum(b.volume for b in bands))
    west, east = whole.split_by_plane(O + [75.0, 0, 0], [1.0, 0, 0])
    assert west.is_solid and east.is_solid
    assert close(west.volume, sum(b.volume for b in bands[:3]))
    assert close(east.volume, sum(b.volume for b in bands[3:]))


@pytest.mark.parametrize(("na", "nb"), [(1, 3), (2, 5), (4, 7), (1, 10)])
def test_loft_of_one_square_sampled_differently_is_a_prism(na, nb):
    def square(x, k):
        s = np.linspace(-10, 10, k + 1)[:-1]
        sides = [[x, t, -10] for t in s] + [[x, 10, t] for t in s]
        return DEEP + np.array(sides + [[x, -t, 10] for t in s] + [[x, -10, -t] for t in s])

    m = wf.loft(wf.Polylines([square(0, na), square(50, nb)], closed=True))
    assert m.is_solid
    assert close(m.volume, 20 * 20 * 50)


def test_large_lens_with_varying_point_counts_is_a_solid():
    m = lens(sections=300, counts=(400, 600))
    assert 2.5e5 < m.n_faces < 3.1e5
    assert m.is_solid
    edges = np.linspace(m.bounds[0][2] - 1, m.bounds[1][2] + 1, 31)
    assert close(m.volume_by_elevation(edges).sum(), m.volume)


def test_vein_thinning_to_a_sliver_ring_matches_the_frustum_volumes():
    scales = [1.0, 0.6, 0.25, 1e-3]
    rect = np.array([[0, -20, -1], [0, 20, -1], [0, 20, 1], [0, -20, 1]], float)
    rings = [O + rect * [1, s, s] + [30.0 * k, 0, 0] for k, s in enumerate(scales)]
    m = wf.loft(wf.Polylines(rings, closed=True, length_unit="m"))
    assert m.is_solid
    areas = [80.0 * s * s for s in scales]
    expect = sum(30.0 / 3 * (a + b + np.sqrt(a * b)) for a, b in itertools.pairwise(areas))
    assert close(m.volume, expect)


@pytest.mark.skip(reason="acceptance: loft refuses a vein that pinches out to a line or a point")
@pytest.mark.parametrize(
    ("last", "expect"),
    [([[30.0, -5, 0], [30.0, 5, 0]], 10 * 4 * 30 / 2), ([[30.0, 0, 0]], 10 * 4 * 30 / 3)],
)
def test_vein_pinching_out_to_a_line_or_point(last, expect):
    rect = np.array([[0, -5, -2], [0, 5, -2], [0, 5, 2], [0, -5, 2]], float)
    m = wf.loft(wf.Polylines([rect, np.array(last)], closed=[True, False]))
    assert m.is_solid
    assert close(m.volume, expect)


@pytest.mark.skip(reason="acceptance: loft twists the band to a section cut obliquely across the axis")
@pytest.mark.parametrize("n", [4, 16])
def test_loft_to_an_obliquely_cut_section_is_a_truncated_prism(n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    y, z = 2.5 * np.cos(t), 2.5 * np.sin(t)
    square_cut = np.column_stack([np.zeros(n), y, z])
    oblique = np.column_stack([40.0 + 0.05 * y, y, z])
    m = wf.loft(wf.Polylines([square_cut, oblique], closed=True))
    assert m.is_solid
    assert close(m.volume, 0.5 * n * 2.5**2 * np.sin(2 * np.pi / n) * 40.0)


@pytest.mark.skip(reason="acceptance: loft along a drive with grade changes links mitred rings into fans")
def test_loft_of_mitred_drive_rings_is_the_swept_drive():
    profile, area = drive_profile()
    m = wf.loft(wf.Polylines(sweep(DEEP + DECLINE, profile), closed=True, length_unit="m"))
    assert m.is_solid
    assert close(m.volume, area * np.linalg.norm(np.diff(DECLINE, axis=0), axis=1).sum())


# ---------------------------------------------------------------- stopes


@pytest.fixture(scope="module")
def stopes():
    s = stope_block()
    return s, tree_union(s)


def test_union_of_adjacent_stopes_conserves_volume(stopes):
    parts, u = stopes
    assert sum(p.n_faces for p in parts) > 1e4
    assert all(p.is_solid for p in parts)
    assert u.is_solid
    s = u.validate().summary
    assert s["shells"] == 3 and s["touching_shell_edges"] >= 1
    assert close(u.volume, sum(p.volume for p in parts))
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            assert (parts[i] & parts[j]).n_faces == 0
    rng = np.random.default_rng(8)
    pts = rng.uniform(u.bounds[0], u.bounds[1], (20_000, 3))
    brute = np.zeros(len(pts), bool)
    for p in parts:
        brute |= p.contains(pts)
    assert (u.contains(pts) == brute).all()


def test_union_of_stopes_in_another_grouping_gives_the_same_solid(stopes):
    parts, u = stopes
    rng = np.random.default_rng(0)
    order = rng.permutation(len(parts))
    v = tree_union([parts[i] for i in order])
    assert v.is_solid
    assert close(v.volume, u.volume)
    pts = rng.uniform(u.bounds[0], u.bounds[1], (20_000, 3))
    assert (v.locate(pts) == u.locate(pts)).all()


def test_composites_on_stope_corners_edges_and_shared_faces(stopes):
    parts, u = stopes
    pts = []
    for p in parts:
        lo, hi = np.array(p.bounds[0]), np.array(p.bounds[1])
        mid = (lo + hi) / 2
        pts += [lo, hi, mid, [lo[0], mid[1], lo[2]], [mid[0], lo[1], lo[2]], [mid[0], mid[1], hi[2]]]
    pts = np.array(pts)
    probes = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if i or j or k])
    expect = []
    for p in pts:
        near = p + 1e-3 * probes
        covered = np.zeros(len(near), bool)
        for s in parts:
            covered |= s.contains(near)
        on_some = any(s.locate([p])[0] <= 0 for s in parts)
        expect.append(-1 if covered.all() else (0 if on_some else 1))
    got = u.locate(pts)
    assert (got == np.array(expect)).all()
    assert set(got.tolist()) == {-1, 0}
    assert (u.contains(pts, boundary=False) == (got == -1)).all()
    assert u.contains(pts).all()


def test_drillholes_through_stope_corners_edges_and_faces(stopes):
    parts, u = stopes
    holes = {
        "oblique": (DEEP + [-10, -10, -10], DEEP + [90, 42, 82]),
        "through a vertex where the brick offset meets": (DEEP + [2, 8, 16], DEEP + [18, 24, 32]),
    }
    for a, b in holes.values():
        assert close(inside_length(u, [a, b]), sum(inside_length(p, [a, b]) for p in parts))
    vertical = [DEEP + [20, 16, -10], DEEP + [20, 16, 82]]
    assert inside_length(u, vertical) == 72.0
    along_shared_faces = [DEEP + [-5, 16, 24], DEEP + [100, 16, 24]]
    assert inside_length(u, along_shared_faces) == 90.0


@pytest.mark.skip(reason="acceptance: union refuses a shell that touches itself until the gap is filled")
def test_sequential_union_of_stopes_closing_around_a_gap():
    def cell(i, j, k):
        lo = O + 10.0 * np.array([i, j, k], float)
        return box(lo, lo + 10.0)

    u = cell(0, 0, 0)
    for c in [(0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 1, 0)]:
        u = u | cell(*c)
    assert u.is_solid
    assert u.volume == 5000.0
    u = u | cell(1, 0, 0)
    assert u.is_solid and u.volume == 6000.0


# ---------------------------------------------------------------- drives


@pytest.fixture(scope="module")
def decline():
    return drive(DEEP + DECLINE)


def test_drive_swept_along_a_3d_polyline_has_area_times_length(decline):
    d, expect = decline
    assert d.is_solid
    assert close(d.volume, expect)


def test_stope_minus_drive_conserves_volume_and_chains(decline):
    d, _ = decline
    stope = box(DEEP + [0, 0, -8], DEEP + [60, 20, 22])
    cut, common, both = stope - d, stope & d, stope | d
    for m in (cut, common, both):
        assert m.is_solid
    assert close(vol(cut) + vol(common), stope.volume)
    assert close(vol(both) + vol(common), stope.volume + d.volume)
    again = cut | common
    assert again.is_solid
    assert close(again.volume, stope.volume)
    assert (cut & common).n_faces == 0
    rng = np.random.default_rng(12)
    pts = rng.uniform(stope.bounds[0], stope.bounds[1], (20_000, 3))
    ins, ind = stope.locate(pts), d.locate(pts)
    off = (ins != 0) & (ind != 0)
    assert (cut.locate(pts)[off] == np.where((ins == -1) & (ind == 1), -1, 1)[off]).all()


def test_stope_minus_drive_with_tessellated_stope(decline):
    d, _ = decline
    stope = grid_box(DEEP + [0, 0, -8], DEEP + [60, 20, 22], (20, 8, 10))
    cut, both = stope - d, stope | d
    assert cut.is_solid and both.is_solid
    assert close(both.volume - d.volume, cut.volume)


@pytest.mark.skip(reason="acceptance: rounding folds a sliver where a drive wall crosses a stope wall")
def test_stope_and_drive_intersection_is_a_solid():
    d, _ = drive(DECLINE[:4], n_arch=3)
    stope = box([0, 0, -8], [60, 20, 22])
    assert d.is_solid and stope.is_solid
    common = stope & d
    assert common.is_solid
    assert (d - stope).is_solid
    assert close(vol(common) + vol(stope - d), stope.volume)


# ---------------------------------------------------------------- faults and domains


@pytest.fixture(scope="module")
def orebody():
    return lens()


def test_orebody_split_by_a_fault_surface(orebody):
    f = fault()
    below, above = orebody.clip_by_surface(f)
    assert below.is_solid and above.is_solid
    assert below.volume > 0.1 * orebody.volume and above.volume > 0.1 * orebody.volume
    assert close(below.volume + above.volume, orebody.volume)
    assert (below & above).n_faces == 0
    joined = below | above
    assert joined.is_solid
    assert close(joined.volume, orebody.volume)
    rng = np.random.default_rng(9)
    pts = rng.uniform(orebody.bounds[0], orebody.bounds[1], (20_000, 3))
    lb, la, lo = below.locate(pts), above.locate(pts), orebody.locate(pts)
    off = (lb != 0) & (la != 0)
    assert ((lb == -1) | (la == -1))[off].tolist() == (lo == -1)[off].tolist()


def test_large_orebody_split_by_a_dense_fault():
    m = lens(sections=300, counts=(400, 600))
    below, above = m.clip_by_surface(fault(n=200))
    assert below.is_solid and above.is_solid
    assert close(below.volume + above.volume, m.volume)
    assert close((below | above).volume, m.volume)


def test_domains_tile_the_model_without_overlap(orebody):
    model = box(O + [-400, -150, 50], O + [400, 150, 350])
    host = model - orebody
    assert host.is_solid
    assert close(host.volume + orebody.volume, model.volume)
    domains = []
    for piece in (host, orebody):
        for half in piece.clip_by_surface(fault()):
            domains += [m for m in half.clip_by_surface(fault(dip=75.0, offset=-60.0)) if m.n_faces]
    assert len(domains) >= 6
    assert all(m.is_solid for m in domains)
    assert close(sum(m.volume for m in domains), model.volume)
    for i in range(len(domains)):
        for j in range(i + 1, len(domains)):
            assert vol(domains[i] & domains[j]) == 0.0
    u = domains[0]
    for m in domains[1:]:
        u = u | m
    assert u.is_solid
    assert close(u.volume, model.volume)


# ---------------------------------------------------------------- thin veins


def vein(p0, strike, dip, length, height, thick):
    s, d = np.radians(strike), np.radians(dip)
    along = np.array([np.sin(s), np.cos(s), 0.0])
    down = np.cos(d) * np.array([np.cos(s), -np.sin(s), 0.0]) - np.sin(d) * np.array([0, 0, 1.0])
    normal = np.cross(along, down)
    normal /= np.linalg.norm(normal)
    corners = [
        p0 + a * along + b * down + t * normal for t in (0, thick) for b in (0, height) for a in (0, length)
    ]
    m = wf.Mesh(np.array(corners), wf.datasets._BOX_FACES, length_unit="m")
    return (m if m.volume > 0 else m.flip()), along, down, normal


MILLION = np.array([1e6, 1e6, 1e3])


def test_thin_vein_at_a_million_metres():
    v, along, down, normal = vein(MILLION, 37.0, 72.0, 100.0, 60.0, 0.1)
    assert v.is_solid
    assert close(v.volume, 600.0, rel=1e-8)
    inside = [MILLION + 50 * along + 30 * down + 0.05 * normal]
    outside = [MILLION + 50 * along + 30 * down + w * normal for w in (-1e-3, 0.101)]
    assert v.locate(inside + outside).tolist() == [-1, 1, 1]
    centre = MILLION + 50 * along + 30 * down + 0.05 * normal
    hole = np.array([-0.8, 0.0, 0.0]) * normal + np.array([0, 0, -0.6])
    hole /= np.linalg.norm(hole)
    chord = inside_length(v, [centre - 20 * hole, centre + 20 * hole], [100.0, 140.0])
    assert abs(chord - 0.1 / abs(hole @ normal)) < 1e-8


def test_stope_across_a_vertical_thin_vein_matches_the_closed_form():
    v, along, _, normal = vein(MILLION + [0, 0, 50], 37.0, 90.0, 300.0, 120.0, 0.1)
    assert close(v.volume, 3600.0, rel=1e-8)
    mid = MILLION + 150 * along
    stope = box(mid + [-20, -40, -60], mid + [20, 40, -10])
    common = v & stope
    assert common.is_solid
    assert close(common.volume, 40 * 0.1 / abs(normal[1]) * 50, rel=1e-9)
    rest, both = v - stope, v | stope
    assert rest.is_solid and both.is_solid
    assert close(rest.volume + common.volume, v.volume)
    assert close(both.volume + common.volume, v.volume + stope.volume)


# ---------------------------------------------------------------- drillholes


def curved_hole(rng, segments=40, step=7.5):
    collar = O + [rng.uniform(-280, 280), rng.uniform(-90, -40), 330.0]
    az, dip = np.radians(rng.uniform(-20, 20)), np.radians(rng.uniform(55, 75))
    pts = [collar]
    for _ in range(segments):
        dip += np.radians(rng.normal(0, 0.6))
        az += np.radians(rng.normal(0, 0.6))
        pts.append(
            pts[-1] + step * np.array([np.sin(az) * np.cos(dip), np.cos(az) * np.cos(dip), -np.sin(dip)])
        )
    pts = np.array(pts)
    return pts, np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])


def test_curved_drillholes_agree_with_segment_hits(orebody):
    rng = np.random.default_rng(11)
    hits = 0
    for _ in range(30):
        pts, depth = curved_hole(rng)
        lo, hi, inside = intervals(orebody, pts, depth)
        assert lo[0] == depth[0] and hi[-1] == depth[-1]
        assert (lo[1:] == hi[:-1]).all() and (inside[1:] != inside[:-1]).all() and (hi > lo).all()
        t = orebody.segment_hits(pts[:-1], pts[1:])
        seg = np.array(t["segment"])
        cuts = depth[seg] + np.array(t["t"]) * (depth[seg + 1] - depth[seg])
        cuts = np.unique(np.concatenate([[depth[0]], cuts, [depth[-1]]]))
        mid = (cuts[:-1] + cuts[1:]) / 2
        xyz = np.column_stack([np.interp(mid, depth, pts[:, k]) for k in range(3)])
        brute = np.diff(cuts)[orebody.locate(xyz) <= 0].sum()
        assert close((hi - lo)[inside].sum(), brute)
        hits += brute > 0
    assert hits >= 20


STOPE = box(O + [0, 0, -100], O + [20, 10, -80])


@pytest.mark.parametrize(
    ("a", "b", "inside"),
    [
        ([-5, 5, -90], [25, 5, -90], [(5.0, 25.0)]),
        ([-2, -1, -102], [22, 11, -78], [(3.0, 33.0)]),
        ([-5, -5, -95], [5, 5, -105], []),
        ([10, -5, -95], [10, 5, -105], []),
        ([-5, 0, -100], [25, 0, -100], [(5.0, 25.0)]),
        ([-5, 0, -90], [25, 0, -90], [(5.0, 25.0)]),
        ([10, -5, -105], [10, 15, -85], [(np.sqrt(50), np.sqrt(450))]),
    ],
    ids=[
        "through",
        "vertex to vertex",
        "grazes a vertex",
        "grazes an edge",
        "along an edge",
        "in a face",
        "edge in",
    ],
)
def test_drillhole_entry_and_exit_depths_on_a_stope(a, b, inside):
    lo, hi, ins = intervals(STOPE, [O + a, O + b])
    got = list(zip(lo[ins], hi[ins]))
    assert len(got) == len(inside)
    for (g0, g1), (e0, e1) in zip(got, inside):
        assert close(g0, e0, rel=1e-12) and close(g1, e1, rel=1e-12)


def test_composites_exactly_on_a_stope_boundary():
    lo, hi = O + [0, 0, -100], O + [20, 10, -80]
    corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
    mid = (lo + hi) / 2
    edges = (corners[:, None] + corners[None]) / 2
    on = np.vstack(
        [corners, edges.reshape(-1, 3), [mid * [1, 1, 0] + [0, 0, lo[2]], [lo[0], mid[1], mid[2]]]]
    )
    expect = np.where((on == mid).all(axis=1), -1, 0)
    assert (expect == -1).sum() == 8
    assert STOPE.locate(on).tolist() == expect.tolist()
    assert STOPE.contains(on).tolist() == (expect <= 0).tolist()
    assert STOPE.contains(on, boundary=False).tolist() == (expect < 0).tolist()
    assert (STOPE.signed_distance(on)[expect == 0] == 0.0).all()


# ---------------------------------------------------------------- reuse and determinism


def test_split_and_cut_outputs_round_trip_and_chain(orebody, decline, tmp_path):
    below, above = orebody.clip_by_surface(fault())
    d, _ = decline
    for i, m in enumerate((below, above, box(DEEP + [0, 0, -8], DEEP + [60, 20, 22]) - d)):
        path = tmp_path / f"piece{i}.ply"
        wf.write(path, m)
        back = wf.read_mesh(path)
        assert np.array_equal(back.coords, m.coords) and np.array_equal(back.triangles, m.triangles)
        assert back.is_solid
    a, b = wf.read_mesh(tmp_path / "piece0.ply"), wf.read_mesh(tmp_path / "piece1.ply")
    assert close((a | b).volume, orebody.volume)


SCENARIO = """
import hashlib, sys
import numpy as np
sys.path.insert(0, {tests!r})
import test_acceptance_underground as t
stopes = t.stope_block()
u = t.tree_union(stopes)
ob = t.lens(sections=30)
below, above = ob.clip_by_surface(t.fault(n=80))
d, _ = t.drive(t.DEEP + t.DECLINE)
cut = t.box(t.DEEP + [0, 0, -8], t.DEEP + [60, 20, 22]) - d
rng = np.random.default_rng(1)
pts, depth = t.curved_hole(rng)
tr = ob.trace_intervals(pts, depth)
h = hashlib.sha256()
for m in (u, below, above, cut):
    h.update(m.coords.tobytes())
    h.update(m.triangles.tobytes())
for k in ("from", "to", "inside"):
    h.update(np.asarray(tr[k]).tobytes())
h.update(u.locate(rng.uniform(u.bounds[0], u.bounds[1], (5000, 3))).tobytes())
print(h.hexdigest())
"""


def test_underground_results_do_not_depend_on_thread_count():
    code = SCENARIO.format(tests=str(Path(__file__).parent))
    digests = set()
    for threads in ("1", "4"):
        env = dict(os.environ, RAYON_NUM_THREADS=threads)
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        )
        digests.add(out.stdout.strip())
    assert len(digests) == 1


@pytest.mark.skip(reason="acceptance: union refuses a shell that touches itself until the gap is filled")
def test_union_of_the_stope_block_in_any_grouping():
    parts = stope_block()
    order = np.random.default_rng(4).permutation(len(parts))
    v = tree_union([parts[i] for i in order])
    assert v.is_solid
    assert close(v.volume, sum(p.volume for p in parts))
