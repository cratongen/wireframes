"""Second round of geological and underground acceptance scenarios: lofts that pinch out at both
ends, share sections or fan around an axis, non-convex sections, stopes at strike and at third-panel
offsets, drives on stope floors and walls, declines through stope blocks, faults through section
rings and past vertical, domains whose block fractions add to one, thin lofted veins at 1e6
coordinates, drillholes through mesh vertices and along mesh edges, and composites on vertices. All
data are seeded and synthetic."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import test_acceptance_underground as ug
import wireframes as wf
from conftest import close

O = wf.datasets.ORIGIN
DEEP = ug.DEEP
MILLION = ug.MILLION
box = wf.datasets.box
vol = ug.vol
grid_box = ug.grid_box
tree_union = ug.tree_union


def solid_or_empty(m):
    return m.n_faces == 0 or m.is_solid


def ellipse(x, n, a=40.0, b=10.0, phase=0.0):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False) + phase
    return np.column_stack([np.full(n, x), a * np.cos(t), b * np.sin(t)])


def loft(rings, closed=True):
    return wf.loft(wf.Polylines(rings, closed=closed, length_unit="m"))


def bands(rings):
    return [loft(rings[k : k + 2], [len(r) > 2 for r in rings[k : k + 2]]) for k in range(len(rings) - 1)]


def rotate_z(m, degrees, about):
    a = np.radians(degrees)
    r = np.eye(4)
    r[:2, :2] = [[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]
    to, back = np.eye(4), np.eye(4)
    to[:3, 3], back[:3, 3] = -np.asarray(about), about
    return m.transform(back @ r @ to)


def plane_grid(centre, u, v, n, m, half_u, half_v):
    a, b = np.meshgrid(np.linspace(-half_u, half_u, n), np.linspace(-half_v, half_v, m))
    xyz = centre + a.ravel()[:, None] * np.asarray(u, float) + b.ravel()[:, None] * np.asarray(v, float)
    faces = []
    for j in range(m - 1):
        for i in range(n - 1):
            k = j * n + i
            faces += [[k, k + 1, k + n + 1], [k, k + n + 1, k + n]]
    return wf.Mesh(xyz, np.array(faces), length_unit="m")


def lofted_vein(sections=40, thick=0.1, dip=70.0, seed=7):
    """Vein 0.1 m thick at 1e6 coordinates, lofted from elliptical sections with 6 to 29 points."""
    rng = np.random.default_rng(seed)
    d = np.radians(dip)
    rings = []
    for k in range(sections):
        n = int(rng.integers(6, 30))
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        y, z = 30 * np.cos(t), thick / 2 * np.sin(t)
        rings.append(
            MILLION
            + np.column_stack(
                [
                    np.full(n, 5.0 * k),
                    y * np.cos(d) - z * np.sin(d),
                    y * np.sin(d) + z * np.cos(d) + 0.5 * np.sin(k / 5),
                ]
            )
        )
    return rings


def brute_inside_length(m, a, b):
    """Inside length of segment a-b from its exact crossings, classifying each piece by its midpoint."""
    length = np.linalg.norm(b - a)
    hits = m.segment_hits([a], [b])
    t = np.unique(np.concatenate([[0.0], np.array(hits["t"]), [1.0]]))
    mid = (t[:-1] + t[1:]) / 2
    return (np.diff(t) * length)[m.locate(a + mid[:, None] * (b - a)) <= 0].sum()


# ---------------------------------------------------------------- lofting


def test_vein_pinching_out_at_both_ends_adds_up_band_by_band():
    rings = [np.array([[0.0, 0, 0]])]
    rings += [ellipse(10.0 * k, 7 + 5 * k, a=4.0 * k, b=float(k)) for k in range(1, 5)]
    rings += [np.array([[50.0, -3, 0], [50.0, 3, 0]])]
    m = loft([O + r for r in rings], [False, True, True, True, True, False])
    assert m.is_solid
    parts = bands([O + r for r in rings])
    assert all(p.is_solid for p in parts)
    assert close(m.volume, sum(p.volume for p in parts))


def test_lofts_sharing_a_section_union_to_the_whole_loft():
    rng = np.random.default_rng(2)
    rings = [
        O
        + ellipse(
            20.0 * k, int(rng.integers(15, 60)), a=30 + 3 * np.sin(k), b=8.0 + k % 3, phase=rng.uniform()
        )
        for k in range(9)
    ]
    whole, west, east = loft(rings), loft(rings[:5]), loft(rings[4:])
    assert whole.is_solid and west.is_solid and east.is_solid
    assert close(west.volume + east.volume, whole.volume)
    assert (west & east).n_faces == 0
    u = west | east
    assert u.is_solid
    assert u.validate().summary["shells"] == 1
    assert close(u.volume, whole.volume)
    pts = np.random.default_rng(3).uniform(whole.bounds[0], whole.bounds[1], (20_000, 3))
    lw, lu = whole.locate(pts), u.locate(pts)
    off = (lw != 0) & (lu != 0)
    assert (lw[off] == lu[off]).all()


def test_loft_of_radial_sections_fanned_about_a_vertical_axis():
    rings = []
    for k, ang in enumerate(np.radians([0.0, 10, 20, 30, 40])):
        n = 24 + 9 * k
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        r, z = 100 + 10 * np.cos(t), 8 * np.sin(t)
        rings.append(O + np.column_stack([r * np.cos(ang), r * np.sin(ang), z]))
    m = loft(rings)
    assert m.is_solid
    assert close(m.volume, sum(b.volume for b in bands(rings)))
    assert 0.97 < m.volume / (np.pi * 10 * 8 * 2 * np.pi * 100 * 40 / 360) < 1.0


def test_thin_lofted_vein_at_a_million_metres_adds_up_band_by_band():
    rings = lofted_vein()
    m = loft(rings)
    assert m.is_solid
    assert close(m.volume, sum(b.volume for b in bands(rings)))


def test_loft_of_similar_l_shaped_sections_is_a_frustum():
    scale = 0.5
    corners = np.array([[0, 0], [10, 0], [10, 2], [2, 2], [2, 10], [0, 10]], float)
    corners -= corners.mean(axis=0)
    big = np.column_stack([np.zeros(6), corners])
    small = np.column_stack([np.full(6, 10.0), scale * corners])
    m = loft([big, small])
    assert m.is_solid
    a, b = 36.0, 36.0 * scale**2
    assert close(m.volume, 10.0 / 3 * (a + b + np.sqrt(a * b)))


def test_loft_does_not_depend_on_start_vertex_or_digitizing_direction():
    a, b = ellipse(0, 8), ellipse(25, 11, phase=0.37)
    ref = loft([a, b]).volume
    for shift in range(11):
        for step in (1, -1):
            assert loft([a, np.roll(b, shift, axis=0)[::step]]).volume == ref


# ---------------------------------------------------------------- stopes


@pytest.fixture(scope="module")
def stopes():
    s = ug.stope_block()
    return s, tree_union(s)


def test_union_with_each_of_its_stopes_is_idempotent(stopes):
    parts, u = stopes
    for p in parts:
        again, rest, common = u | p, u - p, u & p
        assert again.is_solid and rest.is_solid and common.is_solid
        assert close(again.volume, u.volume)
        assert close(rest.volume, u.volume - p.volume)
        assert close(common.volume, p.volume)
        assert (rest & p).n_faces == 0


def test_symmetric_difference_of_face_sharing_stopes_is_their_union(stopes):
    parts, _ = stopes
    a, b = parts[0], parts[1]
    x = a ^ b
    assert x.is_solid
    assert x.volume == a.volume + b.volume


def test_stopes_offset_by_a_third_of_a_panel_at_utm_coordinates():
    rng = np.random.default_rng(5)
    size = np.array([20.0, 16.0, 24.0])
    parts = []
    for lev in range(3):
        for panel in range(4):
            lo = DEEP + [panel * size[0] + lev * size[0] / 3, 0, lev * size[2]]
            parts.append(grid_box(lo, lo + size, tuple(int(v) for v in rng.integers(2, 7, 3))))
    u = tree_union(parts)
    assert u.is_solid
    assert close(u.volume, sum(p.volume for p in parts))


def test_union_of_face_sharing_stopes_rotated_to_strike():
    a = rotate_z(grid_box(DEEP, DEEP + [20, 16, 24], (8, 1, 2)), 37.0, DEEP)
    b = rotate_z(grid_box(DEEP + [0, 16, 0], DEEP + [20, 32, 24], (1, 3, 2)), 37.0, DEEP)
    assert a.is_solid and b.is_solid
    u = a | b
    assert u.is_solid
    assert close(u.volume, a.volume + b.volume)


# ---------------------------------------------------------------- drives


def test_drive_floor_on_the_stope_floor():
    profile, area = ug.drive_profile()
    path = DEEP + np.array([[-30.0, 10, 0], [80.0, 10, 0]]) - [0, 0, profile[:, 1].min()]
    d = ug.tube(ug.sweep(path, profile))
    assert d.bounds[0][2] == DEEP[2]
    stope = box(DEEP, DEEP + [60, 20, 24])
    cut, common, both = stope - d, stope & d, stope | d
    assert cut.is_solid and common.is_solid and both.is_solid
    assert close(common.volume, area * 60.0)
    assert close(cut.volume + common.volume, stope.volume)
    assert close(both.volume + common.volume, stope.volume + d.volume)


def test_drive_along_a_tessellated_stope_wall_and_floor():
    d = box(DEEP + [-20, 0, 0], DEEP + [80, 5, 4])
    stope = grid_box(DEEP, DEEP + [60, 20, 24], (7, 5, 3))
    cut, common, both = stope - d, stope & d, stope | d
    assert cut.is_solid and common.is_solid and both.is_solid
    assert common.volume == 1200.0
    assert close(cut.volume, stope.volume - 1200.0)
    assert close(both.volume, stope.volume + 800.0)


def test_decline_through_the_stope_block_matches_stope_by_stope(stopes):
    parts, u = stopes
    d, _ = ug.drive(DEEP + np.array([[-30, -5, 3.0], [20, 8, 10.0], [60, 20, 30.0], [100, 40, 50.0]]))
    cut, common = u - d, u & d
    assert cut.is_solid and common.is_solid
    assert close(cut.volume, sum(vol(p - d) for p in parts))
    assert close(common.volume, sum(vol(p & d) for p in parts))
    assert close(cut.volume + common.volume, u.volume)


def test_thin_vein_minus_a_straight_drive_is_a_solid():
    d, _ = ug.drive(np.array([[20, -40, 0.0], [60, 30, 0.0]]))
    v, *_ = ug.vein(np.array([0.0, -20, 20]), 90.0, 70.0, 150.0, 50.0, 0.1)
    cut, common = v - d, v & d
    assert cut.is_solid and common.is_solid
    assert close(cut.volume + common.volume, v.volume)


# ---------------------------------------------------------------- faults and domains


@pytest.fixture(scope="module")
def lens31():
    return ug.lens(sections=31)


def test_fault_exactly_through_a_section_ring_matches_the_plane_split(lens31):
    x = np.linspace(-300, 300, 31)[12]
    f = plane_grid(O + [x, 0, 200], [0, 1.0, 0], [0, 0, 1.0], 9, 7, 200, 100)
    below, above = lens31.clip_by_surface(f)
    west, east = lens31.split_by_plane(O + [x, 0, 0], [1.0, 0, 0])
    for m in (below, above, west, east):
        assert m.is_solid
    assert close(below.volume + above.volume, lens31.volume)
    assert close(min(below.volume, above.volume), min(west.volume, east.volume))
    assert close(max(below.volume, above.volume), max(west.volume, east.volume))


def test_orebody_split_by_a_fault_that_turns_past_vertical(lens31):
    n, m = 60, 40
    u, w = np.meshgrid(np.linspace(-400, 400, n), np.linspace(0, 1, m))
    ang = np.radians(40 + 100 * w)
    xyz = np.column_stack([u.ravel(), (150 * np.cos(ang) - 20).ravel(), (50 + 150 * np.sin(ang)).ravel()])
    faces = []
    for j in range(m - 1):
        for i in range(n - 1):
            k = j * n + i
            faces += [[k, k + 1, k + n + 1], [k, k + n + 1, k + n]]
    below, above = lens31.clip_by_surface(wf.Mesh(O + xyz, np.array(faces), length_unit="m"))
    assert below.is_solid and above.is_solid
    assert min(below.volume, above.volume) > 0.2 * lens31.volume
    assert close(below.volume + above.volume, lens31.volume)
    assert close((below | above).volume, lens31.volume)


def test_fault_that_ends_inside_the_orebody_is_refused(lens31):
    f = plane_grid(O + [0, 0, 200], [1.0, 0, 0], [0, 0.5, 1.0], 10, 6, 100, 60)
    with pytest.raises(wf.errors.InvalidInput, match="does not cut through"):
        lens31.clip_by_surface(f)


def test_lofted_orebody_clipped_by_an_inclined_planar_fault(lens31):
    f = plane_grid(O + [37.0, 0, 200], [0, 1.0, 0], [0.2, 0, 1.0], 30, 20, 300, 200)
    below, above = lens31.clip_by_surface(f)
    assert below.is_solid and above.is_solid
    assert close(below.volume + above.volume, lens31.volume)


def test_clip_of_a_stope_by_an_inclined_planar_fault_is_closed():
    stope = box(DEEP + [60, 0, 0], DEEP + [80, 16, 24])
    f = plane_grid(DEEP + [40, 16, 40], [1.0, 0.3, 0], [0.1, 0.2, 1.0], 25, 25, 200, 200)
    lo, hi = np.array(stope.bounds[0]) - 15, np.array(stope.bounds[1]) + 15
    c = f.face_centroids
    f = f.filter_faces(((c > lo) & (c < hi)).all(axis=1))
    assert f.n_faces == 17
    below, above = stope.clip_by_surface(f)
    assert below.is_closed and above.is_closed
    assert below.is_solid and above.is_solid
    assert close(below.volume + above.volume, stope.volume)


def test_lofted_orebody_split_by_random_planes(lens31):
    rng = np.random.default_rng(2)
    for _ in range(30):
        origin = O + [rng.uniform(-250, 250), rng.uniform(-30, 30), 200 + rng.uniform(-10, 10)]
        a, b = lens31.split_by_plane(origin, rng.normal(size=3))
        assert solid_or_empty(a) and solid_or_empty(b)
        assert close(vol(a) + vol(b), lens31.volume)


def test_domain_block_fractions_add_to_one(lens31):
    model = box(O + [-400, -150, 50], O + [400, 150, 350])
    host = model - lens31
    domains = [m for piece in (host, lens31) for m in piece.clip_by_surface(ug.fault()) if m.n_faces]
    assert len(domains) == 4 and all(m.is_solid for m in domains)
    centroids = O + np.mgrid[-395:400:10, -145:150:10, 55:350:10].reshape(3, -1).T
    fractions = np.array([m.block_fraction(centroids, [10.0, 10, 10]) for m in domains])
    assert np.abs(fractions.sum(axis=0) - 1).max() < 1e-12
    for m, f in zip(domains, fractions):
        assert close(1000.0 * f.sum(), m.volume)
    rng = np.random.default_rng(6)
    pts = rng.uniform(model.bounds[0], model.bounds[1], (20_000, 3))
    loc = np.array([m.locate(pts) for m in domains])
    off = (loc != 0).all(axis=0)
    assert ((loc[:, off] == -1).sum(axis=0) == 1).all()


def test_rotated_block_fractions_add_to_the_orebody_volume(lens31):
    a = np.radians(30.0)
    r = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    g = (np.mgrid[-30:30, -20:20, -6:6].reshape(3, -1).T + 0.5) * 12.5
    f = lens31.block_fraction(O + [0, 0, 200] + g @ r.T, [12.5, 12.5, 12.5], rotation=r)
    assert ((f >= 0) & (f <= 1)).all()
    assert close(12.5**3 * f.sum(), lens31.volume)


def test_stope_block_and_host_tile_the_model(stopes):
    parts, u = stopes
    model = box(DEEP + [-30, -20, -10], DEEP + [100, 60, 110])
    host = model - u
    assert host.is_solid
    assert close(host.volume + u.volume, model.volume)
    for p in parts:
        assert vol(host & p) == 0.0
    centroids = DEEP + np.mgrid[-27.5:100:5, -17.5:60:5, -7.5:110:5].reshape(3, -1).T
    total = host.block_fraction(centroids, [5.0, 5, 5])
    for p in parts:
        total += p.block_fraction(centroids, [5.0, 5, 5])
    assert np.abs(total - 1).max() < 1e-12


# ---------------------------------------------------------------- thin veins


def test_lofted_vein_against_a_tessellated_stope_block():
    v = loft(lofted_vein())
    parts = []
    for i in range(5):
        for j in range(2):
            for k in range(2):
                lo = MILLION + [-10 + 45 * i, -15 + 15 * j, -30 + 30 * k]
                parts.append(grid_box(lo, lo + [45, 15, 30], (3 + i, 2 + j, 2 + k)))
    u = tree_union(parts)
    common = [v & p for p in parts]
    assert all(solid_or_empty(c) for c in common)
    vu = v & u
    assert vu.is_solid
    assert close(sum(vol(c) for c in common), vu.volume)
    assert close(vu.volume + vol(v - u), v.volume)


def test_lofted_vein_split_by_random_planes_at_a_million_metres():
    v = loft(lofted_vein())
    rng = np.random.default_rng(1)
    for _ in range(40):
        a, b = v.split_by_plane(
            MILLION + [rng.uniform(0, 190), rng.uniform(-5, 5), rng.uniform(-5, 5)], rng.normal(size=3)
        )
        assert solid_or_empty(a) and solid_or_empty(b)
        assert close(vol(a) + vol(b), v.volume)


# ---------------------------------------------------------------- drillholes and composites


def test_every_vertex_of_a_solid_is_on_its_boundary(stopes, lens31):
    for m in (lens31, stopes[1]):
        assert (m.locate(m.coords) == 0).all()
        assert (m.signed_distance(m.coords) == 0.0).all()
        assert m.contains(m.coords).all()
        assert not m.contains(m.coords, boundary=False).any()


def test_drillholes_through_lofted_mesh_vertices(lens31):
    rng = np.random.default_rng(3)
    for k in rng.choice(lens31.n_vertices, 200, replace=False):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        a, b = lens31.coords[k] - 130 * d, lens31.coords[k] + 170 * d
        lo, hi, inside = ug.intervals(lens31, [a, b])
        assert lo[0] == 0.0 and (lo[1:] == hi[:-1]).all() and (inside[1:] != inside[:-1]).all()
        got, brute = (hi - lo)[inside].sum(), brute_inside_length(lens31, a, b)
        # A hole within an ulp of a vertex can graze a chord a few ulps of the UTM coordinates
        # long; the brute midpoint rounds out of it, so the oracle resolves only that spacing.
        assert close(got, brute) or abs(got - brute) <= 8 * np.spacing(np.abs(a).max())


def test_drillholes_along_lofted_mesh_edges_are_inside(lens31):
    tri = lens31.triangles
    for f in np.random.default_rng(4).choice(len(tri), 100, replace=False):
        p, q = lens31.coords[tri[f, 0]], lens31.coords[tri[f, 1]]
        lo, hi, inside = ug.intervals(lens31, [p, q])
        assert inside.all() and len(lo) == 1
        assert lo[0] == 0.0 and close(hi[0], np.linalg.norm(q - p))


@pytest.mark.parametrize(
    ("a", "b", "inside"),
    [
        ([5, 5, -80], [5, 5, -120], [(0.0, 20.0)]),
        ([5, 5, -60], [5, 5, -100], [(20.0, 40.0)]),
        ([0, 0, -80], [20, 10, -100], [(0.0, 30.0)]),
        ([0, 5, -90], [20, 5, -90], [(0.0, 20.0)]),
        ([0, 0, -80], [40, 0, -80], [(0.0, 20.0)]),
    ],
    ids=[
        "collar on a face",
        "toe on a face",
        "collar on a corner",
        "both ends on faces",
        "along an edge then out",
    ],
)
def test_holes_collared_or_ending_on_a_stope_boundary(a, b, inside):
    lo, hi, ins = ug.intervals(ug.STOPE, [O + a, O + b])
    assert list(zip(lo[ins].tolist(), hi[ins].tolist())) == inside


def test_repeated_survey_station_and_depths_not_starting_at_zero():
    hole = [O + [5, 5, -70], O + [5, 5, -85], O + [5, 5, -85], O + [5, 5, -110]]
    lo, hi, ins = ug.intervals(ug.STOPE, hole, [0.0, 15, 15, 40])
    assert (lo.tolist(), hi.tolist(), ins.tolist()) == ([0.0, 10, 30], [10.0, 30, 40], [False, True, False])
    lo, hi, ins = ug.intervals(ug.STOPE, hole[::3], [123.5, 163.5])
    assert (lo.tolist(), hi.tolist(), ins.tolist()) == (
        [123.5, 133.5, 153.5],
        [133.5, 153.5, 163.5],
        [False, True, False],
    )


def test_drillhole_campaign_on_a_large_orebody():
    m = ug.lens(sections=300, counts=(400, 600))
    assert m.n_faces > 2.5e5
    rng = np.random.default_rng(13)
    total = 0.0
    for _ in range(40):
        pts, depth = ug.curved_hole(rng)
        lo, hi, inside = ug.intervals(m, pts, depth)
        assert lo[0] == depth[0] and hi[-1] == depth[-1] and (lo[1:] == hi[:-1]).all()
        got = (hi - lo)[inside].sum()
        brute = sum(brute_inside_length(m, pts[k], pts[k + 1]) for k in range(len(pts) - 1))
        assert close(got, brute)
        total += got
    assert total > 500.0


# ---------------------------------------------------------------- reuse and determinism


def test_fault_blocks_round_trip_through_datamine_csv(lens31, tmp_path):
    pieces = [
        q
        for p in lens31.clip_by_surface(ug.fault())
        for q in p.clip_by_surface(ug.fault(dip=75.0, offset=-60.0))
    ]
    pieces = [q for q in pieces if q.n_faces]
    assert len(pieces) == 3
    for k, q in enumerate(pieces):
        path = tmp_path / f"block{k}_tr.csv"
        wf.write(path, q)
        back = wf.read_mesh(path)
        assert np.array_equal(back.coords, q.coords) and np.array_equal(back.triangles, q.triangles)
        assert back.is_solid and back.volume == q.volume
    again = tree_union([wf.read_mesh(tmp_path / f"block{k}_tr.csv") for k in range(3)])
    assert close(again.volume, lens31.volume)


SCENARIO = """
import hashlib, sys
import numpy as np
sys.path.insert(0, {tests!r})
import test_acceptance_underground as ug
import test_acceptance_underground_deep as t
h = hashlib.sha256()
def add(m):
    h.update(m.coords.tobytes())
    h.update(m.triangles.tobytes())
v = t.loft(t.lofted_vein())
ob = ug.lens(sections=31)
model = t.box(ug.O + [-400, -150, 50], ug.O + [400, 150, 350])
for m in (v, ob, model - ob, *ob.split_by_plane(ug.O + [10.0, 0, 200], [1.0, 0.3, 0.2])):
    add(m)
c = ug.O + np.mgrid[-395:400:20, -145:150:20, 55:350:20].reshape(3, -1).T
h.update(ob.block_fraction(c, [20.0, 20, 20]).tobytes())
rng = np.random.default_rng(1)
for _ in range(5):
    pts, depth = ug.curved_hole(rng)
    tr = ob.trace_intervals(pts, depth)
    for k in ("from", "to", "inside"):
        h.update(np.asarray(tr[k]).tobytes())
print(h.hexdigest())
"""


def test_deep_scenarios_do_not_depend_on_thread_count():
    code = SCENARIO.format(tests=str(Path(__file__).parent))
    digests = set()
    for threads in ("1", "3", "8"):
        env = dict(os.environ, RAYON_NUM_THREADS=threads)
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        )
        digests.add(out.stdout.strip())
    assert len(digests) == 1
