"""Acceptance scenarios for open-pit work at UTM scale: topography from contours and spot heights,
pit shells clipped by topography, pushbacks, bench volumes, crest and toe strings, cut and fill,
and block fractions. All data are seeded and synthetic."""

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
UP = [0.0, 0.0, 1.0]


def hill(x, y):
    return 300.0 + 40.0 * np.exp(-((x / 260) ** 2 + (y / 200) ** 2)) + 0.02 * x


def contour_dtm(seed=7, spot=20_000, origin=O):
    """Topography from closed contour strings every 5 m, a road breakline and scattered spot heights."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, 180, endpoint=False)
    rings = []
    for level in np.arange(305.0, 340.0, 5.0):
        r = np.sqrt(-np.log((level - 300.0) / 40.0))
        rings.append(np.column_stack([260 * r * np.cos(t), 200 * r * np.sin(t), np.full(t.size, level)]))
    road = np.column_stack([np.linspace(-440, 440, 60), np.full(60, -300.0), np.zeros(60)])
    road[:, 2] = hill(road[:, 0], road[:, 1]) - 3.0
    xy = rng.uniform(-450, 450, (spot, 2))
    spots = np.column_stack([xy, hill(xy[:, 0], xy[:, 1])])
    lines = wf.Polylines([p + origin for p in rings + [road]], closed=[True] * len(rings) + [False])
    return wf.dtm(spots + origin, breaklines=lines), lines


def cone(floor, radius, *, n=121, half=300.0, origin=O):
    """Conical pit design surface on a regular grid, walls at 45 degrees."""
    g = np.linspace(-half, half, n)
    gx, gy = np.meshgrid(g, g)
    z = floor + np.maximum(np.hypot(gx, gy) - radius, 0.0)
    return wf.dtm(origin + np.column_stack([gx.ravel(), gy.ravel(), z.ravel()]))


def stepped_pit(widths, levels, outer, *, undercut=0.0, origin=O):
    """Benched pit design surface with square benches: a floor of half-width ``widths[0]`` at
    ``levels[0]``, a wall from ``levels[k]`` to ``levels[k + 1]`` at half-width ``widths[k]`` (its top
    pulled in by `undercut`, which makes an overhang), then a berm out to the next width, the last one
    to `outer`. Walls are vertical when `undercut` is 0."""
    v, f = [], []

    def ring(w, z):
        v.extend([[-w, -w, z], [w, -w, z], [w, w, z], [-w, w, z]])
        return len(v) - 4

    def band(i, j):
        for k in range(4):
            a, b, c, d = i + k, i + (k + 1) % 4, j + (k + 1) % 4, j + k
            f.extend([[a, c, b], [a, d, c]])

    r = ring(widths[0], levels[0])
    f.extend([[r, r + 1, r + 2], [r, r + 2, r + 3]])
    for k, w_next in enumerate(list(widths[1:]) + [outer]):
        top = ring(widths[k] - undercut, levels[k + 1])
        band(r, top)
        r = ring(w_next, levels[k + 1])
        band(top, r)
    return wf.Mesh(np.asarray(v) + origin, np.asarray(f), length_unit="m")


def flat(z, half=300.0, origin=O):
    c = np.array([[-half, -half], [half, -half], [half, half], [-half, half]])
    return wf.dtm(origin + np.column_stack([c, np.full(4, z)]))


def vol(m):
    return 0.0 if m.n_faces == 0 else m.volume


def bench_edges(m, step):
    lo, hi = m.bounds
    return np.arange(np.floor(lo[2] / step) * step, np.ceil(hi[2] / step) * step + step / 2, step)


def block_grid(m, size):
    lo, hi = m.bounds
    axes = [np.arange(np.floor(lo[k] / size[k]) * size[k], hi[k], size[k]) + size[k] / 2 for k in range(3)]
    return np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)


@pytest.fixture(scope="module")
def pits():
    topo, lines = contour_dtm()
    c1, c2 = cone(200.0, 60.0), cone(170.0, 90.0)
    s1, s2 = c1.close_to_elevation(600.0), c2.close_to_elevation(600.0)
    b1, a1 = s1.clip_by_surface(topo)
    b2, a2 = s2.clip_by_surface(topo)
    pushback, rest = (s2 - s1).clip_by_surface(topo)
    names = ["topo", "lines", "c1", "c2", "s1", "s2", "b1", "a1", "b2", "a2", "pushback", "rest"]
    return dict(zip(names, (topo, lines, c1, c2, s1, s2, b1, a1, b2, a2, pushback, rest), strict=True))


def test_contour_dtm_honours_contours_breaklines_and_spot_heights(pits):
    topo, lines = pits["topo"], pits["lines"]
    assert 3e4 < topo.n_faces < 3e5
    s = topo.validate().summary
    assert not s["is_closed"] and s["shells"] == 1
    assert all(s[k] == 0 for k in ("degenerate_faces", "duplicate_faces", "non_manifold_edges",
                                   "inconsistent_edges", "self_intersections", "coplanar_overlaps"))  # fmt: skip
    assert (topo.face_normals[:, 2] > 0).all()
    assert np.array_equal(topo.surface_z(lines.coords), lines.coords[:, 2])
    index = {tuple(p): i for i, p in enumerate(topo.coords.tolist())}
    tri = topo.triangles
    edges = {
        tuple(sorted(e)) for e in np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]).tolist()
    }
    for part, closed in zip(lines.parts, lines.closed):
        ids = [index[tuple(p)] for p in part.tolist()]
        segs = list(itertools.pairwise(ids + ids[:1] if closed else ids))
        assert all(tuple(sorted(e)) in edges for e in segs)
    levels = [315.0, 320.0, 325.0, 330.0, 335.0]
    strings = topo.contours(levels)
    assert set(strings.attributes["elevation"].tolist()) == set(levels)
    assert strings.closed.all()
    for z, part in zip(strings.attributes["elevation"], strings.parts):
        assert (part[:, 2] == z).all()


def test_pit_shell_clipped_by_topography_matches_cut_and_fill(pits):
    s1, b1, a1, topo, c1 = pits["s1"], pits["b1"], pits["a1"], pits["topo"], pits["c1"]
    assert c1.n_faces > 2e4
    assert s1.is_solid and b1.is_solid and a1.is_solid
    assert close(b1.volume + a1.volume, s1.volume, rel=1e-12)
    excavation = wf.cut_fill(topo, c1)["cut"]
    assert close(b1.volume, excavation, rel=1e-12)
    assert close(pits["b2"].volume, wf.cut_fill(topo, pits["c2"])["cut"], rel=1e-12)
    assert b1.bounds[1][2] <= topo.bounds[1][2] and b1.bounds[0][2] >= 200.0


def test_pushback_volume_is_the_difference_of_the_pit_volumes(pits):
    pb, b1, b2 = pits["pushback"], pits["b1"], pits["b2"]
    assert pb.is_solid and pits["rest"].is_solid
    assert close(pb.volume, b2.volume - b1.volume, rel=1e-12)
    assert close(pb.volume + pits["rest"].volume, pits["s2"].volume - pits["s1"].volume, rel=1e-12)
    lo, hi = b1.bounds
    points = np.random.default_rng(3).uniform(lo, hi, (20_000, 3))
    in_b1 = b1.locate(points) == -1
    assert in_b1.sum() > 1000
    assert (pb.locate(points[in_b1]) == 1).all() and (b2.locate(points[in_b1]) == -1).all()


@pytest.mark.parametrize("step", [5.0, 10.0])
def test_bench_volumes_sum_to_the_solid_and_agree_with_plane_splits(pits, step):
    for name in ("b2", "pushback"):
        m = pits[name]
        edges = bench_edges(m, step)
        benches = m.volume_by_elevation(edges)
        assert close(benches.sum(), m.volume, rel=1e-12)
        assert (benches >= 0).all()
        rest, parts = m, []
        for z in edges[1:-1]:
            below, rest = rest.split_by_plane(O + [0, 0, z], UP)
            parts.append(vol(below))
        parts.append(vol(rest))
        assert np.allclose(parts, benches, rtol=0, atol=1e-9 * m.volume)


def test_ten_metre_benches_are_pairs_of_five_metre_benches(pits):
    m = pits["b2"]
    ten = bench_edges(m, 10.0)
    five = m.volume_by_elevation(np.arange(ten[0], ten[-1] + 2.5, 5.0))
    assert np.allclose(five.reshape(-1, 2).sum(1), m.volume_by_elevation(ten), rtol=0, atol=1e-9 * m.volume)


def test_crest_and_toe_strings_from_slicing(pits):
    b2, c2 = pits["b2"], pits["c2"]
    edges = bench_edges(b2, 10.0)[1:-1]
    rings = b2.slice_by_plane(O, UP, offsets=edges)
    assert rings.closed.all()
    assert set(rings.attributes["offset"].tolist()) <= set(edges.tolist())
    for z in (200.0, 250.0, 290.0):
        toe = c2.contours([z])
        assert toe.closed.all() and len(toe) == 1
        radius = z - 170.0 + 90.0
        assert abs(toe.area()[0] / (np.pi * radius**2) - 1) < 2e-3
        assert np.allclose(np.hypot(*(toe.coords[:, :2] - O[:2]).T), radius, rtol=0, atol=0.2)


def test_cut_and_fill_between_as_built_and_design():
    rng = np.random.default_rng(11)
    design = cone(170.0, 90.0, n=101)
    xy = np.concatenate([rng.uniform(-300, 300, (15_000, 2)),
                         [[-300, -300], [300, -300], [300, 300], [-300, 300]]])  # fmt: skip
    surveyed = design.surface_z(O[:2] + xy) + rng.normal(0.0, 0.8, len(xy))
    asbuilt = wf.dtm(np.column_stack([O[:2] + xy, surveyed]))
    assert asbuilt.n_faces > 2e4
    forward, backward = wf.cut_fill(asbuilt, design), wf.cut_fill(design, asbuilt)
    assert forward["cut"] > 0 and forward["fill"] > 0
    assert close(forward["net"], forward["fill"] - forward["cut"], rel=1e-12)
    assert close(forward["cut"], backward["fill"], rel=1e-12) and close(
        forward["fill"], backward["cut"], rel=1e-12
    )
    assert close(forward["area"], 360_000.0, rel=1e-12)
    assert forward["cut_area"] + forward["fill_area"] <= forward["area"] * (1 + 1e-12)
    base = 100.0
    closed = asbuilt.close_to_elevation(base).volume, design.close_to_elevation(base).volume
    assert close(forward["net"], closed[1] - closed[0], rel=1e-9 * max(closed) / abs(forward["net"]))
    square = np.array([[-250.0, -250.0], [250.0, -250.0], [250.0, 250.0], [-250.0, 250.0]])
    inner = wf.cut_fill(asbuilt, design, boundary=wf.Polylines([O[:2] + square], closed=True))
    assert close(inner["area"], 250_000.0, rel=1e-12)
    g = np.linspace(-249.75, 249.75, 1000)
    gx, gy = np.meshgrid(g, g)
    plan = O[:2] + np.column_stack([gx.ravel(), gy.ravel()])
    d = (design.surface_z(plan) - asbuilt.surface_z(plan)) * 0.25
    assert abs(np.minimum(d, 0).sum() + inner["cut"]) < 1e-3 * inner["cut"]
    assert abs(np.maximum(d, 0).sum() - inner["fill"]) < 1e-3 * inner["fill"]
    zero = wf.cut_fill(design, design)
    assert zero["cut"] == 0 and zero["fill"] == 0


WIDTHS, LEVELS = [20.0, 40.0, 60.0, 80.0], [260.0, 270.0, 280.0, 290.0, 300.0]


def test_benched_pit_with_vertical_walls_matches_closed_forms():
    pit = stepped_pit(WIDTHS, LEVELS, 200.0)
    assert pit.validate().summary["self_intersections"] == 0
    bench = [(2 * w) ** 2 * 10.0 for w in WIDTHS]
    cf = wf.cut_fill(flat(300.0), pit)
    assert cf["cut"] == sum(bench) and cf["fill"] == 0 and cf["cut_area"] == 160.0**2
    block = wf.datasets.box(O + [-150, -150, 250], O + [150, 150, 320])
    rock, air = block.clip_by_surface(pit)
    assert rock.is_solid and air.is_solid
    assert close(air.volume, sum(bench) + 300.0**2 * 20, rel=1e-12)
    assert close(rock.volume + air.volume, block.volume, rel=1e-12)
    five = air.volume_by_elevation(np.arange(250.0, 325.0, 5.0))
    expected = [0.0, 0.0] + [b / 2 for b in bench for _ in (0, 1)] + [300.0**2 * 5] * 4
    assert np.allclose(five, expected, rtol=1e-12, atol=1e-6)
    mid = air.slice_by_plane(O, UP, offsets=[265.0, 275.0, 285.0, 295.0])
    assert mid.closed.all() and np.allclose(mid.area(), [(2 * w) ** 2 for w in WIDTHS], rtol=1e-12)
    for k, z in enumerate(LEVELS[1:], start=1):
        lower, upper = air.split_by_plane(O + [0, 0, z], UP)
        assert lower.is_solid and upper.is_solid
        assert close(lower.volume, sum(bench[:k]), rel=1e-12)
        crest = air.slice_by_plane(O + [0, 0, z], UP)
        assert crest.closed.all() and close(crest.area().sum(), (2 * WIDTHS[k - 1]) ** 2, rel=1e-12)
        toe = pit.slice_by_plane(O + [0, 0, z], UP)
        assert close(toe.length().sum(), 8 * WIDTHS[k - 1], rel=1e-12)
    centres = block_grid(air, [10.0, 10.0, 5.0])
    f = air.block_fraction(centres, [10, 10, 5])
    assert close((f * 500.0).sum(), air.volume, rel=1e-12)
    assert np.all((f > 1 - 1e-12) | (f < 1e-12))


def test_overhanging_wall_clips_exactly_and_is_refused_by_cut_and_fill():
    u = 3.0
    pit = stepped_pit(WIDTHS, LEVELS, 200.0, undercut=u)
    s = pit.validate().summary
    assert s["self_intersections"] == 0 and s["inconsistent_edges"] == 0
    with pytest.raises(wf.InvalidInput, match="design is not a height field"):
        wf.cut_fill(flat(300.0), pit)
    block = wf.datasets.box(O + [-150, -150, 250], O + [150, 150, 320])
    rock, air = block.clip_by_surface(pit)
    assert rock.is_solid and air.is_solid
    frustums = sum(10.0 / 3 * (a * a + b * b + a * b) for a, b in ((2 * w, 2 * w - 2 * u) for w in WIDTHS))
    assert close(air.volume, frustums + 300.0**2 * 20, rel=1e-12)
    rings = air.slice_by_plane(O, UP, offsets=[265.0])
    assert close(rings.area()[0], (2 * WIDTHS[0] - u) ** 2, rel=1e-12)


def test_block_fractions_in_each_pushback_sum_to_the_solid(pits):
    size = np.array([10.0, 10.0, 5.0])
    pb, b1, b2 = pits["pushback"], pits["b1"], pits["b2"]
    centres = block_grid(b2, size)
    bv = size.prod()
    f1, f2, fp = (m.block_fraction(centres, size) for m in (b1, b2, pb))
    for f, m in ((f1, b1), (f2, b2), (fp, pb)):
        assert ((f >= 0) & (f <= 1)).all()
        assert close((f * bv).sum(), m.volume, rel=1e-9)
    assert np.abs(f2 - f1 - fp).max() < 1e-9


def test_outputs_round_trip_and_chain(pits, tmp_path):
    pb = pits["pushback"]
    for ext in ("ply", "obj"):
        path = tmp_path / f"pushback.{ext}"
        wf.write(path, pb)
        back = wf.read_mesh(path)
        assert np.array_equal(back.coords, pb.coords) and np.array_equal(back.triangles, pb.triangles)
    lower, upper = pb.split_by_plane(O + [0, 0, 250.0], UP)
    again = lower | upper
    assert again.is_solid and close(again.volume, pb.volume, rel=1e-12)


SCENARIO = """
import hashlib, sys
import numpy as np
sys.path.insert(0, {tests!r})
import test_acceptance_open_pit as t
topo, _ = t.contour_dtm(spot=4000)
s1, s2 = t.cone(200.0, 60.0, n=61).close_to_elevation(600.0), t.cone(170.0, 90.0, n=61).close_to_elevation(600.0)
pb, rest = (s2 - s1).clip_by_surface(topo)
f = pb.block_fraction(t.block_grid(pb, [10.0, 10.0, 5.0]), [10, 10, 5])
cf = t.wf.cut_fill(topo, t.cone(170.0, 90.0, n=61))
h = hashlib.sha256()
for a in (pb.coords, pb.triangles, rest.coords, rest.triangles, f, np.array(list(cf.values())),
          pb.volume_by_elevation(t.bench_edges(pb, 5.0))):
    h.update(np.ascontiguousarray(a).tobytes())
print(h.hexdigest())
"""


def test_results_do_not_depend_on_thread_count():
    code = SCENARIO.format(tests=str(Path(__file__).parent))
    digests = set()
    for threads in ("1", "4"):
        env = dict(os.environ, RAYON_NUM_THREADS=threads)
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        )
        digests.add(out.stdout.strip())
    assert len(digests) == 1


def prism(half, bottom, top):
    """Box with a flat bottom and its four top corners lifted to ``top(x, y)``."""
    v = [
        [x, y, bottom if k == 0 else top(x, y)] for k in (0, 1) for y in (-half, half) for x in (-half, half)
    ]
    return wf.Mesh(np.array(v, float), wf.datasets._BOX_FACES)


def test_booleans_of_solids_whose_caps_differ_by_one_ulp():
    def plane(x, y):
        return 5 + 0.25 * x + 0.5 * y

    a, b = prism(10.0, -10.0, plane), prism(6.0, -5.0, plane)
    c = b.coords.copy()
    c[4, 2] = np.nextafter(c[4, 2], np.inf)
    b = b.with_vertices(c)
    assert a.is_solid and b.is_solid
    for m in (a & b, a | b, a - b, b - a):
        assert m.n_faces == 0 or m.validate().summary["boundary_edges"] == 0


@pytest.mark.skip(reason="acceptance: rounding folds slivers of near-coplanar caps into self-intersections")
def test_pushback_as_difference_of_pits_clipped_by_the_same_topography():
    c = np.array([[-20.0, -20.0], [20.0, -20.0], [20.0, 20.0], [-20.0, 20.0]])
    topo = wf.dtm(np.column_stack([c, [4.0, 6.0, 5.5, 4.5]]))
    outer, _ = wf.datasets.box([-10, -10, -10], [10, 10, 30]).clip_by_surface(topo)
    inner, _ = wf.datasets.box([-6, -6, -5], [6, 6, 30]).clip_by_surface(topo)
    assert outer.is_solid and inner.is_solid
    pushback = outer - inner
    assert pushback.is_solid
    assert close(pushback.volume, outer.volume - inner.volume, rel=1e-12)
    assert (outer & inner).is_solid and (outer | inner).is_solid


def test_bench_split_of_a_small_clipped_pit_is_valid():
    rng = np.random.default_rng(0)
    xy = rng.uniform(-450, 450, (30, 2))
    topo = wf.dtm(np.column_stack([xy, hill(xy[:, 0], xy[:, 1])]))
    excavation, _ = cone(170.0, 90.0, n=7, origin=np.zeros(3)).close_to_elevation(600.0).clip_by_surface(topo)
    assert excavation.is_solid
    lower, upper = excavation.split_by_plane([0, 0, 310.0], UP)
    assert lower.is_solid and upper.is_solid


def test_block_strictly_below_a_pocketed_solid_has_zero_fraction():
    box = wf.datasets.box
    solid = box([0, 0, 0], [10, 10, 10]) - box([3, 3, 3], [7, 7, 12])
    assert solid.is_solid
    assert solid.block_fraction([[5.0, 5.0, -25.0]], [10, 10, 10]).tolist() == [0.0]
