"""Second round of open-pit acceptance scenarios at UTM scale: benched polygon pits against frustum
closed forms, pushbacks that share a wall with the previous stage, bench and block reconciliation per
pushback, cut and fill over a partition of the plan, topography with a vertical cliff and an
overhang, translation exactness, large meshes and determinism. All data are seeded and synthetic."""

import itertools
import os
import resource
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import wireframes as wf
from conftest import close
from test_acceptance_open_pit import bench_edges, block_grid, cone, contour_dtm, flat

O = wf.datasets.ORIGIN
UP = [0.0, 0.0, 1.0]
DATA = Path(__file__).parent / "data"
SHAPE = 1e-10

SQUARE = np.array([[1.0, 0.0], [6.123233995736766e-17, 1.0], [-1.0, 1.2246467991473532e-16],
                   [-1.8369701987210297e-16, -1.0]])  # fmt: skip
OCTAGON = np.array([[1.0, 0.0], [0.7071067811865476, 0.7071067811865475], [6.123233995736766e-17, 1.0],
                    [-0.7071067811865475, 0.7071067811865476], [-1.0, 1.2246467991473532e-16],
                    [-0.7071067811865477, -0.7071067811865475], [-1.8369701987210297e-16, -1.0],
                    [0.7071067811865474, -0.7071067811865477]])  # fmt: skip


def ring_dirs(n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([np.cos(t), np.sin(t)])


def benched_pit(dirs, benches, *, floor=200.0, toe=40.0, last_berm=30.0, east=0.0, origin=O):
    """Pit design surface with regular polygon rings along `dirs`: a floor of circumradius `toe` at
    `floor`, then `benches` batters 10 m high and 5 m wide in plan, each followed by a 6 m berm (the
    last one `last_berm` wide). Every ring is stretched by ``1 + east * max(cos, 0)``, so a later stage
    moves only the east wall and shares the west wall with the previous stage exactly. Returns the
    mesh and the ``(radius, z)`` of each ring."""
    n = len(dirs)
    k = 1.0 + east * np.maximum(dirs[:, 0], 0.0)
    rings, r, z = [(toe, floor)], toe, floor
    for b in range(benches):
        r, z = r + 5.0, z + 10.0
        rings.append((r, z))
        r += 6.0 if b < benches - 1 else last_berm
        rings.append((r, z))
    v = [[0.0, 0.0, floor]] + [
        [r * k[i] * dirs[i, 0], r * k[i] * dirs[i, 1], z] for r, z in rings for i in range(n)
    ]
    f = [[0, 1 + i, 1 + (i + 1) % n] for i in range(n)]
    for j in range(len(rings) - 1):
        a, b = 1 + j * n, 1 + (j + 1) * n
        for i in range(n):
            i1 = (i + 1) % n
            f.extend([[a + i, b + i, b + i1], [a + i, b + i1, a + i1]])
    return wf.Mesh(np.asarray(v) + origin, np.asarray(f), length_unit="m"), rings


def polygon_area(n, r):
    """Area of a regular polygon of circumradius `r`. Closed forms are compared at `SHAPE`, not 1e-12:
    ring corners at x = 7e5 round to the 1.2e-10 m spacing of f64 there, which moves each wall."""
    return n / 2 * r * r * np.sin(2 * np.pi / n)


def pit_volume(n, rings, lo, hi):
    """Closed-form volume of an unstretched benched pit between elevations `lo` and `hi`: a sum of
    frustums of similar regular polygons."""
    total = 0.0
    for (ra, za), (rb, zb) in itertools.pairwise(rings):
        a, b = max(lo, za), min(hi, zb)
        if zb > za and b > a:
            r0, r1 = (ra + (rb - ra) * (z - za) / (zb - za) for z in (a, b))
            p, q = polygon_area(n, r0), polygon_area(n, r1)
            total += (b - a) / 3 * (p + q + np.sqrt(p * q))
    return total


def is_solid(m):
    return m.n_faces == 0 or m.validate().summary["is_solid"]


def test_benched_polygon_pit_matches_frustum_closed_forms():
    n = 256
    pit, rings = benched_pit(ring_dirs(n), 10)
    assert pit.n_faces > 1e4
    solid = pit.close_to_elevation(600.0)
    ground = flat(300.0, half=2000.0)
    ex, air = solid.clip_by_surface(ground)
    assert ex.is_solid and air.is_solid
    assert close(ex.volume, pit_volume(n, rings, 200.0, 300.0), rel=SHAPE)
    assert close(ex.volume + air.volume, solid.volume, rel=1e-12)
    for step in (5.0, 10.0):
        edges = np.arange(200.0, 300.0 + step / 2, step)
        benches = ex.volume_by_elevation(edges)
        expected = [pit_volume(n, rings, a, b) for a, b in itertools.pairwise(edges)]
        assert np.allclose(benches, expected, rtol=SHAPE, atol=0)
        assert close(benches.sum(), ex.volume, rel=1e-12)
    for z in np.arange(210.0, 300.0, 10.0):
        lower, upper = ex.split_by_plane(O + [0, 0, z], UP)
        assert lower.is_solid and upper.is_solid
        assert close(lower.volume, pit_volume(n, rings, 200.0, z), rel=SHAPE)
        assert close(lower.volume + upper.volume, ex.volume, rel=1e-12)
    crest_levels = np.arange(210.0, 300.0, 10.0)
    crests = ex.slice_by_plane(O, UP, offsets=crest_levels)
    assert crests.closed.all() and len(crests) == len(crest_levels)
    design = {tuple(p) for p in pit.coords.tolist()}
    for part in crests.parts:
        assert all(tuple(p) in design for p in part.tolist())
    radii = [r for r, z in rings[1::2]][: len(crest_levels)]
    assert np.allclose(crests.area(), [polygon_area(n, r) for r in radii], rtol=SHAPE, atol=0)
    mids = ex.slice_by_plane(O, UP, offsets=crest_levels - 5.0)
    assert np.allclose(mids.area(), [polygon_area(n, r - 2.5) for r in radii], rtol=SHAPE, atol=0)
    toes = pit.slice_by_plane(O, UP, offsets=crest_levels)
    assert toes.closed.all()
    assert np.allclose(toes.length(), [2 * n * r * np.sin(np.pi / n) for r in radii], rtol=SHAPE, atol=0)


@pytest.fixture(scope="module")
def stages():
    """Three pit stages whose pushbacks move only the east wall, under flat ground at 305 m."""
    ground = flat(305.0, half=2000.0)
    solids = [benched_pit(ring_dirs(96), 12, east=e)[0].close_to_elevation(600.0) for e in (0.0, 0.3, 0.6)]
    cuts = [solids[0]] + [solids[k] - solids[k - 1] for k in (1, 2)]
    pushbacks = [c.clip_by_surface(ground)[0] for c in cuts]
    pits = [s.clip_by_surface(ground)[0] for s in solids]
    return {"solids": solids, "cuts": cuts, "pushbacks": pushbacks, "pits": pits}


def test_pushbacks_sharing_a_wall_add_up_to_each_stage(stages):
    solids, cuts, pbs, pits = stages["solids"], stages["cuts"], stages["pushbacks"], stages["pits"]
    for m in cuts + pbs + pits:
        assert m.is_solid
    for k in (1, 2):
        assert close(cuts[k].volume, solids[k].volume - solids[k - 1].volume, rel=1e-12)
    for k in range(3):
        assert close(sum(p.volume for p in pbs[: k + 1]), pits[k].volume, rel=1e-12)
    edges = np.arange(200.0, 310.0, 5.0)
    by_bench = sum(p.volume_by_elevation(edges) for p in pbs)
    assert np.allclose(by_bench, pits[2].volume_by_elevation(edges), rtol=0, atol=1e-9 * pits[2].volume)
    lo, hi = pits[0].bounds
    points = np.random.default_rng(5).uniform(lo, hi, (20_000, 3))
    in_first = pits[0].locate(points) == -1
    assert in_first.sum() > 1000
    assert (pbs[1].locate(points[in_first]) == 1).all() and (pbs[2].locate(points[in_first]) == 1).all()


def test_block_fractions_of_each_pushback_add_up_block_by_block(stages):
    pbs, final = stages["pushbacks"], stages["pits"][2]
    size = np.array([10.0, 10.0, 5.0])
    centres = block_grid(final, size)
    total = final.block_fraction(centres, size)
    parts = [p.block_fraction(centres, size) for p in pbs]
    assert np.abs(total - sum(parts)).max() < 1e-9
    for f, m in zip(parts + [total], pbs + [final], strict=True):
        assert ((f >= 0) & (f <= 1)).all()
        assert close((f * size.prod()).sum(), m.volume, rel=1e-9)


def test_pushback_round_trips_through_mining_formats_and_chains(stages, tmp_path):
    pb = stages["pushbacks"][1]
    edges = np.arange(200.0, 310.0, 10.0)
    for name in ("pushback_tr.csv", "pushback.dxf", "pushback.ply"):
        wf.write(tmp_path / name, pb)
        back = wf.read_mesh(tmp_path / name)
        assert np.array_equal(back.coords, pb.coords) and np.array_equal(back.triangles, pb.triangles)
        assert back.is_solid
        assert np.array_equal(back.volume_by_elevation(edges), pb.volume_by_elevation(edges))
        lower, upper = back.split_by_plane(O + [0, 0, 250.0], UP)
        assert lower.is_solid and upper.is_solid
        assert close(lower.volume + upper.volume, pb.volume, rel=1e-12)


def test_cut_and_fill_adds_up_over_a_partition_of_the_plan():
    rng = np.random.default_rng(11)
    design = cone(170.0, 90.0, n=101)
    corners = np.array([[-300.0, -300.0], [300.0, -300.0], [300.0, 300.0], [-300.0, 300.0]])
    xy = np.concatenate([rng.uniform(-300, 300, (15_000, 2)), corners])
    asbuilt = wf.dtm(
        np.column_stack([O[:2] + xy, design.surface_z(O[:2] + xy) + rng.normal(0, 0.8, len(xy))])
    )
    whole = wf.cut_fill(asbuilt, design)
    xs = np.linspace(-300.0, 300.0, 13)
    seam = np.column_stack([xs, np.where(np.arange(13) % 2 == 0, 7.0, -11.0)])
    north = np.vstack([seam, [[300.0, 300.0], [-300.0, 300.0]]])
    south = np.vstack([seam[::-1], [[-300.0, -300.0], [300.0, -300.0]]])
    hole = np.array([[-50.0, -50.0], [60.0, -40.0], [50.0, 70.0], [-40.0, 50.0]])

    def within(*rings):
        return wf.cut_fill(asbuilt, design, boundary=wf.Polylines([O[:2] + r for r in rings], closed=True))

    halves = [within(north), within(south)]
    holed, plug = within(corners, hole), within(hole)
    shoelace = 0.5 * abs(
        np.dot(hole[:, 0], np.roll(hole[:, 1], -1)) - np.dot(hole[:, 1], np.roll(hole[:, 0], -1))
    )
    assert close(plug["area"], shoelace, rel=1e-12)
    for key in ("cut", "fill", "net", "cut_area", "fill_area", "area"):
        scale = max(abs(whole[key]), whole["area"] * 1e-3)
        assert abs(sum(h[key] for h in halves) - whole[key]) <= 1e-10 * scale
        assert abs(holed[key] + plug[key] - whole[key]) <= 1e-10 * scale


def test_clip_at_utm_scale_matches_the_same_clip_at_the_local_origin():
    topo, _ = contour_dtm(spot=4000)
    solid = cone(170.0, 90.0, n=61).close_to_elevation(600.0)
    local = [m.with_vertices(m.coords - O) for m in (solid, topo)]
    for m, near in zip((solid, topo), local, strict=True):
        assert np.array_equal(near.coords + O, m.coords)
    far = solid.clip_by_surface(topo)
    here = local[0].clip_by_surface(local[1])
    for a, b in zip(far, here, strict=True):
        assert a.is_solid and b.is_solid
        assert a.n_faces == b.n_faces and a.n_vertices == b.n_vertices
        assert close(a.volume, b.volume, rel=1e-12)
        assert np.abs(np.sort(a.coords - O, axis=0) - np.sort(b.coords, axis=0)).max() < 1e-9


def test_block_fractions_refine_into_sub_blocks_and_rotate():
    topo, _ = contour_dtm(spot=4000)
    pit, _ = cone(170.0, 90.0, n=61).close_to_elevation(600.0).clip_by_surface(topo)
    size = np.array([10.0, 10.0, 10.0])
    centres = block_grid(pit, size)
    parent = pit.block_fraction(centres, size)
    offsets = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)]) * size / 4
    children = pit.block_fraction((centres[:, None] + offsets).reshape(-1, 3), size / 2).reshape(-1, 8)
    assert np.abs(parent * size.prod() - children.sum(1) * size.prod() / 8).max() < 1e-9 * size.prod()
    assert close((parent * size.prod()).sum(), pit.volume, rel=1e-9)
    theta = np.radians(30.0)
    rot = np.array(
        [[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]]
    )
    lo, hi = (np.asarray(b) for b in pit.bounds)
    mid = (lo + hi) / 2
    g = np.arange(-40, 41) * 10.0
    axes = np.meshgrid(g, g, np.arange(lo[2] - mid[2] - 5.0, hi[2] - mid[2] + 10.0, 10.0), indexing="ij")
    rotated = mid + np.column_stack([a.ravel() for a in axes]) @ rot.T
    f = pit.block_fraction(rotated, size, rotation=rot)
    assert ((f >= 0) & (f <= 1)).all()
    assert close((f * size.prod()).sum(), pit.volume, rel=1e-9)


def cliff(undercut, low=300.0, high=330.0, half=300.0, step=10.0, origin=O):
    """Topography on a 10 m grid: `low` west of x = 0 and `high` east of it, joined by a cliff that is
    vertical when `undercut` is 0 and overhangs the low ground by `undercut` otherwise."""
    ys = np.arange(-half, half + step / 2, step)
    ny, v, f = len(ys), [], []

    def plate(xs, z):
        base = len(v)
        v.extend([x, y, z] for x in xs for y in ys)
        for i in range(len(xs) - 1):
            for j in range(ny - 1):
                a = base + i * ny + j
                f.extend([[a, a + ny, a + ny + 1], [a, a + ny + 1, a + 1]])
        return base

    west = np.arange(-half, step / 2, step)
    bottom = plate(west, low) + (len(west) - 1) * ny
    east = np.arange(0.0, half + step / 2, step)
    east[0] = -undercut
    top = plate(east, high)
    for j in range(ny - 1):
        f.extend([[bottom + j, top + j + 1, bottom + j + 1], [bottom + j, top + j, top + j + 1]])
    return wf.Mesh(np.asarray(v) + origin, np.asarray(f), length_unit="m")


@pytest.mark.parametrize("undercut", [0.0, 4.0])
def test_topography_with_a_cliff_clips_benches_exactly(undercut):
    topo = cliff(undercut)
    s = topo.validate().summary
    assert s["shells"] == 1 and s["self_intersections"] == 0 and s["non_manifold_edges"] == 0
    block = wf.datasets.box(O + [-100, -100, 250], O + [100, 100, 400])
    below, above = block.clip_by_surface(topo)
    assert below.is_solid and above.is_solid
    expected = 200 * 100 * 50 + 200 * 100 * 80 + undercut * 30 / 2 * 200
    assert close(below.volume, expected, rel=1e-12)
    assert close(below.volume + above.volume, block.volume, rel=1e-12)
    edges = np.arange(250.0, 405.0, 5.0)
    benches = below.volume_by_elevation(edges)
    slab = [200 * 200 * 5 if z < 300 else 200 * 100 * 5 + undercut * (z + 2.5 - 300) / 30 * 200 * 5
            if z < 330 else 0.0 for z in edges[:-1]]  # fmt: skip
    assert np.allclose(benches, slab, rtol=1e-12, atol=1e-6)
    lower, upper = below.split_by_plane(O + [0, 0, 300.0], UP)
    assert lower.is_solid and upper.is_solid and close(lower.volume, 200 * 200 * 50, rel=1e-12)
    face = above.slice_by_plane(O + [0, 0, 315.0], UP)
    assert face.closed.all() and close(face.area().sum(), (100 - undercut / 2) * 200, rel=1e-12)
    probe = O[:2] + np.array([[0.0, 5.0], [-undercut / 2, 5.0], [-1e-9, 5.0], [-undercut - 1e-9, 5.0]])
    assert topo.surface_z(probe).tolist() == [330.0, 330.0, 330.0 if undercut else 300.0, 300.0]
    if undercut:
        with pytest.raises(wf.InvalidInput, match="not a height field"):
            wf.cut_fill(flat(300.0), topo)
    else:
        cf = wf.cut_fill(flat(300.0), topo)
        assert cf["fill"] == 30.0 * 300 * 600 and cf["cut"] == 0 and cf["fill_area"] == 300 * 600


def test_large_pit_shell_clipped_by_a_large_topography():
    topo, _ = contour_dtm(spot=60_000)
    design = cone(170.0, 90.0, n=301)
    assert topo.n_faces > 1e5 and design.n_faces > 1e5
    solid = design.close_to_elevation(600.0)
    below, above = solid.clip_by_surface(topo)
    assert below.is_solid and above.is_solid
    assert close(below.volume + above.volume, solid.volume, rel=1e-12)
    assert close(below.volume, wf.cut_fill(topo, design)["cut"], rel=1e-12)
    assert close(below.volume_by_elevation(bench_edges(below, 5.0)).sum(), below.volume, rel=1e-12)


SCENARIO = """
import hashlib, sys
import numpy as np
sys.path.insert(0, {tests!r})
import test_acceptance_open_pit_stages as t
ground = t.flat(305.0, half=2000.0)
solids = [t.benched_pit(t.ring_dirs(48), 12, east=e)[0].close_to_elevation(600.0) for e in (0.0, 0.3)]
pb = (solids[1] - solids[0]).clip_by_surface(ground)[0]
pit = solids[1].clip_by_surface(ground)[0]
size = np.array([10.0, 10.0, 5.0])
f = pb.block_fraction(t.block_grid(pit, size), size)
h = hashlib.sha256()
for a in (pb.coords, pb.triangles, pit.coords, pit.triangles, f, pb.volume_by_elevation(np.arange(200.0, 310.0, 5.0)),
          pit.slice_by_plane(t.O, t.UP, offsets=[250.0]).coords):
    h.update(np.ascontiguousarray(a).tobytes())
print(h.hexdigest())
"""


def test_stage_results_do_not_depend_on_thread_count():
    code = SCENARIO.format(tests=str(Path(__file__).parent))
    digests = set()
    for threads in ("1", "3", "4"):
        env = dict(os.environ, RAYON_NUM_THREADS=threads)
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        )
        digests.add(out.stdout.strip())
    assert len(digests) == 1


@pytest.mark.skip(reason="acceptance: clipping an octagonal pit by flat ground folds the cap at utm scale")
def test_octagonal_pit_clipped_by_flat_ground_is_a_solid():
    pit, rings = benched_pit(OCTAGON, 1)
    solid = pit.close_to_elevation(600.0)
    lower, upper = solid.split_by_plane(O + [0, 0, 205.0], UP)
    assert lower.is_solid and upper.is_solid
    below, above = solid.clip_by_surface(flat(205.0, half=2000.0))
    assert below.is_solid and above.is_solid
    assert close(below.volume, pit_volume(8, rings, 200.0, 205.0), rel=SHAPE)


@pytest.mark.skip(reason="acceptance: union of a previous pit and its pushback is not a solid")
def test_previous_pit_and_pushback_unite_into_the_next_pit():
    ground = flat(215.0, half=2000.0, origin=np.zeros(3))
    first, second = (
        benched_pit(SQUARE, 2, east=e, origin=np.zeros(3))[0].close_to_elevation(600.0) for e in (0.0, 0.3)
    )
    pit = first.clip_by_surface(ground)[0]
    pushback = (second - first).clip_by_surface(ground)[0]
    assert pit.is_solid and pushback.is_solid
    union = pit | pushback
    assert union.is_solid
    assert close(union.volume, second.clip_by_surface(ground)[0].volume, rel=1e-12)


def load_pair(name):
    z = np.load(DATA / name)
    return wf.Mesh(z["av"], z["at"]), wf.Mesh(z["bv"], z["bt"])


@pytest.mark.skip(reason="acceptance: booleans of touching pushback slivers are not solids")
def test_booleans_of_touching_pushback_slivers_are_solids():
    pit, pushback = load_pair("pushback_pair_touching.npz")
    assert pit.is_solid and pushback.is_solid
    both = pit.volume + pushback.volume
    meet = pit & pushback
    assert meet.n_faces == 0 or (meet.is_solid and meet.volume < 1e-9 * both)
    for m, v in ((pit | pushback, both), (pit - pushback, pit.volume), (pushback - pit, pushback.volume)):
        assert m.is_solid and close(m.volume, v, rel=1e-9)


RUNAWAY = """
import sys
import numpy as np
import wireframes as wf
z = np.load({path!r})
a, b = wf.Mesh(z["av"], z["at"]), wf.Mesh(z["bv"], z["bt"])
m = a & b
solid = m.n_faces == 0 or m.is_solid
print(m.n_faces, solid, 0.0 if m.n_faces == 0 else m.volume if solid else float("nan"))
"""


def limit_memory():
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))


def intersect_runaway_pair():
    code = RUNAWAY.format(path=str(DATA / "pushback_pair_runaway.npz"))
    env = dict(os.environ, RAYON_NUM_THREADS="2", MALLOC_ARENA_MAX="2")
    out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=300, check=False,
                         preexec_fn=limit_memory)  # fmt: skip
    assert out.returncode == 0, out.stderr[-2000:]
    return out.stdout.split()


def test_intersection_of_previous_pit_and_pushback_pieces_stays_small():
    intersect_runaway_pair()


@pytest.mark.skip(reason="acceptance: rounded slivers between walls that agree to an ulp self-intersect")
def test_intersection_of_previous_pit_and_pushback_pieces_is_a_thin_solid():
    _, solid, volume = intersect_runaway_pair()
    assert solid == "True" and float(volume) < 1e-6


@pytest.mark.skip(reason="acceptance: union and difference of previous pit and pushback pieces fail rounding")
def test_union_and_difference_of_previous_pit_and_pushback_pieces():
    pit, pushback = load_pair("pushback_pair_runaway.npz")
    assert pit.is_solid and pushback.is_solid
    for m, v in ((pit | pushback, pit.volume + pushback.volume), (pit - pushback, pit.volume),
                 (pushback - pit, pushback.volume)):  # fmt: skip
        assert m.is_solid and close(m.volume, v, rel=1e-9)


def test_dump_beyond_the_topography_is_not_called_below_it():
    c = np.array([[-20.0, -20.0], [20.0, -20.0], [20.0, 20.0], [-20.0, 20.0]])
    ground = wf.dtm(np.column_stack([c, np.zeros(4)]))
    dump = wf.datasets.box([99, -1, 49], [101, 1, 51])
    try:
        below, above = dump.clip_by_surface(ground)
    except wf.InvalidInput:
        return
    assert below.n_faces == 0 and close(above.volume, dump.volume)


def test_an_empty_pushback_chains_into_further_operations():
    box = wf.datasets.box
    pit = box(O + [0, 0, 200], O + [10, 10, 210])
    empty = pit & box(O + [50, 50, 200], O + [60, 60, 210])
    assert empty.n_faces == 0
    assert empty.volume == 0.0
    assert empty.volume_by_elevation([200.0, 205.0, 210.0]).tolist() == [0.0, 0.0]
    assert empty.block_fraction(O + [[5.0, 5.0, 205.0]], [10, 10, 10]).tolist() == [0.0]
    assert (pit | empty).volume == pit.volume and (pit - empty).volume == pit.volume
    assert (pit & empty).n_faces == 0 and (empty - pit).n_faces == 0


def test_dtm_never_silently_drops_a_conflicting_height():
    corners = np.array([[0.0, 0.0, 300.0], [10.0, 0.0, 300.0], [0.0, 10.0, 300.0], [10.0, 10.0, 300.0]])
    crest = wf.Polylines([O + [[2.0, 2.0, 305.0], [8.0, 8.0, 305.0]]], closed=[False])
    spots = np.vstack([corners, [[2.0, 2.0, 300.0]]]) + O
    try:
        topo = wf.dtm(spots, breaklines=crest)
    except wf.InvalidInput:
        pass
    else:
        assert np.array_equal(topo.surface_z(crest.coords[:, :2]), crest.coords[:, 2])
    with pytest.raises(wf.InvalidInput):
        wf.dtm(np.vstack([corners, [[5.0, 5.0, 301.0], [5.0, 5.0, 309.0]]]) + O)


def test_contour_at_the_pit_floor_traces_the_floor_toe():
    design = cone(170.0, 90.0, n=41, half=200.0)
    toe = design.contours([170.0])
    assert len(toe) == 1 and toe.closed.all()
    corners = design.coords[design.triangles]
    flat = (corners[:, :, 2] == 170.0).all(axis=1)
    a, b, c = (corners[flat, k, :2] for k in range(3))
    u, w = b - a, c - a
    floor = 0.5 * np.abs(u[:, 0] * w[:, 1] - u[:, 1] * w[:, 0]).sum()
    assert close(toe.area()[0], floor, rel=1e-12)
