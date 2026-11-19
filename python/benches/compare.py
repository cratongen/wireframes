"""Time wireframes against trimesh, manifold3d and scipy on the same seeded inputs.

Every round runs each implementation once, in turn, so warm-up and drift favour none of
them. The table gives the median over rounds and the range, in milliseconds, with the
1-minute load average before and after each kernel. Before timing, each kernel checks
that the implementations agree on its result.

Run one kernel at a time with ``--only``, for example::

    mise run bench:compare -- --only union --rounds 5
"""

from __future__ import annotations

import argparse
import importlib.metadata as md
import os
import platform
import statistics
import sys
import time
from pathlib import Path

import manifold3d as mf
import numpy as np
import scipy.spatial
import trimesh
import wireframes as wf

SIZES = {"1k": (25, 20), "10k": (100, 50), "100k": (500, 100)}
SHIFT = (0.37, 0.29, 0.11)
OFFSETS = -0.49 + 0.02 * np.arange(50)
TRIMESH_BATCH = 1000


def torus(major, minor):
    """The vertices and triangles of ``wireframes_core::fixtures::torus(major, minor)``."""
    th = np.arange(major) / major * (2 * np.pi)
    ph = np.arange(minor) / minor * (2 * np.pi)
    r = 2.0 + 0.5 * np.cos(ph)[None, :]
    x, y = r * np.cos(th)[:, None], r * np.sin(th)[:, None]
    v = np.stack([x, y, np.broadcast_to(0.5 * np.sin(ph), x.shape)], axis=-1).reshape(-1, 3)
    i, j = np.meshgrid(np.arange(major), np.arange(minor), indexing="ij")

    def vid(i, j):
        return (i % major) * minor + (j % minor)

    a, b, c, d = vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)
    f = np.stack([np.stack([a, b, c], -1), np.stack([a, c, d], -1)], axis=2).reshape(-1, 3)
    return np.ascontiguousarray(v), f.astype(np.int64)


def to_manifold(v, f):
    return mf.Manifold(mf.Mesh64(np.ascontiguousarray(v, np.float64), np.ascontiguousarray(f, np.uint64)))


def to_trimesh(v, f):
    return trimesh.Trimesh(v, f, process=False)


def points(seed, n, lo, hi):
    return np.random.default_rng(seed).uniform(lo, hi, size=(n, 3))


def close(x, y, rel=1e-9):
    return abs(x - y) <= rel * max(abs(x), abs(y))


def boolean_kernels(op):
    for size, (major, minor) in SIZES.items():
        v, f = torus(major, minor)
        w = v + SHIFT
        a, b = wf.Mesh(v, f), wf.Mesh(w, f)
        ma, mb = to_manifold(v, f), to_manifold(w, f)
        ta, tb = to_trimesh(v, f), to_trimesh(w, f)
        ours = getattr(wf, op)
        theirs = {
            "union": lambda ma=ma, mb=mb: (ma + mb).to_mesh64(),
            "difference": lambda ma=ma, mb=mb: (ma - mb).to_mesh64(),
        }[op]
        engine = {"union": trimesh.boolean.union, "difference": trimesh.boolean.difference}[op]

        def check(ours=ours, engine=engine, a=a, b=b, ma=ma, mb=mb, ta=ta, tb=tb):
            vol = ours(a, b).volume
            vm = (ma + mb if op == "union" else ma - mb).volume()
            vt = engine([ta, tb], engine="manifold").volume
            assert close(vol, vm, 1e-6) and close(vol, vt, 1e-6), (vol, vm, vt)
            return f"volume {vol:.6f}, manifold3d {vm:.6f}, trimesh {vt:.6f}"

        yield (
            f"{op}",
            f"two {size}-triangle tori",
            check,
            {
                "wireframes": lambda ours=ours, a=a, b=b: ours(a, b),
                "manifold3d": theirs,
                "trimesh": lambda engine=engine, ta=ta, tb=tb: engine([ta, tb], engine="manifold"),
            },
        )


def validate_kernels():
    for size, (major, minor) in SIZES.items():
        v, f = torus(major, minor)
        a = wf.Mesh(v, f)

        def check(a=a, v=v, f=f):
            assert a.is_solid
            assert to_manifold(v, f).status() == mf.Error.NoError
            assert to_trimesh(v, f).is_volume
            return "all three accept the torus"

        yield (
            "validate (manifold3d: construct; trimesh: is_volume)",
            f"{size}-triangle torus",
            check,
            {
                "wireframes": lambda a=a: a.validate(),
                "manifold3d": lambda v=v, f=f: to_manifold(v, f).status(),
                "trimesh": lambda v=v, f=f: to_trimesh(v, f).is_volume,
            },
        )


def slice_kernel():
    v, f = torus(*SIZES["100k"])
    a, m, t = wf.Mesh(v, f), to_manifold(v, f), to_trimesh(v, f)

    def check():
        parts = len(a.slice_by_plane([0, 0, 0], [0, 0, 1], offsets=OFFSETS).parts)
        loops = sum(len(m.slice(h).to_polygons()) for h in OFFSETS)
        lines = trimesh.intersections.mesh_multiplane(t, [0, 0, 0], [0, 0, 1], OFFSETS)[0]
        assert parts == loops == 100 and len(lines) == 50, (parts, loops, len(lines))
        return f"{parts} closed loops from wireframes and manifold3d"

    yield (
        "slice by 50 planes",
        "100k-triangle torus",
        check,
        {
            "wireframes": lambda: a.slice_by_plane([0, 0, 0], [0, 0, 1], offsets=OFFSETS),
            "manifold3d": lambda: [m.slice(h).to_polygons() for h in OFFSETS],
            "trimesh": lambda: trimesh.intersections.mesh_multiplane(t, [0, 0, 0], [0, 0, 1], OFFSETS),
        },
    )


def query_kernels():
    v, f = torus(*SIZES["100k"])
    a = wf.Mesh(v, f)
    p = points(1, 100_000, [-2.6, -2.6, -0.6], [2.6, 2.6, 0.6])

    def batched(fn, t):
        return np.concatenate([fn(t, p[i : i + TRIMESH_BATCH]) for i in range(0, len(p), TRIMESH_BATCH)])

    def tm_contains():
        return batched(lambda t, q: t.contains(q), to_trimesh(v, f))

    def tm_distance():
        return batched(trimesh.proximity.signed_distance, to_trimesh(v, f))

    def check_contains():
        ours, theirs = a.contains(p), tm_contains()
        differ = int((ours != theirs).sum())
        assert differ <= 10, differ
        return f"{int(ours.sum())} inside; trimesh differs on {differ}"

    yield (
        f"contains (fresh index each call; trimesh in batches of {TRIMESH_BATCH})",
        "100k points, 100k-triangle torus",
        check_contains,
        {"wireframes": lambda: wf.Mesh(v, f).contains(p), "trimesh": tm_contains},
    )

    def check_distance():
        ours, theirs = a.signed_distance(p), -tm_distance()
        inside = a.contains(p)
        tri = v[f]
        worst = 0.0
        for k in np.argsort(-np.abs(ours - theirs))[:50]:
            cp = trimesh.triangles.closest_point(tri, np.repeat(p[k][None], len(f), 0))
            brute = np.linalg.norm(cp - p[k], axis=1).min() * (-1 if inside[k] else 1)
            worst = max(worst, abs(ours[k] - brute))
        assert worst <= 1e-12, worst
        return (
            f"on the 50 points where they differ most, wireframes is within {worst:.1e} of brute force; "
            f"trimesh differs from wireframes by up to {np.abs(ours - theirs).max():.1e}"
        )

    yield (
        f"signed distance (fresh index each call; trimesh in batches of {TRIMESH_BATCH})",
        "100k points, 100k-triangle torus",
        check_distance,
        {"wireframes": lambda: wf.Mesh(v, f).signed_distance(p), "trimesh": tm_distance},
    )

    n = (100, 100, 50)
    lo, hi = np.array([-2.5, -2.5, -0.5]), np.array([2.5, 2.5, 0.5])
    size = (hi - lo) / n
    k, j, i = np.meshgrid(*(np.arange(c) for c in n[::-1]), indexing="ij")
    centroids = lo + (np.stack([i, j, k], -1).reshape(-1, 3) + 0.5) * size

    def check_fraction():
        total = a.block_fraction(centroids, size).sum() * size.prod()
        assert close(total, a.volume), (total, a.volume)
        return f"fractions sum to the torus volume {a.volume:.6f}"

    yield (
        "block fraction (fresh index each call)",
        "100 x 100 x 50 blocks, 100k-triangle torus",
        check_fraction,
        {"wireframes": lambda: wf.Mesh(v, f).block_fraction(centroids, size)},
    )


def dtm_kernel():
    p = points(4, 100_000, [0, 0, 0], [1000, 1000, 1])
    p[:, 2] += 100 + 0.05 * p[:, 0] - 0.02 * p[:, 1]

    def check():
        ours = wf.dtm(p)
        theirs = scipy.spatial.Delaunay(p[:, :2])
        assert len(ours.triangles) == len(theirs.simplices), (len(ours.triangles), len(theirs.simplices))
        return f"{len(ours.triangles)} triangles from both"

    yield (
        "dtm (scipy: Delaunay in plan)",
        "100k points",
        check,
        {"wireframes": lambda: wf.dtm(p), "scipy": lambda: scipy.spatial.Delaunay(p[:, :2])},
    )


def io_kernels(tmp):
    v, f = torus(*SIZES["100k"])
    a, t = wf.Mesh(v, f), to_trimesh(v, f)
    names = {
        "obj": "m.obj",
        "ply": "m.ply",
        "off": "m.off",
        "stl": "m.stl",
        "dxf": "m.dxf",
        "vtk": "m.vtk",
        "vtu": "m.vtu",
        "datamine csv": "m_tr.csv",
    }
    for fmt, name in names.items():
        ours, theirs = tmp / "wf" / name, tmp / "trimesh" / name
        ours.parent.mkdir(parents=True, exist_ok=True)
        theirs.parent.mkdir(parents=True, exist_ok=True)
        has_trimesh = fmt in ("obj", "ply", "off", "stl")

        def check(ours=ours, theirs=theirs, has_trimesh=has_trimesh):
            wf.write(ours, a)
            assert len(wf.read_mesh(ours).triangles) == len(f)
            if has_trimesh:
                t.export(theirs)
                assert len(trimesh.load(theirs, process=False).faces) == len(f)
            return f"{ours.stat().st_size} bytes written by wireframes"

        write = {"wireframes": lambda ours=ours: wf.write(ours, a)}
        read = {"wireframes": lambda ours=ours: wf.read_mesh(ours)}
        if has_trimesh:
            write["trimesh"] = lambda theirs=theirs: t.export(theirs)
            read["trimesh"] = lambda theirs=theirs: trimesh.load(theirs, process=False)
        yield f"write {fmt}", "100k-triangle torus", check, write
        yield f"read {fmt}", "100k-triangle torus", lambda: "written above", read


def run(kernel, rounds):
    name, dataset, check, impls = kernel
    load0 = os.getloadavg()[0]
    note = check()
    times = {k: [] for k in impls}
    for k, fn in impls.items():
        fn()
    for _ in range(rounds):
        for k, fn in impls.items():
            t0 = time.perf_counter()
            fn()
            times[k].append(1e3 * (time.perf_counter() - t0))
    load1 = os.getloadavg()[0]
    return name, dataset, times, note, (load0, load1)


def cell(ms):
    if not ms:
        return "n/a"
    med, lo, hi = statistics.median(ms), min(ms), max(ms)
    digits = 3 if med < 10 else 1 if med < 1000 else 0
    return f"{med:.{digits}f} ({lo:.{digits}f} to {hi:.{digits}f})"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--rounds", type=int, default=5, help="timed rounds per kernel (default 5)")
    parser.add_argument("--only", default="", help="run kernels whose name contains this text")
    parser.add_argument("--tmp", type=Path, default=Path("target/compare"), help="folder for I/O files")
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error("--rounds must be at least 1")
    kernels = [
        *boolean_kernels("union"),
        *boolean_kernels("difference"),
        *validate_kernels(),
        *slice_kernel(),
        *query_kernels(),
        *dtm_kernel(),
        *io_kernels(args.tmp),
    ]
    kernels = [k for k in kernels if args.only in k[0]]
    if not kernels:
        parser.error(f"no kernel name contains {args.only!r}")
    versions = ", ".join(f"{p} {md.version(p)}" for p in ("numpy", "trimesh", "manifold3d", "scipy", "rtree"))
    cpu = (
        next(
            (
                s.split(":", 1)[1].strip()
                for s in Path("/proc/cpuinfo").read_text().splitlines()
                if "model name" in s
            ),
            platform.machine(),
        )
        if Path("/proc/cpuinfo").exists()
        else platform.machine()
    )
    print(f"# {platform.platform()}, {cpu}, {os.cpu_count()} logical CPUs")
    print(f"# python {sys.version.split()[0]}, wireframes {wf.__version__}, {versions}")
    print(f"# RAYON_NUM_THREADS={os.environ.get('RAYON_NUM_THREADS', 'unset')}, rounds={args.rounds}")
    print("| kernel | dataset | wireframes | trimesh | manifold3d | scipy | load before, after | check |")
    print("|---|---|---|---|---|---|---|---|")
    for k in kernels:
        name, dataset, times, note, load = run(k, args.rounds)
        cols = " | ".join(cell(times.get(i, [])) for i in ("wireframes", "trimesh", "manifold3d", "scipy"))
        print(f"| {name} | {dataset} | {cols} | {load[0]:.1f}, {load[1]:.1f} | {note} |", flush=True)


if __name__ == "__main__":
    main()
