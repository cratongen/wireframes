# Performance

Every number on this page comes from benchmarks run on 2026-11-19 on one machine, in
one session, one command at a time. Cells marked "not measured" have no benchmark
behind them.

## Machine and versions

| Item | Value |
|---|---|
| CPU | 13th Gen Intel Core i9-13980HX, 32 logical CPUs |
| Memory | 31 GiB; every run capped at 4 GiB with `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0` |
| OS | Linux 7.2.15 (Omarchy), glibc 2.44, laptop with a desktop session open |
| Rust | rustc 1.97.1, release profile, criterion 0.8.2, rayon 1.12.0 |
| Python | 3.13.16, numpy 2.5.5, trimesh 5.2.1 (with rtree 1.4.1), manifold3d 3.5.5, scipy 1.18.1 |
| wireframes | 0.1.0, code at commit `5d5ab30` |
| Threads | `RAYON_NUM_THREADS` unset, so wireframes uses 32 threads; trimesh, manifold3d and scipy run with their defaults |

Results do not depend on the thread count, but timings do.

## Inputs

All inputs are seeded and synthetic, and Rust and Python build the same ones:

- **Torus**: `fixtures::torus(major, minor)`, radii 2 and 0.5, `2 * major * minor`
  triangles: 1k = `torus(25, 20)`, 10k = `torus(100, 50)`, 100k = `torus(500, 100)`.
  `python/benches/compare.py` builds the same grid with numpy.
- **Boolean pairs**: a torus and a copy translated by `(0.37, 0.29, 0.11)`.
- **Points**: 1e5 uniform points in `[-2.6, 2.6] x [-2.6, 2.6] x [-0.6, 0.6]`
  (SplitMix64 seed 1 in Rust, numpy `default_rng(1)` in Python).
- **Block model**: 100 x 100 x 50 blocks of 0.05 x 0.05 x 0.02 over the torus bounds.
- **Slices**: 50 horizontal planes at z = -0.49, -0.47, ..., 0.49.
- **DTM**: 1e5 uniform points over a 1000 x 1000 square, z a plane plus noise of 1.

## Commands

```sh
export CARGO_BUILD_JOBS=8 RUST_TEST_THREADS=4
systemd-run --user --scope -q -p MemoryMax=4G -p MemorySwapMax=0 mise run -j 1 bench
for k in union difference validate slice "block fraction" dtm contains "signed distance" "write " "read "; do
  systemd-run --user --scope -q -p MemoryMax=4G -p MemorySwapMax=0 \
    mise run -j 1 bench:compare -- --only "$k" --rounds 5
done
```

## Rust kernels (criterion)

Median time per call with criterion's 95% confidence interval, 10 samples each. The
1-minute load average was 4.2 before the run and 13.5 after it, most of it from the
benchmarks' own 32 threads.

| Kernel | Input | 1k | 10k | 100k |
|---|---|---|---|---|
| Arrangement (`arrangement_size`) | two tori | 13.7 ms (13.6 to 13.8) | 55.5 ms (54.7 to 55.9) | 218 ms (216 to 220) |
| `union` | two tori | 27.7 ms (27.5 to 28.1) | 139 ms (138 to 141) | 958 ms (950 to 968) |
| `difference` | two tori | 27.9 ms (27.6 to 28.1) | 135 ms (134 to 135) | 1.06 s (1.04 to 1.08) |
| `validate` | one torus | 2.53 ms (2.52 to 2.54) | 14.5 ms (14.4 to 14.7) | 126 ms (125 to 128) |

| Kernel | Input | Median (95% CI) |
|---|---|---|
| `slice_by_plane`, 50 planes | 100k torus | 86.3 ms (85.9 to 86.8) |
| `SolidIndex::new` (validates, then builds the BVH) | 100k torus | 156 ms (155 to 157) |
| `locate_all` | 1e5 points, 100k torus | 4.68 ms (4.66 to 4.72) |
| `signed_distance` | 1e5 points, 100k torus | 56.9 ms (56.4 to 57.2) |
| `block_fraction` | 100 x 100 x 50 blocks, 100k torus | 127 ms (126 to 127) |
| `trace_intervals` | 200-vertex trace, 100k torus | 66.5 ms (66.5 to 66.6) |
| `dtm` | 1e5 points | 61.6 ms (61.5 to 61.8) |

Readers and writers work on bytes in memory, on the 100k torus, with default options
(binary PLY, STL, VTK and VTU; text for the rest):

| Format | Write | Read |
|---|---|---|
| OBJ | 15.0 ms | 13.9 ms |
| PLY | 0.720 ms | 2.86 ms |
| OFF | 14.1 ms | 12.3 ms |
| STL | 3.29 ms | 6.94 ms |
| DXF (`3DFACE`) | 144 ms | 134 ms |
| VTK legacy | 4.12 ms | 1.91 ms |
| VTU | 4.40 ms | 2.06 ms |
| Datamine `_tr.csv`/`_pt.csv` | 24.6 ms | 56.5 ms |

## Against trimesh, manifold3d and scipy (Python)

`python/benches/compare.py` calls every library from Python on the same arrays. Each
round runs every implementation once, in turn. Cells give the median of 5 rounds in
milliseconds and the range, after one untimed warm-up call. Before timing, each kernel
checks that the libraries agree: equal volumes to 1e-6 for booleans, the same 100 loops
for slices, the same inside count for `contains`, the same triangle count for the DTM.

| Kernel | Input | wireframes | trimesh | manifold3d | scipy | Load before, after |
|---|---|---|---|---|---|---|
| union | two 1k tori | 27.9 (25.2 to 39.1) | 3.82 (2.59 to 4.16) | 2.32 (1.72 to 2.65) | n/a | 9.8, 9.8 |
| union | two 10k tori | 127 (123 to 134) | 13.7 (13.5 to 17.9) | 5.15 (4.83 to 5.38) | n/a | 9.8, 9.8 |
| union | two 100k tori | 941 (923 to 967) | 45.6 (44.1 to 49.9) | 17.2 (16.1 to 17.5) | n/a | 9.8, 9.1 |
| difference | two 1k tori | 27.4 (24.8 to 31.6) | 2.97 (2.47 to 3.75) | 1.77 (1.56 to 2.02) | n/a | 8.4, 8.4 |
| difference | two 10k tori | 130 (127 to 135) | 13.5 (13.2 to 15.1) | 4.92 (4.62 to 5.66) | n/a | 8.4, 8.4 |
| difference | two 100k tori | 1049 (1028 to 1069) | 46.2 (45.4 to 51.2) | 15.2 (14.9 to 17.0) | n/a | 8.4, 7.9 |
| validate | 1k torus | 2.32 (2.03 to 2.61) | 1.16 (0.65 to 1.29) | 0.764 (0.469 to 0.785) | n/a | 7.3, 7.3 |
| validate | 10k torus | 12.8 (11.7 to 13.7) | 3.91 (3.84 to 4.64) | 4.57 (4.06 to 5.01) | n/a | 7.3, 7.3 |
| validate | 100k torus | 117 (115 to 128) | 43.3 (42.0 to 45.9) | 19.3 (18.9 to 19.9) | n/a | 7.3, 7.3 |
| slice by 50 planes | 100k torus | 84.5 (77.1 to 90.2) | 123 (110 to 125) | 20.0 (19.0 to 22.7) | n/a | 7.3, 7.3 |
| contains | 1e5 points, 100k torus | 200 (195 to 210) | 52373 (50338 to 52528) | n/a | n/a | 5.4, 21.4 |
| signed distance | 1e5 points, 100k torus | 221 (212 to 223) | 47672 (47358 to 48153) | n/a | n/a | 22.3, 26.5 |
| block fraction | 100 x 100 x 50 blocks, 100k torus | 288 (278 to 299) | n/a | n/a | n/a | 6.9, 6.9 |
| dtm | 1e5 points | 61.4 (60.3 to 61.6) | n/a | n/a | 254 (253 to 259) | 6.9, 6.4 |

Files on disk (`target/compare`, page cache warm), 100k torus:

| Format | wireframes write | trimesh write | wireframes read | trimesh read |
|---|---|---|---|---|
| OBJ | 18.7 (17.6 to 19.2) | 98.9 (94.5 to 104) | 14.5 (14.3 to 15.1) | 67.1 (66.6 to 67.3) |
| PLY | 1.35 (1.24 to 1.41) | 1.03 (0.95 to 1.07) | 3.13 (3.08 to 3.29) | 0.684 (0.668 to 0.723) |
| OFF | 19.2 (18.3 to 20.7) | 100 (90.7 to 104) | 13.9 (13.6 to 14.3) | 104 (97.9 to 108) |
| STL | 5.21 (4.77 to 5.39) | 2.64 (2.39 to 2.67) | 8.23 (8.18 to 8.68) | 3.07 (2.50 to 3.44) |
| DXF | 166 (152 to 177) | n/a | 157 (151 to 159) | n/a |
| VTK legacy | 5.80 (4.77 to 7.64) | n/a | 2.36 (2.33 to 2.65) | n/a |
| VTU | 5.40 (5.12 to 5.68) | n/a | 2.60 (2.45 to 3.22) | n/a |
| Datamine CSV pair | 32.7 (30.0 to 34.1) | n/a | 64.2 (64.2 to 64.6) | n/a |

The load average was between 20.8 and 24.3 during the file runs, falling from the
multithreaded trimesh query runs before them.

What each column times:

- **Booleans**: wireframes validates both inputs and the result and computes the
  arrangement with exact predicates. manifold3d times `(a + b).to_mesh64()` on
  prebuilt `Manifold` objects, so its construction is left out. trimesh calls
  `trimesh.boolean.union` with `engine="manifold"`, which converts to manifold3d
  inside the call.
- **validate**: the three columns do different work. wireframes looks for every
  problem in its report, self-intersections included. manifold3d times
  `Manifold(Mesh64(...)).status()`, which checks that the mesh is manifold. trimesh
  times `Trimesh(...).is_volume`, which checks watertightness, winding and volume
  sign. Neither manifold3d nor trimesh looks for self-intersections.
- **contains, signed distance, block fraction**: every call builds a fresh
  `wf.Mesh`, so the time includes validation and the BVH (156 ms of it, from the
  criterion row for `SolidIndex::new`). trimesh also builds its index in each call.
  trimesh runs in batches of 1000 points, because one call on 1e4 points peaked at
  3.9 GB. Its signed distance differs from wireframes by up to 8.1e-6; brute force
  over all faces agrees with wireframes to within 2.7e-16 on the 50 points where
  the two differ most.
- **dtm**: scipy's `Delaunay` triangulates in plan and builds no mesh; it returns
  the same 199968 triangles.

## Two fixes measured before and after

The first comparison put wireframes booleans about 100 times behind manifold3d at
100k triangles. Stage timers on the 100k union showed validation taking most of the
time: a boolean validates both inputs and twice validates its own rounded result.
In `validate`, 640k candidate face pairs each ran the exact triangle-triangle test,
although for a pair that shares an edge or a corner and is not coplanar the answer
never depends on it. `validate` now runs that test only where it can change the
answer (commit `5b49623`). The arrangement also sorted 209k points with an exact
comparator on one thread; it now sorts them in parallel (commit `5d5ab30`).

Medians of 3 criterion runs per binary, alternating before and after, 10 samples per
run (`--bench --discard-baseline` on saved bench binaries):

| Kernel | Before | After |
|---|---|---|
| `validate`, 100k torus | 361 ms | 136 ms |
| `union`, two 100k tori | 1.89 s | 1.20 s |
| `difference`, two 100k tori | 2.31 s | 1.25 s |

The parallel sort, measured the same way against the first fix alone:

| Kernel | Before | After |
|---|---|---|
| arrangement, two 100k tori | 331 ms | 217 ms |
| `union`, two 100k tori | 1.06 s | 0.973 s |

The 1-minute load average ranged from 2.6 to 20.2 during these runs.

The remaining gap to manifold3d is structural. wireframes decides topology with exact
predicates that always evaluate with expansions (the floating-point filter in front of
them is not written yet), and every boolean validates three meshes before it returns.
Speed with the filter: not measured.

## Not measured

- PyVista or VTK as a reference: not measured.
- Read and write speed of files larger than the page cache: not measured.
- Block fraction, arrangement, trace intervals and the Datamine, DXF, VTK and VTU
  formats against another library: not measured; trimesh and manifold3d have no
  equivalent.
