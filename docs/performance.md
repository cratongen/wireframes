# Performance

Every number on this page comes from benchmarks run on one machine, one command at a
time: on 2026-11-19 at commit `5d5ab30`, and on 2026-11-20 for the rounding section. Cells marked "not measured" have no benchmark
behind them.

## Machine and versions

| Item | Value |
|---|---|
| CPU | 13th Gen Intel Core i9-13980HX, 32 logical CPUs |
| Memory | 31 GiB; every run capped at 4 GiB with `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0` |
| OS | Linux 7.2.15 (Omarchy), glibc 2.44, laptop with a desktop session open |
| Rust | rustc 1.97.1, release profile, criterion 0.8.2, rayon 1.12.0 |
| wireframes | 0.1.0, code at commit `5d5ab30` (rounding section: `ba8a660` and `cf79bfb`) |
| Threads | `RAYON_NUM_THREADS` unset, so wireframes uses 32 threads |

Results do not depend on the thread count, but timings do.

## Inputs

All inputs are seeded and synthetic:

- **Torus**: `fixtures::torus(major, minor)`, radii 2 and 0.5, `2 * major * minor`
  triangles: 1k = `torus(25, 20)`, 10k = `torus(100, 50)`, 100k = `torus(500, 100)`.
- **Boolean pairs**: a torus and a copy translated by `(0.37, 0.29, 0.11)`.
- **Points**: 1e5 uniform points in `[-2.6, 2.6] x [-2.6, 2.6] x [-0.6, 0.6]`
  (SplitMix64 seed 1).
- **Block model**: 100 x 100 x 50 blocks of 0.05 x 0.05 x 0.02 over the torus bounds.
- **Slices**: 50 horizontal planes at z = -0.49, -0.47, ..., 0.49.
- **DTM**: 1e5 uniform points over a 1000 x 1000 square, z a plane plus noise of 1.

## Commands

```sh
export CARGO_BUILD_JOBS=8 RUST_TEST_THREADS=4
systemd-run --user --scope -q -p MemoryMax=4G -p MemorySwapMax=0 mise run -j 1 bench
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

## Two fixes measured before and after

Stage timers on the 100k union showed validation taking most of the time: a boolean validates both inputs and twice validates its own rounded result.
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

## Rounding to f64, measured before and after

The rounding work of `docs/design/rounding.md` (method C, conforming of inputs, flat-run
removal, descents and slot closing) landed in commits `257bbba` to `cf79bfb`. These
benches ran on 2026-11-20 on the same machine against bench binaries built from
`ba8a660`, the commit before the first kernel change, and from `cf79bfb`. Both commits
already run `validate` on each boolean result, so the boolean rows below include that
call. The bench sources are identical at the two commits.

Each binary ran three times, alternating before and after, with
`--bench --discard-baseline --noplot 'torus_pair/(union|difference|validate|slice)'`,
10 samples per run, under `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0`.
The table gives the median of the three criterion medians. The 1-minute load average
ranged from 7.9 to 16.5 across the six runs, most of it from the benchmarks' own 32
threads.

| Kernel | Input | `ba8a660` | `cf79bfb` |
|---|---|---|---|
| `union` | two 1k tori | 31.4 ms | 41.1 ms |
| `union` | two 10k tori | 155 ms | 200 ms |
| `union` | two 100k tori | 1.09 s | 1.47 s |
| `difference` | two 1k tori | 31.3 ms | 40.9 ms |
| `difference` | two 10k tori | 149 ms | 194 ms |
| `difference` | two 100k tori | 1.17 s | 1.51 s |
| `validate` | 1k torus | 2.47 ms | 2.46 ms |
| `validate` | 10k torus | 14.4 ms | 14.5 ms |
| `validate` | 100k torus | 126 ms | 125 ms |
| `slice_by_plane`, 50 planes | 100k torus | 86.5 ms | 86.3 ms |

On these torus pairs `union` and `difference` take 29 to 35 percent longer at `cf79bfb` (ratios of the medians above).
Which stage costs that time: not measured. The cost of the retries (conforming,
descents, slot closing) on inputs that need them: not measured.

wireframes decides topology with exact
predicates that always evaluate with expansions (the floating-point filter in front of
them is not written yet), and every boolean validates three meshes before it returns.
Speed with the filter: not measured.

## Not measured

- Other libraries as a reference: not measured.
- The share of the `validate` call that booleans, `clip_by_surface` and
  `split_by_plane` run on each result: not measured on its own. The `5d5ab30` boolean
  rows in the first table predate that call; the rows in the rounding table include it.
- `clip_by_surface` and `split_by_plane`: no bench exists, not measured.
- Read and write speed of files larger than the page cache: not measured.
