# wireframes

`wireframes` is a Rust-core, Python-exposed library to read, write, plot, measure,
and combine 3D solids, surfaces, and strings for mining, geology, and petroleum.

## Scope and relationship to boitata

`wireframes` is a standalone geometry package with zero dependency on `boitata`.
Once `wireframes` reaches Milestone 7, `boitata` will depend on `wireframes` and delete
its internal mesh type (`crates/core/src/mesh.rs`) and mesh readers
(`crates/io`: PLY, VTK, Gmsh `.msh`, Medit `.mesh`). Grids, block models,
drillholes, variograms, kriging, and simulation stay in `boitata`.

## Stack

- **Rust workspace**: `crates/core`, `crates/io`, `crates/ops`,
  `crates/boolean`, `crates/query`, `crates/build`, and `crates/py`.
- **Exact geometry**: the `robust` crate for adaptive-precision orientation and
  in-circle predicates, paired with implicit intersection points for mesh
  arrangements.
- **Python bindings**: `pyo3` (`abi3-py311`) built with `maturin`.
- **Runtime dependencies**: `numpy` is the only required Python dependency.
  `matplotlib` is an optional extra (`pip install "wireframes[plot]"`).
- **Tooling and testing**: `mise` for task execution and toolchain pinning,
  `uv` for Python environments, `proptest` for property-based invariants,
  `criterion` for benchmarks, and `mkdocs` for documentation.

## First-class types and Python API

The Python package is imported as `import wireframes as wf` and follows `boitata`'s
container conventions:

- **`wf.Mesh`**: triangle surface mesh representing either an open surface
  (topography, fault, horizon, pit shell) or a closed solid (orebody, stope,
  lithological wireframe), carrying per-vertex and per-face attribute tables,
  `crs`, and `length_unit`.
- **`wf.Polylines`**: 3D strings (open paths and closed polygons) with per-line
  and per-vertex attributes for pit crests, toes, section interpretations, and
  survey breaklines.
- **`wf.PointSet`**: 3D point cloud with per-point attributes (collar surveys,
  LiDAR pickups, blast-hole pickups).
- **VTU cell data**: tetrahedra and hexahedra in `.vtu` files are read and
  written as cell arrays and attributes so finite-element and reservoir meshes
  round-trip without loss, without turning volume meshes into a first-class
  modeling type.

Every container exposes `.coords` as a `(n, 3)` `float64` NumPy array, `.x`,
`.y`, and `.z` 1D views, `.bounds`, `.column_names`, `to_table()`, and zero-copy
or thin conversions via `.to_polars()`, `.to_pandas()`, and `.to_pyarrow()`.
Plotting functions in `wf.plot` take `ax=None`, create a figure when `ax` is
omitted, and return `(fig, ax)`.

## Supported formats

| Format | Extension | Read | Write | Notes |
|---|---|---|---|---|
| Wavefront OBJ | `.obj` | M1 | M1 | Triangle faces and groups |
| Polygon File Format | `.ply` | M1 | M1 | ASCII and binary little/big endian, vertex and face properties |
| Object File Format | `.off` | M1 | M1 | Triangle faces |
| Stereolithography | `.stl` | M1 | M1 | ASCII and binary |
| AutoCAD DXF | `.dxf` | M2 | M2 | `3DFACE`, `MESH`, `POLYFACE`, `POLYLINE`, `LWPOLYLINE`, `POINT`, layers |
| VTK legacy | `.vtk` | M2 | M2 | `POLYDATA` and `UNSTRUCTURED_GRID` |
| VTK XML unstructured | `.vtu` | M2 | M2 | Triangle, quad, tetra, and hexa cells with point/cell data |
| VTK XML polydata | `.vtp` | M2 | M2 | Points, lines, and polygons (`Verts`, `Lines`, `Polys`, `Strips`) with point/cell data |
| Datamine wireframe CSV | `_tr.csv` / `_pt.csv` | M2 | M2 | Paired triangle (`TRISOL`/`PID1`/`PID2`/`PID3`) and point (`PID`/`XP`/`YP`/`ZP`) tables |

## Non-goals

- **Proprietary binary formats**: AutoCAD `.dwg` and binary Datamine `.dm` files
  are closed specifications and will not be read or written. Export DXF or
  `_tr.csv`/`_pt.csv` from the source package instead.
- **Interactive 3D viewers**: `wireframes` renders static publication and check plots
  through `matplotlib` (`wireframes[plot]`). Interactive WebGL or desktop viewers are
  out of scope.
- **Grids and block models**: regular, masked, and sub-blocked grids live in
  `boitata.BlockModel`; `wireframes` provides the geometric queries (`block_fraction`,
  `contains`, `signed_distance`) that `boitata` calls.

## Exact booleans and validation contract

Mine planning and resource estimation depend on conservation of volume. A silent
vertex shift or an unreported hole fill during a boolean operation changes
tonnage without warning. `wireframes` enforces three rules:

1. **Exact predicates and implicit points**: triangle-triangle intersections and
   arrangement cell classifications use exact adaptive-precision predicates
   (`robust::orient3d`, `robust::orient2d`) and exact rational or symbolic
   implicit coordinates for intersection vertices. Topological decisions never use
   floating-point epsilon comparisons; coordinates round to `f64` only when
   materializing the output `Mesh`.
2. **Refuse dirty solids by default**: boolean operations (`union`,
   `intersection`, `difference`, `clip_by_surface`) validate their inputs first.
   Every result is a solid: the exact result is rounded to `f64` with exact
   tests that collapse or snap what rounding cannot place, and where two inputs
   agree only to rounding their faces are first made to coincide exactly
   (`docs/design/rounding.md`). No vertex ends more than four ulps of its largest
   coordinate from its exact position, and `rounding_report=True` returns how
   far vertices moved and a proven bound on the volume change. A result these
   operations cannot certify raises `wf.KernelError` with its `MeshReport` and
   mesh (spec tolerance policy, item 5); no case known to the tests does. If an input has boundary
   edges, non-manifold edges or vertices, degenerate or duplicate faces,
   inconsistent winding, inward shells, or self-intersections, the function
   raises `wf.InvalidMesh` carrying the full `MeshReport` (problem counts and
   exact face/vertex indices).
3. **Explicit, audited repair**: nothing is modified silently. Calling
   `mesh.repair(...)` or passing `repair=True` to a boolean operation runs only
   the requested repair steps (welding duplicate vertices, dropping degenerate
   and duplicate faces, orienting shells, filling holes) and records every change
   in the returned report.
4. **Standing volume identity**: every boolean test suite and property test checks
   `vol(A union B) + vol(A intersect B) = vol(A) + vol(B)` within `f64`
   materialization rounding, alongside `Mesh.validate().summary.is_closed`.

## Roadmap

### M1: Core types, validation report, and OBJ/PLY/OFF/STL I/O

- **Scope**: initialize the Rust workspace (`crates/core`, `crates/io`,
  `crates/py`), `mise.toml`, `pyproject.toml`, and `python/wireframes/`. Implement
  `Mesh`, `Polylines`, `PointSet`, attribute tables, and `Mesh.validate()`
  producing `MeshReport` (`MeshSummary` and `MeshProblem` list covering
  degenerate faces, duplicate faces, duplicate vertices, boundary edges,
  non-manifold edges, non-manifold vertices, inconsistent winding, inward shells,
  and self-intersections). Implement ASCII and binary readers and writers for
  OBJ, PLY, OFF, and STL with round-trip property tests.
- **Acceptance check**:
  ```sh
  mise run verify && mise exec -- uv run --no-sync pytest python/tests/test_io_m1.py -k roundtrip
  ```
- **Status**: done. The acceptance check passes.

### M2: Mining and subsurface formats (DXF, VTK/VTU, Datamine CSV pairs)

- **Scope**: read and write DXF (`3DFACE`, polyface meshes, 3D polylines, points,
  and layer attributes), legacy VTK (`POLYDATA` and `UNSTRUCTURED_GRID`), XML VTU
  (surfaces plus data-only round-trip of tetrahedral and hexahedral cells), XML VTP
  polydata (meshes, strings, and points), and
  Datamine `_tr.csv`/`_pt.csv` wireframe pairs.
- **Acceptance check**:
  ```sh
  mise exec -- uv run --no-sync pytest python/tests/test_io_m2.py -k "dxf or vtk or vtu or datamine"
  ```
- **Status**: partial. Every format reads and writes, and the acceptance check
  passes. Still missing: a test against a real Datamine Studio export (spec
  decision 14). Only synthetic `_tr.csv`/`_pt.csv` fixtures exist so far.

### M3: Measure, validate, repair, and matplotlib plots

- **Scope**: exact signed and unsigned volume, surface area, face areas, face
  normals, bounding boxes, and area/volume centroids. Explicit `Mesh.repair()`
  (vertex welding within tolerance, degenerate/duplicate face removal, consistent
  winding and outward shell orientation, planar and minimal-area hole filling)
  returning `(repaired_mesh, repair_log)`. Implement `wf.plot.plan`,
  `wf.plot.section`, and `wf.plot.surface3d` on `matplotlib`.
- **Acceptance check**:
  ```sh
  mise exec -- uv run --no-sync pytest python/tests/test_measure_repair_plot.py
  ```
- **Status**: done. Hole filling offers `fill_holes="fan"` and `"min_area"`, and
  `repair()` returns `(mesh, RepairReport)`. The acceptance check passes.

### M4: Exact booleans, surface cutting, and plane slicing

- **Scope**: native Rust mesh arrangements with exact predicates and implicit
  intersection points. Implement `union`, `intersection`, `difference`, and
  `symmetric_difference` on closed solids; `clip_by_surface` to split or trim a
  solid against an open surface (for example, clipping an orebody or stope by
  topography or a pit shell); `slice_by_plane` producing closed/open `Polylines`
  for plan and section views; and `cut_and_fill` volume and thickness maps
  between two surfaces.
- **Acceptance check**:
  ```sh
  mise exec -- uv run --no-sync pytest python/tests/test_booleans.py -k "volume_identity and closed_output"
  ```
- **Status**: partial. Union, intersection, difference, symmetric difference,
  `clip_by_surface`, `split_by_plane`, `slice_by_plane` and `cut_fill` work, and
  the acceptance check passes. Solids may hold shells that touch along edges or
  at vertices, and a shell may touch itself the same way; booleans accept and
  return them. Known gaps:
  - A CI pipeline does not exist yet.
  - The semi-static filter in front of the exact predicates is not written.
  - Neighbouring solids that share a wall cannot sit in one valid mesh. Two
    shells that share all or part of a face, or part of an edge, or where a
    vertex of one rests inside a face of the other, are reported as coplanar
    overlaps or self-intersections, so `is_solid` is false and booleans raise
    `wf.InvalidMesh` (example 08, section 7). Stope files digitized against a
    common wall must be joined with `union` (one shell, stope attributes kept
    per face) or kept as separate meshes before any boolean. Accepting
    conforming face contacts would change the solid definition in the
    glossary and every boolean's classification; not started.
  - Rounding exact results to `f64` follows `docs/design/rounding.md` (method C
    with conforming of inputs that agree only to rounding). Every reduced case
    from the old refusals (pits and pushbacks clipped by one topography, prisms
    whose caps differ by one ulp, face-sharing stopes rotated to strike, the
    crosscut chain of example 07, the 648 inclined splits of `sphere(2)`) now
    gives solids within proven bounds, and the acceptance tests that were
    skipped or ignored for it run. So do the adversarial cases in
    `crates/boolean/tests/attacks.rs`: pits clipped by one topography and
    united, chains of lids and caps tilted by 1e-13 rad, repeated union then
    difference of one body, a stope split along its own face plane, a rotated
    stope block minus a raise, development turned to strike, and from a second
    round nine face-plane splits of grid stopes, repeated union then difference
    of turned bodies, stopes dipping 62.5 degrees minus a drive along strike or
    on their floor, and a sixteen-stope chain against a common wall minus a
    drive along strike. One case a property test found later still raises
    `wf.KernelError` (`uncertified`) and is ignored: stacked stopes at dip 89.6
    minus a drive on their floor. A third round of mining scenarios
    (`crates/boolean/tests/attacks_2.rs`) added thirteen ignored cases; five
    splits of stope blocks through their own vertices or a slot a few ulps
    wide now give solids, while two level splits of one block still raise
    `uncertified`, and a repeated union then
    difference of a raise in a worked orebody and the running union of
    pushbacks clipped by one ground raise `uncertified` or take minutes on a
    few thousand faces. Intersecting the block with a box whose top lies on
    the level plane certified in every case the sweep tried. The method has no general
    3D validity theorem, so `uncertified` stays possible for inputs the tests
    have not met.

### M5: Spatial queries for drillholes and block models

- **Scope**: BVH-accelerated exact `contains` (point-in-solid with boundary
  classification), `signed_distance` and closest point on surface, ray and
  desurveyed drillhole trace intersections returning entry/exit depths and pierced
  face attributes, and exact `block_fraction` (proportion of each axis-aligned or
  rotated 3D block inside a solid).
- **Acceptance check**:
  ```sh
  mise exec -- uv run --no-sync pytest python/tests/test_queries.py -k "block_fraction_sums_to_solid_volume"
  ```
- **Status**: done. `contains`, `locate`, `signed_distance`, `closest_point`,
  `ray_hits`, `segment_hits`, `trace_intervals` and `block_fraction` (with
  `rotation`) run on a BVH. The acceptance check passes.

### M6: Surface and solid construction and mesh editing

- **Scope**: 2.5D constrained Delaunay triangulation (`dtm`) from points and 3D
  breakline strings; solid construction from sequential section strings
  (`loft`) with branching and end-capping; polyline extrusion along a
  vector or path; quadric error decimation; volume-preserving Taubin smoothing;
  isotropic surface offset; and isotropic remeshing.
- **Acceptance check**:
  ```sh
  mise exec -- uv run --no-sync pytest python/tests/test_construct.py
  ```
- **Status**: partial. `dtm` with breaklines, `triangulate`, `extrude`, `loft`,
  `close_to_elevation`, `contours`, `decimate`, `smooth`, `offset` and `remesh`
  exist, and the acceptance check passes. Missing against the scope above:
  - `loft` links one ring to one ring. A section with several rings raises
    `wf.InvalidInput` (spec decision 15). Not started.
  - `extrude` sweeps along a vector only, not along a path. Not started; a
    drive along a curved centreline needs lofted rings.
  - `offset` moves vertices by `distance` along angle-weighted normals, so
    faces move less where they meet at an angle (by `distance / sqrt(3)` at a
    box corner). Neither a mitred nor an isotropic (Minkowski) offset exists
    yet (spec decision 16). Not started.
  - `smooth` runs Taubin smoothing, which changes a solid's volume slightly
    (by -0.29 percent for the orebody after 5 iterations in example 05).
    `preserve_volume=True` rescales the result about its volume centroid so
    the volume matches the input to rounding; it is off by default and needs
    a closed solid.

### M7: Boitata adopts wireframes

- **Scope**: publish `wireframes` wheels, add `wireframes` as a dependency of `boitata`,
  replace `boitata`'s `Mesh` (`crates/core/src/mesh.rs`) and mesh readers
  (`crates/io`: PLY, VTK, `.msh`, `.mesh`) with `wireframes`, and pass `boitata`'s
  entire verification and example suite.
- **Acceptance check**:
  ```sh
  mise -C ../boitata run verify   # boitata checked out next to wireframes
  ```
- **Status**: partial. The pre-release `0.1.0a1` is on PyPI: install it with
  `pip install --pre wireframes`. It has wheels for Linux x86_64 and aarch64
  (glibc 2.28 or newer), macOS (Intel and Apple silicon) and Windows x64, and
  pip builds from the sdist on other platforms. The release workflow builds and
  tests each wheel when a `v*` tag is pushed. Pull requests have no CI yet, and
  boitata does not depend on wireframes yet. `mise run compat` checks that a
  locally built wheel installs with pip, uv, pixi and conda.

## Performance

`mise run bench` runs the criterion benches. The table below gives criterion medians
measured on a 13th Gen Intel Core i9-13980HX (32 logical CPUs, wireframes on 32
threads). The `validate`, boolean and slicing rows are medians of three runs at commit
`cf79bfb` on 2026-11-20, interleaved with three runs of `ba8a660`, the commit before
the rounding work (1-minute load 7.9 to 16.5); the boolean rows include the `validate`
call each boolean runs on its result. The other rows come from commit `5d5ab30` on
2026-11-19. The [performance page](https://github.com/cratongen/wireframes/blob/main/docs/performance.md) lists the commands, versions,
confidence intervals, load, the before and after rows and the results for every kernel
and reader. Comparison with other libraries: not measured.

| Operation | Dataset | Median |
|---|---|---|
| STL binary write, read (in memory) | 100k-triangle torus | 3.29 ms, 6.94 ms |
| PLY binary write, read (in memory) | 100k-triangle torus | 0.720 ms, 2.86 ms |
| Mesh validation (`validate`) | 100k-triangle torus | 125 ms |
| Exact boolean (`union`) | Two 100k-triangle tori | 1.47 s (1.09 s before rounding) |
| Exact boolean (`difference`) | Two 100k-triangle tori | 1.51 s (1.17 s before rounding) |
| Plane slicing (`slice_by_plane`) | 100k-triangle torus across 50 sections | 86.3 ms |
| Point location (`locate_all`), index prebuilt | 1e5 points in a 100k-triangle torus | 4.68 ms |
| Signed distance, index prebuilt | 1e5 points, 100k-triangle torus | 56.9 ms |
| Exact block fraction (`block_fraction`) | 100 x 100 x 50 blocks, 100k-triangle torus | 127 ms |
| DTM (`dtm`) | 1e5 points | 61.6 ms |

## License

MIT
