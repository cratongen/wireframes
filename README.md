# drusa

`drusa` is a Rust-core, Python-exposed library to read, write, plot, measure,
and combine 3D solids, surfaces, and strings for mining, geology, and petroleum.

> A druse (*drusa* in Portuguese) is a rock cavity lined with crystal faces.

## Scope and relationship to boitata

`drusa` is a standalone geometry package with zero dependency on `boitata`.
Once `drusa` reaches Milestone 7, `boitata` will depend on `drusa` and delete
its internal mesh type (`crates/core/src/mesh.rs`) and mesh readers
(`crates/io`: PLY, VTK, Gmsh `.msh`, Medit `.mesh`). Grids, block models,
drillholes, variograms, kriging, and simulation stay in `boitata`.

## Stack

- **Rust workspace**: `crates/core`, `crates/io`, `crates/boolean`,
  `crates/query`, `crates/construct`, and `crates/py`.
- **Exact geometry**: the `robust` crate for adaptive-precision orientation and
  in-circle predicates, paired with implicit intersection points for mesh
  arrangements.
- **Python bindings**: `pyo3` (`abi3-py311`) built with `maturin`.
- **Runtime dependencies**: `numpy` is the only required Python dependency.
  `matplotlib` is an optional extra (`pip install "drusa[plot]"`).
- **Tooling and testing**: `mise` for task execution and toolchain pinning,
  `uv` for Python environments, `proptest` for property-based invariants,
  `criterion` for benchmarks, and `mkdocs` for documentation.

## First-class types and Python API

The Python package is imported as `import drusa as dr` and follows `boitata`'s
container conventions:

- **`dr.Mesh`**: triangle surface mesh representing either an open surface
  (topography, fault, horizon, pit shell) or a closed solid (orebody, stope,
  lithological wireframe), carrying per-vertex and per-face attribute tables,
  `crs`, and `length_unit`.
- **`dr.Polylines`**: 3D strings (open paths and closed polygons) with per-line
  and per-vertex attributes for pit crests, toes, section interpretations, and
  survey breaklines.
- **`dr.PointSet`**: 3D point cloud with per-point attributes (collar surveys,
  LiDAR pickups, blast-hole pickups).
- **VTU cell data**: tetrahedra and hexahedra in `.vtu` files are read and
  written as cell arrays and attributes so finite-element and reservoir meshes
  round-trip without loss, without turning volume meshes into a first-class
  modeling type.

Every container exposes `.coords` as a `(n, 3)` `float64` NumPy array, `.x`,
`.y`, and `.z` 1D views, `.bounds`, `.column_names`, `to_table()`, and zero-copy
or thin conversions via `.to_polars()`, `.to_pandas()`, and `.to_pyarrow()`.
Plotting functions in `dr.plot` take `ax=None`, create a figure when `ax` is
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
| Datamine wireframe CSV | `_tr.csv` / `_pt.csv` | M2 | M2 | Paired triangle (`TRISOL`/`PID1`/`PID2`/`PID3`) and point (`PID`/`XP`/`YP`/`ZP`) tables |

## Non-goals

- **Proprietary binary formats**: AutoCAD `.dwg` and binary Datamine `.dm` files
  are closed specifications and will not be read or written. Export DXF or
  `_tr.csv`/`_pt.csv` from the source package instead.
- **Interactive 3D viewers**: `drusa` renders static publication and check plots
  through `matplotlib` (`drusa[plot]`). Interactive WebGL or desktop viewers are
  out of scope.
- **Grids and block models**: regular, masked, and sub-blocked grids live in
  `boitata.BlockModel`; `drusa` provides the geometric queries (`block_fraction`,
  `contains`, `signed_distance`) that `boitata` calls.

## Exact booleans and validation contract

Mine planning and resource estimation depend on conservation of volume. A silent
vertex shift or an unreported hole fill during a boolean operation changes
tonnage without warning. `drusa` enforces three rules:

1. **Exact predicates and implicit points**: triangle-triangle intersections and
   arrangement cell classifications use exact adaptive-precision predicates
   (`robust::orient3d`, `robust::orient2d`) and exact rational or symbolic
   implicit coordinates for intersection vertices. Topological decisions never use
   floating-point epsilon comparisons; coordinates round to `f64` only when
   materializing the output `Mesh`.
2. **Refuse dirty solids by default**: boolean operations (`union`,
   `intersection`, `difference`, `clip_by_surface`) validate their inputs first.
   Valid closed inputs produce a valid closed output. If an input has boundary
   edges, non-manifold edges or vertices, degenerate or duplicate faces,
   inconsistent winding, inward shells, or self-intersections, the function
   raises `dr.InvalidMesh` carrying the full `MeshReport` (problem counts and
   exact face/vertex indices).
3. **Explicit, audited repair**: nothing is modified silently. Calling
   `mesh.repair(...)` or passing `repair=True` to a boolean operation runs only
   the requested repair steps (welding duplicate vertices, dropping degenerate
   and duplicate faces, orienting shells, filling holes) and records every change
   in the returned report.
4. **Standing volume identity**: every boolean test suite and property test checks
   `vol(A union B) + vol(A intersect B) = vol(A) + vol(B)` within `f64`
   materialization rounding, alongside `Mesh.validate().summary.is_closed`. The
   C++ `Manifold` library is used in CI as an external test oracle and benchmark
   baseline, never as a dependency of `drusa`.

## Roadmap

### M1: Core types, validation report, and OBJ/PLY/OFF/STL I/O

- **Scope**: initialize the Rust workspace (`crates/core`, `crates/io`,
  `crates/py`), `mise.toml`, `pyproject.toml`, and `python/drusa/`. Implement
  `Mesh`, `Polylines`, `PointSet`, attribute tables, and `Mesh.validate()`
  producing `MeshReport` (`MeshSummary` and `MeshProblem` list covering
  degenerate faces, duplicate faces, duplicate vertices, boundary edges,
  non-manifold edges, non-manifold vertices, inconsistent winding, inward shells,
  and self-intersections). Implement ASCII and binary readers and writers for
  OBJ, PLY, OFF, and STL with round-trip property tests.
- **Acceptance check**:
  ```sh
  mise run verify && mise exec -- pytest python/tests/test_io_m1.py -k roundtrip
  ```

### M2: Mining and subsurface formats (DXF, VTK/VTU, Datamine CSV pairs)

- **Scope**: read and write DXF (`3DFACE`, polyface meshes, 3D polylines, points,
  and layer attributes), legacy VTK (`POLYDATA` and `UNSTRUCTURED_GRID`), XML VTU
  (surfaces plus data-only round-trip of tetrahedral and hexahedral cells), and
  Datamine `_tr.csv`/`_pt.csv` wireframe pairs.
- **Acceptance check**:
  ```sh
  mise exec -- pytest python/tests/test_io_m2.py -k "dxf or vtk or vtu or datamine"
  ```

### M3: Measure, validate, repair, and matplotlib plots

- **Scope**: exact signed and unsigned volume, surface area, face areas, face
  normals, bounding boxes, and area/volume centroids. Explicit `Mesh.repair()`
  (vertex welding within tolerance, degenerate/duplicate face removal, consistent
  winding and outward shell orientation, planar and minimal-area hole filling)
  returning `(repaired_mesh, repair_log)`. Implement `dr.plot.plan`,
  `dr.plot.section`, and `dr.plot.surface3d` on `matplotlib`.
- **Acceptance check**:
  ```sh
  mise exec -- pytest python/tests/test_measure_repair_plot.py
  ```

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
  mise exec -- pytest python/tests/test_booleans.py -k "volume_identity and closed_output"
  ```

### M5: Spatial queries for drillholes and block models

- **Scope**: BVH-accelerated exact `contains` (point-in-solid with boundary
  classification), `signed_distance` and closest point on surface, ray and
  desurveyed drillhole trace intersections returning entry/exit depths and pierced
  face attributes, and exact `block_fraction` (proportion of each axis-aligned or
  rotated 3D block inside a solid).
- **Acceptance check**:
  ```sh
  mise exec -- pytest python/tests/test_queries.py -k "block_fraction_sums_to_solid_volume"
  ```

### M6: Surface and solid construction and mesh editing

- **Scope**: 2.5D constrained Delaunay triangulation (`dtm`) from points and 3D
  breakline strings; solid construction from sequential section strings
  (`link_strings`) with branching and end-capping; polyline extrusion along a
  vector or path; quadric error decimation; volume-preserving Taubin smoothing;
  isotropic surface offset; and isotropic remeshing.
- **Acceptance check**:
  ```sh
  mise exec -- pytest python/tests/test_construct.py
  ```

### M7: Boitata adopts drusa

- **Scope**: publish `drusa` wheels, add `drusa` as a dependency of `boitata`,
  replace `boitata`'s `Mesh` (`crates/core/src/mesh.rs`) and mesh readers
  (`crates/io`: PLY, VTK, `.msh`, `.mesh`) with `drusa`, and pass `boitata`'s
  entire verification and example suite.
- **Acceptance check**:
  ```sh
  mise -C /home/user/Projects/boitata run verify
  ```

## Performance

No performance numbers are stated before the code exists. Every number in this
table will come from reproducible benchmarks we run on the same machine and in
the same session using `criterion` (Rust) and interleaved multi-round timings
reporting medians (Python), with system load recorded, comparing `drusa` against
`trimesh`, `PyVista`/VTK, and `Manifold`:

| Operation | Dataset | drusa | trimesh | PyVista / VTK | Manifold |
|---|---|---|---|---|---|
| STL / PLY binary read & write | 500k-triangle pit shell | not measured | not measured | not measured | not measured |
| Mesh validation (`validate`) | 250k-triangle geological solid | not measured | not measured | not measured | not measured |
| Exact boolean (`union`, `difference`) | Two 50k-triangle intersecting solids | not measured | not measured | not measured | not measured |
| Plane slicing (`slice_by_plane`) | 100k-triangle solid across 50 sections | not measured | not measured | not measured | not measured |
| Point-in-solid (`contains`) | 1M block centroids in 50k-triangle solid | not measured | not measured | not measured | not measured |
| Exact block fraction (`block_fraction`) | 250k blocks intersecting a wireframe | not measured | not measured | not measured | not measured |

## License

MIT
