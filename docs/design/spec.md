# wireframes API contract

Status: binding contract for Milestones 1 to 6. Every parallel implementer follows
the names, signatures and semantics in this file. If the code needs a change here,
the lead edits this file first, then the code follows.

Sources: `AGENTS.md` and `README.md` (charter), a read-only survey of
`/home/user/Projects/boitata`, the Meshes.jl reference catalog, and format notes
for OBJ, PLY, OFF, STL, DXF, VTK legacy, VTK XML and Datamine CSV pairs.

## 0. Decisions where the sources disagree

Each line: the decision, then why.

1. Rust crates are named by this contract: `wireframes-core`, `wireframes-io`, `wireframes-ops`,
   `wireframes-boolean`, `wireframes-query`, `wireframes-build`, `wireframes-py`. The charter's
   `construct` crate becomes `wireframes-build`, and `measure`/`repair` move from core to
   `wireframes-ops`, so core stays a dependency-light kernel that every other crate needs.
2. Attributes use a wireframes-owned `Table` (typed columns with validity masks), not
   Arrow `RecordBatch`: it keeps numpy the only Python runtime dependency and keeps
   wireframes free of an Arrow version that boitata would then have to match. boitata
   converts with a small adapter in M7.
3. The Python `Table` drops boitata's `__arrow_c_*` capsule methods: they need an
   Arrow implementation. `to_pyarrow`, `to_polars` and `to_pandas` build from numpy
   and import their library lazily.
4. `Mesh` stays triangle-only (`u32` triangles). Readers fan-triangulate polygons
   from their first corner and record a face column `polygon`, so writers that
   support polygons (OBJ, PLY, OFF, VTK, VTU) rebuild the original polygons and
   round trips stay exact.
5. Volume cells (tetra, hexa, wedge, pyramid, polyhedron) round-trip as data in
   `wireframes_io::UnstructuredGrid`, which is not a modeling type (charter).
6. Boolean kernel: Cherchi et al. 2020 ("Fast and robust mesh arrangements using
   floating-point arithmetic") style implicit points (LPI, TPI) with a floating
   filter and Shewchuk expansion fallback, chosen over a bignum rational kernel
   (slower, needs a dependency) and over CGAL corefinement (C++).
7. Classification of arrangement patches uses exact ray parity plus label
   propagation across intersection curves, chosen over generalized winding numbers
   (floating sums cannot decide exactly) because inputs must be valid solids.
8. Point-in-solid casts a vertical ray (+z) and breaks ties by symbolic
   perturbation, chosen over an arbitrary ray direction because the axis ray reduces
   every crossing test to exact `orient2d` plus one `orient3d`.
9. Block fraction and bench volumes use analytic clamped-column integration
   (divergence theorem with the field `(0, 0, clamp(z) - z0)`), chosen over block
   discretization (boitata `proportion`) because it has no sampling error.
10. CDT for `dtm` and polygon triangulation uses the `spade` crate (pure Rust,
    exact predicates through `robust`, constraint splitting), chosen over writing a
    CDT, because a correct CDT is thousands of lines.
11. Self-intersection in validation is exact and includes touching between faces
    that share no vertex. boitata ignores touching, but touching solids break
    booleans, so wireframes reports it.
12. wireframes carries `crs` and `length_unit` as labels and checks them for equality.
    CRS transforms stay in boitata, which calls `Mesh::with_vertices`.
13. DXF is read and written by hand: a group-code reader (ASCII and binary) and an
    ASCII R2000 writer with shortest round-trip floats. No DXF crate is used; the
    `dxf` crate's writer neither guarantees round-trip float text nor writes
    polyface vertex flags correctly.
14. Datamine triangle id aliases are `TRIANGLE`, `TRIANG`, `TRI` and `TRISOL`
    (README says `TRISOL`, the survey says `TRIANGLE`; accept both, write
    `TRIANGLE`). The id is the first alias in that priority order that the header
    holds, whatever the column order; the other aliases are face attributes (a
    `TRISOL` surface id beside `TRIANGLE` is kept), so only `TRIANGLE` is a
    reserved face column name on write. Verify against a real export before M2
    closes.
15. `loft` (solid from sections) links one ring to one ring in v1. A section with
    several rings errors with `Unsupported`. Branching is deferred and listed in
    section 12.
16. `offset` moves vertices along angle-weighted normals and validates the output.
    It refuses an offset whose result self-intersects. A true Minkowski offset is
    deferred.
17. Gmsh `.msh` and Medit `.mesh` are not in this contract. The README says boitata
    will delete those readers in M7, so M7 must either keep them in boitata or add
    them here. Deferred, listed in section 12.
18. A valid solid may consist of several closed, outward-oriented, manifold shells
    that touch each other along edges or at vertices, as long as their interiors do
    not overlap (owner decision, 2026-11-19). Contact between shells is information,
    not a defect: `validate` counts it, booleans accept it and return it. A shell
    may also touch itself along edges or at vertices on the same terms (revised
    2026-11-19 after the underground acceptance round: refusing it made unions
    depend on their order, since an intermediate union of stopes closing around a
    gap touches itself until the gap is filled).

## 1. Workspace

```
wireframes/
  Cargo.toml                  workspace, members = ["crates/*"], resolver = "3", edition 2024
  rust-toolchain.toml         stable, pinned
  mise.toml                   tools and tasks (section 11)
  pyproject.toml              maturin, module-name = "wireframes._wireframes", python-source = "python"
  crates/
    core/      package wireframes-core
    io/        package wireframes-io
    ops/       package wireframes-ops
    boolean/   package wireframes-boolean
    query/     package wireframes-query
    build/     package wireframes-build
    py/        package wireframes-py (cdylib, lib name _wireframes)
  python/wireframes/
  python/tests/
  examples/
  docs/
```

Dependency graph (no cycles):

```
core  <- io
core  <- ops
core, ops <- query
core, ops <- boolean
core, ops, query, boolean <- build
all <- py
```

External dependencies (workspace `[workspace.dependencies]`, exact versions pinned
in `Cargo.lock`):

| crate | used by | why |
|---|---|---|
| `robust = "1"` | core | adaptive exact `orient2d`, `orient3d`, `incircle`, `insphere` |
| `rayon = "1"` | core, ops, query, boolean; io and build (dev) | data parallel loops; thread-count determinism tests |
| `thiserror = "2"` | core | `Error` derive |
| `flate2 = "1"` | io | zlib in VTU |
| `spade = "2"` | build | constrained Delaunay |
| `pyo3` (feature `abi3-py311`, `multiple-pymethods`), `numpy` | py | bindings, matching minor pair |
| dev: `proptest = "1"`, `num-bigint = "0.4"`, `num-rational = "0.4"`, `criterion = "0.8"` | all | property tests, exact oracle, benchmarks |

No other dependency without a line in this table.

## 2. Global rules

- Coordinates are `f64`, finite. Triangle indices are `u32` and in range.
  Constructors reject anything else with `Error::InvalidInput`.
- Orientation: a solid's triangles wind counter-clockwise seen from outside; face
  normal is `(b - a) x (c - a)`; signed volume is positive for an outward solid.
- **Solid**: a mesh whose `MeshSummary::is_solid()` is true (section 3.12): one or
  more closed, outward, manifold shells whose interiors do not overlap; shells may
  touch each other, and a shell may touch itself, at exactly coincident edges and
  vertices.
- Determinism: results never depend on thread count. Parallel maps collect in input
  order. Parallel sums split work into fixed chunks of 4096 items (independent of
  thread count), sum each chunk with Neumaier compensated summation, then sum the
  chunk partials sequentially in chunk order. Ties in priority queues break by the
  smaller index.
- Memory: no fixed caps. Algorithms size buffers from the input.
- Errors: every fallible function returns `wireframes_core::Result<T>`. Messages are
  lowercase and name the offending argument, for example
  `"triangles must index existing vertices: face 12 uses 9 of 8"`.
- Immutability: every operation returns a new object; inputs are never modified.
  No function repairs its input unless the caller passes `repair = true` or calls
  `repair` directly.
- CRS and unit: a combining operation (boolean, cut and fill, thickness, any query
  that takes two geometries or points with labels) calls `units::check_same`. An
  undeclared value matches anything; two declared values that differ after
  normalization give `Error::CrsMismatch` or `Error::UnitMismatch`. The result takes
  the declared value.

### Tolerance policy

1. Topological and side-of decisions are exact: `orient2d`, `orient3d`,
   `incircle`, and implicit-point predicates. No epsilon decides orientation, side
   of plane, inside or outside, crossing or degeneracy.
2. User tolerances exist only where the user asks for a geometric merge: repair
   welding (`weld`), validation duplicate-vertex detection (`tolerance`), `dtm`
   duplicate removal (`min_spacing`), polyline joining, decimation error bound,
   remesh target length. Each is absolute, in the data's length unit, finite and
   `>= 0`; `tolerance::check(name, value)` returns `Error::InvalidInput` otherwise.
   Zero means exact equality (bitwise after mapping `-0.0` to `+0.0`).
3. Measures (area, volume, centroid, fractions, distances) are `f64` results.
   Before accumulation, implementations subtract a reference point (the mesh's
   bounding-box minimum, or the query box minimum) so UTM-sized coordinates do not
   lose digits.
4. Tests compare measures with `close(a, b) := |a - b| <= 1e-9 * max(1, |a|, |b|)`
   unless a test states a different bound and the reason.
5. Materialized intersection points round each exact coordinate to the nearest
   `f64`. Repeated and collinear corners that rounding makes are removed per
   output shell (section 6.2 item 7), and so is a sliver that rounding folds onto
   the coplanar neighbour across its longest edge, such as the kink where a stope
   wall cuts a drive wall quad that is planar only to rounding, a sliver that
   rounding turns away from its own plane, and a sliver that rounding tilts into a
   face of the other input when its neighbour comes from the same input face.
   Rounding can still leave a self-intersection that the exact arrangement did not
   have between faces of different inputs. This is not rare: it
   happens whenever two input surfaces agree only to an ulp, so the exact result
   holds slivers thinner than `f64` can place. Stopes tessellated on different grids
   and rotated to strike are one case: their shared walls are planar only to
   rounding, so the exact union keeps pockets an ulp thick between them. Open-pit work hits it often: booleans
   of a pit and its pushback after each was clipped by the topography (their shared
   walls and caps were computed separately and differ by ulps), and
   `clip_by_surface` of polygonal pits with many wall quads by a surface whose
   faces cross the wall quads near their diagonals (the open-pit acceptance tests
   measured 127 invalid clips in 200 seeded tilted planes for an octagonal pit at
   UTM coordinates). Prefer operations on the exact inputs: clip the difference
   (`(s2 - s1).clip_by_surface(topo)`) rather than taking the difference of clipped
   pieces. The result is closed and manifold, and its volume is right to rounding,
   but `validate` reports the self-intersections, so `is_solid` is false; a fix
   needs snap rounding of the arrangement, which is not implemented. Output tests check
   combinatorial closure and manifoldness, and `validate` on outputs is run in the
   property tests with the failure rate reported, never hidden.

## 3. wireframes-core

Files: `crates/core/src/{lib.rs, error.rs, table.rs, mesh.rs, polylines.rs,
points.rs, units.rs, tolerance.rs, predicates.rs, exact.rs, intersect.rs, aabb.rs,
bvh.rs, topology.rs, validate.rs}`.

`lib.rs` re-exports: `Error, Result, Table, Column, ColumnData, Mesh, Polylines,
PointSet, Aabb, Bvh, EdgeMap, HalfEdges, MeshProblemKind, MeshProblem,
MeshSummary, MeshReport, validate`, and the modules `predicates, exact, intersect,
units, tolerance`.

### 3.1 error.rs

```rust
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0}")]
    InvalidInput(String),
    #[error("invalid mesh: {context}: {defects}; the error's report lists each problem ..., and mesh.repair() or repair=True fixes what repair can")]
    InvalidMesh { context: String, report: Box<MeshReport> },
    #[error("{}", format_location(path, line, message))]
    Format { path: Option<std::path::PathBuf>, line: Option<u64>, message: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Unsupported(String),
    #[error("crs mismatch: {left} vs {right}")]
    CrsMismatch { left: String, right: String },
    #[error("length unit mismatch: {left} vs {right}")]
    UnitMismatch { left: String, right: String },
    #[error("missing column {name:?}; columns are {available:?}")]
    MissingColumn { name: String, available: Vec<String> },
}
pub type Result<T> = std::result::Result<T, Error>;
```

`InvalidMesh` is the only way a boolean, clip, query or build function reports an
input that fails validation; `context` names the argument (`"a"`, `"b"`,
`"solid"`, `"surface"`). The message names the cause and the way out: `{defects}`
is `MeshSummary::defects()`, the nonzero defect counts (`boundary_edges=18,
self_intersections=29`, or `not closed`), followed by the pointer to the report
and to repair, so `wf.union(topography, orebody)` says why it refuses.

### 3.2 table.rs

```rust
#[derive(Clone, Debug, PartialEq)]
pub enum ColumnData { Float64(Vec<f64>), Int64(Vec<i64>), Bool(Vec<bool>), Utf8(Vec<String>) }

#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    pub data: ColumnData,
    pub valid: Option<Vec<bool>>,   // None means every row is valid
    pub unit: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct Table { /* names: Vec<String>, columns: Vec<Column>, rows: usize */ }
```

`Column`: `len()`, `is_valid(i) -> bool`, `null_count()`, `type_name() -> &'static str`
(`"float64" | "int64" | "bool" | "utf8"`), `nulls(kind, n) -> Column` (all-null column
of a type), constructors `Column::f64(Vec<f64>)`, `Column::i64`, `Column::bool`,
`Column::utf8`. Null payload values are `NaN`, `0`, `false`, `""`.

`Table`:
- `new(rows: usize) -> Table` (no columns), `rows()`, `num_columns()`,
  `names() -> &[String]`, `column(name) -> Option<&Column>`,
  `require(name) -> Result<&Column>` (`MissingColumn`).
- `with_column(self, name, Column) -> Result<Table>`: add or replace in place of the
  old position; length must equal `rows`.
- `without_column(self, name) -> Table`.
- `take(&self, rows: &[u32]) -> Table`: gather; indices must be `< rows` (panics are
  bugs; callers pass valid indices).
- `take_or_null(&self, rows: &[Option<u32>]) -> Table`: `None` gives a null row.
- `filter(&self, keep: &[bool]) -> Result<Table>`.
- `concat(a: &Table, b: &Table) -> Result<Table>`: rows of `a` then `b`; columns are
  the union in first-seen order (`a` first); a column missing on one side is null
  there; the same name with different types gives `InvalidInput`.

Row-count invariants: vertex tables have one row per vertex, face tables one per
triangle, polyline feature tables one per feature, polyline vertex tables one per
vertex, point tables one per point.

### 3.3 mesh.rs

```rust
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    vertex_attributes: Table,
    face_attributes: Table,
    pub crs: Option<String>,
    pub length_unit: Option<String>,
}
```

- `Mesh::new(vertices: Vec<[f64;3]>, triangles: Vec<[u32;3]>) -> Result<Mesh>`:
  checks finiteness and index range; attributes start empty with correct row counts.
- `Mesh::empty() -> Mesh`.
- `vertices() -> &[[f64;3]]`, `triangles() -> &[[u32;3]]`, `n_vertices()`,
  `n_faces()`, `is_empty()` (no faces), `corners(f: u32) -> [[f64;3];3]`.
- `vertex_attributes() -> &Table`, `face_attributes() -> &Table`.
- `with_vertex_column(self, name, Column) -> Result<Mesh>`,
  `with_face_column(self, name, Column) -> Result<Mesh>`,
  `with_attributes(self, vertex: Table, face: Table) -> Result<Mesh>` (row counts
  checked).
- `with_vertices(self, Vec<[f64;3]>) -> Result<Mesh>`: same count, finite; topology
  and attributes kept (used by boitata for CRS and local grid transforms).
- `with_crs(self, Option<String>) -> Mesh`, `with_length_unit(self, Option<String>) -> Mesh`.
- `bounds() -> Option<Aabb>` over used and unused vertices; `None` when no vertices.
- `face_aabbs() -> Vec<Aabb>`.
- `filter_faces(&self, keep: &[bool]) -> Result<Mesh>`: keeps listed faces and their
  attributes, drops vertices no kept face uses, keeps vertex order.
- `concat(meshes: &[&Mesh]) -> Result<Mesh>`: appends vertices and faces in order,
  offsets indices, `Table::concat` on attributes, CRS and unit checked with
  `units::check_same`; no welding.
- `validate(&self, tolerance: f64) -> Result<MeshReport>` (calls `validate::validate`).
- `is_closed(&self) -> bool`: at least one non-degenerate face and every edge of the
  non-degenerate faces used by exactly two faces, or by more that `Shells` pairs into
  shells in contact (combinatorial plus exact radial order; equals
  `MeshSummary::is_closed`).

### 3.4 polylines.rs

```rust
#[derive(Clone, Debug, PartialEq)]
pub struct Polylines {
    vertices: Vec<[f64; 3]>,
    parts: Vec<u32>,      // offsets into vertices, len = n_parts + 1, starts at 0
    features: Vec<u32>,   // offsets into parts, len = n_features + 1, starts at 0
    closed: Vec<bool>,    // one per part
    vertex_attributes: Table,
    attributes: Table,    // one row per feature
    pub crs: Option<String>,
    pub length_unit: Option<String>,
}
```

Rules (boitata compatible): a closed part does not repeat its first vertex; a closed
part has `>= 3` vertices, an open part `>= 1` (one vertex is a point, for example
the pinch-out of a `loft`); a feature has `>= 0` parts; holes and
multipart features are several parts of one feature; plan inside is even-odd over
the feature's closed parts.

- `Polylines::new(vertices, parts, features, closed, attributes: Table) -> Result<Polylines>`;
  `with_vertex_attributes(self, Table) -> Result<Polylines>`.
- `from_parts(parts: Vec<Vec<[f64;3]>>, closed: Vec<bool>, feature_of_part: Option<Vec<u32>>) -> Result<Polylines>`
  (`feature_of_part` non-decreasing from 0; `None` means one feature per part).
- `len()` (features), `is_empty()`, `n_parts()`, `n_vertices()`, `vertices()`,
  `part(p) -> &[[f64;3]]`, `part_range(p) -> Range<usize>`,
  `feature_parts(f) -> Range<usize>`, `closed() -> &[bool]`, `rings(f) -> Vec<usize>`.
- `attributes()`, `vertex_attributes()`, `with_column(self, name, Column)` (feature),
  `with_vertex_column(self, name, Column)`, `with_attributes(self, Table)`.
- `filter(&self, keep: &[bool]) -> Result<Polylines>` (whole features).
- `bounds() -> Option<Aabb>`.
- `concat(items: &[&Polylines]) -> Result<Polylines>`.
- `length() -> Vec<f64>` (3D per feature, includes closing segments),
  `area() -> Vec<f64>` (plan, even-odd, per feature), `feature_contains(f, [f64;2]) -> bool`
  (plan, even-odd; a point on a ring edge counts as inside; exact via `orient2d`),
  `locate([f64;2]) -> Option<usize>`,
  `distance([f64;2], feature: Option<usize>, signed: bool) -> Option<f64>`
  (plan distance to the nearest ring edge; negative inside when `signed`).
- `with_crs`, `with_length_unit`, `with_vertices` like `Mesh`.

### 3.5 points.rs

```rust
#[derive(Clone, Debug, PartialEq)]
pub struct PointSet {
    coords: Vec<[f64; 3]>,
    attributes: Table,
    pub crs: Option<String>,
    pub length_unit: Option<String>,
}
```

`new(coords) -> Result`, `with_attributes`, `with_column`, `coords()`, `attributes()`,
`len()`, `is_empty()`, `bounds()`, `filter(&[bool])`, `with_crs`, `with_length_unit`,
`with_vertices`.

### 3.6 units.rs, tolerance.rs

- `units::normalize_crs(&str) -> String`: trims; `epsg:4326` and WKT ending in
  `AUTHORITY["EPSG","4326"]` or `ID["EPSG",4326]` become `EPSG:4326`; anything else
  is kept verbatim.
- `units::length_factor(unit: &str) -> Result<f64>`: metres per unit for `m, km, cm,
  mm, ft, us_ft, in, yd` (aliases `metre, meter, foot, feet, us_survey_foot`);
  unknown gives `InvalidInput`.
- `units::check_same(what: &str, a: (&Option<String>, &Option<String>), b: (...)) -> Result<(Option<String>, Option<String>)>`:
  returns the merged `(crs, length_unit)`.
- `tolerance::check(name: &str, value: f64) -> Result<f64>`.
- `tolerance::canonical(x: f64) -> f64`: maps `-0.0` to `+0.0` (used by every exact
  weld and hash).

### 3.7 predicates.rs

Wrappers over `robust` returning `f64` with the exact sign (magnitude is only
meaningful for its sign). `robust` stays exact while every coordinate is zero or
has a magnitude in `[2^-120, 2^120]`; outside that range the wrappers evaluate the
same determinant with `Expansion`, so the sign is exact for any finite input:

```rust
pub fn orient2d(a: [f64;2], b: [f64;2], c: [f64;2]) -> f64   // > 0 when c is left of a->b
pub fn orient3d(a: [f64;3], b: [f64;3], c: [f64;3], d: [f64;3]) -> f64 // > 0 when d is below plane abc (robust convention, documented)
pub fn incircle(a: [f64;2], b: [f64;2], c: [f64;2], d: [f64;2]) -> f64
pub fn insphere(a: [f64;3], b: [f64;3], c: [f64;3], d: [f64;3], e: [f64;3]) -> f64
pub fn orient2d_axis(a: [f64;3], b: [f64;3], c: [f64;3], drop: usize) -> f64 // projection dropping axis `drop`
pub fn side_of_plane(p: [f64;3], a: [f64;3], b: [f64;3], c: [f64;3]) -> Sign  // Positive = side the normal (b-a)x(c-a) points to
pub enum Sign { Negative, Zero, Positive }
pub fn collinear3(a: [f64;3], b: [f64;3], c: [f64;3]) -> bool // all three orient2d_axis are zero
```

`side_of_plane` hides the `robust` sign convention; every other module uses it, not
raw `orient3d`.

### 3.8 exact.rs

Exact arithmetic for implicit points (decision 6):

```rust
pub struct Expansion(/* Vec<f64> nonoverlapping, increasing magnitude; or a signed integer times 2^e */);
impl Expansion {
    pub fn from_f64(x: f64) -> Self;
    pub fn sum(&self, o: &Self) -> Self; pub fn diff(&self, o: &Self) -> Self;
    pub fn product(&self, o: &Self) -> Self; pub fn scale(&self, b: f64) -> Self;
    pub fn sign(&self) -> Sign;
    pub fn estimate(&self) -> f64;               // nearest-ish f64 with the exact sign, clamped to the finite range
    pub fn components(&self) -> Option<&[f64]>;  // None for the integer form
    pub fn to_dyadic(&self) -> (bool, Vec<u64>, i64); // exact value: ±limbs * 2^exp
}
pub fn quotient_estimate(n: &Expansion, d: &Expansion) -> f64; // n / d within a few ulps at any magnitude
pub fn round_quotient(n: &Expansion, d: &Expansion) -> f64;    // nearest f64, ties to even, ±inf past f64::MAX
pub fn shell_volume_sign(v: &[[f64;3]], faces: &[[u32;3]]) -> Sign; // exact sign of a closed shell's volume

#[derive(Clone, Debug)]
pub enum ImplicitPoint {
    Explicit([f64; 3]),
    /// Line through p, q meets plane through r, s, t.
    Lpi { p: [f64;3], q: [f64;3], r: [f64;3], s: [f64;3], t: [f64;3] },
    /// Three planes (each through three explicit points) meet.
    Tpi { planes: [[[f64;3];3];3] },
}
impl ImplicitPoint {
    pub fn to_f64(&self) -> [f64; 3];            // each coordinate rounded to nearest
    pub fn is_valid(&self) -> bool;              // denominator nonzero (exact)
}
pub fn orient3d_i(a: &ImplicitPoint, b: [f64;3], c: [f64;3], d: [f64;3]) -> Sign;
pub fn orient2d_i(a: &ImplicitPoint, b: &ImplicitPoint, c: &ImplicitPoint, drop: usize) -> Sign;
pub fn cmp_axis(a: &ImplicitPoint, b: &ImplicitPoint, axis: usize) -> std::cmp::Ordering;
pub fn same_point(a: &ImplicitPoint, b: &ImplicitPoint) -> bool;
pub fn plane_side(p: [f64;3], origin: [f64;3], normal: [f64;3]) -> Sign; // exact sign of (p - origin) . normal
```

Every `Expansion` operation is exact for any finite input. An operation whose
`f64` components would underflow (a product below `2^-968`) or overflow redoes
itself on a signed integer times a power of two, and its result keeps that form.
Inputs in the usual range never leave the `f64` form.

Each predicate first evaluates with a semi-static interval filter; on an uncertain
sign it evaluates in homogeneous coordinates with `Expansion`. `to_f64` rounds the
exact quotient to nearest (`estimate` of numerator and denominator expansions, then
one correction step verified with an expansion comparison).

### 3.9 intersect.rs

```rust
pub enum TriTri { Disjoint, Touch, Cross, CoplanarOverlap }
pub fn tri_tri(a: [[f64;3];3], b: [[f64;3];3]) -> TriTri;
pub enum SegTri { Miss, Interior, Edge(u8), Vertex(u8), Coplanar }
pub fn seg_tri(p: [f64;3], q: [f64;3], t: [[f64;3];3]) -> (SegTri, Option<ImplicitPoint>);
pub fn point_in_triangle(p: [f64;3], t: [[f64;3];3]) -> bool; // closed triangle, exact, coplanar required
```

All decisions use `predicates` and `exact`. `Touch` means the closed triangles meet
but their relative interiors do not cross; `Cross` means an edge of one meets the
other's relative interior or the interiors cross; `CoplanarOverlap` means coplanar
with a 2D overlap of positive area.

### 3.10 aabb.rs, bvh.rs

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb { pub min: [f64; 3], pub max: [f64; 3] }
impl Aabb {
    pub fn empty() -> Aabb; pub fn from_points(&[[f64;3]]) -> Aabb;
    pub fn union(&self, o: &Aabb) -> Aabb; pub fn intersects(&self, o: &Aabb) -> bool; // closed boxes
    pub fn contains(&self, p: [f64;3]) -> bool; pub fn center(&self) -> [f64;3];
    pub fn distance2(&self, p: [f64;3]) -> f64; pub fn expand(&self, d: f64) -> Aabb;
    pub fn is_empty(&self) -> bool;
}

pub struct Bvh { /* nodes, leaf item order */ }
impl Bvh {
    pub fn new(boxes: &[Aabb]) -> Bvh; // median split on the longest axis of centroids, leaf size 4, ties by index
    pub fn len(&self) -> usize;
    pub fn query(&self, b: &Aabb, f: impl FnMut(u32));
    pub fn pairs(&self, other: &Bvh, f: impl FnMut(u32, u32));
    pub fn self_pairs(&self, f: impl FnMut(u32, u32)); // i < j
    pub fn ray(&self, origin: [f64;3], dir: [f64;3], t_max: f64, f: impl FnMut(u32));
    pub fn nearest(&self, p: [f64;3], dist2: impl FnMut(u32) -> f64) -> Option<(u32, f64)>; // smaller index on ties
}
```

Traversal order is fixed by the tree, so callbacks run in a deterministic order.
Box tests use floating comparisons and only prune; exact predicates make every final
decision.

### 3.11 topology.rs

```rust
pub struct EdgeMap {
    pub edges: Vec<[u32; 2]>,          // sorted endpoint pair, edges ordered by (min, max)
    pub face_edges: Vec<[u32; 3]>,     // edge id of local edge e = corner e -> (e+1)%3
    pub edge_face_offsets: Vec<u32>,   // CSR into edge_faces, len = edges + 1
    pub edge_faces: Vec<u32>,          // faces of each edge, ascending
    pub vertex_face_offsets: Vec<u32>,
    pub vertex_faces: Vec<u32>,
}
impl EdgeMap {
    pub fn new(mesh: &Mesh, skip: Option<&[bool]>) -> EdgeMap; // skip marks faces to leave out (degenerates)
    pub fn faces_of_edge(&self, e: u32) -> &[u32];
    pub fn faces_of_vertex(&self, v: u32) -> &[u32];
    pub fn boundary_edges(&self) -> Vec<u32>;
    pub fn non_manifold_edges(&self) -> Vec<u32>;
}

pub struct HalfEdges { pub twin: Vec<u32> } // half-edge h = 3*f + e; twin u32::MAX on boundary
impl HalfEdges {
    pub fn new(mesh: &Mesh) -> Result<HalfEdges>; // InvalidInput if an edge has more than 2 faces
    pub fn next(h: u32) -> u32; pub fn prev(h: u32) -> u32;
    pub fn face(h: u32) -> u32; pub fn origin(&self, mesh: &Mesh, h: u32) -> u32;
    pub fn boundary_loops(&self, mesh: &Mesh) -> Vec<Vec<u32>>; // vertex loops, deterministic start at smallest half-edge
    pub fn components(&self, mesh: &Mesh) -> Vec<u32>;          // edge-connected component label per face, labels by first face
}

pub enum Radial { Paired(Vec<[usize; 2]>), Coincident, Unbalanced }
pub fn radial_pairs(forward: &[bool], turn: impl Fn(usize, usize) -> Sign, same: impl Fn(usize, usize) -> bool) -> Radial;
pub fn radial_pairs_at(mesh: &Mesh, faces: &[u32], p: [f64; 3], q: [f64; 3]) -> Radial;

pub struct Shells { pub of: Vec<u32>, pub count: usize, pub non_manifold: Vec<u32>, pub contacts: Vec<u32>, /* private pairs */ }
impl Shells {
    pub fn new(mesh: &Mesh, edges: &EdgeMap) -> Shells;
    pub fn links(&self, edges: &EdgeMap, id: u32) -> Vec<[u32; 2]>; // face pairs joined across edge id
}
```

Radial pairing: the faces around one edge `p -> q` sort by exact angle about the
edge direction (`turn(a, b)` is the sign of `det[q - p, r_a - p, r_b - p]` for third
corners `r`, `same` decides coplanar half-planes by `orient2d_axis`). A face running
`q -> p` bounds a wedge of solid with the next face counter-clockwise, which must run
`p -> q`. Two faces in one half-plane give `Coincident`; directions that do not
alternate give `Unbalanced`. `Shells` joins faces across two-face edges and across
radial pairs at edges with more faces; an edge whose faces do not pair, or whose
pairs meet one shell twice, is in `non_manifold`, the rest of the multi-face edges
are `contacts`. Shell labels are numbered by first face.

### 3.12 validate.rs

Port of boitata plus two kinds and exact geometry. Declaration order is the sort
order, and boitata's kinds keep their relative order.

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MeshProblemKind {
    DegenerateFace, DuplicateFace, DuplicateVertex, BoundaryEdge, NonManifoldEdge,
    NonManifoldVertex, InconsistentWinding, InwardShell, SelfIntersection,
    CoplanarOverlap, UnusedVertex,
}
impl MeshProblemKind { pub fn name(&self) -> &'static str; } // snake_case: ..., "coplanar_overlap", "unused_vertex"

#[derive(Clone, Debug, PartialEq)]
pub struct MeshProblem {
    pub kind: MeshProblemKind,
    pub face: Option<u32>,
    pub vertex: Option<u32>,
    pub edge: Option<u8>,    // local corner e -> (e+1)%3 of `face`
    pub other: Option<u32>,  // earlier duplicate, other face of an edge, earlier face of a pair
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct MeshSummary {
    pub degenerate_faces: usize, pub duplicate_faces: usize, pub duplicate_vertices: usize,
    pub boundary_edges: usize, pub non_manifold_edges: usize, pub non_manifold_vertices: usize,
    pub inconsistent_edges: usize, pub shells: usize, pub inward_shells: usize,
    pub self_intersections: usize, pub coplanar_overlaps: usize, pub unused_vertices: usize,
    pub touching_shell_edges: usize, pub touching_shell_vertices: usize, // information, not defects
    pub is_closed: bool,
}
impl MeshSummary {
    pub fn is_solid(&self) -> bool; // is_closed and every count except shells, unused_vertices and the touching counts is 0
    pub fn to_pairs(&self) -> Vec<(&'static str, i64)>; // stable order, used by Python dict
}

#[derive(Clone, Debug, PartialEq)]
pub struct MeshReport { pub summary: MeshSummary, pub problems: Vec<MeshProblem> }
pub fn validate(mesh: &Mesh, tolerance: f64) -> Result<MeshReport>;
```

Semantics:
- DegenerateFace: a repeated corner index, or the three corners exactly collinear
  (`collinear3`). boitata's `AREA_EPS` is dropped (tolerance policy rule 1).
- DuplicateFace: same corner set as an earlier face, any order; `other` = earlier.
- DuplicateVertex: within `tolerance` (Euclidean) of an earlier vertex; at
  `tolerance = 0` exact equality. Grid hash with cell `max(tolerance, tiny)` and
  neighbour cells; `other` = the earliest such vertex. Exempt: an earlier vertex
  with exactly equal coordinates when the shells of both are all closed (a contact
  between shells, or of a shell with itself).
- UnusedVertex: no face uses the vertex.
- BoundaryEdge (1 face), NonManifoldEdge (3 or more faces, except a contact in
  `Shells::contacts`), NonManifoldVertex (faces around the vertex form more than one
  fan, joined through `Shells::links`, unless every fan lies in a closed shell),
  InconsistentWinding (two faces run the shared edge in the same direction).
- Shells: `Shells::new` on the non-degenerate faces. A shell is closed when none of
  its faces has a boundary, non-manifold or inconsistent edge.
- Contacts of closed shells: a contact point is an exact coordinate where more than
  one vertex umbrella (fan) of closed shells meets, by shared index or coincident
  copies. Every edge location with more than two faces of closed shells that has a
  contact point at one end, or is in `Shells::contacts`, is a candidate:
  `radial_pairs_at` on all its faces must pair each face with a face it is linked to
  across that edge (`Shells::links`); then the edge counts once in
  `touching_shell_edges`, else the first unlinked pair (or, when `Unbalanced`, the
  first two faces) is a `SelfIntersection` (the interiors overlap there).
  `touching_shell_vertices` counts the contact points and the ends of the counted
  edges. A shell touching itself along an edge can have ends where its faces form
  one umbrella; that edge is then shared by index. `Coincident` adds nothing: the coplanar overlap or duplicate face is
  reported already. Contacts are only at exactly coincident vertices and edges; a
  vertex of one shell on the interior of another shell's edge or face stays a
  `SelfIntersection`.
- InwardShell: per closed, consistently wound shell, nesting depth = number of other
  closed shells containing it (`topology::shell_location`: the first of its vertices
  off the other shell, exact point-in-shell by the section 7.1 rule; when every
  vertex lies on the other shell, as for a cavity touching it at all its vertices,
  the first exact point inside one of its edges that is off it). Even depth must
  have positive signed volume, odd depth negative, by the exact sign of
  `exact::shell_volume_sign` (an `f64` sum with a proven error bound, exact
  arithmetic when the bound does not settle the sign). `repair` orients by the same
  depth.
- SelfIntersection: a pair of non-degenerate faces with `tri_tri` = `Cross` or
  `Touch` where the faces share no vertex, or share exactly one vertex and meet
  elsewhere too, or share an edge and are not coplanar but meet elsewhere (cannot
  happen for planar triangles, kept for completeness). Candidates from
  `Bvh::self_pairs`.
- CoplanarOverlap: a pair with `tri_tri` = `CoplanarOverlap`, any shared vertices.
- Degenerate faces are left out of edge, vertex, shell and intersection checks.
- `problems` sorted by kind, then face, then vertex, then other.
- `tolerance` invalid gives `InvalidInput`.

## 4. wireframes-io

Files: `crates/io/src/{lib.rs, common.rs, grid.rs, obj.rs, ply.rs, off.rs, stl.rs,
dxf.rs, vtk.rs, vtu.rs, datamine_csv.rs}`; fixtures in `crates/io/tests/fixtures/`
(each under 1 MB, synthetic).

### 4.1 lib.rs

```rust
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Data {
    pub mesh: Option<Mesh>,
    pub polylines: Option<Polylines>,
    pub points: Option<PointSet>,
    pub grid: Option<UnstructuredGrid>,
    pub length_unit: Option<String>,  // from the file (DXF $INSUNITS), never applied as a rescale
    pub warnings: Vec<String>,        // e.g. "skipped 3 INSERT entities"
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format { Obj, Ply, Off, Stl, Dxf, Vtk, Vtu, Vtp, DatamineCsv }

#[derive(Clone, Debug, Default)]
pub struct ReadOptions { pub length_unit: Option<String>, pub crs: Option<String> } // labels set on every returned object

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum DxfEntity { #[default] Face3d, Polyface }

#[derive(Clone, Debug, Default)]
pub struct WriteOptions {
    pub ascii: bool,                 // OBJ/OFF/DXF/Datamine are always text
    pub dxf_entity: DxfEntity,
    pub zlib: bool,                  // VTU only
    pub skip_unsupported_columns: bool, // else a column the format cannot store errors
}

pub fn detect(path: &Path) -> Result<Format>; // by extension, case-insensitive; *_tr.csv / *_pt.csv
pub fn read(path: impl AsRef<Path>) -> Result<Data>;
pub fn read_with(path: impl AsRef<Path>, opts: &ReadOptions) -> Result<Data>;
pub fn write(path: impl AsRef<Path>, data: &Data) -> Result<()>;
pub fn write_with(path: impl AsRef<Path>, data: &Data, opts: &WriteOptions) -> Result<()>;
pub fn read_mesh(path: impl AsRef<Path>) -> Result<Mesh>;           // InvalidInput "file holds no mesh"
pub fn write_mesh(path: impl AsRef<Path>, mesh: &Mesh, opts: &WriteOptions) -> Result<()>;
pub fn read_polylines(path: impl AsRef<Path>) -> Result<Polylines>;
pub fn read_points(path: impl AsRef<Path>) -> Result<PointSet>;
```

Each module exposes `pub fn read(r: &[u8]) -> Result<Data>` and
`pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>>`, except
`datamine_csv` which takes and returns the pair `(pt: &[u8], tr: &[u8])`. `lib.rs`
does all file access, so modules are pure and testable on bytes. Writing a part of
`Data` that a format cannot hold (polylines to STL, a grid to OBJ) gives
`Unsupported`.

`grid.rs`:

```rust
#[derive(Clone, Debug, PartialEq)]
pub struct UnstructuredGrid {
    pub points: Vec<[f64; 3]>,
    pub cell_types: Vec<u8>,          // VTK cell type ids, verbatim
    pub offsets: Vec<i64>,            // len = cells + 1, starts at 0
    pub connectivity: Vec<i64>,
    pub faces: Option<(Vec<i64>, Vec<i64>)>, // polyhedron face stream (offsets, ids), verbatim
    pub point_data: Table,
    pub cell_data: Table,
}
impl UnstructuredGrid {
    pub fn surface(&self) -> Result<Mesh>; // cells of types 5..9 (strips alternate winding, pixel 0,1,3,2), fan from first corner
    pub fn n_cells(&self) -> usize;
}
```

`common.rs`: shortest round-trip float formatting (`format!("{}", x)` / `{:?}` for
`f64`), a whitespace tokenizer with line numbers, the exact `weld` (bitwise after
`canonical`, first-appearance order), the polygon fan with `polygon` column, the
binary scalar reader (ASCII token, LE, BE) shared by PLY, OFF and VTK, and the
cp1252 fallback decoder.

### 4.2 Format contract

Rows list what each format reads into `Data` and writes back. "Exact" means write
then read returns equal `Data` (vertex `f64` bits, triangle indices, attribute
values and nulls, part and feature structure) for a canonical input as defined in
4.3.

| module | extensions | reads | writes | attributes | exact round trip |
|---|---|---|---|---|---|
| obj | .obj | `v` (first 3 values), `f` any size (first index field, 1-based and negative relative at parse time, forward references allowed), `o`, `g` (several names joined by one space), `l` as polylines, `p` as points, `\` continuation, CRLF; ignores `vt vn s mtllib usemtl` | `v`, `o`/`g` on change, `f` 1-based, polygons rebuilt from `polygon` | face `object`, `group` (Utf8); faces before the first `g` get `"default"` only if some `g` exists | yes |
| ply | .ply | ascii, binary LE/BE; all 8 scalar types and aliases; `vertex_indices` or `vertex_index`; unknown elements read through; `tristrips` with -1 restart; int64/uint64 refused; a binary body error gives its byte offset from the start of the file | `double x y z`, `list uchar int vertex_indices` (`list int int` above 255 corners), ascii or binary LE | numeric vertex and face properties; integer types read as Int64, float types as Float64; Float64 NaN reads as null; nulls write as NaN (Float64) or the column writes as double if Int64 has nulls; Bool writes as uchar and reads back as Int64; Utf8 unsupported | yes for Float64 without NaN and Int64 within i32 |
| off | .off | `[ST][C][N][4][n]OFF`, counts on keyword line or next, a file that ends early errors with the line it ends on, comments, `4OFF` divide by w, `nOFF` dim 3 only, binary OFF refused with `Unsupported`; face colours ignored | `OFF`, then polygons | none | yes (vertices and faces) |
| stl | .stl | binary first (`len >= 84 && 84 + 50*count == len`), fallback rules from the format notes with a warning, then ASCII tokens case-insensitive, several solids, each `solid` closed by `endsolid` before the next `solid` or the end of the file (else a truncation error); a file starting with `solid` that fails as ASCII but declares fewer triangles than its length holds reads as padded binary with the trailing-bytes warning | binary header `wireframes` padded, or ASCII `solid wireframes` | face `solid` (Utf8, ASCII only), face `stl_attribute` (Int64, binary, only when any value is nonzero) | only for f32-representable vertices (`stl::f32_rounding(mesh)` counts the vertices binary STL moves; Python `write` warns with `UserWarning` when any move); corners are exact-welded on read |
| dxf | .dxf | a file must end with the `0`/`EOF` group (else a truncation error naming the last line read); a leading UTF-8 byte order mark is skipped; an entity error names the line of the group that holds the fault; a float group that is not finite errors; paper space entities (group 67 = 1) are skipped with a count in `warnings`; `3DFACE` (every corner needs its x, y and z groups; tri, or quad fanned with `polygon`), polyface `POLYLINE` (flag 64), polygon mesh (flag 16, M x N with wrap flags, fitted vertices when 75 set), `MESH` level-0 cage (hand-parsed from ASCII; a smoothing level (91) above 0 adds a `warnings` entry that only the cage is read; binary DXF reports skipped MESH in `warnings`), 3D `POLYLINE`, 2D `POLYLINE` and `LWPOLYLINE` (OCS to WCS by the arbitrary axis algorithm, bulge ignored with a warning when nonzero), `LINE`, `POINT`; `INSERT` skipped with a count in `warnings`; `$INSUNITS` to `length_unit` (1 in, 2 ft, 4 mm, 5 cm, 6 m, 7 km, else none); cp1252 fallback before AC1021 | hand-written ASCII R2000 with HEADER `$INSUNITS`, LAYER table for every layer used, faces as 3DFACE (default) or POLYFACE per layer split at 32767 vertices or faces; polylines as 3D POLYLINE (flag 8, closed flag 1); points as POINT | face, feature and point `layer` (Utf8, default `"0"`), `color` (Int64 ACI, absent when BYLAYER everywhere) | yes for mesh (canonical), polylines and points |
| vtk | .vtk | legacy 2.0 to 5.1, ASCII and BINARY (big-endian), POLYDATA (VERTICES, LINES, POLYGONS, TRIANGLE_STRIPS in that cell order) and UNSTRUCTURED_GRID; SCALARS (with lookup table skipped), VECTORS, NORMALS, TENSORS, FIELD; `METADATA` skipped; `%20` names | version 5.1 OFFSETS/CONNECTIVITY vtktypeint64, BINARY unless `ascii`; POLYDATA for Mesh, Polylines and PointSet; UNSTRUCTURED_GRID for `grid` | numeric point and cell data; component arrays as `<name>_0..`; nulls as `<name>_valid` UInt8 companion | yes |
| vtu | .vtu, .vtp | XML; ascii, inline base64, appended raw or base64, zlib, UInt32/UInt64 headers, BigEndian; several Pieces concatenated with id offsets; vtp cell order Verts, Lines, Polys, Strips (fixes boitata bug 1); integer arrays read as i64/u64 (fixes bug 2); polyhedron faces (both legacy and VTK 9.4 names) kept in `faces`; LZ4/LZMA refused; an array whose data ends early errors with the file byte offset of its data; a file not ending in `</VTKFile>` is truncated and errors | appended raw, UInt64 headers, LittleEndian; zlib when `opts.zlib`; `.vtu` writes a mesh as triangle cells or the `grid`; `.vtp` writes mesh, polylines, points | as vtk | yes |
| datamine_csv | `<name>_tr.csv` + `<name>_pt.csv` | pair found from either name; delimiter sniffed (`,` or `;`, with `,` decimal under `;`); headers trimmed, unquoted, case-insensitive; point `XP YP ZP PID`; triangle id the first of `TRIANGLE`/`TRIANG`/`TRI`/`TRISOL` present (by that priority, not column order; the rest are face attributes), `PID1 PID2 PID3`; a format error names the file of the pair that holds it (`pit_pt.csv:3: ...`); ids parsed as f64 then required integral; `-1.0E30` is null; PID map via HashMap (gaps allowed); missing PID errors with the row; duplicate PID errors; a non-numeric type/unit row after the header skipped | `XP,YP,ZP,PID` (PID 1-based, vertex order) and `TRIANGLE,PID1,PID2,PID3` (1-based), then attributes; `,` and `.` | extra point columns as vertex attributes, extra triangle columns as face attributes; numeric to Float64, else Utf8; `SID` kept as a face column | yes |

### 4.3 Canonical meshes for round-trip tests

A mesh is canonical for a format when: (STL, DXF 3DFACE) no two vertices have equal
coordinates, no vertex is unused, and vertices appear in first-use order over
faces; (STL) every coordinate is f32-representable; (all) column types are ones the
format stores. `wireframes_io::canonicalize(mesh: &Mesh, format: Format) -> Mesh` exists
in `common.rs` as a test helper (`pub` under `#[doc(hidden)]`) so property tests
generate any mesh and compare against its canonical form.

## 5. wireframes-ops

Files: `crates/ops/src/{lib.rs, measure.rs, slices.rs, thickness.rs, repair.rs,
components.rs, transform.rs}`.

### 5.1 measure.rs

```rust
pub fn area(mesh: &Mesh) -> f64;
pub fn face_areas(mesh: &Mesh) -> Vec<f64>;
pub fn face_normals(mesh: &Mesh) -> Vec<[f64; 3]>;      // unit; [0,0,0] for degenerate faces
pub fn face_centroids(mesh: &Mesh) -> Vec<[f64; 3]>;
pub fn vertex_normals(mesh: &Mesh) -> Vec<[f64; 3]>;    // angle-weighted, unit; [0,0,0] for isolated vertices
pub fn dip(mesh: &Mesh) -> Vec<(f64, f64)>;             // per face (dip degrees 0..90, dip direction degrees 0..360 from +y clockwise); NaN dip direction for horizontal faces
pub fn volume(mesh: &Mesh) -> Result<f64>;              // signed; InvalidMesh unless closed and consistently wound (EdgeMap check, not full validate; an edge with more than two faces passes when `Shells` pairs it as a contact)
pub fn centroid(mesh: &Mesh) -> Result<[f64; 3]>;       // volume centroid, same precondition, InvalidInput when volume is 0
pub fn area_centroid(mesh: &Mesh) -> Option<[f64; 3]>;  // None when area is 0
```

### 5.2 slices.rs

```rust
/// Volume of the solid inside an axis-aligned box, by clamped-column integration.
pub fn box_volume(mesh: &Mesh, faces: &[u32], b: &Aabb) -> f64;
pub fn volume_by_elevation(mesh: &Mesh, edges: &[f64]) -> Result<Vec<f64>>;
```

`box_volume` clips each listed face's plan projection to the box's xy rectangle
(Sutherland-Hodgman on the 3D triangle against the four vertical planes), splits the
pieces at `z = b.min.z` and `z = b.max.z`, and integrates
`sign(n_z) * (clamp(z, zmin, zmax) - zmin)` over each piece's plan area. Coordinates
are shifted by `b.min` first. Callers pass every face whose plan bounding box
meets the box's xy rectangle and whose max z is above `b.min.z`. Faces entirely
above the box count: each contributes `zmax - zmin` over its clipped footprint. A
face entirely below contributes 0 and may be left out. Precondition: closed,
consistently outward mesh (unchecked here; public callers check).

`volume_by_elevation`: `edges` strictly increasing, at least 2; result `i` is the
solid volume with `edges[i] <= z < edges[i+1]`; InvalidMesh unless closed and
consistently wound.

### 5.3 thickness.rs

```rust
pub enum ThicknessMode { Vertical, Normal }
pub fn surface_z(surface: &Mesh, xy: &[[f64; 2]]) -> Vec<f64>; // highest surface z at (x,y); NaN where no face covers; vertical faces ignored; on-edge counts as covered (exact orient2d)
pub fn vertical_distance(surface: &Mesh, points: &[[f64; 3]]) -> Vec<f64>; // z - surface_z
pub fn thickness(top: &Mesh, bottom: &Mesh, xy: &[[f64; 2]], mode: ThicknessMode) -> Result<Vec<f64>>;
```

`Vertical`: `surface_z(top) - surface_z(bottom)`, negative where the surfaces cross,
NaN where either is missing. `Normal`: unsigned distance from the top point
`(x, y, z_top)` to the closest point of `bottom`, signed like `Vertical`. CRS and
unit checked.

### 5.4 repair.rs

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HoleFill { Fan, MinArea }

#[derive(Clone, Debug, PartialEq)]
pub struct RepairOptions {
    pub weld: Option<f64>,          // None: no weld; Some(0.0): exact equal coordinates
    pub remove_degenerate: bool,    // default true
    pub remove_duplicate_faces: bool, // default true
    pub remove_unused_vertices: bool, // default true
    pub orient: bool,               // default true: consistent winding, outward by nesting parity
    pub fill_holes: Option<HoleFill>, // default None
    pub weld_touching_shells: bool, // default false: closed shells in contact keep their own vertices
}
impl Default for RepairOptions { /* weld Some(0.0), the four bools true, fill_holes None, weld_touching_shells false */ }

#[derive(Clone, Debug, PartialEq)]
pub enum RepairAction {
    VertexMerged { vertex: u32, into: u32 },
    VertexRemoved { vertex: u32 },
    VertexSplit { vertex: u32, copy: u32 },            // input vertex, output index of its extra copy
    FaceRemoved { face: u32, reason: MeshProblemKind }, // DegenerateFace or DuplicateFace
    FaceFlipped { face: u32 },
    FaceAdded { face: u32, hole: u32 },                 // face index in the output
}

#[derive(Clone, Debug, PartialEq)]
pub struct RepairReport {
    pub actions: Vec<RepairAction>,     // in step order: weld (merges, then splits), degenerate, duplicate, unused, orient, fill
    pub vertex_map: Vec<Option<u32>>,   // input vertex -> output vertex (first copy of a split vertex)
    pub face_map: Vec<Option<u32>>,     // input face -> output face
    pub before: MeshSummary,
    pub after: MeshSummary,
    pub holes_left: usize,              // loops not filled (pinched at a vertex used by 2 boundary edges)
}

pub fn repair(mesh: &Mesh, opts: &RepairOptions) -> Result<(Mesh, RepairReport)>;
```

Semantics: steps run in the listed order. Weld keeps the earliest vertex of each
cluster and its attributes; clusters are transitive closures within `weld`. Unless
`weld_touching_shells`, weld then keeps closed shells in contact apart: `Shells`
(section 3.11) on the welded mesh, after degenerate and duplicate removal, finds
the closed shells (no boundary, non-manifold or inconsistent edge), and a cluster
used by several closed shells becomes one output vertex per shell (open faces share
one), each with the attributes of the earliest input vertex that shell used there.
The first copy keeps the cluster's earliest vertex, as without contacts. An input
vertex shared by several closed shells gets a `VertexSplit` per extra copy. Without
`weld`, nothing is merged or split. Splitting needs a consistent winding around the
contact, as `radial_pairs` does; a contact it cannot pair stays welded and is
reported as non-manifold.
Degenerate and duplicate removal keeps the earliest face. Orient flips the fewest
faces per edge-connected component to agree with that component's first face, then
flips whole closed components to match nesting parity (exact containment);
non-orientable components are left as they were and counted in `after`. `Fan` fills
each loop from its first vertex; `MinArea` uses the O(n^3) minimum-area
triangulation (Liepa) and falls back to `Fan` above 256 loop vertices with a
recorded action per added face. New faces go last with null face attributes and the
neighbouring winding. Self-intersections are never repaired. `repair` is
idempotent: repairing its own output gives zero actions.

### 5.5 components.rs

```rust
pub fn component_labels(mesh: &Mesh) -> Vec<u32>;       // `Shells` of the non-degenerate faces (touching shells are separate), degenerate faces join their edge neighbours; label = order of first face
pub fn components(mesh: &Mesh) -> Result<Vec<Mesh>>;    // one mesh per label, faces and vertices in input order
```

### 5.6 transform.rs

```rust
pub fn transform(mesh: &Mesh, m: [[f64; 4]; 4]) -> Result<Mesh>; // affine, row-major, applied to column vectors; flips faces when det < 0
pub fn translate(mesh: &Mesh, d: [f64; 3]) -> Mesh;
pub fn to_length_unit(mesh: &Mesh, unit: &str) -> Result<Mesh>; // InvalidInput when no unit declared
pub fn flip(mesh: &Mesh) -> Mesh;
```

## 6. wireframes-boolean

Files: `crates/boolean/src/{lib.rs, arrangement.rs, triangulate.rs, classify.rs,
ops.rs, clip.rs, slice.rs, cutfill.rs}`.

### 6.1 Types and functions

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op { Union, Intersection, Difference, SymmetricDifference }

#[derive(Clone, Debug, Default)]
pub struct BooleanOptions { pub repair: bool } // true: repair each input with RepairOptions::default() first

#[derive(Clone, Debug)]
pub struct BooleanResult {
    pub mesh: Mesh,
    pub face_source: Vec<(u8, u32)>,   // per output face: (0 = a, 1 = b, input face index after repair)
    pub repairs: [Option<RepairReport>; 2],
}

pub fn boolean(a: &Mesh, b: &Mesh, op: Op, opts: &BooleanOptions) -> Result<BooleanResult>;
pub fn union(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult>;
pub fn intersection(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult>;
pub fn difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult>;
pub fn symmetric_difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult>;

pub struct ClipResult { pub below: Mesh, pub above: Mesh, pub repairs: [Option<RepairReport>; 2] }
pub fn clip_by_surface(solid: &Mesh, surface: &Mesh, opts: &BooleanOptions) -> Result<ClipResult>;
pub fn split_by_plane(mesh: &Mesh, origin: [f64; 3], normal: [f64; 3]) -> Result<(Mesh, Mesh)>;
pub fn slice_by_plane(mesh: &Mesh, origin: [f64; 3], normal: [f64; 3], offsets: &[f64]) -> Result<Polylines>;

pub struct CutFill { pub cut: f64, pub fill: f64, pub net: f64, pub cut_area: f64, pub fill_area: f64, pub area: f64 }
pub fn cut_fill(base: &Mesh, design: &Mesh, boundary: Option<&Polylines>) -> Result<CutFill>;
```

### 6.2 Semantics

Input checks (booleans and `clip_by_surface`): CRS and unit via `check_same`.
Without `repair`, each solid input must pass `validate(0.0)` with
`summary.is_solid()`, else `InvalidMesh { context: "a" | "b" | "solid", report }`.
A mesh with no faces (the empty result of a boolean) is the empty solid and is
accepted wherever a solid is: by booleans, `clip_by_surface` (both parts empty),
`volume` (0), `volume_by_elevation` (zeros) and the solid queries (every point
outside, every block fraction 0). `validate` still reports it as not closed, so
`is_solid` stays false for it.
With `repair`, the input is repaired first, then validated; failing still errors.
The `surface` of `clip_by_surface` must have no degenerate, duplicate, non-manifold,
inconsistent, self-intersecting or coplanar-overlap problems (boundary edges allowed).

Arrangement (`arrangement.rs`):
1. Candidate pairs from `Bvh::pairs` on face boxes.
2. `tri_tri` per pair, producing intersection segments whose endpoints are
   `ImplicitPoint`s (explicit input vertices, LPI of an edge of one triangle with the
   other's plane, TPI for coplanar configurations).
3. Points on a face are merged with `exact::same_point`, never by distance.
4. `triangulate.rs`: per face, constrained triangulation of its corners, points
   and segments in the face's dominant-axis projection, using `orient2d_i` only.
5. `classify.rs`: patches are maximal face regions bounded by intersection curves.
   A patch with an input vertex not on the other mesh takes that vertex's exact
   classification (section 7.1). Other patches take the label of a neighbour across
   an intersection edge, flipped or kept by the exact side test of the neighbour's
   opposite vertex against the crossing face plane.
6. Keep rules (B faces flipped when kept by `Difference`):
   - Union: A outside B, B outside A. Intersection: A inside B, B inside A.
   - Difference: A outside B, B inside A (flipped).
   - SymmetricDifference: A outside B and B outside A as they are, A inside B and
     B inside A flipped; every coplanar patch dropped.
   - Coplanar same-facing patches: keep A's once for Union and Intersection, drop
     both for Difference and SymmetricDifference. Coplanar opposite-facing: drop both
     for Union, Intersection and SymmetricDifference, keep A's once for Difference.
   - Shells (`ops::shells`): kept sub-triangles join across edges with two of them;
     at an edge with more, `radial_pairs` orders them exactly about the edge
     (`turn` is the sign of `det[d, m_a, m_b]` and `same` the sign of `m_a . m_b`,
     with `d = q - p` from the implicit endpoints and `m` the exact input-face normal,
     negated when the face is flipped or runs `q -> p`) and joins the pairs. Each
     shell gets its own vertex copies, so shells in contact share coordinates and no
     indices.
7. Output vertices: input vertices keep their `f64`; implicit points round with
   `to_f64`. Topology comes from the exact arrangement, then each output shell is
   cleaned of what rounding alone breaks, without moving any point: vertices of
   one shell whose rounded coordinates are equal become one vertex (keeping an
   input vertex's attribute row when one of them has it), faces that then repeat a
   vertex are dropped, opposite faces on the same three vertices cancel, and a face
   whose rounded corners are exactly collinear (`collinear3`) is removed by
   splitting the face across its longest edge at its middle corner (the new faces
   keep that face's source). A face whose middle corner rounds to strictly inside the
   coplanar face across its longest edge (exact `orient3d` and `orient2d_axis` on the
   rounded points) is a fold and is removed the same way. Each such split flips the
   longest edge to the edge from the middle corner to the neighbour's far corner; it
   is skipped where that edge already exists or where the longest edge was itself
   made by a split, so every split removes an edge the shell came in with and the
   cleanup ends (flips on rounded points can otherwise cycle without bound). A face
   that rounding turned away from the exact plane it lies in (the sign of
   `n . ((b - a) x (c - a))` on the rounded corners, exact, with `n` the input
   face's normal rounded to `f64` and negated when the face is flipped, is not
   positive) is split away the same way, whether or not its neighbour is coplanar.
   Where the split would make an edge that already exists, a flat or turned face
   loses its shortest edge instead: one end merges into the other (an input vertex
   is kept) when both are different exact points, the link condition holds, and no
   face moved by the merge is turned or flat. Last, `validate` runs on each shell;
   the faces of its degenerate, self-intersection and coplanar-overlap problems and
   their edge neighbours are split across their longest edge where the neighbour
   there comes from the same input face and both new faces face its plane. These
   rounds repeat while the shell's problem count drops; a round that does not lower
   it is undone. Merges and these splits move the surface by at most the width of a
   face that rounding turned, folded or flattened. Every step keeps each edge used
   once in each direction. Equal coordinates in different shells stay apart
   (touching shells), and so do the copies of one exact point within a shell.
   Rounding can still leave self-intersections between faces of different inputs
   (tolerance policy rule 5).
8. Attributes: output face attributes are `Table::concat` of A's and B's face
   tables taken by `face_source`; vertex attributes keep input rows, null for new
   vertices. An empty result (nothing kept) is `Mesh::empty()` with the merged
   schema. Inputs may be solids with touching shells. Result is a solid when inputs
   are solids (exact topology), and a valid input to further booleans: parts that
   touch along an edge or at a vertex, including the `a - b` and `b - a` parts of
   `SymmetricDifference` along the intersection curve, are separate shells. A result
   shell that touches itself gets one vertex copy per umbrella of faces about each
   exact point, the umbrellas found through the same radial pairs; where it touches
   itself along an edge whose ends it does not touch, that edge keeps its four faces
   on shared vertices. `clip_by_surface` uses the same split per piece. A result that
   is not closed and manifold after rounding, or that keeps a face whose rounded
   corners are collinear (which `validate` leaves out, so its edges are open), errors
   with `InvalidInput("internal: rounding the result to f64 left it open or
   non-manifold")` (`"... the clip result ..."` for `clip_by_surface`), never
   returns it.

`clip_by_surface`: `above` is the part of the solid on the side the surface normals
point to, `below` the other; both are closed and outward, capped by the pieces of
the surface inside the solid (cap faces carry the surface's face attributes; the
`below` cap is flipped). The surface boundary must lie outside the solid: if any
boundary edge of the surface has a point strictly inside the solid,
`InvalidInput("surface does not cut through the solid")`. A surface that misses
the solid returns the whole solid on the side given by the exact classification of
one solid vertex, and an empty mesh on the other. A part of the solid that the
surface does not cut is classified by the nearest surface crossing on the line
through one of its points along z, then x, then y (x and y only when the point is
an input vertex). When none of those lines meets the surface, the part lies
beyond the surface's extent and has no side: `InvalidInput("surface does not reach
part of the solid ...")`, never a guess.

`split_by_plane`: returns `(negative, positive)` where positive means
`(p - origin) . normal > 0`; points on the plane go to positive (boitata rule).
Crossing points on edges are computed exactly and rounded once, before the cap is
built; a crossing that rounds exactly onto an end of its edge is that input vertex.
Section points of a closed mesh whose rounded coordinates project onto the same
point along the dominant axis of the normal (shells touching on the plane, or
exact points that round within a few ulps of each other, such as the crossings
next to a vertex that lies an ulp off an inclined plane) are one point of the cap
triangulation, represented by an input vertex when one of them is; each output
shell then gets its own copy of such a point, so touching shells stay apart in
both parts. Each part then gets the cleanup of section 6.2 item 7 (equal
coordinates merge, repeated and cancelling faces go, collinear corners and folds
are split away), repeated until nothing changes. A part of a closed mesh that is
still not closed and manifold errors with `InvalidInput("internal: rounding the
split to f64 left it open or non-manifold")`. Side tests use
`exact::plane_side`, which evaluates the sign of `(p - origin) . normal` exactly with
expansions, so `origin` and `normal` are taken as given f64 values. Closed meshes are
capped: the cap is a CDT of the rounded section points in the dominant-axis
projection, with exact predicates on those `f64` points, so the cap triangles written
out never overlap or flip (a CDT of the exact points can hold slivers that rounding
flips). Open meshes are trimmed. Cap faces get null attributes. `slice_by_plane` uses
the same side rule. Where the mesh touches the plane from one side (a vertex or an
edge in the plane with every incident face on the positive side), the exact positive
part touches itself there and is not a solid; the cap does not pass through that
contact as a vertex or edge, so `validate` reports self-intersections on that part.

`slice_by_plane`: boitata `section` semantics: one feature per line per offset;
feature columns `offset` (Float64) and the face attributes of the faces it crosses;
a line splits where face attributes change; closed meshes give closed rings wound
counter-clockwise around the inside seen from the normal side (clockwise around
cavities); open surfaces give open lines; a vertex on the plane counts on the normal
side. `offsets` empty means `[0.0]`. Error when a face column is named `offset`.

`cut_fill`: both surfaces must be height fields (no two faces with positive plan
area overlap in plan, decided exactly), else `InvalidInput` naming the surface and a
face. The integration domain is the plan overlap of both surfaces intersected with
the even-odd region of `boundary` (feature 0 and all its rings when one feature,
all features otherwise). The plan overlay of both triangulations (and boundary
edges) is built with exact `orient2d` on implicit 2D points; on each overlay cell
`d = z_design - z_base` is linear, cells are split on `d = 0`; `cut` integrates
`-d` where `d < 0`, `fill` integrates `d` where `d > 0`, `net = fill - cut`, areas are
the plan areas of each sign and of the whole domain.

## 7. wireframes-query

Files: `crates/query/src/{lib.rs, solid.rs, distance.rs, ray.rs, trace.rs,
fraction.rs}`.

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i8)]
pub enum Location { Inside = -1, Boundary = 0, Outside = 1 } // sign agrees with signed distance

pub struct SolidIndex { /* mesh: Mesh, bvh: Bvh, plan bvh */ }
impl SolidIndex {
    pub fn new(solid: &Mesh) -> Result<SolidIndex>; // validate(0.0) and is_solid, else InvalidMesh { context: "solid" }
    pub fn mesh(&self) -> &Mesh;
    pub fn locate(&self, p: [f64; 3]) -> Location;
    pub fn locate_all(&self, points: &[[f64; 3]]) -> Vec<Location>;       // rayon, input order
    pub fn contains_all(&self, points: &[[f64; 3]], boundary_inside: bool) -> Vec<bool>;
    pub fn signed_distance(&self, points: &[[f64; 3]]) -> Vec<f64>;       // negative inside, 0 on boundary
    pub fn block_fraction(&self, blocks: &Blocks) -> Result<Vec<f64>>;
    pub fn trace_intervals(&self, trace: &Trace) -> Result<Vec<Interval>>;
}

pub struct SurfaceIndex { /* mesh, bvh */ }
impl SurfaceIndex {
    pub fn new(mesh: &Mesh) -> SurfaceIndex; // any mesh
    pub fn closest(&self, points: &[[f64; 3]]) -> Vec<Closest>;
    pub fn distance(&self, points: &[[f64; 3]]) -> Vec<f64>;
    pub fn ray_hits(&self, origins: &[[f64; 3]], dirs: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>>;
    pub fn segment_hits(&self, a: &[[f64; 3]], b: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>>;
}

pub struct Closest { pub point: [f64; 3], pub face: u32, pub distance: f64 }
pub struct Hit { pub t: f64, pub face: u32, pub point: [f64; 3], pub kind: HitKind, pub entering: Option<bool> }
pub enum HitKind { Interior, Edge, Vertex, Coplanar }

pub struct Blocks {
    pub centroids: Vec<[f64; 3]>,
    pub size: [f64; 3],                    // positive, finite
    pub rotation: Option<[[f64; 3]; 3]>,   // columns are the block axes in world; orthonormal, det +1 (checked to 1e-12)
}

pub struct Trace { pub xyz: Vec<[f64; 3]>, pub depth: Vec<f64> } // desurveyed vertices and their measured depths, non-decreasing
pub struct Interval { pub from: f64, pub to: f64, pub inside: bool, pub entry_face: Option<u32>, pub exit_face: Option<u32> }
```

### 7.1 Semantics

- `locate`: boundary first: `Boundary` when the point lies in a closed face
  (`orient3d == 0` and `point_in_triangle`). Otherwise cast the ray `+z` from `p`; a
  face crosses when `(x, y)` lies in its plan projection and the face is above `p`
  (exact `side_of_plane`). Plan-edge and plan-vertex ties use symbolic perturbation
  of the query point by `(e, e^2)` in `(x, y)`: `orient2d(a, b, p) == 0` is
  replaced by the sign of the perturbed determinant (first `b.y - a.y` term sign,
  then `a.x - b.x` term sign), so each crossing is counted once. Vertical faces
  never cross. Odd count is `Inside`. Each point costs one BVH ray walk over a plan
  BVH.
- Boundary policy: `contains_all(.., boundary_inside = true)` treats the solid as a
  closed set (boitata `contains` behaviour for points off the surface);
  `false` treats it as open.
- `signed_distance`: magnitude from `SurfaceIndex::closest`, sign from `locate`.
- `closest`: exact nearest face by BVH with ties broken by the smaller face index;
  point computed with Ericson's region test in f64.
- `ray_hits`: `dirs` nonzero; hits with `t >= 0` sorted by `t` then face; a ray
  crossing an edge reports one hit per distinct exact point (duplicates across faces
  sharing the edge collapse to the smallest face); `entering` is `Some` for closed
  meshes (sign of `dir . normal`, exact by `side_of_plane`), else `None`.
  `segment_hits` same with `0 <= t <= 1`.
- `trace_intervals`: intervals cover `[depth[0], depth[last]]` in order, merge
  adjacent equal labels, and alternate `inside`. Crossing parameters come from exact
  segment-face tests; depths interpolate linearly along each trace segment. A run
  lying in a face counts as inside (closed set). Zero-length intervals are dropped.
- `block_fraction`: per block, `ops::slices::box_volume` over candidate faces divided
  by block volume, clamped to `[0, 1]`. A block that no candidate face meets (no
  part of a face lies in the closed block) is wholly inside or outside, so its
  fraction is exactly 1 or 0 from `locate` of its centroid, with no integration
  residue from the faces above it. For rotated blocks the solid is transformed
  once into the block frame (rotation transpose about the first centroid), then the
  same axis-aligned integration runs. Parallel per block, input order.
- Free functions mirroring the methods for one-off calls: `locate(solid, points)`,
  `signed_distance(solid, points)`, `block_fraction(solid, blocks)`,
  `trace_intervals(solid, trace)`, `closest(mesh, points)`, `ray_hits(mesh, o, d)`.

## 8. wireframes-build

Files: `crates/build/src/{lib.rs, dtm.rs, polygon.rs, loft.rs, extrude.rs, close.rs,
decimate.rs, smooth.rs, offset.rs, remesh.rs, contour.rs}`.

```rust
#[derive(Clone, Debug, Default)]
pub struct DtmOptions<'a> {
    pub breaklines: Option<&'a Polylines>,
    pub boundary: Option<&'a Polylines>,  // closed rings, even-odd; triangles outside are dropped
    pub min_spacing: f64,                 // plan duplicate merge distance, default 0 (exact)
    pub max_edge: Option<f64>,            // drop boundary triangles with a plan edge longer than this, repeatedly
}
pub fn dtm(points: &[[f64; 3]], opts: &DtmOptions) -> Result<Mesh>;
pub fn triangulate_polygons(polygons: &Polylines) -> Result<Mesh>;
pub fn loft(sections: &Polylines, opts: &LoftOptions) -> Result<Mesh>;
pub struct LoftOptions { pub cap_ends: bool } // default true
pub fn extrude(polygons: &Polylines, vector: [f64; 3]) -> Result<Mesh>;
pub fn close_to_elevation(surface: &Mesh, z: f64) -> Result<Mesh>;
pub struct DecimateOptions { pub target_faces: Option<usize>, pub max_error: Option<f64>, pub preserve_boundary: bool }
pub fn decimate(mesh: &Mesh, opts: &DecimateOptions) -> Result<Mesh>;
pub fn smooth_taubin(mesh: &Mesh, iterations: usize, lambda: f64, mu: f64) -> Result<Mesh>;
pub fn offset(mesh: &Mesh, distance: f64) -> Result<Mesh>;
pub fn remesh(mesh: &Mesh, edge_length: f64, iterations: usize) -> Result<Mesh>;
pub fn contours(surface: &Mesh, levels: &[f64]) -> Result<Polylines>;
```

Semantics:
- `dtm`: plan constrained Delaunay (spade) of the points, breakline vertices and
  boundary vertices; breaklines and boundary rings are constraint edges; crossing
  constraints split at the plan intersection with z interpolated along the
  breakline that is inserted first (feature order); points closer than
  `min_spacing` in plan to an earlier point are dropped, keeping the earlier one;
  at `min_spacing = 0` only exact 3D duplicates merge, and two input vertices
  (points, then breakline vertices, then boundary vertices) at one plan position
  with different z give `InvalidInput` naming both, never a silent choice;
  faces wind counter-clockwise in plan (normals up). Errors: fewer than 3 points
  not all collinear in plan, non-finite input.
- `triangulate_polygons`: per feature, best-fit plane (Newell normal), CDT of its
  closed parts with even-odd holes, lifted back to the original vertices (no new
  vertices unless rings cross); face column `feature` (Int64). Normals follow the
  Newell normal of the outer ring.
- `loft`: features in order are consecutive sections; each must have exactly one
  part, a closed ring, except that the first and last sections may be an open part
  (a line, or a single point) where the body pinches out. An open part links as the
  degenerate ring that runs along it and back, and gets no cap. Rings are oriented consistently (Newell normals agree with the
  direction to the next section's centre) and each closed ring then starts at its
  lexicographically smallest rotation, so the result (vertices, triangles and
  volume, bit for bit) does not depend on where or in which direction a ring was
  digitized. Each band is linked in its overlay
  frame: both rings are centred on their own area centroid, the component along
  the axis between the two centres is removed, so the band length does not enter
  the cost, and when both rings have area the second is scaled about its centre to
  the area of the first, so similar sections (a shrinking non-convex outline) link
  vertex to vertex. Start vertices are aligned by the closest overlay pair, and dynamic
  programming minimizes the total overlay triangle area over the band (O(m n)).
  Measured in 3D, a long band's double fan and its straight walls have nearly the
  same area, so a 3D cost cannot separate them; the overlay cost can. End caps by `triangulate_polygons`. Output is a solid
  (validated; `InvalidInput` with the report summary if linking produced
  self-intersections). Face column `section` (Int64, band index, caps -1 and -2).
- `extrude`: each closed ring triangulated as caps at 0 and `vector`, walls as
  quads split on the shorter diagonal; output is an outward solid; `vector` must not
  be parallel to the ring plane.
- `close_to_elevation`: walls from the surface boundary loop straight down (or up)
  to `z`, a flat base triangulated from the projected boundary; output solid; the
  surface must be a height field with one boundary loop entirely above `z` or
  entirely below.
- `decimate`: Garland-Heckbert quadric collapse, optimal placement solved in f64
  with fallback to edge endpoints or midpoint; a collapse is rejected when it breaks
  the link condition, flips a face normal (exact `side_of_plane` test), or (with
  `preserve_boundary`) moves a boundary vertex off its boundary line; stops at
  `target_faces` or when the next cost exceeds `max_error` squared; ties by edge
  index; face and vertex attributes of survivors are kept.
- `smooth_taubin`: uniform Laplacian weights, `lambda > 0`, `mu < -lambda`,
  boundary vertices fixed; defaults in Python `lambda=0.5, mu=-0.53`.
- `offset`: vertices move by `distance` along angle-weighted vertex normals
  (outward positive); output validated and refused with `InvalidMesh` when it
  self-intersects.
- `remesh`: Botsch-Kobbelt isotropic remeshing (split above 4/3 L, collapse below
  4/5 L, valence flips, tangential relaxation, projection to the input by
  `SurfaceIndex::closest`); boundary edges split only; face attributes copied from
  the closest input face, vertex attributes dropped.
- `contours`: per level, `slice_by_plane` with origin `(0,0,0)`, normal `(0,0,1)`
  and `offsets = [level]` (a vertex on the level counts above, so a flat area at
  the level is traced where the surface drops below it), followed by each line of
  `slice_by_plane` with normal `(0,0,-1)` and `offsets = [-level]` that shares no
  segment (exact coordinates, either direction) with the first trace, reversed to
  its winding: the toe of a flat area at the level that the surface does not cross
  there, such as a pit floor or the inside of a berm. A line that shares a segment
  is the same contour going round a flat triangle on its other side and is left
  out. The `offset` column is renamed `elevation`.

## 9. Determinism and parallelism summary

Parallel: `validate` (pair tests), `SolidIndex::locate_all`, `contains_all`,
`signed_distance`, `block_fraction`, `SurfaceIndex::closest`, `ray_hits`,
`surface_z`, `face_areas`. Each collects in input order; any sum follows section 2.
Booleans may run tri-tri tests in parallel and sort their outputs by face index
before use. Every crate has a test that runs one function with
`rayon::ThreadPoolBuilder::new().num_threads(1)` and with 4 threads and asserts
bitwise equal results.

## 10. wireframes-py and the Python API

### 10.1 Rust binding layout

`crates/py/Cargo.toml`: `crate-type = ["cdylib"]`, lib name `_wireframes`, pyo3
features `abi3-py311`, `extension-module`, `multiple-pymethods`. Each area adds its
methods to the shared pyclasses in its own file with a separate `#[pymethods]`
block. `lib.rs` declares every module and registers every class and function.

Files: `crates/py/src/{lib.rs, errors.rs, containers.rs, validate.rs, io.rs, ops.rs,
boolean.rs, query.rs, build.rs}`.

- `containers.rs`: `PyTable` (`wireframes.Table`), `PyMesh` (`wireframes.Mesh`),
  `PyPolylines`, `PyPointSet`, `PyUnstructuredGrid`, numpy conversion helpers
  (`points_arg` accepts `(n, 3)` or `(n, 2)` with z = 0 where stated, a `PointSet`, or
  a `Mesh` meaning its vertices). Caches `SolidIndex` and `SurfaceIndex` lazily on
  `PyMesh` (`OnceLock`).
- `errors.rs`: maps `Error` to Python exceptions created in Rust and re-exported by
  `python/wireframes/errors.py`.
- All long-running calls release the GIL (`py.detach` or `allow_threads`, whichever
  the pinned pyo3 provides).

### 10.2 Python package

Files: `python/wireframes/{__init__.py, _wireframes.pyi, errors.py, _containers.py, plot.py,
py.typed}`. `import wireframes as wf`; flat exports listed in `__all__`.

`errors.py`:

```python
class WireframesError(Exception): ...


class InvalidInput(WireframesError, ValueError): ...


class InvalidMesh(InvalidInput):
    report: MeshReport


class FileError(WireframesError, OSError): ...  # Error::Io and Error::Format


class MissingColumn(WireframesError, KeyError): ...
```

`Unsupported`, `CrsMismatch` and `UnitMismatch` raise `InvalidInput`.

Conventions (boitata): required geometry args positional, options keyword-only;
properties for cheap derived values, methods for anything taking points; objects
immutable; choices as string literals; numpy arrays out (`float64`, `int64`,
`bool`); arrays that copy a container's state (`coords`, `x`/`y`/`z`,
`triangles`, `parts`, `closed`, `feature`, columns from `[]` or `column`, grid
`points`, `cell_types`, `offsets`, `connectivity`) are read-only, so a write
raises `ValueError` instead of vanishing; `.copy()` gives a writeable array, and
`with_vertices` or `with_*_column` build a changed container; NumPy-style docstrings (a `Parameters` section whenever the callable takes arguments, and an `Examples` section that runs as a doctest with `wf` and `np` in scope), present tense, first line "A copy with ..." or a
noun phrase; messages lowercase.

Type aliases in the stub:

```python
ArrayLike = npt.ArrayLike
Path: TypeAlias = str | PathLike[str]
Points: TypeAlias = ArrayLike | PointSet | Mesh
Values: TypeAlias = ArrayLike | Sequence[str | None]
Floats = npt.NDArray[np.float64]
```

#### Table

```python
class Table:
    def __init__(self, data: Mapping[str, Values] | None = None) -> None: ...

    num_rows: int  # property
    column_names: list[str]
    units: dict[str, str]

    def column(self, name: str) -> np.ndarray: ...  # float64 (NaN null), int64, bool, object (None null)
    def __getitem__(self, name: str) -> np.ndarray: ...
    def mask(self, name: str) -> npt.NDArray[np.bool_]: ...  # True where valid
    def filter(self, mask: ArrayLike) -> Table: ...
    def drop_null(self, *columns: str) -> Table: ...
    def with_column(self, name: str, values: Values, *, unit: str | None = None) -> Table: ...
    def with_columns(self, data: Mapping[str, Values]) -> Table: ...
    def to_dict(self) -> dict[str, np.ndarray]: ...
    def to_pyarrow(self) -> Any: ...
    def to_polars(self) -> Any: ...
    def to_pandas(self) -> Any: ...
    def __len__(self) -> int: ...
```

Int64 and Bool columns with nulls come out as their numpy dtype with
`mask(name)` giving validity; `to_pyarrow`, `to_polars` and `to_pandas` apply the
mask as nulls.

#### Mesh

```python
class Mesh:
    def __init__(
        self,
        vertices: ArrayLike,
        triangles: ArrayLike,
        *,
        vertex_attributes: Mapping[str, Values] | Table | None = None,
        face_attributes: Mapping[str, Values] | Table | None = None,
        crs: str | None = None,
        length_unit: str | None = None,
    ) -> None: ...

    # properties
    coords: Floats  # (n, 3)
    x: Floats
    y: Floats
    z: Floats
    triangles: npt.NDArray[np.int64]  # (m, 3)
    n_vertices: int
    n_faces: int
    bounds: tuple[list[float], list[float]] | None
    column_names: list[str]  # vertex columns, then face columns
    vertex_attributes: Table
    face_attributes: Table
    crs: str | None
    length_unit: str | None
    is_closed: bool
    is_solid: bool  # validate(tolerance=0).summary["is_solid"]
    area: float
    volume: float  # raises InvalidMesh when not closed and consistently wound
    centroid: Floats  # volume centroid, same rule
    face_areas: Floats
    face_normals: Floats
    face_centroids: Floats
    vertex_normals: Floats
    dip: Floats  # (m, 2) dip, dip direction in degrees

    def __len__(self) -> int: ...  # number of vertices (boitata)
    def __getitem__(
        self, name: str
    ) -> np.ndarray: ...  # vertex column, else face column; a name in both raises InvalidInput
    # attributes and editing (containers.rs, ops.rs)
    def with_vertex_column(self, name: str, values: Values, *, unit: str | None = None) -> Mesh: ...
    def with_face_column(self, name: str, values: Values, *, unit: str | None = None) -> Mesh: ...
    def with_attributes(
        self, *, vertex: Table | Mapping | None = None, face: Table | Mapping | None = None
    ) -> Mesh: ...
    def with_vertices(self, vertices: ArrayLike) -> Mesh: ...
    def filter_faces(self, mask: ArrayLike) -> Mesh: ...
    def to_table(
        self, *, on: Literal["vertices", "faces"] = "vertices"
    ) -> Table: ...  # vertices: x, y, z + vertex columns; faces: v0, v1, v2 + face columns
    def to_polars(self, *, on: Literal["vertices", "faces"] = "vertices") -> Any: ...
    def to_pandas(self, *, on: Literal["vertices", "faces"] = "vertices") -> Any: ...
    def to_pyarrow(self, *, on: Literal["vertices", "faces"] = "vertices") -> Any: ...
    def transform(self, matrix: ArrayLike) -> Mesh: ...  # (4, 4)
    def translate(self, offset: ArrayLike) -> Mesh: ...
    def flip(self) -> Mesh: ...
    def to_length_unit(self, unit: str) -> Mesh: ...
    # validation and repair (validate.rs, ops.rs)
    def validate(self, *, tolerance: float = 0.0) -> MeshReport: ...
    def repair(
        self,
        *,
        weld: float | None = 0.0,
        remove_degenerate: bool = True,
        remove_duplicate_faces: bool = True,
        remove_unused_vertices: bool = True,
        orient: bool = True,
        fill_holes: Literal["fan", "min_area"] | None = None,
        weld_touching_shells: bool = False,
    ) -> tuple[Mesh, RepairReport]: ...
    def components(self) -> list[Mesh]: ...
    def component_labels(self) -> npt.NDArray[np.int64]: ...
    # measures (ops.rs)
    def volume_by_elevation(self, edges: ArrayLike) -> Floats: ...
    def surface_z(self, points: Points) -> Floats: ...  # (n,2) or (n,3), z ignored
    def vertical_distance(self, points: Points) -> Floats: ...
    def thickness(
        self, bottom: Mesh, points: Points, *, mode: Literal["vertical", "normal"] = "vertical"
    ) -> Floats: ...
    # booleans (boolean.rs)
    @overload
    def union(self, other: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
    @overload
    def union(
        self, other: Mesh, *, repair: Literal[True]
    ) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
    # intersection, difference, symmetric_difference: same two overloads
    def clip_by_surface(
        self, surface: Mesh, *, repair: bool = False
    ) -> tuple[Mesh, Mesh]: ...  # (below, above)
    def split_by_plane(
        self, origin: ArrayLike, normal: ArrayLike
    ) -> tuple[Mesh, Mesh]: ...  # (negative, positive)
    def slice_by_plane(
        self, origin: ArrayLike, normal: ArrayLike, *, offsets: ArrayLike | None = None
    ) -> Polylines: ...
    # queries (query.rs)
    def locate(self, points: Points) -> npt.NDArray[np.int8]: ...  # -1 inside, 0 boundary, 1 outside
    def contains(self, points: Points, *, boundary: bool = True) -> npt.NDArray[np.bool_]: ...
    def signed_distance(self, points: Points) -> Floats: ...
    def distance(self, points: Points) -> Floats: ...
    def closest_point(
        self, points: Points
    ) -> tuple[Floats, npt.NDArray[np.int64], Floats]: ...  # points, face, distance
    def ray_hits(
        self, origins: ArrayLike, directions: ArrayLike
    ) -> Table: ...  # ray, t, face, x, y, z, kind, entering (null when open)
    def segment_hits(
        self, starts: ArrayLike, ends: ArrayLike
    ) -> Table: ...  # segment, t, face, x, y, z, kind, entering
    def trace_intervals(
        self, trace: ArrayLike, depths: ArrayLike
    ) -> Table: ...  # from, to, inside, entry_face, exit_face
    def block_fraction(
        self, centroids: Points, size: ArrayLike, *, rotation: ArrayLike | None = None
    ) -> Floats: ...
    # construction and editing (build.rs)
    def close_to_elevation(self, z: float) -> Mesh: ...
    def decimate(
        self,
        *,
        target_faces: int | None = None,
        max_error: float | None = None,
        preserve_boundary: bool = True,
    ) -> Mesh: ...
    def smooth(self, *, iterations: int = 10, lambda_: float = 0.5, mu: float = -0.53) -> Mesh: ...
    def offset(self, distance: float) -> Mesh: ...
    def remesh(self, edge_length: float, *, iterations: int = 5) -> Mesh: ...
    def contours(self, levels: ArrayLike) -> Polylines: ...
```

`MeshReport`: `.summary -> dict[str, int | bool]` (keys from `to_pairs`, plus
`"is_solid"`), `.problems -> Table` (columns `kind` (str), `face`, `vertex`, `edge`,
`other` as int64 with nulls), `.is_solid -> bool`, `__len__` = number of problems.

`RepairReport`: `.actions -> Table` (`action` str in `vertex_merged, vertex_split,
vertex_removed, face_removed, face_flipped, face_added`, `vertex`, `face`, `other`
(merge target, split copy or hole), `reason`),
`.vertex_map`, `.face_map` (int64, -1 removed), `.before`, `.after` (dicts),
`.holes_left`, `.summary -> dict[str, int]` (count per action).

#### Polylines

```python
class Polylines:
    def __init__(self, parts: Sequence[ArrayLike], *, closed: bool | Sequence[bool] = False,
                 features: Sequence[int] | None = None,
                 attributes: Mapping[str, Values] | Table | None = None,
                 vertex_attributes: Mapping[str, Values] | Table | None = None,
                 crs: str | None = None, length_unit: str | None = None) -> None: ...
    coords: Floats; x: Floats; y: Floats; z: Floats
    parts: list[Floats]; closed: npt.NDArray[np.bool_]; feature: npt.NDArray[np.int64]  # per part
    attributes: Table; vertex_attributes: Table
    bounds: tuple[list[float], list[float]] | None
    column_names: list[str]; crs: str | None; length_unit: str | None
    def __len__(self) -> int: ...                 # features
    def filter(self, mask: ArrayLike) -> Polylines: ...
    def with_column(self, name: str, values: Values, *, unit: str | None = None) -> Polylines: ...
    def with_vertex_column(self, name: str, values: Values, *, unit: str | None = None) -> Polylines: ...
    def length(self) -> Floats: ...
    def area(self) -> Floats: ...
    def contains(self, points: Points, *, feature: int | None = None) -> npt.NDArray[np.bool_]: ...
    def locate(self, points: Points) -> npt.NDArray[np.int64]: ...  # -1 when no feature holds the point
    def distance(self, points: Points, *, feature: int | None = None, signed: bool = False) -> Floats: ...
    def to_table(self) -> Table: ...              # long: feature, part, vertex, x, y, z, closed + vertex columns + feature columns
    def to_polars(self) -> Any: ...; def to_pandas(self) -> Any: ...; def to_pyarrow(self) -> Any: ...
    @staticmethod
    def from_table(table: Table | Mapping[str, Values], *, feature: str = "ID", x: str = "X", y: str = "Y",
                   z: str | None = None, part: str | None = None, closed: bool = False,
                   crs: str | None = None, length_unit: str | None = None) -> Polylines: ...
    def to_points(self) -> PointSet: ...           # adds feature and part columns
```

#### PointSet, UnstructuredGrid, Dataset

```python
class PointSet:
    def __init__(self, coords: ArrayLike, attributes: Mapping[str, Values] | Table | None = None, *,
                 crs: str | None = None, length_unit: str | None = None) -> None: ...
    coords: Floats; x: Floats; y: Floats; z: Floats; attributes: Table
    bounds: tuple[list[float], list[float]] | None; column_names: list[str]
    crs: str | None; length_unit: str | None
    def __len__(self) -> int: ...
    def __getitem__(self, name: str) -> np.ndarray: ...
    def filter(self, mask: ArrayLike) -> PointSet: ...
    def with_column(self, name: str, values: Values, *, unit: str | None = None) -> PointSet: ...
    def to_table(self) -> Table: ...
    def to_polars(self) -> Any: ...; def to_pandas(self) -> Any: ...; def to_pyarrow(self) -> Any: ...

class UnstructuredGrid:
    points: Floats; cell_types: npt.NDArray[np.uint8]; offsets: npt.NDArray[np.int64]
    connectivity: npt.NDArray[np.int64]; point_data: Table; cell_data: Table
    n_cells: int
    def surface(self) -> Mesh: ...

class Dataset(NamedTuple):                        # _containers.py
    mesh: Mesh | None; polylines: Polylines | None; points: PointSet | None
    grid: UnstructuredGrid | None; length_unit: str | None; warnings: list[str]
```

#### Free functions

```python
def read(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> Dataset: ...
def read_mesh(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> Mesh: ...
def read_polylines(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> Polylines: ...
def read_points(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> PointSet: ...
def write(path: Path, data: Mesh | Polylines | PointSet | UnstructuredGrid | Dataset, *,
          ascii: bool = False, dxf_entity: Literal["3dface", "polyface"] = "3dface",
          zlib: bool = False, skip_unsupported_columns: bool = False) -> None: ...
def concat(meshes: Sequence[Mesh]) -> Mesh: ...
def union(a: Mesh, b: Mesh, *, repair: bool = False) -> Mesh | tuple[Mesh, tuple[RepairReport, RepairReport]]: ...  # same overloads as the method
def intersection(...); def difference(...); def symmetric_difference(...)   # same as union
def cut_fill(base: Mesh, design: Mesh, *, boundary: Polylines | None = None) -> dict[str, float]: ...  # cut, fill, net, cut_area, fill_area, area
def dtm(points: Points, *, breaklines: Polylines | None = None, boundary: Polylines | None = None,
        min_spacing: float = 0.0, max_edge: float | None = None) -> Mesh: ...
def triangulate(polygons: Polylines) -> Mesh: ...
def loft(sections: Polylines, *, cap_ends: bool = True) -> Mesh: ...
def extrude(polygons: Polylines, vector: ArrayLike) -> Mesh: ...
```

`read*` set `crs` and `length_unit` on every returned object; a file unit (DXF) is
used when the argument is `None`, and an argument that disagrees with the file unit
raises `InvalidInput`. `write` to binary STL warns (`UserWarning`) when float32
rounding moves any vertex, with the count, the largest change and the way out
(`ascii=True` or a local origin); the file is still written.

#### plot.py (matplotlib, gated behind `wireframes[plot]`)

```python
def plan(
    data: Mesh | Polylines | PointSet | Sequence[Mesh | Polylines | PointSet],
    values: str | ArrayLike | None = None,
    *,
    ax=None,
    cmap: str = "viridis",
    edges: bool = False,
    label: bool = False,
) -> tuple[Figure, Axes]: ...
def section(
    data: Mesh | Polylines | PointSet | Sequence[...],
    origin: ArrayLike,
    normal: ArrayLike,
    *,
    thickness: float | None = None,
    values: str | None = None,
    ax=None,
) -> tuple[Figure, Axes]: ...
def surface3d(
    data: Mesh | Polylines | PointSet | Sequence[...],
    values: str | ArrayLike | None = None,
    *,
    ax=None,
    cmap: str = "viridis",
    problems: MeshReport | None = None,
    azim: float = -60.0,
    elev: float = 30.0,
) -> tuple[Figure, Axes]: ...
```

- `plan`: meshes as `tripcolor` (face column or vertex column, else z), polylines as
  lines (rings closed), points as scatter; equal aspect; xy axes in the data unit.
- `section`: meshes drawn as `slice_by_plane(origin, normal)` traces; polylines and
  points within `thickness / 2` of the plane projected onto it; axes are the
  in-plane horizontal distance from `origin` and elevation (vertical sections) or
  the plane basis (others).
- `surface3d`: `plot_trisurf`-style `Poly3DCollection`, lines and points in 3D;
  `problems` draws boundary and non-manifold edges in red and problem faces in
  orange; equal box aspect.
- Each returns `(fig, ax)`; with `ax` given, `fig = ax.figure`. Importing
  `wireframes.plot` without matplotlib raises `ImportError` with the install hint
  `pip install "wireframes[plot]"`.

`__init__.py` exports: `Table, Mesh, Polylines, PointSet, UnstructuredGrid,
Dataset, MeshReport, RepairReport, WireframesError, InvalidInput, InvalidMesh,
FileError, MissingColumn, read, read_mesh, read_polylines, read_points, write,
concat, union, intersection, difference, symmetric_difference, cut_fill, dtm,
triangulate, loft, extrude, plot` (lazy attribute), `__version__`.

## 11. Tooling

`mise.toml` tasks (lead owns):

| task | runs |
|---|---|
| `build` | `cargo build --workspace` |
| `py:sync`, `py:build` | `uv sync`; `maturin develop --release` inside the uv venv |
| `test` | `py:build` first, then `cargo test --workspace --exclude wireframes-py`, then `pytest python/tests` |
| `lint` | `cargo fmt --check`, `cargo clippy --workspace --all-targets -- -D warnings`, `ruff check`, `ruff format --check` |
| `verify` | `lint` then `test` |
| `bench` | `cargo bench` (criterion), interleaved rounds, medians |
| `examples`, `docs`, `docs-build`, `compat` | as in AGENTS.md |

Proptest cases: 256 per property by default (`PROPTEST_CASES` overrides). The
`[env]` table of `mise.toml` sets `PROPTEST_RNG_SEED`, so every run through `mise`
generates the same cases; committed `proptest-regressions` files replay earlier
failures on top of them. To explore new cases, override the seed for one run
(`PROPTEST_RNG_SEED=<n> mise run test`).

`meshio` and `trimesh` are dev-only dependencies used inside tests as independent
writers and readers; tests that need them skip when they are missing.

## 12. Invariants and theory checks

Standing tests, by crate. Each is a Rust test (`tests/invariants.rs` per crate) and,
where a Python API exists, a pytest twin. "P" marks a proptest.

wireframes-core
- P `Table::take` then `concat` of the parts in order equals the original table.
- P `Mesh::new` rejects every out-of-range index and every non-finite coordinate.
- P `predicates` signs agree with a `num-rational` exact oracle on random and
  near-degenerate inputs (points snapped to a coarse grid, then perturbed by one ulp).
- P `exact::orient3d_i` and `cmp_axis` agree with the rational oracle for LPI and TPI
  points; `to_f64` returns the nearest f64 of the rational value.
- P `tri_tri` is symmetric in its arguments and invariant under corner rotation.
- P `Bvh::query` returns exactly the boxes a brute-force scan returns.
- `validate` on boitata's fixture set gives boitata's counts (ported fixtures), plus
  new fixtures for touching faces and coplanar overlaps. Exception by section 0
  item 18: tetrahedra sharing one vertex are a solid with one shell contact, where
  boitata reports a non-manifold vertex.
- Closed shells touching at a vertex or along an edge, by shared index or by
  coincident copies, are solids with the expected contact counts; a contact whose
  interiors overlap is a `SelfIntersection`. The same holds for one shell touching
  itself (`fixtures::voxels` of `EDGE_LOOP` and `VERTEX_CHAIN`).
- P `validate` is invariant under face order permutation (counts), vertex index
  relabelling (counts), and rigid translation by integers (counts).
- A unit cube and an icosphere report `is_solid`; removing one face gives exactly 3
  boundary edges.

wireframes-io
- P write then read is identity for canonical meshes, per format, ascii and binary.
- P read then write then read is identity on every fixture.
- Fixtures: binary STL whose header starts with "solid"; VTP with both Polys and
  Strips and cell data; multi-piece VTU; polyface with negative indices; Datamine
  pair with gapped PIDs and `;` delimiter; PLY with an unknown element between
  vertex and face; OBJ with negative indices and forward references; DXF with cp1252
  layer names; VTU with a polyhedron cell.
- `UnstructuredGrid` tetra and hexa cells round-trip through `.vtu` and `.vtk`
  with identical types, offsets, connectivity and cell data.
- Polygon faces read from OBJ/PLY/OFF write back as the same polygons.

wireframes-ops
- P `volume` of a closed mesh is invariant under translation (within `close`) and
  changes sign under `flip`.
- P `volume` of an axis-aligned box mesh equals `dx * dy * dz`; `area` equals the
  analytic area; `centroid` equals the box centre.
- P `sum(volume_by_elevation(edges)) == volume` when edges span the bounds.
- P `box_volume` summed over a grid of boxes covering the bounds equals `volume`.
- P `repair` is idempotent; `repair` output of a shuffled, duplicated, partly flipped
  copy of a solid is a solid with the original volume.
- P `vertex_map` and `face_map` are consistent with the output (every mapped face's
  corners map to the output face's corners).
- P `components` faces sum to the input's face count and their volumes sum to the
  input volume.
- P `thickness(top, bottom, Vertical)` equals the analytic difference for two planes.

wireframes-boolean
- P `vol(A union B) + vol(A intersection B) == vol(A) + vol(B)`.
- P `vol(A difference B) == vol(A) - vol(A intersection B)`.
- P `vol(A symmetric_difference B) == vol(A union B) - vol(A intersection B)`.
- `vol(A xor B) == vol(A) + vol(B) - 2 vol(A and B)` and `vol(A or B) + vol(A and B)
  == vol(A) + vol(B)` on 160 seeded random grid-box pairs, with inputs of touching
  shells (edge and vertex contact) and nested pairs; every result validates as a
  solid (`tests/touching.rs`).
- Edge, vertex and face contact pairs, a difference that leaves bodies joined at
  edges and vertices, and a chain of ten mixed booleans each fed its previous result,
  checked against an exact voxel count.
- Four shells on one edge and four octant boxes meeting at one vertex; a cavity
  touching its outer shell at every vertex it has; exact edge and vertex contacts at
  UTM scale and contacts 1e-6 apart; shells touching along a face are a coplanar
  overlap, not a contact. Inputs with permuted vertices, faces and corners give the
  same faces (as coordinate triples) and the same summary or the same error, and 1
  and 8 threads give identical results (`tests/touching.rs`).
- `split_by_plane` of a closed mesh that meets the plane only at a vertex, with and
  without touching shells, puts it whole on one side.
- P closed manifold in gives closed manifold out (`is_closed`, no non-manifold edges
  or vertices, consistent winding) for every op.
- P `A union A` and `A intersection A` have the volume of A; `A difference A` is empty.
- P commutativity of union and intersection in volume and in `locate` of random
  points.
- P points sampled at random: `locate` on the result agrees with the boolean of
  `locate` on the inputs, for points not on either boundary.
- Coplanar cases: two unit cubes sharing a face, shifted by half a cell, nested,
  and one face-on-face with opposite facing; each gets an exact expected volume.
- UTM-scale coordinates (x near 7e5, y near 7.5e6) give the same volumes, within
  `close`, as the same solids at the origin.
- P `clip_by_surface`: `vol(below) + vol(above) == vol(solid)`, both closed.
- P `split_by_plane`: volumes add up; open meshes keep their area.
- P `slice_by_plane` of a closed mesh gives only closed rings; ring areas of a
  horizontal slice of a box equal its footprint.
- P `cut_fill(S, S)` gives zero cut and fill; `cut_fill(base, base shifted up by h)`
  gives `fill == h * area`; swapping arguments swaps cut and fill.
- Invalid inputs (open mesh, self-intersecting mesh) raise `InvalidMesh` with the
  report; `repair=True` on a mesh with duplicate vertices succeeds and returns the
  report.

wireframes-query
- P sign of `signed_distance(p)` equals `locate(p)` value.
- P `locate` agrees with the winding number (rational oracle on small meshes) for
  points off the surface, including points sharing x or y with vertices.
- P mesh vertices locate as `Boundary`.
- P `sum(block_fraction * block_volume)` over a block grid covering the solid equals
  the solid volume (README M5 acceptance `block_fraction_sums_to_solid_volume`), for
  axis-aligned and rotated grids.
- P `block_fraction` is 1 for blocks strictly inside and 0 for blocks strictly
  outside (decided by `locate` of corners on a convex solid).
- P `trace_intervals` inside lengths along a straight line through a convex solid
  equal the chord length.
- P `closest` distance is never larger than the distance to any vertex, and equals
  a brute-force scan.
- P `ray_hits` count parity from an outside origin matches `locate` of the ray end.
- Thread count: `locate_all` and `block_fraction` equal bitwise for 1 and 4 threads.

wireframes-build
- P `dtm` output is Delaunay away from constraints (exact `incircle` on every
  unconstrained edge) and contains every breakline segment as a chain of edges.
- P `dtm` vertices are exactly the input points minus merged duplicates; faces wind
  up in plan.
- P `triangulate_polygons` area equals `Polylines.area()` for planar xy polygons.
- P `extrude` volume equals plan area times `vector.z` for vertical extrusions.
- P `loft` of two identical squares offset by h is a box of volume `a * h`.
- P `close_to_elevation` volume equals `cut_fill` against a flat plane at `z`.
- P `decimate` keeps a solid closed and manifold; its volume error is reported, not
  bounded, except for planar inputs, where volume and boundary stay exact.
- P `smooth_taubin` keeps boundary vertices fixed and preserves topology.
- P `remesh` keeps a closed input closed; Hausdorff distance (sampled) to the input
  stays below `edge_length`.
- P `contours` of a plane `z = ax + by + c` at level L are straight lines lying on
  the plane.

wireframes-py
- Every name in `_wireframes.pyi` exists at runtime and every public runtime name is in
  the stub (`python/tests/test_api.py`).
- Error mapping: each Rust `Error` variant raises the mapped Python class.
- numpy dtypes and shapes of every property match the stub.
- Plot functions return `(fig, ax)` and reuse a passed `ax`.

## 13. Ownership

Each area edits only its files. The lead creates the skeleton first (workspace,
empty crates with the module files below declared, `pyproject.toml`, `mise.toml`,
`python/wireframes/__init__.py` exporting every name, `crates/py/src/lib.rs` registering
every module) so areas never touch shared files. A needed change to a shared file
or to this spec goes to the lead.

| area | owner of | files |
|---|---|---|
| L lead | workspace and contract | `Cargo.toml`, `Cargo.lock`, `rust-toolchain.toml`, `mise.toml`, `pyproject.toml`, `README.md`, `AGENTS.md`, `CLAUDE.md`, `docs/design/**`, `docs/mkdocs.yml`, `crates/*/Cargo.toml`, `crates/py/src/lib.rs`, `python/wireframes/__init__.py` |
| A core | kernel | `crates/core/src/**`, `crates/core/tests/**`, `crates/core/benches/**`, `crates/py/src/validate.rs`, `python/tests/test_validate.py` |
| B io | formats | `crates/io/src/**`, `crates/io/tests/**` (fixtures included), `crates/io/benches/**`, `crates/py/src/io.rs`, `python/tests/test_io_m1.py`, `python/tests/test_io_m2.py`, `examples/01-files/**` |
| C ops | measure, repair | `crates/ops/src/**`, `crates/ops/tests/**`, `crates/ops/benches/**`, `crates/py/src/ops.rs`, `python/tests/test_measure_repair_plot.py`, `examples/02-measure-repair/**` |
| D boolean | arrangements | `crates/boolean/src/**`, `crates/boolean/tests/**`, `crates/boolean/benches/**`, `crates/py/src/boolean.rs`, `python/tests/test_booleans.py`, `examples/03-booleans/**` |
| E query | queries | `crates/query/src/**`, `crates/query/tests/**`, `crates/query/benches/**`, `crates/py/src/query.rs`, `python/tests/test_queries.py`, `examples/04-queries/**` |
| F build | construction | `crates/build/src/**`, `crates/build/tests/**`, `crates/build/benches/**`, `crates/py/src/build.rs`, `python/tests/test_construct.py`, `examples/05-construct/**` |
| G python | bindings shell and Python package | `crates/py/src/{errors.rs, containers.rs}`, `python/wireframes/{_wireframes.pyi, errors.py, _containers.py, plot.py, py.typed}`, `python/tests/{test_containers.py, test_api.py, conftest.py}`, `examples/00-basics/**`, `docs/**` except `docs/design/**` |

Rules for cross-area work:
- `test_measure_repair_plot.py` (C) calls `wf.plot.*` through the API in 10.2; G
  implements it.
- G writes the full `_wireframes.pyi` from this spec at the start, docstrings included,
  so areas never edit it; an area that finds a stub mismatch reports it to G.
- Areas depend only on the published signatures of lower crates. Until a lower crate
  lands, an area may write a local test double inside its own `tests/` folder, never
  inside another crate.
- Order of landing: L skeleton, then A (core) and B (io) and G in parallel, then C,
  then D and E in parallel, then F.

## 14. Skeleton notes (lead, after the skeleton landed)

These refine the sections above; where they differ, this section wins.

- Bindings: each `crates/py/src/<area>.rs` has `pub(crate) fn register(m)`, which
  `lib.rs` calls; an area adds its functions and classes there, so `lib.rs` stays
  untouched. `lib.rs` carries `#![allow(dead_code)]` until G lands the containers;
  G removes it. `errors::to_py` maps `Error` to the Python classes (done).
- `wireframes.read` is `_containers.read` (G), which wraps `_wireframes._read` (B) returning
  `(mesh, polylines, points, grid, length_unit, warnings)` in a `Dataset`.
- Unimplemented Rust bodies are `todo!()`; unimplemented Python functions raise
  `NotImplementedError`. Stub crates carry `#![allow(unused_variables)]` in
  `lib.rs`, which the owner removes when the crate is complete.
- `mkdocs.yml` lives at the repository root (mkdocs refuses a config inside
  `docs_dir`); `docs/design/` is excluded from the site.
- Core additions beyond section 3, all public: `Error::{invalid, format,
  format_at, with_path}`; `Column::with_valid`, `Column::with_unit`,
  `Table::columns`; `ColumnData` equality compares floats by bits with every NaN
  equal (so null rows compare equal after a round trip); `predicates::{Sign::of,
  to_i8, flip, times}`, `kept_axes(drop)` (cyclic axes, so `orient2d_axis` has the
  sign of normal component `drop`), `projection_axis(t)`;
  `exact::{Expansion::zero, diff_f64, neg, components}`,
  `ImplicitPoint::homogeneous`, `cmp_lex`; `Polylines::{part_offsets,
  feature_offsets}`; `Bvh::{is_empty, item_box}`; `validate::{degenerate_mask,
  earliest_within}`; `topology::{locate_point, ray_crosses}` (the section 7.1
  ray rule, brute force over a face list: core needs it for `InwardShell`, ops for
  orientation in `repair`; query adds the BVH); `topology::locate_implicit` (the
  same rule for an exact `ImplicitPoint`) and `topology::shell_location(vertices,
  inner, outer) -> Option<i8>` (nesting of one closed shell in another; a point
  inside an edge is the `Lpi` of the edge with the axis plane through a float
  strictly between two of its coordinates; `None` when every vertex and edge lies
  on `outer`).
- `orient3d_i(a, b, c, d)` has the sign of `det[a - d; b - d; c - d]`, the
  convention of `predicates::orient3d`.
- `tri_tri`: `Touch` means the intersection lies in both boundaries; `Cross`
  means some common point is in the relative interior of either triangle (a
  corner resting on the other's interior is `Cross`). Both triangles must be
  non-degenerate. `seg_tri` with `Coplanar` returns no point.
- Validation pairs: two corners count as shared when they have the same index or
  bitwise equal coordinates, so a triangle soup reports duplicate vertices and
  boundary edges but no self-intersections. Coplanar overlaps are reported for
  any pair except two faces with the same corner indices (`DuplicateFace`
  covers those), so coincident faces of touching solids are caught.
- `wireframes_core::fixtures` (public, used by every crate's tests): `Rng`
  (SplitMix64), `box_mesh`, `cube`, `tetrahedron`, `sphere(subdivisions)`,
  `torus(major, minor)`, `flipped`, `translated`, `nested_shells`, `soup`,
  `open_cube`, `pit_surface(n, seed)`, `random_points`, `coplanar_pairs` (the
  four section 12 boolean cases), and `broken::{degenerate_face, duplicate_face,
  duplicate_vertex, unused_vertex, non_manifold_fin, bowtie, edge_overlapping_shells,
  flipped_face, inward_cube, overlapping_cubes, touching_cubes, folded_face}`.
  Valid solids with shell contacts: `touching_tetrahedra` (one shared vertex index),
  `edge_touching_cubes(welded)`, `vertex_touching_cubes` and
  `vertex_touching_cavity` (a tetrahedron with its faces split at their centroids
  and the inward cavity spanned by those four centroids).
- Not done in core: the semi-static floating filter in front of the `exact`
  predicates (they always evaluate with expansions; results are exact, speed
  not measured).

## 15. Deferred (not in this contract)

Gmsh `.msh` and Medit `.mesh` I/O (needed before M7, see decision 17); loft
branching; true Minkowski offset; self-intersection repair; mesh-mesh intersection
lines as Polylines (falls out of the arrangement; add as `intersection_lines` in a
later revision); drape of polylines onto a surface; polygon booleans on
Polylines; polyline simplify, resample, buffer; GOCAD, OMF, LAS/LAZ, RESQML,
GRDECL, ZMAP readers; tetrahedralization; marching tetrahedra (stays in boitata
until it can return a wireframes `Mesh`).
