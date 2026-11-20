# Changelog

All notable changes to wireframes are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Booleans, `clip_by_surface` and `split_by_plane` round exact results with
  exact tests instead of refusing them (`docs/design/rounding.md`): vertices
  that round together merge, vertices within rounding of a face are inserted
  into it, flat runs and sub-ulp shells go, and inputs that agree only to
  rounding are conformed so their faces coincide exactly. Each result carries a
  `RoundingReport` (`wf.RoundingReport`, `rounding_report=True`) with the moved
  vertices, the move bound and a proven volume bound. The 188 refused inclined
  splits of `sphere(2)`, the pushback booleans, the rotated stopes and union 11
  of example 07 now give solids, and the six skipped Python acceptance tests and
  the ignored Rust one run.

- Rust workspace with `core` (mesh, polylines, point sets, typed tables, exact
  predicates, validation), `io`, `ops` (measures, repair, components,
  transforms, thickness), `boolean`, `query`, `build` and `py` crates.
- Readers and writers for OBJ, PLY, OFF, STL, DXF, VTK legacy, VTU/VTP and
  Datamine `_tr.csv`/`_pt.csv` pairs, with `UnstructuredGrid` for volume cells.
- Exact booleans (`union`, `intersection`, `difference`, `symmetric_difference`),
  `clip_by_surface`, `split_by_plane`, `slice_by_plane` and `cut_fill`.
- Queries: `contains`, `locate`, `signed_distance`, `closest_point`, `ray_hits`,
  `segment_hits`, `trace_intervals` and exact `block_fraction`.
- Construction: `dtm` with breaklines, `triangulate`, `extrude`, `loft`,
  `close_to_elevation`, `contours`, `decimate`, `smooth` (with `preserve_volume`
  to keep a solid's volume), `offset`, `remesh`.
- Python package `wireframes` with a full `.pyi` stub, seeded datasets, matplotlib
  plots and `to_polars`/`to_pandas`/`to_pyarrow` conversions.
- Numbered examples for every public capability, and an mkdocs guide.
- Criterion benches for I/O, booleans, slicing, validation and queries.
  `docs/performance.md` records the timings
  taken at commit `5d5ab30`; the cost of the `validate` call that booleans
  now run on each result is not measured.
- Tests: boitata's validate and repair fixtures ported to `wireframes-ops`,
  thread-count determinism tests for `io` and `build`, commutativity of union
  and intersection in `locate`, a Datamine triangle-id alias test, and
  DataFrame conversion tests that check null handling.
- `mise run compat`, which installs the wheel with pip, uv, pixi and conda.
- `wf.KernelError` (`Error::Kernel` in Rust), a `RuntimeError` for limits of
  the kernel on valid input. `report` and `mesh` hold the failed result when
  there is one.
- Solids of touching shells: `validate` counts contacts between closed shells
  in `touching_shell_edges` and `touching_shell_vertices` and keeps them
  solid; booleans and `clip_by_surface` accept and return such solids, so
  results chain. `repair` keeps touching shells apart (`vertex_split` actions)
  unless `weld_touching_shells=True`.

### Changed

- Booleans, `clip_by_surface` and `split_by_plane` of a solid return only
  solids. A result that rounding to f64 left self-intersecting raises
  `wf.KernelError` with its report and mesh; before, it came back closed but
  not solid, and only `validate` showed it.
- Kernel failures (an arrangement triangulation, a classification, rounding a
  result) raise `wf.KernelError` instead of `wf.InvalidInput` with an
  `internal:` message, so they no longer blame the input.
- Test names now match the README milestone acceptance filters.
- `loft` takes an open line or a single point as its first or last section,
  so a vein can pinch out to a wedge or a pyramid.
- An open polyline part may hold a single vertex (a point). OBJ `l`, DXF
  polylines and VTK/VTP line cells with one vertex read as such parts instead
  of being skipped or refused.
- `symmetric_difference` of crossing solids is a valid solid: `a - b` and
  `b - a` are separate shells touching along the intersection curve.
- Union of solids that share only an edge or a vertex, and differences that
  leave bodies joined at an edge or vertex, return separate shells instead of
  raising `InvalidInput`.
- `components` splits shells that touch along an edge.
- A shell may touch itself along an edge or at a vertex: `validate` counts the
  contact in `touching_shell_edges` and `touching_shell_vertices`, and booleans
  and `clip_by_surface` return such shells instead of raising `InvalidInput`.
  Unions of stopes no longer depend on their order.

### Removed

- Unused `dxf` and `regex` dependencies; DXF reading is a hand-written
  group-code reader.

### Fixed

- Exact arithmetic no longer underflows or overflows: an expansion operation
  that would leave the `f64` range continues on an integer times a power of
  two, and `orient2d`, `orient3d`, `incircle` and `insphere` fall back to it
  outside `[2^-120, 2^120]`. A split whose section is narrower than about
  1e-154 keeps its cap, and booleans, `clip_by_surface`, `split_by_plane` and
  `slice_by_plane` work on coordinates near 1e150 and in the subnormal range.
- `validate` and `repair` decide shell orientation by the exact sign of the
  enclosed volume. The `f64` volume underflowed to zero for tiny solids and
  overflowed for coordinates near 1e150, which flagged an outward shell as
  inward.
- `validate` and `repair` no longer take a cavity that touches its outer shell
  at every vertex it has for an inward shell; nesting falls back to an exact
  point inside one of its edges.
- `validate` accepts a cavity that touches its outer shell along an edge, and
  two cavities that touch each other along an edge, as contacts. It reported
  them as self-intersections because each wedge of solid around the edge runs
  from one shell to the other, while the same contact at a vertex was valid.
  An edge contact is now valid whenever the faces around it alternate in
  direction, so no wedge is inside twice.
- `split_by_plane` no longer panics when a closed mesh meets the plane only at
  a vertex.
- `split_by_plane` of a closed mesh at a plane it touches from above, such
  as a closed DTM split at the elevation of a local minimum, returns two
  solids. The cap passes through the contact vertex or edge, so the upper
  part touches itself there; before, it came back with 4 to 7
  self-intersections.
- `split_by_plane` triangulates the cap from the rounded section points, so
  cap triangles no longer overlap where the section has long, nearly collinear
  runs of points (a bench split of a clipped pit).
- Boolean and `clip_by_surface` results no longer keep degenerate faces,
  duplicate vertices or open edges where distinct exact points round to the
  same float: each output shell merges such vertices and removes faces with
  repeated or collinear corners.
- Booleans remove a sliver that rounding folds onto the coplanar face across
  its longest edge, where a stope wall cuts a drive wall quad that is planar
  only to rounding. In a seeded scan of 240 stope and drive booleans, 29
  results were invalid before and none after.
- Rounding no longer cycles on a flat run it removes and puts back: the union
  of two pit shells clipped by one topography ran for more than 17 minutes
  and now returns a solid. A run's removal is kept only when the faces it
  makes hold fewer flat faces than the run.
- Rounding certifies the adversarial cases that raised `uncertified`: a stope
  split along its own face plane (a sub-ulp sliver part is removed within its
  volume bound), repeated union then difference of one body, a rotated stope
  block minus a raise, development turned to strike, and chains of lids and
  caps tilted by 1e-13 rad. When the other operations stall, the rounding
  descends: it keeps the local flip, join, collapse or insertion that lowers
  the problem count most (`docs/design/rounding.md` section 9).
- `block_fraction` is exactly 0 or 1 for a block that no face meets, instead
  of picking up `f64` residue from the faces above it.
- `loft` links each band with both rings centred and the band axis removed.
  The 3D area cost could not tell a long band's straight walls from a double
  fan, so sections cut obliquely across the axis and mitred drive rings along
  a bent path came back refused or with the wrong volume.

### Known issues

- Rounding has no general 3D validity theorem: a boolean, clip or split the
  rounding operations cannot certify raises `wf.KernelError` whose message
  starts with `uncertified`. No case in the tests raises it, the adversarial
  cases of `crates/boolean/tests/attacks.rs` included.
- `loft` does not branch, `extrude` has no path sweep, and `offset` is not a
  Minkowski offset.
- No real Datamine Studio export is tested yet.
- Solids that share a wall (all or part of a face) are not valid in one mesh;
  join them with `union` first.
- No release, published wheels or CI exist yet; install by building the wheel
  from a clone.
