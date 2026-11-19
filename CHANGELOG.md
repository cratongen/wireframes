# Changelog

All notable changes to wireframes are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

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
  `close_to_elevation`, `contours`, `decimate`, `smooth`, `offset`, `remesh`.
- Python package `wireframes` with a full `.pyi` stub, seeded datasets, matplotlib
  plots and `to_polars`/`to_pandas`/`to_pyarrow` conversions.
- Numbered examples for every public capability, and an mkdocs guide.
- Criterion benches for I/O, booleans, slicing, validation and queries (not
  timed yet).
- Tests: boitata's validate and repair fixtures ported to `wireframes-ops`,
  thread-count determinism tests for `io` and `build`, commutativity of union
  and intersection in `locate`, a Datamine triangle-id alias test, and
  DataFrame conversion tests that check null handling.
- `mise run compat`, which installs the wheel with pip, uv, pixi and conda.
- Solids of touching shells: `validate` counts contacts between closed shells
  in `touching_shell_edges` and `touching_shell_vertices` and keeps them
  solid; booleans and `clip_by_surface` accept and return such solids, so
  results chain. `repair` keeps touching shells apart (`vertex_split` actions)
  unless `weld_touching_shells=True`.

### Changed

- Test names now match the README milestone acceptance filters.
- `symmetric_difference` of crossing solids is a valid solid: `a - b` and
  `b - a` are separate shells touching along the intersection curve.
- Union of solids that share only an edge or a vertex, and differences that
  leave bodies joined at an edge or vertex, return separate shells instead of
  raising `InvalidInput`.
- `components` splits shells that touch along an edge.

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
- `split_by_plane` no longer panics when a closed mesh meets the plane only at
  a vertex.
- `split_by_plane` triangulates the cap from the rounded section points, so
  cap triangles no longer overlap where the section has long, nearly collinear
  runs of points (a bench split of a clipped pit).
- Boolean and `clip_by_surface` results no longer keep degenerate faces,
  duplicate vertices or open edges where distinct exact points round to the
  same float: each output shell merges such vertices and removes faces with
  repeated or collinear corners.
- `block_fraction` is exactly 0 or 1 for a block that no face meets, instead
  of picking up `f64` residue from the faces above it.

### Known issues

- A boolean result shell that touches itself raises `InvalidInput`.
- A boolean result can fail validation after its exact vertices are rounded to
  f64: rounding can fold a sliver thinner than the float spacing onto a
  neighbouring face. Seen in chains of grid-box booleans and in booleans of
  pits clipped by the same topography (their caps differ by about one ulp);
  for pushbacks, clip the difference of the unclipped pits instead.
- `split_by_plane` at a plane the mesh touches from above (a DTM contour
  level, for example) returns a positive part that touches itself and is not
  a solid.
- A cavity that touches its outer shell along an edge is refused as a
  self-intersection, while one that touches it at a vertex is a valid contact.
- `loft` does not branch, `extrude` has no path sweep, and `offset` is not a
  Minkowski offset.
- No real Datamine Studio export is tested yet.
