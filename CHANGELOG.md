# Changelog

All notable changes to drusa are listed here. The format follows
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
- Python package `drusa` with a full `.pyi` stub, seeded datasets, matplotlib
  plots and `to_polars`/`to_pandas`/`to_pyarrow` conversions.
- Numbered examples for every public capability, and an mkdocs guide.
- Criterion benches for I/O, booleans, slicing, validation and queries (not
  timed yet).
- Tests: boitata's validate and repair fixtures ported to `drusa-ops`,
  thread-count determinism tests for `io` and `build`, commutativity of union
  and intersection in `locate`, a Datamine triangle-id alias test, and
  DataFrame conversion tests that check null handling.
- `mise run compat`, which installs the wheel with pip, uv, pixi and conda.

### Changed

- Test names now match the README milestone acceptance filters.
- Spec 6.2 item 8 now states that `symmetric_difference` of crossing solids is
  not a solid.

### Removed

- Unused `dxf` and `regex` dependencies; DXF reading is a hand-written
  group-code reader.

### Known issues

- `symmetric_difference` of solids whose surfaces cross keeps duplicate
  vertices along the intersection curve.
- A split whose section is narrower than about 1e-154 loses its cap.
- `loft` does not branch, `extrude` has no path sweep, and `offset` is not a
  Minkowski offset.
- No real Datamine Studio export is tested yet.
