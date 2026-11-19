# Wireframes

## Purpose

Wireframes is a Rust-core, Python-exposed library to read, write, plot, measure, and
combine solids, surfaces, and strings for mining, geology, and petroleum. It has
zero dependency on `boitata`; once mature, `boitata` will depend on `wireframes` and
delete its own mesh type and mesh readers. The primary user is the author for
resource modeling and mine planning work, and the package is published as open
source under the MIT license.

## Commands

Pinned via `mise.toml` (Python, uv) and `rust-toolchain.toml` (Rust):

```sh
mise install        # one-time toolchain setup
mise run build      # cargo build --workspace
mise run test       # rust tests, build python module with maturin, run pytest
mise run lint       # rustfmt --check, clippy -D warnings, ruff check + format --check
mise run verify     # lint + test; must pass before finishing any task
mise run bench      # criterion benchmarks; pass a filter to run a subset
mise run oracle     # compare with Manifold, trimesh and scipy (dev-only oracles)
mise run examples   # rerun examples/NN-*/example_NN.py, rewrite their README.md pages and figures
mise run docs       # mkdocs serve
mise run docs-build # mkdocs build --strict
mise run compat     # install check across pip, uv, pixi, and conda
```

## Layout

- `crates/` - Rust workspace:
  - `core`: `Mesh`, `Polylines`, `PointSet`, attribute tables, validation report,
    and exact geometric predicates (`robust` crate).
  - `io`: readers and writers for OBJ, PLY, OFF, STL, DXF, VTK legacy, VTU,
    VTP, and Datamine `_tr.csv`/`_pt.csv` pairs.
  - `boolean`: native mesh arrangements, exact intersection predicates, implicit
    intersection points, boolean operations, plane slicing, and surface cutting.
  - `query`: point-in-solid, signed distance, ray/drillhole intersections, and
    exact block volume fractions inside solids.
  - `ops`: measures (volume, area, centroid, thickness, slices), repair,
    connected components, and transforms.
  - `build`: construction and editing: constrained Delaunay DTM with breaklines,
    lofting solids from section strings, extrusion, contours, decimation,
    smoothing, offset, and remeshing.
  - `py`: `pyo3` (`abi3-py311`) bindings built with `maturin`.
- `python/wireframes/` - Python package imported as `import wireframes as wf` (thin
  wrappers, `_wireframes.pyi` stub, `plot.py`, `errors.py`, `_containers.py`).
- `python/tests/` - pytest suite and property tests.
- `examples/` - numbered script examples rendered into documentation.
- `docs/` - `mkdocs` documentation site.

## Conventions

- Follow Conventional Commits (`feat:`, `fix:`, `docs:`, `refactor:`, `test:`,
  `perf:`, `chore:`), one logical change per PR, and seeded synthetic or public
  datasets only.
- Every feature ships together: Rust implementation, `pyo3` binding, `.pyi`
  stub, NumPy-style docstring, and an example. Every geometric algorithm also
  ships with a rerunnable invariant check (for example, the boolean volume
  identity `vol(A union B) + vol(A intersect B) = vol(A) + vol(B)`).
- Match `boitata`'s Python API style: flat top-level exports, `.x`/`.y`/`.z`
  coordinate views, `.column_names`, `to_polars`/`to_pandas`/`to_pyarrow` on
  tabular containers, `wf.plot` functions accepting `ax=None` and returning
  `(fig, ax)`, `numpy` as the sole required runtime dependency, and `matplotlib`
  gated behind `wireframes[plot]`.
- Booleans refuse invalid meshes by default and raise `InvalidMesh` carrying the
  `MeshReport`. Repair is opt-in via `mesh.repair()` or `repair=True`, and every
  modification is recorded in the returned report.
- Results must not depend on thread count. Parallel reductions use deterministic
  ordering.
- Impose no fixed RAM caps in algorithms or spatial indexes; size working memory
  from available system RAM.
- No code comments except where an invariant, coordinate convention, or
  predicate degeneracy rule is non-obvious from the code.

## Domain glossary

- **Solid**: one or more closed, manifold, outward-oriented triangle shells
  (`is_closed` and zero self-intersections, non-manifold edges, inconsistent
  windings, or inward shells) representing a 3D geological domain, stope, or
  excavation. Shells may touch each other along edges or at vertices when their
  interiors do not overlap; `validate` counts those contacts
  (`touching_shell_edges`, `touching_shell_vertices`) as information, not
  defects. A shell may touch itself along edges or at vertices on the same
  terms.
- **Surface / DTM**: open triangle mesh representing topography, a fault, an
  unconformity, a pit shell, or a stratigraphic horizon.
- **String**: a 3D polyline (`Polylines`) representing pit crests and toes,
  grade-control polygons, fault traces, or geological section interpretations.
- **Breakline**: a 3D string enforced as constrained edges during Delaunay
  triangulation of a surface.
- **Arrangement**: the subdivision of triangle faces induced by exact pairwise
  intersections before classifying cells as inside or outside.
- **Cut and fill**: net and gross volumes trapped above and below a reference
  surface when compared against a design or updated topography surface.

## Never do

- Never add AI attribution or co-author trailers to commits or PR descriptions
  (no `Co-Authored-By: Claude...`, no "Generated with Claude Code" footers).
- Never use floating-point tolerances to decide topological orientation or
  side-of-plane tests in booleans; use exact predicates (`robust`) and implicit
  intersection points.
- Never silently repair or alter input meshes inside boolean, cutting, or volume
  functions.
- Never add the C++ `Manifold` library as a build or runtime dependency; use it
  strictly as an external test and benchmark oracle.
- Never support closed proprietary formats (DWG, binary Datamine `.dm`) or commit
  client data or files over 1 MB.
- Never guess or extrapolate performance numbers; write `"not measured"` until a
  benchmark is run.
- Never push or force-push without asking first, or call a code task done
  without running `mise run verify`.

## Definition of done

- `mise run verify` passes (once Milestone 1 lands the build skeleton).
- New or changed Rust functionality includes its `pyo3` binding, `.pyi` stub,
  NumPy-style docstring, invariant test, and example.
- `mise run docs-build` succeeds when documentation or examples change.
