# wireframes

`wireframes` reads, writes, plots, measures and combines 3D solids, surfaces and strings for mining, geology and
petroleum. A Rust core does the geometry with exact predicates; you call it from Python as `import wireframes as wf`.

The pre-release `0.1.0a1` is on PyPI. pip installs a pre-release only when you ask for one:

```sh
pip install --pre wireframes           # numpy is the only required dependency
pip install --pre "wireframes[plot]"   # adds matplotlib, for wf.plot
```

PyPI has wheels for Linux x86_64 and aarch64 (glibc 2.28 or newer), macOS (Intel and Apple silicon) and
Windows x64, for CPython 3.11 and newer. On another platform pip builds wireframes from the sdist, which needs a
Rust toolchain.

## A first look

```python
import wireframes as wf

orebody = wf.datasets.orebody()  # seeded synthetic lens at UTM scale
block = wf.datasets.box(wf.datasets.ORIGIN + [-50, -80, 150], wf.datasets.ORIGIN + [50, 80, 250])

print(orebody.volume, orebody.is_solid)
inside = orebody & block  # exact boolean intersection
print(inside.volume + (orebody - block).volume - orebody.volume)  # 0 to rounding

wf.write("lens.dxf", orebody, skip_unsupported_columns=True)
fig, ax = wf.plot.section([orebody, inside], wf.datasets.ORIGIN, [0, 1, 0])
```

## What the guide covers

| Page | Capability |
|---|---|
| [Containers and tables](guide/containers.md) | `Mesh`, `Polylines`, `PointSet`, `Table`, attributes, labels, errors |
| [Files](guide/files.md) | OBJ, PLY, OFF, STL, DXF, VTK, VTU/VTP, Datamine CSV pairs |
| [Measure, validate and repair](guide/measure-repair.md) | volume, area, slabs, thickness, `validate`, `repair`, transforms |
| [Booleans and cutting](guide/booleans.md) | `a \| b`, `a & b`, `a - b`, `a ^ b`, surface and plane cuts, cut and fill |
| [Queries](guide/queries.md) | point in solid, distances, rays, drillhole intervals, block fractions |
| [Construction](guide/construction.md) | DTM, polygons, extrusion, lofting, closing, decimation, smoothing, remeshing |
| [Plotting](guide/plotting.md) | plan, section and 3D views |
| [Open pit, worked example](examples/06-open-pit/README.md) | topography from contours, pit stages, bench tonnes, survey against design |
| [Underground, worked example](examples/07-underground/README.md) | orebody from section strings, fault blocks, stopes and development, drillhole intercepts, block tonnage report |
| [Files and repair, worked example](examples/08-files-repair/README.md) | Datamine pair and DXF as other packages write them, validation and repair reports, format conversion with attributes, the rules for a valid solid |

Each page ends with a numbered example from `examples/`. `mise run examples` reruns them all and rewrites their
pages and figures; every example asserts the invariant it shows.

## Rules the library keeps

- Topological decisions use exact predicates; no epsilon decides inside, outside or crossing.
- Booleans, cuts and volume functions refuse invalid meshes and raise `wf.InvalidMesh` with the full report.
  Repair runs only when you call `mesh.repair()` or pass `repair=True`, and it records every change.
- Results do not depend on the number of threads.
- Every object is immutable; methods return new objects.

The [performance page](performance.md) gives the criterion medians, with the machine,
versions and commands behind each number.
