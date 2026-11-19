# wireframes

`wireframes` reads, writes, plots, measures and combines 3D solids, surfaces and strings for mining, geology and
petroleum. A Rust core does the geometry with exact predicates; you call it from Python as `import wireframes as wf`.

```sh
pip install wireframes            # numpy is the only required dependency
pip install "wireframes[plot]"    # adds matplotlib for wf.plot
```

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

Each page ends with a numbered example from `examples/`. `mise run examples` reruns them all and rewrites their
figures; every example asserts the invariant it shows.

## Rules the library keeps

- Topological decisions use exact predicates; no epsilon decides inside, outside or crossing.
- Booleans, cuts and volume functions refuse invalid meshes and raise `wf.InvalidMesh` with the full report.
  Repair runs only when you call `mesh.repair()` or pass `repair=True`, and it records every change.
- Results do not depend on the number of threads.
- Every object is immutable; methods return new objects.

Speed has not been measured against other packages yet.
