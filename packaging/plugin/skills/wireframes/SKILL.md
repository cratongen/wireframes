---
name: wireframes
description: How to write correct code with wireframes (`import wireframes as wf`): solids, surfaces and strings, files, booleans, measures and block fractions. Use whenever code imports wireframes, or the user asks for solids, surfaces or strings in Python and wireframes is installed or wanted.
---

# wireframes

wireframes is imported as `wf`. It is not on PyPI. From a checkout of
https://github.com/cratongen/wireframes, install it with `mise run py:build`. Python >= 3.11. numpy is the
only required dependency. `wf.plot` needs matplotlib.

```python
import wireframes as wf
```

## Look it up before you write it

Do not write a wireframes call from memory. The installed package carries its own reference, and it matches the
installed version exactly:

```sh
python -c "import wireframes, pathlib; print(wireframes.__version__, pathlib.Path(wireframes.__file__).with_name('_wireframes.pyi'))"
```

- Search that `_wireframes.pyi` stub for `Mesh`, `Polylines`, `PointSet`, `Table`, `UnstructuredGrid` and the
  module functions. `datasets`, `plot`, `read` and the error classes live in `datasets.py`, `plot.py`,
  `_containers.py` and `errors.py` next to the stub.
- `python -c "import wireframes as wf; help(wf.Mesh.block_fraction)"` prints the NumPy-style docstring, with
  every argument and its default.
- `wf.__all__` lists every public name. `dir(wf.datasets)` and `dir(wf.plot)` list the datasets and the plots.
- The docs at https://cratongen.github.io/wireframes/ are the guide: containers, files, measure and repair,
  booleans, queries, construction, plotting, and the performance page. In a checkout the same pages are under
  `docs/`.

If a name is not in `wf.__all__`, it does not exist. Say so and offer the closest one that does.

## Kinds of object

Ask which kind you hold, and the calls follow.

| Kind | Examples | You call | Saves with |
| --- | --- | --- | --- |
| Container | `Mesh`, `Polylines`, `PointSet`, `Table` | `with_vertex_column`, `with_face_column`, `with_column`, `filter`, `filter_faces`, `coords` | `wf.write` |
| Grid data | `UnstructuredGrid` | `surface`, `points`, `cell_data` | `wf.write` |
| Report | `MeshReport`, `RepairReport`, `RoundingReport` | `summary`, `is_solid`, `volume_bound` | keep the object |
| Function | `union`, `dtm`, `loft`, `extrude`, `cut_fill`, `block_fraction` | one call on the whole input | its result |

`UnstructuredGrid` has no constructor. `wf.read` of a `.vtu` or `.vtk` file returns one on `Dataset.grid`, and
`grid.surface()` is the triangle `Mesh`.

One run through the kinds you can build directly:

```python
import numpy as np
import wireframes as wf

solid = wf.datasets.box([0, 0, 0], [2, 2, 2])
other = wf.datasets.box([1, 1, 1], [3, 3, 3])
both = wf.union(solid, other)
report = solid.validate()
points = wf.PointSet([[1, 1, 1], [2, 1, 1], [5, 5, 5]])
where = solid.locate(points)
fractions = solid.block_fraction([[2, 1, 1], [1, 1, 1]], [1, 1, 1])
ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=True)
prism = wf.extrude(ring, [0, 0, 2])
base_ring = [[0, 0, 0], [20, 0, 0], [20, 10, 0], [0, 10, 0]]
top_ring = [[5, 2.5, 10], [15, 2.5, 10], [15, 7.5, 10], [5, 7.5, 10]]
frustum = wf.loft(wf.Polylines([base_ring, top_ring], closed=True, features=[0, 1]))
ground = wf.dtm([[0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]])
volumes = wf.cut_fill(ground, ground.translate([0, 0, -2]))
table = wf.Table({"depth": [12.0, np.nan]}).drop_null()
fig, ax = wf.plot.plan(solid)
```

`where` is int8: -1 inside, 0 on the boundary, 1 outside. `fractions` is the exact fraction of each block.

## Rules that hold everywhere

- **Options are keyword-only.** `wf.Polylines(parts, True)` raises `TypeError`. Write `closed=True`.
  `Mesh` takes vertices and triangles first, then `vertex_attributes`, `face_attributes`, `crs` and
  `length_unit` by keyword.
- **Containers never change in place.** `with_vertex_column`, `with_face_column`, `with_column`, `filter`,
  `filter_faces`, `translate` and `repair` return a new object. Assign it.
- **`Mesh` has no `with_column` and no `filter`.** Vertex columns go through `with_vertex_column`, face
  columns through `with_face_column`, and a subset of triangles through `filter_faces`. `mesh["grade"]` is
  the vertex column of that name, or the face column when there is no vertex column.
- **`Polylines`, `PointSet` and `Table` take `with_column` and `filter`.** A polyline feature groups parts.
  A ring does not repeat its first vertex. `closed` is one flag or one flag per part.
- **Booleans need solids.** `a | b`, `a & b`, `a - b` and `a ^ b` are `wf.union`, `wf.intersection`,
  `wf.difference` and `wf.symmetric_difference`. An input that is not a solid raises `wf.InvalidMesh` with
  the report. `repair=True` repairs both inputs first and returns the repair reports. `rounding_report=True`
  adds a `wf.RoundingReport` last.
- **Shells may touch along an edge or at a vertex.** They may not share a face, part of an edge, or a vertex
  sitting inside another shell's face. `validate` reports those as coplanar overlaps or self-intersections,
  and the boolean raises `wf.InvalidMesh`. Join face-sharing solids with `union`, or keep them as separate
  meshes.
- **A result the rounding cannot certify raises `wf.KernelError`.** `error.report` is the `MeshReport` and
  `error.mesh` is the rounded mesh, when there is one. Leave that mesh out of the answer you report.
- **`volume` does not look for self-intersections.** Overlapping shells count the overlap twice. Read
  `report.is_solid` from one `validate()` before you report a volume.
- **Axes.** Coordinates are easting, northing and elevation, z up. A solid's triangles wind counter-clockwise
  seen from outside, so the signed volume is positive. `dip` is `(m, 2)` dip and dip direction in degrees,
  dip direction clockwise from north.
- **`crs` and `length_unit` are labels.** A call that combines two objects raises `wf.InvalidInput` when both
  declare a label and the labels differ. An undeclared label matches either. `to_length_unit` rescales.
  wireframes does not transform a CRS.
- **`loft` links one ring to one ring.** A section with two or more rings raises `wf.InvalidInput`. The first
  and last sections may be a line or a point. `extrude` sweeps rings along a vector. `offset(distance)` moves
  each vertex by `distance` along its angle-weighted normal, so a face at a box corner moves less than
  `distance`.
- **`smooth(preserve_volume=True)` rescales a closed solid** so its volume matches the input. It is off by
  default. `decimate` and `remesh` return a new mesh.
- **Results do not depend on the number of threads.**
- **These are not scikit-learn estimators.** There is no `fit` and no `Pipeline`.

## Missing values

Float nulls are `NaN`. Text and integer nulls are `None`, and `table.mask(name)` is false on those rows.
The Datamine reader turns `-1.0E30` into a null. `drop_null` drops the rows where a named column, or every
column, is null.

`locate`, `contains`, `signed_distance`, `distance`, `closest_point`, `block_fraction`, `ray_hits` and
`segment_hits` raise `wf.InvalidInput` when a point is not finite. Drop those rows before the call.

```python
import numpy as np
import wireframes as wf

assayed = wf.Table({"depth": [12.0, np.nan], "au": [1.2, 0.4]}).drop_null("depth")
```

## Errors

wireframes errors derive from `wf.WireframesError` and from a builtin.

| Class | Also a | When |
| --- | --- | --- |
| `wf.InvalidInput` | `ValueError` | bad shape or value, unsupported column, CRS or unit mismatch |
| `wf.InvalidMesh` | `wf.InvalidInput` | validation failed; `.report` is the `MeshReport` |
| `wf.KernelError` | `RuntimeError` | the kernel stopped on valid input; `.report` and `.mesh` |
| `wf.FileError` | `OSError` | unreadable file or malformed content |
| `wf.MissingColumn` | `KeyError` | the message lists the columns the object has |

Read the message before changing the call. It names the argument.

## Saving

`wf.write(path, data)` picks the format from the extension. `data` is a `Mesh`, `Polylines`, `PointSet`,
`UnstructuredGrid` or `Dataset`. `wf.read(path)` returns a `wf.Dataset` with `mesh`, `polylines`, `points`,
`grid`, `length_unit` and `warnings`. `wf.read_mesh`, `wf.read_polylines` and `wf.read_points` return that
one part and raise `wf.InvalidInput` when the file has none.

There is no `to_json` on these objects. Readers do not repair and do not rescale.

`write` raises `wf.InvalidInput` when a column cannot be stored. `skip_unsupported_columns=True` drops those
columns. `ascii=True` writes text PLY, STL and VTK. `dxf_entity` is `"3dface"` or `"polyface"`. `zlib=True`
compresses VTU arrays.

PLY, VTK and VTU keep numeric columns. Datamine `_tr.csv` / `_pt.csv` keeps text and numeric columns. DXF
keeps face columns `layer` and `color`. STL stores 32-bit floats: at a northing of 7.5 million their spacing
is half a metre, so keep STL for local coordinates.

```python
import tempfile
from pathlib import Path

import wireframes as wf

solid = wf.datasets.box([0, 0, 0], [2, 2, 2])
path = Path(tempfile.mkdtemp()) / "box.ply"
wf.write(path, solid)
back = wf.read_mesh(path)
```

## Data

`wf.datasets` builds every sample from a seed, offline. Use one instead of inventing coordinates.
`box(lo, hi)`, `orebody(seed=0)`, `topography(seed=0)`, `pit(ground)` and `drillholes(seed=0)`. `box` and
`orebody` are solids. `topography` and `pit` are surfaces. `drillholes` is open `Polylines` with feature
columns `hole` and `depth`. Coordinates sit at UTM scale around `wf.datasets.ORIGIN`.

Grouping and joining belong to Polars or pandas: `mesh.to_polars()`, `mesh.to_pandas()`, `mesh.to_pyarrow()`.
Pass `on="faces"` for one row per triangle.

## Write it fast

Correct wireframes code still does the expensive work more than once when it loops in Python or rebuilds the
mesh between calls. The performance page records the measured times. The paths below are what the code does.

### Validate once, then reuse the mesh

`locate`, `contains`, `signed_distance`, `trace_intervals` and `block_fraction` need a solid. The first of
those calls on a mesh runs `validate` and stores the index on that object. Later calls on the same object
reuse the index.

`mesh.is_solid` and `mesh.validate()` run the checks and do not store the index. A loop that reads
`is_solid` pays them on every pass. Call `report = mesh.validate()` once and read `report.is_solid` and
`report.summary`.

`with_vertex_column`, `with_face_column`, `translate`, `filter_faces`, `repair` and every boolean return a
new mesh. The index stays on the old one.

A boolean validates each input and its result. Calling `validate()` yourself just before the boolean repeats
that work. `repair=True` repairs both inputs on that call. Pass it when the report says the mesh needs
repair, not on a solid you already checked.

### One call for the whole array

`block_fraction(centroids, size)` takes every centroid in one `(n, 3)` array and walks the blocks in
parallel. One call for the model. A Python loop that calls it once per block is the slow shape, and a new
mesh on each pass validates again.

The same holds for `locate`, `contains`, `signed_distance`, `distance`, `closest_point`, `ray_hits` and
`segment_hits`: pass the whole array. `trace_intervals(trace, depths)` takes the whole trace.

`block_fraction` takes one `size` for every centroid in the call. On a regular or masked boitata model that
size is `model.size`, and only when `model.extents is None` and `model.rotation` is `(0, 0, 0)`. A
sub-blocked model stores a different box on each row. This call has no size per row, so group the rows that
share a size and call once per group.

### Read once, write when you leave

`wf.read` or `wf.read_mesh` once, keep the object, and `wf.write` when the result leaves the process. A
write followed by a read builds a new mesh and drops the index. STL and DXF also drop columns or digits, so
the object you get back is a different mesh.

### Pieces that do not touch

`wf.concat` appends meshes and does not weld. One boolean on that mesh validates it once. A Python loop of
unions validates once per piece. The booleans guide puts separated crosscuts in one mesh, then unions once,
and computes each reported volume from the original solids:

```python
import wireframes as wf

drive = wf.datasets.box([0, 0, 0], [100, 5, 4])
crosscuts = [wf.datasets.box([x, 2, 0.5], [x + 4, 30, 3.5]) for x in range(10, 100, 20)]
development = drive | wf.concat(crosscuts)
```

### Do not promise a speed

Give the user a time or a speedup only after you ran the code and timed it. Otherwise say it was not
measured. Numbers that are already on the performance page stay there, with the machine and the command.

## Hand it to a sibling

wireframes is one of a family of standalone libraries for mining and geoscience at https://github.com/cratongen.
None imports another. They pass NumPy arrays, and tables pass through Arrow, so pandas, polars and pyarrow
work in between.

| Library | Import | What it does |
| --- | --- | --- |
| boitata | `import boitata as bt` | Geostatistics: drill holes, compositing, block models, variograms, kriging, simulation |
| segyrs | `import segyrs` | Read, write and slice SEG-Y files, including files larger than memory |
| transiograms | `import transiograms as tg` | Markov chains, hidden Markov models, transiograms and categorical simulation |
| wireframes | `import wireframes as wf` | Solids, surfaces and strings: read, write, measure and combine |

<!-- Below this line: what wireframes hands to each sibling and takes from it, with one short example per neighbour.
Only calls that exist today in both libraries, checked against ../<sibling>/python/<pkg>/*.pyi. -->

boitata has its own `Mesh`, `Polylines` and `BlockModel`. Pass arrays. `bt.Mesh` takes vertices and triangles.
`bt.Polylines` takes a sequence of parts. A face column does not cross with the mesh: copy the array with
`with_face_column`. `bt.BlockModel.coords` are the cell centres and `bt.BlockModel.size` is the block size.
`origin` is the minimum corner of cell `(0, 0, 0)`.

```python
import boitata as bt
import wireframes as wf

solid = wf.datasets.box([0, 0, 0], [10, 20, 5], name="stope")
mesh = bt.Mesh(solid.coords, solid.triangles, length_unit=solid.length_unit)
mesh = mesh.with_face_column("name", solid.face_attributes["name"])

crest = wf.Polylines([[[0, 0, 100], [40, 0, 100], [40, 30, 100], [0, 30, 100]]], closed=True, length_unit="m")
strings = bt.Polylines(crest.parts, closed=True, length_unit=crest.length_unit)

model = bt.BlockModel(origin=(0.0, 0.0, 0.0), size=(5.0, 10.0, 5.0), count=(2, 2, 1), length_unit="m")
fractions = solid.block_fraction(model.coords, model.size)
model = model.with_column("inside", fractions)
```

Use that `block_fraction` call when `model.extents is None` and `model.rotation` is `(0, 0, 0)`. boitata
stores a rotated model's angles as `(azimuth, dip, rake)` and does not export the `(3, 3)` axis matrix
`block_fraction` takes in `rotation=`. For a rotated model, build that matrix yourself. Its columns are the
block axes in world coordinates, orthonormal, determinant +1.

`transiograms.Grid` is an axis-aligned grid. `grid.coords` are the cell centroids and `grid.size` is the cell
size, so the same fraction call takes them:

```python
import transiograms as tg
import wireframes as wf

solid = wf.datasets.box([0, 0, 0], [10, 20, 5])
grid = tg.Grid((0.0, 0.0, 0.0), (5.0, 10.0, 5.0), (2, 2, 1))
fractions = solid.block_fraction(grid.coords, grid.size)
```

segyrs reads and writes SEG-Y. Its calls take a segyrs file. wireframes has no SEG-Y reader, and segyrs has
no mesh or string type. There is no handoff.

## Before you report a result

Run the code. Then check the identity the booleans promise, and the block-fraction sum the queries promise.
`vol(a | b) + vol(a & b)` equals `vol(a) + vol(b)` within rounding. Over blocks that cover a solid,
fraction times block volume sums to `volume`. Report a volume only when `validate().is_solid` is true, and
report the numbers you got, in the data's length unit.

```python
import wireframes as wf

a = wf.datasets.box([0, 0, 0], [2, 2, 2])
b = wf.datasets.box([1, 1, 1], [3, 3, 3])
gap = (a | b).volume + (a & b).volume - a.volume - b.volume
cover = a.block_fraction([[1, 1, 1]], [2, 2, 2])
balance = cover[0] * 8.0 - a.volume
```
