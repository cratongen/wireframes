# Containers and tables

`wireframes` has three geometry containers and one table type.

| Type | Holds | Row tables |
|---|---|---|
| `wf.Mesh` | triangles: an open surface (topography, fault, pit shell) or a closed solid (orebody, stope) | `vertex_attributes`, `face_attributes` |
| `wf.Polylines` | 3D strings: open lines and closed rings, grouped into features | `attributes` (per feature), `vertex_attributes` |
| `wf.PointSet` | points | `attributes` |
| `wf.Table` | typed columns with validity masks and units | |

`wf.UnstructuredGrid` holds the volume cells of a `.vtu` or `.vtk` file so they survive a round trip; it is data,
and `grid.surface()` turns its surface cells into a `Mesh`. Its members are `points`, `cell_types` (VTK type ids),
`connectivity` with `offsets` (`n_cells + 1` entries), `n_cells`, and the `point_data` and `cell_data` tables.

## Shared accessors

Every container has `.coords` as an `(n, 3)` float64 array, `.x`, `.y`, `.z`, `.bounds` as `(min, max)` or
`None`, `.column_names`, `.crs`, `.length_unit`, `to_table()`, and `to_pandas()`, `to_polars()` and
`to_pyarrow()`, which import their library when you call them.

Objects never change. `with_vertex_column`, `with_face_column`, `with_vertices`, `filter_faces`, `translate`
and the rest return new objects, and the arrays you get back are copies.

```python
import wireframes as wf
import numpy as np

vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
triangles = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
mesh = wf.Mesh(
    vertices, triangles, face_attributes={"rock": ["ox"] * len(triangles)}, crs="EPSG:31983", length_unit="m"
)
mesh = mesh.with_vertex_column("grade", np.random.default_rng(0).random(mesh.n_vertices), unit="g/t")
mesh["rock"]  # a vertex column, else a face column
mesh.to_table(on="faces")  # v0, v1, v2 and the face columns
```

A `Mesh` winds a solid's triangles counter-clockwise seen from outside, so its signed volume is positive.
`mesh.with_attributes(vertex=..., face=...)` replaces whole attribute tables (`None` keeps one), and
`wf.concat([a, b])` appends meshes in order without welding; a column missing from one mesh is null there.

## Polylines

A part is one line or ring; a ring does not repeat its first vertex. A feature groups parts, so a polygon with a
hole is one feature with two rings. Inside is even-odd over a feature's rings, and a point on an edge counts as
inside.

```python
outer = [[0, 0, 100], [40, 0, 100], [40, 30, 100], [0, 30, 100]]
hole = [[10, 10, 100], [20, 10, 100], [20, 20, 100], [10, 20, 100]]
points = [[5, 5, 100], [15, 15, 100]]
table = wf.Table({"ID": [1, 1, 1], "X": [0.0, 10, 10], "Y": [0.0, 0, 10], "Z": [100.0, 100, 100]})

pad = wf.Polylines([outer, hole], closed=True, features=[0, 0], attributes={"name": ["pad 3"]})
pad.area(), pad.length()
pad.contains(points)  # exact, in plan
pad.distance(points, signed=True)  # negative inside
wf.Polylines.from_table(table, feature="ID", x="X", y="Y", z="Z")
pad.to_points()  # every vertex, with feature, part and vertex columns
```

## Tables

```python
t = wf.Table({"hole": ["DH1", None], "depth": [120.0, np.nan], "n": [3, None]})
t["depth"]  # float64 with NaN for null
t["n"], t.mask("n")  # int64 payload and validity
t.with_column("au", [1.2, 0.4], unit="g/t").units
t.with_columns({"lab": ["A", "B"], "batch": [1, 2]})
t.drop_null()
```

Floats map to float64, integers to int64, booleans to bool and strings to utf8; `None` (and NaN in a float
column) marks a null.

## Protocol

`Mesh`, `Polylines`, `PointSet` and `Table` implement `__minetypes__`, protocol version 1. The dict holds
plain values and NumPy arrays. A stored float64, int64, bool, coordinate or triangle buffer is a read-only
view. A column with a validity mask is copied: a float null becomes NaN, an int null widens the column to
float64, and a bool or text null becomes int32 codes with `-1` for the missing row. There is no mask array
in the dict. `from_minetypes` copies a dict, or any object whose `__minetypes__` returns one, back into a
wireframes object. An open part needs at least two vertices; one vertex belongs in a `PointSet`.

```python
box = wf.datasets.box([0, 0, 0], [1, 1, 1])
spec = box.__minetypes__()
assert spec["kind"] == "mesh" and spec["triangles"].dtype.name == "uint32"
assert wf.Mesh.from_minetypes(spec).volume == 1.0
column = wf.Table({"au": [1.0, None]}).__minetypes__()["columns"]["au"]
assert column[0] == 1.0 and np.isnan(column[1])
```

## Labels

`crs` and `length_unit` are labels. Operations that combine two geometries check that declared labels match and
raise `InvalidInput` when they differ; an undeclared label matches anything. `mesh.to_length_unit("ft")`
rescales coordinates, and CRS transforms stay outside `wireframes`.

## Errors

| Class | Also a | Raised for |
|---|---|---|
| `wf.WireframesError` | `Exception` | base of every wireframes error |
| `wf.InvalidInput` | `ValueError` | bad shapes or values, unsupported content, CRS or unit mismatch |
| `wf.InvalidMesh` | `InvalidInput` | a mesh that fails validation; `.report` holds the `MeshReport` |
| `wf.KernelError` | `RuntimeError` | a kernel limit on valid input, such as a boolean result that rounding to f64 left self-intersecting; `.report` and `.mesh` hold the failed result, or `None` |
| `wf.FileError` | `OSError` | unreadable files and malformed content |
| `wf.MissingColumn` | `KeyError` | an unknown column name; the message lists the columns |

## Example

![Topography, orebody outline and drillholes](../examples/00-basics/basics.png)

The full example, with the output of each step: [example 00](../examples/00-basics/README.md).
