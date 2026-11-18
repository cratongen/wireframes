# Containers and tables

`drusa` has three geometry containers and one table type.

| Type | Holds | Row tables |
|---|---|---|
| `dr.Mesh` | triangles: an open surface (topography, fault, pit shell) or a closed solid (orebody, stope) | `vertex_attributes`, `face_attributes` |
| `dr.Polylines` | 3D strings: open lines and closed rings, grouped into features | `attributes` (per feature), `vertex_attributes` |
| `dr.PointSet` | points | `attributes` |
| `dr.Table` | typed columns with validity masks and units | |

`dr.UnstructuredGrid` holds the volume cells of a `.vtu` or `.vtk` file so they survive a round trip; it is data,
and `grid.surface()` turns its surface cells into a `Mesh`.

## Shared accessors

Every container has `.coords` as an `(n, 3)` float64 array, `.x`, `.y`, `.z`, `.bounds` as `(min, max)` or
`None`, `.column_names`, `.crs`, `.length_unit`, `to_table()`, and `to_pandas()`, `to_polars()` and
`to_pyarrow()`, which import their library when you call them.

Objects never change. `with_vertex_column`, `with_face_column`, `with_vertices`, `filter_faces`, `translate`
and the rest return new objects, and the arrays you get back are copies.

```python
import drusa as dr
import numpy as np

mesh = dr.Mesh(vertices, triangles, face_attributes={"rock": ["ox"] * len(triangles)}, crs="EPSG:31983", length_unit="m")
mesh = mesh.with_vertex_column("grade", np.random.default_rng(0).random(mesh.n_vertices), unit="g/t")
mesh["rock"]                  # a vertex column, else a face column
mesh.to_table(on="faces")     # v0, v1, v2 and the face columns
```

A `Mesh` winds a solid's triangles counter-clockwise seen from outside, so its signed volume is positive.

## Polylines

A part is one line or ring; a ring does not repeat its first vertex. A feature groups parts, so a polygon with a
hole is one feature with two rings. Inside is even-odd over a feature's rings, and a point on an edge counts as
inside.

```python
pad = dr.Polylines([outer, hole], closed=True, features=[0, 0], attributes={"name": ["pad 3"]})
pad.area(), pad.length()
pad.contains(points)                    # exact, in plan
pad.distance(points, signed=True)       # negative inside
dr.Polylines.from_table(table, feature="ID", x="X", y="Y", z="Z")
```

## Tables

```python
t = dr.Table({"hole": ["DH1", None], "depth": [120.0, np.nan], "n": [3, None]})
t["depth"]          # float64 with NaN for null
t["n"], t.mask("n") # int64 payload and validity
t.with_column("au", [1.2, 0.4], unit="g/t").units
t.drop_null()
```

Floats map to float64, integers to int64, booleans to bool and strings to utf8; `None` (and NaN in a float
column) marks a null.

## Labels

`crs` and `length_unit` are labels. Operations that combine two geometries check that declared labels match and
raise `InvalidInput` when they differ; an undeclared label matches anything. `mesh.to_length_unit("ft")`
rescales coordinates, and CRS transforms stay outside `drusa`.

## Errors

| Class | Also a | Raised for |
|---|---|---|
| `dr.DrusaError` | `Exception` | base of every drusa error |
| `dr.InvalidInput` | `ValueError` | bad shapes or values, unsupported content, CRS or unit mismatch |
| `dr.InvalidMesh` | `InvalidInput` | a mesh that fails validation; `.report` holds the `MeshReport` |
| `dr.FileError` | `OSError` | unreadable files and malformed content |
| `dr.MissingColumn` | `KeyError` | an unknown column name; the message lists the columns |

## Example

![Topography, orebody outline and drillholes](../examples/00-basics/basics.png)

??? example "examples/00-basics/example_00.py"

    ```python
    --8<-- "examples/00-basics/example_00.py"
    ```
