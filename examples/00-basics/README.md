# Meshes, strings, points and tables

`wireframes` has three geometry containers. A `Mesh` holds triangles: an open surface such as topography, or a closed
solid such as an orebody. `Polylines` holds 3D strings grouped into features, and `PointSet` holds points. Each
carries attribute tables, a `crs` label and a `length_unit` label, and never changes once built: every method
returns a new object.

```python
import importlib.util

import numpy as np
import wireframes as wf
from common import save
```

## A mesh from arrays

A mesh takes `(n, 3)` vertices and `(m, 3)` triangles. A solid winds its triangles counter-clockwise seen from
outside, so its signed volume is positive. Here a 20 m cube gets a face column naming its sides.

```python
corners = np.array([[x, y, z] for z in (0, 20) for y in (0, 20) for x in (0, 20)], float)
faces = wf.datasets.box([0, 0, 0], [1, 1, 1]).triangles
sides = ["bottom", "bottom", "top", "top", "south", "south", "north", "north", "west", "west", "east", "east"]
cube = wf.Mesh(corners, faces, face_attributes={"side": sides}, crs="EPSG:31983", length_unit="m")
print(cube)
print(f"volume {cube.volume:.1f} m3, area {cube.area:.1f} m2, solid: {cube.is_solid}")
print(cube.to_table(on="faces").to_dict()["side"][:4])
```

```text
Mesh(8 vertices, 12 triangles, columns ["side"], crs EPSG:31983, length_unit m)
volume 8000.0 m3, area 2400.0 m2, solid: True
['bottom' 'bottom' 'top' 'top']
```

## Tables keep types and nulls

Attribute tables hold float64, int64, bool and string columns. `None` marks a null; numeric nulls come back as
NaN, and `mask` tells which rows hold a value.

```python
table = wf.Table({"hole": ["DH1", "DH2", None], "depth": [120.0, None, 80.5], "assayed": [True, False, True]})
print(table)
print(table["depth"], table.mask("hole"))
print(table.drop_null().num_rows, "rows without nulls")
```

```text
Table(3 rows; hole: utf8, depth: float64, assayed: bool)
[120.    nan  80.5] [ True  True False]
1 rows without nulls
```

## The synthetic deposit

`wf.datasets` builds a seeded synthetic deposit at UTM scale: a lens-shaped orebody, a rolling topography and a
fan of drillholes. The plan view draws the topography by elevation, the orebody outline at 200 m and the holes.

```python
orebody = wf.datasets.orebody()
topography = wf.datasets.topography()
holes = wf.datasets.drillholes()
collars = wf.PointSet(
    np.array([p[0] for p in holes.parts]), {"hole": holes.attributes["hole"]}, length_unit="m"
)
print(orebody, holes, collars, sep="\n")

outline = orebody.slice_by_plane(wf.datasets.ORIGIN + [0, 0, 200], [0, 0, 1])
fig, ax = wf.plot.plan([topography, outline, holes, collars], label=False)
ax.set_title("Topography, orebody at 200 m and drillholes")
save(fig, "basics")
```

```text
Mesh(216 vertices, 428 triangles, columns ["section"], crs None, length_unit m)
Polylines(12 features, 12 parts, 24 vertices, columns ["hole", "depth"])
PointSet(12 points, columns ["hole"])
```

![basics](basics.png)

Every container converts to a table, and from there to pandas, polars or pyarrow when they are installed.

```python
print(holes.to_table().column_names)
assert cube.volume == 8000.0
```

```text
['feature', 'part', 'vertex', 'x', 'y', 'z', 'closed', 'hole', 'depth']
```

## Joining meshes and editing tables

`wf.concat` appends meshes in order without welding, so volumes add. `with_attributes` swaps whole attribute
tables, and `Table.with_columns` adds several columns at once.

```python
second = cube.translate([40, 0, 0])
pair = wf.concat([cube, second])
print(pair.n_faces, "faces,", f"{pair.volume:.0f} m3")
assert pair.volume == cube.volume + second.volume
relabelled = pair.with_attributes(face={"block": [1] * 12 + [2] * 12})
print(relabelled.column_names, table.with_columns({"lab": ["A", "B", "C"]}).column_names)
```

```text
24 faces, 16000 m3
['block'] ['hole', 'depth', 'assayed', 'lab']
```

## Strings from a table

`Polylines.from_table` groups rows into features by an id column, as a CSV export of a pit crest would give.
`length` measures each feature in 3D, `contains` tests points in plan, and `to_points` turns every vertex into a
point with its feature, part and vertex index.

```python
crest = wf.Polylines.from_table(
    {"ID": ["pad"] * 4, "X": [0, 30, 30, 0], "Y": [0, 0, 30, 30], "Z": [5, 5, 5, 5]}, z="Z", closed=True
)
print(crest.length(), crest.contains([[15, 15], [45, 15]]))
assert crest.length()[0] == 120.0 and crest.contains([[15, 15]])[0]
vertices = crest.to_points()
print(vertices.column_names, len(vertices))
```

```text
[120.] [ True False]
['feature', 'part', 'ID'] 4
```

The DataFrame converters keep nulls and need only the library they name.

```python
for lib in ("polars", "pandas", "pyarrow"):
    if importlib.util.find_spec(lib):
        print(
            lib,
            len(getattr(table, f"to_{lib}")()),
            "table rows,",
            len(getattr(collars, f"to_{lib}")()),
            "collars",
        )
```

```text
polars 3 table rows, 12 collars
pandas 3 table rows, 12 collars
pyarrow 3 table rows, 12 collars
```

Full script: [`example_00.py`](example_00.py), run by `mise run examples`.
