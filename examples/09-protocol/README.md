# The mining protocol

`Mesh`, `Polylines`, `PointSet` and `Table` expose `__minetypes__`. The dict is
protocol version 1: plain values and NumPy arrays, with no validity mask.
`from_minetypes` reads that dict, or any object that returns one.
`block_fraction` reads a `grid` or `block_model`. `trace_intervals` reads `drillholes`.

## A mesh crosses as views

`coords` and `triangles` share memory with the mesh. Triangles stay uint32.
A masked column is copied, and a null float becomes NaN.

```python
import numpy as np
import wireframes as wf

box = wf.datasets.box([0, 0, 0], [1, 1, 1])
spec = box.__minetypes__()
print(spec["kind"], spec["triangles"].dtype.name, spec["triangles"].flags.writeable)
assert spec["protocol"] == 1
assert np.shares_memory(spec["coords"], box.__minetypes__()["coords"])
assert wf.Mesh.from_minetypes(spec).volume == 1.0

table = wf.Table({"au": [1.2, None]})
print(table.__minetypes__()["columns"]["au"].tolist())
```

```text
mesh uint32 False
[1.2, nan]
```

## A grid is the only argument

`origin` is the minimum corner of cell `(0, 0, 0)`. The fraction of the one
cell that fills the box is 1. An open part of one vertex is a `PointSet`.

```python
solid = wf.datasets.box([0, 0, 0], [2, 2, 2])
model = {
    "protocol": 1,
    "kind": "grid",
    "columns": {},
    "units": {},
    "categories": {},
    "crs": None,
    "length_unit": None,
    "origin": [0.0, 0.0, 0.0],
    "size": [2.0, 2.0, 2.0],
    "count": [1, 1, 1],
    "rotation": [0.0, 0.0, 0.0],
}
print(solid.block_fraction(model).tolist())
assert solid.block_fraction(model).tolist() == [1.0]

point = wf.Polylines([[[0.0, 0.0, 0.0]]])
try:
    point.__minetypes__()
except wf.InvalidInput as error:
    print(error)
```

```text
[1.0]
part 0 has 1 vertex; open parts need at least 2. A single point belongs in a PointSet
```

Full script: [`example_09.py`](example_09.py), run by `mise run examples`.
