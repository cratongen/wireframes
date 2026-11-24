# Queries

| Method | Returns |
|---|---|
| `locate(points)` | int8: -1 inside, 0 on the surface, 1 outside |
| `contains(points, boundary=True)` | bool; `boundary=False` treats the solid as open |
| `signed_distance(points)` | distance to the surface, negative inside |
| `distance(points)`, `closest_point(points)` | unsigned distance; closest points, faces and distances |
| `ray_hits(origins, directions)`, `segment_hits(starts, ends)` | table of hits: `t`, `face`, `x`, `y`, `z`, `kind`, `entering` |
| `trace_intervals(trace, depths=None)` | table of runs: `from`, `to`, `inside`, `entry_face`, `exit_face` |
| `block_fraction(centroids, size=None, rotation=None)` | fraction of each block inside the solid |

`locate`, `contains`, `signed_distance`, `trace_intervals` and `block_fraction` need a solid; `distance`,
`closest_point` and the hit queries take any mesh. The mesh builds its spatial index on the first query and
keeps it, and every call releases the GIL.

`locate` casts a vertical ray and counts crossings with exact predicates; a ray through an edge or a vertex
counts once, by symbolic perturbation. A vertex of the mesh is always on its surface.

`trace_intervals` takes a desurveyed trace, `(n, 3)` vertices, and the measured depth of each vertex,
non-decreasing. The runs cover `depths[0]` to `depths[-1]` in order and alternate between outside and inside;
the depth of each crossing comes from an exact segment-face test and is interpolated linearly along its trace
segment. A run that lies in a face counts as inside, and zero-length runs are dropped. `entry_face` and
`exit_face` are null at the two ends of the trace, where a run starts or stops without crossing a face; null
reads as 0 in the int64 array, so check `iv.mask("entry_face")`. Depths out of order, or a depth count that
differs from the vertex count, raise `wf.InvalidInput`.

`block_fraction` integrates the solid over each block with the divergence theorem. There is no sampling error:
over blocks that cover the solid, fraction times block volume sums to the solid volume, and the fractions of
solids that split another solid (fault blocks of an orebody, or its mined and unmined parts) add up block by block
to the fraction of the whole. `rotation` takes a `(3, 3)` matrix whose columns are the block axes. A size that is
not positive raises `wf.InvalidInput`.

A protocol object of kind `grid` or `block_model` is the only argument. `origin` is the minimum corner of cell
`(0, 0, 0)`, and `rotation` is `(azimuth, dip, rake)` in degrees. A sub-blocked model passes `parent`, `offset`
and `scale`. A `drillholes` object is the only argument of `trace_intervals`: each hole is the slice of
`path_coords` and `path_depths` between `path_offsets`, and the table adds `hole_index` as int64.

Both need a solid: an open or otherwise invalid mesh raises `wf.InvalidMesh` with the report.

```python
import numpy as np
import wireframes as wf

solid = wf.datasets.orebody()
lo, hi = solid.bounds
axes = [np.arange(a + s / 2, b + s, s) for a, b, s in zip(lo, hi, [10, 10, 5])]
centroids = np.column_stack([g.ravel() for g in np.meshgrid(*axes, indexing="ij")])
trace_xyz = np.array([solid.centroid + [0, 0, 150], solid.centroid - [0, 0, 150]])
depths = np.array([0.0, 300.0])

f = solid.block_fraction(centroids, [10, 10, 5])
iv = solid.trace_intervals(trace_xyz, depths)
inside = iv.filter(iv["inside"])

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
cube = wf.datasets.box([0, 0, 0], [2, 2, 2])
assert cube.block_fraction(model).tolist() == [1.0]
```

## Example

![Signed distance and block fractions](../examples/04-queries/queries.png)

The full example, with the output of each step: [example 04](../examples/04-queries/README.md).
[Example 07](../examples/07-underground/README.md) intersects desurveyed drillholes with an orebody and its fault
blocks, checks that each entry and exit point lies on the surface, and builds a tonnage report from block
fractions with a block-by-block balance.
