# Measure, validate and repair

## Measures

| Property or method | Returns |
|---|---|
| `area`, `face_areas` | total and per-triangle area |
| `volume`, `centroid` | signed volume and volume centroid; `InvalidMesh` unless closed and consistently wound |
| `face_normals`, `face_centroids`, `vertex_normals` | `(m, 3)` or `(n, 3)` arrays |
| `dip` | `(m, 2)` dip and dip direction in degrees, clockwise from north |
| `volume_by_elevation(edges)` | solid volume between consecutive elevations |
| `surface_z(points)`, `vertical_distance(points)` | surface elevation above points, and point height above it |
| `thickness(bottom, points, mode=...)` | vertical or normal thickness between two surfaces |

Measures subtract the bounding-box corner before summing, so UTM coordinates keep their digits.
`volume_by_elevation` integrates clamped columns exactly, and its slabs add back to `volume`.

## Validation

`mesh.validate(tolerance=0.0)` returns a `MeshReport`. `report.summary` counts degenerate and duplicate faces,
duplicate and unused vertices, boundary and non-manifold edges, non-manifold vertices, inconsistent windings,
shells, inward shells, self-intersections and coplanar overlaps, and adds `is_closed` and `is_solid`.
`report.problems` is a table with one row per problem and the face, vertex and edge involved.

A solid is closed, manifold, consistently and outward wound, and free of self-intersections. Its closed shells
may touch each other along edges or at vertices: the summary counts those contacts in `touching_shell_edges`
and `touching_shell_vertices`, and they leave `is_solid` true. Every check uses exact predicates; `tolerance`
only widens the duplicate-vertex test within one shell.

## Repair

```python
fixed, report = mesh.repair(weld=0.0, fill_holes="min_area")
report.actions  # one row per change: action, vertex, face, other, reason
report.vertex_map  # output vertex of each input vertex, -1 when removed
report.face_map  # output face of each input face, -1 when removed
report.holes_left  # boundary loops still open after the run
report.summary  # count per action
```

Steps run in order: weld, remove degenerate faces, remove duplicate faces, remove unused vertices, orient,
fill holes. Welding keeps closed shells that touch apart: each keeps its own copy of a shared vertex, and the
report lists a `vertex_split` row where the input shared one between shells. Pass `weld_touching_shells=True`
to merge those vertices as well. Orienting flips the fewest faces to agree with each shell's first face, then flips whole shells
outward by nesting parity. `repair` never touches self-intersections, and running it on its own output changes
nothing.

## Components and transforms

`components()` splits a mesh into edge-connected pieces and `component_labels()` labels each triangle. Shells
that touch along an edge or at a vertex come out as separate pieces.
`transform(matrix)` applies a `(4, 4)` affine matrix and flips the winding when the determinant is negative;
`translate`, `flip` and `to_length_unit` cover the common cases.

## Example

![A hole drawn from the validation report, and volume by slab](../examples/02-measure-repair/measure_repair.png)

??? example "examples/02-measure-repair/example_02.py"

    ```python
    --8 < --"examples/02-measure-repair/example_02.py"
    ```
