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
`volume_by_elevation` integrates clamped columns exactly, and its slabs add back to `volume`. A slab the solid
does not reach can come out as a rounding residue of either sign, such as `-4.5e-13`, instead of 0.

`volume` refuses a mesh that is not closed or not consistently wound, but it does not look for
self-intersections: two shells that overlap give a number that counts the overlap twice. Check `is_solid`
before you report a volume.

## Validation

`mesh.validate(tolerance=0.0)` returns a `MeshReport`. `report.summary` counts degenerate and duplicate faces,
duplicate and unused vertices, boundary and non-manifold edges, non-manifold vertices, inconsistent windings,
shells, inward shells, self-intersections and coplanar overlaps, and adds `is_closed` and `is_solid`.
`report.problems` is a table with one row per problem: `kind`, `face`, `vertex`, `edge` (the local edge 0, 1 or
2 of `face`) and `other` (the earlier duplicate, or the other face of an edge or a pair), with nulls where a
column does not apply.

Every check uses exact predicates. `tolerance` only widens the duplicate-vertex test: `validate(tolerance=0.002)`
reports vertices within 2 mm of each other, which finds corners that another package wrote a little apart.

### What a valid solid is

`is_solid` is true when all of these hold:

1. Closed: every edge joins two triangles, or an even number that pair into wedges of solid around it.
2. Manifold: around each vertex the triangles of a shell form one fan, unless every fan there belongs to a
   closed shell.
3. Consistent winding: two triangles that share an edge run it in opposite directions.
4. Outward: each closed shell has positive signed volume, and a shell inside another one (a cavity) has negative
   signed volume, by exact nesting.
5. No self-intersections or coplanar overlaps: no two triangles cross, touch away from shared vertices, or
   overlap in one plane.
6. No degenerate or duplicate faces, and no duplicate vertices except exact copies where closed shells touch.

Closed shells may touch each other along edges or at vertices, as primary stopes in a checkerboard do, and one
shell may touch itself the same way, as a ring of stopes around a pillar does. Cavities follow the same rule: a
cavity may touch its outer shell, and two cavities may touch each other, along edges or at vertices. The summary counts those
contacts in `touching_shell_edges` and `touching_shell_vertices`, and they leave `is_solid` true. Contacts are
exact: the edges and vertices must have equal coordinates. Only closed shells touch; where an open shell meets
another one the report lists non-manifold vertices or edges instead. Shells that share part or all of a face,
or meet on part of an edge, overlap: `validate` reports coplanar overlaps or self-intersections, and `union`
merges such pieces into one shell.

## Repair

```python
import wireframes as wf

box = wf.datasets.box([0, 0, 0], [10, 10, 10])
mesh = box.filter_faces([True] * 11 + [False])  # one triangle missing

fixed, report = mesh.repair(weld=0.0, fill_holes="min_area")
report.actions  # one row per change: action, vertex, face, other, reason
report.vertex_map  # output vertex of each input vertex, -1 when removed
report.face_map  # output face of each input face, -1 when removed
report.before, report.after  # validation summaries of input and output
report.holes_left  # boundary loops still open after the run
report.summary  # count per action
```

Steps run in order: weld, remove degenerate faces, remove duplicate faces, remove unused vertices, orient,
fill holes. `repair` returns a copy and leaves the input as it was.

- Weld merges vertices within `weld` (0 for equal coordinates, None to skip), keeping the earliest vertex of
  each group and its attributes (`vertex_merged`). It keeps closed shells that touch apart: each keeps its own
  copy of a shared vertex, and the report lists a `vertex_split` row where the input shared one between shells,
  as a mesh read from DXF or STL does. Pass `weld_touching_shells=True` to merge those vertices instead.
- Degenerate and duplicate faces go (`face_removed`, with the problem kind in `reason`), keeping the earliest
  of each duplicate.
- Orienting flips the fewest faces to agree with each shell's first face, then flips whole shells outward by
  nesting parity (`face_flipped`).
- Holes stay open unless `fill_holes` is `"fan"` or `"min_area"` (`face_added`). Added faces go last with null
  face columns: `repair` never invents attribute values. `holes_left` counts loops left open.

`repair` never moves a vertex except by welding, never touches self-intersections or coplanar overlaps, and
running it on its own output changes nothing. Its counts balance: input vertices minus `vertex_merged` and
`vertex_removed`, plus `vertex_split`, are the output vertices, and input faces minus `face_removed` plus
`face_added` are the output faces. A script should check `report.after["is_solid"]` and `report.holes_left`
before using the result.

## Components and transforms

`components()` splits a mesh into edge-connected pieces and `component_labels()` labels each triangle. Shells
that touch along an edge or at a vertex come out as separate pieces.
`transform(matrix)` applies a `(4, 4)` affine matrix and flips the winding when the determinant is negative;
`translate`, `flip` and `to_length_unit` cover the common cases.

## Examples

![A hole drawn from the validation report, and volume by slab](../examples/02-measure-repair/measure_repair.png)

The full examples, with the output of each step: [example 02](../examples/02-measure-repair/README.md) measures,
validates and repairs a solid, and [example 08](../examples/08-files-repair/README.md) repairs a Datamine pair
and a DXF from other packages, reads each repair report and sets out the rules above case by case.
