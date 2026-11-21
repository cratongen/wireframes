# Booleans and cutting

| Operator | Function | Keeps |
|---|---|---|
| `a \| b` | `wf.union(a, b)` | inside `a` or `b` |
| `a & b` | `wf.intersection(a, b)` | inside both |
| `a - b` | `wf.difference(a, b)` | inside `a`, outside `b` |
| `a ^ b` | `wf.symmetric_difference(a, b)` | inside exactly one |

The methods `a.union(b)` and so on do the same. Each input must be a solid (`is_solid`); otherwise the call
raises `wf.InvalidMesh` with the report. Pass `repair=True` to repair both inputs first; the call then returns
`(mesh, (report_a, report_b))`.

The arrangement splits every pair of crossing triangles along their exact intersection, using implicit
intersection points and exact orientation tests. Coordinates become f64 only when the result is written out.
Face columns follow the triangles they came from, so a `source` column tells you which side built each face.

The volumes obey inclusion and exclusion to rounding:

- `vol(a | b) + vol(a & b) == vol(a) + vol(b)`
- `vol(a - b) == vol(a) - vol(a & b)`
- `vol(a ^ b) == vol(a | b) - vol(a & b) == vol(a) + vol(b) - 2 vol(a & b)`

## Touching shells

A solid may hold several closed shells that touch each other along edges or at vertices, as long as their
interiors do not overlap. Two stopes that share a corner, or `a - b` and `b - a` meeting along their
intersection curve, are one valid solid. Every boolean accepts such solids and returns them:

```python
import wireframes as wf

a = wf.datasets.box([0, 0, 0], [1, 1, 1])
b = wf.datasets.box([1, 1, 0], [2, 2, 1])  # shares the edge x = y = 1 with a
two = a | b
two.validate().summary["touching_shell_edges"]  # 1: information, not a defect
two.is_solid  # True
(two | wf.datasets.box([0.5, 0.5, 0], [1.5, 1.5, 1])).volume  # 2.5: results chain
```

Each shell of a result keeps its own vertices where it touches another, and `validate` reports the contacts in
`touching_shell_edges` and `touching_shell_vertices`. Inputs may share those vertices by index instead; both
forms validate. `a ^ b` holds `a - b` and `b - a` as separate shells along the intersection curve.

Validation still refuses contacts that are not between separate closed shells: duplicate vertices inside one
shell, a vertex of one shell resting inside another shell's face, coplanar overlapping faces (two solids that
share a face; their union merges them instead), and contacts where the interiors overlap. One shell may touch
itself on the same terms, such as an L-shaped body whose two arms meet along an edge, or a ring of stopes that
closes around a gap: the result is a solid with one shell and the contact in `touching_shell_edges`, so a
stope-by-stope union does not depend on the order of the stopes.

## Rounding

Exact predicates decide every crossing, and the result is exact until its new vertices are written as f64. Where
two inputs agree only to the last digit, the exact result can hold slivers thinner than f64 can place, and
rounding them alone would fold them into their neighbours. Booleans, `clip_by_surface` and `split_by_plane`
therefore round with exact tests: vertices that round together become one, a vertex within rounding of a face it
does not touch is inserted into that face, and pieces whose rounded volume is zero go. Where two inputs agree to
rounding with different triangulations, such as two pits clipped by the same topography, each input's vertices
are inserted into the other's faces there and edges are flipped until both match, before the exact arrangement
runs. No vertex ends more than four ulps of its largest coordinate from its exact position, and the volume of the
result differs from the exact result's by at most a bound computed with it. Every result is then checked with
`validate`.

Open-pit work meets this when it combines pieces that were each clipped by the same topography: their caps come
from the same topography triangles cut at different places.

```python
topography = wf.datasets.topography()
sky = topography.translate([0, 0, 200])  # pit walls daylight above the ground
closed1 = wf.datasets.pit(sky, radius=60).close_to_elevation(600)
closed2 = wf.datasets.pit(sky, radius=90).close_to_elevation(600)
pit1 = closed1.clip_by_surface(topography)[0]
pit2 = closed2.clip_by_surface(topography)[0]
pushback, rounding = wf.difference(pit2, pit1, rounding_report=True)
pushback.is_solid  # True
rounding.moved_input_vertices, rounding.volume_bound  # 0 and about 7e-07 m3
```

`rounding_report=True` on any boolean, `clip_by_surface` or `split_by_plane` returns a `wf.RoundingReport`
last: how many vertices moved and how far at most (`max_move`, `move_bound`), what each operation did, and
`volume_bound`, a proven bound on the difference from the exact result's volume.

Clipping the difference of the closed, unclipped solids once gives the pushback too:

```python
pushback = (closed2 - closed1).clip_by_surface(topography)[0]
pushback.is_solid  # True
```

A result the rounding cannot certify raises `wf.KernelError` whose message starts with `uncertified`;
`error.report` holds the `MeshReport` of the rounded result and `error.mesh` the result itself, closed, so you can
still read its volume. Ignored tests in `crates/boolean/tests/attacks.rs` still raise it: stacked stopes at dip 89.6 minus a
drive on their floor, splits of rotated stope blocks through one of their own vertices, a split through a slot a few ulps
wide, repeated union then difference of a raise in a worked orebody, and the running union of ten pushbacks. Two of those
splits panic instead of raising. To cut a stope block at a level that passes through its vertices, intersect it with a box
whose top lies on the level; the tests certified that in every case they tried. The method has no general validity
theorem, so inputs the tests have not met can raise it too. An exact coordinate beyond the f64 range raises it
with `non-finite`.

Underground work meets it in long chains of booleans. Each result carries rounded vertices along its intersection
lines, and a boolean that takes it as input treats them as exact; long, thin triangles near them would fold over
each other once rounded again. Two habits keep the number of moved vertices down. Put pieces that do not touch into one mesh with
`wf.concat` and run one boolean, instead of one boolean per piece. And build each reported quantity from the
original solids rather than from an earlier result:

```python
drive = wf.datasets.box([0, 0, 0], [100, 5, 4])
crosscuts = [wf.datasets.box([x, 2, 0.5], [x + 4, 30, 3.5]) for x in range(10, 100, 20)]
stopes = wf.concat([wf.datasets.box([x - 3, 20, 0], [x + 7, 40, 20]) for x in range(10, 100, 20)])
orebody = wf.datasets.box([0, 25, 0], [100, 45, 20])
development = drive | wf.concat(crosscuts)

level = drive | wf.concat(crosscuts)  # one union; the crosscuts do not touch each other
dilution = (stopes - orebody).volume - ((stopes & development) - orebody).volume  # from original solids
```

[Example 06](../examples/06-open-pit/README.md) shows both orders of clipping, and
[example 07](../examples/07-underground/README.md) shows a chain of twelve unions that used to fail at its
eleventh step.

## Cutting

```python
solid = wf.datasets.orebody()
fault = wf.datasets.topography().translate([0, 0, solid.centroid[2] - 300])  # spans the orebody in plan
origin, normal = solid.centroid, [1, 0, 0]

below, above = solid.clip_by_surface(fault)  # both capped; above is where the normals point
negative, positive = solid.split_by_plane(origin, normal)  # exact side test; on-plane points go positive
lines = solid.slice_by_plane(origin, normal, offsets=[0, 10, 20])
```

`clip_by_surface` needs the surface boundary outside the solid: a fault or topography that ends inside the solid
raises `wf.InvalidInput`, so extend the surface past the solid first. Both inputs must pass validation. Cap faces
take the surface's face columns, so a `surface` column on the topography marks which faces of a clipped pit came
from the ground and which from the design. `slice_by_plane` returns one feature per line with an `offset` column
and the face columns of the triangles it crosses; a line splits where those columns change.

## Cut and fill

`wf.cut_fill(base, design, boundary=None)` overlays two height-field surfaces in plan and integrates
`design - base` exactly over each overlay cell. It returns `cut`, `fill`, `net`, `cut_area`, `fill_area` and
`area`, in the data's units. Swapping the surfaces swaps cut and fill.

The domain is where both surfaces exist in plan, inside the closed rings of `boundary` when given. A surface with
two faces over one plan point, such as an overhanging wall, raises `wf.InvalidInput`. Over one domain, nets add
up: `net(a, c) == net(a, b) + net(b, c)` to rounding, which checks a survey against a design and the topography
in one line.

## Example

![Sections of the orebody against a mining block](../examples/03-booleans/booleans.png)

The full example, with the output of each step: [example 03](../examples/03-booleans/README.md).
