# Booleans and cutting

| Operator | Function | Keeps |
|---|---|---|
| `a \| b` | `dr.union(a, b)` | inside `a` or `b` |
| `a & b` | `dr.intersection(a, b)` | inside both |
| `a - b` | `dr.difference(a, b)` | inside `a`, outside `b` |
| `a ^ b` | `dr.symmetric_difference(a, b)` | inside exactly one |

The methods `a.union(b)` and so on do the same. Each input must be a solid (`is_solid`); otherwise the call
raises `dr.InvalidMesh` with the report. Pass `repair=True` to repair both inputs first; the call then returns
`(mesh, (report_a, report_b))`.

The arrangement splits every pair of crossing triangles along their exact intersection, using implicit
intersection points and exact orientation tests. Coordinates become f64 only when the result is written out.
Face columns follow the triangles they came from, so a `source` column tells you which side built each face.

The volumes obey inclusion and exclusion to rounding:

- `vol(a | b) + vol(a & b) == vol(a) + vol(b)`
- `vol(a - b) == vol(a) - vol(a & b)`
- `vol(a ^ b) == vol(a | b) - vol(a & b)`

!!! note "Symmetric difference"
    `a ^ b` holds `a - b` and `b - a` in one mesh. The two pieces meet along the intersection curve, so
    `validate` reports duplicate vertices there and the result is not a solid you can feed to another boolean
    or query. Compute `a - b` and `b - a` separately when you need solids.

Solids that touch only along an edge or at a vertex give a non-manifold result, which raises `InvalidInput`.

## Cutting

```python
below, above = solid.clip_by_surface(fault)                 # both capped; above is where the normals point
negative, positive = mesh.split_by_plane(origin, normal)    # exact side test; on-plane points go positive
lines = mesh.slice_by_plane(origin, normal, offsets=[0, 10, 20])
```

`clip_by_surface` needs the surface boundary outside the solid. `slice_by_plane` returns one feature per line
with an `offset` column and the face columns of the triangles it crosses; a line splits where those columns
change.

## Cut and fill

`dr.cut_fill(base, design, boundary=None)` overlays two height-field surfaces in plan and integrates
`design - base` exactly over each overlay cell. It returns `cut`, `fill`, `net`, `cut_area`, `fill_area` and
`area`. Swapping the surfaces swaps cut and fill.

## Example

![Sections of the orebody against a mining block](../examples/03-booleans/booleans.png)

??? example "examples/03-booleans/example_03.py"

    ```python
    --8<-- "examples/03-booleans/example_03.py"
    ```
