# Construction

| Call | Builds |
|---|---|
| `wf.dtm(points, breaklines=None, boundary=None, min_spacing=0.0, max_edge=None)` | surface by constrained Delaunay in plan |
| `wf.triangulate(polygons)` | each feature's rings triangulated in their best-fit plane, holes by even-odd |
| `wf.extrude(polygons, vector)` | solid swept by closed rings |
| `wf.loft(sections, cap_ends=True)` | solid linking one ring per feature to the next |
| `surface.close_to_elevation(z)` | solid between a height-field surface and a flat floor |
| `mesh.contours(levels)` | contour lines with an `elevation` column |
| `mesh.decimate(target_faces=None, max_error=None)` | fewer triangles by quadric edge collapse |
| `mesh.smooth(iterations=10, lambda_=0.5, mu=-0.53)` | Taubin smoothing, boundary fixed |
| `mesh.remesh(edge_length, iterations=5)` | isotropic remeshing |
| `mesh.offset(distance)` | vertices moved along their normals, outward when positive |

`dtm` keeps breaklines and boundary rings as edges, drops points closer than `min_spacing` in plan to an earlier
one, and winds faces upward. It keeps the `crs` and `length_unit` of a `PointSet` or `Mesh` you pass it.

`loft` orients the rings, aligns their start vertices, and picks each band's triangles by dynamic programming to
minimise their area with both rings centred and seen along the band axis, so long bands and sections cut obliquely
across the axis link vertex to vertex instead of in fans. It takes one closed ring per section, except that the first and last sections may be an open line or a single point
where a vein pinches out; branching is not supported yet. `offset` refuses a
result that intersects itself.

Every solid builder validates its output. The checks the test suite runs on them:

- the volume of an extrusion equals ring area times height;
- a loft of two squares gives the box or frustum volume;
- the volume of `close_to_elevation(z)` equals the fill of `cut_fill` against a flat plane at `z`;
- contours of a plane are straight lines on that plane.

## Example

![A DTM with a breakline, and edited meshes in section](../examples/05-construct/construct.png)

??? example "examples/05-construct/example_05.py"

    ```python
    --8 < --"examples/05-construct/example_05.py"
    ```
