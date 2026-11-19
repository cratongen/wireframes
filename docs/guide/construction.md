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
| `mesh.smooth(iterations=10, lambda_=0.5, mu=-0.53, preserve_volume=False)` | Taubin smoothing, boundary fixed; `preserve_volume=True` keeps a solid's volume |
| `mesh.remesh(edge_length, iterations=5)` | isotropic remeshing |
| `mesh.offset(distance)` | vertices moved along their normals, outward when positive |

`dtm` keeps breaklines and boundary rings as edges, drops points closer than `min_spacing` in plan to an earlier
one, and winds faces upward. Where two breaklines cross, the new vertex takes its elevation from the one inserted
first, in feature order, so put contours with exact elevations before strings that interpolate. With the default
`min_spacing=0`, two points at one plan position with different elevations raise `wf.InvalidInput` naming both.
It keeps the `crs` and `length_unit` of a `PointSet` or `Mesh` you pass it.

`loft` orients the rings, aligns their start vertices, and picks each band's triangles by dynamic programming to
minimise their area with both rings centred and seen along the band axis, so long bands and sections cut obliquely
across the axis link vertex to vertex instead of in fans. It takes one closed ring per section, except that the first and last sections may be an open line or a single point
where a vein pinches out. Branching is not supported yet: a section with two or more rings raises
`wf.InvalidInput` naming the section, so loft each branch on its own and unite the solids. The loft keeps every
string vertex as a vertex of the solid, so the solid passes through the interpretation on each section.

`extrude` sweeps each closed ring along a straight vector. The volume of the prism is the ring's area in its own
plane times the component of the vector along the ring normal; note that `Polylines.area()` is the plan area,
which is zero for a ring in a vertical plane such as a drive profile. `extrude` has no path sweep: build a drive
along a curved centreline by lofting rings placed along it. `offset` refuses a result that intersects itself.
It moves vertices by `distance`, not faces: at a box corner each face moves by `distance / sqrt(3)`, so it is not
a skin of constant thickness. `loft` links one ring per section; a section with two or more rings, as where a
lens splits, raises `InvalidInput`, so loft each branch on its own and unite the pieces.

Every solid builder validates its output. The checks the test suite runs on them:

- the volume of an extrusion equals ring area times height, and a drive profile swept along strike gives profile
  area times length ([example 07](../examples/07-underground/README.md));
- a loft of two squares gives the box or frustum volume;
- the volume of `close_to_elevation(z)` equals the fill of `cut_fill` against a flat plane at `z`;
- a pit design closed above the ground and clipped by the topography has the volume of the design's cut against
  the topography ([example 06](../examples/06-open-pit/README.md));
- contours of a plane are straight lines on that plane.

## Example

![A DTM with a breakline, and edited meshes in section](../examples/05-construct/construct.png)

The full example, with the output of each step: [example 05](../examples/05-construct/README.md).
