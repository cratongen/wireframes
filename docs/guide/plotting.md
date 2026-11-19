# Plotting

`wf.plot` draws static matplotlib figures. It needs matplotlib (the `wireframes[plot]` extra); `import wireframes` alone never
imports matplotlib. Each function takes one container or a list of them, an optional `ax`, and returns
`(fig, ax)`; with `ax` given, it draws there.

| Function | Draws |
|---|---|
| `wf.plot.plan(data, values=None, edges=False, label=False)` | meshes coloured by a column or by z, polylines as lines, points as dots |
| `wf.plot.section(data, origin, normal, thickness=None, values=None)` | mesh section lines, plus polylines and points within `thickness / 2` of the plane |
| `wf.plot.surface3d(data, values=None, problems=None, azim=-60, elev=30)` | shaded triangles, lines and points in 3D |

`values` names a face or vertex column, or gives one value per face or vertex. In `section`, `values` names a
column of the section lines that picks each line's colour. A vertical section plots the distance along the
section against elevation; any other plane uses an in-plane basis.

`surface3d(mesh, problems=mesh.validate())` draws boundary and non-manifold edges in red and other problem
faces in orange, which makes holes and folds easy to find.

```python
import wireframes as wf

topography = wf.datasets.topography()
orebody = wf.datasets.orebody()
holes = wf.datasets.drillholes()
origin = orebody.centroid
outline = orebody.slice_by_plane(origin, [0, 0, 1])
holed = orebody.filter_faces(orebody.face_centroids[:, 2] < origin[2] + 40)

fig, ax = wf.plot.plan([topography, outline, holes], edges=True)
fig, ax = wf.plot.section([orebody, holes], origin, [1, 0, 0], thickness=40)
fig, ax = wf.plot.surface3d(holed, problems=holed.validate())
```

![Plan view](../examples/00-basics/basics.png)

![3D view with a hole marked](../examples/02-measure-repair/measure_repair.png)
