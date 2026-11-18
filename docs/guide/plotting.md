# Plotting

`dr.plot` draws static matplotlib figures. It needs `pip install "drusa[plot]"`; `import drusa` alone never
imports matplotlib. Each function takes one container or a list of them, an optional `ax`, and returns
`(fig, ax)`; with `ax` given, it draws there.

| Function | Draws |
|---|---|
| `dr.plot.plan(data, values=None, edges=False, label=False)` | meshes coloured by a column or by z, polylines as lines, points as dots |
| `dr.plot.section(data, origin, normal, thickness=None, values=None)` | mesh section lines, plus polylines and points within `thickness / 2` of the plane |
| `dr.plot.surface3d(data, values=None, problems=None, azim=-60, elev=30)` | shaded triangles, lines and points in 3D |

`values` names a face or vertex column, or gives one value per face or vertex. In `section`, `values` names a
column of the section lines that picks each line's colour. A vertical section plots the distance along the
section against elevation; any other plane uses an in-plane basis.

`surface3d(mesh, problems=mesh.validate())` draws boundary and non-manifold edges in red and other problem
faces in orange, which makes holes and folds easy to find.

```python
fig, ax = dr.plot.plan([topography, outline, holes], edges=True)
fig, ax = dr.plot.section([orebody, holes], origin, [1, 0, 0], thickness=40)
fig, ax = dr.plot.surface3d(holed, problems=holed.validate())
```

![Plan view](../examples/00-basics/basics.png)

![3D view with a hole marked](../examples/02-measure-repair/measure_repair.png)
