# Files

`dr.read(path)` returns a `Dataset` with `mesh`, `polylines`, `points`, `grid`, the `length_unit` the file
declares and the reader's `warnings`. `dr.read_mesh`, `dr.read_polylines` and `dr.read_points` return one part
and raise `InvalidInput` when the file has none. `dr.write(path, data)` takes a container, a grid or a `Dataset`;
the extension picks the format.

| Format | Extensions | Holds | Columns it stores |
|---|---|---|---|
| Wavefront OBJ | `.obj` | faces, lines, points | face `object`, `group` |
| Polygon File Format | `.ply` | faces, ASCII or binary | numeric vertex and face columns |
| Object File Format | `.off` | faces | none |
| Stereolithography | `.stl` | triangles, ASCII or binary, 32-bit floats | face `solid`, `stl_attribute` |
| AutoCAD DXF | `.dxf` | 3DFACE, POLYFACE, MESH, polylines, points, layers | `layer`, `color` |
| VTK legacy | `.vtk` | polydata and unstructured grids | numeric point and cell data |
| VTK XML | `.vtu`, `.vtp` | as VTK, with zlib compression | numeric point and cell data |
| Datamine wireframe | `<name>_tr.csv` and `<name>_pt.csv` | triangles | extra point and triangle columns |

Text formats write floats in their shortest round-trip form, so a write and a read return the same bits. STL
stores 32-bit floats: at a northing of 7.5 million their spacing is half a metre, so keep STL for local
coordinates. Readers fan polygons from their first corner and record a face column `polygon`; writers that know
polygons (OBJ, PLY, OFF, VTK, VTU) rebuild them.

```python
data = dr.read("pit.dxf", crs="EPSG:31983")      # unit from $INSUNITS
mesh = dr.read_mesh("lens_tr.csv")               # either file of a Datamine pair
dr.write("lens.vtu", mesh, zlib=True)
dr.write("lens.dxf", mesh, dxf_entity="polyface", skip_unsupported_columns=True)
```

Options of `write`:

- `ascii=True` writes text PLY, STL and VTK.
- `dxf_entity` is `"3dface"` or `"polyface"`.
- `zlib=True` compresses VTU arrays.
- `skip_unsupported_columns=True` leaves out columns the format cannot store; without it they raise.

A `length_unit` argument that disagrees with the unit a DXF declares raises `InvalidInput`. Closed formats
(DWG, binary Datamine `.dm`) are out of scope: export DXF or a CSV pair from the source package.

## Example

![Orebody read back from a Datamine pair](../examples/01-files/files.png)

??? example "examples/01-files/example_01.py"

    ```python
    --8<-- "examples/01-files/example_01.py"
    ```
