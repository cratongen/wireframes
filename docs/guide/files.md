# Files

`wf.read(path)` returns a `Dataset` with `mesh`, `polylines`, `points`, `grid`, the `length_unit` the file
declares and the reader's `warnings`. `wf.read_mesh`, `wf.read_polylines` and `wf.read_points` return one part
and raise `InvalidInput` when the file has none. `wf.write(path, data)` takes a container, a grid or a `Dataset`;
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
import wireframes as wf

mesh = wf.datasets.orebody()
wf.write("lens.vtu", mesh, zlib=True)
wf.write("lens.dxf", mesh, dxf_entity="polyface", skip_unsupported_columns=True)
wf.write("lens_tr.csv", mesh)  # writes lens_tr.csv and lens_pt.csv
data = wf.read("lens.dxf", crs="EPSG:31983")  # unit from $INSUNITS
mesh = wf.read_mesh("lens_pt.csv")  # either file of a Datamine pair
```

Options of `write`:

- `ascii=True` writes text PLY, STL and VTK.
- `dxf_entity` is `"3dface"` or `"polyface"`.
- `zlib=True` compresses VTU arrays.
- `skip_unsupported_columns=True` leaves out columns the format cannot store; without it they raise.

## Files from other packages

Readers never repair or rescale what they read. Check `warnings`, then run `validate` on the mesh before you
measure it; [Measure, validate and repair](measure-repair.md) describes the report.

A DXF reads model space only. Paper space entities (title blocks, viewports) and `INSERT` block references are
skipped and counted in `warnings`, and so is the smoothing of a MESH, whose level-0 cage is what you get. A
`3DFACE` whose fourth corner repeats the third is a triangle; a true quad becomes two triangles with a face
column `polygon`. Formats that store corners per face (`3DFACE`, STL) are joined into a mesh by merging corners
with equal coordinates, so touching solids share their contact vertices, and corners a few millimetres apart
stay apart and show in `validate` as boundary edges. A DXF without a `HEADER` declares no unit: pass
`length_unit`.

A Datamine pair may use `,` or `;` between fields, with `,` as the decimal mark under `;`. The reader accepts
Windows line ends, a trailing delimiter, a row of field types under the header, gaps in `PID` and `-1.0E30` for
an absent value, which reads as null. Columns after `PID1 PID2 PID3` become face columns; a column reads as
text when any value is quoted or not numeric, so the writer quotes every text value and codes such as `"007"`
come back as text. The writer always uses `,` and `.`, and numbers `PID` from 1 in vertex order.

A truncated DXF (no `0`/`EOF` group), a triangle that names a `PID` missing from the point file, and a
non-finite coordinate raise `FileError` with the file and line. A `length_unit` argument that disagrees with
the unit a DXF declares raises `InvalidInput`; convert with `to_length_unit` after reading. Closed formats
(DWG, binary Datamine `.dm`) are out of scope: export DXF or a CSV pair from the source package.

## Keeping attributes across formats

| Target | Text columns | Numeric columns | Nulls |
|---|---|---|---|
| Datamine pair | kept | kept, read back as float | kept (`-1.0E30`) |
| DXF | `layer` only | `color` only | none |
| VTK, VTU, PLY | no | kept | kept (NaN or a `_valid` array) |
| OBJ | `object`, `group` | no | none |
| STL | no | no | none |

`write` raises `InvalidInput` naming a column the format cannot store. Copy what you need into a column it
does store, for example `mesh.with_face_column("layer", mesh["ZONE"])` before writing DXF, and pass
`skip_unsupported_columns=True` to drop the rest on purpose.

## Examples

![Orebody read back from a Datamine pair](../examples/01-files/files.png)

The full examples, with the output of each step: [example 01](../examples/01-files/README.md) writes and reads
every format, and [example 08](../examples/08-files-repair/README.md) reads a Datamine pair and a DXF as other
packages write them, repairs both and converts them while keeping their attributes.
