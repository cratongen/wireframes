"""
# Reading and writing files

`dr.read` returns everything a file holds as a `Dataset`: a mesh, polylines, points, a volume grid, the unit the
file declares and the reader's warnings. `dr.write` picks the format from the extension. Text formats print
floats in their shortest round-trip form, so coordinates survive a write and a read bit for bit.
"""

# %% [hidden]
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
sys.path.insert(0, str(HERE.parent))

# %%
import drusa as dr
import numpy as np
from common import save

# %% [markdown]
# ## One mesh, every format
#
# The orebody goes out to each format and comes back. STL and DXF 3DFACE store loose triangles, so their readers
# weld equal corners and number vertices in first-use order. Every format but STL returns the triangles bit for
# bit; STL stores 32-bit floats, whose spacing at a northing of 7.5 million is half a metre, so a vertex moves by
# up to 0.25 m. Keep STL for local coordinates. Formats that cannot store a column skip it when asked to.

# %%
orebody = dr.datasets.orebody()
original = orebody.coords[orebody.triangles]
folder = Path(tempfile.mkdtemp())
for name in [
    "lens.obj",
    "lens.ply",
    "lens.off",
    "lens.stl",
    "lens.dxf",
    "lens.vtk",
    "lens.vtu",
    "lens_tr.csv",
]:
    path = folder / name
    dr.write(path, orebody, skip_unsupported_columns=True)
    back = dr.read_mesh(path)
    shift = np.abs(back.coords[back.triangles] - original).max()
    print(
        f"{name:12} {path.stat().st_size:8,} bytes  largest shift {shift:.1e} m  columns {back.column_names}"
    )
    assert shift == 0.0 or (name.endswith(".stl") and shift <= 0.25)

# %% [markdown]
# ## Strings, points and labels
#
# Drillhole traces go to DXF as 3D polylines on a layer. `crs` and `length_unit` are labels set at read time; a
# DXF that declares `$INSUNITS` sets the unit itself, and a disagreeing argument raises `InvalidInput`.

# %%
holes = dr.datasets.drillholes()
dr.write(
    folder / "holes.dxf", holes.with_column("layer", holes.attributes["hole"]), skip_unsupported_columns=True
)
data = dr.read(folder / "holes.dxf", crs="EPSG:31983", length_unit="m")
print(data.polylines, data.polylines.crs, data.polylines.length_unit, data.warnings)
print(data.polylines.attributes["layer"][:3])

# %% [markdown]
# ## A Datamine pair
#
# A Datamine wireframe is two CSV files, `<name>_tr.csv` and `<name>_pt.csv`; either name reads the pair.

# %%
pair = dr.read_mesh(folder / "lens_pt.csv")
print((folder / "lens_tr.csv").read_text().splitlines()[:2])
print(pair)

# %% [markdown]
# The figure shows the mesh read back from the Datamine pair, coloured by the `section` column the loft wrote.

# %%
fig, ax = dr.plot.surface3d(pair, "section", azim=-50, elev=25)
ax.set_title("Orebody read back from lens_tr.csv")
save(fig, "files")
