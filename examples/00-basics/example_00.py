"""
# Meshes, strings, points and tables

`drusa` has three geometry containers. A `Mesh` holds triangles: an open surface such as topography, or a closed
solid such as an orebody. `Polylines` holds 3D strings grouped into features, and `PointSet` holds points. Each
carries attribute tables, a `crs` label and a `length_unit` label, and never changes once built: every method
returns a new object.
"""

# %% [hidden]
import sys
from pathlib import Path

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
sys.path.insert(0, str(HERE.parent))

# %%
import drusa as dr
import numpy as np
from common import save

# %% [markdown]
# ## A mesh from arrays
#
# A mesh takes `(n, 3)` vertices and `(m, 3)` triangles. A solid winds its triangles counter-clockwise seen from
# outside, so its signed volume is positive. Here a 20 m cube gets a face column naming its sides.

# %%
corners = np.array([[x, y, z] for z in (0, 20) for y in (0, 20) for x in (0, 20)], float)
faces = dr.datasets.box([0, 0, 0], [1, 1, 1]).triangles
sides = ["bottom", "bottom", "top", "top", "south", "south", "north", "north", "west", "west", "east", "east"]
cube = dr.Mesh(corners, faces, face_attributes={"side": sides}, crs="EPSG:31983", length_unit="m")
print(cube)
print(f"volume {cube.volume:.1f} m3, area {cube.area:.1f} m2, solid: {cube.is_solid}")
print(cube.to_table(on="faces").to_dict()["side"][:4])

# %% [markdown]
# ## Tables keep types and nulls
#
# Attribute tables hold float64, int64, bool and string columns. `None` marks a null; numeric nulls come back as
# NaN, and `mask` tells which rows hold a value.

# %%
table = dr.Table({"hole": ["DH1", "DH2", None], "depth": [120.0, None, 80.5], "assayed": [True, False, True]})
print(table)
print(table["depth"], table.mask("hole"))
print(table.drop_null().num_rows, "rows without nulls")

# %% [markdown]
# ## The synthetic deposit
#
# `dr.datasets` builds a seeded synthetic deposit at UTM scale: a lens-shaped orebody, a rolling topography and a
# fan of drillholes. The plan view draws the topography by elevation, the orebody outline at 200 m and the holes.

# %%
orebody = dr.datasets.orebody()
topography = dr.datasets.topography()
holes = dr.datasets.drillholes()
collars = dr.PointSet(
    np.array([p[0] for p in holes.parts]), {"hole": holes.attributes["hole"]}, length_unit="m"
)
print(orebody, holes, collars, sep="\n")

outline = orebody.slice_by_plane(dr.datasets.ORIGIN + [0, 0, 200], [0, 0, 1])
fig, ax = dr.plot.plan([topography, outline, holes, collars], label=False)
ax.set_title("Topography, orebody at 200 m and drillholes")
save(fig, "basics")

# %% [markdown]
# Every container converts to a table, and from there to pandas, polars or pyarrow when they are installed.

# %%
print(holes.to_table().column_names)
assert cube.volume == 8000.0
