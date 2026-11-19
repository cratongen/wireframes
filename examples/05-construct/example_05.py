"""
# Building surfaces and solids

`wf.dtm` triangulates survey points into a surface with constrained Delaunay, keeping breaklines as edges.
`wf.extrude`, `wf.loft` and `Mesh.close_to_elevation` turn strings and surfaces into solids, and `decimate`,
`smooth`, `remesh` and `offset` edit meshes. Each solid comes out validated.
"""

# %% [hidden]
import sys
from pathlib import Path

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
sys.path.insert(0, str(HERE.parent))

# %%
import matplotlib.pyplot as plt
import numpy as np
import wireframes as wf
from common import ACCENT, HIGHLIGHT, save

# %% [markdown]
# ## A surface with a breakline
#
# A road crest surveyed as a string must stay an edge of the surface. Without the breakline, triangles cross it
# and the surface cuts the corner of the road.

# %%
o = wf.datasets.ORIGIN
rng = np.random.default_rng(11)
xy = rng.uniform(-200, 200, (300, 2))
z = 300 + 0.05 * xy[:, 0] + 8 * np.tanh(xy[:, 1] / 30)
survey = wf.PointSet(o + np.column_stack([xy, z]), length_unit="m")
crest_x = np.linspace(-200, 200, 9)
crest = wf.Polylines([o + np.column_stack([crest_x, np.full(9, 10.0), 300 + 0.05 * crest_x + 3])])
plain = wf.dtm(survey)
road = wf.dtm(survey, breaklines=crest)
print(plain, road, sep="\n")

# %% [markdown]
# ## Contours and a solid down to a floor
#
# `contours` slices the surface at each level. `close_to_elevation` drops walls from the boundary of a
# height-field surface to a flat floor; its volume equals the fill of the surface over a flat plane at that
# elevation, which `cut_fill` computes independently.

# %%
lines = road.contours(np.arange(290.0, 315.0, 2.5))
print(lines)
solid = road.close_to_elevation(250.0)
lo, hi = road.bounds
corners = [[lo[0], lo[1], 250.0], [hi[0], lo[1], 250.0], [hi[0], hi[1], 250.0], [lo[0], hi[1], 250.0]]
floor = wf.dtm(np.array(corners))
fill = wf.cut_fill(floor, road)["fill"]
print(f"solid {solid.volume:,.1f} m3, cut_fill against the floor {fill:,.1f} m3")
assert solid.is_solid and abs(solid.volume - fill) <= 1e-8 * fill

# %% [markdown]
# ## Solids from strings
#
# `extrude` sweeps closed rings along a vector; `loft` links one section ring to the next, choosing the
# triangles of each band by dynamic programming to minimise their area.

# %%
pad = wf.Polylines([[[0, 0], [40, 0], [40, 25], [15, 40], [0, 25]]], closed=True, length_unit="m")
stockpile = wf.extrude(pad, [0, 0, 12])
print(f"extruded {stockpile.volume:,.1f} m3 = area {pad.area()[0]:,.1f} m2 x 12 m")
assert abs(stockpile.volume - pad.area()[0] * 12) <= 1e-9 * stockpile.volume
square = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], float)
column = wf.loft(wf.Polylines([square, square * [0.6, 0.6, 1] + [2, 2, 20]], closed=True))
frustum = 20 / 3 * (100 + 36 + np.sqrt(100 * 36))
print(f"lofted frustum {column.volume:,.4f} m3, analytic {frustum:,.4f} m3, solid: {column.is_solid}")
assert abs(column.volume - frustum) <= 1e-9 * frustum

# %% [markdown]
# `triangulate` fills each feature's closed rings with triangles, treating inner rings as holes by even-odd, so a
# pad with a sump becomes a flat mesh whose area is the polygon area.

# %%
sump = [[15, 10], [25, 10], [25, 20], [15, 20]]
holed_pad = wf.Polylines([pad.parts[0][:, :2], sump], closed=True, features=[0, 0])
floor = wf.triangulate(holed_pad)
print(f"triangulated {floor.n_faces} triangles, {floor.area:,.1f} m2 = {holed_pad.area()[0]:,.1f} m2")
assert abs(floor.area - holed_pad.area()[0]) <= 1e-9 * floor.area

# %% [markdown]
# ## Editing meshes
#
# `decimate` collapses edges by quadric error and keeps a solid closed; `smooth` applies Taubin's
# shrink-free smoothing with the boundary fixed; `remesh` evens out edge lengths; `offset` moves every vertex
# along its normal and refuses a result that intersects itself.

# %%
orebody = wf.datasets.orebody()
light = orebody.decimate(target_faces=150)
smooth = orebody.smooth(iterations=5)
even = orebody.remesh(15.0, iterations=3)
grown = orebody.offset(2.0)
for name, m in [("decimate", light), ("smooth", smooth), ("remesh", even), ("offset 2 m", grown)]:
    print(f"{name:11} {m.n_faces:5} triangles, volume {m.volume:12,.0f} m3, closed: {m.is_closed}")

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
wf.plot.plan(road, ax=ax1)
for part in crest.parts:
    ax1.plot(part[:, 0], part[:, 1], color=HIGHLIGHT, linewidth=1.6)
for part in lines.parts:
    ax1.plot(part[:, 0], part[:, 1], color="w", linewidth=0.6)
ax1.set_title("DTM with a breakline (orange) and contours")
wf.plot.section([orebody, light], o, [1, 0, 0], ax=ax2)
ring = even.slice_by_plane(o, [1, 0, 0])
for part, closed in zip(ring.parts, ring.closed, strict=True):
    part = np.vstack([part, part[:1]]) if closed else part
    ax2.plot(part[:, 1] - o[1], part[:, 2], color=ACCENT, linewidth=1.0)
ax2.set_title("Cross section: orebody, decimated (black) and remeshed (blue)")
save(fig, "construct")
