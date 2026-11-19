"""
# Booleans, clipping and cut and fill

`a | b`, `a & b`, `a - b` and `a ^ b` combine two solids by an exact mesh arrangement: every crossing between
triangles is decided by exact predicates on implicit intersection points, and coordinates are rounded to f64 only
when the result is written out. Inputs that are not solids are refused with the validation report.
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
from common import ACCENT, HIGHLIGHT, LIGHT, save

# %% [markdown]
# ## The orebody and a mining block
#
# The block is 100 m along strike, 160 m across and 100 m high. Each mesh carries a face column naming it, and
# the results keep it on the faces that came from each side.

# %%
o = wf.datasets.ORIGIN
orebody = wf.datasets.orebody().with_face_column("source", ["orebody"] * 428)
block = wf.datasets.box(o + [-50, -80, 150], o + [50, 80, 250]).with_face_column("source", ["block"] * 12)
parts = {
    "a | b": orebody | block,
    "a & b": orebody & block,
    "a - b": orebody - block,
    "a ^ b": orebody ^ block,
}
for name, part in parts.items():
    print(f"{name}: {part.volume:12,.1f} m3, {part.n_faces:5} triangles")

# %% [markdown]
# The volumes obey the inclusion-exclusion identities to rounding. These are the checks the test suite runs on
# random solids.

# %%
a, b = orebody.volume, block.volume
union, inter, diff, sym = (parts[k].volume for k in ("a | b", "a & b", "a - b", "a ^ b"))
print(f"union + intersection - a - b = {union + inter - a - b:+.2e} m3")
print(f"difference - (a - intersection) = {diff - (a - inter):+.2e} m3")
print(f"symmetric difference - (union - intersection) = {sym - (union - inter):+.2e} m3")
assert abs(union + inter - a - b) <= 1e-9 * (a + b)

# %% [markdown]
# Each operator has a method and a function form that give the same mesh. Where the two surfaces cross, the
# symmetric difference keeps the `a - b` and `b - a` pieces apart along the intersection curve, so `validate`
# reports duplicate vertices there; compute `a - b` and `b - a` separately when you need solids.

# %%
same = [orebody.symmetric_difference(block), wf.symmetric_difference(orebody, block)]
assert all(np.array_equal(m.triangles, parts["a ^ b"].triangles) for m in same)
print(parts["a ^ b"].validate().summary["duplicate_vertices"], "duplicate vertices along the curve")

# %% [markdown]
# A section along strike shows where each part of the outline came from.

# %%
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained", sharey=True)
colors = {"orebody": ACCENT, "block": HIGHLIGHT}
along = [0.0, 1.0, 0.0]
for ax, name in zip(axes, ["a & b", "a - b", "a | b"], strict=True):
    for ring in orebody.slice_by_plane(o, along).parts:
        ax.fill(ring[:, 0] - o[0], ring[:, 2], color=LIGHT)
    lines = parts[name].slice_by_plane(o, along)
    for f, (ring, closed) in enumerate(zip(lines.parts, lines.closed, strict=True)):
        ring = np.vstack([ring, ring[:1]]) if closed else ring
        ax.plot(ring[:, 0] - o[0], ring[:, 2], color=colors[lines.attributes["source"][f]], linewidth=1.4)
    ax.set_title(f"orebody {name.replace('a', '').replace('b', 'block').strip()}")
    ax.set_aspect("equal")
    ax.set_xlabel("Easting from origin (m)")
axes[0].set_ylabel("Elevation (m)")
save(fig, "booleans")

# %% [markdown]
# ## Cutting with a surface and a plane
#
# `clip_by_surface` splits a solid by an open surface and caps both parts with the pieces of the surface inside
# it. A dipping fault plane, built here as a small DTM, splits the orebody; `split_by_plane` does the same for an
# exact plane.

# %%
g = np.linspace(-300, 300, 7)
gx, gy = np.meshgrid(g, g)
fault = wf.dtm(o + np.column_stack([gx.ravel(), gy.ravel(), (200 + 0.4 * gx + 0.1 * gy).ravel()]))
below, above = orebody.clip_by_surface(fault)
print(
    f"below {below.volume:,.0f} + above {above.volume:,.0f} - orebody = {below.volume + above.volume - a:+.2e} m3"
)
west, east = orebody.split_by_plane(o + [25.0, 0, 0], [1.0, 0.0, 0.0])
print(
    f"west {west.volume:,.0f} + east {east.volume:,.0f} - orebody = {west.volume + east.volume - a:+.2e} m3"
)
assert below.is_solid and above.is_solid

# %% [markdown]
# ## Cut and fill
#
# `cut_fill` overlays two height-field surfaces in plan and integrates their difference exactly per overlay cell.
# A pit design against the topography gives the excavated volume as cut. The design is a 20 m grid that follows
# the topography outside the pit, so where it bridges a hollow of the topography it sits above it, and that shows
# as fill.

# %%
topography = wf.datasets.topography()
pit = wf.datasets.pit(topography)
volumes = wf.cut_fill(topography, pit)
print({k: round(v) for k, v in volumes.items()})
assert abs(volumes["net"] - (volumes["fill"] - volumes["cut"])) <= 1e-6 * volumes["cut"]
