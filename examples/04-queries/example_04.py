"""
# Points, rays, drillholes and blocks against a solid

`locate` tells whether each point is inside, on or outside a solid, exactly: it counts crossings of a vertical
ray with exact predicates and breaks ties by symbolic perturbation. `signed_distance`, `ray_hits`,
`trace_intervals` and `block_fraction` build on the same index, which the mesh builds once and keeps.
"""

# %% [hidden]
import sys
from pathlib import Path

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
sys.path.insert(0, str(HERE.parent))

# %%
import drusa as dr
import matplotlib.pyplot as plt
import numpy as np
from common import HIGHLIGHT, save

# %% [markdown]
# ## Inside, on or outside
#
# The sign of `locate` (-1 inside, 0 on the surface, 1 outside) agrees with the sign of `signed_distance`. A
# vertex of the mesh is on its surface.

# %%
o = dr.datasets.ORIGIN
orebody = dr.datasets.orebody()
rng = np.random.default_rng(3)
points = o + rng.uniform([-220, -80, 160], [220, 80, 240], (2000, 3))
where = orebody.locate(points)
distance = orebody.signed_distance(points)
print(np.unique(where, return_counts=True))
assert np.array_equal(np.sign(distance).astype(np.int8), where)
assert (orebody.locate(orebody.coords) == 0).all()

# %% [markdown]
# ## Drillhole intervals
#
# `trace_intervals` splits a desurveyed trace into runs inside and outside the solid, with the faces where it
# enters and leaves. Depths interpolate along each trace segment.

# %%
holes = dr.datasets.drillholes()
rows = []
for name, trace in zip(holes.attributes["hole"], holes.parts, strict=True):
    iv = orebody.trace_intervals(trace, [0.0, 350.0])
    for start, end in zip(iv["from"][iv["inside"]], iv["to"][iv["inside"]], strict=True):
        rows.append((name, start, end))
        print(f"{name}: inside from {start:6.2f} m to {end:6.2f} m ({end - start:5.2f} m)")

# %% [markdown]
# ## Block fractions
#
# `block_fraction` integrates the solid over each block with the divergence theorem, so there is no sampling
# error: over blocks that cover the solid, fraction times block volume sums to the solid volume.

# %%
size = np.array([20.0, 10.0, 5.0])
lo, hi = (np.array(b) for b in orebody.bounds)
axes = [np.arange(lo[k] + size[k] / 2 - size[k], hi[k] + size[k], size[k]) for k in range(3)]
centres = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
fraction = orebody.block_fraction(centres, size)
print(
    f"{len(centres):,} blocks, {np.count_nonzero(fraction):,} touch the solid, {np.sum(fraction == 1):,} full"
)
print(f"sum of fraction x volume - solid volume = {fraction.sum() * size.prod() - orebody.volume:+.2e} m3")
assert abs(fraction.sum() * size.prod() - orebody.volume) <= 1e-8 * orebody.volume

# %% [markdown]
# The top panel shows the signed distance on a horizontal plane at 200 m, the bottom panel the block fractions
# of the 200 m bench with the drillhole intercepts.

# %%
gx, gy = np.meshgrid(np.linspace(-240, 240, 241), np.linspace(-90, 90, 91))
plane = o + np.column_stack([gx.ravel(), gy.ravel(), np.full(gx.size, 200.0)])
field = orebody.signed_distance(plane).reshape(gx.shape)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6.4), layout="constrained", sharex=True)
image = ax1.pcolormesh(gx, gy, field, cmap="RdBu", vmin=-30, vmax=30, shading="auto")
ax1.contour(gx, gy, field, [0.0], colors="k", linewidths=1.0)
fig.colorbar(image, ax=ax1, label="Signed distance (m)")
ax1.set_title("Signed distance at 200 m")
bench = np.isclose(centres[:, 2], axes[2][np.argmin(np.abs(axes[2] - (o[2] + 200)))])
c = centres[bench] - o
art = ax2.scatter(c[:, 0], c[:, 1], c=fraction[bench], marker="s", s=24, cmap="cividis")
fig.colorbar(art, ax=ax2, label="Fraction inside")
for trace in holes.parts:
    t = trace - o
    ax2.plot(t[:, 0], t[:, 1], color=HIGHLIGHT, linewidth=0.8)
ax2.set_title("Block fractions of the 200 m bench and drillhole traces")
for ax in (ax1, ax2):
    ax.set_aspect("equal")
    ax.set_ylabel("Northing from origin (m)")
ax2.set_xlabel("Easting from origin (m)")
save(fig, "queries")
