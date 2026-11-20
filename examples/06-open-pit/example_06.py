"""
# Open pit: topography, pit stages, bench tonnes and a survey

This example follows one open-pit workflow from survey data to a reconciliation, with every number checked by an
independent route. It builds a topographic surface from contours and a breakline, closes a pit design into a
solid and clips it by the topography, splits the pit into mining stages by boolean difference, reports volumes and
tonnes per bench, and compares an end-of-stage survey with the design.

All data are seeded and synthetic, at UTM scale: easting near 700 000 m and northing near 7 500 000 m. Lengths are
in metres (m), areas in square metres (m2), volumes in cubic metres (m3), densities in tonnes per cubic metre
(t/m3) and masses in tonnes (t), with kt for 1 000 t and Mt for 1 000 000 t.
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
from common import ACCENT, GRAY, HIGHLIGHT, INK, LIGHT, SERIES, save
from matplotlib.patches import Patch

O = wf.datasets.ORIGIN
UP = [0.0, 0.0, 1.0]
NORTH = [0.0, 1.0, 0.0]

# %% [markdown]
# ## 1. Topography from contours and a breakline
#
# A digital terrain model (DTM) is a surface of triangles that stands for the ground. Each point of the plan has
# one elevation on it, which makes the DTM a *height field*. Surveyors rarely hand over a DTM. They hand over
# *contours*, strings that join points of equal elevation (here every 5 m), *spot heights*, single surveyed
# points, and *breaklines*, strings along a sharp change of slope such as a creek bed, a road edge or a crest.
#
# The terrain below is a hill crossed in the north by a creek. The function `ground` plays the part of the true
# ground: the example takes its contours, samples the creek bed along its thalweg (the line of lowest points),
# draws 60 spot heights from it plus one at each corner of the 1 km square, and then leaves it alone until the
# accuracy check.

# %%
HALF = 500.0
PIT_CENTRE = O + [40.0, -30.0, 0.0]


def ground(x, y):
    """Elevation in m of the synthetic terrain at plan offsets x, y in m from the origin."""
    hill = 45.0 * np.exp(-(((x - 40) / 280) ** 2 + ((y + 30) / 220) ** 2))
    creek = 14.0 * np.exp(-(((y - 0.35 * x - 300) / 45) ** 2))
    return 290.0 + hill + 0.02 * x - creek


g = np.linspace(-HALF, HALF, 201)
gx, gy = np.meshgrid(g, g)
reference = wf.dtm(O + np.column_stack([gx.ravel(), gy.ravel(), ground(gx, gy).ravel()]))
contours = reference.contours(np.arange(275.0, 340.0, 5.0))
creek_x = np.linspace(-HALF, 230.0, 30)
creek_xy = np.column_stack([creek_x, 0.35 * creek_x + 300])
creek = wf.Polylines([O + np.column_stack([creek_xy, ground(*creek_xy.T)])], length_unit="m")
rng = np.random.default_rng(7)
corners = np.array([[-HALF, -HALF], [HALF, -HALF], [HALF, HALF], [-HALF, HALF]])
spot_xy = np.vstack([rng.uniform(-HALF, HALF, (60, 2)), corners])
spots = wf.PointSet(O + np.column_stack([spot_xy, ground(*spot_xy.T)]), length_unit="m")
print(contours)
print(f"{len(spots)} spot heights, creek breakline {creek.length()[0]:,.0f} m long")

# %% [markdown]
# `wf.dtm` triangulates the spot heights by constrained Delaunay triangulation in plan, and passes the contours
# and the creek as breaklines, so every segment of every string becomes an edge of the surface. Where the creek
# crosses a contour, the new vertex takes its elevation from the string inserted first; the contours come first
# because they carry exact elevations. A second DTM without the creek shows what the breakline buys.

# %%
strings = wf.Polylines(contours.parts + creek.parts, closed=[*contours.closed, False], length_unit="m")
topography = wf.dtm(spots, breaklines=strings)
no_creek = wf.dtm(spots, breaklines=contours)
print(topography)
summary = topography.validate().summary
print(f"closed: {summary['is_closed']}, self-intersections: {summary['self_intersections']}")

# %% [markdown]
# Two checks test the surface. The DTM must pass through every vertex of every string, with no interpolation,
# and it must stay close to the true ground away from the strings. The second check samples 5 000 points.

# %%
on_strings = topography.surface_z(contours.coords) - contours.coords[:, 2]
print(f"largest misfit at contour vertices: {np.abs(on_strings).max()} m")
assert np.array_equal(topography.surface_z(contours.coords), contours.coords[:, 2])
check_xy = rng.uniform(-0.9 * HALF, 0.9 * HALF, (5000, 2))
for name, surface in [("with creek", topography), ("without creek", no_creek)]:
    error = surface.surface_z(O[:2] + check_xy) - ground(*check_xy.T)
    print(f"DTM {name:13}: RMS error {np.sqrt(np.mean(error**2)):.2f} m, largest {np.abs(error).max():.2f} m")

# %% [markdown]
# The plan shows the DTM coloured by elevation with its source strings. The section cuts across the creek: without
# the breakline, triangles bridge the creek bed between the two contours on its banks, and the surface floats
# above the true thalweg.

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8), layout="constrained")
wf.plot.plan(topography, ax=ax1)
for part in contours.parts:
    ax1.plot(part[:, 0], part[:, 1], color="w", linewidth=0.5)
ax1.plot(creek.x, creek.y, color=HIGHLIGHT, linewidth=1.6, label="creek breakline")
ax1.scatter(spots.x, spots.y, s=6, color=INK, label="spot heights")
ax1.legend(loc="lower left")
ax1.set_title("DTM coloured by elevation, contours every 5 m")
ax1.set_xlabel("Easting (m)")
ax1.set_ylabel("Northing (m)")

across = np.array([1.0, 0.35, 0.0]) / np.hypot(1.0, 0.35)
station = O + [0.0, 300.0, 0.0]
side = np.array([-across[1], across[0], 0.0])
for surface, colour, label in [(topography, ACCENT, "DTM with creek"), (no_creek, HIGHLIGHT, "DTM without")]:
    for part in surface.slice_by_plane(station, across).parts:
        ax2.plot((part - station) @ side, part[:, 2], color=colour, linewidth=1.3, label=label)
s = np.linspace(-150, 150, 301)
true_xy = (station[:2] - O[:2]) + np.outer(s, side[:2])
ax2.plot(s, ground(*true_xy.T), color=GRAY, linestyle="--", linewidth=1.0, label="true ground")
ax2.set_xlim(-150, 150)
ax2.legend(loc="lower right")
ax2.set_title("Section across the creek at easting 700 000 m")
ax2.set_xlabel("Distance across the creek, positive to the north-north-west (m)")
ax2.set_ylabel("Elevation (m)")
save(fig, "open_pit_topography")

# %% [markdown]
# ## 2. The pit design and its stages
#
# A *pit shell* is the surface of the final excavation. A pit wall rises in *benches*: each bench is a 10 m
# *batter* (the steep face, here 5 m wide in plan, about 63 degrees) whose foot is the *toe* and whose top edge is
# the *crest*, followed by a 6 m *berm*, the flat ledge that catches rock falling from above. Batter and berm
# together give an overall wall slope of about 42 degrees. The design below has a 40 m floor radius at 220 m
# elevation, 12 benches, and a last 60 m wide ledge above the ground, so the shell runs past the line where it
# meets the topography, the *daylight line*.
#
# The pit is mined in three *stages*: a starter pit, then two *pushbacks* (also called cutbacks), each of which
# widens the pit to the east. Each stage's shell is the starter shell stretched east, so all three share the west
# wall vertex for vertex.


# %%
def benched_shell(stretch, *, floor=220.0, toe=40.0, benches=12, n=64):
    """Pit design surface: a polygonal floor at `floor` of radius `toe`, then `benches` 10 m benches with 5 m
    batters and 6 m berms. Each ring is stretched east by ``1 + stretch * max(cos, 0)``."""
    angle = np.linspace(0, 2 * np.pi, n, endpoint=False)
    unit = np.column_stack([np.cos(angle), np.sin(angle)])
    k = 1.0 + stretch * np.maximum(unit[:, 0], 0.0)
    rings, r, z = [(toe, floor)], toe, floor
    for b in range(benches):
        r, z = r + 5.0, z + 10.0
        rings.append((r, z))
        r += 6.0 if b < benches - 1 else 60.0
        rings.append((r, z))
    v = [[0.0, 0.0, floor]] + [
        [r * k[i] * unit[i, 0], r * k[i] * unit[i, 1], z] for r, z in rings for i in range(n)
    ]
    f = [[0, 1 + i, 1 + (i + 1) % n] for i in range(n)]
    for j in range(len(rings) - 1):
        a, b = 1 + j * n, 1 + (j + 1) * n
        for i in range(n):
            i1 = (i + 1) % n
            f.extend([[a + i, b + i, b + i1], [a + i, b + i1, a + i1]])
    faces = {"surface": ["design"] * len(f)}
    return wf.Mesh(np.asarray(v) + PIT_CENTRE, np.asarray(f), face_attributes=faces, length_unit="m")


shells = [benched_shell(s) for s in (0.0, 0.3, 0.6)]
topography = topography.with_face_column("surface", ["topography"] * topography.n_faces)
for k, shell in enumerate(shells, start=1):
    print(f"stage {k} shell: {shell.n_faces} triangles, rim at {shell.bounds[1][2]:.0f} m")
print(f"highest ground: {topography.bounds[1][2]:.1f} m")

# %% [markdown]
# A shell is an open surface and has no volume. `close_to_elevation(400)` turns it into a *solid*, a closed and
# consistently wound set of triangles: it drops walls from the shell's rim to 400 m and closes them with a flat
# top. The solid holds the air above the shell up to 400 m, higher than any ground.
#
# `clip_by_surface(topography)` then splits that solid along the topography into the part below the ground, which
# is the rock to excavate, and the part above it. Both parts come back capped by the piece of topography inside
# the solid, and each cap face keeps the topography's `surface` column.

# %%
closed = [shell.close_to_elevation(400.0) for shell in shells]
pits = [solid.clip_by_surface(topography)[0] for solid in closed]
for k, pit in enumerate(pits, start=1):
    print(f"stage {k} pit: {pit.volume:12,.0f} m3, solid: {pit.is_solid}")

# %% [markdown]
# `wf.cut_fill(base, design)` gives an independent check. It overlays two height fields in plan and integrates
# their difference over every cell of the overlay; *cut* is the volume where the design lies below the base,
# *fill* where it lies above. The cut of each shell against the topography must equal the clipped pit's volume.
# The two routes are independent: one builds a 3D mesh arrangement and measures a closed solid, the other
# integrates heights in plan.

# %%
for k, (shell, pit) in enumerate(zip(shells, pits, strict=True), start=1):
    cut = wf.cut_fill(topography, shell)["cut"]
    print(
        f"stage {k}: clipped pit {pit.volume:,.1f} m3, cut {cut:,.1f} m3, difference {pit.volume - cut:+.1e} m3"
    )
    assert pit.is_solid and abs(pit.volume - cut) <= 1e-9 * cut

# %% [markdown]
# Each stage's own volume is the difference between consecutive closed shells, clipped by the topography.
# Subtracting the clipped pits gives the same stage as a solid too: section 5 shows what the rounding does there.

# %%
stages = [closed[0]] + [closed[k] - closed[k - 1] for k in (1, 2)]
stages = [solid.clip_by_surface(topography)[0] for solid in stages]
for k, stage in enumerate(stages, start=1):
    print(f"stage {k}: {stage.volume:12,.0f} m3, solid: {stage.is_solid}")
total = sum(stage.volume for stage in stages)
print(
    f"sum of stages {total:,.1f} m3, final pit {pits[2].volume:,.1f} m3, difference {total - pits[2].volume:+.1e}"
)
assert all(stage.is_solid for stage in stages)
assert abs(total - pits[2].volume) <= 1e-9 * pits[2].volume

# %% [markdown]
# A plan, an east-west section through the pit centre and a 3D view show the three stages. The plan draws each
# stage's daylight line. It is the outer edge of the pit's `design` faces, and the function `daylight` walks the
# edges that only one triangle uses and returns them as one closed ring.


# %%
def daylight(surface):
    """The boundary of an open surface with one boundary loop, as a closed ring."""
    t = surface.triangles
    edges = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    _, index, count = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    following = dict(edges[count[index] == 1].tolist())
    ring = [next(iter(following))]
    while following[ring[-1]] != ring[0]:
        ring.append(following[ring[-1]])
    assert len(ring) == len(following), "the surface has more than one boundary loop"
    return wf.Polylines([surface.coords[ring]], closed=True, length_unit="m")


designs = [pit.filter_faces(pit["surface"] == "design") for pit in pits]
outlines = [daylight(d) for d in designs]
for k, line in enumerate(outlines, start=1):
    print(f"stage {k} daylight line: {line.length()[0]:,.0f} m long, enclosing {line.area()[0]:,.0f} m2")


# %%
def plain(mesh):
    """The mesh without its face columns, so a section through it comes out as whole closed rings."""
    return wf.Mesh(mesh.coords, mesh.triangles, length_unit="m")


fig = plt.figure(figsize=(12, 9.2), layout="constrained")
top, bottom = fig.subfigures(2, 1, height_ratios=[1.0, 1.0])
ax1, ax2 = top.subplots(1, 2, width_ratios=[1.0, 1.3])
for part in topography.contours(np.arange(270.0, 340.0, 5.0)).parts:
    ax1.plot(part[:, 0], part[:, 1], color=LIGHT, linewidth=0.6)
for part in plain(pits[2]).slice_by_plane(O, UP, offsets=np.arange(225.0, 340.0, 10.0)).parts:
    ax1.plot(*np.vstack([part, part[:1]])[:, :2].T, color=GRAY, linewidth=0.4)
for k, line in enumerate(outlines):
    ring = np.vstack([line.parts[0], line.parts[0][:1]])
    ax1.plot(ring[:, 0], ring[:, 1], color=SERIES[k], linewidth=1.8, label=f"stage {k + 1} daylight")
ax1.axhline(PIT_CENTRE[1], color=INK, linewidth=0.6, linestyle=":")
ax1.set_aspect("equal")
ax1.set_xlim(PIT_CENTRE[0] - 260, PIT_CENTRE[0] + 330)
ax1.set_ylim(PIT_CENTRE[1] - 230, PIT_CENTRE[1] + 230)
ax1.ticklabel_format(useOffset=False, style="plain")
ax1.locator_params(nbins=5)
ax1.legend(loc="lower left", fontsize=7)
ax1.set_title("Plan: final pit cut halfway up each batter")
ax1.set_xlabel("Easting (m)")
ax1.set_ylabel("Northing (m)")

for k, stage in enumerate(stages):
    for ring in plain(stage).slice_by_plane(PIT_CENTRE, NORTH).parts:
        ax2.fill(ring[:, 0] - PIT_CENTRE[0], ring[:, 2], color=SERIES[k], label=f"stage {k + 1}")
for part in topography.slice_by_plane(PIT_CENTRE, NORTH).parts:
    ax2.plot(part[:, 0] - PIT_CENTRE[0], part[:, 2], color=INK, linewidth=1.0, label="topography")
handles, labels = ax2.get_legend_handles_labels()
unique = dict(zip(labels, handles, strict=True))
ax2.legend(unique.values(), unique.keys(), loc="lower right", fontsize=7)
ax2.set_aspect("equal")
ax2.set_xlim(-200, 280)
ax2.set_ylim(200, 350)
ax2.set_title("Section along the dotted line, looking north")
ax2.set_xlabel("Easting from the pit centre (m)")
ax2.set_ylabel("Elevation (m)")

ax3 = bottom.add_subplot(projection="3d")
labelled = [
    plain(stage).translate(-PIT_CENTRE).with_face_column("stage", [f"stage {k + 1}"] * stage.n_faces)
    for k, stage in enumerate(stages)
]
wf.plot.surface3d(wf.concat(labelled), "stage", ax=ax3, azim=-60, elev=-25)
swatches = [
    Patch(color=plt.get_cmap("viridis")(c), label=f"stage {k + 1}")
    for k, c in enumerate([1 / 6, 1 / 2, 5 / 6])
]
ax3.legend(handles=swatches, loc="upper left", fontsize=7)
ax3.set_title("The stages seen from below, so the benches show; stage 2 is the thin band between the others")
save(fig, "open_pit_stages")

# %% [markdown]
# ## 3. Volumes and tonnes by bench
#
# Rock near the surface is weathered. The *base of oxidation* is the surface below which rock is fresh; above it
# lies *oxide*, softer and lighter. Here it sits 12 m to 20 m below the topography. The *in situ dry bulk density*
# converts volume to mass: 2.0 t/m3 for oxide and 2.7 t/m3 for fresh rock (synthetic values).
#
# `clip_by_surface` splits each stage at the base of oxidation, and `volume_by_elevation` integrates each piece
# exactly between the bench elevations, every 10 m from the pit floor. Bench 220 holds the rock from 220 m up to
# 230 m. The result is a table with one row per stage, material and bench.

# %%
DENSITY = {"oxide": 2.0, "fresh": 2.7}
x, y = (topography.coords[:, :2] - PIT_CENTRE[:2]).T
depth = 16.0 + 4.0 * np.sin(x / 90.0) * np.cos(y / 120.0)
oxidation = wf.Mesh(
    topography.coords - np.column_stack([0 * x, 0 * y, depth]), topography.triangles, length_unit="m"
)
edges = np.arange(220.0, 350.0, 10.0)
rows = []
for k, stage in enumerate(stages, start=1):
    fresh, oxide = stage.clip_by_surface(oxidation)
    assert oxide.is_solid and fresh.is_solid
    assert abs(oxide.volume + fresh.volume - stage.volume) <= 1e-9 * stage.volume
    for material, piece in [("oxide", oxide), ("fresh", fresh)]:
        volumes = piece.volume_by_elevation(edges)
        for bench, volume in zip(edges[:-1], volumes, strict=True):
            rows.append((k, material, bench, volume, volume * DENSITY[material]))
table = wf.Table(
    {
        "stage": [r[0] for r in rows],
        "material": [r[1] for r in rows],
        "bench": [r[2] for r in rows],
        "volume": [r[3] for r in rows],
        "tonnes": [r[4] for r in rows],
    }
)
print(table)

# %% [markdown]
# The bench report lists kilotonnes per bench and stage, highest bench first. A bench that a piece does not reach
# can show `-0.0`: the exact integral is zero there, and rounding leaves a few 1e-13 m3 of either sign.
#
# The balance check below the report is the point of the exercise. Per stage, the benches must add up to the
# stage volume; the stages must add up to the final pit; and the final pit must match the cut of the final shell
# against the topography.

# %%
stage_col, bench_col, mat_col = table["stage"], table["bench"], table["material"]
volume_col, tonnes_col = table["volume"], table["tonnes"]
print(
    f"{'bench (m)':>11} "
    + " ".join(f"{f'stage {k} oxide':>14} {f'stage {k} fresh':>14}" for k in (1, 2, 3))
    + f" {'total':>9}"
)
for bench in edges[-2::-1]:
    cells = []
    for k in (1, 2, 3):
        for material in ("oxide", "fresh"):
            pick = (stage_col == k) & (bench_col == bench) & (mat_col == material)
            cells.append(tonnes_col[pick].sum() / 1e3)
    if sum(cells) > 0:
        print(
            f"{bench:5.0f}-{bench + 10:3.0f}  "
            + " ".join(f"{c:14,.1f}" for c in cells)
            + f" {sum(cells):9,.1f}"
        )
print("kt per bench and stage")

for k, stage in enumerate(stages, start=1):
    pick = stage_col == k
    print(
        f"stage {k}: benches {volume_col[pick].sum():12,.1f} m3, solid {stage.volume:12,.1f} m3, "
        f"{tonnes_col[pick].sum() / 1e6:6.3f} Mt"
    )
    assert abs(volume_col[pick].sum() - stage.volume) <= 1e-9 * stage.volume
    assert volume_col[pick].min() >= -1e-9 * stage.volume
final_cut = wf.cut_fill(topography, shells[2])["cut"]
print(f"all benches {volume_col.sum():,.1f} m3, final shell cut {final_cut:,.1f} m3")
print(
    f"total {tonnes_col.sum() / 1e6:.3f} Mt, of which oxide {tonnes_col[mat_col == 'oxide'].sum() / 1e6:.3f} Mt"
)
assert abs(volume_col.sum() - final_cut) <= 1e-9 * final_cut

# %%
fig, ax = plt.subplots(figsize=(8, 4.6), layout="constrained")
bench_mid = edges[:-1] + 5.0
left = np.zeros(len(bench_mid))
for k in (1, 2, 3):
    for material, hatch in [("oxide", "///"), ("fresh", None)]:
        pick = (stage_col == k) & (mat_col == material)
        kt = tonnes_col[pick] / 1e3
        ax.barh(bench_mid, kt, height=8.5, left=left, color=SERIES[k - 1], hatch=hatch, edgecolor="w", linewidth=0.3,
                label=f"stage {k} {material}")  # fmt: skip
        left += kt
ax.set_yticks(edges)
ax.set_xlabel("Mass per bench (kt)")
ax.set_ylabel("Bench elevation (m)")
ax.set_title("Tonnes by bench, stage and material (oxide hatched)")
ax.legend(loc="lower right", ncols=3, fontsize=7)
save(fig, "open_pit_benches")

# %% [markdown]
# ## 4. Cut and fill against a survey
#
# When the starter pit is finished, the mine surveyors pick up the pit floor and walls. Comparing that survey with
# the design gives two numbers: *overbreak*, rock taken below the design surface, and *underbreak*, rock left
# above it. In `cut_fill(design, survey)` they are the cut and the fill.
#
# The design surface to compare is the part of the stage 1 pit that came from the shell, its `design` faces from
# section 2. Its outer edge, the stage 1 daylight line, bounds the survey.

# %%
design, outline = designs[0], outlines[0]
print(design)

# %% [markdown]
# The synthetic survey has 6 000 pickups inside the daylight line, each on the design with 0.3 m of survey noise.
# Two deviations are built in: the floor blast in the west broke 0.8 m below grade, and the last bench in the
# north-east is still standing, 10 m above the floor. The survey surface is a DTM of the pickups bounded by the
# daylight line.

# %%
rng = np.random.default_rng(42)
lo, hi = outline.bounds
xy = rng.uniform(lo[:2], hi[:2], (12_000, 2))
xy = xy[outline.contains(xy)][:6000]
z = design.surface_z(xy) + rng.normal(0.0, 0.3, len(xy))
dx, dy = (xy - PIT_CENTRE[:2]).T
z = np.where((dx < -10) & (z < 222.0), z - 0.8, z)
z = np.where((dx > 15) & (dy > 15), np.maximum(z, 230.0), z)
survey = wf.dtm(wf.PointSet(np.column_stack([xy, z]), length_unit="m"), boundary=outline)
print(survey)

# %% [markdown]
# Three overlays, all over the same plan domain inside the daylight line:
#
# - topography to survey: the volume mined, as cut;
# - topography to design: the volume designed, which must equal the stage 1 pit solid;
# - design to survey: overbreak (cut) and underbreak (fill).
#
# Net volume is fill minus cut. Integration is linear, so the nets must add up:
# net(topography to survey) = net(topography to design) + net(design to survey).

# %%
mined = wf.cut_fill(topography, survey, boundary=outline)
designed = wf.cut_fill(topography, design, boundary=outline)
versus = wf.cut_fill(design, survey, boundary=outline)
print(f"mined      {mined['cut']:12,.0f} m3 cut, {mined['fill']:8,.0f} m3 fill")
print(f"designed   {designed['cut']:12,.0f} m3 cut; stage 1 pit {pits[0].volume:,.0f} m3")
print(f"overbreak  {versus['cut']:12,.0f} m3 over {versus['cut_area']:,.0f} m2")
print(f"underbreak {versus['fill']:12,.0f} m3 over {versus['fill_area']:,.0f} m2")
balance = mined["net"] - designed["net"] - versus["net"]
print(f"net balance {balance:+.1e} m3 against {designed['cut']:,.0f} m3 designed")
assert abs(designed["cut"] - pits[0].volume) <= 1e-9 * pits[0].volume
assert abs(balance) <= 1e-9 * designed["cut"]

# %% [markdown]
# The small fill of the mined overlay is survey noise near the crest, where a few pickups land above the ground.
# The map colours each survey vertex by its height above the design: red for underbreak, blue for
# overbreak.

# %%
survey = survey.with_vertex_column("survey_minus_design", design.vertical_distance(survey.coords), unit="m")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), layout="constrained", width_ratios=[1.0, 1.3])
wf.plot.plan([survey, outline], "survey_minus_design", ax=ax1, cmap="RdBu_r")
ax1.collections[0].set_clim(-2, 2)
ax1.set_title("Survey minus design (m)")
ax1.set_xlabel("Easting (m)")
ax1.set_ylabel("Northing (m)")
for surface, colour, label, width in [(design, GRAY, "design", 2.5), (survey, ACCENT, "survey", 1.0)]:
    for part in surface.slice_by_plane(PIT_CENTRE + [0.0, 25.0, 0.0], NORTH).parts:
        ax2.plot(part[:, 0] - PIT_CENTRE[0], part[:, 2], color=colour, linewidth=width, label=label)
ax2.legend(loc="upper left")
ax2.set_title("Section 25 m north of the pit centre, looking north")
ax2.set_xlabel("Easting from the pit centre (m)")
ax2.set_ylabel("Elevation (m)")
save(fig, "open_pit_survey")

# %% [markdown]
# ## 5. What can go wrong
#
# The library refuses bad geometry instead of guessing. Two such cases come up in open-pit work, and a third case,
# surfaces that agree only to rounding, shows what the library does where the exact result is thinner than a
# 64-bit number can place.
#
# **A hole in a solid.** A pit solid that lost triangles on export, here three faces of the stage 2 pushback, is
# open. `validate()` returns a report with a count per problem and one row per offending edge or face, and every
# volume, boolean or clipping function raises `wf.InvalidMesh` carrying that report. `repair(fill_holes="fan")`
# closes simple holes, but only when you call it, and its report lists each change.

# %%
holed = stages[1].filter_faces(np.arange(stages[1].n_faces) >= 3)
report = holed.validate()
print({k: v for k, v in report.summary.items() if v and k != "shells"})
hole = report.problems.filter(report.problems["kind"] == "boundary_edge")
print(f"faces along the hole: {hole['face'].tolist()}")
try:
    holed.clip_by_surface(oxidation)
except wf.InvalidMesh as error:
    print(f"{type(error).__name__}: {error}")
    assert error.report.summary["boundary_edges"] == report.summary["boundary_edges"] > 0
else:
    raise AssertionError("an open pit solid must be refused")

# %% [markdown]
# **A spot height that disagrees with a contour.** A spot height at the plan position of a contour vertex but at
# another elevation gives two heights for one point of the plan. `wf.dtm` raises `wf.InvalidInput` and names both
# points instead of picking one. Correct or remove the spot height; the contour carries the mapped elevation.

# %%
clash = contours.coords[0] + [0.0, 0.0, 1.5]
try:
    wf.dtm(np.vstack([spots.coords, clash]), breaklines=strings)
except wf.InvalidInput as error:
    print(f"{type(error).__name__}: {error}")
else:
    raise AssertionError("a spot height off its contour must be refused")

# %% [markdown]
# **Surfaces that agree only to rounding.** Exact predicates decide every crossing, but the output vertices are
# written as 64-bit floating-point numbers. Subtracting the clipped stage 1 pit from the clipped stage 2 pit meets
# the limit of that rounding: both pits carry a cap cut from the same topography triangles, and where the caps
# overlap they agree only to the last digit. The exact difference holds slivers thinner than a 64-bit number can
# place, and rounding them would fold them into their neighbours. So the boolean first inserts each pit's cap
# vertices into the other pit's cap where they lie within rounding of it, and flips edges until the two caps match
# there. The caps then coincide exactly, the slivers vanish, and the difference is a solid. No input vertex moves.
# `rounding_report=True` returns what the rounding did and the proven bound on how far the volume can be from the
# exact result's.

# %%
pushback, rounding = wf.difference(pits[1], pits[0], rounding_report=True)
assert isinstance(rounding, wf.RoundingReport) and rounding.moved_input_vertices == 0
assert pushback.is_solid
assert abs(pushback.volume - stages[1].volume) <= 1e-9 * stages[1].volume
print(f"solid: {pushback.is_solid}, {pushback.volume:,.1f} m3 against stage 2 at {stages[1].volume:,.1f} m3")
print(
    f"{rounding.inserted_points} points inserted, {rounding.flipped_edges} edges flipped, "
    f"{rounding.moved_vertices} vertices moved; volume within {rounding.volume_bound:.1e} m3 of the exact result"
)
print(
    f"{rounding.rounded_vertices} rounded vertices, {rounding.merged_vertices} merged, "
    f"{rounding.collapsed_edges} collapsed, {rounding.removed_shells} shells removed, "
    f"{rounding.cancelled_faces} faces cancelled; largest move {max(rounding.max_move):.1e} m "
    f"within {max(rounding.move_bound):.1e} m; {rounding.steps} steps of at most {rounding.step_bound}"
)
