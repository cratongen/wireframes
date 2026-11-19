import numpy as np
import pytest
import wireframes as wf
from conftest import close


def test_box_measures(box, rng):
    for _ in range(10):
        lo = rng.uniform(-5, 5, 3)
        size = rng.uniform(0.1, 4, 3)
        b = box(lo, lo + size)
        dx, dy, dz = size
        assert close(b.volume, dx * dy * dz)
        assert close(b.area, 2 * (dx * dy + dy * dz + dx * dz))
        assert np.allclose(b.centroid, lo + size / 2, rtol=1e-12, atol=1e-12)
        assert close(b.flip().volume, -b.volume)
        assert close(b.translate([7e5, 7.5e6, 100]).volume, b.volume, rel=1e-8)


def test_face_measures(box):
    b = box([0, 0, 0], [1, 1, 1])
    assert b.face_areas.sum() == pytest.approx(6.0)
    n = b.face_normals
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)
    assert np.allclose(((b.face_centroids - 0.5) * n).sum(axis=1), 0.5)
    assert np.allclose(np.linalg.norm(b.vertex_normals, axis=1), 1.0)
    dip = b.dip
    assert set(np.round(dip[:, 0]).tolist()) == {0.0, 90.0}
    assert np.isnan(dip[dip[:, 0] == 0, 1]).all()


def test_volume_by_elevation_sums_to_volume(orebody):
    lo, hi = orebody.bounds
    edges = np.linspace(lo[2] - 1, hi[2] + 1, 13)
    slabs = orebody.volume_by_elevation(edges)
    assert len(slabs) == 12 and (slabs >= -1e-6).all()
    assert slabs.sum() == pytest.approx(orebody.volume, rel=1e-9)


def test_transform_flips_with_negative_determinant(box):
    b = box([0, 0, 0], [1, 2, 3])
    mirror = np.diag([-1.0, 1.0, 1.0, 1.0])
    m = b.transform(mirror)
    assert m.volume == pytest.approx(6.0) and m.is_solid
    scale = np.diag([2.0, 2.0, 2.0, 1.0])
    assert b.transform(scale).volume == pytest.approx(48.0)
    with pytest.raises(wf.InvalidInput):
        b.transform(np.eye(3))


def test_to_length_unit(box):
    b = box([0, 0, 0], [1, 1, 1])
    ft = b.to_length_unit("ft")
    assert ft.length_unit == "ft" and ft.volume == pytest.approx(1 / 0.3048**3)
    with pytest.raises(wf.InvalidInput):
        wf.Mesh(b.coords, b.triangles).to_length_unit("ft")


def test_repair_restores_a_shuffled_duplicated_flipped_solid(orebody, rng):
    soup_v = orebody.coords[orebody.triangles].reshape(-1, 3)
    tris = np.arange(len(soup_v)).reshape(-1, 3)
    flip = rng.random(len(tris)) < 0.3
    tris[flip] = tris[flip][:, ::-1]
    tris = np.vstack([tris, tris[:5]])
    soup = wf.Mesh(soup_v, tris[rng.permutation(len(tris))])
    assert not soup.is_solid
    fixed, report = soup.repair()
    assert fixed.is_solid
    assert fixed.volume == pytest.approx(orebody.volume, rel=1e-9)
    s = report.summary
    assert s["vertex_merged"] > 0 and s["face_removed"] == 5 and s["face_flipped"] > 0
    assert report.before["is_solid"] is False and report.after["is_solid"] is True
    assert len(report) == len(report.actions)
    vmap, fmap = report.vertex_map, report.face_map
    assert len(vmap) == soup.n_vertices and len(fmap) == soup.n_faces
    kept = fmap >= 0
    assert np.array_equal(
        np.sort(vmap[soup.triangles[kept]], axis=1), np.sort(fixed.triangles[fmap[kept]], axis=1)
    )
    again, second = fixed.repair()
    assert len(second) == 0 and np.array_equal(again.triangles, fixed.triangles)


def test_repair_fills_holes(box):
    open_box = box([0, 0, 0], [1, 1, 1]).filter_faces(np.arange(12) > 1)
    for fill in ("fan", "min_area"):
        fixed, report = open_box.repair(fill_holes=fill)
        assert fixed.is_solid and fixed.volume == pytest.approx(1.0)
        assert report.summary["face_added"] == 2 and report.holes_left == 0
    _, report = open_box.repair()
    assert report.holes_left == 1
    with pytest.raises(wf.InvalidInput):
        open_box.repair(fill_holes="nope")


def test_components_add_up(box):
    m = wf.concat([box([0, 0, 0], [1, 1, 1]), box([3, 0, 0], [5, 1, 1])])
    parts = m.components()
    assert [p.n_faces for p in parts] == [12, 12]
    assert sum(p.volume for p in parts) == pytest.approx(m.volume)
    assert m.component_labels().tolist() == [0] * 12 + [1] * 12


def test_thickness_of_two_planes():
    g = np.linspace(0, 10, 6)
    gx, gy = np.meshgrid(g, g)
    xy = np.column_stack([gx.ravel(), gy.ravel()])
    top = wf.dtm(np.column_stack([xy, 5 + 0.1 * xy[:, 0]]))
    bottom = wf.dtm(np.column_stack([xy, 1 + 0.05 * xy[:, 1]]))
    q = np.array([[2.5, 3.5], [7.0, 1.0], [20.0, 20.0]])
    t = top.thickness(bottom, q)
    expected = 5 + 0.1 * q[:, 0] - 1 - 0.05 * q[:, 1]
    assert np.allclose(t[:2], expected[:2]) and np.isnan(t[2])
    assert np.allclose(top.surface_z(q)[:2], 5 + 0.1 * q[:2, 0])
    p3 = np.column_stack([q[:2], [10.0, 0.0]])
    assert np.allclose(top.vertical_distance(p3), p3[:, 2] - (5 + 0.1 * q[:2, 0]))
    normal = top.thickness(bottom, q[:2], mode="normal")
    assert (normal <= t[:2] + 1e-12).all() and (normal > 0).all()
    with pytest.raises(wf.InvalidInput):
        top.thickness(bottom, q, mode="sideways")


plt = pytest.importorskip("matplotlib.pyplot")


@pytest.mark.parametrize("view", ["plan", "section", "surface3d"])
def test_plot_returns_fig_ax_and_reuses_ax(view, orebody):
    holes = wf.datasets.drillholes(n=3)
    points = holes.to_points()
    fn = getattr(wf.plot, view)
    args = (wf.datasets.ORIGIN, [1, 0, 0]) if view == "section" else ()
    fig, ax = fn([orebody, holes, points], *args)
    assert ax.figure is fig
    projection = "3d" if view == "surface3d" else None
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(projection=projection)
    got = fn(orebody, *args, ax=ax2)
    assert got == (fig2, ax2)
    plt.close("all")


def test_plot_values_and_problems(orebody):
    _, ax = wf.plot.plan(orebody, "section", edges=True, label=True)
    assert ax.get_xlabel() == "x (m)"
    _, ax = wf.plot.section(orebody, wf.datasets.ORIGIN, [0, 1, 0], values="section", thickness=10)
    assert len(ax.lines) > 0
    _, ax = wf.plot.section(orebody, wf.datasets.ORIGIN + [0, 0, 200], [0.2, 0.1, 1])
    assert ax.get_xlabel() == "u (m)"
    holed = orebody.filter_faces(np.arange(orebody.n_faces) != 3)
    _, ax = wf.plot.surface3d(holed, np.arange(holed.n_faces), problems=holed.validate())
    assert len(ax.collections) == 2
    with pytest.raises(wf.InvalidInput):
        wf.plot.plan(orebody, np.ones(5))
    with pytest.raises(wf.InvalidInput):
        wf.plot.plan("not geometry")
    plt.close("all")


def test_repair_keeps_touching_shells_apart(box):
    a = box([0, 0, 0], [1, 1, 1])
    two = a | box([1, 1, 0], [2, 2, 1])
    assert two.is_solid and two.n_vertices == 16
    fixed, report = two.repair()
    assert len(report) == 0 and fixed.n_vertices == 16
    welded, report = two.repair(weld_touching_shells=True)
    assert report.summary["vertex_merged"] == 2 and welded.n_vertices == 14
    assert welded.is_solid
    split, report = welded.repair()
    assert report.summary["vertex_split"] == 2 and split.n_vertices == 16
    actions = report.actions
    assert set(actions["action"].tolist()) == {"vertex_split"}
    assert split.validate().summary["touching_shell_edges"] == 1
    assert len(split.repair()[1]) == 0


def test_cavities_touching_along_an_edge_are_solids_like_touching_bodies(box):
    holes = [box([1, 1, 1], [2, 2, 2]), box([2, 2, 1], [3, 3, 2])]
    bodies = wf.concat(holes)
    cavities = wf.concat([box([0, 0, 0], [4, 4, 4])] + [h.flip() for h in holes])
    for m in (bodies, cavities):
        s = m.validate().summary
        assert m.is_solid and s["touching_shell_edges"] == 1 and s["self_intersections"] == 0
    assert close(cavities.volume, 62.0)
    assert close((box([0, 0, 0], [4, 4, 4]) - bodies).volume, 62.0)


def test_plot_colours_text_face_and_point_columns_by_category():
    from wireframes import plot

    m = wf.datasets.box([0, 0, 0], [1, 1, 1]).with_face_column("layer", ["ORE", "WASTE", None] * 4)
    fig, ax = plot.plan(m, "layer")
    bar = fig.axes[-1]
    assert [t.get_text() for t in bar.get_yticklabels()] == ["ORE", "WASTE"]
    plt.close(fig)
    fig, _ = plot.surface3d(m, "layer")
    plt.close(fig)
    points = wf.PointSet([[0, 0, 0], [1, 0, 0], [2, 0, 0]], {"hole": ["DH2", "DH1", "DH2"]})
    fig, ax = plot.plan(points, "hole")
    assert ax.collections[0].get_array().tolist() == [0.0, 1.0, 0.0]
    plt.close(fig)
    rock = m.with_vertex_column("rock", ["a", "b"] * 4)
    with pytest.raises(wf.InvalidInput, match="'rock' cannot shade triangles"):
        plot.plan(rock, "rock")


def test_plot_missing_column_names_vertex_and_face_columns():
    from wireframes import plot

    m = (
        wf.datasets.box([0, 0, 0], [1, 1, 1])
        .with_face_column("rock", [1.0] * 12)
        .with_vertex_column("au", [0.5] * 8)
    )
    for view in (plot.plan, plot.surface3d):
        with pytest.raises(wf.MissingColumn, match=r'"gold"; columns are \["au", "rock"\]'):
            view(m, "gold")
