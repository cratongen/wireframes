"""Static matplotlib views: plan, section and 3D surface.

Every function takes ``ax=None``, creates a figure when it is omitted, and returns ``(fig, ax)``.
Needs matplotlib: ``pip install "drusa[plot]"``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib import colormaps
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
except ImportError as e:
    raise ImportError('drusa.plot needs matplotlib: pip install "drusa[plot]"') from e

from ._drusa import Mesh, MeshReport, PointSet, Polylines
from .errors import InvalidInput

__all__ = ["plan", "section", "surface3d"]


def _items(data: Any) -> list[Mesh | Polylines | PointSet]:
    items = list(data) if isinstance(data, Sequence) and not isinstance(data, str) else [data]
    for item in items:
        if not isinstance(item, (Mesh, Polylines, PointSet)):
            raise InvalidInput(f"plot data must be Mesh, Polylines or PointSet, not {type(item).__name__}")
    return items


def _figure(ax: Any, projection: str | None = None) -> tuple[Any, Any]:
    if ax is not None:
        return ax.figure, ax
    fig = plt.figure(figsize=(6.4, 5.2), layout="constrained")
    return fig, fig.add_subplot(projection=projection)


def _unit_label(items: list[Any], axis: str) -> str:
    unit = next((i.length_unit for i in items if i.length_unit), None)
    return f"{axis} ({unit})" if unit else axis


def _mesh_values(mesh: Mesh, values: Any) -> tuple[np.ndarray, bool]:
    """``(values, per_face)``; vertex z when `values` is None."""
    if values is None:
        return mesh.z, False
    if isinstance(values, str):
        if values in mesh.face_attributes:
            return np.asarray(mesh.face_attributes[values], dtype=float), True
        return np.asarray(mesh.vertex_attributes[values], dtype=float), False
    v = np.asarray(values, dtype=float)
    if len(v) == mesh.n_faces:
        return v, True
    if len(v) == mesh.n_vertices:
        return v, False
    raise InvalidInput(
        f"values must have one value per face ({mesh.n_faces}) or per vertex ({mesh.n_vertices})"
    )


def _face_values(mesh: Mesh, values: Any) -> np.ndarray:
    v, per_face = _mesh_values(mesh, values)
    return v if per_face else v[mesh.triangles].mean(axis=1)


def _point_values(points: PointSet, values: Any) -> np.ndarray | None:
    if values is None:
        return None
    if isinstance(values, str):
        return np.asarray(points[values], dtype=float) if values in points.attributes else None
    v = np.asarray(values, dtype=float)
    return v if len(v) == len(points) else None


def _ring(part: np.ndarray, closed: bool) -> np.ndarray:
    return np.vstack([part, part[:1]]) if closed else part


def plan(
    data: Mesh | Polylines | PointSet | Sequence[Mesh | Polylines | PointSet],
    values: str | Any | None = None,
    *,
    ax: Any = None,
    cmap: str = "viridis",
    edges: bool = False,
    label: bool = False,
) -> tuple[Any, Any]:
    """Plan view of meshes, polylines and points.

    Parameters
    ----------
    data : Mesh, Polylines, PointSet or sequence of them
        What to draw, bottom first.
    values : str or array_like, optional
        Colour of meshes and points: a face or vertex column name, or one value per face or vertex. None
        colours meshes by elevation.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; None creates a figure.
    cmap : str
        Colormap name.
    edges : bool
        Draw triangle edges.
    label : bool
        Write each polyline feature's index at its first vertex.

    Returns
    -------
    tuple of (Figure, Axes)
    """
    items = _items(data)
    fig, ax = _figure(ax)
    for item in items:
        if isinstance(item, Mesh):
            if item.n_faces == 0:
                continue
            v, per_face = _mesh_values(item, values)
            style = {"cmap": cmap, "edgecolors": "0.3" if edges else "none", "linewidth": 0.2}
            if per_face:
                art = ax.tripcolor(item.x, item.y, item.triangles, facecolors=v, **style)
            else:
                art = ax.tripcolor(item.x, item.y, item.triangles, v, shading="gouraud", **style)
            if edges and not per_face:
                ax.triplot(item.x, item.y, item.triangles, color="0.3", linewidth=0.2)
            fig.colorbar(art, ax=ax, shrink=0.8, label=values if isinstance(values, str) else "z")
        elif isinstance(item, Polylines):
            features = item.feature
            for p, (part, closed) in enumerate(zip(item.parts, item.closed, strict=True)):
                line = _ring(part, closed)
                ax.plot(line[:, 0], line[:, 1], color="k", linewidth=1.0)
                first = p == 0 or features[p] != features[p - 1]
                if label and first:
                    ax.annotate(str(features[p]), (part[0, 0], part[0, 1]), fontsize=7)
        else:
            c = _point_values(item, values)
            if c is None:
                ax.scatter(item.x, item.y, color="k", s=8)
            else:
                ax.scatter(item.x, item.y, c=c, s=8, cmap=cmap)
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_xlabel(_unit_label(items, "x"))
    ax.set_ylabel(_unit_label(items, "y"))
    return fig, ax


def _basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    """``(u, v, vertical)``: in-plane axes; for a vertical plane ``u`` is horizontal and ``v`` is +z."""
    n = normal / np.linalg.norm(normal)
    if abs(n[2]) < 1e-12:
        u = np.array([-n[1], n[0], 0.0])
        return u / np.linalg.norm(u), np.array([0.0, 0.0, 1.0]), True
    axis = np.eye(3)[np.argmin(np.abs(n))]
    u = np.cross(n, axis)
    u /= np.linalg.norm(u)
    return u, np.cross(n, u), False


def section(
    data: Mesh | Polylines | PointSet | Sequence[Mesh | Polylines | PointSet],
    origin: Any,
    normal: Any,
    *,
    thickness: float | None = None,
    values: str | None = None,
    ax: Any = None,
) -> tuple[Any, Any]:
    """Section view along the plane through `origin` with `normal`.

    Meshes are drawn as their ``slice_by_plane`` lines; polylines and points within ``thickness / 2`` of
    the plane are projected onto it. A vertical plane has the horizontal distance from `origin` and the
    elevation as axes; any other plane uses an in-plane basis.

    Parameters
    ----------
    data : Mesh, Polylines, PointSet or sequence of them
        What to draw.
    origin, normal : array_like
        The section plane.
    thickness : float, optional
        Width of the slab of polylines and points shown; None shows all of them.
    values : str, optional
        Feature column of the section lines (a face column of the mesh) that picks each line's colour.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; None creates a figure.

    Returns
    -------
    tuple of (Figure, Axes)
    """
    items = _items(data)
    o = np.asarray(origin, dtype=float)
    n = np.asarray(normal, dtype=float)
    if o.shape != (3,) or n.shape != (3,) or not np.linalg.norm(n) > 0:
        raise InvalidInput("origin and normal must have 3 values and normal must be nonzero")
    u, v, vertical = _basis(n)
    unit_n = n / np.linalg.norm(n)
    fig, ax = _figure(ax)

    def project(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        d = p - o
        return d @ u, (p[:, 2] if vertical else d @ v)

    def near(p: np.ndarray) -> np.ndarray:
        if thickness is None:
            return np.ones(len(p), dtype=bool)
        return np.abs((p - o) @ unit_n) <= thickness / 2

    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for item in items:
        if isinstance(item, Mesh):
            lines = item.slice_by_plane(o, n)
            colors = ["k"] * len(lines)
            if values is not None and len(lines):
                labels = lines.attributes[values]
                keys = list(dict.fromkeys(labels.tolist()))
                colors = [cycle[keys.index(x) % len(cycle)] for x in labels.tolist()]
            for p, (part, closed) in enumerate(zip(lines.parts, lines.closed, strict=True)):
                s, z = project(_ring(part, closed))
                ax.plot(s, z, color=colors[lines.feature[p]], linewidth=1.2)
        elif isinstance(item, Polylines):
            for part, closed in zip(item.parts, item.closed, strict=True):
                line = _ring(part, closed)
                s, z = project(line)
                keep = near(line)
                ax.plot(np.where(keep, s, np.nan), np.where(keep, z, np.nan), color="0.3", linewidth=0.8)
        else:
            keep = near(item.coords)
            s, z = project(item.coords[keep])
            ax.scatter(s, z, s=8, color="0.2")
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_xlabel(_unit_label(items, "distance along section" if vertical else "u"))
    ax.set_ylabel(_unit_label(items, "z" if vertical else "v"))
    return fig, ax


def surface3d(
    data: Mesh | Polylines | PointSet | Sequence[Mesh | Polylines | PointSet],
    values: str | Any | None = None,
    *,
    ax: Any = None,
    cmap: str = "viridis",
    problems: MeshReport | None = None,
    azim: float = -60.0,
    elev: float = 30.0,
) -> tuple[Any, Any]:
    """3D view of meshes as shaded triangles, with polylines and points.

    Parameters
    ----------
    data : Mesh, Polylines, PointSet or sequence of them
        What to draw.
    values : str or array_like, optional
        Colour of the meshes: a face or vertex column, or one value per face or vertex; None colours by
        elevation.
    ax : mpl_toolkits.mplot3d.Axes3D, optional
        3D axes to draw on; None creates a figure.
    cmap : str
        Colormap name.
    problems : MeshReport, optional
        Validation of the first mesh: boundary and non-manifold edges are drawn in red, other problem
        faces in orange.
    azim, elev : float
        View angles in degrees.

    Returns
    -------
    tuple of (Figure, Axes3D)
    """
    items = _items(data)
    fig, ax = _figure(ax, projection="3d")
    lo = np.full(3, np.inf)
    hi = np.full(3, -np.inf)
    first_mesh = True
    for item in items:
        if item.bounds is not None:
            lo = np.minimum(lo, item.bounds[0])
            hi = np.maximum(hi, item.bounds[1])
        if isinstance(item, Mesh):
            if item.n_faces == 0:
                continue
            tris = item.coords[item.triangles]
            fv = _face_values(item, values)
            norm = Normalize(np.nanmin(fv), np.nanmax(fv))
            faces = colormaps[cmap](norm(fv))
            shade = 0.6 + 0.4 * np.abs(item.face_normals[:, 2])
            faces[:, :3] *= shade[:, None]
            if problems is not None and first_mesh:
                _mark_problems(ax, item, problems, faces)
            ax.add_collection3d(Poly3DCollection(tris, facecolors=faces, edgecolors="none"))
            first_mesh = False
        elif isinstance(item, Polylines):
            for part, closed in zip(item.parts, item.closed, strict=True):
                line = _ring(part, closed)
                ax.plot(line[:, 0], line[:, 1], line[:, 2], color="k", linewidth=1.0)
        else:
            ax.scatter(item.x, item.y, item.z, s=6, color="k")
    if np.all(np.isfinite(lo)):
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(np.maximum(hi - lo, 1e-9 * max(np.max(hi - lo), 1.0)))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel(_unit_label(items, "x"))
    ax.set_ylabel(_unit_label(items, "y"))
    ax.set_zlabel(_unit_label(items, "z"))
    return fig, ax


def _mark_problems(ax: Any, mesh: Mesh, report: MeshReport, faces: np.ndarray) -> None:
    table = report.problems
    kind, face, edge = table["kind"], table["face"], table["edge"]
    has_face, has_edge = table.mask("face"), table.mask("edge")
    segments = []
    for k, f, e, hf, he in zip(kind, face, edge, has_face, has_edge, strict=True):
        if not hf or f >= mesh.n_faces:
            continue
        if k in ("boundary_edge", "non_manifold_edge") and he:
            t = mesh.triangles[f]
            segments.append(mesh.coords[[t[e], t[(e + 1) % 3]]])
        else:
            faces[f] = (1.0, 0.55, 0.0, 1.0)
    if segments:
        ax.add_collection3d(Line3DCollection(segments, colors="red", linewidths=2.0))
