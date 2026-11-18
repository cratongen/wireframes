"""Static matplotlib plots: plan, section and 3D surface views."""

from __future__ import annotations

try:
    import matplotlib  # noqa: F401
except ImportError as e:
    raise ImportError('drusa.plot needs matplotlib: pip install "drusa[plot]"') from e


def plan(data, values=None, *, ax=None, cmap="viridis", edges=False, label=False):
    """Plan view of meshes, polylines and points; returns ``(fig, ax)``."""
    raise NotImplementedError("plan is not implemented yet")


def section(data, origin, normal, *, thickness=None, values=None, ax=None):
    """Section view along the plane through `origin` with `normal`; returns ``(fig, ax)``."""
    raise NotImplementedError("section is not implemented yet")


def surface3d(data, values=None, *, ax=None, cmap="viridis", problems=None, azim=-60.0, elev=30.0):
    """3D view of meshes, polylines and points; returns ``(fig, ax)``."""
    raise NotImplementedError("surface3d is not implemented yet")
