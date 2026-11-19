"""Read, write, plot, measure and combine solids, surfaces and strings."""

from . import _containers, errors
from ._containers import Dataset, read
from ._wireframes import (
    Mesh,
    MeshReport,
    PointSet,
    Polylines,
    RepairReport,
    Table,
    UnstructuredGrid,
    __version__,
    concat,
    cut_fill,
    difference,
    dtm,
    extrude,
    intersection,
    loft,
    read_mesh,
    read_points,
    read_polylines,
    symmetric_difference,
    triangulate,
    union,
    write,
)
from .errors import FileError, InvalidInput, InvalidMesh, MissingColumn, WireframesError

del _containers

__all__ = [
    "Dataset",
    "FileError",
    "InvalidInput",
    "InvalidMesh",
    "Mesh",
    "MeshReport",
    "MissingColumn",
    "PointSet",
    "Polylines",
    "RepairReport",
    "Table",
    "UnstructuredGrid",
    "WireframesError",
    "__version__",
    "concat",
    "cut_fill",
    "datasets",
    "difference",
    "dtm",
    "errors",
    "extrude",
    "intersection",
    "loft",
    "plot",
    "read",
    "read_mesh",
    "read_points",
    "read_polylines",
    "symmetric_difference",
    "triangulate",
    "union",
    "write",
]


def __getattr__(name: str):
    if name in ("plot", "datasets"):
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module 'wireframes' has no attribute {name!r}")
