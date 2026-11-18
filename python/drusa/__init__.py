"""Read, write, plot, measure and combine solids, surfaces and strings."""

from . import errors
from ._containers import Dataset, read
from ._drusa import (
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
from .errors import DrusaError, FileError, InvalidInput, InvalidMesh, MissingColumn

__all__ = [
    "Dataset",
    "DrusaError",
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
    "__version__",
    "concat",
    "cut_fill",
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
    if name == "plot":
        import importlib

        return importlib.import_module(f"{__name__}.plot")
    raise AttributeError(f"module 'drusa' has no attribute {name!r}")
