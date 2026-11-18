"""Pure-Python containers and wrappers around the compiled module."""

from __future__ import annotations

from os import PathLike
from typing import NamedTuple

from . import _drusa


class Dataset(NamedTuple):
    """Everything one file holds."""

    mesh: _drusa.Mesh | None
    polylines: _drusa.Polylines | None
    points: _drusa.PointSet | None
    grid: _drusa.UnstructuredGrid | None
    length_unit: str | None
    warnings: list[str]


def read(path: str | PathLike[str], *, crs: str | None = None, length_unit: str | None = None) -> Dataset:
    """Everything in the file at `path`, labelled with `crs` and `length_unit`."""
    return Dataset(*_drusa._read(path, crs=crs, length_unit=length_unit))
