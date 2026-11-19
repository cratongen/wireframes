"""Pure-Python accessors and conversions shared by the compiled containers."""

from __future__ import annotations

import math
from collections.abc import Mapping
from itertools import pairwise
from os import PathLike
from typing import Any, NamedTuple

import numpy as np

from . import _wireframes
from ._wireframes import Mesh, PointSet, Polylines, Table
from .errors import InvalidInput


class Dataset(NamedTuple):
    """Everything one file holds.

    Attributes
    ----------
    mesh : Mesh or None
        Triangles of the file, polygons fanned from their first corner.
    polylines : Polylines or None
        Lines and rings of the file.
    points : PointSet or None
        Isolated points of the file.
    grid : UnstructuredGrid or None
        Volume cells of a ``.vtu`` or ``.vtk`` file.
    length_unit : str or None
        Unit the file declares (DXF ``$INSUNITS``), never applied as a rescale.
    warnings : list of str
        What the reader skipped, such as ``"skipped 3 INSERT entities"``.
    """

    mesh: Mesh | None
    polylines: Polylines | None
    points: PointSet | None
    grid: _wireframes.UnstructuredGrid | None
    length_unit: str | None
    warnings: list[str]


def read(path: str | PathLike[str], *, crs: str | None = None, length_unit: str | None = None) -> Dataset:
    """Everything in the file at `path`.

    Parameters
    ----------
    path : str or PathLike
        File to read; the extension picks the format (``.obj``, ``.ply``, ``.off``, ``.stl``, ``.dxf``,
        ``.vtk``, ``.vtu``, ``.vtp``, or a Datamine ``_tr.csv``/``_pt.csv`` pair named by either file).
    crs : str, optional
        Coordinate reference system label set on every returned object.
    length_unit : str, optional
        Length unit set on every returned object; None takes the unit the file declares.

    Returns
    -------
    Dataset
        Mesh, polylines, points, grid, the file's length unit and reader warnings.

    Raises
    ------
    FileError
        When the file cannot be read or its content is malformed.
    InvalidInput
        When `length_unit` disagrees with the unit the file declares.
    """
    return Dataset(*_wireframes._read(path, crs=crs, length_unit=length_unit))


def _is_null(x: Any) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def _column(values: Any) -> tuple[str, Any, Any]:
    """``(kind, data, valid)`` of a column given as an array or sequence; ``valid`` is None when no row is null."""
    a = np.asarray(values)
    if a.ndim != 1:
        raise InvalidInput(f"column values must be one-dimensional, got shape {a.shape}")
    if a.dtype.kind == "b":
        return "bool", np.ascontiguousarray(a, dtype=bool), None
    if a.dtype.kind in "iu":
        return "int64", np.ascontiguousarray(a, dtype=np.int64), None
    if a.dtype.kind == "f":
        data = np.ascontiguousarray(a, dtype=np.float64)
        valid = ~np.isnan(data)
        return "float64", data, None if valid.all() else valid
    if a.dtype.kind in "US":
        return "utf8", [str(x) for x in a.tolist()], None
    items = a.tolist()
    valid = np.array([not _is_null(x) for x in items], dtype=bool)
    present = [x for x, ok in zip(items, valid, strict=True) if ok]
    mask = None if valid.all() else valid
    if all(isinstance(x, str) for x in present) and present:
        return "utf8", ["" if not ok else x for x, ok in zip(items, valid, strict=True)], mask
    if all(isinstance(x, (bool, np.bool_)) for x in present) and present:
        return "bool", np.array([bool(x) if ok else False for x, ok in zip(items, valid, strict=True)]), mask
    if all(isinstance(x, (int, np.integer)) and not isinstance(x, bool) for x in present) and present:
        data = np.array([int(x) if ok else 0 for x, ok in zip(items, valid, strict=True)], dtype=np.int64)
        return "int64", data, mask
    try:
        data = np.array([float(x) if ok else np.nan for x, ok in zip(items, valid, strict=True)])
    except (TypeError, ValueError):
        raise InvalidInput("column values must be all numbers, all booleans or all strings") from None
    return "float64", data, mask


def _require(module: str, method: str) -> Any:
    import importlib

    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(f"{method} needs {module}: pip install {module}") from e


def _to_pyarrow(self: Table) -> Any:
    pa = _require("pyarrow", "to_pyarrow")
    arrays = {}
    for name in self.column_names:
        values, valid = self.column(name), self.mask(name)
        mask = None if valid.all() else ~valid
        if values.dtype == object:
            arrays[name] = pa.array(
                [v if ok else None for v, ok in zip(values, valid, strict=True)], pa.string()
            )
        else:
            arrays[name] = pa.array(values, mask=mask)
    return pa.table(arrays) if arrays else pa.table({}).slice(0, 0)


def _to_polars(self: Table) -> Any:
    pl = _require("polars", "to_polars")
    series = []
    for name in self.column_names:
        values, valid = self.column(name), self.mask(name)
        if valid.all() and values.dtype != object:
            series.append(pl.Series(name, values))
        else:
            dtype = {"f": pl.Float64, "i": pl.Int64, "b": pl.Boolean, "O": pl.String}[values.dtype.kind]
            items = [v if ok else None for v, ok in zip(values.tolist(), valid, strict=True)]
            series.append(pl.Series(name, items, dtype=dtype))
    return pl.DataFrame(series)


def _to_pandas(self: Table) -> Any:
    pd = _require("pandas", "to_pandas")
    columns = {}
    for name in self.column_names:
        values, valid = self.column(name), self.mask(name)
        if valid.all() or values.dtype.kind in "fO":
            columns[name] = values
        elif values.dtype.kind == "i":
            columns[name] = pd.arrays.IntegerArray(values, ~valid)
        else:
            columns[name] = pd.arrays.BooleanArray(values, ~valid)
    return pd.DataFrame(columns, index=pd.RangeIndex(self.num_rows))


def _drop_null(self: Table, *columns: str) -> Table:
    """A copy without the rows where any of `columns` is null; no names means every column.

    Parameters
    ----------
    *columns : str
        Columns to check.

    Returns
    -------
    Table
        The rows where every checked column is valid.
    """
    keep = np.ones(self.num_rows, dtype=bool)
    for name in columns or self.column_names:
        keep &= self.mask(name)
    return self.filter(keep)


def _from_table(
    table: Table | Mapping[str, Any],
    *,
    feature: str = "ID",
    x: str = "X",
    y: str = "Y",
    z: str | None = None,
    part: str | None = None,
    closed: bool = False,
    crs: str | None = None,
    length_unit: str | None = None,
) -> Polylines:
    """Polylines from a long table of vertices, one row per vertex.

    Consecutive rows with the same `feature` value form one feature, and within it consecutive rows with
    the same `part` value form one part. The other columns become vertex columns, and the `feature`
    column a feature column.

    Parameters
    ----------
    table : Table or mapping of str to array_like
        Vertex rows in drawing order.
    feature, x, y : str
        Columns holding the feature id and the plan coordinates.
    z : str, optional
        Column holding elevation; None sets z to 0.
    part : str, optional
        Column holding the part id within a feature; None makes each feature one part.
    closed : bool
        Whether every part is a ring; a ring does not repeat its first vertex.
    crs, length_unit : str, optional
        Labels of the result.

    Returns
    -------
    Polylines
        One feature per run of equal `feature` values.
    """
    if not isinstance(table, Table):
        table = Table(dict(table))
    ids = table[feature]
    xs, ys = table[x].astype(float), table[y].astype(float)
    zs = table[z].astype(float) if z is not None else np.zeros(len(xs))
    parts_of = table[part] if part is not None else np.zeros(len(xs), dtype=np.int64)
    coords = np.column_stack([xs, ys, zs])
    breaks = [0]
    features = []
    first_rows = []
    f = -1
    for i in range(len(xs)):
        new_feature = i == 0 or ids[i] != ids[i - 1]
        if new_feature or parts_of[i] != parts_of[i - 1]:
            if i > 0:
                breaks.append(i)
            if new_feature:
                f += 1
                first_rows.append(i)
            features.append(f)
    breaks.append(len(xs))
    pieces = [coords[a:b] for a, b in pairwise(breaks)]
    if not len(xs):
        pieces, features = [], []
    skip = {feature, x, y, z, part}
    vertex = {n: table[n] for n in table.column_names if n not in skip}
    return Polylines(
        pieces,
        closed=closed,
        features=features,
        attributes={feature: ids[first_rows]},
        vertex_attributes=vertex,
        crs=crs,
        length_unit=length_unit,
    )


def _to_points(self: Polylines) -> PointSet:
    """Every vertex as a point, with ``feature`` and ``part`` columns, the vertex columns and the feature columns.

    Returns
    -------
    PointSet
        One point per vertex, in vertex order.
    """
    table = self.to_table()
    attributes = {n: table[n] for n in table.column_names if n not in ("x", "y", "z", "vertex", "closed")}
    return PointSet(self.coords, attributes, crs=self.crs, length_unit=self.length_unit)


def _axis(a: int) -> property:
    names = "xyz"
    return property(lambda self: self.coords[:, a], doc=f"``(n,)`` {names[a]} coordinates.")


for _cls in (Mesh, Polylines, PointSet):
    _cls.x, _cls.y, _cls.z = (_axis(a) for a in range(3))

for _cls in (Mesh, Polylines, PointSet):
    _cls.to_polars = lambda self, **kw: self.to_table(**kw).to_polars()
    _cls.to_pandas = lambda self, **kw: self.to_table(**kw).to_pandas()
    _cls.to_pyarrow = lambda self, **kw: self.to_table(**kw).to_pyarrow()
    _cls.to_polars.__doc__ = "``to_table()`` as a polars DataFrame; needs polars."
    _cls.to_pandas.__doc__ = "``to_table()`` as a pandas DataFrame; needs pandas."
    _cls.to_pyarrow.__doc__ = "``to_table()`` as a pyarrow Table; needs pyarrow."

Table.to_pyarrow = _to_pyarrow
Table.to_polars = _to_polars
Table.to_pandas = _to_pandas
Table.to_pyarrow.__doc__ = "The table as a pyarrow Table, nulls kept; needs pyarrow."
Table.to_polars.__doc__ = "The table as a polars DataFrame, nulls kept; needs polars."
Table.to_pandas.__doc__ = "The table as a pandas DataFrame; int and bool nulls become nullable dtypes."
Table.drop_null = _drop_null
Polylines.from_table = staticmethod(_from_table)
Polylines.to_points = _to_points
