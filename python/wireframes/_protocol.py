"""Read and write the ``__minetypes__`` protocol, version 1.

The producer side lives in the extension. This module checks dicts by hand and
turns them into wireframes objects, block fractions and drillhole runs.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from ._wireframes import Mesh, PointSet, Polylines, Table
from .errors import InvalidInput

_OWN = (Mesh, PointSet, Polylines, Table)
_COLUMN_DTYPES = (
    np.dtype(np.float64),
    np.dtype(np.float32),
    np.dtype(np.int64),
    np.dtype(np.int32),
    np.dtype(np.bool_),
)
_KINDS = (
    "table",
    "categories",
    "grid",
    "block_model",
    "point_set",
    "mesh",
    "polylines",
    "drillholes",
    "intervals",
    "composites",
    "seismic_volume",
    "prestack_volume",
)

_array_block_fraction = None
_array_trace_intervals = None


def _bad(message: str) -> InvalidInput:
    return InvalidInput(message)


def _native(array: np.ndarray) -> bool:
    order = array.dtype.byteorder
    if order in "=|":
        return True
    native = "<" if np.little_endian else ">"
    return order == native


def _is_int(value: Any) -> bool:
    return isinstance(value, (int, np.integer)) and not isinstance(value, bool)


def _require(d: Mapping[str, Any], key: str) -> Any:
    if key not in d:
        raise _bad(f"key {key!r} is missing")
    return d[key]


def _version(d: Mapping[str, Any]) -> None:
    value = _require(d, "protocol")
    if not _is_int(value) or int(value) != 1:
        raise _bad(f"key 'protocol' must be 1, not {value!r}")


def _label(d: Mapping[str, Any], key: str) -> str | None:
    value = _require(d, key)
    if value is None:
        return None
    if not isinstance(value, str) or value.strip() == "":
        raise _bad(f"key {key!r} is blank; use None when it is undeclared")
    return value


def _as_dict(obj: Any) -> dict[str, Any]:
    method = getattr(obj, "__minetypes__", None)
    if callable(method) and not isinstance(obj, dict):
        spec = method()
    elif isinstance(obj, dict):
        spec = obj
    else:
        raise _bad("expected a protocol dict or an object with __minetypes__")
    if not isinstance(spec, dict):
        raise _bad("__minetypes__ must return a dict")
    return spec


def _foreign(obj: Any) -> dict[str, Any] | None:
    """Protocol dict of something that is not already a wireframes container."""
    if isinstance(obj, _OWN):
        return None
    method = getattr(obj, "__minetypes__", None)
    if callable(method):
        spec = method()
        if not isinstance(spec, dict):
            raise _bad("__minetypes__ must return a dict")
        return spec
    if isinstance(obj, dict) and "protocol" in obj and "kind" in obj:
        return obj
    return None


def _header(d: Mapping[str, Any], kind: str) -> None:
    _version(d)
    got = _require(d, "kind")
    if got != kind:
        if got not in _KINDS:
            raise _bad(f"key 'kind' is {got!r}, which is not a protocol kind")
        raise _bad(f"key 'kind' must be {kind!r}, not {got!r}")
    _label(d, "crs")
    _label(d, "length_unit")


def _check_array(
    array: Any, key: str, dtype: np.dtype | None, ndim: int, cols: int | None, length: int | None
) -> np.ndarray:
    if not isinstance(array, np.ndarray):
        raise _bad(f"key {key!r} must be a numpy array")
    if dtype is not None and array.dtype != dtype:
        raise _bad(f"key {key!r} must be dtype {np.dtype(dtype).name}, not {array.dtype}")
    if array.ndim != ndim:
        raise _bad(f"key {key!r} must have {ndim} dimensions, not shape {array.shape}")
    if cols is not None and (array.shape[-1] != cols):
        raise _bad(f"key {key!r} must have shape (..., {cols}), not {array.shape}")
    if length is not None and array.shape[0] != length:
        raise _bad(f"key {key!r} must have length {length}, not {array.shape[0]}")
    if not array.flags.c_contiguous:
        raise _bad(f"key {key!r} must be C-contiguous")
    if not _native(array):
        raise _bad(f"key {key!r} must be native-endian")
    return array


def _finite(array: np.ndarray, key: str) -> None:
    if not np.isfinite(array).all():
        raise _bad(f"key {key!r} must be finite")


def _dict(d: Mapping[str, Any], key: str) -> dict[str, Any]:
    value = _require(d, key)
    if not isinstance(value, dict):
        raise _bad(f"key {key!r} must be a dict")
    return value


def _category(name: str, entry: Any) -> list[str]:
    if not isinstance(entry, dict):
        raise _bad(f"category {name!r} must be a dict")
    for key in ("names", "colors", "mapping", "other"):
        if key not in entry:
            raise _bad(f"category {name!r} is missing key {key!r}")
    names = list(entry["names"])
    if (
        not names
        or len(set(names)) != len(names)
        or any(not isinstance(item, str) or item == "" for item in names)
    ):
        raise _bad(f"category {name!r} names must be distinct and non-empty")
    colors = entry["colors"]
    if colors is not None:
        colors = list(colors)
        if len(colors) != len(names) or any(not isinstance(item, str) or item == "" for item in colors):
            raise _bad(f"category {name!r} colors must be non-empty strings, one per name")
    mapping = entry["mapping"]
    if not isinstance(mapping, dict):
        raise _bad(f"category {name!r} mapping must be a dict")
    for key, value in mapping.items():
        if not isinstance(key, str) or not isinstance(value, str) or key in names or value not in names:
            raise _bad(f"category {name!r} mapping must send an alias to a name")
    other = entry["other"]
    if other is not None and other != names[-1]:
        raise _bad(f"category {name!r} other must equal the last name")
    return names


def _column_values(name: str, array: np.ndarray, entry: Any) -> Any:
    """Values ``Table`` can store. Category columns become text, or bool for false/true."""
    if entry is not None:
        names = _category(name, entry)
        _check_array(array, name, np.dtype(np.int32), 1, None, None)
        if np.any((array < -1) | (array >= len(names))):
            raise _bad(f"key {name!r} has a code outside -1 .. {len(names) - 1}")
        if (
            names == ["false", "true"]
            and entry["colors"] is None
            and entry["mapping"] == {}
            and entry["other"] is None
        ):
            return [None if int(code) < 0 else bool(int(code)) for code in array]
        return [None if int(code) < 0 else names[int(code)] for code in array]
    _check_array(array, name, None, 1, None, None)
    if array.dtype.kind in "UT":
        return [str(item) for item in array.tolist()]
    if array.dtype not in _COLUMN_DTYPES:
        raise _bad(f"key {name!r} has dtype {array.dtype}, which a column does not use")
    if array.dtype.kind == "f":
        if np.isinf(array).any():
            raise _bad(f"key {name!r} must not contain infinity")
        if array.dtype == np.float32:
            return np.array(array, dtype=np.float64)
        return array
    if array.dtype == np.int32:
        return np.array(array, dtype=np.int64)
    return array


def _columns_table(columns: Any, units: Any, categories: Any, rows: int) -> Table:
    if not isinstance(columns, dict) or not isinstance(units, dict) or not isinstance(categories, dict):
        raise _bad("columns, units and categories must be dicts")
    for name in categories:
        if name not in columns:
            raise _bad(f"key 'categories' names {name!r}, which is not a column")
    built: dict[str, Any] = {}
    for name, array in columns.items():
        if not isinstance(name, str) or name == "":
            raise _bad("a column name is empty")
        if not isinstance(array, np.ndarray) or array.shape != (rows,):
            shape = getattr(array, "shape", None)
            raise _bad(f"key {name!r} must have shape ({rows},), not {shape}")
        built[name] = _column_values(name, array, categories.get(name))
    for name, unit in units.items():
        if name not in built:
            raise _bad(f"key 'units' names {name!r}, which is not a column")
        if not isinstance(unit, str) or unit == "":
            raise _bad(f"column {name!r} has a blank unit")
    table = Table(built) if built else Table._with_rows(rows)
    if table.num_rows != rows:
        raise _bad(f"key 'rows' is {rows}, but the columns have length {table.num_rows}")
    for name, unit in units.items():
        table = table.with_column(name, built[name], unit=unit)
    return table


def _table_from(obj: Any) -> Table:
    d = _as_dict(obj)
    _header(d, "table")
    if d["crs"] is not None or d["length_unit"] is not None:
        raise _bad("a wireframes Table stores no crs or length_unit; both must be None")
    rows = _require(d, "rows")
    if not _is_int(rows) or int(rows) < 0:
        raise _bad(f"key 'rows' must be an integer >= 0, not {rows!r}")
    return _columns_table(_require(d, "columns"), _require(d, "units"), _require(d, "categories"), int(rows))


def _mesh_from(cls: type[Mesh], obj: Any) -> Mesh:
    d = _as_dict(obj)
    _header(d, "mesh")
    for key in ("columns", "units", "categories"):
        if _dict(d, key):
            raise _bad(f"key {key!r} on a mesh must be empty")
    coords = _check_array(_require(d, "coords"), "coords", np.dtype(np.float64), 2, 3, None)
    _finite(coords, "coords")
    if coords.shape[0] > 2**32 - 1:
        raise _bad(f"key 'coords' has {coords.shape[0]} rows; a mesh allows at most {2**32 - 1}")
    triangles = _check_array(_require(d, "triangles"), "triangles", np.dtype(np.uint32), 2, 3, None)
    if triangles.size and int(triangles.max()) >= coords.shape[0]:
        raise _bad("key 'triangles' indexes a vertex that coords does not have")
    vertex = _columns_table(
        _require(d, "vertex_columns"),
        _require(d, "vertex_units"),
        _require(d, "vertex_categories"),
        coords.shape[0],
    )
    face = _columns_table(
        _require(d, "face_columns"),
        _require(d, "face_units"),
        _require(d, "face_categories"),
        triangles.shape[0],
    )
    return cls(
        coords,
        triangles,
        vertex_attributes=vertex,
        face_attributes=face,
        crs=_label(d, "crs"),
        length_unit=_label(d, "length_unit"),
    )


def _offsets(d: Mapping[str, Any], key: str, total: int) -> np.ndarray:
    offsets = _check_array(_require(d, key), key, np.dtype(np.uint32), 1, None, None)
    if offsets.size == 0 or int(offsets[0]) != 0 or int(offsets[-1]) != total:
        raise _bad(f"key {key!r} must start at 0 and end at {total}")
    if np.any(np.diff(offsets.astype(np.int64)) < 0):
        raise _bad(f"key {key!r} must be non-decreasing")
    return offsets


def _polylines_from(cls: type[Polylines], obj: Any) -> Polylines:
    d = _as_dict(obj)
    _header(d, "polylines")
    coords = _check_array(_require(d, "coords"), "coords", np.dtype(np.float64), 2, 3, None)
    _finite(coords, "coords")
    if coords.shape[0] > 2**31 - 1:
        raise _bad(f"key 'coords' has {coords.shape[0]} rows; polylines allow at most {2**31 - 1}")
    parts = _offsets(d, "part_offsets", coords.shape[0])
    n_parts = len(parts) - 1
    features = _offsets(d, "feature_offsets", n_parts)
    closed = _check_array(_require(d, "closed"), "closed", np.dtype(np.bool_), 1, None, n_parts)
    spans = np.diff(parts.astype(np.int64))
    for p, (count, ring) in enumerate(zip(spans, closed, strict=True)):
        need = 3 if ring else 2
        if int(count) < need:
            noun = "vertex" if int(count) == 1 else "vertices"
            which = "closed" if ring else "open"
            raise _bad(
                f"part {p} has {int(count)} {noun}; {which} parts need at least {need}. "
                "A single point belongs in a PointSet"
            )
    feature_of_part: list[int] = []
    for f, (a, b) in enumerate(itertools.pairwise(features)):
        feature_of_part.extend([f] * (int(b) - int(a)))
    pieces = [coords[int(a) : int(b)] for a, b in itertools.pairwise(parts)]
    attributes = _columns_table(
        _require(d, "columns"), _require(d, "units"), _require(d, "categories"), len(features) - 1
    )
    vertex = _columns_table(
        _require(d, "vertex_columns"),
        _require(d, "vertex_units"),
        _require(d, "vertex_categories"),
        coords.shape[0],
    )
    return cls(
        pieces,
        closed=closed.tolist(),
        features=feature_of_part,
        attributes=attributes,
        vertex_attributes=vertex,
        crs=_label(d, "crs"),
        length_unit=_label(d, "length_unit"),
    )


def _points_from(cls: type[PointSet], obj: Any) -> PointSet:
    d = _as_dict(obj)
    _header(d, "point_set")
    coords = _check_array(_require(d, "coords"), "coords", np.dtype(np.float64), 2, 3, None)
    _finite(coords, "coords")
    attributes = _columns_table(
        _require(d, "columns"), _require(d, "units"), _require(d, "categories"), coords.shape[0]
    )
    return cls(coords, attributes, crs=_label(d, "crs"), length_unit=_label(d, "length_unit"))


def _triple(value: Any, key: str, *, positive: bool) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (3,) or not np.isfinite(array).all():
        raise _bad(f"key {key!r} must be 3 finite floats")
    if positive and np.any(array <= 0):
        raise _bad(f"key {key!r} must be greater than 0")
    return array


def _count(value: Any) -> tuple[int, int, int]:
    if isinstance(value, np.ndarray):
        items = value.tolist()
    else:
        items = list(value)
    if len(items) != 3 or any(not _is_int(item) or int(item) < 1 for item in items):
        raise _bad(f"key 'count' must be 3 integers >= 1, not {value!r}")
    nx, ny, nz = (int(item) for item in items)
    if nx * ny * nz > 2**63 - 1:
        raise _bad("key 'count' has more cells than a signed 64-bit integer holds")
    return nx, ny, nz


def _frame(rotation: np.ndarray) -> np.ndarray:
    azimuth, dip, rake = (math.radians(float(x)) for x in rotation)
    sa, ca = math.sin(azimuth), math.cos(azimuth)
    sd, cd = math.sin(dip), math.cos(dip)
    sr, cr = math.sin(rake), math.cos(rake)
    r1 = np.array([[sa, ca, 0.0], [-ca, sa, 0.0], [0.0, 0.0, 1.0]])
    r2 = np.array([[cd, 0.0, -sd], [0.0, 1.0, 0.0], [sd, 0.0, cd]])
    r3 = np.array([[1.0, 0.0, 0.0], [0.0, cr, sr], [0.0, -sr, cr]])
    aniso = r3 @ r2 @ r1
    return np.stack((-aniso[1], aniso[0], aniso[2]))


def _centroids(
    origin: np.ndarray,
    size: np.ndarray,
    count: tuple[int, int, int],
    index: np.ndarray,
    centre: np.ndarray,
    frame: np.ndarray,
) -> np.ndarray:
    nx, ny, _nz = count
    plane = np.uint64(nx * ny)
    ix = (index % np.uint64(nx)).astype(np.float64)
    iy = ((index // np.uint64(nx)) % np.uint64(ny)).astype(np.float64)
    iz = (index // plane).astype(np.float64)
    local = np.column_stack(
        (
            (ix + centre[:, 0]) * size[0],
            (iy + centre[:, 1]) * size[1],
            (iz + centre[:, 2]) * size[2],
        )
    )
    return origin + local @ frame


def _check_subgrid(subgrid: Any, offset: np.ndarray, scale: np.ndarray) -> None:
    if subgrid is None:
        return
    items = subgrid.tolist() if isinstance(subgrid, np.ndarray) else list(subgrid)
    if len(items) != 3 or any(not _is_int(item) or not 1 <= int(item) <= 2**32 - 1 for item in items):
        raise _bad("key 'subgrid' must be 3 integers in 1 .. 2**32 - 1, or None")
    for axis, item in enumerate(items):
        m = np.float32(int(item))
        for edge in (
            (offset[:, axis] - np.float32(0.5) * scale[:, axis]) * m,
            (offset[:, axis] + np.float32(0.5) * scale[:, axis]) * m,
        ):
            nearest = np.rint(edge)
            if np.any(np.abs(edge - nearest) > np.float32(1e-5)) or np.any((nearest < 0) | (nearest > m)):
                raise _bad("key 'subgrid' does not match offset and scale")


def _blocks(d: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    kind = d["kind"]
    _header(d, kind)
    origin = _triple(_require(d, "origin"), "origin", positive=False)
    size = _triple(_require(d, "size"), "size", positive=True)
    count = _count(_require(d, "count"))
    rotation = _triple(_require(d, "rotation"), "rotation", positive=False)
    n_parent = count[0] * count[1] * count[2]
    if kind == "grid":
        index = np.arange(n_parent, dtype=np.uint64)
        centre = np.full((n_parent, 3), 0.5, dtype=np.float64)
        scale = np.ones((n_parent, 3), dtype=np.float64)
    else:
        layout = _require(d, "layout")
        if layout not in ("regular", "masked", "subblocked"):
            raise _bad(f"key 'layout' must be 'regular', 'masked' or 'subblocked', not {layout!r}")
        if layout == "regular":
            for key in ("index", "parent", "offset", "scale", "subgrid"):
                if _require(d, key) is not None:
                    raise _bad(f"key {key!r} must be None on a regular block model")
            index = np.arange(n_parent, dtype=np.uint64)
            centre = np.full((n_parent, 3), 0.5, dtype=np.float64)
            scale = np.ones((n_parent, 3), dtype=np.float64)
        elif layout == "masked":
            for key in ("parent", "offset", "scale", "subgrid"):
                if _require(d, key) is not None:
                    raise _bad(f"key {key!r} must be None on a masked block model")
            index = _check_array(_require(d, "index"), "index", np.dtype(np.uint64), 1, None, None)
            if index.size and (np.any(index[1:] <= index[:-1]) or np.any(index >= np.uint64(n_parent))):
                raise _bad("key 'index' must be strictly increasing and address a parent cell")
            centre = np.full((len(index), 3), 0.5, dtype=np.float64)
            scale = np.ones((len(index), 3), dtype=np.float64)
        else:
            if _require(d, "index") is not None:
                raise _bad("key 'index' must be None on a sub-blocked model")
            index = _check_array(_require(d, "parent"), "parent", np.dtype(np.uint64), 1, None, None)
            if index.size and (np.any(index[1:] < index[:-1]) or np.any(index >= np.uint64(n_parent))):
                raise _bad("key 'parent' must be non-decreasing and address a parent cell")
            offset = _check_array(_require(d, "offset"), "offset", np.dtype(np.float32), 2, 3, len(index))
            raw_scale = _check_array(_require(d, "scale"), "scale", np.dtype(np.float32), 2, 3, len(index))
            if not np.isfinite(offset).all() or np.any((offset <= 0) | (offset >= 1)):
                raise _bad("key 'offset' must be finite and inside (0, 1)")
            if not np.isfinite(raw_scale).all() or np.any((raw_scale <= 0) | (raw_scale > 1)):
                raise _bad("key 'scale' must be finite and inside (0, 1]")
            lo = offset - np.float32(0.5) * raw_scale
            hi = offset + np.float32(0.5) * raw_scale
            if np.any(lo < np.float32(-1e-6)) or np.any(hi > np.float32(1 + 1e-6)):
                raise _bad("key 'offset' and 'scale' put a sub-block outside its parent cell")
            _check_subgrid(_require(d, "subgrid"), offset, raw_scale)
            centre = offset.astype(np.float64)
            scale = raw_scale.astype(np.float64)
    columns = _dict(d, "columns")
    for name, array in columns.items():
        _check_array(array, name, None, 1, None, len(index))
        if array.dtype.kind == "f" and np.isinf(array).any():
            raise _bad(f"key {name!r} must not contain infinity")
        if array.dtype.kind not in "fiuUbT" and array.dtype != np.dtype(np.bool_):
            raise _bad(f"key {name!r} has dtype {array.dtype}, which a column does not use")
        if array.dtype.kind in "iu" and array.dtype not in (np.dtype(np.int32), np.dtype(np.int64)):
            raise _bad(f"key {name!r} has dtype {array.dtype}, which a column does not use")
    _dict(d, "units")
    _dict(d, "categories")
    frame = _frame(rotation)
    centroids = _centroids(origin, size, count, index, centre, frame)
    world = scale * size
    axes = np.ascontiguousarray(frame.T)
    if np.array_equal(axes, np.eye(3)):
        axes = None
    return centroids, world, axes


def _fractions(mesh: Mesh, spec: Mapping[str, Any]) -> np.ndarray:
    centroids, world, axes = _blocks(spec)
    n = len(centroids)
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    order = np.arange(n)
    if n > 1:
        keys = np.ascontiguousarray(world).view(np.uint8).reshape(n, -1)
        order = np.lexsort(tuple(keys[:, k] for k in range(keys.shape[1])))
    out = np.empty(n, dtype=np.float64)
    start = 0
    while start < n:
        end = start + 1
        while end < n and np.array_equal(world[order[end]], world[order[start]]):
            end += 1
        chosen = order[start:end]
        out[chosen] = _array_block_fraction(mesh, centroids[chosen], world[chosen[0]], rotation=axes)
        start = end
    return out


def block_fraction(self: Mesh, centroids: Any, size: Any = None, *, rotation: Any = None) -> np.ndarray:
    """Fraction of each block inside the solid.

    Parameters
    ----------
    centroids : array_like, PointSet, Mesh, or protocol object
        ``(n, 3)`` block centres, or one object of kind ``grid`` or
        ``block_model``. A protocol object is the only argument: its
        ``origin``, ``size``, ``count`` and ``rotation`` place the blocks,
        and a sub-blocked model uses ``parent``, ``offset`` and ``scale``.
    size : array_like, optional
        Block size along its three axes, positive. Omit it for a protocol
        object. One size applies to every centre in an array call.
    rotation : array_like, optional
        ``(3, 3)`` matrix whose columns are the block axes in world
        coordinates (orthonormal, determinant +1). Omit it for a protocol
        object. None is axis-aligned.

    Returns
    -------
    numpy.ndarray
        Values in ``[0, 1]``, one per centre or per protocol row. Over blocks
        that cover the solid, fraction times block volume sums to the solid
        volume. A sub-block volume is ``dx * dy * dz * su * sv * sw`` with
        the float32 scale widened before the product.

    Raises
    ------
    InvalidMesh
        When the mesh is not a solid; ``report`` holds the validation.
    InvalidInput
        When ``size`` is not three positive finite values, ``rotation`` is
        not a rotation matrix, or the protocol dict breaks a key, dtype,
        shape, contiguity or endian rule. The message names the key or the
        kind.

    Examples
    --------
    >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    >>> m.block_fraction([[2, 1, 1], [1, 1, 1]], [1, 1, 1]).tolist()
    [0.5, 1.0]
    """
    spec = _foreign(centroids)
    if spec is not None:
        if size is not None or rotation is not None:
            raise _bad(
                "block_fraction reads size and rotation from a grid or block_model; pass that object alone"
            )
        kind = spec.get("kind", None)
        if kind not in ("grid", "block_model"):
            shown = kind if "kind" in spec else None
            raise _bad(f"block_fraction accepts kind 'grid' or 'block_model', not {shown!r}")
        return _fractions(self, spec)
    if size is None:
        raise _bad("block_fraction needs a size, or one grid or block_model protocol object")
    return _array_block_fraction(self, centroids, size, rotation=rotation)


def _sequence_str(value: Any, key: str, *, distinct: bool) -> list[str]:
    if isinstance(value, np.ndarray):
        items = value.tolist()
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        items = list(value)
    else:
        raise _bad(f"key {key!r} must be a sequence of strings")
    if any(not isinstance(item, str) or item == "" for item in items):
        raise _bad(f"key {key!r} must be non-empty strings")
    if distinct and len(set(items)) != len(items):
        raise _bad(f"key {key!r} has a repeated id")
    return items


def _check_survey(survey: Any, holes: int) -> None:
    if survey is None:
        return
    if not isinstance(survey, dict):
        raise _bad("key 'survey' must be a dict or None")
    hole = _check_array(_require(survey, "hole_index"), "hole_index", np.dtype(np.int32), 1, None, None)
    depth = _check_array(_require(survey, "depth"), "depth", np.dtype(np.float64), 1, None, len(hole))
    azimuth = _check_array(_require(survey, "azimuth"), "azimuth", np.dtype(np.float64), 1, None, len(hole))
    dip = _check_array(_require(survey, "dip"), "dip", np.dtype(np.float64), 1, None, len(hole))
    if len(hole) and (np.any((hole < 0) | (hole >= holes)) or np.any(hole[1:] < hole[:-1])):
        raise _bad("key 'hole_index' must be sorted and address a hole")
    _finite(depth, "depth")
    _finite(azimuth, "azimuth")
    _finite(dip, "dip")
    if np.any(depth < 0) or np.any((azimuth < 0) | (azimuth >= 360)) or np.any((dip < -90) | (dip > 90)):
        raise _bad("key 'survey' has a depth, azimuth or dip outside its range")
    for start in range(holes):
        rows = hole == start
        taken = depth[rows]
        if len(taken) > 1 and np.any(np.diff(taken) <= 0):
            raise _bad("key 'depth' must be strictly increasing inside a hole")


def _check_intervals(spec: Any, holes: int, key: str) -> None:
    if spec is None:
        return
    if not isinstance(spec, dict):
        raise _bad(f"key {key!r} must be a dict or None")
    hole = _check_array(_require(spec, "hole_index"), "hole_index", np.dtype(np.int32), 1, None, None)
    start = _check_array(_require(spec, "from_depth"), "from_depth", np.dtype(np.float64), 1, None, len(hole))
    end = _check_array(_require(spec, "to_depth"), "to_depth", np.dtype(np.float64), 1, None, len(hole))
    _dict(spec, "units")
    _dict(spec, "categories")
    columns = _dict(spec, "columns")
    for name, array in columns.items():
        _check_array(array, name, None, 1, None, len(hole))
        if array.dtype.kind == "f" and np.isinf(array).any():
            raise _bad(f"key {name!r} must not contain infinity")
    if len(hole) == 0:
        return
    if np.any((hole < 0) | (hole >= holes)) or np.any(hole[1:] < hole[:-1]):
        raise _bad(f"key {key!r} hole_index must be sorted and address a hole")
    _finite(start, "from_depth")
    _finite(end, "to_depth")
    if np.any(start < 0) or np.any(end <= start):
        raise _bad(f"key {key!r} needs 0 <= from_depth < to_depth")
    for h in range(holes):
        rows = np.flatnonzero(hole == h)
        if len(rows) > 1 and np.any(start[rows[1:]] < end[rows[:-1]] - 1e-9):
            raise _bad(f"key {key!r} has overlapping intervals")


def _drillholes(d: Mapping[str, Any]) -> list[tuple[int, np.ndarray, np.ndarray]]:
    _header(d, "drillholes")
    holes = _sequence_str(_require(d, "holes"), "holes", distinct=True)
    _check_array(_require(d, "coords"), "coords", np.dtype(np.float64), 2, 3, len(holes))
    _finite(d["coords"], "coords")
    for name, array in _dict(d, "columns").items():
        _check_array(array, name, None, 1, None, len(holes))
    _dict(d, "units")
    _dict(d, "categories")
    _check_survey(_require(d, "survey"), len(holes))
    _check_intervals(_require(d, "intervals"), len(holes), "intervals")
    tables = _dict(d, "tables")
    for name, table in tables.items():
        _check_intervals(table, len(holes), name)
    method = _require(d, "method")
    if method not in ("minimum_curvature", "tangential", "balanced_tangential"):
        raise _bad(f"key 'method' must name a desurvey method, not {method!r}")
    path = _check_array(_require(d, "path_coords"), "path_coords", np.dtype(np.float64), 2, 3, None)
    depth = _check_array(_require(d, "path_depths"), "path_depths", np.dtype(np.float64), 1, None, len(path))
    offsets = _offsets(d, "path_offsets", len(path))
    if len(offsets) != len(holes) + 1:
        raise _bad("key 'path_offsets' must have one entry more than holes")
    _finite(path, "path_coords")
    _finite(depth, "path_depths")
    runs = []
    for h, (a, b) in enumerate(itertools.pairwise(offsets)):
        a, b = int(a), int(b)
        if b - a < 2:
            raise _bad(f"hole {h} has {b - a} path vertices; a hole needs at least 2")
        taken = depth[a:b]
        if np.any(np.diff(taken) <= 0):
            raise _bad("key 'path_depths' must be strictly increasing inside a hole")
        runs.append((h, path[a:b], taken))
    return runs


def _runs_table(parts: list[tuple[int, Table]]) -> Table:
    if not parts:
        return Table(
            {
                "from": np.zeros(0),
                "to": np.zeros(0),
                "inside": np.zeros(0, dtype=bool),
                "entry_face": np.zeros(0, dtype=np.int64),
                "exit_face": np.zeros(0, dtype=np.int64),
                "hole_index": np.zeros(0, dtype=np.int64),
            }
        )
    rows: dict[str, list[Any]] = {
        name: [] for name in ("from", "to", "inside", "entry_face", "exit_face", "hole_index")
    }
    for hole, iv in parts:
        rows["from"].extend(iv["from"].tolist())
        rows["to"].extend(iv["to"].tolist())
        rows["inside"].extend(bool(v) for v in iv["inside"].tolist())
        for name in ("entry_face", "exit_face"):
            mask = iv.mask(name)
            rows[name].extend(int(v) if ok else None for v, ok in zip(iv[name], mask, strict=True))
        rows["hole_index"].extend([hole] * len(iv))
    return Table(rows)


def trace_intervals(self: Mesh, trace: Any, depths: Any = None) -> Table:
    """Inside and outside runs of a trace through the solid.

    Parameters
    ----------
    trace : array_like or protocol object
        ``(n, 3)`` desurveyed vertices, or one object of kind ``drillholes``.
        A drillholes object is the only argument. Each hole is the slice of
        ``path_coords`` and ``path_depths`` between ``path_offsets``.
    depths : array_like, optional
        Measured depth of each vertex, non-decreasing. Omit it for a
        drillholes object.

    Returns
    -------
    Table
        ``from``, ``to``, ``inside`` (alternating), ``entry_face`` and
        ``exit_face`` (int64, null where a run starts or ends at an end of
        the trace instead of at a face; null reads as 0, so check
        ``Table.mask``). A drillholes object also returns ``hole_index``,
        the position in ``holes``, stored as int64. The rows of one hole
        stay together, and holes stay in order. The runs of one trace cover
        ``depths[0]`` to ``depths[-1]``.

    Raises
    ------
    InvalidMesh
        When the mesh is not a solid; ``report`` holds the validation.
    InvalidInput
        When the depths decrease, their count differs from the vertex count,
        or the protocol dict breaks a key, dtype, shape, contiguity or endian
        rule. The message names the key or the kind.

    Examples
    --------
    >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    >>> runs = m.trace_intervals([[1, 1, -1], [1, 1, 3]], [0.0, 4.0])
    >>> runs["from"].tolist(), runs["to"].tolist(), runs["inside"].tolist()
    ([0.0, 1.0, 3.0], [1.0, 3.0, 4.0], [False, True, False])
    """
    spec = _foreign(trace)
    if spec is not None:
        if depths is not None:
            raise _bad("trace_intervals reads depths from a drillholes object; pass that object alone")
        kind = spec.get("kind", None)
        if kind != "drillholes":
            raise _bad(f"trace_intervals accepts kind 'drillholes', not {kind!r}")
        runs = _drillholes(spec)
        parts = [(hole, _array_trace_intervals(self, xyz, depth)) for hole, xyz, depth in runs]
        return _runs_table(parts)
    if depths is None:
        raise _bad("trace_intervals needs depths, or one drillholes protocol object")
    return _array_trace_intervals(self, trace, depths)


def _from_table(obj: Any) -> Table:
    return _table_from(obj)


def _from_mesh(obj: Any) -> Mesh:
    return _mesh_from(Mesh, obj)


def _from_polylines(obj: Any) -> Polylines:
    return _polylines_from(Polylines, obj)


def _from_points(obj: Any) -> PointSet:
    return _points_from(PointSet, obj)


_from_table.__name__ = "from_minetypes"
_from_mesh.__name__ = "from_minetypes"
_from_polylines.__name__ = "from_minetypes"
_from_points.__name__ = "from_minetypes"

_FROM_DOC = """\
Build this object from a protocol dict or from an object with ``__minetypes__``.

Parameters
----------
obj : dict or object
    A dict of protocol version 1, or an object whose ``__minetypes__``
    returns one, of this kind.

Returns
-------
{kind}
    Geometry and columns copied into a wireframes object. The input arrays
    are read and not written. A float32 column widens to float64. An int32
    column that is not a category widens to int64. A category column becomes
    the category name as text, with a null where the code is -1. The pair
    ``["false", "true"]`` with no colors and no mapping becomes a bool
    column. Colors and the mapping are not stored.

Raises
------
InvalidInput
    When the kind, a key, a dtype, a shape, contiguity or endianness does
    not match the protocol. The message names the key or the kind. An open
    part with one vertex names ``PointSet``.

Examples
--------
>>> box = wf.datasets.box([0, 0, 0], [1, 1, 1])
>>> wf.Mesh.from_minetypes(box.__minetypes__()).volume
1.0
"""


def install() -> None:
    global _array_block_fraction, _array_trace_intervals
    if _array_block_fraction is not None:
        return
    _array_block_fraction = Mesh.block_fraction
    _array_trace_intervals = Mesh.trace_intervals
    Mesh.block_fraction = block_fraction
    Mesh.trace_intervals = trace_intervals
    for cls, func, kind in (
        (Table, _from_table, "Table"),
        (Mesh, _from_mesh, "Mesh"),
        (Polylines, _from_polylines, "Polylines"),
        (PointSet, _from_points, "PointSet"),
    ):
        func.__doc__ = _FROM_DOC.format(kind=kind)
        func.__qualname__ = f"{kind}.from_minetypes"
        cls.from_minetypes = staticmethod(func)

    def __arrow_c_stream__(self: Table, requested_schema: Any = None) -> Any:
        """Arrow C stream of this table; nulls stay null. Needs pyarrow.

        Examples
        --------
        >>> t = wf.Table({"au": [1.0, None]})
        >>> import pyarrow as pa
        >>> pa.RecordBatchReader.from_stream(t).read_all().column("au").null_count
        1
        """
        return self.to_pyarrow().__arrow_c_stream__(requested_schema)

    Table.__arrow_c_stream__ = __arrow_c_stream__
