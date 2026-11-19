from collections.abc import Mapping, Sequence
from os import PathLike
from typing import Any, Literal, TypeAlias, overload

import numpy as np
import numpy.typing as npt

ArrayLike = npt.ArrayLike
Path: TypeAlias = str | PathLike[str]
Points: TypeAlias = ArrayLike | PointSet | Mesh
Values: TypeAlias = ArrayLike | Sequence[str | None]
Floats: TypeAlias = npt.NDArray[np.float64]
Ints: TypeAlias = npt.NDArray[np.int64]
Bools: TypeAlias = npt.NDArray[np.bool_]
Bounds: TypeAlias = tuple[list[float], list[float]] | None

__version__: str

class Table:
    """Columnar attribute table: float64, int64, bool and string columns with validity masks and units."""

    def __init__(self, data: Mapping[str, Values] | None = None) -> None: ...
    @property
    def num_rows(self) -> int:
        """Number of rows."""
    @property
    def column_names(self) -> list[str]:
        """Column names in order."""
    @property
    def units(self) -> dict[str, str]:
        """Units of the columns that have one."""
    def column(self, name: str) -> np.ndarray:
        """The values of column `name`: float64 (NaN null), int64, bool or object (None null)."""
    def __getitem__(self, name: str) -> np.ndarray: ...
    def __contains__(self, name: str) -> bool: ...
    def mask(self, name: str) -> Bools:
        """Validity of column `name`: True where the row holds a value."""
    def filter(self, mask: ArrayLike) -> Table:
        """A copy with the rows where `mask` is True."""
    def drop_null(self, *columns: str) -> Table:
        """A copy without the rows where any of `columns` (default: every column) is null."""
    def with_column(self, name: str, values: Values, *, unit: str | None = None) -> Table:
        """A copy with column `name` added, or replaced in place."""
    def with_columns(self, data: Mapping[str, Values]) -> Table:
        """A copy with every column of `data` added or replaced."""
    def to_dict(self) -> dict[str, np.ndarray]:
        """Column name to numpy values."""
    def to_pyarrow(self) -> Any:
        """The table as a pyarrow Table; needs pyarrow."""
    def to_polars(self) -> Any:
        """The table as a polars DataFrame; needs polars."""
    def to_pandas(self) -> Any:
        """The table as a pandas DataFrame; needs pandas."""
    def __len__(self) -> int: ...

class MeshReport:
    """What ``Mesh.validate`` found."""

    @property
    def summary(self) -> dict[str, int | bool]:
        """Problem counts by kind, plus ``shells``, the contact counts ``touching_shell_edges`` and
        ``touching_shell_vertices`` (information, not defects), ``is_closed`` and ``is_solid``."""
    @property
    def problems(self) -> Table:
        """One row per problem: ``kind``, ``face``, ``vertex``, ``edge``, ``other``."""
    @property
    def is_solid(self) -> bool:
        """Closed, manifold, outward and without self-intersections; shells may touch along edges or at vertices."""
    def __len__(self) -> int: ...

class RepairReport:
    """What ``Mesh.repair`` changed."""

    @property
    def actions(self) -> Table:
        """One row per change: ``action``, ``vertex``, ``face``, ``other``, ``reason``."""
    @property
    def vertex_map(self) -> Ints:
        """Output vertex of each input vertex, -1 when removed."""
    @property
    def face_map(self) -> Ints:
        """Output face of each input face, -1 when removed."""
    @property
    def before(self) -> dict[str, int | bool]:
        """Validation summary of the input."""
    @property
    def after(self) -> dict[str, int | bool]:
        """Validation summary of the output."""
    @property
    def holes_left(self) -> int:
        """Boundary loops left open."""
    @property
    def summary(self) -> dict[str, int]:
        """Number of actions of each kind."""
    def __len__(self) -> int: ...

class Mesh:
    """Triangle mesh: an open surface or a closed solid, with vertex and face attribute tables."""

    def __init__(
        self,
        vertices: ArrayLike,
        triangles: ArrayLike,
        *,
        vertex_attributes: Mapping[str, Values] | Table | None = None,
        face_attributes: Mapping[str, Values] | Table | None = None,
        crs: str | None = None,
        length_unit: str | None = None,
    ) -> None: ...
    @property
    def coords(self) -> Floats:
        """``(n, 3)`` vertex coordinates."""
    @property
    def x(self) -> Floats:
        """``(n,)`` x coordinates."""
    @property
    def y(self) -> Floats:
        """``(n,)`` y coordinates."""
    @property
    def z(self) -> Floats:
        """``(n,)`` z coordinates."""
    @property
    def triangles(self) -> Ints:
        """``(m, 3)`` vertex indices of each triangle."""
    @property
    def n_vertices(self) -> int:
        """Number of vertices."""
    @property
    def n_faces(self) -> int:
        """Number of triangles."""
    @property
    def bounds(self) -> Bounds:
        """``(min, max)`` corners, None without vertices."""
    @property
    def column_names(self) -> list[str]:
        """Names of the vertex columns, then the face columns."""
    @property
    def vertex_attributes(self) -> Table:
        """One row per vertex."""
    @property
    def face_attributes(self) -> Table:
        """One row per triangle."""
    @property
    def crs(self) -> str | None:
        """Coordinate reference system label."""
    @property
    def length_unit(self) -> str | None:
        """Length unit label."""
    @property
    def is_closed(self) -> bool:
        """Whether every edge of the non-degenerate faces joins exactly two faces, or more that pair into wedges of solid."""
    @property
    def is_solid(self) -> bool:
        """Whether ``validate()`` finds a solid."""
    @property
    def area(self) -> float:
        """Total triangle area."""
    @property
    def volume(self) -> float:
        """Signed enclosed volume; raises InvalidMesh unless closed and consistently wound; overlaps count twice."""
    @property
    def centroid(self) -> Floats:
        """Volume centroid, under the same rule as ``volume``."""
    @property
    def face_areas(self) -> Floats:
        """Area of each triangle."""
    @property
    def face_normals(self) -> Floats:
        """``(m, 3)`` unit normals, zero for degenerate triangles."""
    @property
    def face_centroids(self) -> Floats:
        """``(m, 3)`` triangle centroids."""
    @property
    def vertex_normals(self) -> Floats:
        """``(n, 3)`` angle-weighted unit normals."""
    @property
    def dip(self) -> Floats:
        """``(m, 2)`` dip and dip direction in degrees."""
    def __len__(self) -> int: ...
    def __getitem__(self, name: str) -> np.ndarray: ...
    def with_vertex_column(self, name: str, values: Values, *, unit: str | None = None) -> Mesh:
        """A copy with vertex column `name` added or replaced."""
    def with_face_column(self, name: str, values: Values, *, unit: str | None = None) -> Mesh:
        """A copy with face column `name` added or replaced."""
    def with_attributes(
        self,
        *,
        vertex: Table | Mapping[str, Values] | None = None,
        face: Table | Mapping[str, Values] | None = None,
    ) -> Mesh:
        """A copy with new attribute tables; None keeps a table."""
    def with_vertices(self, vertices: ArrayLike) -> Mesh:
        """A copy with the vertices moved; triangles and attributes kept."""
    def filter_faces(self, mask: ArrayLike) -> Mesh:
        """A copy with the triangles where `mask` is True and the vertices they use."""
    def to_table(self, *, on: Literal["vertices", "faces"] = "vertices") -> Table:
        """Vertices (x, y, z and vertex columns) or faces (v0, v1, v2 and face columns)."""
    def to_polars(self, *, on: Literal["vertices", "faces"] = "vertices") -> Any:
        """``to_table(on=on)`` as a polars DataFrame."""
    def to_pandas(self, *, on: Literal["vertices", "faces"] = "vertices") -> Any:
        """``to_table(on=on)`` as a pandas DataFrame."""
    def to_pyarrow(self, *, on: Literal["vertices", "faces"] = "vertices") -> Any:
        """``to_table(on=on)`` as a pyarrow Table."""
    def transform(self, matrix: ArrayLike) -> Mesh:
        """A copy with the ``(4, 4)`` affine `matrix` applied."""
    def translate(self, offset: ArrayLike) -> Mesh:
        """A copy moved by `offset`."""
    def flip(self) -> Mesh:
        """A copy with every triangle wound the other way."""
    def to_length_unit(self, unit: str) -> Mesh:
        """A copy rescaled to `unit`."""
    def validate(self, *, tolerance: float = 0.0) -> MeshReport:
        """Every topological and geometric problem, decided exactly; touching closed shells are not problems."""
    def repair(
        self,
        *,
        weld: float | None = 0.0,
        remove_degenerate: bool = True,
        remove_duplicate_faces: bool = True,
        remove_unused_vertices: bool = True,
        orient: bool = True,
        fill_holes: Literal["fan", "min_area"] | None = None,
        weld_touching_shells: bool = False,
    ) -> tuple[Mesh, RepairReport]:
        """A repaired copy and the report of every change; touching closed shells keep their own vertices.

        Self-intersections are never repaired and added faces have null face columns.
        """
    def components(self) -> list[Mesh]:
        """One mesh per edge-connected component; touching shells are separate components."""
    def component_labels(self) -> Ints:
        """Edge-connected component of each triangle."""
    def volume_by_elevation(self, edges: ArrayLike) -> Floats:
        """Solid volume between consecutive elevations."""
    def surface_z(self, points: Points) -> Floats:
        """Highest surface elevation above each point, NaN where none."""
    def vertical_distance(self, points: Points) -> Floats:
        """Point elevation minus ``surface_z``."""
    def thickness(
        self, bottom: Mesh, points: Points, *, mode: Literal["vertical", "normal"] = "vertical"
    ) -> Floats:
        """Thickness between this surface and `bottom` at each point."""
    @overload
    def union(self, other: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
    @overload
    def union(
        self, other: Mesh, *, repair: Literal[True]
    ) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
    @overload
    def intersection(self, other: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
    @overload
    def intersection(
        self, other: Mesh, *, repair: Literal[True]
    ) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
    @overload
    def difference(self, other: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
    @overload
    def difference(
        self, other: Mesh, *, repair: Literal[True]
    ) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
    @overload
    def symmetric_difference(self, other: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
    @overload
    def symmetric_difference(
        self, other: Mesh, *, repair: Literal[True]
    ) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
    def __or__(self, other: Mesh) -> Mesh: ...
    def __and__(self, other: Mesh) -> Mesh: ...
    def __sub__(self, other: Mesh) -> Mesh: ...
    def __xor__(self, other: Mesh) -> Mesh: ...
    def clip_by_surface(self, surface: Mesh, *, repair: bool = False) -> tuple[Mesh, Mesh]:
        """``(below, above)`` parts of the solid cut by an open surface, both capped."""
    def split_by_plane(self, origin: ArrayLike, normal: ArrayLike) -> tuple[Mesh, Mesh]:
        """``(negative, positive)`` parts split by a plane."""
    def slice_by_plane(
        self, origin: ArrayLike, normal: ArrayLike, *, offsets: ArrayLike | None = None
    ) -> Polylines:
        """Section lines at each offset along the normal."""
    def locate(self, points: Points) -> npt.NDArray[np.int8]:
        """-1 inside, 0 on the boundary, 1 outside."""
    def contains(self, points: Points, *, boundary: bool = True) -> Bools:
        """Whether each point lies in the solid."""
    def signed_distance(self, points: Points) -> Floats:
        """Distance to the surface, negative inside."""
    def distance(self, points: Points) -> Floats:
        """Unsigned distance to the closest face."""
    def closest_point(self, points: Points) -> tuple[Floats, Ints, Floats]:
        """Closest points, their faces and distances."""
    def ray_hits(self, origins: ArrayLike, directions: ArrayLike) -> Table:
        """Hits of each ray: ray, t, face, x, y, z, kind, entering."""
    def segment_hits(self, starts: ArrayLike, ends: ArrayLike) -> Table:
        """Hits of each segment: segment, t, face, x, y, z, kind, entering."""
    def trace_intervals(self, trace: ArrayLike, depths: ArrayLike) -> Table:
        """Inside and outside runs of a trace: from, to, inside, entry_face, exit_face."""
    def block_fraction(
        self, centroids: Points, size: ArrayLike, *, rotation: ArrayLike | None = None
    ) -> Floats:
        """Exact fraction of each block inside the solid."""
    def close_to_elevation(self, z: float) -> Mesh:
        """The solid between this surface and the plane at `z`."""
    def decimate(
        self,
        *,
        target_faces: int | None = None,
        max_error: float | None = None,
        preserve_boundary: bool = True,
    ) -> Mesh:
        """A copy with fewer triangles."""
    def smooth(self, *, iterations: int = 10, lambda_: float = 0.5, mu: float = -0.53) -> Mesh:
        """A Taubin-smoothed copy."""
    def offset(self, distance: float) -> Mesh:
        """A copy with each vertex moved by `distance` along its normal; faces at corners move less."""
    def remesh(self, edge_length: float, *, iterations: int = 5) -> Mesh:
        """An isotropically remeshed copy."""
    def contours(self, levels: ArrayLike) -> Polylines:
        """Contour lines at each level."""

class Polylines:
    """3D strings: open lines and closed rings grouped into features."""

    def __init__(
        self,
        parts: Sequence[ArrayLike],
        *,
        closed: bool | Sequence[bool] = False,
        features: Sequence[int] | None = None,
        attributes: Mapping[str, Values] | Table | None = None,
        vertex_attributes: Mapping[str, Values] | Table | None = None,
        crs: str | None = None,
        length_unit: str | None = None,
    ) -> None: ...
    @property
    def coords(self) -> Floats:
        """``(n, 3)`` coordinates of every vertex."""
    @property
    def x(self) -> Floats:
        """``(n,)`` x coordinates."""
    @property
    def y(self) -> Floats:
        """``(n,)`` y coordinates."""
    @property
    def z(self) -> Floats:
        """``(n,)`` z coordinates."""
    @property
    def parts(self) -> list[Floats]:
        """``(k, 3)`` vertices of each part."""
    @property
    def closed(self) -> Bools:
        """Whether each part is a ring."""
    @property
    def feature(self) -> Ints:
        """Feature index of each part."""
    @property
    def attributes(self) -> Table:
        """One row per feature."""
    @property
    def vertex_attributes(self) -> Table:
        """One row per vertex."""
    @property
    def bounds(self) -> Bounds:
        """``(min, max)`` corners."""
    @property
    def column_names(self) -> list[str]:
        """Names of the feature columns, then the vertex columns."""
    @property
    def crs(self) -> str | None:
        """Coordinate reference system label."""
    @property
    def length_unit(self) -> str | None:
        """Length unit label."""
    def __len__(self) -> int: ...
    def __getitem__(self, name: str) -> np.ndarray:
        """The feature or vertex column `name`."""
    def filter(self, mask: ArrayLike) -> Polylines:
        """A copy with the features where `mask` is True."""
    def with_column(self, name: str, values: Values, *, unit: str | None = None) -> Polylines:
        """A copy with feature column `name` added or replaced."""
    def with_vertex_column(self, name: str, values: Values, *, unit: str | None = None) -> Polylines:
        """A copy with vertex column `name` added or replaced."""
    def length(self) -> Floats:
        """3D length of each feature."""
    def area(self) -> Floats:
        """Plan area of each feature, even-odd."""
    def contains(self, points: Points, *, feature: int | None = None) -> Bools:
        """Whether each point lies inside in plan."""
    def locate(self, points: Points) -> Ints:
        """First feature holding each point in plan, -1 for none."""
    def distance(self, points: Points, *, feature: int | None = None, signed: bool = False) -> Floats:
        """Plan distance to the nearest ring edge."""
    def to_table(self) -> Table:
        """One row per vertex with feature, part, vertex, x, y, z, closed and the columns."""
    def to_polars(self) -> Any:
        """``to_table()`` as a polars DataFrame."""
    def to_pandas(self) -> Any:
        """``to_table()`` as a pandas DataFrame."""
    def to_pyarrow(self) -> Any:
        """``to_table()`` as a pyarrow Table."""
    @staticmethod
    def from_table(
        table: Table | Mapping[str, Values],
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
        """Polylines from a long table of vertices."""
    def to_points(self) -> PointSet:
        """Every vertex as a point with feature and part columns."""

class PointSet:
    """3D point cloud with an attribute table."""

    def __init__(
        self,
        coords: ArrayLike,
        attributes: Mapping[str, Values] | Table | None = None,
        *,
        crs: str | None = None,
        length_unit: str | None = None,
    ) -> None: ...
    @property
    def coords(self) -> Floats:
        """``(n, 3)`` coordinates."""
    @property
    def x(self) -> Floats:
        """``(n,)`` x coordinates."""
    @property
    def y(self) -> Floats:
        """``(n,)`` y coordinates."""
    @property
    def z(self) -> Floats:
        """``(n,)`` z coordinates."""
    @property
    def attributes(self) -> Table:
        """One row per point."""
    @property
    def bounds(self) -> Bounds:
        """``(min, max)`` corners."""
    @property
    def column_names(self) -> list[str]:
        """Names of the attribute columns."""
    @property
    def crs(self) -> str | None:
        """Coordinate reference system label."""
    @property
    def length_unit(self) -> str | None:
        """Length unit label."""
    def __len__(self) -> int: ...
    def __getitem__(self, name: str) -> np.ndarray: ...
    def filter(self, mask: ArrayLike) -> PointSet:
        """A copy with the points where `mask` is True."""
    def with_column(self, name: str, values: Values, *, unit: str | None = None) -> PointSet:
        """A copy with column `name` added or replaced."""
    def to_table(self) -> Table:
        """x, y, z and the attribute columns."""
    def to_polars(self) -> Any:
        """``to_table()`` as a polars DataFrame."""
    def to_pandas(self) -> Any:
        """``to_table()`` as a pandas DataFrame."""
    def to_pyarrow(self) -> Any:
        """``to_table()`` as a pyarrow Table."""

class UnstructuredGrid:
    """Volume and mixed cells kept as data for round trips."""

    @property
    def points(self) -> Floats:
        """``(n, 3)`` point coordinates."""
    @property
    def cell_types(self) -> npt.NDArray[np.uint8]:
        """VTK cell type id of each cell."""
    @property
    def offsets(self) -> Ints:
        """``n_cells + 1`` offsets into ``connectivity``."""
    @property
    def connectivity(self) -> Ints:
        """Point ids of every cell."""
    @property
    def point_data(self) -> Table:
        """One row per point."""
    @property
    def cell_data(self) -> Table:
        """One row per cell."""
    @property
    def n_cells(self) -> int:
        """Number of cells."""
    def surface(self) -> Mesh:
        """Triangles of the surface cells."""

def _read(
    path: Path, *, crs: str | None = None, length_unit: str | None = None
) -> tuple[
    Mesh | None, Polylines | None, PointSet | None, UnstructuredGrid | None, str | None, list[str]
]: ...
def read_mesh(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> Mesh:
    """The mesh in the file at `path`."""

def read_polylines(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> Polylines:
    """The polylines in the file at `path`."""

def read_points(path: Path, *, crs: str | None = None, length_unit: str | None = None) -> PointSet:
    """The points in the file at `path`."""

def write(
    path: Path,
    data: Mesh | Polylines | PointSet | UnstructuredGrid | Any,
    *,
    ascii: bool = False,
    dxf_entity: Literal["3dface", "polyface"] = "3dface",
    zlib: bool = False,
    skip_unsupported_columns: bool = False,
) -> None:
    """Writes `data` to `path`; the extension picks the format."""

def concat(meshes: Sequence[Mesh]) -> Mesh:
    """Meshes appended in order, no welding."""

@overload
def union(a: Mesh, b: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
@overload
def union(a: Mesh, b: Mesh, *, repair: Literal[True]) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
@overload
def intersection(a: Mesh, b: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
@overload
def intersection(
    a: Mesh, b: Mesh, *, repair: Literal[True]
) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
@overload
def difference(a: Mesh, b: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
@overload
def difference(
    a: Mesh, b: Mesh, *, repair: Literal[True]
) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
@overload
def symmetric_difference(a: Mesh, b: Mesh, *, repair: Literal[False] = False) -> Mesh: ...
@overload
def symmetric_difference(
    a: Mesh, b: Mesh, *, repair: Literal[True]
) -> tuple[Mesh, tuple[RepairReport, RepairReport]]: ...
def cut_fill(base: Mesh, design: Mesh, *, boundary: Polylines | None = None) -> dict[str, float]:
    """Cut, fill, net and areas of `design` against `base`."""

def dtm(
    points: Points,
    *,
    breaklines: Polylines | None = None,
    boundary: Polylines | None = None,
    min_spacing: float = 0.0,
    max_edge: float | None = None,
) -> Mesh:
    """Surface by constrained Delaunay triangulation in plan."""

def triangulate(polygons: Polylines) -> Mesh:
    """Each feature's rings triangulated, holes by even-odd."""

def loft(sections: Polylines, *, cap_ends: bool = True) -> Mesh:
    """Solid linking consecutive section rings; the end sections may be a line or a point."""

def extrude(polygons: Polylines, vector: ArrayLike) -> Mesh:
    """Solid swept by each ring along `vector`."""
