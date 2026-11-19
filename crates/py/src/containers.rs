use std::sync::OnceLock;

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyMapping};
use wireframes_core::{Column, ColumnData, Error, Mesh, PointSet, Polylines, Table};
use wireframes_io::UnstructuredGrid;
use wireframes_query::{SolidIndex, SurfaceIndex};

use crate::errors::{invalid, to_py};

pub(crate) fn numpy(py: Python) -> PyResult<Bound<PyModule>> {
    py.import("numpy")
}

fn asarray<'py>(obj: &Bound<'py, PyAny>, dtype: &str) -> PyResult<Bound<'py, PyAny>> {
    numpy(obj.py())?.call_method1("ascontiguousarray", (obj, dtype))
}

/// `asarray`, with an empty input as `(0, 3)` so `[]` means no rows.
fn asarray2<'py>(obj: &Bound<'py, PyAny>, dtype: &str) -> PyResult<Bound<'py, PyAny>> {
    let a = asarray(obj, dtype)?;
    if a.getattr("size")?.extract::<usize>()? == 0 {
        return a.call_method1("reshape", ((0, 3),));
    }
    Ok(a)
}

/// A 1-D array-like of floats.
pub(crate) fn floats(obj: &Bound<PyAny>, what: &str) -> PyResult<Vec<f64>> {
    let a: PyReadonlyArray1<f64> = asarray(obj, "float64")
        .and_then(|a| Ok(a.extract()?))
        .map_err(|_| invalid(format!("{what} must be a 1-D array of numbers")))?;
    Ok(a.as_array().to_vec())
}

/// A 1-D array-like of booleans.
pub(crate) fn bools(obj: &Bound<PyAny>, what: &str) -> PyResult<Vec<bool>> {
    let a: PyReadonlyArray1<bool> = asarray(obj, "bool")
        .and_then(|a| Ok(a.extract()?))
        .map_err(|_| invalid(format!("{what} must be a 1-D array of booleans")))?;
    Ok(a.as_array().to_vec())
}

/// Rows of a 2-D array-like of floats with `cols` columns, or 2 columns with
/// z = 0 when `allow2`.
pub(crate) fn rows3(obj: &Bound<PyAny>, what: &str, allow2: bool) -> PyResult<Vec<[f64; 3]>> {
    let shape = if allow2 { "(n, 3) or (n, 2)" } else { "(n, 3)" };
    let bad = || invalid(format!("{what} must be a {shape} array of numbers"));
    let a: PyReadonlyArray2<f64> = asarray2(obj, "float64")
        .and_then(|a| Ok(a.extract()?))
        .map_err(|_| bad())?;
    let a = a.as_array();
    match a.ncols() {
        3 => Ok(a.rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect()),
        2 if allow2 => Ok(a.rows().into_iter().map(|r| [r[0], r[1], 0.0]).collect()),
        _ if a.nrows() == 0 => Ok(vec![]),
        _ => Err(bad()),
    }
}

/// Points given as an array-like, a `PointSet` or a `Mesh` (its vertices).
pub(crate) fn points_arg(obj: &Bound<PyAny>, what: &str, allow2: bool) -> PyResult<Vec<[f64; 3]>> {
    if let Ok(m) = obj.cast::<PyMesh>() {
        return Ok(m.get().inner.vertices().to_vec());
    }
    if let Ok(p) = obj.cast::<PyPointSet>() {
        return Ok(p.get().inner.coords().to_vec());
    }
    rows3(obj, what, allow2)
}

/// `(crs, length_unit)` of a `PointSet` or `Mesh`, none for anything else.
pub(crate) fn labels_of(obj: &Bound<PyAny>) -> (Option<String>, Option<String>) {
    if let Ok(m) = obj.cast::<PyMesh>() {
        let m = &m.get().inner;
        return (m.crs.clone(), m.length_unit.clone());
    }
    if let Ok(p) = obj.cast::<PyPointSet>() {
        let p = &p.get().inner;
        return (p.crs.clone(), p.length_unit.clone());
    }
    (None, None)
}

/// Plan coordinates of points given as `(n, 2)` or `(n, 3)` (z ignored).
pub(crate) fn plan_arg(obj: &Bound<PyAny>, what: &str) -> PyResult<Vec<[f64; 2]>> {
    Ok(points_arg(obj, what, true)?
        .into_iter()
        .map(|p| [p[0], p[1]])
        .collect())
}

pub(crate) fn vec3(obj: &Bound<PyAny>, what: &str) -> PyResult<[f64; 3]> {
    let v = floats(obj, what)?;
    <[f64; 3]>::try_from(v).map_err(|_| invalid(format!("{what} must have 3 values")))
}

/// `a` with numpy's writeable flag cleared: it is a copy of immutable state, so
/// a write would be lost.
pub(crate) fn frozen<'py, T>(a: Bound<'py, T>) -> PyResult<Bound<'py, T>> {
    let kw = PyDict::new(a.py());
    kw.set_item("write", false)?;
    a.as_any().call_method("setflags", (), Some(&kw))?;
    Ok(a)
}

pub(crate) fn array1<'py>(py: Python<'py>, v: Vec<f64>) -> Bound<'py, PyArray1<f64>> {
    PyArray1::from_vec(py, v)
}

pub(crate) fn array_n<'py, const N: usize>(
    py: Python<'py>,
    rows: &[[f64; N]],
) -> Bound<'py, PyArray2<f64>> {
    Array2::from_shape_vec((rows.len(), N), rows.concat())
        .expect("rectangular rows")
        .into_pyarray(py)
}

fn bounds_of(b: Option<wireframes_core::Aabb>) -> Option<(Vec<f64>, Vec<f64>)> {
    b.map(|b| (b.min.to_vec(), b.max.to_vec()))
}

pub(crate) fn missing(name: &str, available: &[String]) -> PyErr {
    to_py(Error::MissingColumn {
        name: name.to_string(),
        available: available.to_vec(),
    })
}

/// `values` as a column, typed by `wireframes._containers._column`.
pub(crate) fn column_arg(values: &Bound<PyAny>, unit: Option<String>) -> PyResult<Column> {
    let py = values.py();
    let parsed = py
        .import("wireframes._containers")?
        .call_method1("_column", (values,))?;
    let (kind, data, valid): (String, Bound<PyAny>, Option<Bound<PyAny>>) = parsed.extract()?;
    let column = match kind.as_str() {
        "float64" => Column::f64(data.extract::<PyReadonlyArray1<f64>>()?.as_array().to_vec()),
        "int64" => Column::i64(data.extract::<PyReadonlyArray1<i64>>()?.as_array().to_vec()),
        "bool" => Column::bool(
            data.extract::<PyReadonlyArray1<bool>>()?
                .as_array()
                .to_vec(),
        ),
        _ => Column::utf8(data.extract::<Vec<String>>()?),
    };
    let column = match valid {
        Some(v) => column.with_valid(bools(&v, "valid")?).map_err(to_py)?,
        None => column,
    };
    Ok(column.with_unit(unit))
}

/// The values of `c` as a numpy array: float64 (NaN null), int64, bool, or
/// object (None null).
pub(crate) fn column_py<'py>(py: Python<'py>, c: &Column) -> PyResult<Bound<'py, PyAny>> {
    let valid = |i: usize| c.is_valid(i);
    Ok(match &c.data {
        ColumnData::Float64(v) => {
            let v: Vec<f64> = (0..v.len())
                .map(|i| if valid(i) { v[i] } else { f64::NAN })
                .collect();
            PyArray1::from_vec(py, v).into_any()
        }
        ColumnData::Int64(v) => PyArray1::from_slice(py, v).into_any(),
        ColumnData::Bool(v) => PyArray1::from_slice(py, v).into_any(),
        ColumnData::Utf8(v) => {
            let items: Vec<Option<&str>> = (0..v.len())
                .map(|i| valid(i).then_some(v[i].as_str()))
                .collect();
            let list = PyList::new(py, items)?;
            let np = numpy(py)?;
            let out = np.call_method1("empty", (v.len(), "object"))?;
            out.set_item(pyo3::types::PySlice::full(py), list)?;
            out
        }
    })
}

/// A table given as a `Table`, a mapping of name to values, or None (no
/// columns), with `rows` rows.
pub(crate) fn table_arg(obj: Option<&Bound<PyAny>>, rows: usize, what: &str) -> PyResult<Table> {
    let Some(obj) = obj.filter(|o| !o.is_none()) else {
        return Ok(Table::new(rows));
    };
    let table = match obj.cast::<PyTable>() {
        Ok(t) => t.get().inner.clone(),
        Err(_) => mapping_table(obj, Some(rows), what)?,
    };
    if table.rows() != rows {
        return Err(invalid(format!(
            "{what} has {} rows, expected {rows}",
            table.rows()
        )));
    }
    Ok(table)
}

fn mapping_table(obj: &Bound<PyAny>, rows: Option<usize>, what: &str) -> PyResult<Table> {
    let map = obj.cast::<PyMapping>().map_err(|_| {
        invalid(format!(
            "{what} must be a Table or a mapping of name to values"
        ))
    })?;
    let mut columns = vec![];
    for item in map.items()?.iter() {
        let (name, values): (String, Bound<PyAny>) = item.extract()?;
        columns.push((name, column_arg(&values, None)?));
    }
    let n = rows.unwrap_or_else(|| columns.first().map_or(0, |(_, c)| c.len()));
    let mut table = Table::new(n);
    for (name, c) in columns {
        if c.len() != n {
            return Err(invalid(format!(
                "{what}: column {name:?} has {} values, expected {n}",
                c.len()
            )));
        }
        table = table.with_column(&name, c).map_err(to_py)?;
    }
    Ok(table)
}

/// Columnar attribute table: float64, int64, bool and string columns, each
/// with a validity mask and an optional unit.
///
/// Parameters
/// ----------
/// data : mapping of str to array_like, optional
///     Column name to values. Floats become float64 (NaN is null), integers
///     int64, booleans bool, strings utf8; None marks a null.
///
/// Examples
/// --------
/// >>> t = wf.Table({"au": [0.5, None], "rock": ["ox", "fr"]})
/// >>> t.column_names, len(t), t["rock"].tolist()
/// (['au', 'rock'], 2, ['ox', 'fr'])
#[pyclass(name = "Table", module = "wireframes", frozen)]
pub(crate) struct PyTable {
    pub(crate) inner: Table,
}

impl From<Table> for PyTable {
    fn from(inner: Table) -> PyTable {
        PyTable { inner }
    }
}

impl PyTable {
    fn get_column(&self, name: &str) -> PyResult<&Column> {
        self.inner
            .column(name)
            .ok_or_else(|| missing(name, self.inner.names()))
    }
}

#[pymethods]
impl PyTable {
    #[new]
    #[pyo3(signature = (data=None))]
    fn new(data: Option<&Bound<PyAny>>) -> PyResult<PyTable> {
        match data.filter(|d| !d.is_none()) {
            None => Ok(Table::new(0).into()),
            Some(d) => Ok(mapping_table(d, None, "data")?.into()),
        }
    }

    /// Number of rows.
    #[getter]
    fn num_rows(&self) -> usize {
        self.inner.rows()
    }

    /// Column names in order.
    #[getter]
    fn column_names(&self) -> Vec<String> {
        self.inner.names().to_vec()
    }

    /// Units of the columns that have one.
    #[getter]
    fn units<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (n, c) in self.inner.names().iter().zip(self.inner.columns()) {
            if let Some(u) = &c.unit {
                d.set_item(n, u)?;
            }
        }
        Ok(d)
    }

    /// The values of column `name`.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     float64 with NaN for nulls, int64, bool, or object with None for
    ///     nulls; ``mask(name)`` tells which int64 and bool rows are null.
    ///
    /// Raises
    /// ------
    /// MissingColumn
    ///     When the table has no such column.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"hole": ["A", "B", "C"], "au": [1.0, None, 2.5]})
    /// >>> t.column("au").tolist()
    /// [1.0, nan, 2.5]
    fn column<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        frozen(column_py(py, self.get_column(name)?)?)
    }

    fn __getitem__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        self.column(py, name)
    }

    fn __contains__(&self, name: &str) -> bool {
        self.inner.column(name).is_some()
    }

    /// Validity of column `name`: True where the row holds a value.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"hole": ["A", "B", "C"], "au": [1.0, None, 2.5]})
    /// >>> t.mask("au").tolist()
    /// [True, False, True]
    fn mask<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyArray1<bool>>> {
        let c = self.get_column(name)?;
        Ok(PyArray1::from_vec(
            py,
            (0..c.len()).map(|i| c.is_valid(i)).collect(),
        ))
    }

    /// A copy with the rows where `mask` is True.
    ///
    /// Parameters
    /// ----------
    /// mask : array_like of bool
    ///     One value per row.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"hole": ["A", "B", "C"], "au": [1.0, None, 2.5]})
    /// >>> t.filter([True, False, True])["hole"].tolist()
    /// ['A', 'C']
    fn filter(&self, mask: &Bound<PyAny>) -> PyResult<PyTable> {
        let keep = bools(mask, "mask")?;
        Ok(self.inner.filter(&keep).map_err(to_py)?.into())
    }

    /// A copy with column `name` added, or replaced in place.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    /// values : array_like or sequence of str or None
    ///     One value per row; None or NaN marks a null.
    /// unit : str, optional
    ///     Unit of the column.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"hole": ["A", "B", "C"], "au": [1.0, None, 2.5]})
    /// >>> t.with_column("cu", [0.1, 0.2, 0.3], unit="%").column_names
    /// ['hole', 'au', 'cu']
    #[pyo3(signature = (name, values, *, unit=None))]
    fn with_column(
        &self,
        name: &str,
        values: &Bound<PyAny>,
        unit: Option<String>,
    ) -> PyResult<PyTable> {
        let c = column_arg(values, unit)?;
        let table = if self.inner.num_columns() == 0 && self.inner.rows() == 0 {
            Table::new(c.len())
        } else {
            self.inner.clone()
        };
        Ok(table.with_column(name, c).map_err(to_py)?.into())
    }

    /// A copy with every column of `data` added or replaced.
    ///
    /// Parameters
    /// ----------
    /// data : Table or mapping of str to array_like
    ///     Columns with one value per row.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"hole": ["A", "B", "C"], "au": [1.0, None, 2.5]})
    /// >>> t.with_columns({"cu": [0.1, 0.2, 0.3], "au": [1.0, 2.0, 3.0]})["au"].tolist()
    /// [1.0, 2.0, 3.0]
    fn with_columns(&self, data: &Bound<PyAny>) -> PyResult<PyTable> {
        let more = mapping_table(data, None, "data")?;
        let mut table = if self.inner.num_columns() == 0 && self.inner.rows() == 0 {
            Table::new(more.rows())
        } else {
            self.inner.clone()
        };
        for (name, c) in more.names().iter().zip(more.columns()) {
            table = table.with_column(name, c.clone()).map_err(to_py)?;
        }
        Ok(table.into())
    }

    /// Column name to numpy values, as ``column`` gives them.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"hole": ["A", "B", "C"], "au": [1.0, None, 2.5]})
    /// >>> t.to_dict()["hole"].tolist()
    /// ['A', 'B', 'C']
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (name, c) in self.inner.names().iter().zip(self.inner.columns()) {
            d.set_item(name, column_py(py, c)?)?;
        }
        Ok(d)
    }

    fn __len__(&self) -> usize {
        self.inner.rows()
    }

    fn __repr__(&self) -> String {
        let cols: Vec<String> = self
            .inner
            .names()
            .iter()
            .zip(self.inner.columns())
            .map(|(n, c)| format!("{n}: {}", c.type_name()))
            .collect();
        format!("Table({} rows; {})", self.inner.rows(), cols.join(", "))
    }
}

/// Triangle mesh: an open surface (topography, fault, pit shell) or a closed
/// solid (orebody, stope), with vertex and face attribute tables.
///
/// Parameters
/// ----------
/// vertices : array_like
///     ``(n, 3)`` finite coordinates.
/// triangles : array_like
///     ``(m, 3)`` vertex indices; a solid winds counter-clockwise seen from
///     outside.
/// vertex_attributes, face_attributes : Table or mapping, optional
///     One row per vertex and one per triangle.
/// crs : str, optional
///     Coordinate reference system label, checked for equality when two
///     geometries combine.
/// length_unit : str, optional
///     Length unit label, such as ``"m"`` or ``"ft"``.
///
/// Examples
/// --------
/// >>> corners = [[0, 0, 0], [6, 0, 0], [0, 6, 0], [0, 0, 6]]
/// >>> tet = wf.Mesh(corners, [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
/// >>> tet.is_solid, tet.volume
/// (True, 36.0)
#[pyclass(name = "Mesh", module = "wireframes", frozen)]
pub(crate) struct PyMesh {
    pub(crate) inner: Mesh,
    solid: OnceLock<SolidIndex>,
    surface: OnceLock<SurfaceIndex>,
}

impl From<Mesh> for PyMesh {
    fn from(inner: Mesh) -> PyMesh {
        PyMesh {
            inner,
            solid: OnceLock::new(),
            surface: OnceLock::new(),
        }
    }
}

impl PyMesh {
    /// The cached point-in-solid index; InvalidMesh when the mesh is not a solid.
    pub(crate) fn solid_index(&self, py: Python) -> PyResult<&SolidIndex> {
        if let Some(index) = self.solid.get() {
            return Ok(index);
        }
        let index = py.detach(|| SolidIndex::new(&self.inner)).map_err(to_py)?;
        Ok(self.solid.get_or_init(|| index))
    }

    pub(crate) fn surface_index(&self, py: Python) -> &SurfaceIndex {
        if let Some(index) = self.surface.get() {
            return index;
        }
        let index = py.detach(|| SurfaceIndex::new(&self.inner));
        self.surface.get_or_init(|| index)
    }

    fn with(&self, f: impl FnOnce(Mesh) -> wireframes_core::Result<Mesh>) -> PyResult<PyMesh> {
        Ok(f(self.inner.clone()).map_err(to_py)?.into())
    }
}

fn triangles_arg(obj: &Bound<PyAny>) -> PyResult<Vec<[u32; 3]>> {
    let bad = || invalid("triangles must be an (m, 3) array of vertex indices");
    let dtype = numpy(obj.py())?
        .call_method1("asarray", (obj,))
        .and_then(|a| a.getattr("dtype"))
        .map_err(|_| bad())?;
    let a: Array2<i64> = match dtype.getattr("kind")?.extract::<String>()?.as_str() {
        "i" | "u" => {
            let a: PyReadonlyArray2<i64> = asarray2(obj, "int64")
                .and_then(|a| Ok(a.extract()?))
                .map_err(|_| bad())?;
            a.as_array().to_owned()
        }
        "f" => {
            let a: PyReadonlyArray2<f64> = asarray2(obj, "float64")
                .and_then(|a| Ok(a.extract()?))
                .map_err(|_| bad())?;
            if let Some(x) = a.as_array().iter().find(|x| x.fract() != 0.0) {
                return Err(invalid(format!(
                    "triangles must be whole vertex indices, but one is {x}"
                )));
            }
            a.as_array().mapv(|x| x as i64)
        }
        _ => {
            return Err(invalid(format!(
                "triangles must be integer vertex indices, not dtype {}",
                dtype.str()?
            )));
        }
    };
    if a.nrows() == 0 {
        return Ok(vec![]);
    }
    if a.ncols() != 3 {
        return Err(bad());
    }
    a.rows()
        .into_iter()
        .enumerate()
        .map(|(f, r)| {
            let mut t = [0u32; 3];
            for k in 0..3 {
                t[k] = u32::try_from(r[k]).map_err(|_| {
                    invalid(format!(
                        "triangles must index existing vertices: face {f} uses {}",
                        r[k]
                    ))
                })?;
            }
            Ok(t)
        })
        .collect()
}

fn on_arg(on: &str) -> PyResult<bool> {
    match on {
        "vertices" => Ok(true),
        "faces" => Ok(false),
        _ => Err(invalid(format!(
            "on must be \"vertices\" or \"faces\", not {on:?}"
        ))),
    }
}

fn append(mut table: Table, more: &Table) -> PyResult<Table> {
    for (name, c) in more.names().iter().zip(more.columns()) {
        table = table.with_column(name, c.clone()).map_err(to_py)?;
    }
    Ok(table)
}

fn coord_table(points: &[[f64; 3]]) -> PyResult<Table> {
    let mut t = Table::new(points.len());
    for (k, name) in ["x", "y", "z"].iter().enumerate() {
        t = t
            .with_column(name, Column::f64(points.iter().map(|p| p[k]).collect()))
            .map_err(to_py)?;
    }
    Ok(t)
}

#[pymethods]
impl PyMesh {
    #[new]
    #[pyo3(signature = (vertices, triangles, *, vertex_attributes=None, face_attributes=None, crs=None, length_unit=None))]
    fn new(
        vertices: &Bound<PyAny>,
        triangles: &Bound<PyAny>,
        vertex_attributes: Option<&Bound<PyAny>>,
        face_attributes: Option<&Bound<PyAny>>,
        crs: Option<String>,
        length_unit: Option<String>,
    ) -> PyResult<PyMesh> {
        let v = rows3(vertices, "vertices", false)?;
        let t = triangles_arg(triangles)?;
        let (nv, nf) = (v.len(), t.len());
        let vt = table_arg(vertex_attributes, nv, "vertex_attributes")?;
        let ft = table_arg(face_attributes, nf, "face_attributes")?;
        let mesh = Mesh::new(v, t)
            .and_then(|m| m.with_attributes(vt, ft))
            .map_err(to_py)?
            .with_crs(crs)
            .with_length_unit(length_unit);
        Ok(mesh.into())
    }

    /// ``(n, 3)`` vertex coordinates.
    #[getter]
    fn coords<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        frozen(array_n(py, self.inner.vertices()))
    }

    /// ``(m, 3)`` vertex indices of each triangle.
    #[getter]
    fn triangles<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<i64>>> {
        let flat: Vec<i64> = self
            .inner
            .triangles()
            .iter()
            .flat_map(|t| t.map(i64::from))
            .collect();
        frozen(
            Array2::from_shape_vec((self.inner.n_faces(), 3), flat)
                .expect("rectangular rows")
                .into_pyarray(py),
        )
    }

    /// Number of vertices.
    #[getter]
    fn n_vertices(&self) -> usize {
        self.inner.n_vertices()
    }

    /// Number of triangles.
    #[getter]
    fn n_faces(&self) -> usize {
        self.inner.n_faces()
    }

    /// ``(min, max)`` corners over every vertex, None without vertices.
    #[getter]
    fn bounds(&self) -> Option<(Vec<f64>, Vec<f64>)> {
        bounds_of(self.inner.bounds())
    }

    /// Names of the vertex columns, then the face columns.
    #[getter]
    fn column_names(&self) -> Vec<String> {
        let mut names = self.inner.vertex_attributes().names().to_vec();
        names.extend_from_slice(self.inner.face_attributes().names());
        names
    }

    /// One row per vertex.
    #[getter]
    fn vertex_attributes(&self) -> PyTable {
        self.inner.vertex_attributes().clone().into()
    }

    /// One row per triangle.
    #[getter]
    fn face_attributes(&self) -> PyTable {
        self.inner.face_attributes().clone().into()
    }

    /// Coordinate reference system label.
    #[getter]
    fn crs(&self) -> Option<String> {
        self.inner.crs.clone()
    }

    /// Length unit label.
    #[getter]
    fn length_unit(&self) -> Option<String> {
        self.inner.length_unit.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.n_vertices()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        let v = self.inner.vertex_attributes().column(name);
        let f = self.inner.face_attributes().column(name);
        match (v, f) {
            (Some(_), Some(_)) => Err(invalid(format!(
                "{name:?} is both a vertex and a face column; use vertex_attributes or face_attributes"
            ))),
            (Some(c), None) | (None, Some(c)) => frozen(column_py(py, c)?),
            (None, None) => Err(missing(name, &self.column_names())),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Mesh({} vertices, {} triangles, columns {:?}, crs {}, length_unit {})",
            self.inner.n_vertices(),
            self.inner.n_faces(),
            self.column_names(),
            self.inner.crs.as_deref().unwrap_or("None"),
            self.inner.length_unit.as_deref().unwrap_or("None")
        )
    }

    /// A copy with vertex column `name` added or replaced.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    /// values : array_like or sequence of str or None
    ///     One value per vertex.
    /// unit : str, optional
    ///     Unit of the column.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.with_vertex_column("grade", np.arange(8.0), unit="g/t").column_names
    /// ['grade']
    #[pyo3(signature = (name, values, *, unit=None))]
    fn with_vertex_column(
        &self,
        name: &str,
        values: &Bound<PyAny>,
        unit: Option<String>,
    ) -> PyResult<PyMesh> {
        let c = column_arg(values, unit)?;
        self.with(|m| m.with_vertex_column(name, c))
    }

    /// A copy with face column `name` added or replaced.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    /// values : array_like or sequence of str or None
    ///     One value per triangle.
    /// unit : str, optional
    ///     Unit of the column.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.with_face_column("domain", np.ones(12, dtype=np.int64)).column_names
    /// ['domain']
    #[pyo3(signature = (name, values, *, unit=None))]
    fn with_face_column(
        &self,
        name: &str,
        values: &Bound<PyAny>,
        unit: Option<String>,
    ) -> PyResult<PyMesh> {
        let c = column_arg(values, unit)?;
        self.with(|m| m.with_face_column(name, c))
    }

    /// A copy with new attribute tables; None keeps a table as it is.
    ///
    /// Parameters
    /// ----------
    /// vertex : Table or mapping of str to array_like, optional
    ///     One row per vertex.
    /// face : Table or mapping of str to array_like, optional
    ///     One row per triangle.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.with_attributes(face={"domain": np.ones(12)}).column_names
    /// ['domain']
    #[pyo3(signature = (*, vertex=None, face=None))]
    fn with_attributes(
        &self,
        vertex: Option<&Bound<PyAny>>,
        face: Option<&Bound<PyAny>>,
    ) -> PyResult<PyMesh> {
        let vt = match vertex {
            Some(v) => table_arg(Some(v), self.inner.n_vertices(), "vertex")?,
            None => self.inner.vertex_attributes().clone(),
        };
        let ft = match face {
            Some(f) => table_arg(Some(f), self.inner.n_faces(), "face")?,
            None => self.inner.face_attributes().clone(),
        };
        self.with(|m| m.with_attributes(vt, ft))
    }

    /// A copy with the vertices moved to `vertices`; triangles and attributes kept.
    ///
    /// Parameters
    /// ----------
    /// vertices : array_like
    ///     ``(n, 3)`` coordinates, as many as the mesh has.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.with_vertices(m.coords * 2).volume
    /// 64.0
    fn with_vertices(&self, vertices: &Bound<PyAny>) -> PyResult<PyMesh> {
        let v = rows3(vertices, "vertices", false)?;
        self.with(|m| m.with_vertices(v))
    }

    /// A copy with the triangles where `mask` is True, and only the vertices they use.
    ///
    /// Parameters
    /// ----------
    /// mask : array_like of bool
    ///     One value per triangle.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.filter_faces(m.face_normals[:, 2] > 0).n_faces
    /// 2
    fn filter_faces(&self, mask: &Bound<PyAny>) -> PyResult<PyMesh> {
        let keep = bools(mask, "mask")?;
        Ok(self.inner.filter_faces(&keep).map_err(to_py)?.into())
    }

    /// Vertices (``x``, ``y``, ``z`` and the vertex columns) or faces (``v0``,
    /// ``v1``, ``v2`` and the face columns) as a table.
    ///
    /// Parameters
    /// ----------
    /// on : {"vertices", "faces"}
    ///     Which rows.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.to_table().column_names, m.to_table(on="faces").num_rows
    /// (['x', 'y', 'z'], 12)
    #[pyo3(signature = (*, on="vertices"))]
    fn to_table(&self, on: &str) -> PyResult<PyTable> {
        let m = &self.inner;
        let table = if on_arg(on)? {
            append(coord_table(m.vertices())?, m.vertex_attributes())?
        } else {
            let mut t = Table::new(m.n_faces());
            for k in 0..3 {
                let c = Column::i64(m.triangles().iter().map(|f| i64::from(f[k])).collect());
                t = t.with_column(&format!("v{k}"), c).map_err(to_py)?;
            }
            append(t, m.face_attributes())?
        };
        Ok(table.into())
    }
}

/// 3D strings: open lines and closed rings grouped into features, with
/// feature and vertex attribute tables.
///
/// Parameters
/// ----------
/// parts : sequence of array_like
///     ``(k, 3)`` or ``(k, 2)`` (z = 0) vertices of each part; a ring does not
///     repeat its first vertex.
/// closed : bool or sequence of bool, default False
///     Whether each part is a ring.
/// features : sequence of int, optional
///     Feature index of each part, starting at 0 and increasing by at most 1;
///     None makes each part its own feature. Holes and multipart features are
///     several parts of one feature.
/// attributes : Table or mapping, optional
///     One row per feature.
/// vertex_attributes : Table or mapping, optional
///     One row per vertex.
/// crs, length_unit : str, optional
///     Labels.
///
/// Examples
/// --------
/// >>> crest = [[0, 0, 100], [4, 0, 100], [4, 3, 100], [0, 3, 100]]
/// >>> rings = wf.Polylines([crest], closed=True, attributes={"name": ["crest"]})
/// >>> rings["name"].tolist(), rings.area().tolist()
/// (['crest'], [12.0])
#[pyclass(name = "Polylines", module = "wireframes", frozen)]
pub(crate) struct PyPolylines {
    pub(crate) inner: Polylines,
}

impl From<Polylines> for PyPolylines {
    fn from(inner: Polylines) -> PyPolylines {
        PyPolylines { inner }
    }
}

#[pymethods]
impl PyPolylines {
    #[new]
    #[pyo3(
        signature = (parts, *, closed=None, features=None, attributes=None, vertex_attributes=None, crs=None, length_unit=None),
        text_signature = "(parts, *, closed=False, features=None, attributes=None, vertex_attributes=None, crs=None, length_unit=None)"
    )]
    fn new(
        parts: &Bound<PyAny>,
        closed: Option<&Bound<PyAny>>,
        features: Option<Vec<u32>>,
        attributes: Option<&Bound<PyAny>>,
        vertex_attributes: Option<&Bound<PyAny>>,
        crs: Option<String>,
        length_unit: Option<String>,
    ) -> PyResult<PyPolylines> {
        let parts: Vec<Vec<[f64; 3]>> = parts
            .try_iter()?
            .enumerate()
            .map(|(i, p)| rows3(&p?, &format!("part {i}"), true))
            .collect::<PyResult<_>>()?;
        let n = parts.len();
        let closed = match closed {
            None => vec![false; n],
            Some(c) => match c.extract::<bool>() {
                Ok(b) => vec![b; n],
                Err(_) => bools(c, "closed")?,
            },
        };
        let p = Polylines::from_parts(parts, closed, features).map_err(to_py)?;
        let at = table_arg(attributes, p.len(), "attributes")?;
        let vt = table_arg(vertex_attributes, p.n_vertices(), "vertex_attributes")?;
        let p = p
            .with_attributes(at)
            .and_then(|p| p.with_vertex_attributes(vt))
            .map_err(to_py)?
            .with_crs(crs)
            .with_length_unit(length_unit);
        Ok(p.into())
    }

    /// ``(n, 3)`` coordinates of every vertex, part after part.
    #[getter]
    fn coords<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        frozen(array_n(py, self.inner.vertices()))
    }

    /// ``(k, 3)`` vertices of each part.
    #[getter]
    fn parts<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyArray2<f64>>>> {
        (0..self.inner.n_parts())
            .map(|p| frozen(array_n(py, self.inner.part(p))))
            .collect()
    }

    /// Whether each part is a ring.
    #[getter]
    fn closed<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<bool>>> {
        frozen(PyArray1::from_slice(py, self.inner.closed()))
    }

    /// Feature index of each part.
    #[getter]
    fn feature<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<i64>>> {
        let mut out = Vec::with_capacity(self.inner.n_parts());
        for f in 0..self.inner.len() {
            out.extend(self.inner.feature_parts(f).map(|_| f as i64));
        }
        frozen(PyArray1::from_vec(py, out))
    }

    /// One row per feature.
    #[getter]
    fn attributes(&self) -> PyTable {
        self.inner.attributes().clone().into()
    }

    /// One row per vertex.
    #[getter]
    fn vertex_attributes(&self) -> PyTable {
        self.inner.vertex_attributes().clone().into()
    }

    /// ``(min, max)`` corners, None without vertices.
    #[getter]
    fn bounds(&self) -> Option<(Vec<f64>, Vec<f64>)> {
        bounds_of(self.inner.bounds())
    }

    /// Names of the feature columns, then the vertex columns.
    #[getter]
    fn column_names(&self) -> Vec<String> {
        let mut names = self.inner.attributes().names().to_vec();
        names.extend_from_slice(self.inner.vertex_attributes().names());
        names
    }

    /// Coordinate reference system label.
    #[getter]
    fn crs(&self) -> Option<String> {
        self.inner.crs.clone()
    }

    /// Length unit label.
    #[getter]
    fn length_unit(&self) -> Option<String> {
        self.inner.length_unit.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// The feature or vertex column `name`, as on Mesh and PointSet.
    fn __getitem__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        let f = self.inner.attributes().column(name);
        let v = self.inner.vertex_attributes().column(name);
        match (f, v) {
            (Some(_), Some(_)) => Err(invalid(format!(
                "{name:?} is both a feature and a vertex column; use attributes or vertex_attributes"
            ))),
            (Some(c), None) | (None, Some(c)) => frozen(column_py(py, c)?),
            (None, None) => Err(missing(name, &self.column_names())),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Polylines({} features, {} parts, {} vertices, columns {:?})",
            self.inner.len(),
            self.inner.n_parts(),
            self.inner.n_vertices(),
            self.column_names()
        )
    }

    /// A copy with the features where `mask` is True.
    ///
    /// Parameters
    /// ----------
    /// mask : array_like of bool
    ///     One value per feature.
    ///
    /// Examples
    /// --------
    /// >>> holes = wf.datasets.drillholes(n=3)
    /// >>> len(holes.filter([True, False, True]))
    /// 2
    fn filter(&self, mask: &Bound<PyAny>) -> PyResult<PyPolylines> {
        let keep = bools(mask, "mask")?;
        Ok(self.inner.filter(&keep).map_err(to_py)?.into())
    }

    /// A copy with feature column `name` added or replaced.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    /// values : array_like or sequence of str or None
    ///     One value per feature; None or NaN marks a null.
    /// unit : str, optional
    ///     Unit of the column.
    ///
    /// Examples
    /// --------
    /// >>> holes = wf.datasets.drillholes(n=3)
    /// >>> holes.with_column("hole", ["DH1", "DH2", "DH3"]).attributes["hole"].tolist()
    /// ['DH1', 'DH2', 'DH3']
    #[pyo3(signature = (name, values, *, unit=None))]
    fn with_column(
        &self,
        name: &str,
        values: &Bound<PyAny>,
        unit: Option<String>,
    ) -> PyResult<PyPolylines> {
        let c = column_arg(values, unit)?;
        Ok(self
            .inner
            .clone()
            .with_column(name, c)
            .map_err(to_py)?
            .into())
    }

    /// A copy with vertex column `name` added or replaced.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    /// values : array_like or sequence of str or None
    ///     One value per vertex; None or NaN marks a null.
    /// unit : str, optional
    ///     Unit of the column.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.with_vertex_column("id", [1, 2, 3, 4]).column_names
    /// ['id']
    #[pyo3(signature = (name, values, *, unit=None))]
    fn with_vertex_column(
        &self,
        name: &str,
        values: &Bound<PyAny>,
        unit: Option<String>,
    ) -> PyResult<PyPolylines> {
        let c = column_arg(values, unit)?;
        Ok(self
            .inner
            .clone()
            .with_vertex_column(name, c)
            .map_err(to_py)?
            .into())
    }

    /// 3D length of each feature, closing segments of rings included.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.length().tolist()
    /// [14.0]
    fn length<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        array1(py, self.inner.length())
    }

    /// Plan area of each feature, even-odd over its rings.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.area().tolist()
    /// [12.0]
    fn area<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        array1(py, self.inner.area())
    }

    /// Whether each point lies inside, in plan.
    ///
    /// Inside is even-odd over the rings of one feature; a point on a ring
    /// edge counts as inside. The test is exact.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 2)`` or ``(n, 3)`` points; z is ignored.
    /// feature : int, optional
    ///     Feature to test; None tests every feature.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.contains([[1, 1], [5, 5]]).tolist()
    /// [True, False]
    #[pyo3(signature = (points, *, feature=None))]
    fn contains<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
        feature: Option<usize>,
    ) -> PyResult<Bound<'py, PyArray1<bool>>> {
        let xy = plan_arg(points, "points")?;
        if let Some(f) = feature
            && f >= self.inner.len()
        {
            return Err(invalid(format!(
                "feature {f} does not exist; there are {}",
                self.inner.len()
            )));
        }
        let p = &self.inner;
        let out: Vec<bool> = py.detach(|| {
            xy.iter()
                .map(|&q| match feature {
                    Some(f) => p.feature_contains(f, q),
                    None => p.locate(q).is_some(),
                })
                .collect()
        });
        Ok(PyArray1::from_vec(py, out))
    }

    /// The first feature whose rings hold each point in plan, -1 for none.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 2)`` or ``(n, 3)`` points; z is ignored.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.locate([[1, 1], [5, 5]]).tolist()
    /// [0, -1]
    fn locate<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<i64>>> {
        let xy = plan_arg(points, "points")?;
        let p = &self.inner;
        let out: Vec<i64> = py.detach(|| {
            xy.iter()
                .map(|&q| p.locate(q).map_or(-1, |f| f as i64))
                .collect()
        });
        Ok(PyArray1::from_vec(py, out))
    }

    /// Plan distance from each point to the nearest ring edge.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 2)`` or ``(n, 3)`` points; z is ignored.
    /// feature : int, optional
    ///     Feature to measure to; None measures to every feature.
    /// signed : bool
    ///     Negative inside.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     NaN when there is no ring.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.distance([[1, 1], [2, 5]], signed=True).tolist()
    /// [-1.0, 2.0]
    #[pyo3(signature = (points, *, feature=None, signed=false))]
    fn distance<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
        feature: Option<usize>,
        signed: bool,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let xy = plan_arg(points, "points")?;
        let p = &self.inner;
        let out: Vec<f64> = py.detach(|| {
            xy.iter()
                .map(|&q| p.distance(q, feature, signed).unwrap_or(f64::NAN))
                .collect()
        });
        Ok(array1(py, out))
    }

    /// Every vertex as a row: ``feature``, ``part``, ``vertex``, ``x``, ``y``,
    /// ``z``, ``closed``, the vertex columns and the feature columns.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
    /// >>> ring.to_table().num_rows
    /// 4
    fn to_table(&self) -> PyResult<PyTable> {
        let p = &self.inner;
        let n = p.n_vertices();
        let (mut feature, mut part, mut vertex, mut closed) = (
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
        );
        for f in 0..p.len() {
            for q in p.feature_parts(f) {
                for (k, _) in p.part_range(q).enumerate() {
                    feature.push(f as i64);
                    part.push(q as i64);
                    vertex.push(k as i64);
                    closed.push(p.closed()[q]);
                }
            }
        }
        let rows: Vec<u32> = feature.iter().map(|&f| f as u32).collect();
        let mut t = Table::new(n)
            .with_column("feature", Column::i64(feature))
            .and_then(|t| t.with_column("part", Column::i64(part)))
            .and_then(|t| t.with_column("vertex", Column::i64(vertex)))
            .map_err(to_py)?;
        t = append(t, &coord_table(p.vertices())?)?;
        t = t
            .with_column("closed", Column::bool(closed))
            .map_err(to_py)?;
        t = append(t, p.vertex_attributes())?;
        t = append(t, &p.attributes().take(&rows))?;
        Ok(t.into())
    }
}

/// 3D point cloud with an attribute table: collars, pickups, LiDAR.
///
/// Parameters
/// ----------
/// coords : array_like
///     ``(n, 3)`` finite coordinates.
/// attributes : Table or mapping, optional
///     One row per point.
/// crs, length_unit : str, optional
///     Labels.
///
/// Examples
/// --------
/// >>> collars = wf.PointSet([[10, 20, 300], [15, 20, 298]], {"hole": ["DH1", "DH2"]})
/// >>> len(collars), collars["hole"].tolist(), collars.z.tolist()
/// (2, ['DH1', 'DH2'], [300.0, 298.0])
#[pyclass(name = "PointSet", module = "wireframes", frozen)]
pub(crate) struct PyPointSet {
    pub(crate) inner: PointSet,
}

impl From<PointSet> for PyPointSet {
    fn from(inner: PointSet) -> PyPointSet {
        PyPointSet { inner }
    }
}

#[pymethods]
impl PyPointSet {
    #[new]
    #[pyo3(signature = (coords, attributes=None, *, crs=None, length_unit=None))]
    fn new(
        coords: &Bound<PyAny>,
        attributes: Option<&Bound<PyAny>>,
        crs: Option<String>,
        length_unit: Option<String>,
    ) -> PyResult<PyPointSet> {
        let c = rows3(coords, "coords", false)?;
        let at = table_arg(attributes, c.len(), "attributes")?;
        let p = PointSet::new(c)
            .and_then(|p| p.with_attributes(at))
            .map_err(to_py)?
            .with_crs(crs)
            .with_length_unit(length_unit);
        Ok(p.into())
    }

    /// ``(n, 3)`` coordinates.
    #[getter]
    fn coords<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        frozen(array_n(py, self.inner.coords()))
    }

    /// One row per point.
    #[getter]
    fn attributes(&self) -> PyTable {
        self.inner.attributes().clone().into()
    }

    /// ``(min, max)`` corners, None without points.
    #[getter]
    fn bounds(&self) -> Option<(Vec<f64>, Vec<f64>)> {
        bounds_of(self.inner.bounds())
    }

    /// Names of the attribute columns.
    #[getter]
    fn column_names(&self) -> Vec<String> {
        self.inner.attributes().names().to_vec()
    }

    /// Coordinate reference system label.
    #[getter]
    fn crs(&self) -> Option<String> {
        self.inner.crs.clone()
    }

    /// Length unit label.
    #[getter]
    fn length_unit(&self) -> Option<String> {
        self.inner.length_unit.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        match self.inner.attributes().column(name) {
            Some(c) => frozen(column_py(py, c)?),
            None => Err(missing(name, self.inner.attributes().names())),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PointSet({} points, columns {:?})",
            self.inner.len(),
            self.column_names()
        )
    }

    /// A copy with the points where `mask` is True.
    ///
    /// Parameters
    /// ----------
    /// mask : array_like of bool
    ///     One value per point.
    ///
    /// Examples
    /// --------
    /// >>> p = wf.PointSet([[0, 0, 0], [1, 1, 1]], {"au": [1.5, 2.0]})
    /// >>> len(p.filter([True, False]))
    /// 1
    fn filter(&self, mask: &Bound<PyAny>) -> PyResult<PyPointSet> {
        let keep = bools(mask, "mask")?;
        Ok(self.inner.filter(&keep).map_err(to_py)?.into())
    }

    /// A copy with column `name` added or replaced.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Column name.
    /// values : array_like or sequence of str or None
    ///     One value per point; None or NaN marks a null.
    /// unit : str, optional
    ///     Unit of the column.
    ///
    /// Examples
    /// --------
    /// >>> p = wf.PointSet([[0, 0, 0], [1, 1, 1]], {"au": [1.5, 2.0]})
    /// >>> p.with_column("cu", [0.1, 0.2], unit="%").column_names
    /// ['au', 'cu']
    #[pyo3(signature = (name, values, *, unit=None))]
    fn with_column(
        &self,
        name: &str,
        values: &Bound<PyAny>,
        unit: Option<String>,
    ) -> PyResult<PyPointSet> {
        let c = column_arg(values, unit)?;
        Ok(self
            .inner
            .clone()
            .with_column(name, c)
            .map_err(to_py)?
            .into())
    }

    /// ``x``, ``y``, ``z`` and the attribute columns.
    ///
    /// Examples
    /// --------
    /// >>> p = wf.PointSet([[0, 0, 0], [1, 1, 1]], {"au": [1.5, 2.0]})
    /// >>> p.to_table().column_names
    /// ['x', 'y', 'z', 'au']
    fn to_table(&self) -> PyResult<PyTable> {
        Ok(append(coord_table(self.inner.coords())?, self.inner.attributes())?.into())
    }
}

/// Volume and mixed cells read from ``.vtu`` or ``.vtk``, kept as data so
/// finite-element and reservoir meshes round-trip; not a modeling type.
#[pyclass(name = "UnstructuredGrid", module = "wireframes", frozen)]
pub(crate) struct PyUnstructuredGrid {
    pub(crate) inner: UnstructuredGrid,
}

#[pymethods]
impl PyUnstructuredGrid {
    /// ``(n, 3)`` point coordinates.
    #[getter]
    fn points<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        frozen(array_n(py, &self.inner.points))
    }

    /// VTK cell type id of each cell.
    #[getter]
    fn cell_types<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        frozen(PyArray1::from_slice(py, &self.inner.cell_types))
    }

    /// ``n_cells + 1`` offsets into ``connectivity``.
    #[getter]
    fn offsets<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<i64>>> {
        frozen(PyArray1::from_slice(py, &self.inner.offsets))
    }

    /// Point ids of every cell, cell after cell.
    #[getter]
    fn connectivity<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<i64>>> {
        frozen(PyArray1::from_slice(py, &self.inner.connectivity))
    }

    /// One row per point.
    #[getter]
    fn point_data(&self) -> PyTable {
        self.inner.point_data.clone().into()
    }

    /// One row per cell.
    #[getter]
    fn cell_data(&self) -> PyTable {
        self.inner.cell_data.clone().into()
    }

    /// Number of cells.
    #[getter]
    fn n_cells(&self) -> usize {
        self.inner.n_cells()
    }

    /// Triangles of the surface cells (VTK types 5 to 9), with every point and
    /// the cell data of each triangle's cell.
    ///
    /// Examples
    /// --------
    /// >>> import pathlib, tempfile
    /// >>> folder = pathlib.Path(tempfile.mkdtemp())
    /// >>> vtk = "# vtk DataFile Version 2.0\ncells\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 4 double\n0 0 0 1 0 0 0 1 0 0 0 1\nCELLS 2 9\n4 0 1 2 3\n3 0 1 2\nCELL_TYPES 2\n10\n5\n"
    /// >>> _ = (folder / "cells.vtk").write_text(vtk)
    /// >>> grid = wf.read(folder / "cells.vtk").grid
    /// >>> grid.n_cells, grid.surface().n_faces
    /// (2, 1)
    fn surface(&self) -> PyResult<PyMesh> {
        Ok(self.inner.surface().map_err(to_py)?.into())
    }

    fn __repr__(&self) -> String {
        format!(
            "UnstructuredGrid({} points, {} cells)",
            self.inner.points.len(),
            self.inner.n_cells()
        )
    }
}

/// Meshes appended in order: vertices, triangles (indices offset) and
/// attributes; no welding.
///
/// Parameters
/// ----------
/// meshes : sequence of Mesh
///     Meshes with equal or undeclared ``crs`` and ``length_unit``.
///
/// Returns
/// -------
/// Mesh
///     A column missing from one mesh is null in its rows.
///
/// Examples
/// --------
/// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
/// >>> wf.concat([a, b]).n_faces
/// 24
#[pyfunction]
pub(crate) fn concat(meshes: &Bound<PyAny>) -> PyResult<PyMesh> {
    let items: Vec<Bound<PyAny>> = meshes.try_iter()?.collect::<PyResult<_>>()?;
    let refs: Vec<&Bound<PyMesh>> = items
        .iter()
        .map(|m| m.cast::<PyMesh>())
        .collect::<Result<_, _>>()?;
    let inner: Vec<&Mesh> = refs.iter().map(|m| &m.get().inner).collect();
    Ok(Mesh::concat(&inner).map_err(to_py)?.into())
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<PyTable>()?;
    m.add_class::<PyMesh>()?;
    m.add_class::<PyPolylines>()?;
    m.add_class::<PyPointSet>()?;
    m.add_class::<PyUnstructuredGrid>()?;
    m.add_function(wrap_pyfunction!(concat, m)?)?;
    Ok(())
}
