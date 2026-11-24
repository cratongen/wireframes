//! `__minetypes__` protocol, version 1. Arrays that already match the stored
//! buffer are read-only views of it. A validity mask, a text column, or a
//! widened nullable int is copied; the copy is named in the method docstring.

use std::collections::HashMap;

use numpy::ndarray::{ArrayView1, ArrayView2};
use numpy::{Element, PyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use wireframes_core::{Column, ColumnData, Table};

use crate::containers::{PyMesh, PyPointSet, PyPolylines, PyTable, frozen, numpy};
use crate::errors::invalid;

fn seal<'py, T, D>(array: &Bound<'py, PyArray<T, D>>)
where
    T: Element,
    D: numpy::ndarray::Dimension,
{
    array.readwrite().make_nonwriteable();
}

fn view1<'py, T>(data: &[T], owner: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>>
where
    T: Element,
{
    let view = ArrayView1::from_shape(data.len(), data).map_err(|err| invalid(err.to_string()))?;
    // SAFETY: `data` is a field of the frozen object `owner`. That object does
    // not reallocate the field, and the array holds `owner` for its lifetime.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner.clone()) };
    seal(&array);
    Ok(array.into_any())
}

fn view_rows<'py, T, const N: usize>(
    rows: &[[T; N]],
    owner: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>>
where
    T: Element,
{
    debug_assert_eq!(size_of::<[T; N]>(), N * size_of::<T>());
    let n = rows.len();
    let len = n
        .checked_mul(N)
        .ok_or_else(|| invalid("coords is too large to view"))?;
    let view = ArrayView2::from_shape((n, N), unsafe {
        // SAFETY: `[T; N]` is `N` contiguous `T`s. `rows` lives in `owner`,
        // which the array holds, and the frozen object does not reallocate it.
        std::slice::from_raw_parts(rows.as_ptr().cast::<T>(), len)
    })
    .map_err(|err| invalid(err.to_string()))?;
    let array = unsafe { PyArray2::borrow_from_array(&view, owner.clone()) };
    seal(&array);
    Ok(array.into_any())
}

fn owned<'py>(array: Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    frozen(array)
}

fn unicode<'py>(py: Python<'py>, values: &[String]) -> PyResult<Bound<'py, PyAny>> {
    let list = PyList::new(py, values.iter().map(String::as_str))?;
    let np = numpy(py)?;
    let array = np.call_method1("asarray", (list, np.getattr("str_")?))?;
    owned(array)
}

fn category<'py>(py: Python<'py>, names: &[String]) -> PyResult<Bound<'py, PyDict>> {
    let entry = PyDict::new(py);
    entry.set_item("names", PyList::new(py, names.iter().map(String::as_str))?)?;
    entry.set_item("colors", py.None())?;
    entry.set_item("mapping", PyDict::new(py))?;
    entry.set_item("other", py.None())?;
    Ok(entry)
}

fn reject_inf(name: &str, values: &[f64], valid: impl Fn(usize) -> bool) -> PyResult<()> {
    if let Some(i) = (0..values.len()).find(|&i| valid(i) && values[i].is_infinite()) {
        return Err(invalid(format!(
            "column {name:?} must be finite: row {i} is {}",
            values[i]
        )));
    }
    Ok(())
}

fn label(what: &str, value: &Option<String>) -> PyResult<()> {
    if value.as_ref().is_some_and(|s| s.trim().is_empty()) {
        return Err(invalid(format!(
            "{what} is blank; leave it unset when the object has none"
        )));
    }
    Ok(())
}

struct Exported {
    values: Py<PyAny>,
    unit: Option<String>,
    category: Option<Py<PyDict>>,
}

fn export_column(
    py: Python,
    name: &str,
    column: &Column,
    owner: &Bound<PyAny>,
) -> PyResult<Exported> {
    if name.is_empty() {
        return Err(invalid("a column name is empty"));
    }
    if column.unit.as_ref().is_some_and(|u| u.trim().is_empty()) {
        return Err(invalid(format!("column {name:?} has a blank unit")));
    }
    let n = column.len();
    let masked = column.valid.is_some();
    let (values, category) = match &column.data {
        ColumnData::Float64(v) => {
            reject_inf(name, v, |i| column.is_valid(i))?;
            let values = if masked {
                let out: Vec<f64> = (0..n)
                    .map(|i| if column.is_valid(i) { v[i] } else { f64::NAN })
                    .collect();
                owned(PyArray1::from_vec(py, out).into_any())?
            } else {
                view1(v, owner)?
            };
            (values, None)
        }
        ColumnData::Int64(v) => {
            if masked {
                let out: Vec<f64> = (0..n)
                    .map(|i| {
                        if column.is_valid(i) {
                            v[i] as f64
                        } else {
                            f64::NAN
                        }
                    })
                    .collect();
                (owned(PyArray1::from_vec(py, out).into_any())?, None)
            } else {
                (view1(v, owner)?, None)
            }
        }
        ColumnData::Bool(v) => {
            if masked {
                let codes: Vec<i32> = (0..n)
                    .map(|i| {
                        if !column.is_valid(i) {
                            -1
                        } else if v[i] {
                            1
                        } else {
                            0
                        }
                    })
                    .collect();
                let names = ["false".to_string(), "true".to_string()];
                (
                    owned(PyArray1::from_vec(py, codes).into_any())?,
                    Some(category(py, &names)?),
                )
            } else {
                (view1(v, owner)?, None)
            }
        }
        ColumnData::Utf8(v) => {
            if masked {
                let mut names = Vec::new();
                let mut index: HashMap<&str, i32> = HashMap::new();
                let mut codes = vec![-1; n];
                for i in 0..n {
                    if !column.is_valid(i) {
                        continue;
                    }
                    if v[i].is_empty() {
                        return Err(invalid(format!(
                            "column {name:?} has an empty string; a category name is non-empty"
                        )));
                    }
                    let code = if let Some(&code) = index.get(v[i].as_str()) {
                        code
                    } else {
                        let code = i32::try_from(names.len()).map_err(|_| {
                            invalid(format!("column {name:?} has too many category names"))
                        })?;
                        index.insert(v[i].as_str(), code);
                        names.push(v[i].clone());
                        code
                    };
                    codes[i] = code;
                }
                if names.is_empty() {
                    return Err(invalid(format!(
                        "column {name:?} has only missing values, so it has no category name"
                    )));
                }
                (
                    owned(PyArray1::from_vec(py, codes).into_any())?,
                    Some(category(py, &names)?),
                )
            } else {
                (unicode(py, v)?, None)
            }
        }
    };
    Ok(Exported {
        values: values.unbind(),
        unit: column.unit.clone(),
        category: category.map(|entry| entry.unbind()),
    })
}

fn export_columns<'py>(
    py: Python<'py>,
    table: &Table,
    owner: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, PyDict>, Bound<'py, PyDict>, Bound<'py, PyDict>)> {
    let columns = PyDict::new(py);
    let units = PyDict::new(py);
    let categories = PyDict::new(py);
    for (name, column) in table.names().iter().zip(table.columns()) {
        let exported = export_column(py, name, column, owner)?;
        columns.set_item(name, exported.values)?;
        if let Some(unit) = exported.unit {
            units.set_item(name, unit)?;
        }
        if let Some(entry) = exported.category {
            categories.set_item(name, entry)?;
        }
    }
    Ok((columns, units, categories))
}

fn common<'py>(
    py: Python<'py>,
    kind: &str,
    columns: Bound<'py, PyDict>,
    units: Bound<'py, PyDict>,
    categories: Bound<'py, PyDict>,
    crs: &Option<String>,
    length_unit: &Option<String>,
) -> PyResult<Bound<'py, PyDict>> {
    label("crs", crs)?;
    label("length_unit", length_unit)?;
    let d = PyDict::new(py);
    d.set_item("protocol", 1)?;
    d.set_item("kind", kind)?;
    d.set_item("columns", columns)?;
    d.set_item("units", units)?;
    d.set_item("categories", categories)?;
    d.set_item("crs", crs.as_deref())?;
    d.set_item("length_unit", length_unit.as_deref())?;
    Ok(d)
}

fn empty_dict<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    Ok(PyDict::new(py))
}

#[pymethods]
impl PyTable {
    /// A table of `n` rows and no columns.
    #[staticmethod]
    fn _with_rows(n: usize) -> PyTable {
        Table::new(n).into()
    }

    /// Protocol dict for this table.
    ///
    /// Returns
    /// -------
    /// dict
    ///     ``kind`` is ``"table"`` and ``rows`` is the row count. Each column
    ///     has shape ``(rows,)``. ``crs`` and ``length_unit`` are None: a
    ///     table does not store those labels.
    ///
    /// Notes
    /// -----
    /// An unmasked float64, int64 or bool column is a view of the stored
    /// buffer. A column with a validity mask is copied. A float null becomes
    /// NaN. An int column with a null widens to float64, and a null becomes
    /// NaN; an integer above 2**53 does not survive that widening exactly.
    /// A bool column with a null becomes int32 codes 0 and 1, with -1 for the
    /// null, and categories ``["false", "true"]``. Text with a null becomes
    /// int32 codes plus categories, in first-seen order. Text without a null
    /// is copied into a NumPy unicode array.
    ///
    /// Examples
    /// --------
    /// >>> t = wf.Table({"au": [1.0, None], "rock": ["ox", "fr"]})
    /// >>> d = t.__minetypes__()
    /// >>> d["kind"], d["rows"], d["columns"]["au"].tolist()
    /// ('table', 2, [1.0, nan])
    fn __minetypes__<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let py = this.py();
        let borrowed = this.borrow();
        let table = &borrowed.inner;
        let (columns, units, categories) = export_columns(py, table, this.as_any())?;
        let d = common(py, "table", columns, units, categories, &None, &None)?;
        d.set_item("rows", table.rows())?;
        Ok(d)
    }
}

#[pymethods]
impl PyMesh {
    /// Protocol dict for this mesh.
    ///
    /// Returns
    /// -------
    /// dict
    ///     ``kind`` is ``"mesh"``. ``coords`` is float64 ``(V, 3)`` and
    ///     ``triangles`` is uint32 ``(F, 3)``, both views of the stored
    ///     arrays. ``vertex_columns`` has length ``V`` and ``face_columns``
    ///     length ``F``. ``columns``, ``units`` and ``categories`` are empty.
    ///     ``crs`` and ``length_unit`` are the stored labels, or None.
    ///
    /// Notes
    /// -----
    /// Column copies follow ``Table.__minetypes__``. Triangles stored as
    /// uint32 stay a view. ``V`` is at most 2**32 - 1.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [1, 1, 1])
    /// >>> d = m.__minetypes__()
    /// >>> d["kind"], d["triangles"].dtype.name
    /// ('mesh', 'uint32')
    fn __minetypes__<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let py = this.py();
        let borrowed = this.borrow();
        let mesh = &borrowed.inner;
        let n = mesh.n_vertices();
        if n > u32::MAX as usize {
            return Err(invalid(format!(
                "coords has {n} rows; a mesh allows at most {} rows",
                u32::MAX
            )));
        }
        if mesh.vertices().iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("coords must be finite"));
        }
        let owner = this.as_any();
        let (columns, units, categories) = (empty_dict(py)?, empty_dict(py)?, empty_dict(py)?);
        let d = common(
            py,
            "mesh",
            columns,
            units,
            categories,
            &mesh.crs,
            &mesh.length_unit,
        )?;
        d.set_item("coords", view_rows::<f64, 3>(mesh.vertices(), owner)?)?;
        d.set_item("triangles", view_rows::<u32, 3>(mesh.triangles(), owner)?)?;
        let (vertex_columns, vertex_units, vertex_categories) =
            export_columns(py, mesh.vertex_attributes(), owner)?;
        let (face_columns, face_units, face_categories) =
            export_columns(py, mesh.face_attributes(), owner)?;
        d.set_item("vertex_columns", vertex_columns)?;
        d.set_item("face_columns", face_columns)?;
        d.set_item("vertex_units", vertex_units)?;
        d.set_item("face_units", face_units)?;
        d.set_item("vertex_categories", vertex_categories)?;
        d.set_item("face_categories", face_categories)?;
        Ok(d)
    }
}

#[pymethods]
impl PyPolylines {
    /// Protocol dict for these polylines.
    ///
    /// Returns
    /// -------
    /// dict
    ///     ``kind`` is ``"polylines"``. ``coords`` is float64 ``(V, 3)``.
    ///     ``part_offsets`` is uint32 ``(P + 1,)`` and ``feature_offsets``
    ///     is uint32 ``(M + 1,)``. ``closed`` is bool ``(P,)``. Those four
    ///     arrays are views of the stored buffers. ``columns`` holds the
    ///     feature attributes, one row per feature. ``vertex_columns`` has
    ///     length ``V``.
    ///
    /// Notes
    /// -----
    /// A closed part has at least 3 vertices and does not repeat its first
    /// vertex in storage. An open part has at least 2 vertices. One vertex
    /// is a ``PointSet``, and this method refuses that part.
    /// Column copies follow ``Table.__minetypes__``. ``V`` is at most
    /// 2**31 - 1.
    ///
    /// Examples
    /// --------
    /// >>> ring = wf.Polylines([[[0, 0, 0], [1, 0, 0], [0, 1, 0]]], closed=True)
    /// >>> ring.__minetypes__()["part_offsets"].tolist()
    /// [0, 3]
    fn __minetypes__<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let py = this.py();
        let borrowed = this.borrow();
        let lines = &borrowed.inner;
        let n = lines.n_vertices();
        if n > i32::MAX as usize {
            return Err(invalid(format!(
                "coords has {n} rows; polylines allow at most {} rows",
                i32::MAX
            )));
        }
        if lines.vertices().iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("coords must be finite"));
        }
        let parts = lines.part_offsets();
        let closed = lines.closed();
        for p in 0..lines.n_parts() {
            let count = (parts[p + 1] - parts[p]) as usize;
            let ring = closed[p];
            let need = if ring { 3 } else { 2 };
            if count < need {
                let noun = if count == 1 { "vertex" } else { "vertices" };
                let which = if ring { "closed" } else { "open" };
                return Err(invalid(format!(
                    "part {p} has {count} {noun}; {which} parts need at least {need}. A single point belongs in a PointSet"
                )));
            }
        }
        let owner = this.as_any();
        let (columns, units, categories) = export_columns(py, lines.attributes(), owner)?;
        let d = common(
            py,
            "polylines",
            columns,
            units,
            categories,
            &lines.crs,
            &lines.length_unit,
        )?;
        d.set_item("coords", view_rows::<f64, 3>(lines.vertices(), owner)?)?;
        d.set_item("part_offsets", view1(lines.part_offsets(), owner)?)?;
        d.set_item("feature_offsets", view1(lines.feature_offsets(), owner)?)?;
        d.set_item("closed", view1(lines.closed(), owner)?)?;
        let (vertex_columns, vertex_units, vertex_categories) =
            export_columns(py, lines.vertex_attributes(), owner)?;
        d.set_item("vertex_columns", vertex_columns)?;
        d.set_item("vertex_units", vertex_units)?;
        d.set_item("vertex_categories", vertex_categories)?;
        Ok(d)
    }
}

#[pymethods]
impl PyPointSet {
    /// Protocol dict for this point set.
    ///
    /// Returns
    /// -------
    /// dict
    ///     ``kind`` is ``"point_set"``. ``coords`` is a float64 ``(N, 3)``
    ///     view of the stored points. ``columns`` has length ``N``.
    ///
    /// Notes
    /// -----
    /// Column copies follow ``Table.__minetypes__``. A single polyline vertex
    /// is this kind, not an open part.
    ///
    /// Examples
    /// --------
    /// >>> p = wf.PointSet([[0, 0, 0], [1, 2, 3]], {"au": [0.5, 1.5]})
    /// >>> p.__minetypes__()["columns"]["au"].tolist()
    /// [0.5, 1.5]
    fn __minetypes__<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let py = this.py();
        let borrowed = this.borrow();
        let points = &borrowed.inner;
        if points.coords().iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("coords must be finite"));
        }
        let owner = this.as_any();
        let (columns, units, categories) = export_columns(py, points.attributes(), owner)?;
        let d = common(
            py,
            "point_set",
            columns,
            units,
            categories,
            &points.crs,
            &points.length_unit,
        )?;
        d.set_item("coords", view_rows::<f64, 3>(points.coords(), owner)?)?;
        Ok(d)
    }
}
