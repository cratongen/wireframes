#![allow(unused_variables)]

use std::path::PathBuf;

use pyo3::prelude::*;

use crate::containers::{PyMesh, PyPointSet, PyPolylines};
use crate::errors::not_implemented;

/// `(mesh, polylines, points, grid, length_unit, warnings)`; `drusa.read` wraps it in a `Dataset`.
#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn _read(
    py: Python,
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<Py<PyAny>> {
    Err(not_implemented("read"))
}

#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn read_mesh(path: PathBuf, crs: Option<String>, length_unit: Option<String>) -> PyResult<PyMesh> {
    Err(not_implemented("read_mesh"))
}

#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn read_polylines(
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<PyPolylines> {
    Err(not_implemented("read_polylines"))
}

#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn read_points(
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<PyPointSet> {
    Err(not_implemented("read_points"))
}

#[pyfunction]
#[pyo3(signature = (path, data, *, ascii=false, dxf_entity="3dface", zlib=false, skip_unsupported_columns=false))]
fn write(
    path: PathBuf,
    data: &Bound<PyAny>,
    ascii: bool,
    dxf_entity: &str,
    zlib: bool,
    skip_unsupported_columns: bool,
) -> PyResult<()> {
    Err(not_implemented("write"))
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_read, m)?)?;
    m.add_function(wrap_pyfunction!(read_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(read_polylines, m)?)?;
    m.add_function(wrap_pyfunction!(read_points, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    Ok(())
}
