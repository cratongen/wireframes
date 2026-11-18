#![allow(unused_variables)]

use drusa_core::{Mesh, PointSet, Polylines, Table};
use drusa_io::UnstructuredGrid;
use pyo3::prelude::*;

use crate::errors::not_implemented;

#[pyclass(name = "Table", module = "drusa", frozen)]
pub(crate) struct PyTable {
    pub(crate) inner: Table,
}

#[pyclass(name = "Mesh", module = "drusa", frozen)]
pub(crate) struct PyMesh {
    pub(crate) inner: Mesh,
}

#[pyclass(name = "Polylines", module = "drusa", frozen)]
pub(crate) struct PyPolylines {
    pub(crate) inner: Polylines,
}

#[pyclass(name = "PointSet", module = "drusa", frozen)]
pub(crate) struct PyPointSet {
    pub(crate) inner: PointSet,
}

#[pyclass(name = "UnstructuredGrid", module = "drusa", frozen)]
pub(crate) struct PyUnstructuredGrid {
    pub(crate) inner: UnstructuredGrid,
}

#[pyfunction]
fn concat(meshes: Vec<PyRef<PyMesh>>) -> PyResult<PyMesh> {
    Err(not_implemented("concat"))
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
