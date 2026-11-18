#![allow(unused_variables)]

use pyo3::prelude::*;

use crate::containers::{PyMesh, PyPolylines};
use crate::errors::not_implemented;

#[pyfunction]
#[pyo3(signature = (points, *, breaklines=None, boundary=None, min_spacing=0.0, max_edge=None))]
fn dtm(
    points: &Bound<PyAny>,
    breaklines: Option<PyRef<PyPolylines>>,
    boundary: Option<PyRef<PyPolylines>>,
    min_spacing: f64,
    max_edge: Option<f64>,
) -> PyResult<PyMesh> {
    Err(not_implemented("dtm"))
}

#[pyfunction]
fn triangulate(polygons: PyRef<PyPolylines>) -> PyResult<PyMesh> {
    Err(not_implemented("triangulate"))
}

#[pyfunction]
#[pyo3(signature = (sections, *, cap_ends=true))]
fn loft(sections: PyRef<PyPolylines>, cap_ends: bool) -> PyResult<PyMesh> {
    Err(not_implemented("loft"))
}

#[pyfunction]
fn extrude(polygons: PyRef<PyPolylines>, vector: [f64; 3]) -> PyResult<PyMesh> {
    Err(not_implemented("extrude"))
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dtm, m)?)?;
    m.add_function(wrap_pyfunction!(triangulate, m)?)?;
    m.add_function(wrap_pyfunction!(loft, m)?)?;
    m.add_function(wrap_pyfunction!(extrude, m)?)?;
    Ok(())
}
