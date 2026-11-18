#![allow(unused_variables)]

use pyo3::prelude::*;

use crate::containers::{PyMesh, PyPolylines};
use crate::errors::not_implemented;

macro_rules! boolean_stub {
    ($name:ident) => {
        #[pyfunction]
        #[pyo3(signature = (a, b, *, repair=false))]
        fn $name(a: PyRef<PyMesh>, b: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
            Err(not_implemented(stringify!($name)))
        }
    };
}

boolean_stub!(union);
boolean_stub!(intersection);
boolean_stub!(difference);
boolean_stub!(symmetric_difference);

#[pyfunction]
#[pyo3(signature = (base, design, *, boundary=None))]
fn cut_fill(
    base: PyRef<PyMesh>,
    design: PyRef<PyMesh>,
    boundary: Option<PyRef<PyPolylines>>,
) -> PyResult<Py<PyAny>> {
    Err(not_implemented("cut_fill"))
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(union, m)?)?;
    m.add_function(wrap_pyfunction!(intersection, m)?)?;
    m.add_function(wrap_pyfunction!(difference, m)?)?;
    m.add_function(wrap_pyfunction!(symmetric_difference, m)?)?;
    m.add_function(wrap_pyfunction!(cut_fill, m)?)?;
    Ok(())
}
