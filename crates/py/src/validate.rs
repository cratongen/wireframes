use drusa_core::MeshReport;
use pyo3::prelude::*;

#[pyclass(name = "MeshReport", module = "drusa", frozen)]
pub(crate) struct PyMeshReport {
    pub(crate) inner: MeshReport,
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<PyMeshReport>()?;
    Ok(())
}
