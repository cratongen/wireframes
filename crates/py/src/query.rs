use pyo3::prelude::*;

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    let _ = m;
    Ok(())
}
