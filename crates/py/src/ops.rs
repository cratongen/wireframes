use drusa_ops::RepairReport;
use pyo3::prelude::*;

#[pyclass(name = "RepairReport", module = "drusa", frozen)]
pub(crate) struct PyRepairReport {
    pub(crate) inner: RepairReport,
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<PyRepairReport>()?;
    Ok(())
}
