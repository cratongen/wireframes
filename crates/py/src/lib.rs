use pyo3::prelude::*;

mod boolean;
mod build;
mod containers;
mod errors;
mod io;
mod ops;
mod query;
mod validate;

#[pymodule]
fn _wireframes(m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    containers::register(m)?;
    validate::register(m)?;
    io::register(m)?;
    ops::register(m)?;
    boolean::register(m)?;
    query::register(m)?;
    build::register(m)?;
    Ok(())
}
