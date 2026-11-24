use pyo3::prelude::*;

mod boolean;
mod build;
mod containers;
mod errors;
mod io;
mod ops;
mod protocol;
mod query;
mod validate;

#[pymodule]
fn _wireframes(m: &Bound<PyModule>) -> PyResult<()> {
    let version = env!("CARGO_PKG_VERSION")
        .replace("-alpha.", "a")
        .replace("-beta.", "b")
        .replace("-rc.", "rc");
    m.add("__version__", version)?;
    containers::register(m)?;
    validate::register(m)?;
    io::register(m)?;
    ops::register(m)?;
    boolean::register(m)?;
    query::register(m)?;
    build::register(m)?;
    Ok(())
}
