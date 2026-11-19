use pyo3::prelude::*;
use pyo3::types::PyType;
use wireframes_core::Error;

use crate::containers::PyMesh;
use crate::validate::PyMeshReport;

/// An instance of `wireframes.errors.<kind>`, which subclasses both `WireframesError`
/// and a builtin exception.
pub(crate) fn error(kind: &str, message: impl ToString) -> PyErr {
    Python::attach(|py| {
        match py
            .import("wireframes.errors")
            .and_then(|m| m.getattr(kind))
            .and_then(|c| Ok(c.cast_into::<PyType>()?))
        {
            Ok(class) => PyErr::from_type(class, message.to_string()),
            Err(e) => e,
        }
    })
}

pub(crate) fn invalid(message: impl ToString) -> PyErr {
    error("InvalidInput", message)
}

/// Maps a wireframes error to its Python class; `InvalidMesh` carries `.report`,
/// `KernelError` `.report` and `.mesh` (or `None`).
pub(crate) fn to_py(e: Error) -> PyErr {
    let message = e.to_string();
    match e {
        Error::InvalidMesh { report, .. } => Python::attach(|py| {
            let err = error("InvalidMesh", &message);
            let attached = Py::new(py, PyMeshReport { inner: *report })
                .and_then(|r| err.value(py).setattr("report", r));
            match attached {
                Ok(()) => err,
                Err(e) => e,
            }
        }),
        Error::Kernel { report, mesh, .. } => Python::attach(|py| {
            let err = error("KernelError", &message);
            let value = err.value(py);
            let attached = (|| {
                let report = report
                    .map(|r| Py::new(py, PyMeshReport { inner: *r }))
                    .transpose()?;
                let mesh = mesh.map(|m| Py::new(py, PyMesh::from(*m))).transpose()?;
                value.setattr("report", report)?;
                value.setattr("mesh", mesh)
            })();
            match attached {
                Ok(()) => err,
                Err(e) => e,
            }
        }),
        Error::Format { .. } | Error::Io(_) => error("FileError", message),
        Error::MissingColumn { .. } => error("MissingColumn", message),
        Error::InvalidInput(_)
        | Error::Unsupported(_)
        | Error::CrsMismatch { .. }
        | Error::UnitMismatch { .. } => invalid(message),
    }
}
