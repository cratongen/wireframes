use pyo3::prelude::*;
use pyo3::types::PyDict;
use wireframes_core::{Column, MeshReport, MeshSummary, Table};

use crate::containers::{PyMesh, PyTable};
use crate::errors::to_py;

/// An int64 column with nulls where `values` has None.
pub(crate) fn nullable_i64(values: Vec<Option<i64>>) -> Column {
    let valid: Vec<bool> = values.iter().map(Option::is_some).collect();
    let column = Column::i64(values.into_iter().map(|v| v.unwrap_or(0)).collect());
    if valid.iter().all(|&v| v) {
        column
    } else {
        column.with_valid(valid).expect("same length")
    }
}

/// A utf8 column with nulls where `values` has None.
pub(crate) fn nullable_utf8(values: Vec<Option<String>>) -> Column {
    let valid: Vec<bool> = values.iter().map(Option::is_some).collect();
    let column = Column::utf8(values.into_iter().map(Option::unwrap_or_default).collect());
    if valid.iter().all(|&v| v) {
        column
    } else {
        column.with_valid(valid).expect("same length")
    }
}

/// Problem counts in a stable order, ``is_closed`` and ``is_solid`` as booleans.
pub(crate) fn summary_dict<'py>(py: Python<'py>, s: &MeshSummary) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    for (name, n) in s.to_pairs() {
        if name == "is_closed" {
            d.set_item(name, n != 0)?;
        } else {
            d.set_item(name, n)?;
        }
    }
    d.set_item("is_solid", s.is_solid())?;
    Ok(d)
}

/// What ``Mesh.validate`` found: problem counts and every problem with its
/// face, vertex and edge.
#[pyclass(name = "MeshReport", module = "wireframes", frozen)]
pub(crate) struct PyMeshReport {
    pub(crate) inner: MeshReport,
}

#[pymethods]
impl PyMeshReport {
    /// Problem counts by kind, plus ``shells``, the contact counts
    /// ``touching_shell_edges`` and ``touching_shell_vertices`` (closed shells
    /// that touch each other or themselves along coincident edges or at
    /// coincident vertices without overlapping; information, not defects),
    /// ``is_closed`` and ``is_solid``.
    #[getter]
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        summary_dict(py, &self.inner.summary)
    }

    /// One row per problem, sorted by kind then face: ``kind`` (str), ``face``,
    /// ``vertex``, ``edge`` (local corner ``e`` to ``e + 1``) and ``other``
    /// (earlier duplicate, other face of an edge or of a pair), int64 with
    /// nulls.
    #[getter]
    fn problems(&self) -> PyResult<PyTable> {
        let p = &self.inner.problems;
        let int = |f: &dyn Fn(&wireframes_core::MeshProblem) -> Option<i64>| {
            nullable_i64(p.iter().map(f).collect())
        };
        let table = Table::new(p.len())
            .with_column(
                "kind",
                Column::utf8(p.iter().map(|q| q.kind.name().to_string()).collect()),
            )
            .and_then(|t| t.with_column("face", int(&|q| q.face.map(i64::from))))
            .and_then(|t| t.with_column("vertex", int(&|q| q.vertex.map(i64::from))))
            .and_then(|t| t.with_column("edge", int(&|q| q.edge.map(i64::from))))
            .and_then(|t| t.with_column("other", int(&|q| q.other.map(i64::from))))
            .map_err(to_py)?;
        Ok(table.into())
    }

    /// Closed, manifold, consistently and outward wound, without
    /// self-intersections. Shells may touch each other along edges or at
    /// vertices.
    #[getter]
    fn is_solid(&self) -> bool {
        self.inner.summary.is_solid()
    }

    fn __len__(&self) -> usize {
        self.inner.problems.len()
    }

    fn __repr__(&self) -> String {
        let found: Vec<String> = self
            .inner
            .summary
            .to_pairs()
            .into_iter()
            .filter(|(n, v)| *v != 0 && *n != "is_closed" && *n != "shells")
            .map(|(n, v)| format!("{n}={v}"))
            .collect();
        format!(
            "MeshReport(is_solid={}, shells={}, {})",
            self.inner.summary.is_solid(),
            self.inner.summary.shells,
            if found.is_empty() {
                "no problems".to_string()
            } else {
                found.join(", ")
            }
        )
    }
}

#[pymethods]
impl PyMesh {
    /// Every topological and geometric problem of the mesh, decided exactly.
    ///
    /// The mesh is a solid when it is closed, manifold, consistently and
    /// outward wound (a shell inside another winds inward), and free of
    /// degenerate and duplicate faces, duplicate vertices, self-intersections
    /// and coplanar overlaps. Closed shells may touch each other, and a shell
    /// may touch itself, along exactly coincident edges or at exactly
    /// coincident vertices; the summary counts those contacts as information.
    /// Shells that share part of a face, or meet on part of an edge, overlap
    /// and are reported.
    ///
    /// Parameters
    /// ----------
    /// tolerance : float
    ///     Distance within which two vertices count as duplicates; 0 means
    ///     equal coordinates. Exact copies of a vertex where closed shells
    ///     touch are not duplicates. Every other check is exact and ignores
    ///     `tolerance`.
    ///
    /// Returns
    /// -------
    /// MeshReport
    ///     Counts and every problem.
    ///
    /// Raises
    /// ------
    /// InvalidInput
    ///     When `tolerance` is negative or not finite.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.validate()
    /// MeshReport(is_solid=true, shells=1, no problems)
    /// >>> wf.datasets.topography(n=20).validate().summary["boundary_edges"] > 0
    /// True
    #[pyo3(signature = (*, tolerance=0.0))]
    fn validate(&self, py: Python, tolerance: f64) -> PyResult<PyMeshReport> {
        let report = py
            .detach(|| self.inner.validate(tolerance))
            .map_err(to_py)?;
        Ok(PyMeshReport { inner: report })
    }

    /// Whether every edge of the non-degenerate faces joins exactly two faces,
    /// or more that pair into wedges of solid touching there.
    #[getter]
    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    /// Whether ``validate()`` finds a solid: closed, manifold, outward and
    /// without self-intersections; its shells may touch each other along
    /// edges or at vertices.
    #[getter]
    fn is_solid(&self, py: Python) -> PyResult<bool> {
        Ok(py
            .detach(|| self.inner.validate(0.0))
            .map_err(to_py)?
            .summary
            .is_solid())
    }
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<PyMeshReport>()?;
    Ok(())
}
