use drusa_boolean::{BooleanOptions, Op};
use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::containers::{PyMesh, PyPolylines, floats, vec3};
use crate::errors::to_py;
use crate::ops::PyRepairReport;

fn run(py: Python, a: &PyMesh, b: &PyMesh, op: Op, repair: bool) -> PyResult<Py<PyAny>> {
    let opts = BooleanOptions { repair };
    let result = py
        .detach(|| drusa_boolean::boolean(&a.inner, &b.inner, op, &opts))
        .map_err(to_py)?;
    let mesh = PyMesh::from(result.mesh);
    if !repair {
        return mesh.into_py_any(py);
    }
    let [ra, rb] = result.repairs;
    let report = |r: Option<drusa_ops::RepairReport>| PyRepairReport {
        inner: r.expect("repair=True repairs both inputs"),
    };
    (mesh, (report(ra), report(rb))).into_py_any(py)
}

macro_rules! boolean_fn {
    ($name:ident, $op:expr, $doc:literal) => {
        #[doc = $doc]
        ///
        /// Parameters
        /// ----------
        /// a, b : Mesh
        ///     Solids (``is_solid``) with equal or undeclared ``crs`` and
        ///     ``length_unit``.
        /// repair : bool
        ///     Repair each input with ``Mesh.repair()`` defaults first and
        ///     return the reports.
        ///
        /// Returns
        /// -------
        /// Mesh or tuple of (Mesh, (RepairReport, RepairReport))
        ///     A solid; face attributes follow the faces they came from.
        ///
        /// Raises
        /// ------
        /// InvalidMesh
        ///     When an input is not a solid; ``report`` holds its validation.
        /// InvalidInput
        ///     When the solids touch along an edge or vertex only, so the
        ///     result would not be manifold.
        #[pyfunction]
        #[pyo3(signature = (a, b, *, repair=false))]
        fn $name(
            py: Python,
            a: PyRef<PyMesh>,
            b: PyRef<PyMesh>,
            repair: bool,
        ) -> PyResult<Py<PyAny>> {
            run(py, &a, &b, $op, repair)
        }
    };
}

boolean_fn!(
    union,
    Op::Union,
    "The solid inside `a` or `b` (also ``a | b``), by an exact mesh arrangement."
);
boolean_fn!(
    intersection,
    Op::Intersection,
    "The solid inside both `a` and `b` (also ``a & b``), by an exact mesh arrangement."
);
boolean_fn!(
    difference,
    Op::Difference,
    "The solid inside `a` and outside `b` (also ``a - b``), by an exact mesh arrangement."
);
boolean_fn!(
    symmetric_difference,
    Op::SymmetricDifference,
    "The solid inside exactly one of `a` and `b` (also ``a ^ b``), by an exact mesh arrangement."
);

/// Cut and fill volumes of `design` against `base`, by an exact plan overlay.
///
/// Both surfaces must be height fields. The domain is where both surfaces
/// exist in plan, inside `boundary` when given.
///
/// Parameters
/// ----------
/// base : Mesh
///     Reference surface, such as the current topography.
/// design : Mesh
///     Compared surface, such as a pit design or a later survey.
/// boundary : Polylines, optional
///     Closed rings limiting the domain (even-odd).
///
/// Returns
/// -------
/// dict of str to float
///     ``cut`` (volume where design is below base), ``fill`` (above),
///     ``net`` (fill minus cut), ``cut_area``, ``fill_area`` and ``area``
///     (plan areas).
#[pyfunction]
#[pyo3(signature = (base, design, *, boundary=None))]
fn cut_fill<'py>(
    py: Python<'py>,
    base: PyRef<PyMesh>,
    design: PyRef<PyMesh>,
    boundary: Option<PyRef<PyPolylines>>,
) -> PyResult<Bound<'py, PyDict>> {
    let boundary = boundary.as_ref().map(|b| &b.inner);
    let (base, design) = (&base.inner, &design.inner);
    let r = py
        .detach(|| drusa_boolean::cut_fill(base, design, boundary))
        .map_err(to_py)?;
    let d = PyDict::new(py);
    for (k, v) in [
        ("cut", r.cut),
        ("fill", r.fill),
        ("net", r.net),
        ("cut_area", r.cut_area),
        ("fill_area", r.fill_area),
        ("area", r.area),
    ] {
        d.set_item(k, v)?;
    }
    Ok(d)
}

#[pymethods]
impl PyMesh {
    /// ``drusa.union(self, other, repair=repair)``.
    #[pyo3(signature = (other, *, repair=false))]
    fn union(&self, py: Python, other: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Union, repair)
    }

    /// ``drusa.intersection(self, other, repair=repair)``.
    #[pyo3(signature = (other, *, repair=false))]
    fn intersection(&self, py: Python, other: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Intersection, repair)
    }

    /// ``drusa.difference(self, other, repair=repair)``.
    #[pyo3(signature = (other, *, repair=false))]
    fn difference(&self, py: Python, other: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Difference, repair)
    }

    /// ``drusa.symmetric_difference(self, other, repair=repair)``.
    #[pyo3(signature = (other, *, repair=false))]
    fn symmetric_difference(
        &self,
        py: Python,
        other: PyRef<PyMesh>,
        repair: bool,
    ) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::SymmetricDifference, repair)
    }

    fn __or__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Union, false)
    }

    fn __and__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Intersection, false)
    }

    fn __sub__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Difference, false)
    }

    fn __xor__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::SymmetricDifference, false)
    }

    /// The solid cut in two by an open `surface`, each part capped.
    ///
    /// Parameters
    /// ----------
    /// surface : Mesh
    ///     Open surface whose boundary lies outside the solid; ``above`` is
    ///     the side its normals point to.
    /// repair : bool
    ///     Repair both inputs first.
    ///
    /// Returns
    /// -------
    /// tuple of (Mesh, Mesh)
    ///     ``(below, above)``, both solids; cap faces carry the surface's
    ///     face attributes. A surface that misses the solid gives the whole
    ///     solid on one side and an empty mesh on the other.
    #[pyo3(signature = (surface, *, repair=false))]
    fn clip_by_surface(
        &self,
        py: Python,
        surface: PyRef<PyMesh>,
        repair: bool,
    ) -> PyResult<(PyMesh, PyMesh)> {
        let opts = BooleanOptions { repair };
        let surface = &surface.inner;
        let r = py
            .detach(|| drusa_boolean::clip_by_surface(&self.inner, surface, &opts))
            .map_err(to_py)?;
        Ok((r.below.into(), r.above.into()))
    }

    /// The mesh split by the plane through `origin` with `normal`, decided exactly.
    ///
    /// Closed meshes are capped; open meshes are trimmed. Points on the plane
    /// go to the positive side.
    ///
    /// Returns
    /// -------
    /// tuple of (Mesh, Mesh)
    ///     ``(negative, positive)``, positive where ``(p - origin) . normal > 0``.
    fn split_by_plane(
        &self,
        py: Python,
        origin: &Bound<PyAny>,
        normal: &Bound<PyAny>,
    ) -> PyResult<(PyMesh, PyMesh)> {
        let (o, n) = (vec3(origin, "origin")?, vec3(normal, "normal")?);
        let (neg, pos) = py
            .detach(|| drusa_boolean::split_by_plane(&self.inner, o, n))
            .map_err(to_py)?;
        Ok((neg.into(), pos.into()))
    }

    /// Section lines where planes parallel to the one through `origin` with
    /// `normal` cut the mesh.
    ///
    /// Parameters
    /// ----------
    /// origin, normal : array_like
    ///     The plane at offset 0.
    /// offsets : array_like, optional
    ///     Distances along the unit normal; None means ``[0.0]``.
    ///
    /// Returns
    /// -------
    /// Polylines
    ///     One feature per line, with the feature column ``offset`` and the
    ///     face columns of the faces each line crosses. Closed meshes give
    ///     rings counter-clockwise seen from the normal side.
    #[pyo3(signature = (origin, normal, *, offsets=None))]
    fn slice_by_plane(
        &self,
        py: Python,
        origin: &Bound<PyAny>,
        normal: &Bound<PyAny>,
        offsets: Option<&Bound<PyAny>>,
    ) -> PyResult<PyPolylines> {
        let (o, n) = (vec3(origin, "origin")?, vec3(normal, "normal")?);
        let offsets = match offsets {
            Some(x) => floats(x, "offsets")?,
            None => vec![],
        };
        let lines = py
            .detach(|| drusa_boolean::slice_by_plane(&self.inner, o, n, &offsets))
            .map_err(to_py)?;
        Ok(lines.into())
    }
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(union, m)?)?;
    m.add_function(wrap_pyfunction!(intersection, m)?)?;
    m.add_function(wrap_pyfunction!(difference, m)?)?;
    m.add_function(wrap_pyfunction!(symmetric_difference, m)?)?;
    m.add_function(wrap_pyfunction!(cut_fill, m)?)?;
    Ok(())
}
