use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use wireframes_boolean::{BooleanOptions, Op};

use crate::containers::{PyMesh, PyPolylines, floats, vec3};
use crate::errors::to_py;
use crate::ops::PyRepairReport;

fn run(py: Python, a: &PyMesh, b: &PyMesh, op: Op, repair: bool) -> PyResult<Py<PyAny>> {
    let opts = BooleanOptions { repair };
    let result = py
        .detach(|| wireframes_boolean::boolean(&a.inner, &b.inner, op, &opts))
        .map_err(to_py)?;
    let mesh = PyMesh::from(result.mesh);
    if !repair {
        return mesh.into_py_any(py);
    }
    let [ra, rb] = result.repairs;
    let report = |r: Option<wireframes_ops::RepairReport>| PyRepairReport {
        inner: r.expect("repair=True repairs both inputs"),
    };
    (mesh, (report(ra), report(rb))).into_py_any(py)
}

macro_rules! boolean_fn {
    ($name:ident, $op:expr, $doc:literal, $examples:literal) => {
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
        ///     Parts that touch along an edge or at a vertex are separate
        ///     shells with their own vertices, so the result is valid input
        ///     to another boolean.
        ///
        /// Raises
        /// ------
        /// InvalidMesh
        ///     When an input is not a solid; ``report`` holds its validation.
        /// InvalidInput
        ///     When one shell of the result would touch itself along an
        ///     edge or at a vertex, which no valid solid represents.
        ///
        /// Examples
        /// --------
        #[doc = $examples]
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
    "The solid inside `a` or `b` (also ``a | b``), by an exact mesh arrangement.",
    ">>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])\n>>> wf.union(a, b).volume, (a | b).volume\n(15.0, 15.0)"
);
boolean_fn!(
    intersection,
    Op::Intersection,
    "The solid inside both `a` and `b` (also ``a & b``), by an exact mesh arrangement.",
    ">>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])\n>>> wf.intersection(a, b).volume, (a & b).volume\n(1.0, 1.0)"
);
boolean_fn!(
    difference,
    Op::Difference,
    "The solid inside `a` and outside `b` (also ``a - b``), by an exact mesh arrangement.",
    ">>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])\n>>> wf.difference(a, b).volume, (a - b).volume\n(7.0, 7.0)"
);
boolean_fn!(
    symmetric_difference,
    Op::SymmetricDifference,
    "The solid inside exactly one of `a` and `b` (also ``a ^ b``), by an exact mesh arrangement.",
    ">>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])\n>>> wf.symmetric_difference(a, b).volume, (a ^ b).volume\n(14.0, 14.0)"
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
///
/// Examples
/// --------
/// >>> base = wf.dtm([[0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]])
/// >>> design = base.translate([0, 0, -2])
/// >>> volumes = wf.cut_fill(base, design)
/// >>> volumes["cut"], volumes["fill"], volumes["area"]
/// (200.0, 0.0, 100.0)
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
        .detach(|| wireframes_boolean::cut_fill(base, design, boundary))
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
    /// ``wireframes.union(self, other, repair=repair)``.
    ///
    /// Parameters
    /// ----------
    /// other : Mesh
    ///     The second solid, with equal or undeclared ``crs`` and ``length_unit``.
    /// repair : bool
    ///     Repair both inputs with ``Mesh.repair()`` defaults first and return the
    ///     reports.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.union(b).volume
    /// 15.0
    #[pyo3(signature = (other, *, repair=false))]
    fn union(&self, py: Python, other: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Union, repair)
    }

    /// ``wireframes.intersection(self, other, repair=repair)``.
    ///
    /// Parameters
    /// ----------
    /// other : Mesh
    ///     The second solid, with equal or undeclared ``crs`` and ``length_unit``.
    /// repair : bool
    ///     Repair both inputs with ``Mesh.repair()`` defaults first and return the
    ///     reports.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.intersection(b).volume
    /// 1.0
    #[pyo3(signature = (other, *, repair=false))]
    fn intersection(&self, py: Python, other: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Intersection, repair)
    }

    /// ``wireframes.difference(self, other, repair=repair)``.
    ///
    /// Parameters
    /// ----------
    /// other : Mesh
    ///     The second solid, with equal or undeclared ``crs`` and ``length_unit``.
    /// repair : bool
    ///     Repair both inputs with ``Mesh.repair()`` defaults first and return the
    ///     reports.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.difference(b).volume
    /// 7.0
    #[pyo3(signature = (other, *, repair=false))]
    fn difference(&self, py: Python, other: PyRef<PyMesh>, repair: bool) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Difference, repair)
    }

    /// ``wireframes.symmetric_difference(self, other, repair=repair)``.
    ///
    /// Parameters
    /// ----------
    /// other : Mesh
    ///     The second solid, with equal or undeclared ``crs`` and ``length_unit``.
    /// repair : bool
    ///     Repair both inputs with ``Mesh.repair()`` defaults first and return the
    ///     reports.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.symmetric_difference(b).volume
    /// 14.0
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
    ///
    /// Raises
    /// ------
    /// InvalidInput
    ///     When the surface ends inside the solid, or when a part of the solid
    ///     lies beyond the surface's extent (no line along x, y or z through
    ///     it meets the surface), so it is neither above nor below.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> cut = wf.dtm([[-1, -1, 0.5], [3, -1, 0.5], [3, 3, 0.5], [-1, 3, 0.5]])
    /// >>> below, above = m.clip_by_surface(cut)
    /// >>> below.volume, above.volume
    /// (2.0, 6.0)
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
            .detach(|| wireframes_boolean::clip_by_surface(&self.inner, surface, &opts))
            .map_err(to_py)?;
        Ok((r.below.into(), r.above.into()))
    }

    /// The mesh split by the plane through `origin` with `normal`, decided exactly.
    ///
    /// Closed meshes are capped; open meshes are trimmed. Points on the plane
    /// go to the positive side.
    ///
    /// Parameters
    /// ----------
    /// origin : array_like
    ///     ``(3,)`` point on the plane.
    /// normal : array_like
    ///     ``(3,)`` nonzero plane normal; it need not be unit length.
    ///
    /// Returns
    /// -------
    /// tuple of (Mesh, Mesh)
    ///     ``(negative, positive)``, positive where ``(p - origin) . normal > 0``.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> negative, positive = m.split_by_plane([0, 0, 0.5], [0, 0, 1])
    /// >>> negative.volume, positive.volume
    /// (2.0, 6.0)
    fn split_by_plane(
        &self,
        py: Python,
        origin: &Bound<PyAny>,
        normal: &Bound<PyAny>,
    ) -> PyResult<(PyMesh, PyMesh)> {
        let (o, n) = (vec3(origin, "origin")?, vec3(normal, "normal")?);
        let (neg, pos) = py
            .detach(|| wireframes_boolean::split_by_plane(&self.inner, o, n))
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
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.slice_by_plane([0, 0, 1], [0, 0, 1], offsets=[-0.5, 0.5]).length().tolist()
    /// [8.0, 8.0]
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
            .detach(|| wireframes_boolean::slice_by_plane(&self.inner, o, n, &offsets))
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
