use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use wireframes_boolean::{BooleanOptions, Op, RoundingReport};

use crate::containers::{PyMesh, PyPolylines, floats, vec3};
use crate::errors::to_py;
use crate::ops::PyRepairReport;

/// What rounding the exact result of a boolean, clip or split to f64 did, and
/// the bounds it proves.
///
/// Every output vertex is an input vertex or the rounding of an exact vertex
/// of the result. Rounding moves a vertex by at most half an ulp per
/// coordinate; a vertex that ends anywhere else counts as moved, and no
/// vertex ends farther than ``move_bound`` from its exact position. Lengths
/// are in the mesh's length unit, rounded up; ``max_move`` and
/// ``move_bound`` are ``(x, y, z)``.
///
/// Attributes
/// ----------
/// rounded_vertices : int
///     Output vertices that are rounded intersection points.
/// moved_vertices : int
///     Exact vertices that ended at coordinates other than their own rounding.
/// moved_input_vertices : int
///     The input vertices among them.
/// max_move : tuple of float
///     Largest distance per axis from an exact vertex to where it ended.
/// move_bound : tuple of float
///     The bound on that distance: four ulps of each vertex's largest
///     coordinate.
/// merged_vertices, collapsed_edges, inserted_points, flipped_edges : int
///     Vertices merged with equal coordinates, vertices joined with a
///     neighbour within rounding, points inserted into faces they lie within
///     rounding of (including vertices of one input inserted into the other
///     where the two agree to rounding), and edges flipped to match the other
///     input there.
/// removed_shells, cancelled_faces : int
///     Shells of zero rounded volume removed, and faces removed in opposite
///     pairs.
/// volume_bound : float
///     Proven bound on ``|volume(result) - volume(exact result)|``.
/// steps, step_bound : int
///     Operations applied, and the bound on them proven before the first.
///
/// Examples
/// --------
/// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
/// >>> m, report = wf.union(a, b, rounding_report=True)
/// >>> report.moved_vertices, report.volume_bound < 1e-12, report.steps <= report.step_bound
/// (0, True, True)
#[pyclass(name = "RoundingReport", module = "wireframes", frozen)]
pub(crate) struct PyRoundingReport {
    pub(crate) inner: RoundingReport,
}

#[pymethods]
impl PyRoundingReport {
    #[getter]
    fn rounded_vertices(&self) -> usize {
        self.inner.rounded_vertices
    }
    #[getter]
    fn moved_vertices(&self) -> usize {
        self.inner.moved_vertices
    }
    #[getter]
    fn moved_input_vertices(&self) -> usize {
        self.inner.moved_input_vertices
    }
    #[getter]
    fn max_move(&self) -> (f64, f64, f64) {
        let m = self.inner.max_move;
        (m[0], m[1], m[2])
    }
    #[getter]
    fn move_bound(&self) -> (f64, f64, f64) {
        let m = self.inner.move_bound;
        (m[0], m[1], m[2])
    }
    #[getter]
    fn merged_vertices(&self) -> usize {
        self.inner.merged_vertices
    }
    #[getter]
    fn collapsed_edges(&self) -> usize {
        self.inner.collapsed_edges
    }
    #[getter]
    fn inserted_points(&self) -> usize {
        self.inner.inserted_points
    }
    #[getter]
    fn flipped_edges(&self) -> usize {
        self.inner.flipped_edges
    }
    #[getter]
    fn removed_shells(&self) -> usize {
        self.inner.removed_shells
    }
    #[getter]
    fn cancelled_faces(&self) -> usize {
        self.inner.cancelled_faces
    }
    #[getter]
    fn volume_bound(&self) -> f64 {
        self.inner.volume_bound
    }
    #[getter]
    fn steps(&self) -> usize {
        self.inner.steps
    }
    #[getter]
    fn step_bound(&self) -> usize {
        self.inner.step_bound
    }

    fn __repr__(&self) -> String {
        let r = &self.inner;
        format!(
            "RoundingReport(moved_vertices={}, moved_input_vertices={}, max_move={:?}, \
             volume_bound={:e}, steps={})",
            r.moved_vertices, r.moved_input_vertices, r.max_move, r.volume_bound, r.steps
        )
    }
}

fn run(
    py: Python,
    a: &PyMesh,
    b: &PyMesh,
    op: Op,
    repair: bool,
    rounding_report: bool,
) -> PyResult<Py<PyAny>> {
    let opts = BooleanOptions { repair };
    let result = py
        .detach(|| wireframes_boolean::boolean(&a.inner, &b.inner, op, &opts))
        .map_err(to_py)?;
    let mesh = PyMesh::from(result.mesh);
    let rounding = PyRoundingReport {
        inner: result.rounding,
    };
    let report = |r: Option<wireframes_ops::RepairReport>| PyRepairReport {
        inner: r.expect("repair=True repairs both inputs"),
    };
    let [ra, rb] = result.repairs;
    match (repair, rounding_report) {
        (false, false) => mesh.into_py_any(py),
        (false, true) => (mesh, rounding).into_py_any(py),
        (true, false) => (mesh, (report(ra), report(rb))).into_py_any(py),
        (true, true) => (mesh, (report(ra), report(rb)), rounding).into_py_any(py),
    }
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
        /// rounding_report : bool
        ///     Also return the ``RoundingReport``: how far rounding the exact
        ///     result to f64 moved any vertex, and the proven bound on the
        ///     volume change.
        ///
        /// Returns
        /// -------
        /// Mesh, or a tuple of the Mesh, then ``(RepairReport, RepairReport)`` when ``repair``, then ``RoundingReport`` when ``rounding_report``
        ///     A closed mesh; face attributes follow the faces they came from.
        ///     Parts that touch along an edge or at a vertex are separate
        ///     shells with their own vertices, so the result is valid input
        ///     to another boolean. The result is always a solid.
        ///
        /// Raises
        /// ------
        /// InvalidMesh
        ///     When an input is not a solid; ``report`` holds its validation.
        /// KernelError
        ///     When the rounding operations cannot certify the rounded result
        ///     as a solid (the message starts with ``uncertified``); ``report``
        ///     and ``mesh`` hold the closed rounded result. No case known to
        ///     the tests raises it.
        ///
        /// Examples
        /// --------
        #[doc = $examples]
        #[pyfunction]
        #[pyo3(signature = (a, b, *, repair=false, rounding_report=false))]
        fn $name(
            py: Python,
            a: PyRef<PyMesh>,
            b: PyRef<PyMesh>,
            repair: bool,
            rounding_report: bool,
        ) -> PyResult<Py<PyAny>> {
            run(py, &a, &b, $op, repair, rounding_report)
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
/// Both surfaces must be height fields: no two faces overlap in plan. The
/// domain is where both surfaces exist in plan, inside `boundary` when given.
/// The overlay is integrated exactly per cell, so net volumes over one domain
/// add up: ``net(a, c) == net(a, b) + net(b, c)`` to rounding.
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
/// Raises
/// ------
/// InvalidInput
///     When a surface is not a height field, naming the surface and a face,
///     or when ``boundary`` holds no closed ring.
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
    /// rounding_report : bool
    ///     Also return the ``RoundingReport``, last.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.union(b).volume
    /// 15.0
    #[pyo3(signature = (other, *, repair=false, rounding_report=false))]
    fn union(
        &self,
        py: Python,
        other: PyRef<PyMesh>,
        repair: bool,
        rounding_report: bool,
    ) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Union, repair, rounding_report)
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
    /// rounding_report : bool
    ///     Also return the ``RoundingReport``, last.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.intersection(b).volume
    /// 1.0
    #[pyo3(signature = (other, *, repair=false, rounding_report=false))]
    fn intersection(
        &self,
        py: Python,
        other: PyRef<PyMesh>,
        repair: bool,
        rounding_report: bool,
    ) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Intersection, repair, rounding_report)
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
    /// rounding_report : bool
    ///     Also return the ``RoundingReport``, last.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.difference(b).volume
    /// 7.0
    #[pyo3(signature = (other, *, repair=false, rounding_report=false))]
    fn difference(
        &self,
        py: Python,
        other: PyRef<PyMesh>,
        repair: bool,
        rounding_report: bool,
    ) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Difference, repair, rounding_report)
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
    /// rounding_report : bool
    ///     Also return the ``RoundingReport``, last.
    ///
    /// Examples
    /// --------
    /// >>> a, b = wf.datasets.box([0, 0, 0], [2, 2, 2]), wf.datasets.box([1, 1, 1], [3, 3, 3])
    /// >>> a.symmetric_difference(b).volume
    /// 14.0
    #[pyo3(signature = (other, *, repair=false, rounding_report=false))]
    fn symmetric_difference(
        &self,
        py: Python,
        other: PyRef<PyMesh>,
        repair: bool,
        rounding_report: bool,
    ) -> PyResult<Py<PyAny>> {
        run(
            py,
            self,
            &other,
            Op::SymmetricDifference,
            repair,
            rounding_report,
        )
    }

    fn __or__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Union, false, false)
    }

    fn __and__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Intersection, false, false)
    }

    fn __sub__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::Difference, false, false)
    }

    fn __xor__(&self, py: Python, other: PyRef<PyMesh>) -> PyResult<Py<PyAny>> {
        run(py, self, &other, Op::SymmetricDifference, false, false)
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
    /// rounding_report : bool
    ///     Also return ``(RoundingReport, RoundingReport)`` for ``below`` and
    ///     ``above``.
    ///
    /// Returns
    /// -------
    /// tuple of (Mesh, Mesh), or (Mesh, Mesh, (RoundingReport, RoundingReport))
    ///     ``(below, above)``, both closed; cap faces carry the surface's
    ///     face attributes. A surface that misses the solid gives the whole
    ///     solid on one side and an empty mesh on the other. Both are
    ///     solids.
    ///
    /// Raises
    /// ------
    /// InvalidMesh
    ///     When the solid is not a solid or the surface fails validation;
    ///     ``report`` holds the validation.
    /// InvalidInput
    ///     When the surface ends inside the solid, or when a part of the solid
    ///     lies beyond the surface's extent (no line along x, y or z through
    ///     it meets the surface), so it is neither above nor below.
    /// KernelError
    ///     When the rounding operations cannot certify a part as a solid
    ///     (``uncertified``); ``report`` and ``mesh`` hold that part.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> cut = wf.dtm([[-1, -1, 0.5], [3, -1, 0.5], [3, 3, 0.5], [-1, 3, 0.5]])
    /// >>> below, above = m.clip_by_surface(cut)
    /// >>> below.volume, above.volume
    /// (2.0, 6.0)
    #[pyo3(signature = (surface, *, repair=false, rounding_report=false))]
    fn clip_by_surface(
        &self,
        py: Python,
        surface: PyRef<PyMesh>,
        repair: bool,
        rounding_report: bool,
    ) -> PyResult<Py<PyAny>> {
        let opts = BooleanOptions { repair };
        let surface = &surface.inner;
        let r = py
            .detach(|| wireframes_boolean::clip_by_surface(&self.inner, surface, &opts))
            .map_err(to_py)?;
        let (below, above) = (PyMesh::from(r.below), PyMesh::from(r.above));
        if !rounding_report {
            return (below, above).into_py_any(py);
        }
        let [lo, hi] = r.rounding.map(|inner| PyRoundingReport { inner });
        (below, above, (lo, hi)).into_py_any(py)
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
    /// rounding_report : bool
    ///     Also return ``(RoundingReport, RoundingReport)`` for the two parts.
    ///
    /// Returns
    /// -------
    /// tuple of (Mesh, Mesh), or (Mesh, Mesh, (RoundingReport, RoundingReport))
    ///     ``(negative, positive)``, positive where ``(p - origin) . normal > 0``.
    ///     When the mesh is a solid, both parts are solids; where the mesh
    ///     touches the plane at a vertex or along an edge, such as a local
    ///     minimum of a closed DTM at a contour level, the part on that side
    ///     touches itself there.
    ///
    /// Raises
    /// ------
    /// KernelError
    ///     When the mesh is a solid and the rounding operations cannot
    ///     certify a part as a solid (``uncertified``); ``report`` and
    ///     ``mesh`` hold that part.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> negative, positive = m.split_by_plane([0, 0, 0.5], [0, 0, 1])
    /// >>> negative.volume, positive.volume
    /// (2.0, 6.0)
    #[pyo3(signature = (origin, normal, *, rounding_report=false))]
    fn split_by_plane(
        &self,
        py: Python,
        origin: &Bound<PyAny>,
        normal: &Bound<PyAny>,
        rounding_report: bool,
    ) -> PyResult<Py<PyAny>> {
        let (o, n) = (vec3(origin, "origin")?, vec3(normal, "normal")?);
        let (neg, pos, reports) = py
            .detach(|| wireframes_boolean::split_by_plane_report(&self.inner, o, n))
            .map_err(to_py)?;
        let (neg, pos) = (PyMesh::from(neg), PyMesh::from(pos));
        if !rounding_report {
            return (neg, pos).into_py_any(py);
        }
        let [lo, hi] = reports.map(|inner| PyRoundingReport { inner });
        (neg, pos, (lo, hi)).into_py_any(py)
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
    m.add_class::<PyRoundingReport>()?;
    m.add_function(wrap_pyfunction!(union, m)?)?;
    m.add_function(wrap_pyfunction!(intersection, m)?)?;
    m.add_function(wrap_pyfunction!(difference, m)?)?;
    m.add_function(wrap_pyfunction!(symmetric_difference, m)?)?;
    m.add_function(wrap_pyfunction!(cut_fill, m)?)?;
    Ok(())
}
