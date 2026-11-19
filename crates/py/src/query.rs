use numpy::{PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use wireframes_core::{Column, Table};
use wireframes_query::{Blocks, Hit, HitKind, Trace};

use crate::containers::{PyMesh, PyTable, array_n, array1, floats, numpy, points_arg, rows3};
use crate::errors::{invalid, to_py};
use crate::validate::nullable_i64;

pub(crate) fn finite_points(
    obj: &Bound<PyAny>,
    what: &str,
    allow2: bool,
) -> PyResult<Vec<[f64; 3]>> {
    let p = points_arg(obj, what, allow2)?;
    if let Some(i) = p.iter().position(|q| q.iter().any(|x| !x.is_finite())) {
        return Err(invalid(format!(
            "{what} must be finite: row {i} is {:?}",
            p[i]
        )));
    }
    Ok(p)
}

fn kind_name(k: HitKind) -> &'static str {
    match k {
        HitKind::Interior => "interior",
        HitKind::Edge => "edge",
        HitKind::Vertex => "vertex",
        HitKind::Coplanar => "coplanar",
    }
}

fn hits_table(first: &str, hits: Vec<Vec<Hit>>) -> PyResult<PyTable> {
    let flat: Vec<(usize, Hit)> = hits
        .into_iter()
        .enumerate()
        .flat_map(|(i, hs)| hs.into_iter().map(move |h| (i, h)))
        .collect();
    let n = flat.len();
    let entering_valid: Vec<bool> = flat.iter().map(|(_, h)| h.entering.is_some()).collect();
    let entering = Column::bool(
        flat.iter()
            .map(|(_, h)| h.entering.unwrap_or(false))
            .collect(),
    )
    .with_valid(entering_valid)
    .map_err(to_py)?;
    let f = |g: &dyn Fn(&Hit) -> f64| Column::f64(flat.iter().map(|(_, h)| g(h)).collect());
    let table = Table::new(n)
        .with_column(
            first,
            Column::i64(flat.iter().map(|(i, _)| *i as i64).collect()),
        )
        .and_then(|t| t.with_column("t", f(&|h| h.t)))
        .and_then(|t| {
            t.with_column(
                "face",
                Column::i64(flat.iter().map(|(_, h)| i64::from(h.face)).collect()),
            )
        })
        .and_then(|t| t.with_column("x", f(&|h| h.point[0])))
        .and_then(|t| t.with_column("y", f(&|h| h.point[1])))
        .and_then(|t| t.with_column("z", f(&|h| h.point[2])))
        .and_then(|t| {
            t.with_column(
                "kind",
                Column::utf8(
                    flat.iter()
                        .map(|(_, h)| kind_name(h.kind).to_string())
                        .collect(),
                ),
            )
        })
        .and_then(|t| t.with_column("entering", entering))
        .map_err(to_py)?;
    Ok(table.into())
}

fn rotation_arg(obj: &Bound<PyAny>) -> PyResult<[[f64; 3]; 3]> {
    let bad = || invalid("rotation must be a (3, 3) array of numbers");
    let a: PyReadonlyArray2<f64> = numpy(obj.py())?
        .call_method1("ascontiguousarray", (obj, "float64"))
        .and_then(|a| Ok(a.extract()?))
        .map_err(|_| bad())?;
    let a = a.as_array();
    if a.dim() != (3, 3) {
        return Err(bad());
    }
    Ok(std::array::from_fn(|i| std::array::from_fn(|j| a[[i, j]])))
}

#[pymethods]
impl PyMesh {
    /// Where each point lies relative to the solid, decided exactly.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 3)`` points.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray of int8
    ///     -1 inside, 0 on the boundary, 1 outside; the sign agrees with
    ///     ``signed_distance``.
    ///
    /// Raises
    /// ------
    /// InvalidMesh
    ///     When the mesh is not a solid.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.locate([[1, 1, 1], [2, 1, 1], [5, 5, 5]]).tolist()
    /// [-1, 0, 1]
    fn locate<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<i8>>> {
        let p = finite_points(points, "points", false)?;
        let index = self.solid_index(py)?;
        let loc = py.detach(|| index.locate_all(&p));
        Ok(PyArray1::from_vec(
            py,
            loc.into_iter().map(|l| l as i8).collect(),
        ))
    }

    /// Whether each point lies in the solid.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 3)`` points.
    /// boundary : bool
    ///     Whether a point on the surface counts as inside (closed set).
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.contains([[1, 1, 1], [2, 1, 1], [5, 5, 5]]).tolist()
    /// [True, True, False]
    /// >>> m.contains([[1, 1, 1], [2, 1, 1], [5, 5, 5]], boundary=False).tolist()
    /// [True, False, False]
    #[pyo3(signature = (points, *, boundary=true))]
    fn contains<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
        boundary: bool,
    ) -> PyResult<Bound<'py, PyArray1<bool>>> {
        let p = finite_points(points, "points", false)?;
        let index = self.solid_index(py)?;
        let inside = py.detach(|| index.contains_all(&p, boundary));
        Ok(PyArray1::from_vec(py, inside))
    }

    /// Distance from each point to the surface, negative inside the solid.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 3)`` points.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.signed_distance([[1, 1, 1], [3, 1, 1]]).tolist()
    /// [-1.0, 1.0]
    fn signed_distance<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let p = finite_points(points, "points", false)?;
        let index = self.solid_index(py)?;
        let d = py.detach(|| index.signed_distance(&p));
        Ok(array1(py, d))
    }

    /// Unsigned distance from each point to the closest point of any face.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 3)`` points.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.distance([[1, 1, 1], [3, 1, 1]]).tolist()
    /// [1.0, 1.0]
    fn distance<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let p = finite_points(points, "points", false)?;
        let index = self.surface_index(py);
        let d = py.detach(|| index.distance(&p));
        Ok(array1(py, d))
    }

    /// The closest point of the mesh to each point.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 3)`` points.
    ///
    /// Returns
    /// -------
    /// tuple of (numpy.ndarray, numpy.ndarray, numpy.ndarray)
    ///     ``(n, 3)`` closest points, the face holding each (smaller index on
    ///     ties) and the distances.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> points, faces, distances = m.closest_point([[3, 1, 1]])
    /// >>> points.tolist(), distances.tolist()
    /// ([[2.0, 1.0, 1.0]], [1.0])
    #[allow(clippy::type_complexity)]
    fn closest_point<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<(
        Bound<'py, PyArray2<f64>>,
        Bound<'py, PyArray1<i64>>,
        Bound<'py, PyArray1<f64>>,
    )> {
        let p = finite_points(points, "points", false)?;
        let index = self.surface_index(py);
        let c = py.detach(|| index.closest(&p));
        let pts: Vec<[f64; 3]> = c.iter().map(|c| c.point).collect();
        Ok((
            array_n(py, &pts),
            PyArray1::from_vec(py, c.iter().map(|c| i64::from(c.face)).collect()),
            array1(py, c.iter().map(|c| c.distance).collect()),
        ))
    }

    /// Every face each ray meets, with exact crossing decisions.
    ///
    /// Parameters
    /// ----------
    /// origins, directions : array_like
    ///     ``(n, 3)`` ray starts and nonzero directions.
    ///
    /// Returns
    /// -------
    /// Table
    ///     One row per hit with ``t >= 0``, sorted by ray then ``t``:
    ///     ``ray``, ``t`` (in units of the direction), ``face``, ``x``, ``y``,
    ///     ``z``, ``kind`` (interior, edge, vertex, coplanar) and ``entering``
    ///     (null for open meshes and coplanar hits).
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> hits = m.ray_hits([[1, 1, -1]], [[0, 0, 1]])
    /// >>> hits["t"].tolist(), hits["entering"].tolist()
    /// ([1.0, 3.0], [True, False])
    fn ray_hits(
        &self,
        py: Python,
        origins: &Bound<PyAny>,
        directions: &Bound<PyAny>,
    ) -> PyResult<PyTable> {
        let o = rows3(origins, "origins", false)?;
        let d = rows3(directions, "directions", false)?;
        let index = self.surface_index(py);
        let hits = py.detach(|| index.ray_hits(&o, &d)).map_err(to_py)?;
        hits_table("ray", hits)
    }

    /// Every face each segment meets: ``ray_hits`` with ``0 <= t <= 1``.
    ///
    /// Parameters
    /// ----------
    /// starts, ends : array_like
    ///     ``(n, 3)`` segment end points.
    ///
    /// Returns
    /// -------
    /// Table
    ///     As ``ray_hits``, with the column ``segment`` in place of ``ray``.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.segment_hits([[1, 1, -1]], [[1, 1, 1]])["t"].tolist()
    /// [0.5]
    fn segment_hits(
        &self,
        py: Python,
        starts: &Bound<PyAny>,
        ends: &Bound<PyAny>,
    ) -> PyResult<PyTable> {
        let a = rows3(starts, "starts", false)?;
        let b = rows3(ends, "ends", false)?;
        let index = self.surface_index(py);
        let hits = py.detach(|| index.segment_hits(&a, &b)).map_err(to_py)?;
        hits_table("segment", hits)
    }

    /// Inside and outside runs of a drillhole trace through the solid.
    ///
    /// Parameters
    /// ----------
    /// trace : array_like
    ///     ``(n, 3)`` desurveyed trace vertices.
    /// depths : array_like
    ///     Measured depth of each vertex, non-decreasing.
    ///
    /// Returns
    /// -------
    /// Table
    ///     ``from``, ``to``, ``inside`` (alternating), ``entry_face`` and
    ///     ``exit_face`` (int64, null where the trace meets no face),
    ///     covering ``depths[0]`` to ``depths[-1]``.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> runs = m.trace_intervals([[1, 1, -1], [1, 1, 3]], [0.0, 4.0])
    /// >>> runs["from"].tolist(), runs["to"].tolist(), runs["inside"].tolist()
    /// ([0.0, 1.0, 3.0], [1.0, 3.0, 4.0], [False, True, False])
    fn trace_intervals(
        &self,
        py: Python,
        trace: &Bound<PyAny>,
        depths: &Bound<PyAny>,
    ) -> PyResult<PyTable> {
        let t = Trace {
            xyz: rows3(trace, "trace", false)?,
            depth: floats(depths, "depths")?,
        };
        let index = self.solid_index(py)?;
        let iv = py.detach(|| index.trace_intervals(&t)).map_err(to_py)?;
        let table = Table::new(iv.len())
            .with_column("from", Column::f64(iv.iter().map(|i| i.from).collect()))
            .and_then(|tb| tb.with_column("to", Column::f64(iv.iter().map(|i| i.to).collect())))
            .and_then(|tb| {
                tb.with_column(
                    "inside",
                    Column::bool(iv.iter().map(|i| i.inside).collect()),
                )
            })
            .and_then(|tb| {
                tb.with_column(
                    "entry_face",
                    nullable_i64(iv.iter().map(|i| i.entry_face.map(i64::from)).collect()),
                )
            })
            .and_then(|tb| {
                tb.with_column(
                    "exit_face",
                    nullable_i64(iv.iter().map(|i| i.exit_face.map(i64::from)).collect()),
                )
            })
            .map_err(to_py)?;
        Ok(table.into())
    }

    /// Fraction of each block's volume inside the solid, by exact
    /// clamped-column integration (no sampling).
    ///
    /// Parameters
    /// ----------
    /// centroids : array_like, PointSet or Mesh
    ///     ``(n, 3)`` block centres.
    /// size : array_like
    ///     Block size along its three axes, positive.
    /// rotation : array_like, optional
    ///     ``(3, 3)`` matrix whose columns are the block axes in world
    ///     coordinates (orthonormal, determinant +1); None is axis-aligned.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Values in ``[0, 1]``.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.block_fraction([[2, 1, 1], [1, 1, 1]], [1, 1, 1]).tolist()
    /// [0.5, 1.0]
    #[pyo3(signature = (centroids, size, *, rotation=None))]
    fn block_fraction<'py>(
        &self,
        py: Python<'py>,
        centroids: &Bound<PyAny>,
        size: &Bound<PyAny>,
        rotation: Option<&Bound<PyAny>>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let size = floats(size, "size")?;
        let size: [f64; 3] = size
            .try_into()
            .map_err(|_| invalid("size must have 3 values"))?;
        let blocks = Blocks {
            centroids: points_arg(centroids, "centroids", false)?,
            size,
            rotation: rotation.map(rotation_arg).transpose()?,
        };
        let index = self.solid_index(py)?;
        let f = py.detach(|| index.block_fraction(&blocks)).map_err(to_py)?;
        Ok(array1(py, f))
    }
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    let _ = m;
    Ok(())
}
