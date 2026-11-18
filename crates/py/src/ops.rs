use drusa_core::{Column, Table};
use drusa_ops::{HoleFill, RepairAction, RepairOptions, RepairReport, ThicknessMode};
use numpy::{PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::containers::{PyMesh, PyTable, array_n, array1, floats, plan_arg, points_arg, vec3};
use crate::errors::{invalid, to_py};
use crate::validate::{nullable_i64, nullable_utf8, summary_dict};

/// What ``Mesh.repair`` changed, step by step, and how input vertices and
/// faces map to the output.
#[pyclass(name = "RepairReport", module = "drusa", frozen)]
pub(crate) struct PyRepairReport {
    pub(crate) inner: RepairReport,
}

const ACTIONS: [&str; 5] = [
    "vertex_merged",
    "vertex_removed",
    "face_removed",
    "face_flipped",
    "face_added",
];

fn action_name(a: &RepairAction) -> &'static str {
    ACTIONS[match a {
        RepairAction::VertexMerged { .. } => 0,
        RepairAction::VertexRemoved { .. } => 1,
        RepairAction::FaceRemoved { .. } => 2,
        RepairAction::FaceFlipped { .. } => 3,
        RepairAction::FaceAdded { .. } => 4,
    }]
}

fn index_map<'py>(py: Python<'py>, map: &[Option<u32>]) -> Bound<'py, PyArray1<i64>> {
    PyArray1::from_vec(py, map.iter().map(|m| m.map_or(-1, i64::from)).collect())
}

#[pymethods]
impl PyRepairReport {
    /// One row per change in step order (weld, degenerate, duplicate, unused,
    /// orient, fill): ``action``, ``vertex``, ``face``, ``other`` (the vertex
    /// merged into, or the hole of an added face) and ``reason`` (the problem
    /// kind of a removed face).
    #[getter]
    fn actions(&self) -> PyResult<PyTable> {
        let a = &self.inner.actions;
        let (mut vertex, mut face, mut other, mut reason) = (vec![], vec![], vec![], vec![]);
        for act in a {
            let (v, f, o, r) = match act {
                RepairAction::VertexMerged { vertex, into } => {
                    (Some(*vertex), None, Some(*into), None)
                }
                RepairAction::VertexRemoved { vertex } => (Some(*vertex), None, None, None),
                RepairAction::FaceRemoved { face, reason } => {
                    (None, Some(*face), None, Some(reason.name().to_string()))
                }
                RepairAction::FaceFlipped { face } => (None, Some(*face), None, None),
                RepairAction::FaceAdded { face, hole } => (None, Some(*face), Some(*hole), None),
            };
            vertex.push(v.map(i64::from));
            face.push(f.map(i64::from));
            other.push(o.map(i64::from));
            reason.push(r);
        }
        let table = Table::new(a.len())
            .with_column(
                "action",
                Column::utf8(a.iter().map(|x| action_name(x).to_string()).collect()),
            )
            .and_then(|t| t.with_column("vertex", nullable_i64(vertex)))
            .and_then(|t| t.with_column("face", nullable_i64(face)))
            .and_then(|t| t.with_column("other", nullable_i64(other)))
            .and_then(|t| t.with_column("reason", nullable_utf8(reason)))
            .map_err(to_py)?;
        Ok(table.into())
    }

    /// Output vertex of each input vertex, -1 when removed.
    #[getter]
    fn vertex_map<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        index_map(py, &self.inner.vertex_map)
    }

    /// Output face of each input face, -1 when removed.
    #[getter]
    fn face_map<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        index_map(py, &self.inner.face_map)
    }

    /// Validation summary of the input.
    #[getter]
    fn before<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        summary_dict(py, &self.inner.before)
    }

    /// Validation summary of the output.
    #[getter]
    fn after<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        summary_dict(py, &self.inner.after)
    }

    /// Boundary loops left open.
    #[getter]
    fn holes_left(&self) -> usize {
        self.inner.holes_left
    }

    /// Number of actions of each kind.
    #[getter]
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for name in ACTIONS {
            let n = self
                .inner
                .actions
                .iter()
                .filter(|a| action_name(a) == name)
                .count();
            d.set_item(name, n)?;
        }
        Ok(d)
    }

    fn __len__(&self) -> usize {
        self.inner.actions.len()
    }

    fn __repr__(&self) -> String {
        let counts: Vec<String> = ACTIONS
            .iter()
            .map(|name| {
                let n = self
                    .inner
                    .actions
                    .iter()
                    .filter(|a| action_name(a) == *name)
                    .count();
                format!("{name}={n}")
            })
            .collect();
        format!(
            "RepairReport({}, holes_left={})",
            counts.join(", "),
            self.inner.holes_left
        )
    }
}

fn hole_fill(name: Option<&str>) -> PyResult<Option<HoleFill>> {
    match name {
        None => Ok(None),
        Some("fan") => Ok(Some(HoleFill::Fan)),
        Some("min_area") => Ok(Some(HoleFill::MinArea)),
        Some(other) => Err(invalid(format!(
            "fill_holes must be \"fan\", \"min_area\" or None, not {other:?}"
        ))),
    }
}

fn matrix4(obj: &Bound<PyAny>) -> PyResult<[[f64; 4]; 4]> {
    let bad = || invalid("matrix must be a (4, 4) array of numbers");
    let a: PyReadonlyArray2<f64> = crate::containers::numpy(obj.py())?
        .call_method1("ascontiguousarray", (obj, "float64"))
        .and_then(|a| Ok(a.extract()?))
        .map_err(|_| bad())?;
    let a = a.as_array();
    if a.dim() != (4, 4) {
        return Err(bad());
    }
    Ok(std::array::from_fn(|i| std::array::from_fn(|j| a[[i, j]])))
}

#[pymethods]
impl PyMesh {
    /// Total triangle area.
    #[getter]
    fn area(&self) -> f64 {
        drusa_ops::area(&self.inner)
    }

    /// Signed enclosed volume, positive for an outward solid.
    ///
    /// Raises
    /// ------
    /// InvalidMesh
    ///     When the mesh is not closed and consistently wound.
    #[getter]
    fn volume(&self) -> PyResult<f64> {
        drusa_ops::volume(&self.inner).map_err(to_py)
    }

    /// Volume centroid, under the same rule as ``volume``.
    #[getter]
    fn centroid<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(array1(
            py,
            drusa_ops::centroid(&self.inner).map_err(to_py)?.to_vec(),
        ))
    }

    /// Area of each triangle.
    #[getter]
    fn face_areas<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        array1(py, drusa_ops::face_areas(&self.inner))
    }

    /// ``(m, 3)`` unit normal of each triangle, zero for degenerate ones.
    #[getter]
    fn face_normals<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        array_n(py, &drusa_ops::face_normals(&self.inner))
    }

    /// ``(m, 3)`` centroid of each triangle.
    #[getter]
    fn face_centroids<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        array_n(py, &drusa_ops::face_centroids(&self.inner))
    }

    /// ``(n, 3)`` angle-weighted unit normal of each vertex, zero for isolated ones.
    #[getter]
    fn vertex_normals<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        array_n(py, &drusa_ops::vertex_normals(&self.inner))
    }

    /// ``(m, 2)`` dip (0 to 90 degrees) and dip direction (0 to 360 degrees
    /// clockwise from +y) of each triangle; NaN dip direction when horizontal.
    #[getter]
    fn dip<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let d: Vec<[f64; 2]> = drusa_ops::dip(&self.inner)
            .into_iter()
            .map(|(a, b)| [a, b])
            .collect();
        array_n(py, &d)
    }

    /// A copy with the affine transform `matrix` applied; faces flip when its
    /// determinant is negative, so a solid stays outward.
    ///
    /// Parameters
    /// ----------
    /// matrix : array_like
    ///     ``(4, 4)`` row-major matrix applied to column vectors.
    fn transform(&self, matrix: &Bound<PyAny>) -> PyResult<PyMesh> {
        let m = matrix4(matrix)?;
        Ok(drusa_ops::transform(&self.inner, m).map_err(to_py)?.into())
    }

    /// A copy moved by `offset`.
    fn translate(&self, offset: &Bound<PyAny>) -> PyResult<PyMesh> {
        let d = vec3(offset, "offset")?;
        Ok(drusa_ops::translate(&self.inner, d).map_err(to_py)?.into())
    }

    /// A copy with every triangle wound the other way.
    fn flip(&self) -> PyMesh {
        drusa_ops::flip(&self.inner).into()
    }

    /// A copy with coordinates rescaled to `unit`.
    ///
    /// Raises
    /// ------
    /// InvalidInput
    ///     When the mesh declares no length unit, or `unit` is unknown.
    fn to_length_unit(&self, unit: &str) -> PyResult<PyMesh> {
        Ok(drusa_ops::to_length_unit(&self.inner, unit)
            .map_err(to_py)?
            .into())
    }

    /// A repaired copy and the report of every change made.
    ///
    /// Steps run in order: weld, remove degenerate faces, remove duplicate
    /// faces, remove unused vertices, orient, fill holes.
    /// Self-intersections are never repaired.
    ///
    /// Parameters
    /// ----------
    /// weld : float or None
    ///     Merge vertices within this distance (transitively); 0 merges equal
    ///     coordinates, None skips welding.
    /// remove_degenerate, remove_duplicate_faces, remove_unused_vertices : bool
    ///     Which removals run.
    /// orient : bool
    ///     Make winding consistent and outward by nesting parity.
    /// fill_holes : {"fan", "min_area"} or None
    ///     Close boundary loops.
    ///
    /// Returns
    /// -------
    /// tuple of (Mesh, RepairReport)
    #[pyo3(signature = (*, weld=Some(0.0), remove_degenerate=true, remove_duplicate_faces=true, remove_unused_vertices=true, orient=true, fill_holes=None))]
    #[allow(clippy::too_many_arguments)]
    fn repair(
        &self,
        py: Python,
        weld: Option<f64>,
        remove_degenerate: bool,
        remove_duplicate_faces: bool,
        remove_unused_vertices: bool,
        orient: bool,
        fill_holes: Option<&str>,
    ) -> PyResult<(PyMesh, PyRepairReport)> {
        let opts = RepairOptions {
            weld,
            remove_degenerate,
            remove_duplicate_faces,
            remove_unused_vertices,
            orient,
            fill_holes: hole_fill(fill_holes)?,
        };
        let (mesh, report) = py
            .detach(|| drusa_ops::repair(&self.inner, &opts))
            .map_err(to_py)?;
        Ok((mesh.into(), PyRepairReport { inner: report }))
    }

    /// One mesh per edge-connected component, in order of first face.
    fn components(&self) -> PyResult<Vec<PyMesh>> {
        Ok(drusa_ops::components(&self.inner)
            .map_err(to_py)?
            .into_iter()
            .map(PyMesh::from)
            .collect())
    }

    /// Edge-connected component of each triangle, numbered by first face.
    fn component_labels<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        PyArray1::from_vec(
            py,
            drusa_ops::component_labels(&self.inner)
                .into_iter()
                .map(i64::from)
                .collect(),
        )
    }

    /// Solid volume between consecutive elevations, by exact clamped-column
    /// integration.
    ///
    /// Parameters
    /// ----------
    /// edges : array_like
    ///     Strictly increasing elevations, at least 2.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     ``len(edges) - 1`` volumes; value ``i`` covers ``edges[i] <= z < edges[i + 1]``.
    fn volume_by_elevation<'py>(
        &self,
        py: Python<'py>,
        edges: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let e = floats(edges, "edges")?;
        let v = py
            .detach(|| drusa_ops::volume_by_elevation(&self.inner, &e))
            .map_err(to_py)?;
        Ok(array1(py, v))
    }

    /// Highest surface elevation above each point, NaN where no face covers it.
    ///
    /// Parameters
    /// ----------
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 2)`` or ``(n, 3)``; z is ignored.
    fn surface_z<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let xy = plan_arg(points, "points")?;
        let z = py.detach(|| drusa_ops::surface_z(&self.inner, &xy));
        Ok(array1(py, z))
    }

    /// Point elevation minus ``surface_z``: positive above the surface.
    fn vertical_distance<'py>(
        &self,
        py: Python<'py>,
        points: &Bound<PyAny>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let p = points_arg(points, "points", false)?;
        let d = py.detach(|| drusa_ops::vertical_distance(&self.inner, &p));
        Ok(array1(py, d))
    }

    /// Thickness between this surface (top) and `bottom` at each point.
    ///
    /// Parameters
    /// ----------
    /// bottom : Mesh
    ///     Lower surface.
    /// points : array_like, PointSet or Mesh
    ///     ``(n, 2)`` or ``(n, 3)``; z is ignored.
    /// mode : {"vertical", "normal"}
    ///     ``vertical`` is the elevation difference; ``normal`` the distance
    ///     from the top point to the closest point of `bottom`, signed like
    ///     ``vertical``.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Negative where the surfaces cross, NaN where either is missing.
    #[pyo3(signature = (bottom, points, *, mode="vertical"))]
    fn thickness<'py>(
        &self,
        py: Python<'py>,
        bottom: PyRef<PyMesh>,
        points: &Bound<PyAny>,
        mode: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let mode = match mode {
            "vertical" => ThicknessMode::Vertical,
            "normal" => ThicknessMode::Normal,
            _ => {
                return Err(invalid(format!(
                    "mode must be \"vertical\" or \"normal\", not {mode:?}"
                )));
            }
        };
        let xy = plan_arg(points, "points")?;
        let bottom = &bottom.inner;
        let t = py
            .detach(|| drusa_ops::thickness(&self.inner, bottom, &xy, mode))
            .map_err(to_py)?;
        Ok(array1(py, t))
    }
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<PyRepairReport>()?;
    Ok(())
}
