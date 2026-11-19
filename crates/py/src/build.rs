use pyo3::prelude::*;
use wireframes_build::{DecimateOptions, DtmOptions, LoftOptions};

use crate::containers::{PyMesh, PyPolylines, floats, labels_of, vec3};
use crate::errors::to_py;
use crate::query::finite_points;

/// Surface by constrained Delaunay triangulation in plan.
///
/// Parameters
/// ----------
/// points : array_like, PointSet or Mesh
///     ``(n, 3)`` survey points, at least 3 not collinear in plan.
/// breaklines : Polylines, optional
///     Strings enforced as chains of edges; where two cross, z comes from
///     the one inserted first.
/// boundary : Polylines, optional
///     Closed rings (even-odd); triangles outside are dropped.
/// min_spacing : float
///     Points closer than this in plan to an earlier point are dropped; 0
///     drops only exact duplicates.
/// max_edge : float, optional
///     Drop boundary triangles with a plan edge longer than this, repeatedly.
///
/// Returns
/// -------
/// Mesh
///     Faces wound counter-clockwise in plan (normals up).
#[pyfunction]
#[pyo3(signature = (points, *, breaklines=None, boundary=None, min_spacing=0.0, max_edge=None))]
fn dtm(
    py: Python,
    points: &Bound<PyAny>,
    breaklines: Option<PyRef<PyPolylines>>,
    boundary: Option<PyRef<PyPolylines>>,
    min_spacing: f64,
    max_edge: Option<f64>,
) -> PyResult<PyMesh> {
    let p = finite_points(points, "points", false)?;
    let opts = DtmOptions {
        breaklines: breaklines.as_ref().map(|b| &b.inner),
        boundary: boundary.as_ref().map(|b| &b.inner),
        min_spacing,
        max_edge,
    };
    let (crs, unit) = labels_of(points);
    Ok(py
        .detach(|| wireframes_build::dtm(&p, &opts))
        .map_err(to_py)?
        .with_crs(crs)
        .with_length_unit(unit)
        .into())
}

/// Each feature's rings triangulated in their best-fit plane, holes by even-odd.
///
/// Parameters
/// ----------
/// polygons : Polylines
///     Closed rings; the open parts of a feature are ignored.
///
/// Returns
/// -------
/// Mesh
///     The original vertices, face column ``feature``; normals follow the
///     outer ring's Newell normal.
#[pyfunction]
fn triangulate(py: Python, polygons: PyRef<PyPolylines>) -> PyResult<PyMesh> {
    let p = &polygons.inner;
    Ok(py
        .detach(|| wireframes_build::triangulate_polygons(p))
        .map_err(to_py)?
        .into())
}

/// Solid linking consecutive section rings, one ring per feature.
///
/// Parameters
/// ----------
/// sections : Polylines
///     Features in order, each with exactly one closed part.
/// cap_ends : bool
///     Close the first and last sections.
///
/// Returns
/// -------
/// Mesh
///     A validated solid with face column ``section`` (band index; caps are
///     -1 and -2).
#[pyfunction]
#[pyo3(signature = (sections, *, cap_ends=true))]
fn loft(py: Python, sections: PyRef<PyPolylines>, cap_ends: bool) -> PyResult<PyMesh> {
    let s = &sections.inner;
    let opts = LoftOptions { cap_ends };
    Ok(py
        .detach(|| wireframes_build::loft(s, &opts))
        .map_err(to_py)?
        .into())
}

/// Solid swept by each closed ring along `vector`.
///
/// Parameters
/// ----------
/// polygons : Polylines
///     Closed rings.
/// vector : array_like
///     Sweep direction and length, not parallel to the ring plane.
///
/// Returns
/// -------
/// Mesh
///     An outward solid.
#[pyfunction]
fn extrude(py: Python, polygons: PyRef<PyPolylines>, vector: &Bound<PyAny>) -> PyResult<PyMesh> {
    let v = vec3(vector, "vector")?;
    let p = &polygons.inner;
    Ok(py
        .detach(|| wireframes_build::extrude(p, v))
        .map_err(to_py)?
        .into())
}

#[pymethods]
impl PyMesh {
    /// The solid between this height-field surface and the plane at elevation `z`.
    ///
    /// Walls drop from the boundary loop to `z` and a flat base closes them.
    /// The surface must have one boundary loop, entirely above or entirely
    /// below `z`.
    fn close_to_elevation(&self, py: Python, z: f64) -> PyResult<PyMesh> {
        Ok(py
            .detach(|| wireframes_build::close_to_elevation(&self.inner, z))
            .map_err(to_py)?
            .into())
    }

    /// A copy with fewer triangles, by quadric edge collapse.
    ///
    /// Parameters
    /// ----------
    /// target_faces : int, optional
    ///     Stop at this many faces.
    /// max_error : float, optional
    ///     Stop when the next collapse moves the surface more than this.
    /// preserve_boundary : bool
    ///     Keep boundary vertices on their boundary lines.
    #[pyo3(signature = (*, target_faces=None, max_error=None, preserve_boundary=true))]
    fn decimate(
        &self,
        py: Python,
        target_faces: Option<usize>,
        max_error: Option<f64>,
        preserve_boundary: bool,
    ) -> PyResult<PyMesh> {
        let opts = DecimateOptions {
            target_faces,
            max_error,
            preserve_boundary,
        };
        Ok(py
            .detach(|| wireframes_build::decimate(&self.inner, &opts))
            .map_err(to_py)?
            .into())
    }

    /// A copy smoothed by Taubin's lambda and mu steps, boundary vertices fixed.
    ///
    /// Parameters
    /// ----------
    /// iterations : int
    ///     Number of lambda and mu step pairs.
    /// lambda_ : float
    ///     Shrinking step, positive.
    /// mu : float
    ///     Inflating step, below ``-lambda_``.
    #[pyo3(signature = (*, iterations=10, lambda_=0.5, mu=-0.53))]
    fn smooth(&self, py: Python, iterations: usize, lambda_: f64, mu: f64) -> PyResult<PyMesh> {
        Ok(py
            .detach(|| wireframes_build::smooth_taubin(&self.inner, iterations, lambda_, mu))
            .map_err(to_py)?
            .into())
    }

    /// A copy with every vertex moved by `distance` along its vertex normal,
    /// outward when positive.
    ///
    /// Raises
    /// ------
    /// InvalidMesh
    ///     When the result self-intersects.
    fn offset(&self, py: Python, distance: f64) -> PyResult<PyMesh> {
        Ok(py
            .detach(|| wireframes_build::offset(&self.inner, distance))
            .map_err(to_py)?
            .into())
    }

    /// A copy remeshed towards edges of length `edge_length` (isotropic).
    ///
    /// Face attributes come from the closest input face; vertex attributes
    /// are dropped.
    #[pyo3(signature = (edge_length, *, iterations=5))]
    fn remesh(&self, py: Python, edge_length: f64, iterations: usize) -> PyResult<PyMesh> {
        Ok(py
            .detach(|| wireframes_build::remesh(&self.inner, edge_length, iterations))
            .map_err(to_py)?
            .into())
    }

    /// Contour lines of the surface at each elevation in `levels`.
    ///
    /// Returns
    /// -------
    /// Polylines
    ///     One feature per line, with the feature column ``elevation``.
    fn contours(&self, py: Python, levels: &Bound<PyAny>) -> PyResult<PyPolylines> {
        let l = floats(levels, "levels")?;
        Ok(py
            .detach(|| wireframes_build::contours(&self.inner, &l))
            .map_err(to_py)?
            .into())
    }
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dtm, m)?)?;
    m.add_function(wrap_pyfunction!(triangulate, m)?)?;
    m.add_function(wrap_pyfunction!(loft, m)?)?;
    m.add_function(wrap_pyfunction!(extrude, m)?)?;
    Ok(())
}
