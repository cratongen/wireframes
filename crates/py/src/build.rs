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
///     Points closer than this in plan to an earlier point are dropped. At 0
///     only exact 3D duplicates are dropped, and two points (spot heights,
///     breakline or boundary vertices) at one plan position with different z
///     raise ``InvalidInput``.
/// max_edge : float, optional
///     Drop boundary triangles with a plan edge longer than this, repeatedly.
///
/// Returns
/// -------
/// Mesh
///     Faces wound counter-clockwise in plan (normals up).
///
/// Examples
/// --------
/// >>> s = wf.dtm([[0, 0, 0], [10, 0, 1], [10, 10, 2], [0, 10, 1], [5, 5, 3]])
/// >>> s.n_faces, s.is_solid
/// (4, False)
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
///
/// Examples
/// --------
/// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
/// >>> wf.triangulate(ring).area
/// 12.0
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
///     Features in order, each with exactly one closed part. The first and
///     last features may instead hold one open part, a line or a single
///     point, where the body pinches out; it gets no cap.
/// cap_ends : bool
///     Close the first and last sections.
///
/// Returns
/// -------
/// Mesh
///     A validated solid with face column ``section`` (band index; caps are
///     -1 and -2). Every section vertex is a vertex of the solid.
///
/// Raises
/// ------
/// InvalidInput
///     When a section holds more than one part (branching is not
///     supported), or a middle section is open.
///
/// Examples
/// --------
/// >>> square = [[0, 0], [4, 0], [4, 4], [0, 4]]
/// >>> sections = wf.Polylines([[[x, y, 0] for x, y in square], [[x, y, 10] for x, y in square]], closed=[True, True])
/// >>> wf.loft(sections).volume
/// 160.0
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
///
/// Examples
/// --------
/// >>> ring = wf.Polylines([[[0, 0, 0], [4, 0, 0], [4, 3, 0], [0, 3, 0]]], closed=[True])
/// >>> wf.extrude(ring, [0, 0, 10]).volume
/// 120.0
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
    /// The surface must have one boundary loop and lie strictly above or
    /// strictly below `z`. Closing a pit design surface above the highest
    /// ground and clipping it by the topography gives the excavation.
    ///
    /// Parameters
    /// ----------
    /// z : float
    ///     Elevation of the flat base.
    ///
    /// Returns
    /// -------
    /// Mesh
    ///     A solid whose volume equals the ``cut_fill`` of the surface
    ///     against a flat plane at `z`; walls and base have null face
    ///     attributes.
    ///
    /// Raises
    /// ------
    /// InvalidInput
    ///     When the surface is not a height field with one boundary loop, or
    ///     a vertex lies at or across `z`.
    ///
    /// Examples
    /// --------
    /// >>> s = wf.dtm([[0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]])
    /// >>> s.close_to_elevation(4.0).volume
    /// 600.0
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
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.orebody()
    /// >>> m.decimate(target_faces=100).n_faces <= 100
    /// True
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
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.topography(n=50)
    /// >>> m.smooth(iterations=5).n_faces == m.n_faces
    /// True
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
    /// Parameters
    /// ----------
    /// distance : float
    ///     Move along each vertex normal; positive grows a solid outward.
    ///
    /// Raises
    /// ------
    /// InvalidMesh
    ///     When the result self-intersects.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [2, 2, 2])
    /// >>> m.offset(0.1).volume > m.volume
    /// True
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
    ///
    /// Parameters
    /// ----------
    /// edge_length : float
    ///     Target edge length, positive.
    /// iterations : int
    ///     Number of split, collapse, flip and relax passes.
    ///
    /// Examples
    /// --------
    /// >>> m = wf.datasets.box([0, 0, 0], [4, 4, 4])
    /// >>> r = m.remesh(1.0)
    /// >>> r.n_faces > m.n_faces, r.is_solid
    /// (True, True)
    #[pyo3(signature = (edge_length, *, iterations=5))]
    fn remesh(&self, py: Python, edge_length: f64, iterations: usize) -> PyResult<PyMesh> {
        Ok(py
            .detach(|| wireframes_build::remesh(&self.inner, edge_length, iterations))
            .map_err(to_py)?
            .into())
    }

    /// Contour lines of the surface at each elevation in `levels`.
    ///
    /// Parameters
    /// ----------
    /// levels : array_like
    ///     Elevations to contour.
    ///
    /// Returns
    /// -------
    /// Polylines
    ///     One feature per line, with the feature column ``elevation``. A
    ///     flat area at a level that the surface does not cross there, such
    ///     as a pit floor, gives its toe; a line that only goes round a flat
    ///     triangle on the other side of a contour is left out.
    ///
    /// Examples
    /// --------
    /// >>> s = wf.dtm([[0, 0, 0], [10, 0, 0], [10, 10, 10], [0, 10, 10]])
    /// >>> lines = s.contours([2.5, 5.0])
    /// >>> lines.attributes["elevation"].tolist(), lines.length().tolist()
    /// ([2.5, 5.0], [10.0, 10.0])
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
