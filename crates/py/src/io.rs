use std::path::PathBuf;

use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use wireframes_core::units::length_factor;
use wireframes_io::{Data, DxfEntity, Format, ReadOptions, WriteOptions};

use crate::containers::{PyMesh, PyPointSet, PyPolylines, PyUnstructuredGrid};
use crate::errors::{invalid, to_py};

fn load(
    py: Python,
    path: &PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<Data> {
    let data = py
        .detach(|| wireframes_io::read_with(path, &ReadOptions::default()))
        .map_err(to_py)?;
    let unit = match (length_unit, data.length_unit.clone()) {
        (Some(given), Some(file)) => {
            let (a, b) = (
                length_factor(&given).map_err(to_py)?,
                length_factor(&file).map_err(to_py)?,
            );
            if a != b {
                return Err(invalid(format!(
                    "{}: length_unit {given:?} disagrees with the file unit {file:?}",
                    path.display()
                )));
            }
            Some(given)
        }
        (given, file) => given.or(file),
    };
    Ok(Data {
        mesh: data
            .mesh
            .map(|m| m.with_crs(crs.clone()).with_length_unit(unit.clone())),
        polylines: data
            .polylines
            .map(|p| p.with_crs(crs.clone()).with_length_unit(unit.clone())),
        points: data
            .points
            .map(|p| p.with_crs(crs.clone()).with_length_unit(unit.clone())),
        ..data
    })
}

/// ``(mesh, polylines, points, grid, length_unit, warnings)``; ``wireframes.read``
/// wraps it in a ``Dataset``.
#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn _read(
    py: Python,
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<Py<PyAny>> {
    let d = load(py, &path, crs, length_unit)?;
    (
        d.mesh.map(PyMesh::from),
        d.polylines.map(PyPolylines::from),
        d.points.map(PyPointSet::from),
        d.grid.map(|g| PyUnstructuredGrid { inner: g }),
        d.length_unit,
        d.warnings,
    )
        .into_py_any(py)
}

fn none_in(path: &std::path::Path, what: &str) -> PyErr {
    invalid(format!("{}: file holds no {what}", path.display()))
}

/// The mesh in the file at `path`.
///
/// Parameters
/// ----------
/// path : str or PathLike
///     File to read; see ``wireframes.read`` for the formats.
/// crs, length_unit : str, optional
///     Labels of the result; None takes the unit the file declares.
///
/// Raises
/// ------
/// InvalidInput
///     When the file holds no mesh.
#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn read_mesh(
    py: Python,
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<PyMesh> {
    let d = load(py, &path, crs, length_unit)?;
    Ok(d.mesh.ok_or_else(|| none_in(&path, "mesh"))?.into())
}

/// The polylines in the file at `path`; see ``read_mesh``.
#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn read_polylines(
    py: Python,
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<PyPolylines> {
    let d = load(py, &path, crs, length_unit)?;
    Ok(d.polylines
        .ok_or_else(|| none_in(&path, "polylines"))?
        .into())
}

/// The points in the file at `path`; see ``read_mesh``.
#[pyfunction]
#[pyo3(signature = (path, *, crs=None, length_unit=None))]
fn read_points(
    py: Python,
    path: PathBuf,
    crs: Option<String>,
    length_unit: Option<String>,
) -> PyResult<PyPointSet> {
    let d = load(py, &path, crs, length_unit)?;
    Ok(d.points.ok_or_else(|| none_in(&path, "points"))?.into())
}

fn data_arg(data: &Bound<PyAny>) -> PyResult<Data> {
    if let Ok(m) = data.cast::<PyMesh>() {
        return Ok(Data {
            mesh: Some(m.get().inner.clone()),
            ..Data::default()
        });
    }
    if let Ok(p) = data.cast::<PyPolylines>() {
        return Ok(Data {
            polylines: Some(p.get().inner.clone()),
            ..Data::default()
        });
    }
    if let Ok(p) = data.cast::<PyPointSet>() {
        return Ok(Data {
            points: Some(p.get().inner.clone()),
            ..Data::default()
        });
    }
    if let Ok(g) = data.cast::<PyUnstructuredGrid>() {
        return Ok(Data {
            grid: Some(g.get().inner.clone()),
            ..Data::default()
        });
    }
    let bad = || invalid("data must be a Mesh, Polylines, PointSet, UnstructuredGrid or Dataset");
    if !data.hasattr("warnings")? {
        return Err(bad());
    }
    let part = |name: &str| -> PyResult<Option<Bound<PyAny>>> {
        let v = data.getattr(name)?;
        Ok((!v.is_none()).then_some(v))
    };
    Ok(Data {
        mesh: part("mesh")?
            .map(|m| Ok::<_, PyErr>(m.cast::<PyMesh>().map_err(|_| bad())?.get().inner.clone()))
            .transpose()?,
        polylines: part("polylines")?
            .map(|m| {
                Ok::<_, PyErr>(
                    m.cast::<PyPolylines>()
                        .map_err(|_| bad())?
                        .get()
                        .inner
                        .clone(),
                )
            })
            .transpose()?,
        points: part("points")?
            .map(|m| {
                Ok::<_, PyErr>(
                    m.cast::<PyPointSet>()
                        .map_err(|_| bad())?
                        .get()
                        .inner
                        .clone(),
                )
            })
            .transpose()?,
        grid: part("grid")?
            .map(|m| {
                Ok::<_, PyErr>(
                    m.cast::<PyUnstructuredGrid>()
                        .map_err(|_| bad())?
                        .get()
                        .inner
                        .clone(),
                )
            })
            .transpose()?,
        length_unit: part("length_unit")?.map(|u| u.extract()).transpose()?,
        warnings: vec![],
    })
}

/// Writes `data` to `path`; the extension picks the format.
///
/// Parameters
/// ----------
/// path : str or PathLike
///     Target file; a Datamine pair writes both ``_pt.csv`` and ``_tr.csv``.
/// data : Mesh, Polylines, PointSet, UnstructuredGrid or Dataset
///     What to write.
/// ascii : bool
///     Text instead of binary for PLY, STL and VTK.
/// dxf_entity : {"3dface", "polyface"}
///     How DXF stores triangles.
/// zlib : bool
///     Compress VTU arrays.
/// skip_unsupported_columns : bool
///     Leave out columns the format cannot store instead of raising.
///
/// Raises
/// ------
/// InvalidInput
///     When the format cannot hold part of `data`.
/// FileError
///     When the file cannot be written.
///
/// Warns
/// -----
/// UserWarning
///     When binary STL rounds coordinates to float32 (UTM coordinates do), naming
///     how many vertices moved and by how much.
#[pyfunction]
#[pyo3(signature = (path, data, *, ascii=false, dxf_entity="3dface", zlib=false, skip_unsupported_columns=false))]
fn write(
    py: Python,
    path: PathBuf,
    data: &Bound<PyAny>,
    ascii: bool,
    dxf_entity: &str,
    zlib: bool,
    skip_unsupported_columns: bool,
) -> PyResult<()> {
    let dxf_entity = match dxf_entity {
        "3dface" => DxfEntity::Face3d,
        "polyface" => DxfEntity::Polyface,
        _ => {
            return Err(invalid(format!(
                "dxf_entity must be \"3dface\" or \"polyface\", not {dxf_entity:?}"
            )));
        }
    };
    let data = data_arg(data)?;
    let opts = WriteOptions {
        ascii,
        dxf_entity,
        zlib,
        skip_unsupported_columns,
    };
    py.detach(|| wireframes_io::write_with(&path, &data, &opts))
        .map_err(to_py)?;
    let binary_stl = !ascii && wireframes_io::detect(&path).ok() == Some(Format::Stl);
    if let Some(m) = data.mesh.as_ref().filter(|_| binary_stl) {
        let (moved, worst) = wireframes_io::stl::f32_rounding(m);
        if moved > 0 {
            let msg = format!(
                "{}: binary stl stores coordinates as float32, so {moved} of {} vertices moved \
                 by up to {worst:.3e} and the mesh may no longer be valid; write ascii=True, \
                 or translate to a local origin first",
                path.display(),
                m.n_vertices()
            );
            let category = py.get_type::<pyo3::exceptions::PyUserWarning>();
            PyErr::warn(py, &category, &std::ffi::CString::new(msg)?, 1)?;
        }
    }
    Ok(())
}

pub(crate) fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_read, m)?)?;
    m.add_function(wrap_pyfunction!(read_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(read_polylines, m)?)?;
    m.add_function(wrap_pyfunction!(read_points, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    Ok(())
}
