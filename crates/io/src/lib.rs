//! Readers and writers. Format modules are pure functions on bytes; this file
//! does all file access.

use std::path::{Path, PathBuf};

use wireframes_core::{Error, Mesh, PointSet, Polylines, Result};

pub mod common;
pub mod datamine_csv;
pub mod dxf;
pub mod grid;
pub mod obj;
pub mod off;
pub mod ply;
pub mod stl;
pub mod vtk;
pub mod vtu;

#[doc(hidden)]
pub use common::canonicalize;
pub use grid::UnstructuredGrid;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Data {
    pub mesh: Option<Mesh>,
    pub polylines: Option<Polylines>,
    pub points: Option<PointSet>,
    pub grid: Option<UnstructuredGrid>,
    /// From the file (DXF `$INSUNITS`), never applied as a rescale.
    pub length_unit: Option<String>,
    pub warnings: Vec<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Obj,
    Ply,
    Off,
    Stl,
    Dxf,
    Vtk,
    Vtu,
    Vtp,
    DatamineCsv,
}

/// Labels set on every returned object.
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    pub length_unit: Option<String>,
    pub crs: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum DxfEntity {
    #[default]
    Face3d,
    Polyface,
}

#[derive(Clone, Debug, Default)]
pub struct WriteOptions {
    /// OBJ, OFF, DXF and Datamine are always text.
    pub ascii: bool,
    pub dxf_entity: DxfEntity,
    /// VTU only.
    pub zlib: bool,
    /// Else a column the format cannot store errors.
    pub skip_unsupported_columns: bool,
}

/// By extension, case-insensitive; `*_tr.csv` and `*_pt.csv` are Datamine pairs.
pub fn detect(path: &Path) -> Result<Format> {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let ext = name.rsplit_once('.').map(|(_, e)| e).unwrap_or("");
    Ok(match ext {
        "obj" => Format::Obj,
        "ply" => Format::Ply,
        "off" => Format::Off,
        "stl" => Format::Stl,
        "dxf" => Format::Dxf,
        "vtk" => Format::Vtk,
        "vtu" => Format::Vtu,
        "vtp" => Format::Vtp,
        "csv" if name.ends_with("_tr.csv") || name.ends_with("_pt.csv") => Format::DatamineCsv,
        "csv" => {
            return Err(Error::Unsupported(format!(
                "{}: a csv must be a datamine pair named <name>_tr.csv and <name>_pt.csv",
                path.display()
            )));
        }
        _ => {
            return Err(Error::Unsupported(format!(
                "{}: unknown file extension {ext:?}; wireframes reads and writes .obj, .ply, \
                 .off, .stl, .dxf, .vtk, .vtu, .vtp and datamine pairs named <name>_tr.csv \
                 and <name>_pt.csv",
                path.display()
            )));
        }
    })
}

/// `(points, triangles)` paths of the Datamine pair `path` belongs to.
fn datamine_pair(path: &Path) -> (PathBuf, PathBuf) {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let at = name.len() - "tr.csv".len();
    let kind = &name[at..at + 2];
    let is_pt = kind.eq_ignore_ascii_case("pt");
    let other_kind: String = (if is_pt { "tr" } else { "pt" })
        .chars()
        .zip(kind.chars())
        .map(|(c, o)| {
            if o.is_ascii_uppercase() {
                c.to_ascii_uppercase()
            } else {
                c
            }
        })
        .collect();
    let other = path.with_file_name(format!("{}{other_kind}{}", &name[..at], &name[at + 2..]));
    if is_pt {
        (path.to_path_buf(), other)
    } else {
        (other, path.to_path_buf())
    }
}

fn load(path: &Path) -> Result<Vec<u8>> {
    std::fs::read(path).map_err(|e| {
        Error::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })
}

fn save(path: &Path, bytes: &[u8]) -> Result<()> {
    std::fs::write(path, bytes).map_err(|e| {
        Error::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })
}

pub fn read(path: impl AsRef<Path>) -> Result<Data> {
    read_with(path, &ReadOptions::default())
}

pub fn read_with(path: impl AsRef<Path>, opts: &ReadOptions) -> Result<Data> {
    let path = path.as_ref();
    let format = detect(path)?;
    let parsed = match format {
        Format::DatamineCsv => {
            let (pt, tr) = datamine_pair(path);
            let (pt_bytes, tr_bytes) = (load(&pt)?, load(&tr)?);
            datamine_csv::read_pair(&pt_bytes, &tr_bytes, Some((&pt, &tr)))
        }
        _ => {
            let bytes = load(path)?;
            match format {
                Format::Obj => obj::read(&bytes),
                Format::Ply => ply::read(&bytes),
                Format::Off => off::read(&bytes),
                Format::Stl => stl::read(&bytes),
                Format::Dxf => dxf::read(&bytes),
                Format::Vtk => vtk::read(&bytes),
                Format::Vtu | Format::Vtp => vtu::read(&bytes),
                Format::DatamineCsv => unreachable!("read as a pair above"),
            }
        }
    };
    let mut data = parsed.map_err(|e| e.with_path(path))?;
    if let Some(crs) = &opts.crs {
        data.mesh = data.mesh.map(|m| m.with_crs(Some(crs.clone())));
        data.polylines = data.polylines.map(|p| p.with_crs(Some(crs.clone())));
        data.points = data.points.map(|p| p.with_crs(Some(crs.clone())));
    }
    if let Some(unit) = &opts.length_unit {
        data.mesh = data.mesh.map(|m| m.with_length_unit(Some(unit.clone())));
        data.polylines = data
            .polylines
            .map(|p| p.with_length_unit(Some(unit.clone())));
        data.points = data.points.map(|p| p.with_length_unit(Some(unit.clone())));
    }
    Ok(data)
}

pub fn write(path: impl AsRef<Path>, data: &Data) -> Result<()> {
    write_with(path, data, &WriteOptions::default())
}

pub fn write_with(path: impl AsRef<Path>, data: &Data, opts: &WriteOptions) -> Result<()> {
    let path = path.as_ref();
    let bytes = match detect(path)? {
        Format::Obj => obj::write(data, opts)?,
        Format::Ply => ply::write(data, opts)?,
        Format::Off => off::write(data, opts)?,
        Format::Stl => stl::write(data, opts)?,
        Format::Dxf => dxf::write(data, opts)?,
        Format::Vtk => vtk::write(data, opts)?,
        Format::Vtu => vtu::write(data, opts)?,
        Format::Vtp => vtu::write_vtp(data, opts)?,
        Format::DatamineCsv => {
            let (pt_bytes, tr_bytes) = datamine_csv::write(data, opts)?;
            let (pt, tr) = datamine_pair(path);
            save(&pt, &pt_bytes)?;
            return save(&tr, &tr_bytes);
        }
    };
    save(path, &bytes)
}

/// `InvalidInput("file holds no mesh")` when the file has none.
pub fn read_mesh(path: impl AsRef<Path>) -> Result<Mesh> {
    let path = path.as_ref();
    read(path)?
        .mesh
        .ok_or_else(|| Error::invalid(format!("{}: file holds no mesh", path.display())))
}

pub fn write_mesh(path: impl AsRef<Path>, mesh: &Mesh, opts: &WriteOptions) -> Result<()> {
    let data = Data {
        mesh: Some(mesh.clone()),
        ..Data::default()
    };
    write_with(path, &data, opts)
}

pub fn read_polylines(path: impl AsRef<Path>) -> Result<Polylines> {
    let path = path.as_ref();
    read(path)?
        .polylines
        .ok_or_else(|| Error::invalid(format!("{}: file holds no polylines", path.display())))
}

pub fn read_points(path: impl AsRef<Path>) -> Result<PointSet> {
    let path = path.as_ref();
    read(path)?
        .points
        .ok_or_else(|| Error::invalid(format!("{}: file holds no points", path.display())))
}
