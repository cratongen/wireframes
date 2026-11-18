//! Readers and writers. Format modules are pure functions on bytes; this file
//! does all file access.

#![allow(unused_variables)]

use std::path::Path;

use drusa_core::{Mesh, PointSet, Polylines, Result};

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
    todo!()
}

pub fn read(path: impl AsRef<Path>) -> Result<Data> {
    todo!()
}

pub fn read_with(path: impl AsRef<Path>, opts: &ReadOptions) -> Result<Data> {
    todo!()
}

pub fn write(path: impl AsRef<Path>, data: &Data) -> Result<()> {
    todo!()
}

pub fn write_with(path: impl AsRef<Path>, data: &Data, opts: &WriteOptions) -> Result<()> {
    todo!()
}

/// `InvalidInput("file holds no mesh")` when the file has none.
pub fn read_mesh(path: impl AsRef<Path>) -> Result<Mesh> {
    todo!()
}

pub fn write_mesh(path: impl AsRef<Path>, mesh: &Mesh, opts: &WriteOptions) -> Result<()> {
    todo!()
}

pub fn read_polylines(path: impl AsRef<Path>) -> Result<Polylines> {
    todo!()
}

pub fn read_points(path: impl AsRef<Path>) -> Result<PointSet> {
    todo!()
}
