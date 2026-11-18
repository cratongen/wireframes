#![allow(unused_variables)]

use drusa_core::{Mesh, Result, Table};

/// Volume and mixed cells kept as data for round trips; not a modeling type.
#[derive(Clone, Debug, PartialEq)]
pub struct UnstructuredGrid {
    pub points: Vec<[f64; 3]>,
    /// VTK cell type ids, verbatim.
    pub cell_types: Vec<u8>,
    /// `cells + 1` long, starting at 0.
    pub offsets: Vec<i64>,
    pub connectivity: Vec<i64>,
    /// Polyhedron face stream `(offsets, ids)`, verbatim.
    pub faces: Option<(Vec<i64>, Vec<i64>)>,
    pub point_data: Table,
    pub cell_data: Table,
}

impl UnstructuredGrid {
    /// Triangles of the surface cells (VTK types 5 to 9), fanned from their first corner.
    pub fn surface(&self) -> Result<Mesh> {
        todo!()
    }

    pub fn n_cells(&self) -> usize {
        todo!()
    }
}
