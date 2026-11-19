use wireframes_core::{Error, Mesh, Result, Table};

use crate::vtk::Fan;

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
    /// Triangles of the surface cells (VTK types 5 to 9), fanned from their
    /// first corner; strips alternate winding and pixels run 0, 1, 3, 2. The
    /// mesh keeps every grid point and its point data; face attributes are
    /// the cell data of the source cell.
    pub fn surface(&self) -> Result<Mesh> {
        let n = self.n_cells();
        if self.offsets.len() != n + 1 {
            return Err(Error::invalid(format!(
                "offsets has {} values, there are {n} cells",
                self.offsets.len()
            )));
        }
        let mut fan = Fan::default();
        for c in 0..n {
            let (a, b) = (self.offsets[c], self.offsets[c + 1]);
            if a < 0 || b < a || b as usize > self.connectivity.len() {
                return Err(Error::invalid(format!("cell {c} has offsets {a}..{b}")));
            }
            let ids: Vec<u32> = self.connectivity[a as usize..b as usize]
                .iter()
                .map(|&i| {
                    u32::try_from(i).map_err(|_| Error::invalid(format!("cell {c} uses point {i}")))
                })
                .collect::<Result<_>>()?;
            match self.cell_types[c] {
                5 | 7 | 9 if ids.len() >= 3 => fan.polygon(&ids, c as u32),
                6 if ids.len() >= 3 => fan.strip(&ids, c as u32),
                8 if ids.len() == 4 => fan.polygon(&[ids[0], ids[1], ids[3], ids[2]], c as u32),
                5..=9 => {
                    return Err(Error::invalid(format!(
                        "cell {c} of type {} has {} points",
                        self.cell_types[c],
                        ids.len()
                    )));
                }
                _ => {}
            }
        }
        let faces = self.cell_data.take(&fan.source);
        fan.mesh(self.points.clone(), self.point_data.clone(), faces)
    }

    pub fn n_cells(&self) -> usize {
        self.cell_types.len()
    }
}

#[cfg(test)]
mod tests {
    use wireframes_core::Column;

    use super::*;

    fn grid(cell_types: Vec<u8>, cells: &[&[i64]]) -> UnstructuredGrid {
        let mut offsets = vec![0];
        let mut connectivity = vec![];
        for c in cells {
            connectivity.extend_from_slice(c);
            offsets.push(connectivity.len() as i64);
        }
        let n = cell_types.len() as i64;
        UnstructuredGrid {
            points: (0..6).map(|i| [i as f64, (i * i) as f64, 0.0]).collect(),
            cell_types,
            offsets,
            connectivity,
            faces: None,
            point_data: Table::new(6),
            cell_data: Table::new(n as usize)
                .with_column("id", Column::i64((0..n).collect()))
                .unwrap(),
        }
    }

    #[test]
    fn surface_keeps_surface_cells_and_skips_volumes() {
        let g = grid(
            vec![10, 5, 9, 8],
            &[&[0, 1, 2, 3], &[0, 1, 2], &[0, 1, 2, 3], &[0, 1, 2, 3]],
        );
        assert_eq!(g.n_cells(), 4);
        let m = g.surface().unwrap();
        assert_eq!(m.n_vertices(), 6);
        assert_eq!(
            m.triangles(),
            &[[0, 1, 2], [0, 1, 2], [0, 2, 3], [0, 1, 3], [0, 3, 2]]
        );
        let ids = &m.face_attributes().column("id").unwrap().data;
        assert_eq!(
            *ids,
            wireframes_core::ColumnData::Int64(vec![1, 2, 2, 3, 3])
        );
    }

    #[test]
    fn strips_alternate_winding() {
        let m = grid(vec![6], &[&[0, 1, 2, 3, 4]]).surface().unwrap();
        assert_eq!(m.triangles(), &[[0, 1, 2], [2, 1, 3], [2, 3, 4]]);
    }

    #[test]
    fn short_cells_and_bad_ids_error() {
        assert!(grid(vec![5], &[&[0, 1]]).surface().is_err());
        assert!(grid(vec![5], &[&[0, 1, -1]]).surface().is_err());
        assert!(grid(vec![5], &[&[0, 1, 9]]).surface().is_err());
    }
}
