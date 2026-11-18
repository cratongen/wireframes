use crate::{Aabb, Column, Error, MeshReport, Result, Table, units};

/// Triangle mesh: an open surface or a closed solid, with per-vertex and
/// per-face attributes.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Mesh {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    vertex_attributes: Table,
    face_attributes: Table,
    pub crs: Option<String>,
    pub length_unit: Option<String>,
}

pub(crate) fn check_finite(what: &str, points: &[[f64; 3]]) -> Result<()> {
    match points.iter().position(|p| p.iter().any(|x| !x.is_finite())) {
        None => Ok(()),
        Some(i) => Err(Error::invalid(format!(
            "{what} must be finite: {what} {i} is {:?}",
            points[i]
        ))),
    }
}

fn check_rows(what: &str, table: &Table, rows: usize) -> Result<()> {
    if table.rows() == rows {
        Ok(())
    } else {
        Err(Error::invalid(format!(
            "{what} table has {} rows, expected {rows}",
            table.rows()
        )))
    }
}

impl Mesh {
    /// Checks that coordinates are finite and triangles index existing vertices.
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> Result<Mesh> {
        check_finite("vertices", &vertices)?;
        let n = vertices.len() as u64;
        if let Some((f, t)) = triangles
            .iter()
            .enumerate()
            .find(|(_, t)| t.iter().any(|&i| i as u64 >= n))
        {
            let bad = t.iter().find(|&&i| i as u64 >= n).expect("found above");
            return Err(Error::invalid(format!(
                "triangles must index existing vertices: face {f} uses {bad} of {n}"
            )));
        }
        let (nv, nf) = (vertices.len(), triangles.len());
        Ok(Mesh {
            vertices,
            triangles,
            vertex_attributes: Table::new(nv),
            face_attributes: Table::new(nf),
            crs: None,
            length_unit: None,
        })
    }

    pub fn empty() -> Mesh {
        Mesh::default()
    }

    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }

    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn n_faces(&self) -> usize {
        self.triangles.len()
    }

    /// No faces.
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    pub fn corners(&self, f: u32) -> [[f64; 3]; 3] {
        self.triangles[f as usize].map(|i| self.vertices[i as usize])
    }

    pub fn vertex_attributes(&self) -> &Table {
        &self.vertex_attributes
    }

    pub fn face_attributes(&self) -> &Table {
        &self.face_attributes
    }

    pub fn with_vertex_column(mut self, name: &str, column: Column) -> Result<Mesh> {
        self.vertex_attributes = self.vertex_attributes.with_column(name, column)?;
        Ok(self)
    }

    pub fn with_face_column(mut self, name: &str, column: Column) -> Result<Mesh> {
        self.face_attributes = self.face_attributes.with_column(name, column)?;
        Ok(self)
    }

    pub fn with_attributes(mut self, vertex: Table, face: Table) -> Result<Mesh> {
        check_rows("vertex attribute", &vertex, self.vertices.len())?;
        check_rows("face attribute", &face, self.triangles.len())?;
        self.vertex_attributes = vertex;
        self.face_attributes = face;
        Ok(self)
    }

    /// The same topology and attributes on new coordinates.
    pub fn with_vertices(mut self, vertices: Vec<[f64; 3]>) -> Result<Mesh> {
        if vertices.len() != self.vertices.len() {
            return Err(Error::invalid(format!(
                "vertices has {} rows, mesh has {}",
                vertices.len(),
                self.vertices.len()
            )));
        }
        check_finite("vertices", &vertices)?;
        self.vertices = vertices;
        Ok(self)
    }

    pub fn with_crs(mut self, crs: Option<String>) -> Mesh {
        self.crs = crs;
        self
    }

    pub fn with_length_unit(mut self, unit: Option<String>) -> Mesh {
        self.length_unit = unit;
        self
    }

    /// Box of every vertex, used or not; `None` without vertices.
    pub fn bounds(&self) -> Option<Aabb> {
        (!self.vertices.is_empty()).then(|| Aabb::from_points(&self.vertices))
    }

    pub fn face_aabbs(&self) -> Vec<Aabb> {
        (0..self.triangles.len() as u32)
            .map(|f| Aabb::from_points(&self.corners(f)))
            .collect()
    }

    /// The faces where `keep` is true with their attributes, and the vertices
    /// they use, in input order.
    pub fn filter_faces(&self, keep: &[bool]) -> Result<Mesh> {
        let face_attributes = self.face_attributes.filter(keep)?;
        let mut new_index = vec![u32::MAX; self.vertices.len()];
        for (t, _) in self.triangles.iter().zip(keep).filter(|(_, k)| **k) {
            for &i in t {
                new_index[i as usize] = 0;
            }
        }
        let mut kept = vec![];
        for (v, slot) in new_index.iter_mut().enumerate() {
            if *slot == 0 {
                *slot = kept.len() as u32;
                kept.push(v as u32);
            }
        }
        let triangles: Vec<[u32; 3]> = self
            .triangles
            .iter()
            .zip(keep)
            .filter(|(_, k)| **k)
            .map(|(t, _)| t.map(|i| new_index[i as usize]))
            .collect();
        Ok(Mesh {
            vertices: kept.iter().map(|&v| self.vertices[v as usize]).collect(),
            triangles,
            vertex_attributes: self.vertex_attributes.take(&kept),
            face_attributes,
            crs: self.crs.clone(),
            length_unit: self.length_unit.clone(),
        })
    }

    /// Vertices and faces of every mesh in order, indices offset, no welding.
    pub fn concat(meshes: &[&Mesh]) -> Result<Mesh> {
        let mut out = Mesh::empty();
        let mut labelled = (None, None);
        for (k, m) in meshes.iter().enumerate() {
            labelled = units::check_same(
                &format!("mesh {k}"),
                (&m.crs, &m.length_unit),
                (&labelled.0, &labelled.1),
            )?;
            let offset = u32::try_from(out.vertices.len())
                .ok()
                .filter(|o| o.checked_add(m.vertices.len() as u32).is_some())
                .ok_or_else(|| {
                    Error::invalid("concatenated mesh has more than u32::MAX vertices")
                })?;
            out.vertex_attributes = Table::concat(&out.vertex_attributes, &m.vertex_attributes)?;
            out.face_attributes = Table::concat(&out.face_attributes, &m.face_attributes)?;
            out.vertices.extend_from_slice(&m.vertices);
            out.triangles
                .extend(m.triangles.iter().map(|t| t.map(|i| i + offset)));
        }
        out.crs = labelled.0;
        out.length_unit = labelled.1;
        Ok(out)
    }

    pub fn validate(&self, tolerance: f64) -> Result<MeshReport> {
        crate::validate::validate(self, tolerance)
    }

    /// At least one non-degenerate face, and every edge of the non-degenerate
    /// faces used by exactly two of them.
    pub fn is_closed(&self) -> bool {
        let skip = crate::validate::degenerate_mask(self);
        if skip.iter().all(|&s| s) {
            return false;
        }
        let edges = crate::EdgeMap::new(self, Some(&skip));
        (0..edges.edges.len() as u32).all(|e| edges.faces_of_edge(e).len() == 2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures;

    #[test]
    fn constructor_rejects_bad_input() {
        let e = Mesh::new(vec![[0.0; 3]; 8], vec![[0, 1, 9]])
            .unwrap_err()
            .to_string();
        assert_eq!(
            e,
            "triangles must index existing vertices: face 0 uses 9 of 8"
        );
        assert!(Mesh::new(vec![[f64::NAN, 0.0, 0.0]], vec![]).is_err());
        assert!(Mesh::empty().is_empty() && Mesh::empty().bounds().is_none());
    }

    #[test]
    fn filter_faces_drops_unused_vertices_and_keeps_attributes() {
        let m = fixtures::cube()
            .with_face_column("id", Column::i64((0..12).collect()))
            .unwrap()
            .with_vertex_column("v", Column::i64((0..8).collect()))
            .unwrap();
        let keep: Vec<bool> = (0..12).map(|f| f < 2).collect();
        let f = m.filter_faces(&keep).unwrap();
        assert_eq!(f.n_faces(), 2);
        assert_eq!(
            f.face_attributes().column("id").unwrap().data,
            crate::ColumnData::Int64(vec![0, 1])
        );
        assert_eq!(f.n_vertices(), f.vertex_attributes().rows());
        for t in f.triangles() {
            for k in 0..3 {
                let crate::ColumnData::Int64(v) = &f.vertex_attributes().column("v").unwrap().data
                else {
                    panic!()
                };
                assert_eq!(
                    f.vertices()[t[k] as usize],
                    m.vertices()[v[t[k] as usize] as usize]
                );
            }
        }
        assert!(m.filter_faces(&[true]).is_err());
    }

    #[test]
    fn concat_offsets_and_checks_labels() {
        let a = fixtures::cube().with_length_unit(Some("m".into()));
        let b = fixtures::tetrahedron();
        let c = Mesh::concat(&[&a, &b]).unwrap();
        assert_eq!((c.n_vertices(), c.n_faces()), (12, 16));
        assert_eq!(c.triangles()[12], b.triangles()[0].map(|i| i + 8));
        assert_eq!(c.length_unit.as_deref(), Some("m"));
        let ft = b.with_length_unit(Some("ft".into()));
        assert!(Mesh::concat(&[&a, &ft]).is_err());
    }

    #[test]
    fn closedness_and_bounds() {
        assert!(fixtures::cube().is_closed());
        assert!(!fixtures::open_cube().is_closed());
        assert!(!Mesh::empty().is_closed());
        assert_eq!(
            fixtures::cube().bounds().unwrap(),
            Aabb {
                min: [0.0; 3],
                max: [1.0; 3]
            }
        );
        assert!(fixtures::cube().with_vertices(vec![[0.0; 3]]).is_err());
        let moved = fixtures::cube().with_vertices(vec![[1.0; 3]; 8]).unwrap();
        assert_eq!(moved.n_faces(), 12);
    }
}
