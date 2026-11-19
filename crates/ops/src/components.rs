use wireframes_core::{EdgeMap, Mesh, Result};

/// Edge-connected component label per face, numbered by first face.
pub fn component_labels(mesh: &Mesh) -> Vec<u32> {
    let edges = EdgeMap::new(mesh, None);
    let n = mesh.n_faces();
    let mut label = vec![u32::MAX; n];
    let mut next = 0;
    let mut stack = vec![];
    for f in 0..n {
        if label[f] != u32::MAX {
            continue;
        }
        label[f] = next;
        stack.push(f as u32);
        while let Some(g) = stack.pop() {
            for &e in &edges.face_edges[g as usize] {
                for &h in edges.faces_of_edge(e) {
                    if label[h as usize] == u32::MAX {
                        label[h as usize] = next;
                        stack.push(h);
                    }
                }
            }
        }
        next += 1;
    }
    label
}

/// One mesh per label, faces and vertices in input order.
pub fn components(mesh: &Mesh) -> Result<Vec<Mesh>> {
    let labels = component_labels(mesh);
    let count = labels.iter().max().map_or(0, |&l| l as usize + 1);
    let mut faces: Vec<Vec<u32>> = vec![vec![]; count];
    for (f, &l) in labels.iter().enumerate() {
        faces[l as usize].push(f as u32);
    }
    let mut slot = vec![u32::MAX; mesh.n_vertices()];
    faces
        .iter()
        .map(|fs| {
            let mut verts: Vec<u32> = fs
                .iter()
                .flat_map(|&f| mesh.triangles()[f as usize])
                .collect();
            verts.sort_unstable();
            verts.dedup();
            for (k, &v) in verts.iter().enumerate() {
                slot[v as usize] = k as u32;
            }
            let tris = fs
                .iter()
                .map(|&f| mesh.triangles()[f as usize].map(|i| slot[i as usize]))
                .collect();
            let points = verts.iter().map(|&v| mesh.vertices()[v as usize]).collect();
            Ok(Mesh::new(points, tris)?
                .with_attributes(
                    mesh.vertex_attributes().take(&verts),
                    mesh.face_attributes().take(fs),
                )?
                .with_crs(mesh.crs.clone())
                .with_length_unit(mesh.length_unit.clone()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::measure::volume;
    use wireframes_core::fixtures::{self, broken};

    #[test]
    fn components_split_shells_and_keep_volume() {
        let m = fixtures::nested_shells();
        let labels = component_labels(&m);
        assert_eq!(labels[..12], [0; 12]);
        assert_eq!(labels[12..], [1; 12]);
        let parts = components(&m).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts.iter().map(Mesh::n_faces).sum::<usize>(), m.n_faces());
        let total: f64 = parts.iter().map(|p| volume(p).unwrap()).sum();
        assert_eq!(total, volume(&m).unwrap());
        assert_eq!(parts[1].vertices()[0], [1.0; 3]);
    }

    #[test]
    fn a_shared_vertex_does_not_connect() {
        let parts = components(&broken::touching_tetrahedra()).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[1].n_vertices(), 4);
        assert!(components(&Mesh::empty()).unwrap().is_empty());
    }
}
