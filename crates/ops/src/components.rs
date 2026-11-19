use wireframes_core::topology::Shells;
use wireframes_core::validate::degenerate_mask;
use wireframes_core::{EdgeMap, Mesh, Result};

/// Edge-connected component label per face, numbered by first face. Shells
/// that touch along an edge are separate components ([`Shells`]); a degenerate
/// face joins every face it shares an edge with.
pub fn component_labels(mesh: &Mesh) -> Vec<u32> {
    let n = mesh.n_faces();
    let skip = degenerate_mask(mesh);
    let shells = Shells::new(mesh, &EdgeMap::new(mesh, Some(&skip)));
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn root(parent: &mut [u32], mut i: u32) -> u32 {
        while parent[i as usize] != i {
            parent[i as usize] = parent[parent[i as usize] as usize];
            i = parent[i as usize];
        }
        i
    }
    let mut join = |a: u32, b: u32| {
        let (a, b) = (root(&mut parent, a), root(&mut parent, b));
        parent[a.max(b) as usize] = a.min(b);
    };
    let mut first: Vec<u32> = vec![u32::MAX; shells.count];
    for f in 0..n as u32 {
        let s = shells.of[f as usize];
        if s != u32::MAX {
            if first[s as usize] == u32::MAX {
                first[s as usize] = f;
            }
            join(first[s as usize], f);
        }
    }
    if skip.contains(&true) {
        let all = EdgeMap::new(mesh, None);
        for e in 0..all.edges.len() as u32 {
            let fs = all.faces_of_edge(e);
            if fs.iter().any(|&f| skip[f as usize]) {
                for &g in &fs[1..] {
                    join(fs[0], g);
                }
            }
        }
    }
    let mut label = vec![u32::MAX; n];
    let mut next = 0;
    let mut of_root: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
    for f in 0..n as u32 {
        let r = root(&mut parent, f);
        label[f as usize] = *of_root.entry(r).or_insert_with(|| {
            next += 1;
            next - 1
        });
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
    use wireframes_core::fixtures;

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
        let parts = components(&fixtures::touching_tetrahedra()).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[1].n_vertices(), 4);
        assert!(components(&Mesh::empty()).unwrap().is_empty());
    }
}
