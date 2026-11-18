use crate::Mesh;
use crate::intersect::point_in_triangle;
use crate::predicates::{Sign, orient2d, side_of_plane};
use crate::{Error, Result};

/// Undirected edges of a mesh with face and vertex incidence in CSR form.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeMap {
    /// Sorted endpoint pairs, ordered by `(min, max)`.
    pub edges: Vec<[u32; 2]>,
    /// Edge id of local edge `e` (corner `e` to `(e + 1) % 3`); `u32::MAX` for skipped faces.
    pub face_edges: Vec<[u32; 3]>,
    pub edge_face_offsets: Vec<u32>,
    /// Faces of each edge, ascending.
    pub edge_faces: Vec<u32>,
    pub vertex_face_offsets: Vec<u32>,
    pub vertex_faces: Vec<u32>,
}

fn csr(n: usize, pairs: &mut [(u32, u32)]) -> (Vec<u32>, Vec<u32>) {
    pairs.sort_unstable();
    let mut offsets = vec![0u32; n + 1];
    for &(k, _) in pairs.iter() {
        offsets[k as usize + 1] += 1;
    }
    for i in 0..n {
        offsets[i + 1] += offsets[i];
    }
    (offsets, pairs.iter().map(|&(_, v)| v).collect())
}

impl EdgeMap {
    /// `skip` marks faces to leave out, such as degenerate ones.
    pub fn new(mesh: &Mesh, skip: Option<&[bool]>) -> EdgeMap {
        let tris = mesh.triangles();
        let used = |f: usize| skip.is_none_or(|s| !s[f]);
        let mut keyed: Vec<([u32; 2], u32, u8)> = vec![];
        for (f, t) in tris.iter().enumerate().filter(|(f, _)| used(*f)) {
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                keyed.push(([a.min(b), a.max(b)], f as u32, e as u8));
            }
        }
        keyed.sort_unstable();
        let mut edges: Vec<[u32; 2]> = vec![];
        let mut face_edges = vec![[u32::MAX; 3]; tris.len()];
        let mut pairs = vec![];
        for (key, f, e) in keyed {
            if edges.last() != Some(&key) {
                edges.push(key);
            }
            let id = edges.len() as u32 - 1;
            face_edges[f as usize][e as usize] = id;
            if pairs.last() != Some(&(id, f)) {
                pairs.push((id, f));
            }
        }
        let (edge_face_offsets, edge_faces) = csr(edges.len(), &mut pairs);
        let mut vpairs: Vec<(u32, u32)> = vec![];
        for (f, t) in tris.iter().enumerate().filter(|(f, _)| used(*f)) {
            let mut seen: Vec<u32> = vec![];
            for &v in t {
                if !seen.contains(&v) {
                    seen.push(v);
                    vpairs.push((v, f as u32));
                }
            }
        }
        let (vertex_face_offsets, vertex_faces) = csr(mesh.n_vertices(), &mut vpairs);
        EdgeMap {
            edges,
            face_edges,
            edge_face_offsets,
            edge_faces,
            vertex_face_offsets,
            vertex_faces,
        }
    }

    pub fn faces_of_edge(&self, e: u32) -> &[u32] {
        let e = e as usize;
        &self.edge_faces[self.edge_face_offsets[e] as usize..self.edge_face_offsets[e + 1] as usize]
    }

    pub fn faces_of_vertex(&self, v: u32) -> &[u32] {
        let v = v as usize;
        &self.vertex_faces
            [self.vertex_face_offsets[v] as usize..self.vertex_face_offsets[v + 1] as usize]
    }

    pub fn boundary_edges(&self) -> Vec<u32> {
        (0..self.edges.len() as u32)
            .filter(|&e| self.faces_of_edge(e).len() == 1)
            .collect()
    }

    pub fn non_manifold_edges(&self) -> Vec<u32> {
        (0..self.edges.len() as u32)
            .filter(|&e| self.faces_of_edge(e).len() > 2)
            .collect()
    }
}

/// Twin half-edges: half-edge `h = 3 * f + e` runs from corner `e` to `(e + 1) % 3`
/// of face `f`. Twins pair the two half-edges of an edge whatever their
/// directions; `u32::MAX` on the boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct HalfEdges {
    pub twin: Vec<u32>,
}

impl HalfEdges {
    /// `InvalidInput` when an edge has more than two faces.
    pub fn new(mesh: &Mesh) -> Result<HalfEdges> {
        let map = EdgeMap::new(mesh, None);
        let mut twin = vec![u32::MAX; 3 * mesh.n_faces()];
        let mut halves: Vec<Vec<u32>> = vec![vec![]; map.edges.len()];
        for (f, ids) in map.face_edges.iter().enumerate() {
            for e in 0..3 {
                halves[ids[e] as usize].push(3 * f as u32 + e as u32);
            }
        }
        for (id, h) in halves.iter().enumerate() {
            match h[..] {
                [_] => {}
                [a, b] => {
                    twin[a as usize] = b;
                    twin[b as usize] = a;
                }
                _ => {
                    let [u, v] = map.edges[id];
                    return Err(Error::invalid(format!(
                        "edge {u}-{v} has {} half-edges; half-edges need at most 2",
                        h.len()
                    )));
                }
            }
        }
        Ok(HalfEdges { twin })
    }

    pub fn next(h: u32) -> u32 {
        h - h % 3 + (h + 1) % 3
    }

    pub fn prev(h: u32) -> u32 {
        h - h % 3 + (h + 2) % 3
    }

    pub fn face(h: u32) -> u32 {
        h / 3
    }

    pub fn origin(&self, mesh: &Mesh, h: u32) -> u32 {
        mesh.triangles()[(h / 3) as usize][(h % 3) as usize]
    }

    fn target(&self, mesh: &Mesh, h: u32) -> u32 {
        self.origin(mesh, Self::next(h))
    }

    /// Vertex loops of boundary half-edges. Each loop starts at the smallest
    /// unused boundary half-edge and continues with the smallest unused one
    /// leaving the current vertex; a chain with no way on is returned open.
    pub fn boundary_loops(&self, mesh: &Mesh) -> Vec<Vec<u32>> {
        let boundary: Vec<u32> = (0..self.twin.len() as u32)
            .filter(|&h| self.twin[h as usize] == u32::MAX)
            .collect();
        let mut leaving: Vec<(u32, u32)> = boundary
            .iter()
            .map(|&h| (self.origin(mesh, h), h))
            .collect();
        leaving.sort_unstable();
        let mut used = vec![false; self.twin.len()];
        let mut loops = vec![];
        for &start in &boundary {
            if used[start as usize] {
                continue;
            }
            let mut ring = vec![];
            let mut h = start;
            loop {
                used[h as usize] = true;
                ring.push(self.origin(mesh, h));
                let v = self.target(mesh, h);
                let from = leaving.partition_point(|&(o, _)| o < v);
                match leaving[from..]
                    .iter()
                    .take_while(|&&(o, _)| o == v)
                    .find(|&&(_, g)| !used[g as usize])
                {
                    Some(&(_, g)) => h = g,
                    None => break,
                }
            }
            loops.push(ring);
        }
        loops
    }

    /// Edge-connected component label per face, numbered by first face.
    pub fn components(&self, mesh: &Mesh) -> Vec<u32> {
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
                for e in 0..3 {
                    let t = self.twin[3 * g as usize + e];
                    if t != u32::MAX && label[(t / 3) as usize] == u32::MAX {
                        label[(t / 3) as usize] = next;
                        stack.push(t / 3);
                    }
                }
            }
            next += 1;
        }
        label
    }
}

/// Sign of `orient2d(a, b, p)` with `p` moved to `p + (e, e^2)` for an
/// infinitesimal `e > 0`, so it is never zero when `a != b`.
fn perturbed_orient(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> Sign {
    match Sign::of(orient2d(a, b, p)) {
        Sign::Zero => match Sign::of(a[1] - b[1]) {
            Sign::Zero => Sign::of(b[0] - a[0]),
            s => s,
        },
        s => s,
    }
}

/// Exact location of `p` against the closed triangle set `faces`: -1 inside,
/// 0 on a face, 1 outside. Counts crossings of the ray `+z` from `p` with
/// symbolic perturbation of `p` in plan; meaningful when `faces` bound a closed
/// set. Brute force over every face.
pub fn locate_point(vertices: &[[f64; 3]], faces: &[[u32; 3]], p: [f64; 3]) -> i8 {
    let corners = |t: &[u32; 3]| t.map(|i| vertices[i as usize]);
    if faces.iter().any(|t| {
        let c = corners(t);
        !crate::predicates::collinear3(c[0], c[1], c[2]) && point_in_triangle(p, c)
    }) {
        return 0;
    }
    let crossings = faces.iter().filter(|t| ray_crosses(corners(t), p)).count();
    if crossings % 2 == 1 { -1 } else { 1 }
}

/// Whether the ray `+z` from `p` (perturbed in plan) crosses face `c`. `p` must
/// not lie on the closed face.
pub fn ray_crosses(c: [[f64; 3]; 3], p: [f64; 3]) -> bool {
    let plan = c.map(|v| [v[0], v[1]]);
    let s = Sign::of(orient2d(plan[0], plan[1], plan[2]));
    if s == Sign::Zero {
        return false;
    }
    let q = [p[0], p[1]];
    if (0..3).any(|i| perturbed_orient(plan[i], plan[(i + 1) % 3], q) != s) {
        return false;
    }
    side_of_plane(p, c[0], c[1], c[2]) == s.flip()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures;

    #[test]
    fn edge_map_of_a_cube() {
        let m = fixtures::cube();
        let e = EdgeMap::new(&m, None);
        assert_eq!(e.edges.len(), 18);
        assert!(e.boundary_edges().is_empty() && e.non_manifold_edges().is_empty());
        assert!(e.edges.windows(2).all(|w| w[0] < w[1]));
        for (f, ids) in e.face_edges.iter().enumerate() {
            for &id in ids {
                assert!(e.faces_of_edge(id).contains(&(f as u32)));
            }
        }
        let fan: usize = (0..8).map(|v| e.faces_of_vertex(v).len()).sum();
        assert_eq!(fan, 36);
        let skip: Vec<bool> = (0..12).map(|f| f == 0).collect();
        let e = EdgeMap::new(&m, Some(&skip));
        assert_eq!(e.boundary_edges().len(), 3);
        assert_eq!(e.face_edges[0], [u32::MAX; 3]);
    }

    #[test]
    fn half_edges_twins_loops_and_components() {
        let m = fixtures::open_cube();
        let h = HalfEdges::new(&m).unwrap();
        for (a, &b) in h.twin.iter().enumerate() {
            if b != u32::MAX {
                assert_eq!(h.twin[b as usize], a as u32);
            }
        }
        let loops = h.boundary_loops(&m);
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].len(), 4);
        let two = Mesh::concat(&[&fixtures::cube(), &fixtures::tetrahedron()]).unwrap();
        let labels = HalfEdges::new(&two).unwrap().components(&two);
        assert_eq!(labels.iter().filter(|&&l| l == 0).count(), 12);
        assert_eq!(labels.iter().filter(|&&l| l == 1).count(), 4);
        assert!(HalfEdges::new(&fixtures::broken::non_manifold_fin()).is_err());
        assert_eq!(
            (HalfEdges::next(5), HalfEdges::prev(3), HalfEdges::face(7)),
            (3, 5, 2)
        );
    }

    #[test]
    fn locate_point_in_a_cube() {
        let m = fixtures::cube();
        let loc = |p| locate_point(m.vertices(), m.triangles(), p);
        assert_eq!(loc([0.5, 0.5, 0.5]), -1);
        assert_eq!(loc([0.5, 0.5, 1.5]), 1);
        assert_eq!(loc([1.0, 1.0, 1.0]), 0);
        assert_eq!(loc([0.5, 0.5, 1.0]), 0);
        assert_eq!(loc([0.0, 0.0, 0.5]), 0);
        assert_eq!(loc([0.5, 0.5, -0.5]), 1);
        assert_eq!(loc([1.0, 0.5, -0.5]), 1);
        assert_eq!(loc([0.5, 0.5, 0.5f64.next_up()]), -1);
    }

    proptest::proptest! {
        #[test]
        fn locate_matches_box_test_on_grid_points(g in proptest::array::uniform3(-2i32..=6)) {
            let m = fixtures::cube();
            let p = g.map(|x| x as f64 * 0.25);
            let inside = p.iter().all(|&x| 0.0 < x && x < 1.0);
            let on = p.iter().all(|&x| (0.0..=1.0).contains(&x)) && !inside;
            let want = if inside { -1 } else if on { 0 } else { 1 };
            proptest::prop_assert_eq!(locate_point(m.vertices(), m.triangles(), p), want);
        }
    }
}
