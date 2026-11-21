use crate::Mesh;
use crate::exact::{ImplicitPoint, orient2d_i, orient3d_i};
use crate::intersect::point_in_triangle;
use crate::predicates::{Sign, orient2d, projection_axis, side_of_plane};
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

/// How the oriented faces around one edge bound solid wedges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Radial {
    /// `[i, j]` per wedge of solid: face `i` runs against the edge and `j`,
    /// the next face counter-clockwise about the edge direction, runs along it.
    Paired(Vec<[usize; 2]>),
    /// Two faces lie in the same half-plane.
    Coincident,
    /// Directions do not alternate around the edge.
    Unbalanced,
}

/// Pairs oriented faces around one edge into the wedges of solid they bound.
/// `forward[i]`: face `i` runs the edge along its direction. `turn(a, b)`:
/// sign of the rotation from face `a`'s half-plane to face `b`'s about the
/// edge direction, positive when `b` lies less than a half turn
/// counter-clockwise. `same(a, b)`: whether two half-planes of one plane
/// coincide; asked only when `turn` is zero. Both must be exact.
pub fn radial_pairs(
    forward: &[bool],
    turn: impl Fn(usize, usize) -> Sign,
    same: impl Fn(usize, usize) -> bool,
) -> Radial {
    let n = forward.len();
    if n == 0 || n % 2 == 1 {
        return Radial::Unbalanced;
    }
    let half: Vec<u8> = (0..n)
        .map(|i| match (i, turn(0, i)) {
            (0, _) => 0,
            (_, Sign::Positive) => 1,
            (_, Sign::Negative) => 3,
            (_, Sign::Zero) if same(0, i) => 0,
            _ => 2,
        })
        .collect();
    if half[1..].contains(&0) || half.iter().filter(|&&h| h == 2).count() > 1 {
        return Radial::Coincident;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        half[a].cmp(&half[b]).then_with(|| match turn(a, b) {
            Sign::Positive => std::cmp::Ordering::Less,
            Sign::Negative => std::cmp::Ordering::Greater,
            Sign::Zero => a.cmp(&b),
        })
    });
    if order
        .windows(2)
        .any(|w| half[w[0]] == half[w[1]] && turn(w[0], w[1]) != Sign::Positive)
    {
        return Radial::Coincident;
    }
    if (0..n).any(|k| forward[order[k]] == forward[order[(k + 1) % n]]) {
        return Radial::Unbalanced;
    }
    Radial::Paired(
        (0..n)
            .filter(|&k| !forward[order[k]])
            .map(|k| [order[k], order[(k + 1) % n]])
            .collect(),
    )
}

/// [`radial_pairs`] for faces of `mesh` around the edge whose endpoints have
/// the coordinates `p` and `q`, its direction `q - p`. A face without one
/// corner exactly at `p`, another at `q` and a third off their line bounds no
/// wedge, so the edge is [`Radial::Unbalanced`].
pub fn radial_pairs_at(mesh: &Mesh, faces: &[u32], p: [f64; 3], q: [f64; 3]) -> Radial {
    let key = |x: [f64; 3]| x.map(|c| crate::tolerance::canonical(c).to_bits());
    let (kp, kq) = (key(p), key(q));
    let o = |r: [f64; 3], k: usize| Sign::of(crate::predicates::orient2d_axis(p, q, r, k));
    let mut forward = vec![];
    let mut third = vec![];
    let mut axis = vec![];
    for &f in faces {
        let c = mesh.corners(f);
        let at = |k: [u64; 3]| (0..3).find(|&i| key(c[i]) == k);
        let (Some(i), Some(j)) = (at(kp), at(kq)) else {
            return Radial::Unbalanced;
        };
        if i == j {
            return Radial::Unbalanced;
        }
        let r = c[3 - i - j];
        let Some(k) = (0..3).find(|&k| o(r, k) != Sign::Zero) else {
            return Radial::Unbalanced;
        };
        forward.push((i + 1) % 3 == j);
        third.push(r);
        axis.push(k);
    }
    let turn = |a: usize, b: usize| side_of_plane(third[b], p, q, third[a]);
    let same = |a: usize, b: usize| o(third[a], axis[a]) == o(third[b], axis[a]);
    radial_pairs(&forward, turn, same)
}

/// Shells of a mesh: faces joined across edges with two faces, and across
/// edges with more faces in the pairs [`radial_pairs_at`] finds there.
#[derive(Clone, Debug, PartialEq)]
pub struct Shells {
    /// Shell per face, numbered by first face; `u32::MAX` for faces `edges` skips.
    pub of: Vec<u32>,
    pub count: usize,
    /// Edges with more than two faces that do not pair into wedges of solid.
    pub non_manifold: Vec<u32>,
    /// Edges with more than two faces that pair into wedges of solid: shells
    /// touching each other, or one shell touching itself, along the edge.
    pub contacts: Vec<u32>,
    paired: std::collections::HashMap<u32, Vec<[u32; 2]>>,
}

impl Shells {
    pub fn new(mesh: &Mesh, edges: &EdgeMap) -> Shells {
        let n = mesh.n_faces();
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn root(parent: &mut [u32], mut i: u32) -> u32 {
            while parent[i as usize] != i {
                parent[i as usize] = parent[parent[i as usize] as usize];
                i = parent[i as usize];
            }
            i
        }
        let join = |parent: &mut [u32], a: u32, b: u32| {
            let (a, b) = (root(parent, a), root(parent, b));
            parent[a.max(b) as usize] = a.min(b);
        };
        let mut paired = std::collections::HashMap::new();
        let mut non_manifold = vec![];
        for id in 0..edges.edges.len() as u32 {
            let fs = edges.faces_of_edge(id);
            if fs.len() <= 2 {
                for &g in &fs[1..] {
                    join(&mut parent, fs[0], g);
                }
                continue;
            }
            let [a, b] = edges.edges[id as usize].map(|i| mesh.vertices()[i as usize]);
            match radial_pairs_at(mesh, fs, a, b) {
                Radial::Paired(pairs) => {
                    let pairs: Vec<[u32; 2]> = pairs.iter().map(|p| p.map(|k| fs[k])).collect();
                    for p in &pairs {
                        join(&mut parent, p[0], p[1]);
                    }
                    paired.insert(id, pairs);
                }
                _ => {
                    for &g in &fs[1..] {
                        join(&mut parent, fs[0], g);
                    }
                    non_manifold.push(id);
                }
            }
        }
        let mut contacts: Vec<u32> = paired.keys().copied().collect();
        contacts.sort_unstable();
        non_manifold.sort_unstable();
        let mut of = vec![u32::MAX; n];
        let mut label: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
        for f in 0..n {
            if edges.face_edges[f][0] == u32::MAX {
                continue;
            }
            let r = root(&mut parent, f as u32);
            let next = label.len() as u32;
            of[f] = *label.entry(r).or_insert(next);
        }
        Shells {
            of,
            count: label.len(),
            non_manifold,
            contacts,
            paired,
        }
    }

    /// Face pairs joined across edge `id`: both faces of a two-face edge, the
    /// radial pairs of a contact, and every face with the first otherwise.
    pub fn links(&self, edges: &EdgeMap, id: u32) -> Vec<[u32; 2]> {
        if let Some(p) = self.paired.get(&id) {
            return p.clone();
        }
        let fs = edges.faces_of_edge(id);
        fs[1.min(fs.len())..].iter().map(|&g| [fs[0], g]).collect()
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

/// [`locate_point`] for an exact implicit point, with the same ray and
/// perturbation.
pub fn locate_implicit(vertices: &[[f64; 3]], faces: &[[u32; 3]], p: &ImplicitPoint) -> i8 {
    let x = ImplicitPoint::Explicit;
    let corners = |t: &[u32; 3]| t.map(|i| vertices[i as usize]);
    let side = |c: [[f64; 3]; 3]| orient3d_i(p, c[1], c[2], c[0]);
    let on = |c: [[f64; 3]; 3]| {
        side(c) == Sign::Zero
            && projection_axis(c).is_some_and(|k| {
                let s = orient2d_i(&x(c[0]), &x(c[1]), &x(c[2]), k);
                (0..3).all(|i| orient2d_i(&x(c[i]), &x(c[(i + 1) % 3]), p, k) != s.flip())
            })
    };
    if faces.iter().any(|t| on(corners(t))) {
        return 0;
    }
    let crosses = |c: [[f64; 3]; 3]| {
        let s = orient2d_i(&x(c[0]), &x(c[1]), &x(c[2]), 2);
        let plan = |a: [f64; 3], b: [f64; 3]| match orient2d_i(&x(a), &x(b), p, 2) {
            Sign::Zero => match Sign::of(a[1] - b[1]) {
                Sign::Zero => Sign::of(b[0] - a[0]),
                t => t,
            },
            t => t,
        };
        s != Sign::Zero && (0..3).all(|i| plan(c[i], c[(i + 1) % 3]) == s) && side(c) == s.flip()
    };
    if faces.iter().filter(|t| crosses(corners(t))).count() % 2 == 1 {
        -1
    } else {
        1
    }
}

/// An exact point strictly inside the segment `pq`: where it meets the axis
/// plane through a float strictly between two of its coordinates.
fn inside_segment(p: [f64; 3], q: [f64; 3]) -> Option<ImplicitPoint> {
    (0..3).find_map(|k| {
        let (lo, hi) = (p[k].min(q[k]), p[k].max(q[k]));
        let c = lo * 0.5 + hi * 0.5;
        (lo < c && c < hi).then(|| {
            let mut r = [0.0; 3];
            r[k] = c;
            let (mut s, mut t) = (r, r);
            s[(k + 1) % 3] = 1.0;
            t[(k + 2) % 3] = 1.0;
            ImplicitPoint::Lpi { p, q, r, s, t }
        })
    })
}

/// Location of the closed shell `inner` against the closed shell `outer`: -1
/// inside, 1 outside. Taken from the first vertex of `inner` off `outer`,
/// else from an exact point inside the first edge of `inner` off `outer`
/// (shells may touch at every vertex); `None` when all of those lie on `outer`.
pub fn shell_location(vertices: &[[f64; 3]], inner: &[[u32; 3]], outer: &[[u32; 3]]) -> Option<i8> {
    let v = |i: u32| vertices[i as usize];
    if let Some(l) = inner
        .iter()
        .flatten()
        .map(|&i| locate_point(vertices, outer, v(i)))
        .find(|&l| l != 0)
    {
        return Some(l);
    }
    inner
        .iter()
        .flat_map(|t| (0..3).map(move |e| (t[e], t[(e + 1) % 3])))
        .filter_map(|(a, b)| inside_segment(v(a), v(b)))
        .map(|p| locate_implicit(vertices, outer, &p))
        .find(|&l| l != 0)
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

        #[test]
        fn locate_implicit_matches_locate_point(g in proptest::array::uniform3(-2i32..=14), which in 0usize..3) {
            let m = [fixtures::cube(), fixtures::nested_shells(), fixtures::sphere(1)][which].clone();
            let p = g.map(|x| x as f64 * 0.25 - 0.5);
            let want = locate_point(m.vertices(), m.triangles(), p);
            proptest::prop_assert_eq!(locate_implicit(m.vertices(), m.triangles(), &ImplicitPoint::Explicit(p)), want);
            let q = [p[0] + 0.5, p[1], p[2] + 0.25];
            let x = inside_segment(p, q).expect("coordinates far apart");
            let at = x.to_f64();
            proptest::prop_assert!(at[0] > p[0] && at[0] < q[0]);
        }
    }
}
