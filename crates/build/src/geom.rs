use std::collections::HashMap;

use spade::handles::FixedVertexHandle;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
use wireframes_core::predicates::{kept_axes, orient2d};
use wireframes_core::{Error, Result};

pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

pub fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(sub(a, b))
}

pub fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n > 0.0 { scale(a, 1.0 / n) } else { [0.0; 3] }
}

pub fn face_normal(t: [[f64; 3]; 3]) -> [f64; 3] {
    cross(sub(t[1], t[0]), sub(t[2], t[0]))
}

/// Newell normal of a closed ring, length twice the vector area, about the
/// ring's first vertex.
pub fn newell(ring: &[[f64; 3]]) -> [f64; 3] {
    let o = ring[0];
    let mut n = [0.0; 3];
    for i in 0..ring.len() {
        let a = sub(ring[i], o);
        let b = sub(ring[(i + 1) % ring.len()], o);
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    n
}

pub fn centroid(points: &[[f64; 3]]) -> [f64; 3] {
    let s = points.iter().fold([0.0; 3], |s, &p| add(s, p));
    scale(s, 1.0 / points.len() as f64)
}

pub fn check_finite(what: &str, points: &[[f64; 3]]) -> Result<()> {
    match points.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
        Some(i) => Err(Error::invalid(format!(
            "{what} must be finite: point {i} is {:?}",
            points[i]
        ))),
        None => Ok(()),
    }
}

pub fn check_distinct_corners(triangles: &[[u32; 3]]) -> Result<()> {
    match triangles
        .iter()
        .position(|t| t[0] == t[1] || t[1] == t[2] || t[2] == t[0])
    {
        Some(f) => Err(Error::invalid(format!(
            "mesh triangles must use three distinct vertices: face {f} is {:?}",
            triangles[f]
        ))),
        None => Ok(()),
    }
}

/// Even-odd inside test of `p` against closed rings; exact through `orient2d`.
pub fn even_odd(rings: &[Vec<[f64; 2]>], p: [f64; 2]) -> bool {
    let mut inside = false;
    for ring in rings {
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            if (a[1] > p[1]) != (b[1] > p[1]) {
                let o = orient2d(a, b, p);
                if (o > 0.0) == (b[1] > a[1]) {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

fn key(p: [f64; 2]) -> (u64, u64) {
    let c = wireframes_core::tolerance::canonical;
    (c(p[0]).to_bits(), c(p[1]).to_bits())
}

/// Planar constrained Delaunay triangulation in a frame `(u, v, w)`: `u, v` are
/// the plan coordinates and `w` is carried along, interpolated on split vertices.
pub struct Cdt {
    tri: ConstrainedDelaunayTriangulation<Point2<f64>>,
    pub verts: Vec<[f64; 3]>,
    handle_of: Vec<FixedVertexHandle>,
    id_of: Vec<u32>,
    by_pos: HashMap<(u64, u64), u32>,
    pinned: Vec<bool>,
}

impl Cdt {
    /// Bulk-loads distinct plan positions; `verts[i]` keeps id `i`.
    pub fn new(verts: Vec<[f64; 3]>) -> Result<Cdt> {
        let mut by_pos = HashMap::with_capacity(verts.len());
        for (i, p) in verts.iter().enumerate() {
            if by_pos.insert(key([p[0], p[1]]), i as u32).is_some() {
                return Err(Error::invalid("cdt vertices must be distinct in plan"));
            }
        }
        let pts: Vec<Point2<f64>> = verts.iter().map(|p| Point2::new(p[0], p[1])).collect();
        let tri = ConstrainedDelaunayTriangulation::<Point2<f64>>::bulk_load(pts)
            .map_err(|e| Error::invalid(format!("coordinates out of range: {e:?}")))?;
        let mut handle_of = vec![FixedVertexHandle::from_index(0); verts.len()];
        let mut id_of = vec![u32::MAX; tri.num_vertices()];
        for v in tri.fixed_vertices() {
            let p = tri.vertex(v).position();
            let id = by_pos[&key([p.x, p.y])];
            handle_of[id as usize] = v;
            id_of[v.index()] = id;
        }
        let pinned = vec![false; verts.len()];
        Ok(Cdt {
            tri,
            verts,
            handle_of,
            id_of,
            by_pos,
            pinned,
        })
    }

    /// Keeps `w` of vertex `id` when a later constraint passes through it.
    pub fn pin(&mut self, id: u32) {
        self.pinned[id as usize] = true;
    }

    fn push(&mut self, p: [f64; 3], h: FixedVertexHandle) -> u32 {
        let id = self.verts.len() as u32;
        self.verts.push(p);
        self.pinned.push(true);
        self.handle_of.push(h);
        if self.id_of.len() <= h.index() {
            self.id_of.resize(h.index() + 1, u32::MAX);
        }
        self.id_of[h.index()] = id;
        self.by_pos.insert(key([p[0], p[1]]), id);
        id
    }

    /// Adds the constraint `a -> b`, splitting crossed constraints; a split
    /// vertex takes `w` interpolated along the constraint it splits (the one
    /// inserted first). An unpinned vertex the constraint passes through takes
    /// `w` interpolated along `a -> b`. Every vertex on the constraint is then
    /// pinned. Returns the number of vertices it created.
    pub fn constrain(&mut self, a: u32, b: u32) -> usize {
        if a == b {
            return 0;
        }
        let (ha, hb) = (self.handle_of[a as usize], self.handle_of[b as usize]);
        let crossed: Vec<[u32; 2]> = self
            .tri
            .get_conflicting_edges_between_vertices(ha, hb)
            .map(|e| [e.from().fix(), e.to().fix()].map(|h| self.id_of[h.index()]))
            .collect();
        let before = self.tri.num_vertices();
        let chain = self.tri.add_constraint_and_split(ha, hb, |p| p);
        let on_line: Vec<FixedVertexHandle> = chain
            .iter()
            .flat_map(|&e| self.tri.directed_edge(e).vertices().map(|v| v.fix()))
            .filter(|h| h.index() < before)
            .collect();
        for h in on_line {
            let id = self.id_of[h.index()];
            if id != a && id != b && !self.pinned[id as usize] {
                let w = self.lerp_w([a, b], self.uv(id));
                self.verts[id as usize][2] = w;
            }
            self.pinned[id as usize] = true;
        }
        let after = self.tri.num_vertices();
        for i in before..after {
            let h = FixedVertexHandle::from_index(i);
            let p = self.tri.vertex(h).position();
            let p = [p.x, p.y];
            let along = crossed
                .iter()
                .map(|&[u, v]| (seg_dist2(self.uv(u), self.uv(v), p), [u, v]))
                .min_by(|x, y| x.0.total_cmp(&y.0))
                .map_or([a, b], |(_, e)| e);
            let w = self.lerp_w(along, p);
            self.push([p[0], p[1], w], h);
        }
        after - before
    }

    fn uv(&self, id: u32) -> [f64; 2] {
        let p = self.verts[id as usize];
        [p[0], p[1]]
    }

    fn lerp_w(&self, [u, v]: [u32; 2], p: [f64; 2]) -> f64 {
        let (a, b) = (self.verts[u as usize], self.verts[v as usize]);
        let d = [b[0] - a[0], b[1] - a[1]];
        let l2 = d[0] * d[0] + d[1] * d[1];
        let t = if l2 > 0.0 {
            (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        a[2] + t * (b[2] - a[2])
    }

    /// Counter-clockwise triangles in `(u, v)`, in spade's face order.
    pub fn triangles(&self) -> Vec<[u32; 3]> {
        self.tri
            .inner_faces()
            .map(|f| f.vertices().map(|v| self.id_of[v.fix().index()]))
            .collect()
    }
}

fn seg_dist2(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let l2 = d[0] * d[0] + d[1] * d[1];
    let t = if l2 > 0.0 {
        (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let q = [a[0] + t * d[0] - p[0], a[1] + t * d[1] - p[1]];
    q[0] * q[0] + q[1] * q[1]
}

/// Triangulates closed rings (even-odd holes) in the projection dropping axis
/// `drop`. Returns the vertices (every ring vertex in order, then vertices
/// created where rings cross) and triangles wound so their normal's `drop`
/// component is positive. Vertices sharing a projected position share an id.
pub fn triangulate_rings(
    rings: &[Vec<[f64; 3]>],
    drop: usize,
) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
    let (i, j) = kept_axes(drop);
    let to_uvw = |p: [f64; 3]| [p[i], p[j], p[drop]];
    let mut first: Vec<[f64; 3]> = vec![];
    let mut seen = HashMap::new();
    let mut ids: Vec<Vec<u32>> = vec![];
    let mut out: Vec<[f64; 3]> = vec![];
    let mut alias: Vec<u32> = vec![];
    for ring in rings {
        let mut r = vec![];
        for &p in ring {
            let q = to_uvw(p);
            let id = *seen.entry(key([q[0], q[1]])).or_insert_with(|| {
                first.push(q);
                (first.len() - 1) as u32
            });
            alias.push(id);
            out.push(p);
            r.push(id);
        }
        ids.push(r);
    }
    let mut cdt = Cdt::new(first)?;
    for r in &ids {
        for k in 0..r.len() {
            cdt.constrain(r[k], r[(k + 1) % r.len()]);
        }
    }
    let rings2: Vec<Vec<[f64; 2]>> = rings
        .iter()
        .map(|r| r.iter().map(|p| [p[i], p[j]]).collect())
        .collect();
    let n_first = seen.len();
    let mut out_id: Vec<u32> = vec![u32::MAX; cdt.verts.len()];
    for (k, &a) in alias.iter().enumerate() {
        if out_id[a as usize] == u32::MAX {
            out_id[a as usize] = k as u32;
        }
    }
    for k in n_first..cdt.verts.len() {
        let q = cdt.verts[k];
        let mut p = [0.0; 3];
        p[i] = q[0];
        p[j] = q[1];
        p[drop] = q[2];
        out_id[k] = out.len() as u32;
        out.push(p);
    }
    let tris = cdt
        .triangles()
        .into_iter()
        .filter(|t| {
            let c = t.map(|v| cdt.uv(v));
            let m = [
                (c[0][0] + c[1][0] + c[2][0]) / 3.0,
                (c[0][1] + c[1][1] + c[2][1]) / 3.0,
            ];
            even_odd(&rings2, m)
        })
        .map(|t| t.map(|v| out_id[v as usize]))
        .collect();
    Ok((out, tris))
}

/// Undirected edge key.
pub fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}
