use std::cmp::Reverse;
use std::collections::BinaryHeap;

use drusa_core::predicates::{Sign, collinear3, side_of_plane};
use drusa_core::{Error, HalfEdges, Mesh, Result, tolerance};

use crate::geom::{
    add, check_distinct_corners, check_finite, cross, dot, face_normal, normalize, scale, sub,
};

#[derive(Clone, Debug, Default)]
pub struct DecimateOptions {
    pub target_faces: Option<usize>,
    pub max_error: Option<f64>,
    pub preserve_boundary: bool,
}

type Quadric = [f64; 10];

fn plane_quadric(n: [f64; 3], p: [f64; 3]) -> Quadric {
    let d = -dot(n, p);
    let [a, b, c] = n;
    [
        a * a,
        a * b,
        a * c,
        a * d,
        b * b,
        b * c,
        b * d,
        c * c,
        c * d,
        d * d,
    ]
}

fn qadd(a: &mut Quadric, b: &Quadric) {
    for k in 0..10 {
        a[k] += b[k];
    }
}

fn qeval(q: &Quadric, p: [f64; 3]) -> f64 {
    let [x, y, z] = p;
    let v = q[0] * x * x
        + 2.0 * q[1] * x * y
        + 2.0 * q[2] * x * z
        + 2.0 * q[3] * x
        + q[4] * y * y
        + 2.0 * q[5] * y * z
        + 2.0 * q[6] * y
        + q[7] * z * z
        + 2.0 * q[8] * z
        + q[9];
    v.max(0.0)
}

fn qsolve(q: &Quadric) -> Option<[f64; 3]> {
    let m = [[q[0], q[1], q[2]], [q[1], q[4], q[5]], [q[2], q[5], q[7]]];
    let r = [-q[3], -q[6], -q[8]];
    let det3 = |m: [[f64; 3]; 3]| dot(m[0], cross(m[1], m[2]));
    let det = det3(m);
    let s = m.iter().flatten().fold(0.0f64, |s, x| s.max(x.abs()));
    if det.is_nan() || det.abs() <= 1e-10 * s * s * s {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut mk = m;
        for row in 0..3 {
            mk[row][k] = r[row];
        }
        x[k] = det3(mk) / det;
    }
    x.iter().all(|c| c.is_finite()).then_some(x)
}

struct State {
    pos: Vec<[f64; 3]>,
    origin: [f64; 3],
    q: Vec<Quadric>,
    faces: Vec<[u32; 3]>,
    face_alive: Vec<bool>,
    vf: Vec<Vec<u32>>,
    alive: Vec<bool>,
    stamp: Vec<u32>,
}

struct Collapse {
    keep: u32,
    gone: u32,
    target: [f64; 3],
    cost: f64,
}

impl State {
    fn faces_of(&self, v: u32) -> impl Iterator<Item = u32> + '_ {
        self.vf[v as usize]
            .iter()
            .copied()
            .filter(|&f| self.face_alive[f as usize])
    }

    fn neighbours(&self, v: u32) -> Vec<u32> {
        let mut n: Vec<u32> = self
            .faces_of(v)
            .flat_map(|f| self.faces[f as usize])
            .filter(|&w| w != v)
            .collect();
        n.sort_unstable();
        n.dedup();
        n
    }

    fn edge_faces(&self, a: u32, b: u32) -> Vec<u32> {
        self.faces_of(a)
            .filter(|&f| self.faces[f as usize].contains(&b))
            .collect()
    }

    fn boundary_neighbours(&self, v: u32) -> Vec<u32> {
        self.neighbours(v)
            .into_iter()
            .filter(|&w| self.edge_faces(v, w).len() == 1)
            .collect()
    }

    fn cost(&self, a: u32, b: u32, p: [f64; 3]) -> f64 {
        let mut q = self.q[a as usize];
        qadd(&mut q, &self.q[b as usize]);
        qeval(&q, sub(p, self.origin))
    }

    fn plan(&self, a: u32, b: u32, preserve: bool) -> Option<Collapse> {
        let (pa, pb) = (self.pos[a as usize], self.pos[b as usize]);
        let at = |keep: u32, gone: u32, target: [f64; 3]| Collapse {
            keep,
            gone,
            target,
            cost: self.cost(a, b, target),
        };
        if preserve {
            let (ba, bb) = (self.boundary_neighbours(a), self.boundary_neighbours(b));
            match (ba.is_empty(), bb.is_empty()) {
                (false, true) => return Some(at(a, b, pa)),
                (true, false) => return Some(at(b, a, pb)),
                (false, false) => {
                    let on_line = |v: u32, nb: &[u32], other: u32| {
                        nb.len() == 2
                            && nb.contains(&other)
                            && nb.iter().all(|&w| {
                                w == other
                                    || collinear3(
                                        self.pos[w as usize],
                                        self.pos[v as usize],
                                        self.pos[other as usize],
                                    )
                            })
                    };
                    let mut best: Option<Collapse> = None;
                    for (keep, gone, nb, p) in [(a, b, &bb, pa), (b, a, &ba, pb)] {
                        if on_line(gone, nb, keep) {
                            let c = at(keep, gone, p);
                            if best.as_ref().is_none_or(|x| c.cost < x.cost) {
                                best = Some(c);
                            }
                        }
                    }
                    return best;
                }
                (true, true) => {}
            }
        }
        let mut q = self.q[a as usize];
        qadd(&mut q, &self.q[b as usize]);
        let mut best = at(a, b, pa);
        let other = at(b, a, pb);
        if other.cost < best.cost {
            best = other;
        }
        let d = sub(pa, pb);
        let slack = 1e-12 * (q[0] + q[4] + q[7]) * dot(d, d);
        if let Some(x) = qsolve(&q) {
            let opt = at(a, b, add(x, self.origin));
            return Some(if best.cost <= opt.cost + slack {
                best
            } else {
                opt
            });
        }
        let mid = at(a, b, scale(add(pa, pb), 0.5));
        Some(if best.cost <= mid.cost + slack {
            best
        } else {
            mid
        })
    }

    fn allowed(&self, c: &Collapse) -> bool {
        let (u, v) = (c.keep, c.gone);
        let shared = self.edge_faces(u, v);
        if shared.is_empty() || shared.len() > 2 {
            return false;
        }
        let mut opposite: Vec<u32> = shared
            .iter()
            .map(|&f| {
                *self.faces[f as usize]
                    .iter()
                    .find(|&&w| w != u && w != v)
                    .unwrap()
            })
            .collect();
        opposite.sort_unstable();
        let (nu, nv) = (self.neighbours(u), self.neighbours(v));
        let common: Vec<u32> = nu
            .iter()
            .copied()
            .filter(|w| nv.binary_search(w).is_ok())
            .collect();
        if common != opposite {
            return false;
        }
        let boundary = |x: u32| !self.boundary_neighbours(x).is_empty();
        if shared.len() == 2 && boundary(u) && boundary(v) {
            return false;
        }
        for &o in &opposite {
            let valence = self.neighbours(o).len();
            if valence <= 2 || (valence <= 3 && !boundary(o)) {
                return false;
            }
        }
        for x in [u, v] {
            for f in self.faces_of(x) {
                if shared.contains(&f) {
                    continue;
                }
                let t = self.faces[f as usize];
                let k = t.iter().position(|&w| w == x).unwrap();
                let a = self.pos[t[(k + 1) % 3] as usize];
                let b = self.pos[t[(k + 2) % 3] as usize];
                let moved = self.pos[x as usize];
                if collinear3(c.target, a, b) {
                    return false;
                }
                let n = normalize(face_normal([moved, a, b]));
                let aux = add(a, n);
                let before = side_of_plane(moved, a, b, aux);
                let after = side_of_plane(c.target, a, b, aux);
                if before == Sign::Zero || before != after {
                    return false;
                }
            }
        }
        true
    }

    fn apply(&mut self, c: &Collapse) -> usize {
        let (u, v) = (c.keep, c.gone);
        let shared = self.edge_faces(u, v);
        for &f in &shared {
            self.face_alive[f as usize] = false;
        }
        let moved: Vec<u32> = self.faces_of(v).collect();
        for f in moved {
            for w in self.faces[f as usize].iter_mut() {
                if *w == v {
                    *w = u;
                }
            }
            self.vf[u as usize].push(f);
        }
        self.vf[v as usize].clear();
        self.vf[u as usize].retain(|&f| self.face_alive[f as usize]);
        self.vf[u as usize].sort_unstable();
        self.vf[u as usize].dedup();
        let qv = self.q[v as usize];
        qadd(&mut self.q[u as usize], &qv);
        self.pos[u as usize] = c.target;
        self.alive[v as usize] = false;
        self.stamp[u as usize] += 1;
        shared.len()
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct Entry {
    cost: u64,
    a: u32,
    b: u32,
    sa: u32,
    sb: u32,
}

pub fn decimate(mesh: &Mesh, opts: &DecimateOptions) -> Result<Mesh> {
    check_finite("mesh vertices", mesh.vertices())?;
    let max_error = opts
        .max_error
        .map(|e| tolerance::check("max_error", e))
        .transpose()?;
    if opts.target_faces.is_none() && max_error.is_none() {
        return Err(Error::invalid("decimate needs target_faces or max_error"));
    }
    check_distinct_corners(mesh.triangles())?;
    HalfEdges::new(mesh)?;
    let origin = mesh.bounds().map_or([0.0; 3], |b| b.min);
    let n = mesh.n_vertices();
    let mut st = State {
        pos: mesh.vertices().to_vec(),
        origin,
        q: vec![[0.0; 10]; n],
        faces: mesh.triangles().to_vec(),
        face_alive: vec![true; mesh.n_faces()],
        vf: vec![vec![]; n],
        alive: vec![true; n],
        stamp: vec![0; n],
    };
    let mut edge_count: std::collections::HashMap<(u32, u32), u32> = Default::default();
    for (f, t) in mesh.triangles().iter().enumerate() {
        let c = t.map(|i| sub(st.pos[i as usize], origin));
        let nf = normalize(face_normal(c));
        let pq = plane_quadric(nf, c[0]);
        for k in 0..3 {
            st.vf[t[k] as usize].push(f as u32);
            qadd(&mut st.q[t[k] as usize], &pq);
            *edge_count
                .entry(crate::geom::ekey(t[k], t[(k + 1) % 3]))
                .or_default() += 1;
        }
    }
    for t in mesh.triangles() {
        let c = t.map(|i| sub(st.pos[i as usize], origin));
        let nf = normalize(face_normal(c));
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if edge_count[&crate::geom::ekey(a, b)] == 1 {
                let along = sub(c[(k + 1) % 3], c[k]);
                let side = normalize(cross(along, nf));
                let pq = plane_quadric(side, c[k]);
                qadd(&mut st.q[a as usize], &pq);
                qadd(&mut st.q[b as usize], &pq);
            }
        }
    }
    let mut heap = BinaryHeap::new();
    let push = |st: &State, heap: &mut BinaryHeap<Reverse<Entry>>, a: u32, b: u32| {
        let (a, b) = crate::geom::ekey(a, b);
        if let Some(c) = st.plan(a, b, opts.preserve_boundary) {
            heap.push(Reverse(Entry {
                cost: c.cost.to_bits(),
                a,
                b,
                sa: st.stamp[a as usize],
                sb: st.stamp[b as usize],
            }));
        }
    };
    let mut keys: Vec<(u32, u32)> = edge_count.keys().copied().collect();
    keys.sort_unstable();
    for (a, b) in keys {
        push(&st, &mut heap, a, b);
    }
    let mut faces_left = mesh.n_faces();
    let target = opts.target_faces.unwrap_or(0);
    let limit = max_error.map(|e| e * e);
    while faces_left > target {
        let Some(Reverse(e)) = heap.pop() else { break };
        if !st.alive[e.a as usize]
            || !st.alive[e.b as usize]
            || st.stamp[e.a as usize] != e.sa
            || st.stamp[e.b as usize] != e.sb
        {
            continue;
        }
        if limit.is_some_and(|l| f64::from_bits(e.cost) > l) {
            break;
        }
        let Some(c) = st.plan(e.a, e.b, opts.preserve_boundary) else {
            continue;
        };
        if !st.allowed(&c) {
            continue;
        }
        faces_left -= st.apply(&c);
        for w in st.neighbours(c.keep) {
            push(&st, &mut heap, c.keep, w);
        }
    }
    let mut new_id = vec![u32::MAX; n];
    let mut kept = vec![];
    for i in 0..n {
        if st.alive[i] {
            new_id[i] = kept.len() as u32;
            kept.push(i as u32);
        }
    }
    let verts = kept.iter().map(|&i| st.pos[i as usize]).collect();
    let live: Vec<u32> = (0..st.faces.len() as u32)
        .filter(|&f| st.face_alive[f as usize])
        .collect();
    let tris = live
        .iter()
        .map(|&f| st.faces[f as usize].map(|v| new_id[v as usize]))
        .collect();
    let out = Mesh::new(verts, tris)?.with_attributes(
        mesh.vertex_attributes().take(&kept),
        mesh.face_attributes().take(&live),
    )?;
    Ok(out
        .with_crs(mesh.crs.clone())
        .with_length_unit(mesh.length_unit.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::close_to_elevation;
    use drusa_core::fixtures::{sphere, torus};
    use drusa_ops::area;
    use proptest::prelude::*;

    fn plate(n: u32) -> Mesh {
        let mut v = vec![];
        for j in 0..n {
            for i in 0..n {
                v.push([
                    500_000.0 + i as f64 * 2.0,
                    7_000_000.0 + j as f64 * 3.0,
                    812.5,
                ]);
            }
        }
        let id = |i: u32, j: u32| j * n + i;
        let mut f = vec![];
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        Mesh::new(v, f).unwrap()
    }

    #[test]
    fn sphere_stays_closed_and_reports_volume_change() {
        let m = sphere(3);
        let max_error = 0.02;
        let d = decimate(
            &m,
            &DecimateOptions {
                target_faces: Some(200),
                max_error: Some(max_error),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(d.n_faces() < m.n_faces());
        let report = d.validate(0.0).unwrap();
        assert!(report.summary.is_solid(), "{:?}", report.summary);
        let (v0, v1) = (
            drusa_ops::volume(&m).unwrap(),
            drusa_ops::volume(&d).unwrap(),
        );
        let change = (v1 - v0).abs();
        println!(
            "sphere decimation {} -> {} faces, volume change {change:e} of {v0}",
            m.n_faces(),
            d.n_faces()
        );
        assert!(
            change <= area(&m) * max_error,
            "volume change {change} exceeds area times max_error"
        );
    }

    #[test]
    fn torus_target_reached_and_closed() {
        let m = torus(32, 16);
        let d = decimate(
            &m,
            &DecimateOptions {
                target_faces: Some(300),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(d.n_faces() <= 300 && d.n_faces() >= 298, "{}", d.n_faces());
        assert!(d.is_closed());
        assert!(d.validate(0.0).unwrap().summary.is_solid());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]
        #[test]
        fn planar_inputs_keep_volume_and_boundary(n in 3u32..10, preserve: bool) {
            let top = plate(n);
            let solid = close_to_elevation(&top, 700.0).unwrap();
            let opts = DecimateOptions { max_error: Some(1e-6), preserve_boundary: preserve, ..Default::default() };
            let d = decimate(&solid, &opts).unwrap();
            prop_assert!(d.validate(0.0).unwrap().summary.is_solid());
            let (v0, v1) = (drusa_ops::volume(&solid).unwrap(), drusa_ops::volume(&d).unwrap());
            prop_assert!((v1 - v0).abs() <= 1e-9 * v0, "{} vs {}", v0, v1);
            prop_assert!(d.n_faces() < solid.n_faces() || n == 2);

            let p = decimate(&top, &opts).unwrap();
            prop_assert!((area(&p) - area(&top)).abs() <= 1e-9 * area(&top));
            let corners = [[0, 0], [n - 1, 0], [0, n - 1], [n - 1, n - 1]];
            for [i, j] in corners {
                let c = top.vertices()[(j * n + i) as usize];
                prop_assert!(p.vertices().contains(&c));
            }
            for v in p.vertices() {
                prop_assert!(top.vertices().contains(v), "vertex {:?} moved", v);
            }
        }
    }

    #[test]
    fn needs_a_stop_rule() {
        assert!(decimate(&sphere(1), &DecimateOptions::default()).is_err());
    }
}
