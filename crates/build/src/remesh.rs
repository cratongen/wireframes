use drusa_core::{Error, HalfEdges, Mesh, Result, Table, tolerance};

use drusa_query::SurfaceIndex;

use crate::geom::{add, centroid, check_finite, dist, dot, ekey, face_normal, scale, sub};

struct Work {
    pos: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    face_alive: Vec<bool>,
    vf: Vec<Vec<u32>>,
    alive: Vec<bool>,
}

impl Work {
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

    fn on_boundary(&self, v: u32) -> bool {
        self.neighbours(v)
            .into_iter()
            .any(|w| self.edge_faces(v, w).len() == 1)
    }

    fn edges(&self) -> Vec<(u32, u32)> {
        let mut e: Vec<(u32, u32)> = (0..self.faces.len())
            .filter(|&f| self.face_alive[f])
            .flat_map(|f| {
                let t = self.faces[f];
                [ekey(t[0], t[1]), ekey(t[1], t[2]), ekey(t[2], t[0])]
            })
            .collect();
        e.sort_unstable();
        e.dedup();
        e
    }

    fn len(&self, a: u32, b: u32) -> f64 {
        dist(self.pos[a as usize], self.pos[b as usize])
    }

    fn normal(&self, f: u32) -> [f64; 3] {
        face_normal(self.faces[f as usize].map(|v| self.pos[v as usize]))
    }

    fn split(&mut self, a: u32, b: u32) {
        let m = self.pos.len() as u32;
        self.pos
            .push(scale(add(self.pos[a as usize], self.pos[b as usize]), 0.5));
        self.alive.push(true);
        self.vf.push(vec![]);
        for f in self.edge_faces(a, b) {
            let t = self.faces[f as usize];
            let k = (0..3)
                .find(|&k| {
                    let (x, y) = (t[k], t[(k + 1) % 3]);
                    (x == a && y == b) || (x == b && y == a)
                })
                .unwrap();
            let (x, y, o) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            self.faces[f as usize] = [x, m, o];
            let g = self.faces.len() as u32;
            self.faces.push([m, y, o]);
            self.face_alive.push(true);
            self.vf[y as usize].retain(|&h| h != f);
            self.vf[y as usize].push(g);
            self.vf[o as usize].push(g);
            self.vf[m as usize].extend([f, g]);
        }
    }

    fn flips_ok(&self, x: u32, target: [f64; 3], skip: &[u32]) -> bool {
        self.faces_of(x).filter(|f| !skip.contains(f)).all(|f| {
            let before = self.normal(f);
            let t =
                self.faces[f as usize].map(|v| if v == x { target } else { self.pos[v as usize] });
            let after = face_normal(t);
            dot(before, after) > 0.0 && !drusa_core::predicates::collinear3(t[0], t[1], t[2])
        })
    }

    fn try_collapse(&mut self, a: u32, b: u32, hi: f64) -> bool {
        let shared = self.edge_faces(a, b);
        if shared.len() != 2 {
            return false;
        }
        let (ba, bb) = (self.on_boundary(a), self.on_boundary(b));
        let (keep, gone) = match (ba, bb) {
            (true, true) => return false,
            (true, false) => (a, b),
            (false, true) => (b, a),
            (false, false) => (a.min(b), a.max(b)),
        };
        let target = if ba || bb {
            self.pos[keep as usize]
        } else {
            scale(add(self.pos[a as usize], self.pos[b as usize]), 0.5)
        };
        let mut opposite: Vec<u32> = shared
            .iter()
            .map(|&f| {
                *self.faces[f as usize]
                    .iter()
                    .find(|&&w| w != a && w != b)
                    .unwrap()
            })
            .collect();
        opposite.sort_unstable();
        let (na, nb) = (self.neighbours(a), self.neighbours(b));
        let common: Vec<u32> = na
            .iter()
            .copied()
            .filter(|w| nb.binary_search(w).is_ok())
            .collect();
        if common != opposite {
            return false;
        }
        for &o in &opposite {
            let valence = self.neighbours(o).len();
            if valence <= 3 && !self.on_boundary(o) || valence <= 2 {
                return false;
            }
        }
        if na
            .iter()
            .chain(&nb)
            .any(|&w| w != a && w != b && dist(self.pos[w as usize], target) > hi)
        {
            return false;
        }
        if !self.flips_ok(a, target, &shared) || !self.flips_ok(b, target, &shared) {
            return false;
        }
        for &f in &shared {
            self.face_alive[f as usize] = false;
        }
        let moved: Vec<u32> = self.faces_of(gone).collect();
        for f in moved {
            for w in self.faces[f as usize].iter_mut() {
                if *w == gone {
                    *w = keep;
                }
            }
            self.vf[keep as usize].push(f);
        }
        self.vf[gone as usize].clear();
        let live = &self.face_alive;
        self.vf[keep as usize].retain(|&f| live[f as usize]);
        self.vf[keep as usize].sort_unstable();
        self.vf[keep as usize].dedup();
        self.alive[gone as usize] = false;
        self.pos[keep as usize] = target;
        true
    }

    fn try_flip(&mut self, a: u32, b: u32) -> bool {
        let shared = self.edge_faces(a, b);
        if shared.len() != 2 {
            return false;
        }
        let t = self.faces[shared[0] as usize];
        let (f1, f2) = if (0..3).any(|k| t[k] == a && t[(k + 1) % 3] == b) {
            (shared[0], shared[1])
        } else {
            (shared[1], shared[0])
        };
        let opp = |f: u32| {
            *self.faces[f as usize]
                .iter()
                .find(|&&w| w != a && w != b)
                .unwrap()
        };
        let (c, d) = (opp(f1), opp(f2));
        if c == d || self.neighbours(c).contains(&d) {
            return false;
        }
        let want = |v: u32| if self.on_boundary(v) { 4 } else { 6 };
        let val = |v: u32| self.neighbours(v).len() as i64;
        let dev = |da: i64, db: i64, dc: i64, dd: i64| {
            (val(a) + da - want(a)).abs()
                + (val(b) + db - want(b)).abs()
                + (val(c) + dc - want(c)).abs()
                + (val(d) + dd - want(d)).abs()
        };
        if val(a) <= 3 || val(b) <= 3 || dev(-1, -1, 1, 1) >= dev(0, 0, 0, 0) {
            return false;
        }
        let ref_n = add(self.normal(f1), self.normal(f2));
        let n1 = [c, a, d];
        let n2 = [d, b, c];
        for t in [n1, n2] {
            let p = t.map(|v| self.pos[v as usize]);
            if drusa_core::predicates::collinear3(p[0], p[1], p[2])
                || dot(face_normal(p), ref_n) <= 0.0
            {
                return false;
            }
        }
        self.faces[f1 as usize] = n1;
        self.faces[f2 as usize] = n2;
        self.vf[b as usize].retain(|&f| f != f1);
        self.vf[a as usize].retain(|&f| f != f2);
        self.vf[c as usize].push(f2);
        self.vf[d as usize].push(f1);
        true
    }

    fn vertex_normal(&self, v: u32) -> [f64; 3] {
        crate::geom::normalize(
            self.faces_of(v)
                .fold([0.0; 3], |s, f| add(s, self.normal(f))),
        )
    }
}

pub fn remesh(mesh: &Mesh, edge_length: f64, iterations: usize) -> Result<Mesh> {
    check_finite("mesh vertices", mesh.vertices())?;
    let l = tolerance::check("edge_length", edge_length)?;
    if l == 0.0 {
        return Err(Error::invalid("edge_length must be > 0"));
    }
    if mesh.is_empty() {
        return Err(Error::invalid("mesh must have faces"));
    }
    HalfEdges::new(mesh)?;
    let (hi, lo) = (4.0 / 3.0 * l, 4.0 / 5.0 * l);
    let n = mesh.n_vertices();
    let mut w = Work {
        pos: mesh.vertices().to_vec(),
        faces: mesh.triangles().to_vec(),
        face_alive: vec![true; mesh.n_faces()],
        vf: vec![vec![]; n],
        alive: vec![true; n],
    };
    for (f, t) in mesh.triangles().iter().enumerate() {
        for &v in t {
            w.vf[v as usize].push(f as u32);
        }
    }
    let index = SurfaceIndex::new(mesh);
    for _ in 0..iterations {
        loop {
            let mut long: Vec<(f64, u32, u32)> = w
                .edges()
                .into_iter()
                .map(|(a, b)| (w.len(a, b), a, b))
                .filter(|e| e.0 > hi)
                .collect();
            if long.is_empty() {
                break;
            }
            long.sort_by(|x, y| y.0.total_cmp(&x.0).then((x.1, x.2).cmp(&(y.1, y.2))));
            for (_, a, b) in long {
                if w.edge_faces(a, b).is_empty() || w.len(a, b) <= hi {
                    continue;
                }
                w.split(a, b);
            }
        }
        let mut short: Vec<(f64, u32, u32)> = w
            .edges()
            .into_iter()
            .map(|(a, b)| (w.len(a, b), a, b))
            .filter(|e| e.0 < lo)
            .collect();
        short.sort_by(|x, y| x.0.total_cmp(&y.0).then((x.1, x.2).cmp(&(y.1, y.2))));
        for (_, a, b) in short {
            if w.alive[a as usize] && w.alive[b as usize] && w.len(a, b) < lo {
                w.try_collapse(a, b, hi);
            }
        }
        for (a, b) in w.edges() {
            w.try_flip(a, b);
        }
        let used: Vec<bool> = (0..w.pos.len() as u32)
            .map(|v| w.alive[v as usize] && w.faces_of(v).next().is_some())
            .collect();
        let relaxed: Vec<[f64; 3]> = (0..w.pos.len() as u32)
            .map(|v| {
                let p = w.pos[v as usize];
                if !used[v as usize] || w.on_boundary(v) {
                    return p;
                }
                let nb: Vec<[f64; 3]> =
                    w.neighbours(v).iter().map(|&u| w.pos[u as usize]).collect();
                let d = sub(centroid(&nb), p);
                let nrm = w.vertex_normal(v);
                add(p, sub(d, scale(nrm, dot(nrm, d))))
            })
            .collect();
        let moved: Vec<usize> = (0..w.pos.len())
            .filter(|&v| used[v] && !w.on_boundary(v as u32))
            .collect();
        let targets: Vec<[f64; 3]> = moved.iter().map(|&v| relaxed[v]).collect();
        for (&v, c) in moved.iter().zip(index.closest(&targets)) {
            w.pos[v] = c.point;
        }
    }
    let mut new_id = vec![u32::MAX; w.pos.len()];
    let mut verts = vec![];
    let live: Vec<usize> = (0..w.faces.len()).filter(|&f| w.face_alive[f]).collect();
    for &f in &live {
        for v in w.faces[f] {
            if new_id[v as usize] == u32::MAX {
                new_id[v as usize] = verts.len() as u32;
                verts.push(w.pos[v as usize]);
            }
        }
    }
    let tris: Vec<[u32; 3]> = live
        .iter()
        .map(|&f| w.faces[f].map(|v| new_id[v as usize]))
        .collect();
    let centroids: Vec<[f64; 3]> = live
        .iter()
        .map(|&f| centroid(&w.faces[f].map(|v| w.pos[v as usize])))
        .collect();
    let rows: Vec<u32> = index.closest(&centroids).iter().map(|c| c.face).collect();
    let nv = verts.len();
    let out = Mesh::new(verts, tris)?
        .with_attributes(Table::new(nv), mesh.face_attributes().take(&rows))?;
    Ok(out
        .with_crs(mesh.crs.clone())
        .with_length_unit(mesh.length_unit.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use drusa_core::Column;
    use drusa_core::fixtures::{pit_surface, sphere};

    fn hausdorff_sampled(a: &Mesh, b: &Mesh) -> f64 {
        let ia = SurfaceIndex::new(a);
        let ib = SurfaceIndex::new(b);
        let one = |m: &Mesh, idx: &SurfaceIndex| {
            idx.distance(m.vertices()).into_iter().fold(0.0, f64::max)
        };
        one(a, &ib).max(one(b, &ia))
    }

    #[test]
    fn sphere_stays_closed_and_near_the_input() {
        let m = sphere(2);
        let m = m
            .clone()
            .with_face_column("id", Column::i64((0..m.n_faces() as i64).collect()))
            .unwrap();
        for l in [0.1, 0.3] {
            let r = remesh(&m, l, 5).unwrap();
            assert!(r.is_closed(), "edge length {l}");
            let h = hausdorff_sampled(&m, &r);
            assert!(h < l, "hausdorff {h} for edge length {l}");
            assert_eq!(r.face_attributes().rows(), r.n_faces());
            let lens: Vec<f64> = r
                .triangles()
                .iter()
                .flat_map(|t| (0..3).map(move |k| (t[k], t[(k + 1) % 3])))
                .map(|(a, b)| dist(r.vertices()[a as usize], r.vertices()[b as usize]))
                .collect();
            let mean = lens.iter().sum::<f64>() / lens.len() as f64;
            assert!(mean > 0.5 * l && mean < 1.5 * l, "mean edge {mean} for {l}");
        }
    }

    #[test]
    fn open_surface_keeps_its_corners() {
        let m = pit_surface(8, 1);
        let r = remesh(&m, 6.0, 3).unwrap();
        for c in [[0usize, 0usize], [7, 0], [0, 7], [7, 7]] {
            let p = m.vertices()[c[1] * 8 + c[0]];
            assert!(r.vertices().contains(&p));
        }
        assert!(hausdorff_sampled(&m, &r) < 6.0);
    }

    #[test]
    fn bad_length_errors() {
        assert!(remesh(&sphere(0), 0.0, 1).is_err());
        assert!(remesh(&sphere(0), f64::NAN, 1).is_err());
    }
}
