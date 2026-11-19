//! Conforming two inputs whose surfaces agree only to rounding, before the
//! exact arrangement (section 8 of `docs/design/rounding.md`). A vertex of one
//! input whose hot box meets a face of the other, not at a corner, is inserted
//! into that face; then each face of one input whose three corners are all
//! vertices of the other, lying within reach of its plane, has its edges
//! recovered in the other by flips in the projection of that face. Where the
//! two surfaces agree to rounding they then share vertices and faces exactly,
//! and the exact arrangement treats them as coincident. No input vertex moves;
//! every insertion and flip replaces faces across a tetrahedron, whose exact
//! volume bounds how far the input changed.

use std::collections::{HashMap, VecDeque};

use wireframes_core::intersect::{SegTri, TriTri, seg_tri, tri_tri};
use wireframes_core::predicates::{Sign, collinear3, orient2d_axis, orient3d};
use wireframes_core::{Aabb, Bvh, Mesh, Result};

use crate::round::{box_meets, hot_box, up};

/// One input under conforming.
struct Sheet {
    v: Vec<[f64; 3]>,
    rows: Vec<Option<u32>>,
    t: Vec<[u32; 3]>,
    alive: Vec<bool>,
    /// Input face each face came from.
    src: Vec<u32>,
    edge: HashMap<(u32, u32), usize>,
    at: HashMap<[u64; 3], u32>,
    inc: Vec<Vec<usize>>,
    /// Current faces per input face, and the input faces' boxes.
    kids: Vec<Vec<usize>>,
    boxes: Bvh,
    /// Faces killed and made since the last commit, to undo.
    log: Vec<(usize, bool)>,
    live: usize,
    /// Flip attempts edge recovery has made.
    attempts: usize,
}

/// Whether faces `a` and `b`, with corners equal by coordinates counted as
/// shared, cross or overlap: `validate`'s test for a pair of faces.
fn clash(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> bool {
    let shared: Vec<(usize, usize)> = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .filter(|&(i, j)| key(a[i]) == key(b[j]))
        .collect();
    let coplanar = b.iter().all(|&p| orient3d(a[0], a[1], a[2], p) == 0.0);
    if coplanar && tri_tri(a, b) == TriTri::CoplanarOverlap {
        return shared.len() < 3;
    }
    match shared[..] {
        [] => tri_tri(a, b) != TriTri::Disjoint,
        [(i, j)] => {
            let opposite = |t: [[f64; 3]; 3], c: usize| (t[(c + 1) % 3], t[(c + 2) % 3]);
            let (p, q) = opposite(a, i);
            let (r, s) = opposite(b, j);
            seg_tri(p, q, b).0 != SegTri::Miss || seg_tri(r, s, a).0 != SegTri::Miss
        }
        _ => false,
    }
}

fn key(p: [f64; 3]) -> [u64; 3] {
    p.map(|x| (x + 0.0).to_bits())
}

/// `|det(b - a, c - a, d - a)| / 6`, the volume of tetrahedron `a b c d`,
/// rounded up.
fn tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    up(orient3d(a, b, c, d).abs() / 6.0)
}

impl Sheet {
    fn new(m: &Mesh) -> Sheet {
        let mut s = Sheet {
            v: m.vertices().to_vec(),
            rows: (0..m.n_vertices() as u32).map(Some).collect(),
            t: vec![],
            alive: vec![],
            src: vec![],
            edge: HashMap::new(),
            at: HashMap::new(),
            inc: vec![vec![]; m.n_vertices()],
            kids: vec![vec![]; m.n_faces()],
            boxes: Bvh::new(
                &(0..m.n_faces() as u32)
                    .map(|f| Aabb::from_points(&m.corners(f)))
                    .collect::<Vec<_>>(),
            ),
            log: vec![],
            live: 0,
            attempts: 0,
        };
        for (i, p) in m.vertices().iter().enumerate() {
            s.at.entry(key(*p)).or_insert(i as u32);
        }
        for (f, t) in m.triangles().iter().enumerate() {
            s.push(*t, f as u32);
        }
        s.log.clear();
        s
    }

    /// Keeps the faces made since the last commit when none is flat and none
    /// crosses or overlaps a live face near it; otherwise undoes them.
    fn commit(&mut self) -> bool {
        let made: Vec<usize> = self
            .log
            .iter()
            .filter(|(f, new)| *new && self.alive[*f])
            .map(|(f, _)| *f)
            .collect();
        let flat = |c: [[f64; 3]; 3]| collinear3(c[0], c[1], c[2]);
        let ok = !made.iter().any(|&f| flat(self.corners(f)))
            && made.iter().all(|&f| {
                let c = self.corners(f);
                let mut near = vec![];
                self.boxes
                    .query(&Aabb::from_points(&c), |o| near.push(o as usize));
                near.iter().all(|&o| {
                    self.kids[o].iter().all(|&g| {
                        g == f
                            || !self.alive[g]
                            || flat(self.corners(g))
                            || !clash(c, self.corners(g))
                    })
                })
            });
        if !ok {
            for (f, new) in std::mem::take(&mut self.log).into_iter().rev() {
                if new {
                    self.kill(f);
                } else {
                    self.revive(f);
                }
            }
        }
        self.log.clear();
        ok
    }

    fn revive(&mut self, f: usize) {
        self.alive[f] = true;
        self.live += 1;
        let t = self.t[f];
        for k in 0..3 {
            self.edge.insert((t[k], t[(k + 1) % 3]), f);
        }
    }

    fn push(&mut self, t: [u32; 3], src: u32) -> usize {
        let f = self.t.len();
        self.t.push(t);
        self.alive.push(true);
        self.src.push(src);
        self.kids[src as usize].push(f);
        self.log.push((f, true));
        self.live += 1;
        for k in 0..3 {
            self.edge.insert((t[k], t[(k + 1) % 3]), f);
            self.inc[t[k] as usize].push(f);
        }
        f
    }

    fn kill(&mut self, f: usize) {
        if self.alive[f] {
            self.log.push((f, false));
            self.live -= 1;
        }
        self.alive[f] = false;
        let t = self.t[f];
        for k in 0..3 {
            if self.edge.get(&(t[k], t[(k + 1) % 3])) == Some(&f) {
                self.edge.remove(&(t[k], t[(k + 1) % 3]));
            }
        }
    }

    fn corners(&self, f: usize) -> [[f64; 3]; 3] {
        self.t[f].map(|i| self.v[i as usize])
    }

    fn vertex(&mut self, p: [f64; 3]) -> u32 {
        if let Some(&i) = self.at.get(&key(p)) {
            return i;
        }
        self.v.push(p);
        self.rows.push(None);
        self.inc.push(vec![]);
        let i = self.v.len() as u32 - 1;
        self.at.insert(key(p), i);
        i
    }

    /// [`Sheet::place`], kept only when [`Sheet::commit`] accepts it.
    fn insert(&mut self, p: [f64; 3], f: usize, e: Option<usize>) -> Option<f64> {
        let fresh = !self.at.contains_key(&key(p));
        let vol = self.place(p, f, e);
        if vol.is_some() && self.commit() {
            return vol;
        }
        self.log.clear();
        if fresh {
            self.at.remove(&key(p));
        }
        None
    }

    /// Inserts point `p` into face `f`, on its edge `e` (splitting the face
    /// across that edge too) or inside it; returns the volume of the
    /// tetrahedra between the old and new faces, or `None` when the edge has
    /// no face across it.
    fn place(&mut self, p: [f64; 3], f: usize, e: Option<usize>) -> Option<f64> {
        let t = self.t[f];
        let c = self.corners(f);
        match e {
            None => {
                let h = self.vertex(p);
                self.kill(f);
                for k in 0..3 {
                    self.push([t[k], t[(k + 1) % 3], h], self.src[f]);
                }
                Some(tet(c[0], c[1], c[2], p))
            }
            Some(e) => {
                let (a, b, x) = (t[e], t[(e + 1) % 3], t[(e + 2) % 3]);
                let &g = self.edge.get(&(b, a))?;
                let u = self.t[g];
                let d = u[(0..3).find(|&k| u[k] != a && u[k] != b)?];
                let h = self.vertex(p);
                let vol = tet(c[0], c[1], c[2], p)
                    + tet(
                        self.v[b as usize],
                        self.v[a as usize],
                        self.v[d as usize],
                        p,
                    );
                self.kill(f);
                self.kill(g);
                self.push([a, h, x], self.src[f]);
                self.push([h, b, x], self.src[f]);
                self.push([b, h, d], self.src[g]);
                self.push([h, a, d], self.src[g]);
                Some(up(vol))
            }
        }
    }

    fn live_faces(&self) -> Vec<usize> {
        (0..self.t.len()).filter(|&f| self.alive[f]).collect()
    }

    fn bvh(&self, faces: &[usize]) -> Bvh {
        let boxes: Vec<Aabb> = faces
            .iter()
            .map(|&f| Aabb::from_points(&self.corners(f)))
            .collect();
        Bvh::new(&boxes)
    }

    fn mesh(&self, like: &Mesh) -> Result<(Mesh, Vec<u32>)> {
        let faces = self.live_faces();
        let mut id: Vec<Option<u32>> = vec![None; self.v.len()];
        let (mut v, mut rows) = (vec![], vec![]);
        let tris: Vec<[u32; 3]> = faces
            .iter()
            .map(|&f| {
                self.t[f].map(|i| {
                    *id[i as usize].get_or_insert_with(|| {
                        v.push(self.v[i as usize]);
                        rows.push(self.rows[i as usize]);
                        v.len() as u32 - 1
                    })
                })
            })
            .collect();
        let src: Vec<u32> = faces.iter().map(|&f| self.src[f]).collect();
        let m = Mesh::new(v, tris)?
            .with_attributes(
                like.vertex_attributes().take_or_null(&rows),
                like.face_attributes().take(&src),
            )?
            .with_crs(like.crs.clone())
            .with_length_unit(like.length_unit.clone());
        Ok((m, src))
    }

    /// Recovers edge `u w` by flipping the edges that cross it in the
    /// projection dropping axis `k` (Sloan 1993), each flip only on a strictly
    /// convex quad whose corners lie within `reach` of `plane`. Returns the
    /// volume of the flipped tetrahedra, or `None` when the walk from `u`
    /// meets a vertex on the segment, leaves the plane, or does not end.
    fn recover(
        &mut self,
        u: u32,
        w: u32,
        k: usize,
        plane: [[f64; 3]; 3],
        reach: f64,
    ) -> Option<f64> {
        if self.edge.contains_key(&(u, w)) || self.edge.contains_key(&(w, u)) {
            return Some(0.0);
        }
        let near = |p: [f64; 3]| {
            let n = {
                let a = [0, 1, 2].map(|i| plane[1][i] - plane[0][i]);
                let b = [0, 1, 2].map(|i| plane[2][i] - plane[0][i]);
                [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ]
            };
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            len > 0.0 && orient3d(plane[0], plane[1], plane[2], p).abs() <= reach * len
        };
        let (pu, pw) = (self.v[u as usize], self.v[w as usize]);
        let side = |p: [f64; 3]| Sign::of(orient2d_axis(pu, pw, p, k));
        let mut crossing: VecDeque<(u32, u32)> = VecDeque::new();
        let mut fan: Vec<usize> = self.inc[u as usize]
            .iter()
            .copied()
            .filter(|&f| self.alive[f] && self.t[f].contains(&u))
            .collect();
        fan.sort_unstable();
        fan.dedup();
        let mut first = None;
        for f in fan {
            let t = self.t[f];
            let i = (0..3).find(|&i| t[i] == u)?;
            let (p, q) = (t[(i + 1) % 3], t[(i + 2) % 3]);
            let (sp, sq) = (side(self.v[p as usize]), side(self.v[q as usize]));
            let ahead = |x: [f64; 3]| {
                (0..3)
                    .map(|a| (x[a] - pu[a]) * (pw[a] - pu[a]))
                    .sum::<f64>()
                    > 0.0
            };
            if sp == Sign::Zero && ahead(self.v[p as usize])
                || sq == Sign::Zero && ahead(self.v[q as usize])
            {
                return None;
            }
            if sp != sq
                && sp != Sign::Zero
                && sq != Sign::Zero
                && Sign::of(orient2d_axis(self.v[p as usize], self.v[q as usize], pu, k))
                    != Sign::of(orient2d_axis(self.v[p as usize], self.v[q as usize], pw, k))
            {
                first = Some((p, q));
                break;
            }
        }
        let (mut p, mut q) = first?;
        let limit = 3 * self.live;
        loop {
            crossing.push_back((p, q));
            if crossing.len() > limit {
                return None;
            }
            let &g = self.edge.get(&(q, p))?;
            let t = self.t[g];
            let r = t[(0..3).find(|&i| t[i] != p && t[i] != q)?];
            if r == w {
                break;
            }
            match side(self.v[r as usize]) {
                Sign::Zero => return None,
                s if s == side(self.v[p as usize]) => p = r,
                _ => q = r,
            }
        }
        let mut vol = 0.0;
        let mut budget = flip_budget(crossing.len());
        while let Some((p, q)) = crossing.pop_front() {
            budget = budget.checked_sub(1)?;
            self.attempts += 1;
            let (&f1, &f2) = (self.edge.get(&(p, q))?, self.edge.get(&(q, p))?);
            let (t1, t2) = (self.t[f1], self.t[f2]);
            let r = t1[(0..3).find(|&i| t1[i] != p && t1[i] != q)?];
            let s = t2[(0..3).find(|&i| t2[i] != p && t2[i] != q)?];
            let (vp, vq, vr, vs) = [p, q, r, s].map(|x| self.v[x as usize]).into();
            let convex = Sign::of(orient2d_axis(vr, vs, vp, k))
                .times(Sign::of(orient2d_axis(vr, vs, vq, k)))
                == Sign::Negative;
            if !convex
                || self.edge.contains_key(&(r, s))
                || self.edge.contains_key(&(s, r))
                || ![vp, vq, vr, vs].into_iter().all(near)
            {
                crossing.push_back((p, q));
                continue;
            }
            let (sf1, sf2) = (self.src[f1], self.src[f2]);
            self.kill(f1);
            self.kill(f2);
            self.push([r, p, s], sf2);
            self.push([s, q, r], sf1);
            if !self.commit() {
                crossing.push_back((p, q));
                continue;
            }
            vol += tet(vp, vq, vr, vs);
            let (a, b) = (side(vr), side(vs));
            let crosses = r != u
                && r != w
                && s != u
                && s != w
                && a != Sign::Zero
                && b != Sign::Zero
                && a != b
                && Sign::of(orient2d_axis(vr, vs, pu, k))
                    .times(Sign::of(orient2d_axis(vr, vs, pw, k)))
                    == Sign::Negative;
            if crosses {
                crossing.push_back((r, s));
            }
        }
        Some(up(vol))
    }
}

/// The attempts edge recovery makes on `c` crossing edges: each flip of a
/// convex quad or retry of a reflex one.
fn flip_budget(c: usize) -> usize {
    c.saturating_mul(c).saturating_add(4)
}

/// The bound on conforming inputs with `v` vertices and `f` faces in all: at
/// most `v` insertions (each vertex enters the other input once), each adding
/// two faces, so at most `F = f + 2 v` live faces; at most `6 F` edge
/// recoveries, each crossing at most `3 F` edges and spending at most
/// [`flip_budget`] of them.
fn step_bound(v: usize, f: usize) -> usize {
    let faces = f.saturating_add(v.saturating_mul(2));
    let per = flip_budget(faces.saturating_mul(3));
    v.saturating_add(faces.saturating_mul(6).saturating_mul(per))
        .saturating_add(1)
}

fn dominant(c: [[f64; 3]; 3]) -> usize {
    (0..3)
        .max_by(|&a, &b| {
            orient2d_axis(c[0], c[1], c[2], a)
                .abs()
                .total_cmp(&orient2d_axis(c[0], c[1], c[2], b).abs())
        })
        .expect("three axes")
}

/// Two inputs conformed where they agree to rounding, with the input face
/// each face came from, how many points were inserted and edges flipped, and
/// a bound on the volume between each input and its conformed copy.
pub(crate) struct Conformed {
    pub meshes: [Mesh; 2],
    pub src: [Vec<u32>; 2],
    pub inserted: usize,
    pub flips: usize,
    pub volume: f64,
    pub steps: usize,
    pub step_bound: usize,
}

/// Conforms `a` and `b`; `None` when nothing changes, or when a conformed
/// copy of a closed input is not a solid or one of an open input crosses
/// itself.
pub(crate) fn conform(a: &Mesh, b: &Mesh) -> Result<Option<Conformed>> {
    let mut s = [Sheet::new(a), Sheet::new(b)];
    for sheet in &mut s {
        sheet.live = sheet.alive.len();
    }
    let bound = step_bound(a.n_vertices() + b.n_vertices(), a.n_faces() + b.n_faces());
    let (mut inserted, mut flips, mut volume, mut steps) = (0, 0, 0.0, 0);
    let pts: [Vec<[f64; 3]>; 2] = [a.vertices().to_vec(), b.vertices().to_vec()];
    let rounds = a.n_vertices() + b.n_vertices() + 1;
    for _ in 0..rounds {
        let mut changed = false;
        for into in 0..2 {
            let from = 1 - into;
            let faces = s[into].live_faces();
            let bvh = s[into].bvh(&faces);
            for &p in &pts[from] {
                if s[into].at.contains_key(&key(p)) {
                    continue;
                }
                let mut cand = vec![];
                let (lo, hi) = hot_box(p, &[], false);
                bvh.query(&Aabb { min: lo, max: hi }, |i| cand.push(faces[i as usize]));
                cand.sort_unstable();
                for f in cand {
                    if !s[into].alive[f] {
                        continue;
                    }
                    let c = s[into].corners(f);
                    if !box_meets(&c, lo, hi) || c.iter().any(|&x| box_meets(&[x], lo, hi)) {
                        continue;
                    }
                    let e = (0..3).find(|&e| box_meets(&[c[e], c[(e + 1) % 3]], lo, hi));
                    if let Some(vol) = s[into].insert(p, f, e) {
                        steps += 1;
                        inserted += 1;
                        volume += vol;
                        changed = true;
                        break;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    for into in 0..2 {
        let from = 1 - into;
        let faces = s[from].live_faces();
        for f in faces {
            if !s[from].alive[f] {
                continue;
            }
            let c = s[from].corners(f);
            let Some(ids) = c
                .iter()
                .map(|p| s[into].at.get(&key(*p)).copied())
                .collect::<Option<Vec<u32>>>()
            else {
                continue;
            };
            let k = dominant(c);
            let reach = c
                .iter()
                .flatten()
                .map(|x| 4.0 * crate::round::half_ulp(*x))
                .fold(0.0, f64::max);
            for e in 0..3 {
                steps += 1;
                if let Some(vol) = s[into].recover(ids[e], ids[(e + 1) % 3], k, c, reach) {
                    if vol > 0.0 {
                        flips += 1;
                    }
                    volume += vol;
                }
            }
        }
    }
    steps += s[0].attempts + s[1].attempts;
    assert!(
        steps <= bound,
        "conforming took {steps} steps, past its bound {bound}"
    );
    if inserted == 0 && flips == 0 {
        return Ok(None);
    }
    let (ma, sa) = s[0].mesh(a)?;
    let (mb, sb) = s[1].mesh(b)?;
    for (m, orig) in [(&ma, a), (&mb, b)] {
        let r = m.validate(0.0)?.summary;
        let crossed = r.self_intersections + r.coplanar_overlaps + r.degenerate_faces > 0;
        if (orig.is_closed() && !r.is_solid()) || crossed {
            return Ok(None);
        }
    }
    Ok(Some(Conformed {
        meshes: [ma, mb],
        src: [sa, sb],
        inserted,
        flips,
        volume: up(volume),
        steps,
        step_bound: bound,
    }))
}
