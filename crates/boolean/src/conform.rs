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

use std::collections::{HashMap, HashSet, VecDeque};

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
    /// Side of each face: the two walls of a slot under closing lie on
    /// opposite sides, and edges and clashes are kept per side.
    side: Vec<bool>,
    edge: HashMap<(u32, u32, bool), usize>,
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
    /// Per input face, whether it is a wall of a slot under closing.
    wall: Vec<bool>,
    /// Commits test faces of different sides unless both are walls, and skip
    /// made faces that cancel a live face with the face they cancel.
    strict: bool,
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
    fn new(m: &Mesh, side: &[bool]) -> Sheet {
        let mut s = Sheet {
            v: m.vertices().to_vec(),
            rows: (0..m.n_vertices() as u32).map(Some).collect(),
            t: vec![],
            alive: vec![],
            src: vec![],
            side: vec![],
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
            wall: vec![],
            strict: false,
        };
        for (i, p) in m.vertices().iter().enumerate() {
            s.at.entry(key(*p)).or_insert(i as u32);
        }
        for (f, t) in m.triangles().iter().enumerate() {
            s.push(*t, f as u32, side.get(f).copied().unwrap_or(false));
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
        let cancels: HashSet<usize> = if self.strict {
            made.iter()
                .filter_map(|&f| self.twin(f).map(|g| [f, g]))
                .flatten()
                .collect()
        } else {
            HashSet::new()
        };
        let wall = |s: &Sheet, f: usize| s.wall.get(s.src[f] as usize).copied().unwrap_or(false);
        let apart = |s: &Sheet, f: usize, g: usize| {
            s.side[g] != s.side[f] && (!s.strict || wall(s, f) && wall(s, g))
        };
        let ok = !made.iter().any(|&f| flat(self.corners(f)))
            && made.iter().all(|&f| {
                if cancels.contains(&f) {
                    return true;
                }
                let c = self.corners(f);
                let mut near = vec![];
                self.boxes
                    .query(&Aabb::from_points(&c), |o| near.push(o as usize));
                near.iter().all(|&o| {
                    self.kids[o].iter().all(|&g| {
                        g == f
                            || !self.alive[g]
                            || apart(self, f, g)
                            || cancels.contains(&g)
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

    /// A live face other than `f` on the same three vertices in the opposite
    /// orientation.
    fn twin(&self, f: usize) -> Option<usize> {
        let [a, b, c] = self.t[f];
        let back = [[a, c, b], [c, b, a], [b, a, c]];
        self.inc[a as usize]
            .iter()
            .copied()
            .find(|&g| g != f && self.alive[g] && back.contains(&self.t[g]))
    }

    fn revive(&mut self, f: usize) {
        self.alive[f] = true;
        self.live += 1;
        let t = self.t[f];
        for k in 0..3 {
            self.edge.insert((t[k], t[(k + 1) % 3], self.side[f]), f);
        }
    }

    fn push(&mut self, t: [u32; 3], src: u32, side: bool) -> usize {
        let f = self.t.len();
        self.t.push(t);
        self.alive.push(true);
        self.src.push(src);
        self.side.push(side);
        self.kids[src as usize].push(f);
        self.log.push((f, true));
        self.live += 1;
        for k in 0..3 {
            self.edge.insert((t[k], t[(k + 1) % 3], side), f);
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
            let e = (t[k], t[(k + 1) % 3], self.side[f]);
            if self.edge.get(&e) == Some(&f) {
                self.edge.remove(&e);
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
                    self.push([t[k], t[(k + 1) % 3], h], self.src[f], self.side[f]);
                }
                Some(tet(c[0], c[1], c[2], p))
            }
            Some(e) => {
                let (a, b, x) = (t[e], t[(e + 1) % 3], t[(e + 2) % 3]);
                let sf = self.side[f];
                let &g = self
                    .edge
                    .get(&(b, a, sf))
                    .or_else(|| self.edge.get(&(b, a, !sf)))?;
                let sg = self.side[g];
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
                self.push([a, h, x], self.src[f], sf);
                self.push([h, b, x], self.src[f], sf);
                self.push([b, h, d], self.src[g], sg);
                self.push([h, a, d], self.src[g], sg);
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
    #[allow(clippy::too_many_arguments)]
    fn recover(
        &mut self,
        u: u32,
        w: u32,
        k: usize,
        plane: [[f64; 3]; 3],
        reach: f64,
        on: bool,
    ) -> Option<f64> {
        let has = |s: &Sheet, a: u32, b: u32| s.edge.contains_key(&(a, b, on));
        if has(self, u, w) || has(self, w, u) {
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
            .filter(|&f| self.alive[f] && self.side[f] == on && self.t[f].contains(&u))
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
            let &g = self.edge.get(&(q, p, on))?;
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
            let (&f1, &f2) = (self.edge.get(&(p, q, on))?, self.edge.get(&(q, p, on))?);
            let (t1, t2) = (self.t[f1], self.t[f2]);
            let r = t1[(0..3).find(|&i| t1[i] != p && t1[i] != q)?];
            let s = t2[(0..3).find(|&i| t2[i] != p && t2[i] != q)?];
            let (vp, vq, vr, vs) = [p, q, r, s].map(|x| self.v[x as usize]).into();
            let convex = Sign::of(orient2d_axis(vr, vs, vp, k))
                .times(Sign::of(orient2d_axis(vr, vs, vq, k)))
                == Sign::Negative;
            if !convex
                || has(self, r, s)
                || has(self, s, r)
                || ![vp, vq, vr, vs].into_iter().all(near)
            {
                crossing.push_back((p, q));
                continue;
            }
            let (sf1, sf2) = (self.src[f1], self.src[f2]);
            self.kill(f1);
            self.kill(f2);
            self.push([r, p, s], sf2, on);
            self.push([s, q, r], sf1, on);
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
    /// Faces of slots thinner than rounding cancelled before conforming.
    pub cancelled: usize,
    pub volume: f64,
    pub steps: usize,
    pub step_bound: usize,
}

/// Conforms `a` and `b`; `None` when nothing changes, or when a conformed
/// copy of a closed input is not a solid or one of an open input crosses
/// itself.
pub(crate) fn conform(a: &Mesh, b: &Mesh) -> Result<Option<Conformed>> {
    let mut s = [Sheet::new(a, &[]), Sheet::new(b, &[])];
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
                if let Some(vol) = s[into].recover(ids[e], ids[(e + 1) % 3], k, c, reach, false) {
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
        cancelled: 0,
        volume: up(volume),
        steps,
        step_bound: bound,
    }))
}

/// Closes the slots of a closed input thinner than rounding: two walls of the
/// same input facing each other within an ulp, which earlier rounding left
/// where two bodies agreed only to rounding (section 9 of
/// `docs/design/rounding.md`). Faces facing opposite ways, one with a corner
/// whose narrow hot box meets the other away from its corners, are paired and
/// the pairs two-coloured; the faces of one colour and the rest of the input
/// are conformed as two sheets ([`conform`]), so where the walls agree to
/// rounding they hold the same faces in opposite orientation, and those cancel.
/// Returns the closed mesh, the input face of each face, the faces cancelled,
/// the conforming counts and the bound on the volume between the input and the
/// result; `None` when there is no slot, the pairs do not two-colour, or the
/// result is not a solid.
fn close_slots(m: &Mesh, strict: bool) -> Result<Option<Closed>> {
    if m.is_empty() || !m.is_closed() {
        return Ok(None);
    }
    let mut at: HashMap<[u64; 3], u32> = HashMap::new();
    let mut first: Vec<u32> = vec![];
    let canon: Vec<u32> = m
        .vertices()
        .iter()
        .enumerate()
        .map(|(i, p)| {
            *at.entry(key(*p)).or_insert_with(|| {
                first.push(i as u32);
                first.len() as u32 - 1
            })
        })
        .collect();
    let v: Vec<[f64; 3]> = first.iter().map(|&i| m.vertices()[i as usize]).collect();
    let tris: Vec<[u32; 3]> = m
        .triangles()
        .iter()
        .map(|t| t.map(|i| canon[i as usize]))
        .collect();
    let corners = |f: usize| tris[f].map(|i| v[i as usize]);
    let normal = |c: [[f64; 3]; 3]| {
        let (a, b) = (
            [0, 1, 2].map(|k| c[1][k] - c[0][k]),
            [0, 1, 2].map(|k| c[2][k] - c[0][k]),
        );
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let flush = |f: [[f64; 3]; 3], g: [[f64; 3]; 3]| {
        let n = normal(f);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        g.iter().all(|&p| {
            let reach = p
                .iter()
                .chain(f.iter().flatten())
                .map(|x| 4.0 * crate::round::half_ulp(*x))
                .fold(0.0, f64::max);
            orient3d(f[0], f[1], f[2], p).abs() <= reach * len
        })
    };
    let bvh = Bvh::new(
        &(0..tris.len())
            .map(|f| Aabb::from_points(&corners(f)))
            .collect::<Vec<_>>(),
    );
    let overlap = |f: [[f64; 3]; 3], g: [[f64; 3]; 3]| {
        let k = dominant(f);
        [(f, g), (g, f)].iter().all(|&(t, u)| {
            let s = Sign::of(orient2d_axis(t[0], t[1], t[2], k));
            (0..3).all(|e| {
                u.iter().any(|&p| {
                    Sign::of(orient2d_axis(t[e], t[(e + 1) % 3], p, k)).times(s) == Sign::Positive
                })
            })
        })
    };
    let mut pairs: Vec<Vec<usize>> = vec![vec![]; tris.len()];
    for f in 0..tris.len() {
        let cf = corners(f);
        let n = normal(cf);
        let (lo, hi) = cf.iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(lo, hi), &p| {
                let (a, b) = hot_box(p, &cf, true);
                (
                    [0, 1, 2].map(|k| lo[k].min(a[k])),
                    [0, 1, 2].map(|k| hi[k].max(b[k])),
                )
            },
        );
        let mut near = vec![];
        bvh.query(&Aabb { min: lo, max: hi }, |g| near.push(g as usize));
        near.sort_unstable();
        for g in near {
            let cg = corners(g);
            if g <= f
                || crate::round::facing(cg, n) != Sign::Negative
                || !flush(cf, cg)
                || !flush(cg, cf)
                || !overlap(cf, cg)
            {
                continue;
            }
            pairs[f].push(g);
            pairs[g].push(f);
        }
    }
    let mut colour: Vec<Option<bool>> = vec![None; tris.len()];
    for start in 0..tris.len() {
        if colour[start].is_some() || pairs[start].is_empty() {
            continue;
        }
        colour[start] = Some(false);
        let mut queue = VecDeque::from([start]);
        while let Some(f) = queue.pop_front() {
            let c = colour[f].expect("coloured");
            for &g in &pairs[f] {
                match colour[g] {
                    None => {
                        colour[g] = Some(!c);
                        queue.push_back(g);
                    }
                    Some(x) if x == c => {
                        return Ok(None);
                    }
                    _ => {}
                }
            }
        }
    }
    if colour.iter().all(|c| c.is_none()) {
        return Ok(None);
    }
    let side: Vec<bool> = colour.iter().map(|c| *c == Some(true)).collect();
    let merged = Mesh::new(v.clone(), tris.clone())?;
    let mut sh = Sheet::new(&merged, &side);
    sh.live = sh.alive.len();
    if strict {
        sh.strict = true;
        sh.wall = colour.iter().map(Option::is_some).collect();
    }
    let bound = step_bound(2 * v.len(), tris.len());
    let (mut inserted, mut flips, mut volume, mut steps) = (0, 0, 0.0, 0);
    let mut done: HashSet<(u32, bool)> = HashSet::new();
    for _ in 0..2 * v.len() + 1 {
        let mut changed = false;
        for into in [false, true] {
            let faces: Vec<usize> = sh
                .live_faces()
                .into_iter()
                .filter(|&f| sh.side[f] == into)
                .collect();
            let bvh = sh.bvh(&faces);
            let mut from: Vec<u32> = sh
                .live_faces()
                .into_iter()
                .filter(|&f| sh.side[f] != into)
                .flat_map(|f| sh.t[f])
                .collect();
            from.sort_unstable();
            from.dedup();
            for h in from {
                if done.contains(&(h, into)) {
                    continue;
                }
                let p = sh.v[h as usize];
                let (lo, hi) = hot_box(p, &[p], true);
                let mut cand = vec![];
                bvh.query(&Aabb { min: lo, max: hi }, |i| cand.push(faces[i as usize]));
                cand.sort_unstable();
                for f in cand {
                    if !sh.alive[f] {
                        continue;
                    }
                    let c = sh.corners(f);
                    if !box_meets(&c, lo, hi) || c.iter().any(|&x| box_meets(&[x], lo, hi)) {
                        continue;
                    }
                    let e = (0..3).find(|&e| box_meets(&[c[e], c[(e + 1) % 3]], lo, hi));
                    if let Some(vol) = sh.insert(p, f, e) {
                        done.insert((h, into));
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
    for f in sh.live_faces() {
        if !sh.alive[f] {
            continue;
        }
        let other = !sh.side[f];
        let uses = |sh: &Sheet, x: u32| {
            sh.inc[x as usize]
                .iter()
                .any(|&g| sh.alive[g] && sh.side[g] == other && sh.t[g].contains(&x))
        };
        let ids = sh.t[f];
        if !ids.iter().all(|&x| uses(&sh, x)) {
            continue;
        }
        let c = sh.corners(f);
        let k = dominant(c);
        let reach = c
            .iter()
            .flatten()
            .map(|x| 8.0 * crate::round::half_ulp(*x))
            .fold(0.0, f64::max);
        for e in 0..3 {
            steps += 1;
            if let Some(vol) = sh.recover(ids[e], ids[(e + 1) % 3], k, c, reach, other) {
                if vol > 0.0 {
                    flips += 1;
                }
                volume += vol;
            }
        }
    }
    steps += sh.attempts;
    assert!(
        steps <= bound,
        "closing slots took {steps} steps, past its bound {bound}"
    );
    let rot = |t: [u32; 3]| {
        let k = (0..3).min_by_key(|&k| t[k]).expect("three corners");
        [t[k], t[(k + 1) % 3], t[(k + 2) % 3]]
    };
    let live = sh.live_faces();
    let mut open: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
    let mut gone: HashSet<usize> = HashSet::new();
    for &f in &live {
        let t = sh.t[f];
        if let Some(list) = open.get_mut(&rot([t[0], t[2], t[1]]))
            && let Some(g) = list.pop()
        {
            gone.insert(f);
            gone.insert(g);
            continue;
        }
        open.entry(rot(t)).or_default().push(f);
    }
    if gone.is_empty() {
        return Ok(None);
    }
    let mut kept: Vec<usize> = live.into_iter().filter(|f| !gone.contains(f)).collect();
    let (voids, held) = thin_voids(&sh.v, &kept.iter().map(|&f| sh.t[f]).collect::<Vec<_>>())?;
    volume += held;
    if !voids.is_empty() {
        let drop: HashSet<usize> = voids.iter().map(|&k| kept[k]).collect();
        kept.retain(|f| !drop.contains(f));
    }
    let mut id: Vec<Option<u32>> = vec![None; sh.v.len()];
    let (mut nv, mut rows) = (vec![], vec![]);
    let t: Vec<[u32; 3]> = kept
        .iter()
        .map(|&f| {
            sh.t[f].map(|i| {
                *id[i as usize].get_or_insert_with(|| {
                    nv.push(sh.v[i as usize]);
                    rows.push(sh.rows[i as usize].map(|r| first[r as usize]));
                    nv.len() as u32 - 1
                })
            })
        })
        .collect();
    let src: Vec<u32> = kept.iter().map(|&f| sh.src[f]).collect();
    let out = Mesh::new(nv, t)?
        .with_attributes(
            m.vertex_attributes().take_or_null(&rows),
            m.face_attributes().take(&src),
        )?
        .with_crs(m.crs.clone())
        .with_length_unit(m.length_unit.clone());
    if !out.validate(0.0)?.summary.is_solid() {
        return Ok(None);
    }
    Ok(Some(Closed {
        mesh: out,
        src,
        cancelled: gone.len(),
        inserted,
        flips,
        volume: up(volume),
        steps,
        step_bound: bound,
    }))
}

/// Positions in `tris` of the faces of each shell that is inward (negative
/// volume, exactly) and holds at most the volume of a sheet `4 rho` thick over
/// its faces, `rho` the rounding radius of a face's farthest corner: a void
/// thinner than rounding, which closing a slot can seal off. Also returns the
/// exact volume those voids hold, rounded up.
fn thin_voids(v: &[[f64; 3]], tris: &[[u32; 3]]) -> Result<(Vec<usize>, f64)> {
    use wireframes_core::EdgeMap;
    use wireframes_core::exact::shell_volume_sign;
    use wireframes_core::topology::Shells;
    let mut id: HashMap<u32, u32> = HashMap::new();
    let mut pts = vec![];
    let local: Vec<[u32; 3]> = tris
        .iter()
        .map(|t| {
            t.map(|i| {
                *id.entry(i).or_insert_with(|| {
                    pts.push(v[i as usize]);
                    pts.len() as u32 - 1
                })
            })
        })
        .collect();
    let m = Mesh::new(pts.clone(), local.clone())?;
    let shells = Shells::new(&m, &EdgeMap::new(&m, None));
    let mut groups: HashMap<u32, Vec<usize>> = HashMap::new();
    for (k, &s) in shells.of.iter().enumerate() {
        groups.entry(s).or_default().push(k);
    }
    let mut keys: Vec<u32> = groups.keys().copied().collect();
    keys.sort_unstable();
    let (mut out, mut total) = (vec![], 0.0);
    for g in keys {
        let ks = &groups[&g];
        let faces: Vec<[u32; 3]> = ks.iter().map(|&k| local[k]).collect();
        if ks.len() == local.len() || shell_volume_sign(&pts, &faces) != Sign::Negative {
            continue;
        }
        let corners: Vec<[[f64; 3]; 3]> =
            faces.iter().map(|t| t.map(|i| pts[i as usize])).collect();
        let held = crate::round::volume_change(&corners, &[]);
        let allowed: f64 = corners
            .iter()
            .map(|c| {
                let rho = c
                    .iter()
                    .map(|p| crate::round::norm(p.map(crate::round::half_ulp)))
                    .fold(0.0, f64::max);
                crate::round::tube(*c, 2.0 * up(rho))
            })
            .sum();
        if held <= allowed {
            out.extend(ks);
            total += held;
        }
    }
    Ok((out, up(total)))
}

/// [`close_slots`] repeated on its own result while it closes a slot, at
/// most once per face of `m`: each pass cancels at least two faces of what
/// the last pass left, and the counter bounds the passes whatever conforming
/// adds. The passes' counts, sources and bounds compose.
fn close_all(m: &Mesh, strict: bool) -> Result<Option<Closed>> {
    let mut acc: Option<Closed> = None;
    for _ in 0..m.n_faces() {
        let Some(next) = close_slots(acc.as_ref().map_or(m, |c| &c.mesh), strict)? else {
            break;
        };
        acc = Some(match acc {
            None => next,
            Some(c) => Closed {
                src: next.src.iter().map(|&f| c.src[f as usize]).collect(),
                cancelled: c.cancelled + next.cancelled,
                inserted: c.inserted + next.inserted,
                flips: c.flips + next.flips,
                volume: up(c.volume + next.volume),
                steps: c.steps + next.steps,
                step_bound: c.step_bound.saturating_add(next.step_bound),
                mesh: next.mesh,
            },
        });
    }
    Ok(acc)
}

/// An input with its thin slots closed: the mesh, the input face of each
/// face, the faces cancelled, the points inserted and edges flipped to make
/// the walls agree, a bound on the volume between the input and the mesh, and
/// the steps taken with their bound.
struct Closed {
    mesh: Mesh,
    src: Vec<u32>,
    cancelled: usize,
    inserted: usize,
    flips: usize,
    volume: f64,
    steps: usize,
    step_bound: usize,
}

/// `a` and `b` with their thin slots closed ([`close_slots`]), then conformed
/// where they agree to rounding ([`conform`]); `None` when neither input has a
/// slot to close.
pub(crate) fn close_and_conform(a: &Mesh, b: &Mesh, strict: bool) -> Result<Option<Conformed>> {
    let closed = [close_all(a, strict)?, close_all(b, strict)?];
    if closed.iter().all(Option::is_none) {
        return Ok(None);
    }
    let inputs = [a, b];
    let meshes: [Mesh; 2] = [0, 1].map(|k| {
        closed[k]
            .as_ref()
            .map_or_else(|| inputs[k].clone(), |c| c.mesh.clone())
    });
    let src: [Vec<u32>; 2] = [0, 1].map(|k| {
        closed[k].as_ref().map_or_else(
            || (0..inputs[k].n_faces() as u32).collect(),
            |c| c.src.clone(),
        )
    });
    let pre = closed.iter().flatten();
    let cancelled: usize = pre.clone().map(|c| c.cancelled).sum();
    let inserted: usize = pre.clone().map(|c| c.inserted).sum();
    let flips: usize = pre.clone().map(|c| c.flips).sum();
    let volume: f64 = pre.clone().map(|c| c.volume).sum();
    let steps: usize = pre.clone().map(|c| c.steps).sum();
    let bound = pre.map(|c| c.step_bound).fold(0, usize::saturating_add);
    Ok(Some(match conform(&meshes[0], &meshes[1])? {
        Some(c) => Conformed {
            src: [0, 1].map(|k| c.src[k].iter().map(|&f| src[k][f as usize]).collect()),
            inserted: c.inserted + inserted,
            flips: c.flips + flips,
            cancelled,
            volume: up(c.volume + volume),
            steps: c.steps + steps,
            step_bound: c.step_bound.saturating_add(bound),
            ..c
        },
        None => Conformed {
            meshes,
            src,
            inserted,
            flips,
            cancelled,
            volume: up(volume),
            steps,
            step_bound: bound,
        },
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::BigRational as Q;

    fn first_input(name: &str) -> Mesh {
        let path = format!(
            "{}/tests/data/attacks/{name}.txt",
            env!("CARGO_MANIFEST_DIR")
        );
        let text = std::fs::read_to_string(path).expect("attack case");
        let mut lines = text.lines();
        let head: Vec<usize> = lines
            .next()
            .expect("header")
            .split(' ')
            .skip(1)
            .map(|x| x.parse().expect("count"))
            .collect();
        let mut row = || -> Vec<f64> {
            lines
                .next()
                .expect("row")
                .split(' ')
                .map(|x| x.parse().expect("number"))
                .collect()
        };
        let v: Vec<[f64; 3]> = (0..head[0])
            .map(|_| {
                let p = row();
                [p[0], p[1], p[2]]
            })
            .collect();
        let t: Vec<[u32; 3]> = (0..head[1])
            .map(|_| {
                let p = row();
                [p[0] as u32, p[1] as u32, p[2] as u32]
            })
            .collect();
        Mesh::new(v, t).expect("mesh")
    }

    fn volume6(m: &Mesh) -> Q {
        let q = |x: f64| Q::from_float(x).expect("finite");
        m.triangles().iter().fold(q(0.0), |acc, t| {
            let [a, b, c] = t.map(|i| m.vertices()[i as usize].map(q));
            acc + &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
                + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
                + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
        })
    }

    /// Stopes rotated to strike 118.4 and dip 62.5, united one at a time,
    /// hold slots between their walls thinner than an ulp. Closing them ends
    /// within its step bound, leaves a solid with fewer faces, and changes the
    /// volume by at most its bound, checked exactly.
    #[test]
    fn closing_the_slots_of_rotated_stopes_ends_within_its_bounds() {
        let m = first_input("stopes_dipping_62_minus_drive");
        let c = close_all(&m, false).unwrap().expect("slots to close");
        assert!(c.steps <= c.step_bound);
        assert!(c.cancelled > 0 && c.mesh.n_faces() < m.n_faces() + 2 * c.inserted);
        assert!(c.mesh.validate(0.0).unwrap().summary.is_solid());
        let diff = volume6(&c.mesh) - volume6(&m);
        let bound = Q::from_float(6.0 * c.volume).unwrap();
        assert!(diff.clone() <= bound && -diff <= bound);
        assert_eq!(c.src.len(), c.mesh.n_faces());
    }

    /// A solid with no walls facing each other within rounding has no slot.
    #[test]
    fn a_box_has_no_slot_to_close() {
        let m = wireframes_core::fixtures::box_mesh(
            [700000.0, 7500000.0, -600.0],
            [700010.0, 7500020.0, -580.0],
        );
        assert!(close_all(&m, false).unwrap().is_none());
    }
}
