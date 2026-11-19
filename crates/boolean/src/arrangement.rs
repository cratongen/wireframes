//! Arrangement of two triangle meshes: exact pairwise intersections with
//! implicit points, global exact point identity, and a conforming constrained
//! triangulation of every intersected face.

use std::cmp::Ordering;
use std::collections::HashMap;

use rayon::prelude::*;
use wireframes_core::predicates::{Sign, kept_axes, projection_axis, side_of_plane};
use wireframes_core::tolerance::canonical;
use wireframes_core::{Bvh, Error, Mesh, Result};

use crate::hp::{self, Plane, Pt};
use crate::triangulate::{Geom, cdt};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Loc {
    V(u8),
    E(u8),
    I,
}

pub(crate) struct Face {
    pub c: [[f64; 3]; 3],
    pub v: [u32; 3],
    pub drop: usize,
    pub mesh: u8,
    pub index: u32,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Sub {
    pub v: [u32; 3],
    pub face: u32,
}

pub(crate) struct Arrangement {
    pub pts: Vec<Pt>,
    /// Vertex count of each input; point ids `0..n0` are mesh 0's vertices,
    /// `n0..n0 + n1` mesh 1's, later ids are intersection points.
    pub n_vert: [usize; 2],
    /// The point lies on both inputs.
    pub touched: Vec<bool>,
    pub faces: Vec<Face>,
    pub planes: Vec<Plane>,
    /// Per face, the faces of the other mesh that are coplanar with it and meet it.
    pub coplanar: Vec<Vec<u32>>,
    pub subs: Vec<Sub>,
}

struct PairOut {
    pts: Vec<(Pt, Loc, Loc)>,
    segs: Vec<[usize; 2]>,
    coplanar: bool,
}

fn loc_in(p: &Pt, c: &[[f64; 3]; 3], drop: usize) -> Loc {
    let z = [0, 1, 2]
        .map(|e| hp::orient2d(&Pt::E(c[e]), &Pt::E(c[(e + 1) % 3]), p, drop) == Sign::Zero);
    match (z[0], z[1], z[2]) {
        (true, true, _) => Loc::V(1),
        (_, true, true) => Loc::V(2),
        (true, _, true) => Loc::V(0),
        (true, _, _) => Loc::E(0),
        (_, true, _) => Loc::E(1),
        (_, _, true) => Loc::E(2),
        _ => Loc::I,
    }
}

fn inside_closed(p: &Pt, c: &[[f64; 3]; 3], drop: usize, s: Sign) -> bool {
    (0..3).all(|e| {
        let o = hp::orient2d(&Pt::E(c[e]), &Pt::E(c[(e + 1) % 3]), p, drop);
        o == Sign::Zero || o == s
    })
}

fn orient_of(c: &[[f64; 3]; 3], drop: usize) -> Sign {
    Sign::of(wireframes_core::predicates::orient2d_axis(
        c[0], c[1], c[2], drop,
    ))
}

/// Points of triangle `c` on the plane of the other triangle, given the signs
/// of its corners against that plane.
fn plane_section(c: &[[f64; 3]; 3], s: [Sign; 3], other: &Plane) -> Vec<(Pt, Loc)> {
    let mut out = vec![];
    for i in 0..3 {
        if s[i] == Sign::Zero {
            out.push((Pt::E(c[i]), Loc::V(i as u8)));
        }
        let j = (i + 1) % 3;
        if s[i].times(s[j]) == Sign::Negative {
            out.push((Pt::lpi(c[i], c[j], other), Loc::E(i as u8)));
        }
    }
    out
}

fn pair(t: &Face, g: &Face, pt: &Plane, pg: &Plane) -> Option<PairOut> {
    let st = t.c.map(|p| side_of_plane(p, g.c[0], g.c[1], g.c[2]));
    if st.iter().all(|s| *s == Sign::Zero) {
        return coplanar_pair(t, g);
    }
    if st[0] != Sign::Zero && st.iter().all(|s| *s == st[0]) {
        return None;
    }
    let sg = g.c.map(|p| side_of_plane(p, t.c[0], t.c[1], t.c[2]));
    if sg[0] != Sign::Zero && sg.iter().all(|s| *s == sg[0]) {
        return None;
    }
    let nt = pt.normal_exp();
    let ng = pg.normal_exp();
    let dir: Vec<wireframes_core::exact::Expansion> = (0..3)
        .map(|k| {
            let (i, j) = kept_axes(k);
            nt[i].product(&ng[j]).diff(&nt[j].product(&ng[i]))
        })
        .collect();
    let k = (0..3)
        .filter(|&k| dir[k].sign() != Sign::Zero)
        .max_by(|&a, &b| dir[a].estimate().abs().total_cmp(&dir[b].estimate().abs()))
        .expect("planes are not parallel");
    let sort = |mut v: Vec<(Pt, Loc)>| {
        v.sort_by(|a, b| hp::cmp_axis(&a.0, &b.0, k));
        v
    };
    let a = sort(plane_section(&t.c, st, pg));
    let b = sort(plane_section(&g.c, sg, pt));
    let (a0, a1) = (&a[0], &a[a.len() - 1]);
    let (b0, b1) = (&b[0], &b[b.len() - 1]);
    let full = |x: &(Pt, Loc), from_t: bool, other: Option<&(Pt, Loc)>| -> (Pt, Loc, Loc) {
        match (from_t, other) {
            (true, Some(o)) => (x.0.clone(), x.1, o.1),
            (false, Some(o)) => (x.0.clone(), o.1, x.1),
            (true, None) => (x.0.clone(), x.1, loc_in(&x.0, &g.c, g.drop)),
            (false, None) => (x.0.clone(), loc_in(&x.0, &t.c, t.drop), x.1),
        }
    };
    let lo = match hp::cmp_axis(&a0.0, &b0.0, k) {
        Ordering::Equal => full(a0, true, Some(b0)),
        Ordering::Greater => full(a0, true, None),
        Ordering::Less => full(b0, false, None),
    };
    let hi = match hp::cmp_axis(&a1.0, &b1.0, k) {
        Ordering::Equal => full(a1, true, Some(b1)),
        Ordering::Less => full(a1, true, None),
        Ordering::Greater => full(b1, false, None),
    };
    match hp::cmp_axis(&lo.0, &hi.0, k) {
        Ordering::Greater => None,
        Ordering::Equal => Some(PairOut {
            pts: vec![lo],
            segs: vec![],
            coplanar: false,
        }),
        Ordering::Less => Some(PairOut {
            pts: vec![lo, hi],
            segs: vec![[0, 1]],
            coplanar: false,
        }),
    }
}

fn coplanar_pair(t: &Face, g: &Face) -> Option<PairOut> {
    let drop = t.drop;
    let (s_t, s_g) = (orient_of(&t.c, drop), orient_of(&g.c, drop));
    let mut pts: Vec<(Pt, Loc, Loc)> = vec![];
    for i in 0..3 {
        let p = Pt::E(t.c[i]);
        if inside_closed(&p, &g.c, drop, s_g) {
            let l = loc_in(&p, &g.c, drop);
            pts.push((p, Loc::V(i as u8), l));
        }
    }
    for j in 0..3 {
        let p = Pt::E(g.c[j]);
        if inside_closed(&p, &t.c, drop, s_t) {
            let l = loc_in(&p, &t.c, drop);
            if !matches!(l, Loc::V(_)) {
                pts.push((p, l, Loc::V(j as u8)));
            }
        }
    }
    let o = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| orient_of(&[a, b, c], drop);
    for i in 0..3 {
        let (p0, p1) = (t.c[i], t.c[(i + 1) % 3]);
        for j in 0..3 {
            let (q0, q1) = (g.c[j], g.c[(j + 1) % 3]);
            if o(q0, q1, p0).times(o(q0, q1, p1)) == Sign::Negative
                && o(p0, p1, q0).times(o(p0, p1, q1)) == Sign::Negative
            {
                pts.push((
                    Pt::cross2d(p0, p1, q0, q1, drop),
                    Loc::E(i as u8),
                    Loc::E(j as u8),
                ));
            }
        }
    }
    if pts.is_empty() {
        return None;
    }
    let mut idx: Vec<usize> = (0..pts.len()).collect();
    idx.sort_by(|&a, &b| hp::cmp_proj(&pts[a].0, &pts[b].0, drop));
    let or = |a: usize, b: usize, c: usize| hp::orient2d(&pts[a].0, &pts[b].0, &pts[c].0, drop);
    let mut hull: Vec<usize> = vec![];
    for pass in 0..2 {
        let start = hull.len();
        let seq: Vec<usize> = if pass == 0 {
            idx.clone()
        } else {
            idx.iter().rev().copied().collect()
        };
        for p in seq {
            while hull.len() >= start + 2
                && or(hull[hull.len() - 2], hull[hull.len() - 1], p) != Sign::Positive
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
    }
    let segs = match hull.len() {
        0 | 1 => vec![],
        2 => vec![[hull[0], hull[1]]],
        n => (0..n).map(|i| [hull[i], hull[(i + 1) % n]]).collect(),
    };
    Some(PairOut {
        pts,
        segs,
        coplanar: true,
    })
}

struct FaceGeom<'a> {
    pts: &'a [Pt],
    ids: &'a [u32],
    drop: usize,
}

impl Geom for FaceGeom<'_> {
    fn orient(&self, a: u32, b: u32, c: u32) -> Sign {
        let p = |i: u32| &self.pts[self.ids[i as usize] as usize];
        hp::orient2d(p(a), p(b), p(c), self.drop)
    }

    fn cmp(&self, a: u32, b: u32) -> Ordering {
        let p = |i: u32| &self.pts[self.ids[i as usize] as usize];
        hp::cmp_proj(p(a), p(b), self.drop)
    }

    fn approx(&self, a: u32) -> [f64; 2] {
        let (i, j) = kept_axes(self.drop);
        let p = self.pts[self.ids[a as usize] as usize].approx();
        [p[i], p[j]]
    }
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

fn cmp_pts(a: &Pt, b: &Pt) -> Ordering {
    match (a, b) {
        (Pt::E(x), Pt::E(y)) => x
            .map(canonical)
            .partial_cmp(&y.map(canonical))
            .expect("finite"),
        _ => hp::cmp_lex(a, b),
    }
}

/// Subdivides every face of `a` and `b` along their mutual intersection.
/// Both meshes must have no degenerate faces.
pub(crate) fn build(a: &Mesh, b: &Mesh) -> Result<Arrangement> {
    let n_vert = [a.n_vertices(), b.n_vertices()];
    let mut pts: Vec<Pt> = a
        .vertices()
        .iter()
        .chain(b.vertices())
        .map(|p| Pt::E(*p))
        .collect();
    let mut canon: Vec<u32> = (0..pts.len() as u32).collect();
    let mut first: HashMap<[u64; 3], u32> = HashMap::new();
    for (i, p) in a.vertices().iter().chain(b.vertices()).enumerate() {
        let k = p.map(|x| canonical(x).to_bits());
        canon[i] = *first.entry(k).or_insert(i as u32);
    }
    let mut touched = vec![false; pts.len()];
    for i in n_vert[0]..pts.len() {
        if (canon[i] as usize) < n_vert[0] {
            touched[i] = true;
            touched[canon[i] as usize] = true;
        }
    }
    let mut faces = vec![];
    for (m, mesh) in [a, b].into_iter().enumerate() {
        let off = if m == 0 { 0 } else { n_vert[0] as u32 };
        for (f, t) in mesh.triangles().iter().enumerate() {
            let c = mesh.corners(f as u32);
            let drop = projection_axis(c).ok_or_else(|| {
                Error::invalid(format!(
                    "{} face {f} is degenerate",
                    if m == 0 { "a" } else { "b" }
                ))
            })?;
            faces.push(Face {
                c,
                v: t.map(|i| canon[(i + off) as usize]),
                drop,
                mesh: m as u8,
                index: f as u32,
            });
        }
    }
    let fa = a.n_faces();
    let planes: Vec<Plane> = faces.par_iter().map(|f| Plane::new(f.c)).collect();
    let mut cand = vec![];
    Bvh::new(&a.face_aabbs()).pairs(&Bvh::new(&b.face_aabbs()), |i, j| {
        cand.push((i as usize, fa + j as usize))
    });
    let outs: Vec<(usize, usize, PairOut)> = cand
        .par_iter()
        .filter_map(|&(i, j)| pair(&faces[i], &faces[j], &planes[i], &planes[j]).map(|o| (i, j, o)))
        .collect();

    let mut ids: Vec<Vec<u32>> = Vec::with_capacity(outs.len());
    for (t, g, o) in &outs {
        let v = o
            .pts
            .iter()
            .map(|(p, lt, lg)| match (lt, lg) {
                (Loc::V(k), _) => faces[*t].v[*k as usize],
                (_, Loc::V(k)) => faces[*g].v[*k as usize],
                _ => {
                    pts.push(p.clone());
                    canon.push(canon.len() as u32);
                    touched.push(false);
                    pts.len() as u32 - 1
                }
            })
            .collect();
        ids.push(v);
    }
    let n_input = n_vert[0] + n_vert[1];
    if pts.len() > n_input {
        let mut order: Vec<u32> = (0..n_input as u32)
            .filter(|&i| canon[i as usize] == i)
            .chain(n_input as u32..pts.len() as u32)
            .collect();
        order.sort_by(|&x, &y| cmp_pts(&pts[x as usize], &pts[y as usize]).then(x.cmp(&y)));
        let mut s = 0;
        while s < order.len() {
            let mut e = s + 1;
            while e < order.len()
                && cmp_pts(&pts[order[s] as usize], &pts[order[e] as usize]) == Ordering::Equal
            {
                e += 1;
            }
            let rep = order[s..e].iter().copied().min().expect("non-empty group");
            for &i in &order[s..e] {
                canon[i as usize] = rep;
            }
            s = e;
        }
    }
    for v in ids.iter_mut() {
        for i in v.iter_mut() {
            *i = canon[*i as usize];
            touched[*i as usize] = true;
        }
    }

    let nf = faces.len();
    let mut on_edge: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    let mut extras: Vec<Vec<u32>> = vec![vec![]; nf];
    let mut segs: Vec<Vec<[u32; 2]>> = vec![vec![]; nf];
    let mut coplanar: Vec<Vec<u32>> = vec![vec![]; nf];
    for ((t, g, o), v) in outs.iter().zip(&ids) {
        for (k, (_, lt, lg)) in o.pts.iter().enumerate() {
            for (f, l) in [(*t, lt), (*g, lg)] {
                match l {
                    Loc::V(_) => {}
                    Loc::E(e) => {
                        let fv = faces[f].v;
                        let e = *e as usize;
                        on_edge
                            .entry(ekey(fv[e], fv[(e + 1) % 3]))
                            .or_default()
                            .push(v[k]);
                    }
                    Loc::I => extras[f].push(v[k]),
                }
            }
        }
        for s in &o.segs {
            let (p, q) = (v[s[0]], v[s[1]]);
            if p != q {
                segs[*t].push([p, q]);
                segs[*g].push([p, q]);
            }
        }
        if o.coplanar {
            coplanar[*t].push(*g as u32);
            coplanar[*g].push(*t as u32);
        }
    }
    for (f, face) in faces.iter().enumerate() {
        for e in 0..3 {
            if let Some(list) = on_edge.get(&ekey(face.v[e], face.v[(e + 1) % 3])) {
                extras[f].extend_from_slice(list);
            }
        }
        extras[f].sort_unstable();
        extras[f].dedup();
        extras[f].retain(|p| !face.v.contains(p));
        for s in segs[f].iter_mut() {
            *s = [s[0].min(s[1]), s[0].max(s[1])];
        }
        segs[f].sort_unstable();
        segs[f].dedup();
    }

    let per_face: Vec<Result<Vec<Sub>>> = (0..nf)
        .into_par_iter()
        .map(|f| {
            let face = &faces[f];
            let interior_seg = segs[f]
                .iter()
                .any(|s| !(0..3).any(|e| ekey(face.v[e], face.v[(e + 1) % 3]) == (s[0], s[1])));
            if extras[f].is_empty() && !interior_seg {
                return Ok(vec![Sub {
                    v: face.v,
                    face: f as u32,
                }]);
            }
            let mut local: Vec<u32> = face.v.to_vec();
            local.extend_from_slice(&extras[f]);
            let pos: HashMap<u32, u32> = local
                .iter()
                .enumerate()
                .map(|(i, &g)| (g, i as u32))
                .collect();
            let lsegs: Vec<[u32; 2]> = segs[f].iter().map(|s| [pos[&s[0]], pos[&s[1]]]).collect();
            let geom = FaceGeom {
                pts: &pts,
                ids: &local,
                drop: face.drop,
            };
            let c = cdt(&geom, local.len() as u32, &lsegs).map_err(|e| {
                Error::invalid(format!(
                    "arrangement failed on {} face {}: {e}",
                    if face.mesh == 0 { "a" } else { "b" },
                    face.index
                ))
            })?;
            let flip = orient_of(&face.c, face.drop) == Sign::Negative;
            Ok(c.tris
                .iter()
                .map(|t| {
                    let g = t.map(|i| local[i as usize]);
                    Sub {
                        v: if flip { [g[0], g[2], g[1]] } else { g },
                        face: f as u32,
                    }
                })
                .collect())
        })
        .collect();
    let mut subs = vec![];
    for s in per_face {
        subs.extend(s?);
    }
    Ok(Arrangement {
        pts,
        n_vert,
        touched,
        faces,
        planes,
        coplanar,
        subs,
    })
}

impl Arrangement {
    /// Whether sub-triangle `s` lies inside a coplanar face of the other mesh,
    /// and if so whether that face has the same orientation.
    pub(crate) fn on_other(&self, s: &Sub) -> Option<bool> {
        let f = &self.faces[s.face as usize];
        for &g in &self.coplanar[s.face as usize] {
            let g = &self.faces[g as usize];
            let sg = orient_of(&g.c, f.drop);
            if s.v
                .iter()
                .all(|&p| inside_closed(&self.pts[p as usize], &g.c, f.drop, sg))
            {
                return Some(orient_of(&f.c, f.drop) == sg);
            }
        }
        None
    }

    /// The original input row of a point, as an index into the vertex rows
    /// of `a` followed by those of `b`; `None` for intersection points.
    pub(crate) fn source_row(&self, p: u32) -> Option<u32> {
        ((p as usize) < self.n_vert[0] + self.n_vert[1]).then_some(p)
    }
}
