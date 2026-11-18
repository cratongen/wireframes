//! Patch labelling. A patch is a maximal edge-connected set of sub-triangles
//! of one mesh that does not cross the other mesh's surface. Each patch is
//! labelled from one point strictly inside it by an exact axis-aligned ray.

use std::collections::HashMap;

use drusa_core::predicates::{Sign, kept_axes, orient2d};
use drusa_core::{Aabb, Bvh, Error, Result};

use crate::arrangement::{Arrangement, Sub};
use crate::hp::{self, Plane, Pt};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    In,
    Out,
    OnSame,
    OnOpp,
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

pub(crate) struct Patches {
    /// Patch id per sub-triangle; `u32::MAX` for sub-triangles on the other mesh.
    pub of: Vec<u32>,
    /// Sub-triangle ids of each patch, ascending.
    pub members: Vec<Vec<u32>>,
    /// Coplanar status per sub-triangle.
    pub on: Vec<Option<bool>>,
    /// Undirected edge to sub-triangles, per mesh.
    pub edges: [HashMap<(u32, u32), Vec<u32>>; 2],
}

pub(crate) fn mesh_of(arr: &Arrangement, s: &Sub) -> usize {
    arr.faces[s.face as usize].mesh as usize
}

pub(crate) fn patches(arr: &Arrangement) -> Patches {
    let on: Vec<Option<bool>> = arr.subs.iter().map(|s| arr.on_other(s)).collect();
    let mut edges: [HashMap<(u32, u32), Vec<u32>>; 2] = [HashMap::new(), HashMap::new()];
    for (i, s) in arr.subs.iter().enumerate() {
        let m = mesh_of(arr, s);
        for k in 0..3 {
            edges[m]
                .entry(ekey(s.v[k], s.v[(k + 1) % 3]))
                .or_default()
                .push(i as u32);
        }
    }
    let mut of = vec![u32::MAX; arr.subs.len()];
    let mut members: Vec<Vec<u32>> = vec![];
    for start in 0..arr.subs.len() {
        if on[start].is_some() || of[start] != u32::MAX {
            continue;
        }
        let m = mesh_of(arr, &arr.subs[start]);
        let id = members.len() as u32;
        let mut list = vec![start as u32];
        of[start] = id;
        let mut i = 0;
        while i < list.len() {
            let s = arr.subs[list[i] as usize];
            i += 1;
            for k in 0..3 {
                let e = ekey(s.v[k], s.v[(k + 1) % 3]);
                if edges[1 - m].contains_key(&e) {
                    continue;
                }
                for &n in &edges[m][&e] {
                    if of[n as usize] == u32::MAX && on[n as usize].is_none() {
                        of[n as usize] = id;
                        list.push(n);
                    }
                }
            }
        }
        list.sort_unstable();
        members.push(list);
    }
    Patches {
        of,
        members,
        on,
        edges,
    }
}

/// A point strictly inside a patch, with the axis along which its other two
/// coordinates are explicit `f64` values.
pub(crate) struct Sample {
    pub p: Pt,
    pub axis: usize,
    /// The explicit coordinates of `p` on the kept axes of `axis`.
    pub uv: [f64; 2],
}

impl Sample {
    pub(crate) fn explicit(p: [f64; 3]) -> Sample {
        Sample {
            p: Pt::E(p),
            axis: 2,
            uv: [p[0], p[1]],
        }
    }
}

/// Prefers an input vertex that does not lie on the other mesh, else a point
/// strictly inside one of the patch's sub-triangles, largest first.
pub(crate) fn sample(arr: &Arrangement, members: &[u32]) -> Result<Sample> {
    let n_input = arr.n_vert[0] + arr.n_vert[1];
    for &s in members {
        for &v in &arr.subs[s as usize].v {
            if (v as usize) < n_input && !arr.touched[v as usize] {
                let q = arr.pts[v as usize].explicit().expect("input vertex");
                return Ok(Sample::explicit(q));
            }
        }
    }
    let mut order: Vec<(f64, u32)> = members
        .iter()
        .map(|&s| {
            let sub = &arr.subs[s as usize];
            let drop = arr.faces[sub.face as usize].drop;
            let (i, j) = kept_axes(drop);
            let q = sub.v.map(|v| {
                let a = arr.pts[v as usize].approx();
                [a[i], a[j]]
            });
            (orient2d(q[0], q[1], q[2]).abs(), s)
        })
        .collect();
    order.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    const WEIGHTS: [[f64; 3]; 5] = [
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.5, 0.25, 0.25],
        [0.25, 0.5, 0.25],
        [0.25, 0.25, 0.5],
        [0.6, 0.2, 0.2],
    ];
    for (_, s) in order {
        let sub = &arr.subs[s as usize];
        let face = &arr.faces[sub.face as usize];
        let drop = face.drop;
        let (i, j) = kept_axes(drop);
        let q: Vec<&Pt> = sub.v.iter().map(|&v| &arr.pts[v as usize]).collect();
        let a = q.iter().map(|p| p.approx()).collect::<Vec<_>>();
        let orient = hp::orient2d(q[0], q[1], q[2], drop);
        for w in WEIGHTS {
            let u = w[0] * a[0][i] + w[1] * a[1][i] + w[2] * a[2][i];
            let v = w[0] * a[0][j] + w[1] * a[1][j] + w[2] * a[2][j];
            if !u.is_finite() || !v.is_finite() {
                continue;
            }
            let mut flat = [0.0; 3];
            flat[i] = u;
            flat[j] = v;
            let f = Pt::E(flat);
            if (0..3).all(|k| hp::orient2d(q[k], q[(k + 1) % 3], &f, drop) == orient) {
                return Ok(Sample {
                    p: Pt::on_plane(&arr.planes[sub.face as usize], drop, u, v),
                    axis: drop,
                    uv: [u, v],
                });
            }
        }
    }
    Err(Error::invalid(
        "internal: no f64 point lies strictly inside a patch; it is thinner than f64 resolution",
    ))
}

/// Triangles of one mesh with a box hierarchy, for exact ray queries.
pub(crate) struct Target {
    pub corners: Vec<[[f64; 3]; 3]>,
    bvh: Bvh,
}

pub(crate) struct Hit {
    pub face: u32,
    /// Sign of the face normal's component along the ray axis.
    pub facing: Sign,
}

fn perturbed(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> Sign {
    match Sign::of(orient2d(a, b, p)) {
        Sign::Zero => match Sign::of(a[1] - b[1]) {
            Sign::Zero => Sign::of(b[0] - a[0]),
            s => s,
        },
        s => s,
    }
}

impl Target {
    pub(crate) fn new(corners: Vec<[[f64; 3]; 3]>) -> Target {
        let boxes: Vec<Aabb> = corners.iter().map(|c| Aabb::from_points(c)).collect();
        Target {
            bvh: Bvh::new(&boxes),
            corners,
        }
    }

    /// Faces crossed by the rays `+axis` (first list) and `-axis` (second)
    /// from `p`, with `p` perturbed symbolically in the plane of the other two
    /// axes so that no ray meets an edge or a vertex. `p` must have explicit
    /// coordinates on those two axes and must not lie on any face.
    pub(crate) fn crossings(&self, s: &Sample) -> Result<(Vec<Hit>, Vec<Hit>)> {
        let (p, axis) = (&s.p, s.axis);
        let (i, j) = kept_axes(axis);
        let [pi, pj] = s.uv;
        let mut qb = Aabb {
            min: [0.0; 3],
            max: [0.0; 3],
        };
        qb.min[i] = pi;
        qb.max[i] = pi;
        qb.min[j] = pj;
        qb.max[j] = pj;
        qb.min[axis] = -f64::MAX;
        qb.max[axis] = f64::MAX;
        let mut cand = vec![];
        self.bvh.query(&qb, |f| cand.push(f));
        cand.sort_unstable();
        let (mut up, mut down) = (vec![], vec![]);
        for f in cand {
            let c = self.corners[f as usize];
            let plan = c.map(|v| [v[i], v[j]]);
            let s = Sign::of(orient2d(plan[0], plan[1], plan[2]));
            if s == Sign::Zero {
                continue;
            }
            if (0..3).any(|e| perturbed(plan[e], plan[(e + 1) % 3], [pi, pj]) != s) {
                continue;
            }
            let side = hp::side(p, &Plane::new(c));
            let hit = Hit { face: f, facing: s };
            if side == s.flip() {
                up.push(hit);
            } else if side == s {
                down.push(hit);
            } else {
                return Err(Error::invalid(
                    "internal: classification point lies on the other mesh",
                ));
            }
        }
        Ok((up, down))
    }

    /// Exact containment of a point off the surface of a closed mesh.
    pub(crate) fn inside(&self, s: &Sample) -> Result<bool> {
        Ok(self.crossings(s)?.0.len() % 2 == 1)
    }

    /// The nearest crossed face along `+axis` (`up`) or `-axis`.
    pub(crate) fn nearest(&self, s: &Sample, hits: &[Hit], up: bool) -> Option<usize> {
        let at: Vec<Pt> = hits
            .iter()
            .map(|h| {
                Pt::on_plane(
                    &Plane::new(self.corners[h.face as usize]),
                    s.axis,
                    s.uv[0],
                    s.uv[1],
                )
            })
            .collect();
        (0..hits.len()).reduce(|best, k| {
            let o = hp::cmp_axis(&at[k], &at[best], s.axis);
            let closer = if up { o.is_lt() } else { o.is_gt() };
            if closer { k } else { best }
        })
    }
}

/// Labels every sub-triangle against the other mesh, which must be closed.
pub(crate) fn solid_status(arr: &Arrangement, p: &Patches) -> Result<Vec<Status>> {
    let targets = [0usize, 1].map(|m| {
        Target::new(
            arr.faces
                .iter()
                .filter(|f| f.mesh as usize == m)
                .map(|f| f.c)
                .collect(),
        )
    });
    let labels: Vec<Result<bool>> = {
        use rayon::prelude::*;
        p.members
            .par_iter()
            .map(|members| {
                let m = mesh_of(arr, &arr.subs[members[0] as usize]);
                targets[1 - m].inside(&sample(arr, members)?)
            })
            .collect()
    };
    let labels: Vec<bool> = labels.into_iter().collect::<Result<_>>()?;
    Ok((0..arr.subs.len())
        .map(|s| match p.on[s] {
            Some(true) => Status::OnSame,
            Some(false) => Status::OnOpp,
            None => {
                if labels[p.of[s] as usize] {
                    Status::In
                } else {
                    Status::Out
                }
            }
        })
        .collect())
}
