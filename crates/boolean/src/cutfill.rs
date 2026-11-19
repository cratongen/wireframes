use std::cmp::Ordering;

use rayon::prelude::*;
use wireframes_core::predicates::{Sign, orient2d};
use wireframes_core::{Aabb, Bvh, Error, Mesh, Polylines, Result, units};

use crate::hp::{self, Plane, Pt};
use crate::triangulate::{Geom, cdt, even_odd};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CutFill {
    pub cut: f64,
    pub fill: f64,
    pub net: f64,
    pub cut_area: f64,
    pub fill_area: f64,
    pub area: f64,
}

type P2 = [f64; 2];

fn plan(c: [[f64; 3]; 3]) -> [P2; 3] {
    c.map(|p| [p[0], p[1]])
}

fn plan_sign(t: [P2; 3]) -> Sign {
    Sign::of(orient2d(t[0], t[1], t[2]))
}

fn ccw(t: [P2; 3]) -> [P2; 3] {
    if plan_sign(t) == Sign::Negative {
        [t[0], t[2], t[1]]
    } else {
        t
    }
}

/// Interiors of two non-degenerate plan triangles overlap (exact separating axis test).
fn interiors_overlap(a: [P2; 3], b: [P2; 3]) -> bool {
    let separated = |s: [P2; 3], o: [P2; 3]| {
        (0..3).any(|e| {
            o.iter()
                .all(|&p| Sign::of(orient2d(s[e], s[(e + 1) % 3], p)) != Sign::Positive)
        })
    };
    let (a, b) = (ccw(a), ccw(b));
    !separated(a, b) && !separated(b, a)
}

fn plan_box(t: &[P2]) -> Aabb {
    let p: Vec<[f64; 3]> = t.iter().map(|q| [q[0], q[1], 0.0]).collect();
    Aabb::from_points(&p)
}

fn check_height_field(m: &Mesh, name: &str) -> Result<Vec<bool>> {
    let tris: Vec<[P2; 3]> = (0..m.n_faces() as u32)
        .map(|f| plan(m.corners(f)))
        .collect();
    let flat: Vec<bool> = tris.iter().map(|t| plan_sign(*t) == Sign::Zero).collect();
    let boxes: Vec<Aabb> = tris.iter().map(|t| plan_box(t)).collect();
    let mut pairs = vec![];
    Bvh::new(&boxes).self_pairs(|i, j| pairs.push((i, j)));
    let bad = pairs.par_iter().find_first(|&&(i, j)| {
        !flat[i as usize]
            && !flat[j as usize]
            && interiors_overlap(tris[i as usize], tris[j as usize])
    });
    if let Some((i, j)) = bad {
        return Err(Error::invalid(format!(
            "{name} is not a height field: faces {i} and {j} overlap in plan"
        )));
    }
    Ok(flat)
}

struct Plan2<'a>(&'a [P2]);

impl Geom for Plan2<'_> {
    fn orient(&self, a: u32, b: u32, c: u32) -> Sign {
        Sign::of(orient2d(
            self.0[a as usize],
            self.0[b as usize],
            self.0[c as usize],
        ))
    }

    fn cmp(&self, a: u32, b: u32) -> Ordering {
        let (p, q) = (self.0[a as usize], self.0[b as usize]);
        p.partial_cmp(&q).expect("finite")
    }

    fn approx(&self, a: u32) -> [f64; 2] {
        self.0[a as usize]
    }
}

/// Plan triangles covering the even-odd region of the boundary rings.
fn boundary_triangles(b: &Polylines) -> Result<Vec<[P2; 3]>> {
    let features: Vec<usize> = if b.len() == 1 {
        vec![0]
    } else {
        (0..b.len()).collect()
    };
    let mut pts: Vec<P2> = vec![];
    let mut rings: Vec<Vec<P2>> = vec![];
    for f in features {
        for p in b.feature_parts(f) {
            if b.closed()[p] {
                rings.push(b.part(p).iter().map(|q| [q[0], q[1]]).collect());
            }
        }
    }
    if rings.is_empty() {
        return Err(Error::invalid("boundary has no closed ring"));
    }
    for r in &rings {
        pts.extend_from_slice(r);
    }
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    pts.dedup();
    let id = |p: P2| {
        pts.binary_search_by(|q| q.partial_cmp(&p).expect("finite"))
            .expect("ring vertex") as u32
    };
    let mut segs = vec![];
    for r in &rings {
        for i in 0..r.len() {
            let (a, c) = (id(r[i]), id(r[(i + 1) % r.len()]));
            if a != c {
                segs.push([a, c]);
            }
        }
    }
    let c = cdt(&Plan2(&pts), pts.len() as u32, &segs)
        .map_err(|e| Error::invalid(format!("boundary rings are not simple: {e}")))?;
    Ok(even_odd(&c)
        .into_iter()
        .map(|t| t.map(|i| pts[i as usize]))
        .collect())
}

struct Vert {
    p: Pt,
    /// Supporting line of the edge from this vertex to the next.
    line: (P2, P2),
}

fn p3(p: P2) -> [f64; 3] {
    [p[0], p[1], 0.0]
}

/// Clips a counter-clockwise convex polygon to the closed left side of `a -> b`.
fn clip(poly: Vec<Vert>, a: P2, b: P2) -> Vec<Vert> {
    let n = poly.len();
    if n == 0 {
        return poly;
    }
    let (ea, eb) = (Pt::E(p3(a)), Pt::E(p3(b)));
    let side: Vec<Sign> = poly
        .iter()
        .map(|v| hp::orient2d(&ea, &eb, &v.p, 2))
        .collect();
    let mut out: Vec<Vert> = vec![];
    for i in 0..n {
        let j = (i + 1) % n;
        if side[i] != Sign::Negative {
            let line = if side[i] == Sign::Zero && side[j] == Sign::Negative {
                (a, b)
            } else {
                poly[i].line
            };
            out.push(Vert {
                p: poly[i].p.clone(),
                line,
            });
        }
        if side[i].times(side[j]) == Sign::Negative {
            let (p, q) = poly[i].line;
            out.push(Vert {
                p: Pt::cross2d(p3(p), p3(q), p3(a), p3(b), 2),
                line: if side[i] == Sign::Positive {
                    (a, b)
                } else {
                    poly[i].line
                },
            });
        }
    }
    let mut merged: Vec<Vert> = vec![];
    for v in out {
        match merged.last_mut() {
            Some(last) if hp::cmp_lex(&last.p, &v.p) == Ordering::Equal => last.line = v.line,
            _ => merged.push(v),
        }
    }
    while merged.len() > 1
        && hp::cmp_lex(&merged[0].p, &merged[merged.len() - 1].p) == Ordering::Equal
    {
        merged.pop();
    }
    merged
}

fn tri_poly(t: [P2; 3]) -> Vec<Vert> {
    let t = ccw(t);
    (0..3)
        .map(|k| Vert {
            p: Pt::E(p3(t[k])),
            line: (t[k], t[(k + 1) % 3]),
        })
        .collect()
}

/// Plan area and integral of `d` over a polygon given by f64 vertices and
/// f64 values of `d`, fan-triangulated (the polygon is convex).
fn integrate(xy: &[P2], d: &[f64]) -> (f64, f64) {
    let (mut area, mut int) = (0.0, 0.0);
    for i in 1..xy.len().saturating_sub(1) {
        let (a, b, c) = (xy[0], xy[i], xy[i + 1]);
        let ar = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
        area += ar;
        int += ar * (d[0] + d[i] + d[i + 1]) / 3.0;
    }
    (area, int)
}

fn z_at(c: [[f64; 3]; 3], p: P2) -> f64 {
    let det = (c[1][0] - c[0][0]) * (c[2][1] - c[0][1]) - (c[2][0] - c[0][0]) * (c[1][1] - c[0][1]);
    let l1 =
        ((p[0] - c[0][0]) * (c[2][1] - c[0][1]) - (c[2][0] - c[0][0]) * (p[1] - c[0][1])) / det;
    let l2 =
        ((c[1][0] - c[0][0]) * (p[1] - c[0][1]) - (p[0] - c[0][0]) * (c[1][1] - c[0][1])) / det;
    c[0][2] + l1 * (c[1][2] - c[0][2]) + l2 * (c[2][2] - c[0][2])
}

/// Where `d` reaches 0 between two ends whose exact signs differ. Their f64
/// values can both round to 0, so an equal pair gives the midpoint.
fn crossing(di: f64, dj: f64) -> f64 {
    let den = di - dj;
    if den == 0.0 {
        0.5
    } else {
        (di / den).clamp(0.0, 1.0)
    }
}

/// `[cut, fill, cut_area, fill_area, area]` over one plan cell.
fn cell(
    poly: &[Vert],
    base: &Plane,
    design: &Plane,
    bc: [[f64; 3]; 3],
    dc: [[f64; 3]; 3],
    o: [f64; 3],
) -> [f64; 5] {
    if poly.len() < 3 {
        return [0.0; 5];
    }
    let local = |c: [[f64; 3]; 3]| c.map(|p| [p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
    let (bl, dl) = (local(bc), local(dc));
    let xy: Vec<P2> = poly
        .iter()
        .map(|v| {
            let q = v.p.to_f64();
            [q[0] - o[0], q[1] - o[1]]
        })
        .collect();
    let s: Vec<Sign> = poly
        .iter()
        .map(|v| hp::height_diff_sign(&v.p, design, base))
        .collect();
    let d: Vec<f64> = xy
        .iter()
        .zip(&s)
        .map(|(&p, &sg)| {
            let v = z_at(dl, p) - z_at(bl, p);
            match sg {
                Sign::Zero => 0.0,
                Sign::Positive => v.max(0.0),
                Sign::Negative => v.min(0.0),
            }
        })
        .collect();
    let (total, _) = integrate(&xy, &d);
    let mut out = [0.0, 0.0, 0.0, 0.0, total];
    for (want, slot) in [(Sign::Positive, 1usize), (Sign::Negative, 0usize)] {
        let mut pxy = vec![];
        let mut pd = vec![];
        let n = xy.len();
        for i in 0..n {
            let j = (i + 1) % n;
            if s[i] == want || s[i] == Sign::Zero {
                pxy.push(xy[i]);
                pd.push(d[i]);
            }
            if s[i].times(s[j]) == Sign::Negative {
                let t = crossing(d[i], d[j]);
                pxy.push([
                    xy[i][0] + t * (xy[j][0] - xy[i][0]),
                    xy[i][1] + t * (xy[j][1] - xy[i][1]),
                ]);
                pd.push(0.0);
            }
        }
        let (a, v) = integrate(&pxy, &pd);
        if s.contains(&want) {
            out[slot] = v.abs();
            out[slot + 2] = a;
        }
    }
    out
}

fn neumaier<const K: usize>(items: impl Iterator<Item = [f64; K]>) -> [f64; K] {
    let mut s = [0.0; K];
    let mut c = [0.0; K];
    for x in items {
        for k in 0..K {
            let t = s[k] + x[k];
            c[k] += if s[k].abs() >= x[k].abs() {
                (s[k] - t) + x[k]
            } else {
                (x[k] - t) + s[k]
            };
            s[k] = t;
        }
    }
    std::array::from_fn(|k| s[k] + c[k])
}

/// Volumes between two height-field surfaces over their common plan area,
/// optionally restricted to the even-odd region of `boundary`. `fill` is
/// where `design` is above `base`, `cut` where it is below.
pub fn cut_fill(base: &Mesh, design: &Mesh, boundary: Option<&Polylines>) -> Result<CutFill> {
    let (crs, unit) = units::check_same(
        "design",
        (&base.crs, &base.length_unit),
        (&design.crs, &design.length_unit),
    )?;
    if let Some(b) = boundary {
        units::check_same("boundary", (&crs, &unit), (&b.crs, &b.length_unit))?;
    }
    let base_flat = check_height_field(base, "base")?;
    let design_flat = check_height_field(design, "design")?;
    let clip_tris = boundary.map(boundary_triangles).transpose()?;
    let o = base.bounds().map_or([0.0; 3], |b| b.min);
    let dtris: Vec<[P2; 3]> = (0..design.n_faces() as u32)
        .map(|f| plan(design.corners(f)))
        .collect();
    let dbvh = Bvh::new(&dtris.iter().map(|t| plan_box(t)).collect::<Vec<_>>());
    let cbvh = clip_tris
        .as_ref()
        .map(|c| Bvh::new(&c.iter().map(|t| plan_box(t)).collect::<Vec<_>>()));
    let per_face: Vec<[f64; 5]> = (0..base.n_faces())
        .into_par_iter()
        .map(|f| {
            if base_flat[f] {
                return [0.0; 5];
            }
            let bc = base.corners(f as u32);
            let bt = plan(bc);
            let bp = Plane::new(bc);
            let mut cand = vec![];
            dbvh.query(&plan_box(&bt), |g| cand.push(g));
            cand.sort_unstable();
            let mut parts = vec![];
            for g in cand {
                if design_flat[g as usize] || !interiors_overlap(bt, dtris[g as usize]) {
                    continue;
                }
                let dc = design.corners(g);
                let dt = ccw(dtris[g as usize]);
                let mut poly = tri_poly(bt);
                for e in 0..3 {
                    poly = clip(poly, dt[e], dt[(e + 1) % 3]);
                }
                let dp = Plane::new(dc);
                match (&clip_tris, &cbvh) {
                    (Some(ct), Some(cb)) => {
                        let mut cs = vec![];
                        cb.query(&plan_box(&bt), |k| cs.push(k));
                        cs.sort_unstable();
                        for k in cs {
                            let t = ccw(ct[k as usize]);
                            let mut q: Vec<Vert> = poly
                                .iter()
                                .map(|v| Vert {
                                    p: v.p.clone(),
                                    line: v.line,
                                })
                                .collect();
                            for e in 0..3 {
                                q = clip(q, t[e], t[(e + 1) % 3]);
                            }
                            parts.push(cell(&q, &bp, &dp, bc, dc, o));
                        }
                    }
                    _ => parts.push(cell(&poly, &bp, &dp, bc, dc, o)),
                }
            }
            neumaier(parts.into_iter())
        })
        .collect();
    let [cut, fill, cut_area, fill_area, area] = neumaier(per_face.into_iter());
    Ok(CutFill {
        cut,
        fill,
        net: fill - cut,
        cut_area,
        fill_area,
        area,
    })
}

#[cfg(test)]
mod tests {
    use super::crossing;

    #[test]
    fn crossing_is_finite_when_both_ends_round_to_zero() {
        assert_eq!(crossing(0.0, 0.0), 0.5);
        assert_eq!(crossing(-0.0, 0.0), 0.5);
        assert_eq!(crossing(1.0, -3.0), 0.25);
        assert_eq!(crossing(0.0, -1.0), 0.0);
        assert_eq!(crossing(1e-300, 1e-300), 0.5);
    }
}
