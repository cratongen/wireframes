use std::cmp::Ordering;

use drusa_core::exact::{Expansion, ImplicitPoint};
use drusa_core::predicates::{
    Sign, kept_axes, orient2d_axis, orient3d, projection_axis, side_of_plane,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HitKind {
    Interior,
    Edge,
    Vertex,
    Coplanar,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub t: f64,
    pub face: u32,
    pub point: [f64; 3],
    pub kind: HitKind,
    /// `Some` for closed meshes, except for `Coplanar` hits, which neither enter nor leave.
    pub entering: Option<bool>,
}

type V = [Expansion; 3];

fn ex(p: [f64; 3]) -> V {
    p.map(Expansion::from_f64)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> V {
    [0, 1, 2].map(|k| Expansion::diff_f64(a[k], b[k]))
}

fn cross(u: &V, v: &V) -> V {
    [0, 1, 2].map(|k| {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        u[i].product(&v[j]).diff(&u[j].product(&v[i]))
    })
}

fn dot(u: &V, v: &V) -> Expansion {
    (0..3).fold(Expansion::zero(), |acc, k| acc.sum(&u[k].product(&v[k])))
}

fn normal(t: &[[f64; 3]; 3]) -> V {
    cross(&sub(t[1], t[0]), &sub(t[2], t[0]))
}

/// A ray `origin + t * dir`, `t >= 0`, or a segment `p + t * (q - p)`, `0 <= t <= 1`.
/// Every decision uses the exact direction (`dir`, or `q - p` as an expansion).
#[derive(Clone, Copy, Debug)]
pub(crate) enum Line {
    Ray([f64; 3], [f64; 3]),
    Segment([f64; 3], [f64; 3]),
}

impl Line {
    fn origin(&self) -> [f64; 3] {
        match *self {
            Line::Ray(p, _) | Line::Segment(p, _) => p,
        }
    }

    fn dir(&self) -> V {
        match *self {
            Line::Ray(_, d) => ex(d),
            Line::Segment(p, q) => sub(q, p),
        }
    }

    /// Sign of `det[u - p, v - p, D]`.
    fn triple(&self, u: [f64; 3], v: [f64; 3]) -> Sign {
        match *self {
            Line::Segment(p, q) => Sign::of(orient3d(u, v, q, p)),
            Line::Ray(p, d) => dot(&cross(&sub(u, p), &sub(v, p)), &ex(d)).sign(),
        }
    }

    pub(crate) fn at(&self, t: f64) -> [f64; 3] {
        match *self {
            Line::Ray(p, d) => [0, 1, 2].map(|k| p[k] + t * d[k]),
            Line::Segment(p, q) => [0, 1, 2].map(|k| p[k] + t * (q[k] - p[k])),
        }
    }
}

/// Exact line parameter `num / den` with `den > 0`.
#[derive(Clone, Debug)]
pub(crate) struct Param {
    num: Expansion,
    den: Expansion,
}

impl Param {
    fn new(num: Expansion, den: Expansion) -> Param {
        if den.sign() == Sign::Negative {
            Param {
                num: num.neg(),
                den: den.neg(),
            }
        } else {
            Param { num, den }
        }
    }

    pub(crate) fn zero() -> Param {
        Param::new(Expansion::zero(), Expansion::from_f64(1.0))
    }

    pub(crate) fn one() -> Param {
        Param::new(Expansion::from_f64(1.0), Expansion::from_f64(1.0))
    }

    pub(crate) fn cmp(&self, o: &Param) -> Ordering {
        match self
            .num
            .product(&o.den)
            .diff(&o.num.product(&self.den))
            .sign()
        {
            Sign::Negative => Ordering::Less,
            Sign::Zero => Ordering::Equal,
            Sign::Positive => Ordering::Greater,
        }
    }

    pub(crate) fn value(&self) -> f64 {
        self.num.estimate() / self.den.estimate()
    }
}

/// Where a line meets one closed triangle: the parameter range `[lo, hi]` (a
/// single parameter unless coplanar) and the sign of `D . normal`.
#[derive(Clone, Debug)]
pub(crate) struct RawHit {
    pub lo: Param,
    pub hi: Param,
    pub face: u32,
    pub kind: HitKind,
    pub crossing: Sign,
    pub corner: Option<u8>,
}

/// Exact intersection of `line` with the closed, non-degenerate triangle `t`.
pub(crate) fn face_hit(t: &[[f64; 3]; 3], face: u32, line: &Line) -> Option<RawHit> {
    let p = line.origin();
    let sp = side_of_plane(p, t[0], t[1], t[2]);
    let mut at_end = false;
    let sw = match *line {
        Line::Segment(_, q) => {
            let sq = side_of_plane(q, t[0], t[1], t[2]);
            if sp == sq && sp != Sign::Zero {
                return None;
            }
            at_end = sq == Sign::Zero && sp != Sign::Zero;
            if sq != Sign::Zero { sq } else { sp.flip() }
        }
        Line::Ray(_, d) => {
            let sw = dot(&normal(t), &ex(d)).sign();
            if sw != Sign::Zero && sp == sw {
                return None;
            }
            if sw == Sign::Zero && sp != Sign::Zero {
                return None;
            }
            sw
        }
    };
    if sw == Sign::Zero {
        return coplanar_hit(t, face, line);
    }
    let s = [0, 1, 2].map(|i| line.triple(t[i], t[(i + 1) % 3]));
    if s.contains(&sw.flip()) {
        return None;
    }
    let zeros: Vec<u8> = (0..3u8).filter(|&i| s[i as usize] == Sign::Zero).collect();
    let (kind, corner) = match zeros[..] {
        [] => (HitKind::Interior, None),
        [_] => (HitKind::Edge, None),
        [e, f] => (HitKind::Vertex, Some(if f == e + 1 { f } else { 0 })),
        _ => unreachable!("the line crosses the plane, so not every triple vanishes"),
    };
    let param = if sp == Sign::Zero {
        Param::zero()
    } else if at_end {
        Param::one()
    } else {
        let n = normal(t);
        Param::new(dot(&n, &sub(t[0], p)), dot(&n, &line.dir()))
    };
    Some(RawHit {
        lo: param.clone(),
        hi: param,
        face,
        kind,
        crossing: sw,
        corner,
    })
}

/// The line lies in the plane of `t`: the parameters where it is in the closed
/// triangle are `s * (a_e + t * b_e) >= 0` for each edge `e`.
fn coplanar_hit(t: &[[f64; 3]; 3], face: u32, line: &Line) -> Option<RawHit> {
    let axis = projection_axis(*t)?;
    let (i, j) = kept_axes(axis);
    let s = Sign::of(orient2d_axis(t[0], t[1], t[2], axis));
    let (p, d) = (line.origin(), line.dir());
    let cross2 = |u: &V, v: &V| u[i].product(&v[j]).diff(&u[j].product(&v[i]));
    let mut lo = Param::zero();
    let mut hi = matches!(line, Line::Segment(..)).then(Param::one);
    for e in 0..3 {
        let (u, v) = (t[e], t[(e + 1) % 3]);
        let uv = sub(v, u);
        let a = cross2(&uv, &sub(p, u));
        let b = cross2(&uv, &d);
        let (sa, sb) = (a.sign().times(s), b.sign().times(s));
        if sb == Sign::Zero {
            if sa == Sign::Negative {
                return None;
            }
            continue;
        }
        let bound = Param::new(a.neg(), b);
        if sb == Sign::Positive {
            if bound.cmp(&lo).is_gt() {
                lo = bound;
            }
        } else if hi.as_ref().is_none_or(|h| bound.cmp(h).is_lt()) {
            hi = Some(bound);
        }
    }
    let hi = hi.unwrap_or_else(|| lo.clone());
    if lo.cmp(&hi).is_gt() {
        return None;
    }
    Some(RawHit {
        lo,
        hi,
        face,
        kind: HitKind::Coplanar,
        crossing: Sign::Zero,
        corner: None,
    })
}

/// The `f64` point of a hit: a corner, the origin, the exactly rounded
/// segment-plane point, or `origin + t * dir` rounded in floating point.
pub(crate) fn hit_point(t: &[[f64; 3]; 3], line: &Line, h: &RawHit) -> [f64; 3] {
    if h.lo.num.sign() == Sign::Zero {
        return line.origin();
    }
    if let Some(v) = h.corner {
        return t[v as usize];
    }
    match (h.kind, line) {
        (HitKind::Interior | HitKind::Edge, Line::Segment(p, q)) => ImplicitPoint::Lpi {
            p: *p,
            q: *q,
            r: t[0],
            s: t[1],
            t: t[2],
        }
        .to_f64(),
        _ => line.at(h.lo.value()),
    }
}

/// Sorts hits by exact parameter then face, and keeps one hit per distinct point.
pub(crate) fn sort_and_collapse(hits: &mut Vec<RawHit>) {
    hits.sort_by(|a, b| a.lo.cmp(&b.lo).then(a.face.cmp(&b.face)));
    hits.dedup_by(|b, a| a.lo.cmp(&b.lo).is_eq());
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];

    fn kind(line: Line) -> Option<HitKind> {
        face_hit(&T, 0, &line).map(|h| h.kind)
    }

    #[test]
    fn rays_classify_interior_edge_vertex_and_miss() {
        let down = [0.0, 0.0, -1.0];
        assert_eq!(
            kind(Line::Ray([0.25, 0.25, 1.0], down)),
            Some(HitKind::Interior)
        );
        assert_eq!(kind(Line::Ray([0.5, 0.0, 1.0], down)), Some(HitKind::Edge));
        assert_eq!(kind(Line::Ray([0.5, 0.5, 1.0], down)), Some(HitKind::Edge));
        assert_eq!(
            kind(Line::Ray([1.0, 0.0, 1.0], down)),
            Some(HitKind::Vertex)
        );
        assert_eq!(kind(Line::Ray([0.6, 0.6, 1.0], down)), None);
        assert_eq!(kind(Line::Ray([0.25, 0.25, 1.0], [0.0, 0.0, 1.0])), None);
        let h = face_hit(&T, 0, &Line::Ray([0.25, 0.25, 2.0], [0.0, 0.0, -0.5])).unwrap();
        assert_eq!(h.lo.value(), 4.0);
        assert_eq!(h.crossing, Sign::Negative);
    }

    #[test]
    fn segments_stop_at_their_ends() {
        let s = |a, b| face_hit(&T, 0, &Line::Segment(a, b));
        assert!(s([0.2, 0.2, 1.0], [0.2, 0.2, 0.5]).is_none());
        let h = s([0.2, 0.2, 1.0], [0.2, 0.2, 0.0]).unwrap();
        assert!(h.lo.cmp(&Param::one()).is_eq());
        let h = s([0.2, 0.2, 0.0], [0.2, 0.2, -1.0]).unwrap();
        assert!(h.lo.cmp(&Param::zero()).is_eq());
        assert_eq!(h.kind, HitKind::Interior);
    }

    #[test]
    fn coplanar_lines_report_their_range() {
        let h = face_hit(&T, 0, &Line::Ray([-1.0, 0.25, 0.0], [1.0, 0.0, 0.0])).unwrap();
        assert_eq!(h.kind, HitKind::Coplanar);
        assert_eq!((h.lo.value(), h.hi.value()), (1.0, 1.75));
        let h = face_hit(&T, 0, &Line::Segment([-1.0, 0.25, 0.0], [3.0, 0.25, 0.0])).unwrap();
        assert_eq!((h.lo.value(), h.hi.value()), (0.25, 0.4375));
        assert!(face_hit(&T, 0, &Line::Ray([-1.0, 2.0, 0.0], [1.0, 0.0, 0.0])).is_none());
        assert!(face_hit(&T, 0, &Line::Ray([2.0, 0.25, 0.0], [1.0, 0.0, 0.0])).is_none());
    }
}
