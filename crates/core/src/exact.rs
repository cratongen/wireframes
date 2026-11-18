//! Exact arithmetic on floating-point expansions (Shewchuk 1997) and implicit
//! points (LPI, TPI) defined by explicit `f64` inputs.
//!
//! Expansions assume no intermediate product underflows or overflows; inputs
//! with magnitudes between about `1e-60` and `1e60` satisfy this for every
//! predicate in this module.

use std::cmp::Ordering;

use crate::predicates::{Sign, kept_axes};

/// A sum of nonoverlapping `f64` components in increasing magnitude, zeros removed.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Expansion(Vec<f64>);

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    (x, (a - av) + (b - bv))
}

fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    (x, b - (x - a))
}

fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    (x, a.mul_add(b, -x))
}

impl Expansion {
    pub fn zero() -> Self {
        Expansion(vec![])
    }

    pub fn from_f64(x: f64) -> Self {
        Expansion(if x == 0.0 { vec![] } else { vec![x] })
    }

    /// `a - b` exactly.
    pub fn diff_f64(a: f64, b: f64) -> Self {
        let (x, y) = two_sum(a, -b);
        Expansion([y, x].into_iter().filter(|&c| c != 0.0).collect())
    }

    pub fn components(&self) -> &[f64] {
        &self.0
    }

    fn grow(mut e: Vec<f64>, b: f64) -> Vec<f64> {
        let mut q = b;
        let mut out = Vec::with_capacity(e.len() + 1);
        for c in e.drain(..) {
            let (s, h) = two_sum(q, c);
            if h != 0.0 {
                out.push(h);
            }
            q = s;
        }
        if q != 0.0 {
            out.push(q);
        }
        out
    }

    pub fn sum(&self, o: &Self) -> Self {
        let (big, small) = if self.0.len() >= o.0.len() {
            (self, o)
        } else {
            (o, self)
        };
        let mut e = big.0.clone();
        for &c in &small.0 {
            e = Self::grow(e, c);
        }
        Expansion(e)
    }

    pub fn neg(&self) -> Self {
        Expansion(self.0.iter().map(|c| -c).collect())
    }

    pub fn diff(&self, o: &Self) -> Self {
        self.sum(&o.neg())
    }

    pub fn scale(&self, b: f64) -> Self {
        if b == 0.0 || self.0.is_empty() {
            return Expansion::zero();
        }
        let mut out = Vec::with_capacity(2 * self.0.len());
        let (mut q, h) = two_product(self.0[0], b);
        if h != 0.0 {
            out.push(h);
        }
        for &c in &self.0[1..] {
            let (t1, t0) = two_product(c, b);
            let (s, h) = two_sum(q, t0);
            if h != 0.0 {
                out.push(h);
            }
            let (s2, h) = fast_two_sum(t1, s);
            if h != 0.0 {
                out.push(h);
            }
            q = s2;
        }
        if q != 0.0 {
            out.push(q);
        }
        Expansion(out)
    }

    pub fn product(&self, o: &Self) -> Self {
        let (big, small) = if self.0.len() >= o.0.len() {
            (self, o)
        } else {
            (o, self)
        };
        small
            .0
            .iter()
            .fold(Expansion::zero(), |acc, &c| acc.sum(&big.scale(c)))
    }

    pub fn sign(&self) -> Sign {
        self.0.last().map_or(Sign::Zero, |&c| Sign::of(c))
    }

    /// The components summed in `f64`, within a few ulps of the exact value.
    pub fn estimate(&self) -> f64 {
        self.0.iter().fold(0.0, |a, b| a + b)
    }
}

type V3 = [Expansion; 3];

fn ex(p: [f64; 3]) -> V3 {
    p.map(Expansion::from_f64)
}

fn sub_f(a: [f64; 3], b: [f64; 3]) -> V3 {
    [0, 1, 2].map(|k| Expansion::diff_f64(a[k], b[k]))
}

fn cross(u: &V3, v: &V3) -> V3 {
    [
        u[1].product(&v[2]).diff(&u[2].product(&v[1])),
        u[2].product(&v[0]).diff(&u[0].product(&v[2])),
        u[0].product(&v[1]).diff(&u[1].product(&v[0])),
    ]
}

fn dot(u: &V3, v: &V3) -> Expansion {
    u[0].product(&v[0])
        .sum(&u[1].product(&v[1]))
        .sum(&u[2].product(&v[2]))
}

fn normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> V3 {
    cross(&sub_f(b, a), &sub_f(c, a))
}

#[derive(Clone, Debug)]
pub enum ImplicitPoint {
    Explicit([f64; 3]),
    /// Line through p, q meets plane through r, s, t.
    Lpi {
        p: [f64; 3],
        q: [f64; 3],
        r: [f64; 3],
        s: [f64; 3],
        t: [f64; 3],
    },
    /// Three planes (each through three explicit points) meet.
    Tpi {
        planes: [[[f64; 3]; 3]; 3],
    },
}

impl ImplicitPoint {
    /// Homogeneous coordinates `(X, Y, Z, W)` with the point at `X / W`.
    pub fn homogeneous(&self) -> (V3, Expansion) {
        match self {
            ImplicitPoint::Explicit(p) => (ex(*p), Expansion::from_f64(1.0)),
            ImplicitPoint::Lpi { p, q, r, s, t } => {
                let n = normal(*r, *s, *t);
                let d = sub_f(*q, *p);
                let w = dot(&n, &d);
                let lam = dot(&n, &sub_f(*r, *p));
                let x = [0, 1, 2].map(|k| w.scale(p[k]).sum(&d[k].product(&lam)));
                (x, w)
            }
            ImplicitPoint::Tpi { planes } => {
                let n = planes.map(|p| normal(p[0], p[1], p[2]));
                let d: Vec<Expansion> = (0..3).map(|i| dot(&n[i], &ex(planes[i][0]))).collect();
                let c12 = cross(&n[1], &n[2]);
                let c20 = cross(&n[2], &n[0]);
                let c01 = cross(&n[0], &n[1]);
                let w = dot(&n[0], &c12);
                let x = [0, 1, 2].map(|k| {
                    d[0].product(&c12[k])
                        .sum(&d[1].product(&c20[k]))
                        .sum(&d[2].product(&c01[k]))
                });
                (x, w)
            }
        }
    }

    /// Denominator nonzero, so the point exists.
    pub fn is_valid(&self) -> bool {
        match self {
            ImplicitPoint::Explicit(_) => true,
            _ => self.homogeneous().1.sign() != Sign::Zero,
        }
    }

    /// Each coordinate rounded to the nearest `f64` (ties to even); NaN when invalid.
    pub fn to_f64(&self) -> [f64; 3] {
        if let ImplicitPoint::Explicit(p) = self {
            return *p;
        }
        let (x, w) = self.homogeneous();
        if w.sign() == Sign::Zero {
            return [f64::NAN; 3];
        }
        let (x, w) = if w.sign() == Sign::Negative {
            (x.map(|e| e.neg()), w.neg())
        } else {
            (x, w)
        };
        [0, 1, 2].map(|k| round_quotient(&x[k], &w))
    }
}

/// Nearest `f64` to `n / d` for `d > 0`.
fn round_quotient(n: &Expansion, d: &Expansion) -> f64 {
    let above = |c: f64| n.diff(&d.scale(c)).sign();
    let mut c = n.estimate() / d.estimate() + 0.0;
    if !c.is_finite() {
        c = 0.0;
    }
    let (lo, hi) = match above(c) {
        Sign::Zero => return c,
        Sign::Positive => {
            let mut lo = c;
            loop {
                let hi = lo.next_up();
                match above(hi) {
                    Sign::Zero => return hi,
                    Sign::Negative => break (lo, hi),
                    Sign::Positive => lo = hi,
                }
            }
        }
        Sign::Negative => {
            let mut hi = c;
            loop {
                let lo = hi.next_down();
                match above(lo) {
                    Sign::Zero => return lo,
                    Sign::Positive => break (lo, hi),
                    Sign::Negative => hi = lo,
                }
            }
        }
    };
    match n.scale(2.0).diff(&d.scale(lo)).diff(&d.scale(hi)).sign() {
        Sign::Negative => lo,
        Sign::Positive => hi,
        Sign::Zero => {
            if lo.to_bits() & 1 == 0 {
                lo
            } else {
                hi
            }
        }
    }
}

/// Same convention as [`crate::predicates::orient3d`] with `a` implicit:
/// the sign of `det[a - d; b - d; c - d]`. `a` must be valid.
pub fn orient3d_i(a: &ImplicitPoint, b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Sign {
    let (x, w) = a.homogeneous();
    let m = cross(&sub_f(b, d), &sub_f(c, d));
    let rel = [0, 1, 2].map(|k| x[k].diff(&w.scale(d[k])));
    dot(&rel, &m).sign().times(w.sign())
}

/// Counter-clockwise orientation of the three points projected by dropping axis
/// `drop` (axes kept in cyclic order, see [`crate::predicates::kept_axes`]).
pub fn orient2d_i(a: &ImplicitPoint, b: &ImplicitPoint, c: &ImplicitPoint, drop: usize) -> Sign {
    let (i, j) = kept_axes(drop);
    let (pa, wa) = a.homogeneous();
    let (pb, wb) = b.homogeneous();
    let (pc, wc) = c.homogeneous();
    let m1 = pb[i].product(&wc).diff(&wb.product(&pc[i]));
    let m2 = pb[j].product(&wc).diff(&wb.product(&pc[j]));
    let m3 = pb[i].product(&pc[j]).diff(&pb[j].product(&pc[i]));
    let det = pa[i]
        .product(&m2)
        .diff(&pa[j].product(&m1))
        .sum(&wa.product(&m3));
    det.sign()
        .times(wa.sign())
        .times(wb.sign())
        .times(wc.sign())
}

/// Order of the two points along `axis`.
pub fn cmp_axis(a: &ImplicitPoint, b: &ImplicitPoint, axis: usize) -> Ordering {
    if let (ImplicitPoint::Explicit(p), ImplicitPoint::Explicit(q)) = (a, b) {
        return p[axis].partial_cmp(&q[axis]).unwrap_or(Ordering::Equal);
    }
    let (pa, wa) = a.homogeneous();
    let (pb, wb) = b.homogeneous();
    let s = pa[axis]
        .product(&wb)
        .diff(&pb[axis].product(&wa))
        .sign()
        .times(wa.sign())
        .times(wb.sign());
    match s {
        Sign::Negative => Ordering::Less,
        Sign::Zero => Ordering::Equal,
        Sign::Positive => Ordering::Greater,
    }
}

/// Lexicographic order by x, then y, then z.
pub fn cmp_lex(a: &ImplicitPoint, b: &ImplicitPoint) -> Ordering {
    (0..3)
        .map(|k| cmp_axis(a, b, k))
        .find(|o| o.is_ne())
        .unwrap_or(Ordering::Equal)
}

pub fn same_point(a: &ImplicitPoint, b: &ImplicitPoint) -> bool {
    cmp_lex(a, b) == Ordering::Equal
}

/// Exact sign of `(p - origin) . normal`.
pub fn plane_side(p: [f64; 3], origin: [f64; 3], normal: [f64; 3]) -> Sign {
    let d = sub_f(p, origin);
    (0..3)
        .fold(Expansion::zero(), |acc, k| acc.sum(&d[k].scale(normal[k])))
        .sign()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oracle::{self, Q, near_degenerate, q, qv};
    use proptest::prelude::*;

    fn value(e: &Expansion) -> Q {
        e.components().iter().fold(q(0.0), |acc, &c| acc + q(c))
    }

    fn oracle_point(p: &ImplicitPoint) -> Option<[Q; 3]> {
        match p {
            ImplicitPoint::Explicit(x) => Some(qv(*x)),
            ImplicitPoint::Lpi { p, q, r, s, t } => oracle::lpi(*p, *q, *r, *s, *t),
            ImplicitPoint::Tpi { planes } => oracle::tpi(*planes),
        }
    }

    fn lpi_strategy() -> impl Strategy<Value = ImplicitPoint> {
        proptest::array::uniform5(near_degenerate())
            .prop_map(|[p, q, r, s, t]| ImplicitPoint::Lpi { p, q, r, s, t })
    }

    fn tpi_strategy() -> impl Strategy<Value = ImplicitPoint> {
        proptest::array::uniform3(proptest::array::uniform3(proptest::array::uniform3(
            -1e3..1e3f64,
        )))
        .prop_map(|planes| ImplicitPoint::Tpi { planes })
    }

    fn any_point() -> impl Strategy<Value = ImplicitPoint> {
        prop_oneof![
            near_degenerate().prop_map(ImplicitPoint::Explicit),
            lpi_strategy(),
            tpi_strategy()
        ]
    }

    fn ordering(s: Sign) -> Ordering {
        match s {
            Sign::Negative => Ordering::Less,
            Sign::Zero => Ordering::Equal,
            Sign::Positive => Ordering::Greater,
        }
    }

    #[test]
    fn expansion_arithmetic_is_exact() {
        let a = Expansion::from_f64(1e100).sum(&Expansion::from_f64(1.0));
        let b = a.diff(&Expansion::from_f64(1e100));
        assert_eq!(b.estimate(), 1.0);
        let p = Expansion::from_f64(1.0 + f64::EPSILON)
            .product(&Expansion::from_f64(1.0 - f64::EPSILON));
        assert_eq!(value(&p), q(1.0 + f64::EPSILON) * q(1.0 - f64::EPSILON));
        assert_eq!(Expansion::diff_f64(0.1, 0.1).sign(), Sign::Zero);
        assert_eq!(Expansion::zero().estimate(), 0.0);
    }

    #[test]
    fn lpi_of_axis_line_and_plane() {
        let p = ImplicitPoint::Lpi {
            p: [0.25, 0.5, -1.0],
            q: [0.25, 0.5, 3.0],
            r: [0.0, 0.0, 1.0],
            s: [1.0, 0.0, 1.0],
            t: [0.0, 1.0, 1.0],
        };
        assert!(p.is_valid());
        assert_eq!(p.to_f64(), [0.25, 0.5, 1.0]);
        assert!(same_point(&p, &ImplicitPoint::Explicit([0.25, 0.5, 1.0])));
        let parallel = ImplicitPoint::Lpi {
            p: [0.0; 3],
            q: [1.0, 0.0, 0.0],
            r: [0.0, 0.0, 1.0],
            s: [1.0, 0.0, 1.0],
            t: [0.0, 1.0, 1.0],
        };
        assert!(!parallel.is_valid());
        assert!(parallel.to_f64()[0].is_nan());
    }

    #[test]
    fn tpi_of_three_axis_planes() {
        let p = ImplicitPoint::Tpi {
            planes: [
                [[2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0]],
                [[0.0, 3.0, 0.0], [0.0, 3.0, 1.0], [1.0, 3.0, 0.0]],
                [[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 1.0, -1.0]],
            ],
        };
        assert_eq!(p.to_f64(), [2.0, 3.0, -1.0]);
    }

    #[test]
    fn plane_side_is_exact() {
        let o = [0.1, 0.2, 0.3];
        assert_eq!(plane_side(o, o, [0.3, 0.7, 0.1]), Sign::Zero);
        assert_eq!(
            plane_side([0.1, 0.2, 0.3f64.next_up()], o, [0.0, 0.0, 1.0]),
            Sign::Positive
        );
        assert_eq!(
            plane_side([1e8, 0.0, 0.0], [1e8 - 1.0, 0.0, 0.0], [-1e-300, 0.0, 0.0]),
            Sign::Negative
        );
    }

    proptest! {
        #[test]
        fn sum_and_product_match_rationals(a in -1e30..1e30f64, b in -1e3..1e3f64, c in -1e-3..1e-3f64, d in -1e10..1e10f64) {
            let x = Expansion::from_f64(a).sum(&Expansion::from_f64(b)).sum(&Expansion::from_f64(c));
            let y = Expansion::diff_f64(d, b);
            prop_assert_eq!(value(&x.product(&y)), (q(a) + q(b) + q(c)) * (q(d) - q(b)));
            prop_assert_eq!(value(&x.scale(d)), (q(a) + q(b) + q(c)) * q(d));
            prop_assert_eq!(value(&x.diff(&y)), q(a) + q(b) + q(c) - (q(d) - q(b)));
        }

        #[test]
        fn orient3d_i_matches_oracle(a in any_point(), b in near_degenerate(), c in near_degenerate(), d in near_degenerate()) {
            prop_assume!(a.is_valid());
            let x = oracle_point(&a).unwrap();
            let want = oracle::sign(&oracle::dot(
                &oracle::sub(&x, &qv(d)),
                &oracle::cross(&oracle::sub(&qv(b), &qv(d)), &oracle::sub(&qv(c), &qv(d))),
            ));
            prop_assert_eq!(orient3d_i(&a, b, c, d), want);
        }

        #[test]
        fn orient3d_i_on_explicit_agrees_with_side_of_plane(a in near_degenerate(), b in near_degenerate(), c in near_degenerate(), d in near_degenerate()) {
            let want = Sign::of(crate::predicates::orient3d(a, b, c, d));
            prop_assert_eq!(orient3d_i(&ImplicitPoint::Explicit(a), b, c, d), want);
        }

        #[test]
        fn cmp_axis_and_orient2d_i_match_oracle(a in any_point(), b in any_point(), c in any_point(), axis in 0usize..3) {
            prop_assume!(a.is_valid() && b.is_valid() && c.is_valid());
            let (x, y, z) = (oracle_point(&a).unwrap(), oracle_point(&b).unwrap(), oracle_point(&c).unwrap());
            prop_assert_eq!(cmp_axis(&a, &b, axis), x[axis].cmp(&y[axis]));
            prop_assert_eq!(cmp_axis(&a, &a, axis), Ordering::Equal);
            let (i, j) = kept_axes(axis);
            let det = (&y[i] - &x[i]) * (&z[j] - &x[j]) - (&y[j] - &x[j]) * (&z[i] - &x[i]);
            prop_assert_eq!(orient2d_i(&a, &b, &c, axis), oracle::sign(&det));
        }

        #[test]
        fn to_f64_rounds_to_nearest(a in prop_oneof![lpi_strategy(), tpi_strategy()]) {
            prop_assume!(a.is_valid());
            let x = oracle_point(&a).unwrap();
            let got = a.to_f64();
            for k in 0..3 {
                prop_assert_eq!(got[k].to_bits(), oracle::nearest(&x[k]).to_bits());
            }
        }

        #[test]
        fn plane_side_matches_oracle(p in near_degenerate(), o in near_degenerate(), n in proptest::array::uniform3(-1e3..1e3f64)) {
            let want = oracle::sign(&oracle::dot(&oracle::sub(&qv(p), &qv(o)), &qv(n)));
            prop_assert_eq!(plane_side(p, o, n), want);
            prop_assert_eq!(ordering(plane_side(p, o, n)), ordering(want));
        }
    }
}
