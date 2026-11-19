//! Exact points for the arrangement: explicit `f64` points and homogeneous
//! points `(X, Y, Z) / W` whose coordinates are compressed floating-point
//! expansions (Shewchuk 1997), with `W > 0`. Every predicate here returns the
//! exact sign; nothing is decided on rounded values.

use std::cmp::Ordering;
use std::sync::Arc;

use wireframes_core::exact::{Expansion, round_quotient};
use wireframes_core::predicates::{Sign, kept_axes, orient2d_axis, side_of_plane};

fn ord(s: Sign) -> Ordering {
    match s {
        Sign::Negative => Ordering::Less,
        Sign::Zero => Ordering::Equal,
        Sign::Positive => Ordering::Greater,
    }
}

#[derive(Debug)]
pub(crate) struct Hp {
    x: [Expansion; 3],
    w: Expansion,
    approx: [f64; 3],
}

#[derive(Clone, Debug)]
pub(crate) enum Pt {
    E([f64; 3]),
    H(Arc<Hp>),
}

/// The plane through three explicit points with its exact normal
/// `(b - a) x (c - a)` and offset `n . a`.
#[derive(Clone, Debug)]
pub(crate) struct Plane {
    pub c: [[f64; 3]; 3],
    n: [Expansion; 3],
    d: Expansion,
}

impl Plane {
    pub(crate) fn new(c: [[f64; 3]; 3]) -> Plane {
        let u = [0, 1, 2].map(|k| Expansion::diff_f64(c[1][k], c[0][k]));
        let v = [0, 1, 2].map(|k| Expansion::diff_f64(c[2][k], c[0][k]));
        let n = [0, 1, 2].map(|k| {
            let (i, j) = kept_axes(k);
            u[i].product(&v[j]).diff(&u[j].product(&v[i]))
        });
        let d = (0..3).fold(Expansion::zero(), |acc, k| acc.sum(&n[k].scale(c[0][k])));
        Plane { c, n, d }
    }

    /// Sign of the exact normal's component `k`.
    #[cfg(test)]
    pub(crate) fn normal_sign(&self, k: usize) -> Sign {
        self.n[k].sign()
    }

    pub(crate) fn normal_exp(&self) -> &[Expansion; 3] {
        &self.n
    }
}

impl Hp {
    fn point(mut x: [Expansion; 3], mut w: Expansion) -> Pt {
        match w.sign() {
            Sign::Zero => panic!("homogeneous point with zero weight"),
            Sign::Negative => {
                w = w.neg();
                x = x.map(|c| c.neg());
            }
            Sign::Positive => {}
        }
        let we = w.estimate();
        let approx = [0, 1, 2].map(|k| x[k].estimate() / we);
        Pt::H(Arc::new(Hp { x, w, approx }))
    }
}

impl Pt {
    fn num(&self, k: usize) -> Expansion {
        match self {
            Pt::E(p) => Expansion::from_f64(p[k]),
            Pt::H(h) => h.x[k].clone(),
        }
    }

    fn den(&self) -> Expansion {
        match self {
            Pt::E(_) => Expansion::from_f64(1.0),
            Pt::H(h) => h.w.clone(),
        }
    }

    pub(crate) fn approx(&self) -> [f64; 3] {
        match self {
            Pt::E(p) => *p,
            Pt::H(h) => h.approx,
        }
    }

    pub(crate) fn explicit(&self) -> Option<[f64; 3]> {
        match self {
            Pt::E(p) => Some(*p),
            Pt::H(_) => None,
        }
    }

    /// Each coordinate rounded to the nearest `f64`.
    pub(crate) fn to_f64(&self) -> [f64; 3] {
        match self {
            Pt::E(p) => *p,
            Pt::H(h) => [0, 1, 2].map(|k| round_quotient(&h.x[k], &h.w)),
        }
    }

    /// The point `p + (n / d) (q - p)`, `d != 0`.
    pub(crate) fn on_line(p: [f64; 3], q: [f64; 3], n: &Expansion, d: &Expansion) -> Pt {
        let x = [0, 1, 2].map(|k| {
            d.scale(p[k])
                .sum(&Expansion::diff_f64(q[k], p[k]).product(n))
        });
        Hp::point(x, d.clone())
    }

    /// Line through `p, q` meets `plane`; the line must not be parallel to it.
    pub(crate) fn lpi(p: [f64; 3], q: [f64; 3], plane: &Plane) -> Pt {
        let dot = |a: [f64; 3]| {
            (0..3).fold(Expansion::zero(), |acc, k| {
                acc.sum(&plane.n[k].product(&Expansion::diff_f64(a[k], p[k])))
            })
        };
        let d = dot(q);
        let n = dot(plane.c[0]);
        Pt::on_line(p, q, &n, &d)
    }

    /// Proper crossing of coplanar segments `a0 a1` and `b0 b1`, decided in the
    /// projection dropping `drop`.
    pub(crate) fn cross2d(
        a0: [f64; 3],
        a1: [f64; 3],
        b0: [f64; 3],
        b1: [f64; 3],
        drop: usize,
    ) -> Pt {
        let (i, j) = kept_axes(drop);
        let o = |p: [f64; 3]| {
            Expansion::diff_f64(b1[i], b0[i])
                .product(&Expansion::diff_f64(p[j], b0[j]))
                .diff(&Expansion::diff_f64(b1[j], b0[j]).product(&Expansion::diff_f64(p[i], b0[i])))
        };
        let (o0, o1) = (o(a0), o(a1));
        Pt::on_line(a0, a1, &o0, &o0.diff(&o1))
    }

    /// The point of `plane` whose coordinates on the kept axes of `drop` are
    /// `(u, v)`; the normal's component `drop` must be nonzero.
    pub(crate) fn on_plane(plane: &Plane, drop: usize, u: f64, v: f64) -> Pt {
        let (i, j) = kept_axes(drop);
        let n = &plane.n;
        let mut x = [Expansion::zero(), Expansion::zero(), Expansion::zero()];
        x[i] = n[drop].scale(u);
        x[j] = n[drop].scale(v);
        x[drop] = plane.d.diff(&n[i].scale(u)).diff(&n[j].scale(v));
        Hp::point(x, n[drop].clone())
    }

    /// Segment `a b` meets the plane `(p - origin) . normal = 0`; `a` and `b`
    /// must lie strictly on opposite sides.
    pub(crate) fn plane_cross(a: [f64; 3], b: [f64; 3], origin: [f64; 3], normal: [f64; 3]) -> Pt {
        let dot = |p: [f64; 3]| {
            (0..3).fold(Expansion::zero(), |acc, k| {
                acc.sum(&Expansion::diff_f64(p[k], a[k]).scale(normal[k]))
            })
        };
        Pt::on_line(a, b, &dot(origin), &dot(b))
    }
}

/// Sign of the projected orientation of `a, b, c` dropping axis `drop`: the
/// sign of component `drop` of `(b - a) x (c - a)`.
pub(crate) fn orient2d(a: &Pt, b: &Pt, c: &Pt, drop: usize) -> Sign {
    if let (Pt::E(a), Pt::E(b), Pt::E(c)) = (a, b, c) {
        return Sign::of(orient2d_axis(*a, *b, *c, drop));
    }
    let (i, j) = kept_axes(drop);
    let r = [a, b, c].map(|p| [p.num(i), p.num(j), p.den()]);
    let m = |p: usize, q: usize, s: usize, t: usize| {
        r[1][p].product(&r[2][q]).diff(&r[1][s].product(&r[2][t]))
    };
    r[0][0]
        .product(&m(1, 2, 2, 1))
        .diff(&r[0][1].product(&m(0, 2, 2, 0)))
        .sum(&r[0][2].product(&m(0, 1, 1, 0)))
        .sign()
}

pub(crate) fn cmp_axis(a: &Pt, b: &Pt, k: usize) -> Ordering {
    if let (Pt::E(a), Pt::E(b)) = (a, b) {
        return a[k].partial_cmp(&b[k]).expect("finite coordinates");
    }
    ord(a
        .num(k)
        .product(&b.den())
        .diff(&b.num(k).product(&a.den()))
        .sign())
}

pub(crate) fn cmp_lex(a: &Pt, b: &Pt) -> Ordering {
    (0..3)
        .map(|k| cmp_axis(a, b, k))
        .find(|o| *o != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

/// Lexicographic order on the kept axes of `drop`.
pub(crate) fn cmp_proj(a: &Pt, b: &Pt, drop: usize) -> Ordering {
    let (i, j) = kept_axes(drop);
    cmp_axis(a, b, i).then_with(|| cmp_axis(a, b, j))
}

/// Which side of `plane` holds `p`; `Positive` is the side its normal points to.
pub(crate) fn side(p: &Pt, plane: &Plane) -> Sign {
    match p {
        Pt::E(q) => side_of_plane(*q, plane.c[0], plane.c[1], plane.c[2]),
        Pt::H(h) => {
            let dot = (0..3).fold(Expansion::zero(), |acc, k| {
                acc.sum(&plane.n[k].product(&h.x[k]))
            });
            dot.diff(&plane.d.product(&h.w)).sign()
        }
    }
}

/// Sign of `z_top(x, y) - z_bottom(x, y)` at the plan position of `p`, for
/// two non-vertical planes.
pub(crate) fn height_diff_sign(p: &Pt, top: &Plane, bottom: &Plane) -> Sign {
    let (x, y, w) = (p.num(0), p.num(1), p.den());
    let scaled = |pl: &Plane| {
        pl.d.product(&w)
            .diff(&pl.n[0].product(&x))
            .diff(&pl.n[1].product(&y))
    };
    let s = scaled(top)
        .product(&bottom.n[2])
        .diff(&scaled(bottom).product(&top.n[2]))
        .sign();
    s.times(top.n[2].sign()).times(bottom.n[2].sign())
}

/// `x = m 2^e` with an integer `m`.
fn split_f64(x: f64) -> (bool, u64, i32) {
    let b = x.to_bits();
    let (exp, frac) = (((b >> 52) & 0x7ff) as i32, b & ((1 << 52) - 1));
    let (m, e) = if exp == 0 {
        (frac, -1074)
    } else {
        (frac | (1 << 52), exp - 1075)
    };
    (x < 0.0, m, e)
}

fn add_at(acc: &mut [u64], mut i: usize, v: u128) {
    let mut carry = v;
    while carry != 0 {
        let s = acc[i] as u128 + (carry & u64::MAX as u128);
        acc[i] = s as u64;
        carry = (carry >> 64) + (s >> 64);
        i += 1;
    }
}

/// Exact sign of `sum x_k y_k` over finite `f64` pairs, in integer
/// arithmetic, so no product underflows or overflows.
fn dot_sign_wide(terms: &[(f64, f64)]) -> Sign {
    let parts: Vec<(bool, u128, i32)> = terms
        .iter()
        .filter(|(x, y)| *x != 0.0 && *y != 0.0)
        .map(|&(x, y)| {
            let ((sx, mx, ex), (sy, my, ey)) = (split_f64(x), split_f64(y));
            (sx != sy, mx as u128 * my as u128, ex + ey)
        })
        .collect();
    let Some(emin) = parts.iter().map(|t| t.2).min() else {
        return Sign::Zero;
    };
    let emax = parts.iter().map(|t| t.2).max().expect("nonempty");
    let limbs = (emax - emin) as usize / 64 + 4;
    let mut acc = [vec![0u64; limbs], vec![0u64; limbs]];
    for &(neg, m, e) in &parts {
        let shift = (e - emin) as usize;
        let (i, bit) = (shift / 64, shift % 64);
        let a = &mut acc[usize::from(neg)];
        add_at(a, i, (m & u64::MAX as u128) << bit);
        add_at(a, i + 1, (m >> 64) << bit);
    }
    match acc[0].iter().rev().cmp(acc[1].iter().rev()) {
        Ordering::Greater => Sign::Positive,
        Ordering::Equal => Sign::Zero,
        Ordering::Less => Sign::Negative,
    }
}

/// Exact sign of `(p - origin) . normal` for any finite input.
pub(crate) fn plane_side(p: [f64; 3], origin: [f64; 3], normal: [f64; 3]) -> Sign {
    let safe = |x: f64| x == 0.0 || (x.abs() >= 2f64.powi(-450) && x.abs() <= 2f64.powi(450));
    let terms = [0, 1, 2].map(|k| [(p[k], normal[k]), (-origin[k], normal[k])]);
    let terms = terms.as_flattened();
    if terms.iter().all(|&(x, y)| safe(x) && safe(y)) {
        terms
            .iter()
            .fold(Expansion::zero(), |acc, &(x, y)| {
                acc.sum(&Expansion::from_f64(x).scale(y))
            })
            .sign()
    } else {
        dot_sign_wide(terms)
    }
}

/// `normal` scaled by a power of two so that its largest component lies
/// near 1, when that scaling is exact; else `normal` unchanged. Side tests
/// and plane crossings do not change, and their expansions avoid underflow.
pub(crate) fn unit_scale(normal: [f64; 3]) -> [f64; 3] {
    let big = normal.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    if big == 0.0 {
        return normal;
    }
    let k = -split_f64(big).2 - 52;
    let scale = |x: f64, k: i32| {
        let (half, rest) = (k / 2, k - k / 2);
        x * 2f64.powi(half) * 2f64.powi(rest)
    };
    let n = normal.map(|x| scale(x, k));
    if (0..3).all(|i| n[i].is_finite() && scale(n[i], -k) == normal[i]) {
        n
    } else {
        normal
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use proptest::prelude::*;

    type Q = BigRational;

    fn q(x: f64) -> Q {
        Q::from_float(x).expect("finite")
    }

    fn qz() -> Q {
        Q::from_integer(BigInt::from(0))
    }

    fn val(e: &Expansion) -> Q {
        e.components().iter().fold(qz(), |a, &c| a + q(c))
    }

    fn sgn(x: &Q) -> Sign {
        match x.cmp(&qz()) {
            Ordering::Equal => Sign::Zero,
            Ordering::Greater => Sign::Positive,
            Ordering::Less => Sign::Negative,
        }
    }

    fn abs(x: Q) -> Q {
        if x < qz() { -x } else { x }
    }

    fn exact(p: &Pt) -> [Q; 3] {
        match p {
            Pt::E(v) => v.map(q),
            Pt::H(h) => {
                let w = val(&h.w);
                [0, 1, 2].map(|k| val(&h.x[k]) / &w)
            }
        }
    }

    fn nearest(x: &Q) -> f64 {
        let guess = {
            let n: f64 = x.numer().to_string().parse().unwrap_or(0.0);
            let d: f64 = x.denom().to_string().parse().unwrap_or(1.0);
            n / d
        };
        let mut best = guess;
        let mut c = guess;
        for _ in 0..8 {
            c = c.next_down();
        }
        for _ in 0..17 {
            let (db, dc) = (abs(q(best) - x), abs(q(c) - x));
            if dc < db || (dc == db && c.to_bits() & 1 == 0) {
                best = c;
            }
            c = c.next_up();
        }
        best
    }

    fn coord() -> impl Strategy<Value = f64> {
        prop_oneof![
            (-8i32..8).prop_map(|i| i as f64 * 0.5),
            -1e3..1e3f64,
            (-4i32..4, -3i32..3).prop_map(|(i, u)| i as f64 + u as f64 * f64::EPSILON),
        ]
    }

    fn point() -> impl Strategy<Value = [f64; 3]> {
        proptest::array::uniform3(coord())
    }

    fn any_pt() -> impl Strategy<Value = Pt> {
        prop_oneof![
            point().prop_map(Pt::E),
            (point(), point(), point(), point(), point()).prop_filter_map(
                "parallel",
                |(p, q, a, b, c)| {
                    let pl = Plane::new([a, b, c]);
                    let probe = Pt::lpi_checked(p, q, &pl)?;
                    Some(probe)
                }
            ),
        ]
    }

    impl Pt {
        fn lpi_checked(p: [f64; 3], q: [f64; 3], plane: &Plane) -> Option<Pt> {
            let d = (0..3).fold(Expansion::zero(), |acc, k| {
                acc.sum(&plane.n[k].product(&Expansion::diff_f64(q[k], p[k])))
            });
            (d.sign() != Sign::Zero).then(|| Pt::lpi(p, q, plane))
        }
    }

    #[test]
    fn expansion_ops_are_exact() {
        let a = Expansion::from_f64(1e30)
            .sum(&Expansion::from_f64(1.0))
            .sum(&Expansion::from_f64(1e-30));
        let b = Expansion::diff_f64(3.0, 1e-20);
        assert_eq!(
            val(&a.product(&b)),
            (q(1e30) + q(1.0) + q(1e-30)) * (q(3.0) - q(1e-20))
        );
        assert_eq!(a.diff(&a).sign(), Sign::Zero);
        assert!(
            a.product(&b)
                .components()
                .windows(2)
                .all(|w| w[0].abs() < w[1].abs())
        );
    }

    #[test]
    fn plane_side_survives_underflow_and_overflow() {
        let o = [0.0; 3];
        assert_eq!(
            plane_side([-1e-200, 0.0, 0.0], o, [1e-200, 0.0, 0.0]),
            Sign::Negative
        );
        assert_eq!(
            plane_side([0.0, 1e-200, 0.0], o, [1.0, 1e-300, 0.0]),
            Sign::Positive
        );
        assert_eq!(
            plane_side([1e300, -1e300, 0.0], o, [1e10, 1e10, 0.0]),
            Sign::Zero
        );
        assert_eq!(
            plane_side([1e300, -1e300, 5e-324], o, [1e10, 1e10, 1.0]),
            Sign::Positive
        );
        assert_eq!(
            unit_scale([1e-200, 0.0, -3e-201])[0],
            1e-200 * 2f64.powi(665)
        );
        assert_eq!(unit_scale([5e-324, 0.0, 0.0]), [2f64.powi(-52), 0.0, 0.0]);
    }

    #[test]
    fn lpi_of_axis_line_and_plane_is_exact() {
        let pl = Plane::new([
            [0.0, 0.0, 1.0 / 3.0],
            [1.0, 0.0, 1.0 / 3.0],
            [0.0, 1.0, 1.0 / 3.0],
        ]);
        let p = Pt::lpi([0.25, 0.5, -1.0], [0.25, 0.5, 2.0], &pl);
        assert_eq!(p.to_f64(), [0.25, 0.5, 1.0 / 3.0]);
        assert_eq!(side(&p, &pl), Sign::Zero);
    }

    proptest! {
        #[test]
        fn xp_matches_rationals(a in -1e30..1e30f64, b in -1e3..1e3f64, c in -1e-3..1e-3f64, d in -1e10..1e10f64) {
            let x = Expansion::from_f64(a).sum(&Expansion::from_f64(b)).sum(&Expansion::from_f64(c));
            let y = Expansion::diff_f64(d, b);
            prop_assert_eq!(val(&x.product(&y)), (q(a) + q(b) + q(c)) * (q(d) - q(b)));
            prop_assert_eq!(val(&x.scale(d)), (q(a) + q(b) + q(c)) * q(d));
            prop_assert_eq!(val(&x.diff(&y)), q(a) + q(b) + q(c) - (q(d) - q(b)));
        }

        #[test]
        fn predicates_match_rationals(a in any_pt(), b in any_pt(), c in any_pt(), drop in 0usize..3) {
            let (x, y, z) = (exact(&a), exact(&b), exact(&c));
            let (i, j) = kept_axes(drop);
            let det = (&y[i] - &x[i]) * (&z[j] - &x[j]) - (&y[j] - &x[j]) * (&z[i] - &x[i]);
            prop_assert_eq!(orient2d(&a, &b, &c, drop), sgn(&det));
            prop_assert_eq!(cmp_axis(&a, &b, drop), x[drop].cmp(&y[drop]));
            prop_assert_eq!(cmp_lex(&a, &a), Ordering::Equal);
        }

        #[test]
        fn side_and_on_plane_match_rationals(p in any_pt(), t in proptest::array::uniform3(point()), u in coord(), v in coord()) {
            let pl = Plane::new(t);
            let qv = t.map(|c| c.map(q));
            let e1: Vec<Q> = (0..3).map(|k| &qv[1][k] - &qv[0][k]).collect();
            let e2: Vec<Q> = (0..3).map(|k| &qv[2][k] - &qv[0][k]).collect();
            let n = [
                &e1[1] * &e2[2] - &e1[2] * &e2[1],
                &e1[2] * &e2[0] - &e1[0] * &e2[2],
                &e1[0] * &e2[1] - &e1[1] * &e2[0],
            ];
            let x = exact(&p);
            let s = (0..3).fold(qz(), |acc, k| acc + &n[k] * (&x[k] - &qv[0][k]));
            prop_assert_eq!(side(&p, &pl), sgn(&s));
            for drop in 0..3 {
                if pl.normal_sign(drop) != Sign::Zero {
                    let h = Pt::on_plane(&pl, drop, u, v);
                    prop_assert_eq!(side(&h, &pl), Sign::Zero);
                    let (i, j) = kept_axes(drop);
                    let hx = exact(&h);
                    prop_assert_eq!(&hx[i], &q(u));
                    prop_assert_eq!(&hx[j], &q(v));
                }
            }
        }

        #[test]
        fn to_f64_rounds_to_nearest(p in any_pt()) {
            let x = exact(&p);
            let r = p.to_f64();
            for k in 0..3 {
                prop_assert_eq!(r[k], nearest(&x[k]));
            }
        }

        #[test]
        fn plane_side_matches_rationals_at_any_scale(
            v in proptest::array::uniform9((-4i64..5, -1074i32..1020)),
        ) {
            let x: Vec<f64> = v.iter().map(|&(m, e)| m as f64 * 2f64.powi(e / 2) * 2f64.powi(e - e / 2)).collect();
            prop_assume!(x.iter().all(|c| c.is_finite()));
            let (p, o, n) = ([x[0], x[1], x[2]], [x[3], x[4], x[5]], [x[6], x[7], x[8]]);
            let s = (0..3).fold(qz(), |acc, k| acc + q(n[k]) * (q(p[k]) - q(o[k])));
            prop_assert_eq!(plane_side(p, o, n), sgn(&s));
            let u = unit_scale(n);
            prop_assert_eq!(plane_side(p, o, u), sgn(&s));
            for k in 0..3 {
                prop_assert_eq!(sgn(&q(u[k])), sgn(&q(n[k])));
            }
        }

        #[test]
        fn plane_cross_lies_on_plane(a in point(), b in point(), o in point(), n in proptest::array::uniform3(-4.0..4.0f64)) {
            let (sa, sb) = (plane_side(a, o, n), plane_side(b, o, n));
            prop_assume!(sa.times(sb) == Sign::Negative);
            let x = exact(&Pt::plane_cross(a, b, o, n));
            let s = (0..3).fold(qz(), |acc, k| acc + q(n[k]) * (&x[k] - q(o[k])));
            prop_assert_eq!(s, qz());
        }
    }
}
