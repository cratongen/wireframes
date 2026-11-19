//! Exact rational reference implementations for tests.

use num_rational::BigRational;

use crate::predicates::Sign;

pub type Q = BigRational;

pub fn q(x: f64) -> Q {
    Q::from_float(x).expect("finite")
}

pub fn sign(x: &Q) -> Sign {
    let zero = Q::from_integer(0.into());
    if *x > zero {
        Sign::Positive
    } else if *x < zero {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

pub fn qv(p: [f64; 3]) -> [Q; 3] {
    p.map(q)
}

pub fn sub(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

pub fn cross(u: &[Q; 3], v: &[Q; 3]) -> [Q; 3] {
    [
        &u[1] * &v[2] - &u[2] * &v[1],
        &u[2] * &v[0] - &u[0] * &v[2],
        &u[0] * &v[1] - &u[1] * &v[0],
    ]
}

pub fn dot(u: &[Q; 3], v: &[Q; 3]) -> Q {
    &u[0] * &v[0] + &u[1] * &v[1] + &u[2] * &v[2]
}

/// Side of `p` against the plane through `a, b, c`, positive along `(b - a) x (c - a)`.
pub fn side(p: &[Q; 3], a: &[Q; 3], b: &[Q; 3], c: &[Q; 3]) -> Sign {
    sign(&dot(&cross(&sub(b, a), &sub(c, a)), &sub(p, a)))
}

pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign {
    sign(&((q(b[0]) - q(a[0])) * (q(c[1]) - q(a[1])) - (q(b[1]) - q(a[1])) * (q(c[0]) - q(a[0]))))
}

/// Line through `p, q` meets the plane through `r, s, t`.
pub fn lpi(p: [f64; 3], qq: [f64; 3], r: [f64; 3], s: [f64; 3], t: [f64; 3]) -> Option<[Q; 3]> {
    let (p, qq, r, s, t) = (qv(p), qv(qq), qv(r), qv(s), qv(t));
    let n = cross(&sub(&s, &r), &sub(&t, &r));
    let d = sub(&qq, &p);
    let w = dot(&n, &d);
    if sign(&w) == Sign::Zero {
        return None;
    }
    let lam = dot(&n, &sub(&r, &p)) / w;
    Some([
        &p[0] + &d[0] * &lam,
        &p[1] + &d[1] * &lam,
        &p[2] + &d[2] * &lam,
    ])
}

/// Three planes, each through three points, meet.
pub fn tpi(planes: [[[f64; 3]; 3]; 3]) -> Option<[Q; 3]> {
    let n: Vec<[Q; 3]> = planes
        .iter()
        .map(|p| cross(&sub(&qv(p[1]), &qv(p[0])), &sub(&qv(p[2]), &qv(p[0]))))
        .collect();
    let d: Vec<Q> = (0..3).map(|i| dot(&n[i], &qv(planes[i][0]))).collect();
    let c12 = cross(&n[1], &n[2]);
    let w = dot(&n[0], &c12);
    if sign(&w) == Sign::Zero {
        return None;
    }
    let c20 = cross(&n[2], &n[0]);
    let c01 = cross(&n[0], &n[1]);
    Some([0, 1, 2].map(|k| (&d[0] * &c12[k] + &d[1] * &c20[k] + &d[2] * &c01[k]) / &w))
}

/// Points on a coarse grid near 1e3, each coordinate optionally moved by one ulp.
pub fn near_degenerate() -> impl proptest::strategy::Strategy<Value = [f64; 3]> {
    use proptest::strategy::Strategy;
    (
        proptest::array::uniform3(-4i32..4),
        proptest::array::uniform3(-1i8..=1),
    )
        .prop_map(|(g, u)| {
            [0, 1, 2].map(|k| {
                let x = g[k] as f64 * 0.5 + 1e3;
                match u[k] {
                    -1 => x.next_down(),
                    1 => x.next_up(),
                    _ => x,
                }
            })
        })
}

/// The nearest f64 to `x`, ties to even, by bracketing with rational comparisons.
pub fn nearest(x: &Q) -> f64 {
    let digits = |b: &num_bigint::BigInt| {
        let s = b.to_string();
        let (neg, s) = s
            .strip_prefix('-')
            .map_or((false, s.as_str()), |r| (true, r));
        let keep = s.len().min(17);
        let m: f64 = s[..keep].parse().unwrap();
        (if neg { -m } else { m }, (s.len() - keep) as i32)
    };
    let ((mn, en), (md, ed)) = (digits(x.numer()), digits(x.denom()));
    let mut c = mn / md * 10f64.powi(en - ed);
    if !c.is_finite() {
        c = 0.0;
    }
    while q(c) < *x {
        c = c.next_up();
    }
    while q(c) > *x {
        c = c.next_down();
    }
    let (lo, hi) = (c, c.next_up());
    if q(lo) == *x {
        return lo;
    }
    let mid = (q(lo) + q(hi)) / Q::from_integer(2.into());
    if *x < mid {
        lo
    } else if *x > mid {
        hi
    } else if lo.to_bits() & 1 == 0 {
        lo
    } else {
        hi
    }
}

/// The exact value of an expansion.
pub fn value(e: &crate::exact::Expansion) -> Q {
    let (neg, mag, exp) = e.to_dyadic();
    let n = mag
        .iter()
        .rev()
        .fold(num_bigint::BigInt::from(0), |acc, &l| {
            (acc << 64) + num_bigint::BigInt::from(l)
        });
    let v = Q::from_integer(if neg { -n } else { n });
    v * Q::from_integer(2.into()).pow(i32::try_from(exp).expect("exponent fits i32"))
}
