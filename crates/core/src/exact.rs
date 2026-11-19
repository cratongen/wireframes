//! Exact arithmetic on floating-point expansions (Shewchuk 1997) and implicit
//! points (LPI, TPI) defined by explicit `f64` inputs.
//!
//! Every operation is exact for any finite `f64` inputs. An operation whose
//! floating-point components would underflow or overflow redoes itself on a
//! signed integer times a power of two, and its result keeps that form.

use std::cmp::Ordering;

use crate::predicates::{Sign, kept_axes};

/// A product at least this large keeps its rounding error representable, so
/// `two_product` is exact; smaller nonzero products may lose bits to underflow.
const TINY: f64 = f64::from_bits((1023 - 968) << 52);

/// `x = m 2^e` with an integer `m`, for finite `x`; the flag is the sign.
pub fn split_f64(x: f64) -> (bool, u64, i32) {
    let b = x.to_bits();
    let (exp, frac) = (((b >> 52) & 0x7ff) as i32, b & ((1 << 52) - 1));
    let (m, e) = if exp == 0 {
        (frac, -1074)
    } else {
        (frac | (1 << 52), exp - 1075)
    };
    (x < 0.0, m, e)
}

fn ldexp(mut x: f64, e: i64) -> f64 {
    let mut e = e.clamp(-2400, 2400) as i32;
    while e > 1000 {
        x *= 2f64.powi(1000);
        e -= 1000;
    }
    while e < -1000 {
        x *= 2f64.powi(-1000);
        e += 1000;
    }
    x * 2f64.powi(e)
}

/// `(-1)^neg mag 2^exp`, `mag` little-endian with no high zero limb and odd
/// when nonzero; zero is the empty `mag` with `neg` false.
#[derive(Clone, Debug, Default)]
struct Wide {
    neg: bool,
    mag: Vec<u64>,
    exp: i64,
}

fn cmp_mag(a: &[u64], b: &[u64]) -> Ordering {
    let len = |m: &[u64]| m.iter().rposition(|&l| l != 0).map_or(0, |i| i + 1);
    let (la, lb) = (len(a), len(b));
    la.cmp(&lb)
        .then_with(|| a[..la].iter().rev().cmp(b[..lb].iter().rev()))
}

fn add_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(a.len().max(b.len()) + 1);
    let mut carry = 0u128;
    for i in 0..a.len().max(b.len()) {
        let t = *a.get(i).unwrap_or(&0) as u128 + *b.get(i).unwrap_or(&0) as u128 + carry;
        out.push(t as u64);
        carry = t >> 64;
    }
    out.push(carry as u64);
    out
}

/// `a - b` for `a >= b`.
fn sub_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(a.len());
    let mut borrow = 0u64;
    for (i, &x) in a.iter().enumerate() {
        let (d, b1) = x.overflowing_sub(*b.get(i).unwrap_or(&0));
        let (d, b2) = d.overflowing_sub(borrow);
        out.push(d);
        borrow = u64::from(b1 || b2);
    }
    out
}

impl Wide {
    fn of(x: f64) -> Wide {
        let (neg, m, e) = split_f64(x);
        Wide::norm(neg, vec![m], e as i64)
    }

    fn norm(neg: bool, mut mag: Vec<u64>, mut exp: i64) -> Wide {
        while mag.last() == Some(&0) {
            mag.pop();
        }
        let Some(lo) = mag.iter().position(|&l| l != 0) else {
            return Wide::default();
        };
        mag.drain(..lo);
        exp += 64 * lo as i64;
        let tz = mag[0].trailing_zeros();
        if tz > 0 {
            for i in 0..mag.len() {
                let up = mag.get(i + 1).map_or(0, |&n| n << (64 - tz));
                mag[i] = (mag[i] >> tz) | up;
            }
            exp += tz as i64;
            if mag.last() == Some(&0) {
                mag.pop();
            }
        }
        Wide { neg, mag, exp }
    }

    fn shl(&self, s: i64) -> Vec<u64> {
        let s = usize::try_from(s).expect("exponent of the lower operand");
        let (limbs, bits) = (s / 64, s % 64);
        let mut out = vec![0u64; limbs + self.mag.len() + 1];
        for (i, &l) in self.mag.iter().enumerate() {
            out[limbs + i] |= l << bits;
            if bits > 0 {
                out[limbs + i + 1] |= l >> (64 - bits);
            }
        }
        out
    }

    fn add(&self, o: &Wide) -> Wide {
        if self.mag.is_empty() {
            return o.clone();
        }
        if o.mag.is_empty() {
            return self.clone();
        }
        let e = self.exp.min(o.exp);
        let (a, b) = (self.shl(self.exp - e), o.shl(o.exp - e));
        if self.neg == o.neg {
            return Wide::norm(self.neg, add_mag(&a, &b), e);
        }
        match cmp_mag(&a, &b) {
            Ordering::Equal => Wide::default(),
            Ordering::Greater => Wide::norm(self.neg, sub_mag(&a, &b), e),
            Ordering::Less => Wide::norm(o.neg, sub_mag(&b, &a), e),
        }
    }

    fn neg(&self) -> Wide {
        Wide {
            neg: !self.neg && !self.mag.is_empty(),
            ..self.clone()
        }
    }

    fn mul(&self, o: &Wide) -> Wide {
        if self.mag.is_empty() || o.mag.is_empty() {
            return Wide::default();
        }
        let mut out = vec![0u64; self.mag.len() + o.mag.len()];
        for (i, &a) in self.mag.iter().enumerate() {
            let mut carry = 0u128;
            for (j, &b) in o.mag.iter().enumerate() {
                let t = a as u128 * b as u128 + out[i + j] as u128 + carry;
                out[i + j] = t as u64;
                carry = t >> 64;
            }
            out[i + o.mag.len()] = carry as u64;
        }
        Wide::norm(self.neg != o.neg, out, self.exp + o.exp)
    }

    fn sign(&self) -> Sign {
        match (self.mag.is_empty(), self.neg) {
            (true, _) => Sign::Zero,
            (false, true) => Sign::Negative,
            (false, false) => Sign::Positive,
        }
    }

    /// `(m, e)` with the value within a few ulps of `m 2^e`, `|m|` in `[2^63, 2^64]`.
    fn approx(&self) -> (f64, i64) {
        let Some(&top) = self.mag.last() else {
            return (0.0, 0);
        };
        let n = self.mag.len();
        let lz = top.leading_zeros();
        let mut hi = top << lz;
        if lz > 0 && n > 1 {
            hi |= self.mag[n - 2] >> (64 - lz);
        }
        let m = if self.neg { -(hi as f64) } else { hi as f64 };
        (m, self.exp + 64 * (n as i64 - 1) - lz as i64)
    }
}

/// An exact real number: a sum of nonoverlapping `f64` components in
/// increasing magnitude with zeros removed, or, once those would underflow or
/// overflow, a signed integer times a power of two.
#[derive(Clone, Debug)]
pub struct Expansion(Repr);

#[derive(Clone, Debug)]
enum Repr {
    Fast(Vec<f64>),
    Wide(Wide),
}

impl Default for Expansion {
    fn default() -> Self {
        Expansion::zero()
    }
}

impl PartialEq for Expansion {
    fn eq(&self, o: &Self) -> bool {
        self.diff(o).sign() == Sign::Zero
    }
}

fn compress(mut e: Vec<f64>) -> Vec<f64> {
    let m = e.len();
    if m < 2 {
        return e;
    }
    let mut g = vec![0.0; m];
    let mut bottom = m - 1;
    let mut q = e[m - 1];
    for i in (0..m - 1).rev() {
        let (s, t) = fast_two_sum(q, e[i]);
        if t != 0.0 {
            g[bottom] = s;
            bottom -= 1;
            q = t;
        } else {
            q = s;
        }
    }
    g[bottom] = q;
    e.clear();
    let mut q = g[bottom];
    for &c in &g[bottom + 1..] {
        let (s, t) = fast_two_sum(c, q);
        if t != 0.0 {
            e.push(t);
        }
        q = s;
    }
    if q != 0.0 {
        e.push(q);
    }
    e
}

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

/// A NaN error term marks a product too small for an exact error term; the
/// caller then redoes the operation in [`Wide`].
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    if x.abs() < TINY && a != 0.0 && b != 0.0 {
        return (x, f64::NAN);
    }
    (x, a.mul_add(b, -x))
}

impl Expansion {
    fn checked(v: Vec<f64>, slow: impl FnOnce() -> Wide) -> Self {
        if v.iter().all(|c| c.is_finite()) {
            Expansion(Repr::Fast(v))
        } else {
            Expansion(Repr::Wide(slow()))
        }
    }

    fn wide(&self) -> Wide {
        match &self.0 {
            Repr::Fast(v) => v
                .iter()
                .fold(Wide::default(), |acc, &c| acc.add(&Wide::of(c))),
            Repr::Wide(w) => w.clone(),
        }
    }

    pub fn zero() -> Self {
        Expansion(Repr::Fast(vec![]))
    }

    /// `x` must be finite.
    pub fn from_f64(x: f64) -> Self {
        debug_assert!(x.is_finite());
        Expansion(Repr::Fast(if x == 0.0 { vec![] } else { vec![x] }))
    }

    /// `a - b` exactly.
    pub fn diff_f64(a: f64, b: f64) -> Self {
        let (x, y) = two_sum(a, -b);
        Self::checked([y, x].into_iter().filter(|&c| c != 0.0).collect(), || {
            Wide::of(a).add(&Wide::of(-b))
        })
    }

    /// The floating-point components, or `None` for a value held as an
    /// integer times a power of two.
    pub fn components(&self) -> Option<&[f64]> {
        match &self.0 {
            Repr::Fast(v) => Some(v),
            Repr::Wide(_) => None,
        }
    }

    /// The exact value as `(negative, magnitude, exponent)`: `magnitude` is
    /// little-endian `u64` limbs and the value is `±magnitude 2^exponent`.
    pub fn to_dyadic(&self) -> (bool, Vec<u64>, i64) {
        let w = self.wide();
        (w.neg, w.mag, w.exp)
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
        let (Repr::Fast(a), Repr::Fast(b)) = (&self.0, &o.0) else {
            return Expansion(Repr::Wide(self.wide().add(&o.wide())));
        };
        let (big, small) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        let mut e = big.clone();
        for &c in small {
            e = Self::grow(e, c);
        }
        Self::checked(compress(e), || self.wide().add(&o.wide()))
    }

    pub fn neg(&self) -> Self {
        Expansion(match &self.0 {
            Repr::Fast(v) => Repr::Fast(v.iter().map(|c| -c).collect()),
            Repr::Wide(w) => Repr::Wide(w.neg()),
        })
    }

    pub fn diff(&self, o: &Self) -> Self {
        self.sum(&o.neg())
    }

    /// `b` must be finite.
    pub fn scale(&self, b: f64) -> Self {
        let v = match &self.0 {
            Repr::Wide(w) => return Expansion(Repr::Wide(w.mul(&Wide::of(b)))),
            Repr::Fast(v) => v,
        };
        if b == 0.0 || v.is_empty() {
            return Expansion::zero();
        }
        let mut out = Vec::with_capacity(2 * v.len());
        let (mut q, h) = two_product(v[0], b);
        if h != 0.0 {
            out.push(h);
        }
        for &c in &v[1..] {
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
        Self::checked(compress(out), || self.wide().mul(&Wide::of(b)))
    }

    pub fn product(&self, o: &Self) -> Self {
        let (Repr::Fast(a), Repr::Fast(b)) = (&self.0, &o.0) else {
            return Expansion(Repr::Wide(self.wide().mul(&o.wide())));
        };
        let (big, small) = if a.len() >= b.len() {
            (self, b)
        } else {
            (o, a)
        };
        small
            .iter()
            .fold(Expansion::zero(), |acc, &c| acc.sum(&big.scale(c)))
    }

    pub fn sign(&self) -> Sign {
        match &self.0 {
            Repr::Fast(v) => v.last().map_or(Sign::Zero, |&c| Sign::of(c)),
            Repr::Wide(w) => w.sign(),
        }
    }

    /// Within a few ulps of the exact value, with its exact sign: a nonzero
    /// value below the smallest subnormal gives that subnormal, and one above
    /// `f64::MAX` gives `f64::MAX`.
    pub fn estimate(&self) -> f64 {
        match &self.0 {
            Repr::Fast(v) => v.iter().fold(0.0, |a, b| a + b),
            Repr::Wide(w) => {
                let (m, e) = w.approx();
                let x = ldexp(m, e);
                if m == 0.0 {
                    0.0
                } else if x == 0.0 {
                    f64::from_bits(1).copysign(m)
                } else {
                    x.clamp(-f64::MAX, f64::MAX)
                }
            }
        }
    }
}

/// `n / d` within a few ulps for `d` nonzero, at any magnitude.
pub fn quotient_estimate(n: &Expansion, d: &Expansion) -> f64 {
    if let (Repr::Fast(_), Repr::Fast(_)) = (&n.0, &d.0) {
        let (a, b) = (n.estimate(), d.estimate());
        let q = a / b;
        if (a == 0.0 || a.is_normal()) && b.is_normal() && (q == 0.0 || q.is_normal()) {
            return q;
        }
    }
    let ((mn, en), (md, ed)) = (n.wide().approx(), d.wide().approx());
    ldexp(mn / md, en - ed)
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

/// Sign of the volume enclosed by a closed, consistently wound triangle
/// shell: positive when its normals point outward.
///
/// The `f64` sum relative to one vertex decides when it exceeds its error
/// bound: Shewchuk's orient3d bound (`7u + 56u^2` times the permanent) per
/// face, which covers the rounded differences, plus the recursive summation
/// bound. Differences outside `[2^-300, 2^300]` could underflow or overflow
/// and skip the filter. Otherwise the sum is exact.
pub fn shell_volume_sign(v: &[[f64; 3]], faces: &[[u32; 3]]) -> Sign {
    let o = faces.first().map_or([0.0; 3], |t| v[t[0] as usize]);
    let (lo, hi) = (2f64.powi(-300), 2f64.powi(300));
    let mut filtered = faces.len() as f64 * f64::EPSILON < 0.01;
    let (mut s, mut perm, mut mag) = (0.0f64, 0.0f64, 0.0f64);
    for t in faces {
        if !filtered {
            break;
        }
        let [a, b, c] = t.map(|i| {
            let p = v[i as usize];
            [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
        });
        filtered = [a, b, c]
            .as_flattened()
            .iter()
            .all(|&x| x == 0.0 || (lo..=hi).contains(&x.abs()));
        let m = [
            [b[1] * c[2], b[2] * c[1]],
            [b[2] * c[0], b[0] * c[2]],
            [b[0] * c[1], b[1] * c[0]],
        ];
        let d =
            a[0] * (m[0][0] - m[0][1]) + a[1] * (m[1][0] - m[1][1]) + a[2] * (m[2][0] - m[2][1]);
        s += d;
        mag += d.abs();
        perm += (0..3)
            .map(|k| a[k].abs() * (m[k][0].abs() + m[k][1].abs()))
            .sum::<f64>();
    }
    if filtered {
        let (u, n) = (f64::EPSILON / 2.0, faces.len() as f64);
        let bound = (8.0 * u * perm + (n + 2.0) * u * mag) * (1.0 + 4.0 * n * u);
        if s.abs() > bound {
            return Sign::of(s);
        }
    }
    faces
        .iter()
        .fold(Expansion::zero(), |acc, t| {
            let [a, b, c] = t.map(|i| v[i as usize]);
            acc.sum(&dot(&ex(a), &cross(&ex(b), &ex(c))))
        })
        .sign()
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

/// Nearest `f64` to `n / d` (ties to even) for `d > 0`.
pub fn round_quotient(n: &Expansion, d: &Expansion) -> f64 {
    let times_d = |c: f64| {
        if c.is_finite() {
            d.scale(c)
        } else {
            Expansion::from_f64(2f64.powi(1023).copysign(c))
                .scale(2.0)
                .product(d)
        }
    };
    let above = |c: f64| n.diff(&times_d(c)).sign();
    let c = quotient_estimate(n, d) + 0.0;
    if c.is_nan() {
        return c;
    }
    let (lo, hi) = match above(c) {
        Sign::Zero => return c,
        Sign::Positive => {
            let mut lo = c;
            loop {
                if lo == f64::INFINITY {
                    return lo;
                }
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
                if hi == f64::NEG_INFINITY {
                    return hi;
                }
                let lo = hi.next_down();
                match above(lo) {
                    Sign::Zero => return lo,
                    Sign::Positive => break (lo, hi),
                    Sign::Negative => hi = lo,
                }
            }
        }
    };
    match n.scale(2.0).diff(&times_d(lo)).diff(&times_d(hi)).sign() {
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

type V2 = [Expansion; 2];

fn sub2(a: [f64; 2], b: [f64; 2]) -> V2 {
    [0, 1].map(|k| Expansion::diff_f64(a[k], b[k]))
}

fn cross2(u: &V2, v: &V2) -> Expansion {
    u[0].product(&v[1]).diff(&v[0].product(&u[1]))
}

/// [`crate::predicates::orient2d`] at any magnitude.
pub(crate) fn orient2d_any(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    cross2(&sub2(a, c), &sub2(b, c)).estimate()
}

/// [`crate::predicates::orient3d`] at any magnitude.
pub(crate) fn orient3d_any(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    dot(&sub_f(a, d), &cross(&sub_f(b, d), &sub_f(c, d))).estimate()
}

/// [`crate::predicates::incircle`] at any magnitude.
pub(crate) fn incircle_any(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> f64 {
    let r = [a, b, c].map(|p| sub2(p, d));
    let lift = |v: &V2| v[0].product(&v[0]).sum(&v[1].product(&v[1]));
    lift(&r[0])
        .product(&cross2(&r[1], &r[2]))
        .sum(&lift(&r[1]).product(&cross2(&r[2], &r[0])))
        .sum(&lift(&r[2]).product(&cross2(&r[0], &r[1])))
        .estimate()
}

/// [`crate::predicates::insphere`] at any magnitude.
pub(crate) fn insphere_any(p: [[f64; 3]; 4], e: [f64; 3]) -> f64 {
    let r = p.map(|x| sub_f(x, e));
    let xy = |i: usize, j: usize| r[i][0].product(&r[j][1]).diff(&r[j][0].product(&r[i][1]));
    let (ab, bc, cd, da, ac, bd) = (xy(0, 1), xy(1, 2), xy(2, 3), xy(3, 0), xy(0, 2), xy(1, 3));
    let z = |i: usize| &r[i][2];
    let abc = z(0)
        .product(&bc)
        .diff(&z(1).product(&ac))
        .sum(&z(2).product(&ab));
    let bcd = z(1)
        .product(&cd)
        .diff(&z(2).product(&bd))
        .sum(&z(3).product(&bc));
    let cda = z(2)
        .product(&da)
        .sum(&z(3).product(&ac))
        .sum(&z(0).product(&cd));
    let dab = z(3)
        .product(&ab)
        .sum(&z(0).product(&bd))
        .sum(&z(1).product(&da));
    let lift = |i: usize| dot(&r[i], &r[i]);
    lift(3)
        .product(&abc)
        .diff(&lift(2).product(&dab))
        .sum(&lift(1).product(&cda).diff(&lift(0).product(&bcd)))
        .estimate()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oracle::{self, Q, near_degenerate, q, qv, value};
    use proptest::prelude::*;

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

    fn ldexp2(x: f64, k: i32) -> f64 {
        x * 2f64.powi(k / 2) * 2f64.powi(k - k / 2)
    }

    fn small_point() -> impl Strategy<Value = [f64; 3]> {
        proptest::array::uniform3((-8i32..9).prop_map(f64::from))
    }

    #[test]
    fn quotients_round_to_nearest_at_any_magnitude() {
        let one = Expansion::from_f64(1.0);
        let t = Expansion::from_f64(3.0 * 2f64.powi(-600))
            .product(&Expansion::from_f64(5.0 * 2f64.powi(-600)));
        assert!(t.components().is_none());
        assert_eq!(
            value(&t),
            q(3.0 * 2f64.powi(-600)) * q(5.0 * 2f64.powi(-600))
        );
        assert_eq!(t.estimate(), f64::from_bits(1));
        assert_eq!(t.neg().estimate(), -f64::from_bits(1));
        assert_eq!(round_quotient(&t.scale(0.1), &t), 0.1);
        assert_eq!(round_quotient(&t.scale(-0.1), &t), -0.1);
        let tiny = f64::from_bits(3);
        assert_eq!(round_quotient(&t.scale(tiny), &t), tiny);
        let tie = one.sum(&Expansion::from_f64(2f64.powi(-53)));
        assert_eq!(round_quotient(&tie.product(&t), &t), 1.0);
        let max = Expansion::from_f64(f64::MAX);
        assert_eq!(round_quotient(&max.product(&t), &t), f64::MAX);
        assert_eq!(round_quotient(&max.scale(4.0), &one), f64::INFINITY);
        assert_eq!(round_quotient(&max.scale(-4.0), &one), f64::NEG_INFINITY);
        let ulp = 2f64.powi(970);
        let below = max.sum(&Expansion::from_f64(ulp / 2.0));
        assert_eq!(round_quotient(&below, &one), f64::MAX);
        let mid = max.sum(&Expansion::from_f64(ulp));
        assert_eq!(round_quotient(&mid, &one), f64::INFINITY);
        assert_eq!(max.scale(4.0).estimate(), f64::MAX);
    }

    #[test]
    fn shell_volume_sign_survives_underflow_and_overflow() {
        for k in [-1064, 0, 498] {
            let m = crate::fixtures::cube();
            let v: Vec<[f64; 3]> = m
                .vertices()
                .iter()
                .map(|p| p.map(|x| ldexp2(x, k)))
                .collect();
            assert_eq!(
                shell_volume_sign(&v, m.triangles()),
                Sign::Positive,
                "2^{k}"
            );
            let flipped: Vec<[u32; 3]> = m.triangles().iter().map(|t| [t[0], t[2], t[1]]).collect();
            assert_eq!(shell_volume_sign(&v, &flipped), Sign::Negative, "2^{k}");
        }
    }

    proptest! {
        #[test]
        fn arithmetic_is_exact_at_any_scale(
            v in proptest::array::uniform4((-(1i64 << 53)..(1i64 << 53), -1074i32..971)),
        ) {
            let x: Vec<f64> = v.iter().map(|&(m, e)| ldexp2(m as f64, e)).collect();
            prop_assume!(x.iter().all(|c| c.is_finite()));
            let (a, b, c, d) = (x[0], x[1], x[2], x[3]);
            let s = Expansion::diff_f64(a, b);
            let p = s.product(&Expansion::diff_f64(c, d)).scale(a);
            let want = (q(a) - q(b)) * (q(c) - q(d)) * q(a);
            prop_assert_eq!(value(&p), want.clone());
            prop_assert_eq!(p.sign(), oracle::sign(&want));
            prop_assert_eq!(value(&s.sum(&p)), q(a) - q(b) + &want);
            prop_assert_eq!(value(&p.diff(&p)), q(0.0));
            prop_assert_eq!(Sign::of(p.estimate()), oracle::sign(&want));
        }

        #[test]
        fn predicates_are_exact_at_any_scale(p in proptest::array::uniform5(small_point()), k in -1070i32..1000) {
            let s = |x: [f64; 3]| x.map(|c| ldexp2(c, k));
            let t = |x: [f64; 3]| [x[0], x[1]];
            let [a, b, c, d, e] = p;
            let (sa, sb, sc, sd, se) = (s(a), s(b), s(c), s(d), s(e));
            use crate::predicates::{incircle, insphere, orient2d, orient3d};
            prop_assert_eq!(Sign::of(orient2d(t(sa), t(sb), t(sc))), Sign::of(orient2d(t(a), t(b), t(c))));
            prop_assert_eq!(Sign::of(orient3d(sa, sb, sc, sd)), Sign::of(orient3d(a, b, c, d)));
            prop_assert_eq!(
                Sign::of(incircle(t(sa), t(sb), t(sc), t(sd))),
                Sign::of(incircle(t(a), t(b), t(c), t(d)))
            );
            prop_assert_eq!(Sign::of(insphere(sa, sb, sc, sd, se)), Sign::of(insphere(a, b, c, d, e)));
            let lpi = ImplicitPoint::Lpi { p: a, q: b, r: c, s: d, t: e };
            let slpi = ImplicitPoint::Lpi { p: sa, q: sb, r: sc, s: sd, t: se };
            prop_assume!(lpi.is_valid());
            prop_assert!(slpi.is_valid());
            prop_assert_eq!(orient3d_i(&slpi, sc, sd, se), orient3d_i(&lpi, c, d, e));
            prop_assert_eq!(orient3d_i(&slpi, sa, sd, se), orient3d_i(&lpi, a, d, e));
            prop_assert_eq!(cmp_axis(&slpi, &ImplicitPoint::Explicit(sb), 2), cmp_axis(&lpi, &ImplicitPoint::Explicit(b), 2));
        }

        #[test]
        fn shell_volume_sign_matches_oracle(p in proptest::array::uniform4(near_degenerate()), lift in 0usize..3) {
            let mut p = p;
            p[3][2] = p[0][2] + [0.0, 2f64.powi(-80), -2f64.powi(-80)][lift];
            let faces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
            let want = faces.iter().fold(q(0.0), |acc, t| {
                let [a, b, c] = t.map(|i: u32| qv(p[i as usize]));
                acc + oracle::dot(&a, &oracle::cross(&b, &c))
            });
            prop_assert_eq!(shell_volume_sign(&p, &faces), oracle::sign(&want));
        }

        #[test]
        fn fallback_predicates_match_robust(a in near_degenerate(), b in near_degenerate(), c in near_degenerate(), d in near_degenerate(), e in near_degenerate()) {
            let t = |x: [f64; 3]| [x[0], x[1]];
            let r2 = |x: [f64; 3]| robust::Coord { x: x[0], y: x[1] };
            let r3 = |x: [f64; 3]| robust::Coord3D { x: x[0], y: x[1], z: x[2] };
            prop_assert_eq!(Sign::of(orient2d_any(t(a), t(b), t(c))), Sign::of(robust::orient2d(r2(a), r2(b), r2(c))));
            prop_assert_eq!(Sign::of(orient3d_any(a, b, c, d)), Sign::of(robust::orient3d(r3(a), r3(b), r3(c), r3(d))));
            prop_assert_eq!(
                Sign::of(incircle_any(t(a), t(b), t(c), t(d))),
                Sign::of(robust::incircle(r2(a), r2(b), r2(c), r2(d)))
            );
            prop_assert_eq!(
                Sign::of(insphere_any([a, b, c, d], e)),
                Sign::of(robust::insphere(r3(a), r3(b), r3(c), r3(d), r3(e)))
            );
        }

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
