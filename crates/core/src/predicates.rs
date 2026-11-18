//! Exact orientation predicates. Each returns an `f64` whose sign is exact;
//! its magnitude only approximates the determinant.

use robust::{Coord, Coord3D};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Sign {
    Negative,
    Zero,
    Positive,
}

impl Sign {
    pub fn of(x: f64) -> Sign {
        if x > 0.0 {
            Sign::Positive
        } else if x < 0.0 {
            Sign::Negative
        } else {
            Sign::Zero
        }
    }

    pub fn to_i8(self) -> i8 {
        match self {
            Sign::Negative => -1,
            Sign::Zero => 0,
            Sign::Positive => 1,
        }
    }

    pub fn flip(self) -> Sign {
        match self {
            Sign::Negative => Sign::Positive,
            Sign::Zero => Sign::Zero,
            Sign::Positive => Sign::Negative,
        }
    }

    pub fn times(self, o: Sign) -> Sign {
        Sign::of((self.to_i8() * o.to_i8()) as f64)
    }
}

fn c2(p: [f64; 2]) -> Coord<f64> {
    Coord { x: p[0], y: p[1] }
}

fn c3(p: [f64; 3]) -> Coord3D<f64> {
    Coord3D {
        x: p[0],
        y: p[1],
        z: p[2],
    }
}

/// Positive when `c` is left of `a -> b` (counter-clockwise).
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    robust::orient2d(c2(a), c2(b), c2(c))
}

/// `robust` convention: positive when `d` lies below the plane through `a, b, c`
/// that winds counter-clockwise seen from above, i.e. on the side opposite to the
/// normal `(b - a) x (c - a)`. Prefer [`side_of_plane`].
pub fn orient3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    robust::orient3d(c3(a), c3(b), c3(c), c3(d))
}

/// Positive when `d` is inside the circle through counter-clockwise `a, b, c`.
pub fn incircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> f64 {
    robust::incircle(c2(a), c2(b), c2(c), c2(d))
}

/// Positive when `e` is inside the sphere through `a, b, c, d` (`orient3d(a, b, c, d) > 0`).
pub fn insphere(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3], e: [f64; 3]) -> f64 {
    robust::insphere(c3(a), c3(b), c3(c), c3(d), c3(e))
}

/// The two axes kept when dropping `drop`, in cyclic order, so that the 2D
/// orientation of a projected triangle has the sign of its normal's `drop` component.
pub fn kept_axes(drop: usize) -> (usize, usize) {
    ((drop + 1) % 3, (drop + 2) % 3)
}

/// `orient2d` of the projection dropping axis `drop`; its sign is the sign of
/// component `drop` of `(b - a) x (c - a)`.
pub fn orient2d_axis(a: [f64; 3], b: [f64; 3], c: [f64; 3], drop: usize) -> f64 {
    let (i, j) = kept_axes(drop);
    orient2d([a[i], a[j]], [b[i], b[j]], [c[i], c[j]])
}

/// Which side of the plane through `a, b, c` holds `p`; `Positive` is the side
/// the normal `(b - a) x (c - a)` points to.
pub fn side_of_plane(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Sign {
    Sign::of(orient3d(a, b, c, p)).flip()
}

/// Whether the three points lie on one line (or coincide).
pub fn collinear3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> bool {
    (0..3).all(|k| orient2d_axis(a, b, c, k) == 0.0)
}

/// An axis whose projection keeps the non-degenerate triangle `t` non-degenerate,
/// preferring the largest normal component; `None` when `t` is degenerate.
pub fn projection_axis(t: [[f64; 3]; 3]) -> Option<usize> {
    let o = [0, 1, 2].map(|k| orient2d_axis(t[0], t[1], t[2], k));
    let best = (0..3).max_by(|&i, &j| o[i].abs().total_cmp(&o[j].abs()).then(j.cmp(&i)))?;
    (o[best] != 0.0).then_some(best)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oracle::{self, qv};
    use proptest::prelude::*;

    fn side_oracle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Sign {
        oracle::side(&qv(p), &qv(a), &qv(b), &qv(c))
    }

    fn orient2d_oracle(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign {
        oracle::orient2d(a, b, c)
    }

    use crate::oracle::near_degenerate;

    #[test]
    fn side_of_plane_follows_the_normal() {
        let (a, b, c) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(side_of_plane([0.0, 0.0, 1.0], a, b, c), Sign::Positive);
        assert_eq!(side_of_plane([0.0, 0.0, -1.0], a, b, c), Sign::Negative);
        assert_eq!(side_of_plane([5.0, 5.0, 0.0], a, b, c), Sign::Zero);
        assert!(orient2d([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]) > 0.0);
    }

    #[test]
    fn projection_signs_match_normal_components() {
        let t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.5], [0.0, 1.0, 2.0]];
        let n = [-0.5, -2.0, 1.0];
        for k in 0..3 {
            assert_eq!(Sign::of(orient2d_axis(t[0], t[1], t[2], k)), Sign::of(n[k]));
        }
        assert_eq!(projection_axis(t), Some(1));
        assert!(collinear3([0.0; 3], [1.0, 2.0, 3.0], [2.0, 4.0, 6.0]));
        assert_eq!(
            projection_axis([[0.0; 3], [1.0, 2.0, 3.0], [2.0, 4.0, 6.0]]),
            None
        );
    }

    proptest! {
        #[test]
        fn side_of_plane_agrees_with_rational_oracle(
            p in near_degenerate(), a in near_degenerate(), b in near_degenerate(), c in near_degenerate()
        ) {
            prop_assert_eq!(side_of_plane(p, a, b, c), side_oracle(p, a, b, c));
        }

        #[test]
        fn orient2d_agrees_with_rational_oracle(a in near_degenerate(), b in near_degenerate(), c in near_degenerate()) {
            let p = |x: [f64; 3]| [x[0], x[1]];
            prop_assert_eq!(Sign::of(orient2d(p(a), p(b), p(c))), orient2d_oracle(p(a), p(b), p(c)));
        }

        #[test]
        fn random_side_of_plane_agrees(
            p in proptest::array::uniform3(-1e6..1e6f64), a in proptest::array::uniform3(-1e6..1e6f64),
            b in proptest::array::uniform3(-1e6..1e6f64), c in proptest::array::uniform3(-1e6..1e6f64)
        ) {
            prop_assert_eq!(side_of_plane(p, a, b, c), side_oracle(p, a, b, c));
        }
    }
}
