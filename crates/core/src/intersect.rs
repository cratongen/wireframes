//! Exact triangle-triangle and segment-triangle intersection tests. Triangles
//! must be non-degenerate (corners not collinear).

use std::cmp::Ordering;

use crate::exact::{ImplicitPoint, cmp_lex};
use crate::predicates::{Sign, orient2d_axis, projection_axis, side_of_plane};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TriTri {
    Disjoint,
    /// The closed triangles meet only where both boundaries meet.
    Touch,
    /// Some common point lies in the relative interior of one of them.
    Cross,
    /// Coplanar with an overlap of positive area.
    CoplanarOverlap,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegTri {
    Miss,
    /// The segment meets the triangle's plane once, inside the open triangle.
    Interior,
    /// ... on the open edge from corner `e` to corner `(e + 1) % 3`.
    Edge(u8),
    /// ... at corner `v`.
    Vertex(u8),
    /// The segment lies in the triangle's plane and meets the closed triangle.
    Coplanar,
}

fn sgn(x: f64) -> Sign {
    Sign::of(x)
}

/// Position of `p`, coplanar with `t`, in the closed triangle `t` projected
/// by dropping `axis`: `None` outside, else the `SegTri` location.
fn locate_2d(p: [f64; 3], t: [[f64; 3]; 3], axis: usize) -> Option<SegTri> {
    let s = sgn(orient2d_axis(t[0], t[1], t[2], axis));
    let o = [0, 1, 2].map(|i| sgn(orient2d_axis(t[i], t[(i + 1) % 3], p, axis)));
    if o.iter().any(|&x| x == s.flip()) {
        return None;
    }
    let zeros: Vec<u8> = (0..3u8).filter(|&i| o[i as usize] == Sign::Zero).collect();
    Some(match zeros[..] {
        [] => SegTri::Interior,
        [e] => SegTri::Edge(e),
        [e, f] => SegTri::Vertex(if f == e + 1 { f } else { 0 }),
        _ => unreachable!("t is non-degenerate"),
    })
}

/// Closed point-in-triangle; `p` must be coplanar with `t`, else false.
pub fn point_in_triangle(p: [f64; 3], t: [[f64; 3]; 3]) -> bool {
    if side_of_plane(p, t[0], t[1], t[2]) != Sign::Zero {
        return false;
    }
    projection_axis(t).is_some_and(|axis| locate_2d(p, t, axis).is_some())
}

fn on_segment_2d(p: [f64; 3], a: [f64; 3], b: [f64; 3], axis: usize) -> bool {
    let (i, j) = crate::predicates::kept_axes(axis);
    orient2d_axis(a, b, p, axis) == 0.0
        && p[i] >= a[i].min(b[i])
        && p[i] <= a[i].max(b[i])
        && p[j] >= a[j].min(b[j])
        && p[j] <= a[j].max(b[j])
}

/// Closed segments `ab` and `cd` meet in the projection dropping `axis`.
pub(crate) fn segments_meet_2d(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    axis: usize,
) -> bool {
    let o1 = sgn(orient2d_axis(a, b, c, axis));
    let o2 = sgn(orient2d_axis(a, b, d, axis));
    let o3 = sgn(orient2d_axis(c, d, a, axis));
    let o4 = sgn(orient2d_axis(c, d, b, axis));
    let proper = o1 == o2.flip() && o1 != Sign::Zero && o3 == o4.flip() && o3 != Sign::Zero;
    proper
        || on_segment_2d(c, a, b, axis)
        || on_segment_2d(d, a, b, axis)
        || on_segment_2d(a, c, d, axis)
        || on_segment_2d(b, c, d, axis)
}

/// Where the closed segment `pq` meets the closed triangle `t`, and the point
/// when it is a single one (explicit when it is an endpoint or a corner).
pub fn seg_tri(p: [f64; 3], q: [f64; 3], t: [[f64; 3]; 3]) -> (SegTri, Option<ImplicitPoint>) {
    let sp = side_of_plane(p, t[0], t[1], t[2]);
    let sq = side_of_plane(q, t[0], t[1], t[2]);
    if sp == sq && sp != Sign::Zero {
        return (SegTri::Miss, None);
    }
    let axis = projection_axis(t).expect("triangle must be non-degenerate");
    if sp == Sign::Zero && sq == Sign::Zero {
        let hit = locate_2d(p, t, axis).is_some()
            || locate_2d(q, t, axis).is_some()
            || (0..3).any(|e| segments_meet_2d(p, q, t[e], t[(e + 1) % 3], axis));
        return if hit {
            (SegTri::Coplanar, None)
        } else {
            (SegTri::Miss, None)
        };
    }
    if sp == Sign::Zero || sq == Sign::Zero {
        let x = if sp == Sign::Zero { p } else { q };
        return match locate_2d(x, t, axis) {
            None => (SegTri::Miss, None),
            Some(SegTri::Vertex(v)) => (
                SegTri::Vertex(v),
                Some(ImplicitPoint::Explicit(t[v as usize])),
            ),
            Some(k) => (k, Some(ImplicitPoint::Explicit(x))),
        };
    }
    let o = side_of_plane(t[2], p, t[0], t[1]);
    let s = [0, 1, 2].map(|i| side_of_plane(q, p, t[i], t[(i + 1) % 3]));
    if s.iter().any(|&x| x == o.flip()) {
        return (SegTri::Miss, None);
    }
    let zeros: Vec<u8> = (0..3u8).filter(|&i| s[i as usize] == Sign::Zero).collect();
    match zeros[..] {
        [] => (
            SegTri::Interior,
            Some(ImplicitPoint::Lpi {
                p,
                q,
                r: t[0],
                s: t[1],
                t: t[2],
            }),
        ),
        [e] => (
            SegTri::Edge(e),
            Some(ImplicitPoint::Lpi {
                p,
                q,
                r: t[0],
                s: t[1],
                t: t[2],
            }),
        ),
        [e, f] => {
            let v = if f == e + 1 { f } else { 0 };
            (
                SegTri::Vertex(v),
                Some(ImplicitPoint::Explicit(t[v as usize])),
            )
        }
        _ => unreachable!("p is off the plane of t"),
    }
}

/// Points of triangle `a` on the plane of `b`: corners on it and edge crossings.
/// `sa` are the sides of `a`'s corners against `b`'s plane.
fn plane_section(a: [[f64; 3]; 3], b: [[f64; 3]; 3], sa: [Sign; 3]) -> Vec<ImplicitPoint> {
    let mut pts = vec![];
    for i in 0..3 {
        let j = (i + 1) % 3;
        if sa[i] == Sign::Zero {
            pts.push(ImplicitPoint::Explicit(a[i]));
        }
        if sa[i] != Sign::Zero && sa[j] == sa[i].flip() {
            pts.push(ImplicitPoint::Lpi {
                p: a[i],
                q: a[j],
                r: b[0],
                s: b[1],
                t: b[2],
            });
        }
    }
    pts
}

fn extent(pts: Vec<ImplicitPoint>) -> (ImplicitPoint, ImplicitPoint) {
    let lo = pts
        .iter()
        .min_by(|x, y| cmp_lex(x, y))
        .expect("non-empty")
        .clone();
    let hi = pts
        .iter()
        .max_by(|x, y| cmp_lex(x, y))
        .expect("non-empty")
        .clone();
    (lo, hi)
}

/// Classifies the intersection of two closed, non-degenerate triangles exactly.
pub fn tri_tri(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> TriTri {
    let sb = b.map(|p| side_of_plane(p, a[0], a[1], a[2]));
    if sb[0] != Sign::Zero && sb.iter().all(|&s| s == sb[0]) {
        return TriTri::Disjoint;
    }
    let sa = a.map(|p| side_of_plane(p, b[0], b[1], b[2]));
    if sa[0] != Sign::Zero && sa.iter().all(|&s| s == sa[0]) {
        return TriTri::Disjoint;
    }
    if sb.iter().all(|&s| s == Sign::Zero) {
        return coplanar(a, b);
    }
    let a_straddles = sa.contains(&Sign::Positive) && sa.contains(&Sign::Negative);
    let b_straddles = sb.contains(&Sign::Positive) && sb.contains(&Sign::Negative);
    let (a0, a1) = extent(plane_section(a, b, sa));
    let (b0, b1) = extent(plane_section(b, a, sb));
    let lo = if cmp_lex(&a0, &b0) == Ordering::Less {
        &b0
    } else {
        &a0
    };
    let hi = if cmp_lex(&a1, &b1) == Ordering::Less {
        &a1
    } else {
        &b1
    };
    match cmp_lex(lo, hi) {
        Ordering::Greater => TriTri::Disjoint,
        Ordering::Less if a_straddles || b_straddles => TriTri::Cross,
        Ordering::Less => TriTri::Touch,
        Ordering::Equal => {
            let inside = |p: &ImplicitPoint, x0: &ImplicitPoint, x1: &ImplicitPoint| {
                cmp_lex(x0, p) == Ordering::Less && cmp_lex(p, x1) == Ordering::Less
            };
            if (a_straddles && inside(lo, &a0, &a1)) || (b_straddles && inside(lo, &b0, &b1)) {
                TriTri::Cross
            } else {
                TriTri::Touch
            }
        }
    }
}

fn coplanar(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> TriTri {
    let axis = projection_axis(a).expect("triangle must be non-degenerate");
    let mut weak = false;
    for (t, u) in [(a, b), (b, a)] {
        let s = sgn(orient2d_axis(t[0], t[1], t[2], axis));
        for i in 0..3 {
            let o = u.map(|p| sgn(orient2d_axis(t[i], t[(i + 1) % 3], p, axis)));
            if o.iter().all(|&x| x == s.flip()) {
                return TriTri::Disjoint;
            }
            weak |= o.iter().all(|&x| x != s);
        }
    }
    if weak {
        TriTri::Touch
    } else {
        TriTri::CoplanarOverlap
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const T: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];

    fn shift(t: [[f64; 3]; 3], d: [f64; 3]) -> [[f64; 3]; 3] {
        t.map(|p| [p[0] + d[0], p[1] + d[1], p[2] + d[2]])
    }

    #[test]
    fn crossing_touching_and_disjoint() {
        let vertical = [[0.2, 0.2, -1.0], [0.2, 0.2, 1.0], [0.6, 0.2, 0.5]];
        assert_eq!(tri_tri(T, vertical), TriTri::Cross);
        assert_eq!(tri_tri(T, shift(T, [0.0, 0.0, 1.0])), TriTri::Disjoint);
        let corner_on_face = [[0.2, 0.2, 0.0], [0.2, 0.2, 1.0], [0.6, 0.2, 1.0]];
        assert_eq!(tri_tri(T, corner_on_face), TriTri::Cross);
        let shared_edge = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!(tri_tri(T, shared_edge), TriTri::Touch);
        let shared_vertex = [[1.0, 0.0, 0.0], [2.0, 0.0, 1.0], [2.0, 1.0, 1.0]];
        assert_eq!(tri_tri(T, shared_vertex), TriTri::Touch);
        let edge_on_edge = [[0.5, 0.5, 0.0], [0.5, 0.5, 1.0], [1.0, 1.0, 1.0]];
        assert_eq!(tri_tri(T, edge_on_edge), TriTri::Touch);
    }

    #[test]
    fn interiors_crossing_between_boundary_points_is_cross() {
        let a = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let b = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.0]];
        assert_eq!(tri_tri(a, b), TriTri::Touch);
        let b = [[1.0, -1.0, -1.0], [1.0, 3.0, 1.0], [1.0, 3.0, -1.0]];
        assert_eq!(tri_tri(a, b), TriTri::Cross);
    }

    #[test]
    fn coplanar_cases() {
        assert_eq!(
            tri_tri(T, shift(T, [0.25, 0.25, 0.0])),
            TriTri::CoplanarOverlap
        );
        assert_eq!(tri_tri(T, shift(T, [1.0, 0.0, 0.0])), TriTri::Touch);
        assert_eq!(tri_tri(T, shift(T, [2.0, 0.0, 0.0])), TriTri::Disjoint);
        let mirrored = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        assert_eq!(tri_tri(T, mirrored), TriTri::Touch);
        assert_eq!(tri_tri(T, T), TriTri::CoplanarOverlap);
    }

    #[test]
    fn seg_tri_locations() {
        let up = |x: f64, y: f64| ([x, y, -1.0], [x, y, 1.0]);
        let (p, q) = up(0.25, 0.25);
        let (k, pt) = seg_tri(p, q, T);
        assert_eq!(k, SegTri::Interior);
        assert_eq!(pt.unwrap().to_f64(), [0.25, 0.25, 0.0]);
        let (p, q) = up(0.5, 0.0);
        assert_eq!(seg_tri(p, q, T).0, SegTri::Edge(0));
        let (p, q) = up(0.0, 1.0);
        assert_eq!(seg_tri(p, q, T).0, SegTri::Vertex(2));
        let (p, q) = up(0.0, 0.0);
        assert_eq!(seg_tri(p, q, T).0, SegTri::Vertex(0));
        let (p, q) = up(1.0, 1.0);
        assert_eq!(seg_tri(p, q, T).0, SegTri::Miss);
        assert_eq!(
            seg_tri([0.2, 0.2, 0.0], [0.2, 0.2, 1.0], T).0,
            SegTri::Interior
        );
        assert_eq!(seg_tri([0.2, 0.2, 0.5], [0.2, 0.2, 1.0], T).0, SegTri::Miss);
        assert_eq!(
            seg_tri([-1.0, 0.5, 0.0], [2.0, 0.5, 0.0], T).0,
            SegTri::Coplanar
        );
        assert_eq!(
            seg_tri([-1.0, 2.0, 0.0], [2.0, 2.0, 0.0], T).0,
            SegTri::Miss
        );
        assert!(point_in_triangle([0.5, 0.5, 0.0], T));
        assert!(!point_in_triangle([0.5, 0.5, 1e-300], T));
    }

    fn tri() -> impl Strategy<Value = [[f64; 3]; 3]> {
        proptest::array::uniform3(proptest::array::uniform3(-2i8..=2))
            .prop_map(|t| t.map(|p| p.map(|x| x as f64 * 0.5)))
            .prop_filter("non-degenerate", |t| projection_axis(*t).is_some())
    }

    proptest! {
        #[test]
        fn tri_tri_is_symmetric_and_rotation_invariant(a in tri(), b in tri(), r in 0usize..3, flip in proptest::bool::ANY) {
            let k = tri_tri(a, b);
            prop_assert_eq!(tri_tri(b, a), k);
            let mut rb = [b[r], b[(r + 1) % 3], b[(r + 2) % 3]];
            if flip {
                rb.swap(1, 2);
            }
            prop_assert_eq!(tri_tri(a, rb), k);
        }

        #[test]
        fn disjoint_iff_no_edge_meets_the_other(a in tri(), b in tri()) {
            let hits = (0..3).any(|e| seg_tri(a[e], a[(e + 1) % 3], b).0 != SegTri::Miss)
                || (0..3).any(|e| seg_tri(b[e], b[(e + 1) % 3], a).0 != SegTri::Miss);
            prop_assert_eq!(tri_tri(a, b) == TriTri::Disjoint, !hits);
        }
    }
}
