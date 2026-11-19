//! Exact scaling by a power of two commutes with every operation whose output
//! points are dyadic, so results far from unit scale must equal the unit
//! results scaled. `HUGE` puts coordinates near 1e150; `SUBNORMAL` puts them
//! in the subnormal range, where products underflow.

use wireframes_boolean::{
    BooleanOptions, Op, boolean, clip_by_surface, slice_by_plane, split_by_plane,
};
use wireframes_core::{Mesh, Polylines, fixtures};

mod support;
use support::closed_manifold;

const HUGE: i32 = 498;
const SUBNORMAL: i32 = -1064;

fn ldexp(x: f64, k: i32) -> f64 {
    x * 2f64.powi(k / 2) * 2f64.powi(k - k / 2)
}

fn sp(p: [f64; 3], k: i32) -> [f64; 3] {
    p.map(|x| ldexp(x, k))
}

fn scaled(m: &Mesh, k: i32) -> Mesh {
    let v = m.vertices().iter().map(|&p| sp(p, k)).collect();
    m.clone().with_vertices(v).unwrap()
}

fn lines(p: &Polylines, k: i32) -> (Vec<[f64; 3]>, Vec<u32>) {
    let v = p.vertices().iter().map(|&x| sp(x, k)).collect();
    (v, p.part_offsets().to_vec())
}

/// Two boxes whose faces cross; with power-of-two extents every crossing,
/// including those with face diagonals, is a dyadic point.
fn boxes() -> (Mesh, Mesh) {
    (
        fixtures::box_mesh([0.0; 3], [4.0; 3]),
        fixtures::box_mesh([1.5, 0.5, -1.0], [5.5, 2.5, 3.0]),
    )
}

/// The second box sheared along x, so crossings are no longer dyadic and
/// only scales whose rounding commutes with scaling (normal range) apply.
fn skewed() -> (Mesh, Mesh) {
    let (a, b) = boxes();
    let v = b
        .vertices()
        .iter()
        .map(|p| [p[0] + 0.25 * p[2], p[1], p[2]])
        .collect();
    (a, b.with_vertices(v).unwrap())
}

fn cases() -> [((Mesh, Mesh), &'static [i32]); 2] {
    [(boxes(), &[HUGE, SUBNORMAL]), (skewed(), &[HUGE])]
}

/// An open sheet `z = 1 + x / 2` with normals up, wider than the box.
fn tilted_sheet() -> Mesh {
    let z = |x: f64| 1.0 + x / 2.0;
    Mesh::new(
        vec![
            [-4.0, -4.0, z(-4.0)],
            [8.0, -4.0, z(8.0)],
            [8.0, 8.0, z(8.0)],
            [-4.0, 8.0, z(-4.0)],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap()
}

#[test]
fn booleans_commute_with_extreme_scaling() {
    let opts = BooleanOptions::default();
    for ((a, b), ks) in cases() {
        for op in [Op::Union, Op::Intersection, Op::Difference] {
            let unit = boolean(&a, &b, op, &opts).unwrap().mesh;
            assert!(closed_manifold(&unit), "{op:?}");
            for &k in ks {
                let r = boolean(&scaled(&a, k), &scaled(&b, k), op, &opts).unwrap();
                assert_eq!(r.mesh, scaled(&unit, k), "{op:?} at 2^{k}");
            }
        }
    }
}

#[test]
fn split_and_slice_commute_with_extreme_scaling() {
    let (origin, normal) = ([0.0, 0.0, 1.0], [1.0, 0.0, 1.0]);
    for ((_, b), ks) in cases() {
        let unit = split_by_plane(&b, origin, normal).unwrap();
        assert!(closed_manifold(&unit.0) && closed_manifold(&unit.1));
        let section = slice_by_plane(&b, origin, [0.0, 0.0, 1.0], &[0.5, 1.25]).unwrap();
        assert_eq!(section.n_parts(), 2);
        for &k in ks {
            let s = scaled(&b, k);
            let want = (scaled(&unit.0, k), scaled(&unit.1, k));
            assert_eq!(split_by_plane(&s, sp(origin, k), normal).unwrap(), want);
            assert_eq!(
                split_by_plane(&s, sp(origin, k), sp(normal, k)).unwrap(),
                want
            );
            let offsets = [ldexp(0.5, k), ldexp(1.25, k)];
            let got = slice_by_plane(&s, sp(origin, k), [0.0, 0.0, 1.0], &offsets).unwrap();
            assert_eq!(lines(&got, 0), lines(&section, k), "2^{k}");
        }
    }
}

/// Cap flips are chosen on an `f64` approximation of the points (see
/// `triangulate`), and cocircular cap points tie differently once that
/// approximation lands on the coarse subnormal grid; validity is exact, so
/// vertices, face count and closure still match.
#[test]
fn clip_commutes_with_extreme_scaling() {
    let (a, _) = boxes();
    let sheet = tilted_sheet();
    let opts = BooleanOptions::default();
    let unit = clip_by_surface(&a, &sheet, &opts).unwrap();
    assert!(closed_manifold(&unit.below) && closed_manifold(&unit.above));
    let c = clip_by_surface(&scaled(&a, HUGE), &scaled(&sheet, HUGE), &opts).unwrap();
    assert_eq!(c.below, scaled(&unit.below, HUGE));
    assert_eq!(c.above, scaled(&unit.above, HUGE));
    let c = clip_by_surface(&scaled(&a, SUBNORMAL), &scaled(&sheet, SUBNORMAL), &opts).unwrap();
    for (got, want) in [(&c.below, &unit.below), (&c.above, &unit.above)] {
        let want = scaled(want, SUBNORMAL);
        assert_eq!(got.vertices(), want.vertices());
        assert_eq!(got.n_faces(), want.n_faces());
        assert!(closed_manifold(got));
    }
}
#[test]
fn split_of_a_sliver_near_overflow_is_closed() {
    let tet = Mesh::new(
        vec![
            [-1e150, 0.0, 0.0],
            [1e150, 0.0, 0.0],
            [1e150, 1e150, 0.0],
            [1e150, 0.0, 1e150],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
    .unwrap();
    let (neg, pos) = split_by_plane(&tet, [1e150 - 1e134, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
    assert!(closed_manifold(&neg) && closed_manifold(&pos));
}

#[test]
fn split_with_subnormal_coordinates_is_closed() {
    let d = f64::from_bits(1);
    let tet = Mesh::new(
        vec![
            [-3.0 * d, 0.0, 0.0],
            [5.0 * d, 0.0, 0.0],
            [5.0 * d, 7.0 * d, 0.0],
            [5.0 * d, 0.0, 9.0 * d],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
    .unwrap();
    let (neg, pos) = split_by_plane(&tet, [0.0; 3], [1.0, 0.0, 0.0]).unwrap();
    assert!(closed_manifold(&neg) && closed_manifold(&pos));
    let lines = slice_by_plane(&tet, [0.0; 3], [1.0, 0.0, 0.0], &[]).unwrap();
    assert_eq!(lines.n_parts(), 1);
}
