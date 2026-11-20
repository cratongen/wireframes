//! Failures found by adversarial mining scenarios against the rounding of
//! exact results to `f64`, each reduced to the two inputs of the operation
//! that failed, in `tests/data/attacks/<name>.txt`.

use num_rational::BigRational as Q;
use std::sync::mpsc;
use std::time::Duration;
use wireframes_boolean::{BooleanOptions, BooleanResult, Op, boolean, split_by_plane_report};
use wireframes_core::Mesh;

/// The two inputs in the format of `tests/data/rounding`: per mesh a line
/// `tag vertices faces`, then one vertex and one face per line.
fn pair(name: &str) -> (Mesh, Mesh) {
    let path = format!(
        "{}/tests/data/attacks/{name}.txt",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap();
    let mut lines = text.lines();
    let mut read = || {
        let head: Vec<usize> = lines
            .next()
            .unwrap()
            .split(' ')
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        let v: Vec<[f64; 3]> = (0..head[0])
            .map(|_| {
                let p: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split(' ')
                    .map(|x| x.parse().unwrap())
                    .collect();
                [p[0], p[1], p[2]]
            })
            .collect();
        let t: Vec<[u32; 3]> = (0..head[1])
            .map(|_| {
                let p: Vec<u32> = lines
                    .next()
                    .unwrap()
                    .split(' ')
                    .map(|x| x.parse().unwrap())
                    .collect();
                [p[0], p[1], p[2]]
            })
            .collect();
        Mesh::new(v, t).unwrap()
    };
    let a = read();
    (a, read())
}

/// The results of every attack pair do not depend on the thread count.
#[test]
fn attack_results_do_not_depend_on_thread_count() {
    let cases = [
        ("repeat_overlap", Op::Difference),
        ("repeat_tilted", Op::Difference),
        ("stope_minus_drive", Op::Difference),
        ("development_chain", Op::Union),
    ];
    for (name, op) in cases {
        let (a, b) = pair(name);
        let go = |n: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| boolean(&a, &b, op, &BooleanOptions::default()).unwrap())
        };
        let (r1, r4) = (go(1), go(4));
        assert_eq!(r1.mesh, r4.mesh, "{name}");
        assert_eq!(r1.face_source, r4.face_source, "{name}");
        assert_eq!(r1.rounding, r4.rounding, "{name}");
    }
}

fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}

fn volume6(m: &Mesh) -> Q {
    m.triangles().iter().fold(q(0.0), |acc, t| {
        let [a, b, c] = t.map(|i| m.vertices()[i as usize].map(q));
        acc + &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
            + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
            + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
    })
}

fn solid(m: &Mesh) -> bool {
    m.is_empty() || m.validate(0.0).unwrap().summary.is_solid()
}

/// `op` of the pair `name`, failing after `limit` instead of hanging.
fn run(name: &str, op: Op, limit: Duration) -> (Mesh, Mesh, BooleanResult) {
    let (a, b) = pair(name);
    assert!(solid(&a) && solid(&b), "{name}: inputs");
    let (tx, rx) = mpsc::channel();
    let (x, y) = (a.clone(), b.clone());
    std::thread::spawn(move || {
        let _ = tx.send(boolean(&x, &y, op, &BooleanOptions::default()));
    });
    let r = rx
        .recv_timeout(limit)
        .unwrap_or_else(|_| panic!("{name} {op:?}: no result within {limit:?}"))
        .unwrap_or_else(|e| panic!("{name} {op:?}: {e}"));
    (a, b, r)
}

/// The result is a solid within its move and step bounds.
fn certified(name: &str, r: &BooleanResult) {
    assert!(solid(&r.mesh), "{name}: not a solid");
    assert!(r.rounding.steps <= r.rounding.step_bound, "{name}: steps");
    for k in 0..3 {
        assert!(
            r.rounding.max_move[k] <= r.rounding.move_bound[k],
            "{name}: move"
        );
    }
}

const MINUTE: Duration = Duration::from_secs(60);

/// `op` and the intersection of the pair `name` are solids within their
/// bounds, and `vol(op) + vol(a & b)` equals `vol(a) + vol(b)` for a union or
/// `vol(a)` for a difference, exactly within the sum of the volume bounds.
fn identity_holds(name: &str, op: Op) {
    let (a, b, r) = run(name, op, 5 * MINUTE);
    certified(name, &r);
    let (_, _, i) = run(name, Op::Intersection, 5 * MINUTE);
    certified(name, &i);
    let exact = match op {
        Op::Union => volume6(&a) + volume6(&b),
        _ => volume6(&a),
    };
    let err = volume6(&r.mesh) + volume6(&i.mesh) - exact;
    let bound = q(6.0 * r.rounding.volume_bound) + q(6.0 * i.rounding.volume_bound);
    assert!(
        err.clone() <= bound && -err <= bound,
        "{name}: volume identity"
    );
}

/// A box united with sixteen lids tilted by 1e-13 to 1.6e-12 rad about a
/// horizontal line on its top, one at a time, then with a lid tilted by
/// 1.7e-12 rad.
#[test]
fn union_of_a_chain_of_lids_tilted_by_tiny_angles_finishes() {
    identity_holds("lid_chain", Op::Union);
}

/// A box minus thirteen overlapping caps tilted by alternating angles of
/// 1e-13 rad, one at a time, then minus a fourteenth.
#[test]
fn difference_of_a_chain_of_tilted_caps_finishes() {
    identity_holds("cut_chain", Op::Difference);
}

/// `x = (a | b) - b` for boxes rotated to strike and dip at UTM coordinates,
/// then `(x | b) - b` again: the second difference.
#[test]
fn repeated_union_then_difference_of_one_body_is_a_solid() {
    identity_holds("repeat_overlap", Op::Difference);
}

/// The same loop with `b` turned by a degree in strike and dip: the fifth
/// difference.
#[test]
fn repeated_difference_of_a_tilted_body_finishes() {
    identity_holds("repeat_tilted", Op::Difference);
}

/// Twenty-four face-sharing stopes rotated to strike 271.3 and dip 85 at UTM
/// coordinates, united, minus three drives, then minus a vertical raise.
#[test]
fn stope_block_minus_a_raise_after_drives_is_a_solid() {
    identity_holds("stope_minus_drive", Op::Difference);
}

/// A drive turned to strike 23.7 at UTM coordinates united with six
/// crosscuts and raises one at a time, then with a seventh.
#[test]
fn seventh_union_of_development_turned_to_strike_is_a_solid() {
    identity_holds("development_chain", Op::Union);
}

/// Two conical pit shells clipped by the same contour topography at UTM
/// coordinates: their union.
#[test]
fn union_of_pits_clipped_by_the_same_topography_finishes() {
    identity_holds("pits_sharing_topography", Op::Union);
}

/// A stope rotated to strike 33 and dip 70 at UTM coordinates, split by the
/// plane of one of its own faces (origin and normal on the last line).
#[test]
fn split_of_a_rotated_stope_along_its_own_face_is_two_solids() {
    let (m, _) = pair("split_along_face");
    let path = format!(
        "{}/tests/data/attacks/split_along_face.txt",
        env!("CARGO_MANIFEST_DIR")
    );
    let plane: Vec<f64> = std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .last()
        .unwrap()
        .split(' ')
        .map(|x| x.parse().unwrap())
        .collect();
    let (lo, hi, r) = split_by_plane_report(
        &m,
        [plane[0], plane[1], plane[2]],
        [plane[3], plane[4], plane[5]],
    )
    .unwrap_or_else(|e| panic!("split_along_face: {e}"));
    assert!(solid(&lo) && solid(&hi));
    for x in &r {
        assert!(x.steps <= x.step_bound);
        for k in 0..3 {
            assert!(x.max_move[k] <= x.move_bound[k]);
        }
    }
    let err = volume6(&lo) + volume6(&hi) - volume6(&m);
    let bound = q(6.0 * r[0].volume_bound) + q(6.0 * r[1].volume_bound);
    assert!(err.clone() <= bound && -err <= bound);
}
