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

/// The split of the pair `name` by the plane on its last line is two solids
/// within their bounds whose volumes add up to the input within the summed
/// volume bounds.
fn split_holds(name: &str) {
    let (m, _) = pair(name);
    let path = format!(
        "{}/tests/data/attacks/{name}.txt",
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
    .unwrap_or_else(|e| panic!("{name}: {e}"));
    assert!(solid(&lo) && solid(&hi), "{name}: not solids");
    for x in &r {
        assert!(x.steps <= x.step_bound, "{name}: steps");
        for k in 0..3 {
            assert!(x.max_move[k] <= x.move_bound[k], "{name}: move");
        }
    }
    let err = volume6(&lo) + volume6(&hi) - volume6(&m);
    let bound = q(6.0 * r[0].volume_bound) + q(6.0 * r[1].volume_bound);
    assert!(
        err.clone() <= bound && -err <= bound,
        "{name}: volume identity"
    );
}

/// `x = (a - b)` for boxes rotated to strike 31 and dip 77 at UTM
/// coordinates, `b` turned by 0.7 degrees, then `x = (x | b) - b` four
/// times: the fifth union.
#[test]
fn fifth_union_of_a_repeated_union_then_difference_is_a_solid() {
    identity_holds("repeat_fifth_union", Op::Union);
}

/// A seeded pair of boxes rotated to strike and dip, `b` turned by under a
/// degree, `x = (a - b)`, then `x = (x | b) - b`: the third difference.
#[test]
fn third_difference_of_a_turned_body_is_a_solid() {
    identity_holds("repeat_turned_third_difference", Op::Difference);
}

/// Sixteen grid stopes against a common wall rotated to strike 271.3 and
/// dip 85 at UTM coordinates, united one at a time, minus a drive along
/// strike through all of them.
#[test]
#[ignore = "attack: stope chain minus a drive along strike raises KernelError for union, difference and intersection"]
fn stope_chain_minus_a_drive_along_strike_is_a_solid() {
    identity_holds("stope_chain_minus_drive", Op::Difference);
}

/// Eight grid stopes rotated to strike 118.4 and dip 62.5 at UTM
/// coordinates, united one at a time, minus a drive along strike.
#[test]
#[ignore = "attack: stopes dipping 62.5 minus a drive along strike raises KernelError (self_intersections=231)"]
fn stopes_dipping_62_minus_a_drive_is_a_solid() {
    identity_holds("stopes_dipping_62_minus_drive", Op::Difference);
}

/// The same block minus a drive along strike on the stope floor.
#[test]
#[ignore = "attack: stopes dipping 62.5 minus a drive on their floor raises KernelError (self_intersections=262)"]
fn stopes_dipping_62_minus_a_footwall_drive_is_a_solid() {
    identity_holds("stopes_dipping_62_minus_footwall_drive", Op::Difference);
}

/// A 2 x 2 x 3 grid stope rotated to strike 33 and dip 70 at UTM coordinates,
/// elevation -400, split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_a_is_two_solids() {
    split_holds("split_face_a");
}

/// A seeded grid stope rotated to strike and dip at UTM coordinates near
/// elevation 15, split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_b_is_two_solids() {
    split_holds("split_face_b");
}

/// A seeded grid stope rotated to strike and dip at UTM coordinates near
/// elevation -400, split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_c_is_two_solids() {
    split_holds("split_face_c");
}

/// A seeded grid stope rotated to strike and dip at UTM coordinates near
/// elevation 4800, split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_d_is_two_solids() {
    split_holds("split_face_d");
}

/// The stope of `split_face_d` split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_e_is_two_solids() {
    split_holds("split_face_e");
}

/// Another seeded grid stope rotated to strike and dip at UTM coordinates
/// near elevation 4800, split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_f_is_two_solids() {
    split_holds("split_face_f");
}

/// The stope of `split_face_f` split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_g_is_two_solids() {
    split_holds("split_face_g");
}

/// The stope of `split_face_f` split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_h_is_two_solids() {
    split_holds("split_face_h");
}

/// The stope of `split_face_f` split by the plane through the first vertex of one of its
/// faces with that face's f64 normal.
#[test]
fn split_along_own_face_i_is_two_solids() {
    split_holds("split_face_i");
}
