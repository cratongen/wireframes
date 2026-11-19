mod support;
use support::closed_manifold;

use proptest::prelude::*;
use wireframes_boolean::{BooleanOptions, Op, boolean};
use wireframes_core::fixtures::{self, broken};
use wireframes_core::topology::locate_point;
use wireframes_core::{Error, Mesh};

fn volume(m: &Mesh) -> f64 {
    let o = m.bounds().map_or([0.0; 3], |b| b.min);
    let mut s = 0.0;
    let mut c = 0.0;
    for f in 0..m.n_faces() as u32 {
        let [a, b, d] = m
            .corners(f)
            .map(|p| [p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
        let x = (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0])
            + a[2] * (b[0] * d[1] - b[1] * d[0]))
            / 6.0;
        let t = s + x;
        c += if f64::abs(s) >= f64::abs(x) {
            (s - t) + x
        } else {
            (x - t) + s
        };
        s = t;
    }
    s + c
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * 1f64.max(a.abs()).max(b.abs())
}

fn run(a: &Mesh, b: &Mesh, op: Op) -> Mesh {
    boolean(a, b, op, &BooleanOptions::default())
        .unwrap_or_else(|e| panic!("{op:?} failed: {e}"))
        .mesh
}

fn validates(m: &Mesh) -> bool {
    m.is_empty() || m.validate(0.0).unwrap().summary.is_solid()
}

fn rotate(p: [f64; 3], a: [f64; 3]) -> [f64; 3] {
    let (s, c) = (a.map(f64::sin), a.map(f64::cos));
    let x = [p[0], c[0] * p[1] - s[0] * p[2], s[0] * p[1] + c[0] * p[2]];
    let y = [c[1] * x[0] + s[1] * x[2], x[1], -s[1] * x[0] + c[1] * x[2]];
    [c[2] * y[0] - s[2] * y[1], s[2] * y[0] + c[2] * y[1], y[2]]
}

fn moved(m: &Mesh, angles: [f64; 3], d: [f64; 3]) -> Mesh {
    let v = m
        .vertices()
        .iter()
        .map(|p| {
            let r = rotate(*p, angles);
            [r[0] + d[0], r[1] + d[1], r[2] + d[2]]
        })
        .collect();
    m.clone().with_vertices(v).unwrap()
}

fn check_identities(a: &Mesh, b: &Mesh, tol: f64) {
    let (va, vb) = (volume(a), volume(b));
    let u = run(a, b, Op::Union);
    let i = run(a, b, Op::Intersection);
    let d = run(a, b, Op::Difference);
    let x = run(a, b, Op::SymmetricDifference);
    for m in [&u, &i, &d, &x] {
        assert!(closed_manifold(m), "output not closed manifold");
    }
    let (vu, vi, vd, vx) = (volume(&u), volume(&i), volume(&d), volume(&x));
    assert!(
        close(vu + vi, va + vb, tol),
        "union {vu} + inter {vi} != {va} + {vb}"
    );
    assert!(close(vd, va - vi, tol), "diff {vd} != {va} - {vi}");
    assert!(close(vx, vu - vi, tol), "symdiff {vx} != {vu} - {vi}");
    let u2 = run(b, a, Op::Union);
    let i2 = run(b, a, Op::Intersection);
    assert!(close(volume(&u2), vu, tol) && close(volume(&i2), vi, tol));
}

#[test]
fn validate_failure_rate_is_reported() {
    let mut rng = fixtures::Rng::new(7);
    let (mut runs, mut bad) = (0, 0);
    let (mut sym_runs, mut sym_touch, mut sym_bad) = (0, 0, 0);
    for _ in 0..40 {
        let a = fixtures::sphere(1);
        let ang = [0, 1, 2].map(|_| rng.uniform(-3.2, 3.2));
        let d = [0, 1, 2].map(|_| rng.uniform(-0.8, 0.8));
        let b = moved(&fixtures::cube(), ang, d);
        for op in [Op::Union, Op::Intersection, Op::Difference] {
            let m = run(&a, &b, op);
            assert!(closed_manifold(&m));
            runs += 1;
            bad += usize::from(!validates(&m));
        }
        let x = run(&a, &b, Op::SymmetricDifference);
        assert!(closed_manifold(&x));
        let s = x.validate(0.0).unwrap().summary;
        sym_runs += 1;
        sym_touch += usize::from(s.touching_shell_edges > 0);
        sym_bad += usize::from(!x.is_empty() && !s.is_solid());
    }
    println!("validate(0.0) failed on {bad} of {runs} rounded boolean outputs");
    println!(
        "symmetric difference: {sym_touch} of {sym_runs} have shells touching along the \
         intersection curve; {sym_bad} of {sym_runs} fail validate(0.0)"
    );
}

#[test]
fn overlapping_cubes() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [0.5, 0.25, 0.125]);
    check_identities(&a, &b, 1e-12);
    assert!(close(
        volume(&run(&a, &b, Op::Intersection)),
        0.5 * 0.75 * 0.875,
        1e-12
    ));
}

#[test]
fn coplanar_cases_have_exact_volumes() {
    let want = [
        ("shared_face", 2.0, 0.0),
        ("half_shift", 1.75, 0.25),
        ("nested", 1.0, 0.125),
        ("opposite_face", 1.25, 0.0),
    ];
    for (name, a, b) in fixtures::coplanar_pairs() {
        let (_, u, i) = want.iter().find(|w| w.0 == name).unwrap();
        let mu = run(&a, &b, Op::Union);
        let mi = run(&a, &b, Op::Intersection);
        assert!(closed_manifold(&mu) && closed_manifold(&mi), "{name}");
        assert!(
            close(volume(&mu), *u, 1e-12),
            "{name} union {}",
            volume(&mu)
        );
        assert!(
            close(volume(&mi), *i, 1e-12),
            "{name} intersection {}",
            volume(&mi)
        );
        check_identities(&a, &b, 1e-12);
    }
}

#[test]
fn identical_solids() {
    for a in [fixtures::cube(), fixtures::sphere(1), fixtures::torus(8, 6)] {
        let va = volume(&a);
        assert!(close(volume(&run(&a, &a, Op::Union)), va, 1e-12));
        assert!(close(volume(&run(&a, &a, Op::Intersection)), va, 1e-12));
        assert!(run(&a, &a, Op::Difference).is_empty());
        assert!(run(&a, &a, Op::SymmetricDifference).is_empty());
    }
}

#[test]
fn spheres_and_tori() {
    let a = fixtures::sphere(2);
    let b = moved(&fixtures::torus(12, 8), [0.3, 0.2, 0.1], [0.4, -0.2, 0.1]);
    check_identities(&a, &b, 1e-9);
}

#[test]
fn nested_shells_and_cavities() {
    let a = fixtures::nested_shells();
    let b = fixtures::box_mesh([1.5, 1.5, 1.5], [4.0, 2.5, 2.5]);
    check_identities(&a, &b, 1e-12);
    let i = run(&a, &b, Op::Intersection);
    assert!(close(volume(&i), 1.5 - 0.125, 1e-12));
}

#[test]
fn many_components() {
    let parts: Vec<Mesh> = (0..4)
        .map(|k| fixtures::translated(&fixtures::cube(), [2.0 * k as f64, 0.0, 0.0]))
        .collect();
    let refs: Vec<&Mesh> = parts.iter().collect();
    let a = Mesh::concat(&refs).unwrap();
    let b = fixtures::box_mesh([0.5, 0.5, 0.5], [6.5, 2.0, 0.75]);
    check_identities(&a, &b, 1e-12);
}

#[test]
fn touching_along_an_edge_gives_separate_shells() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [1.0, 1.0, 0.0]);
    let u = run(&a, &b, Op::Union);
    let s = u.validate(0.0).unwrap().summary;
    assert!(closed_manifold(&u) && s.is_solid(), "{s:?}");
    assert_eq!(
        (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
        (2, 1, 2)
    );
    assert!(close(volume(&u), 2.0, 1e-12));
    assert!(run(&a, &b, Op::Intersection).is_empty());
    assert!(close(volume(&run(&a, &b, Op::Difference)), 1.0, 1e-12));
    let x = run(&a, &b, Op::SymmetricDifference);
    assert!(close(volume(&x), 2.0, 1e-12) && validates(&x));
}

#[test]
fn invalid_inputs_are_refused_with_the_report() {
    let cube = fixtures::cube();
    for bad in [fixtures::open_cube(), broken::overlapping_cubes()] {
        match boolean(&bad, &cube, Op::Union, &BooleanOptions::default()) {
            Err(Error::InvalidMesh { context, report }) => {
                assert_eq!(context, "a");
                assert!(!report.summary.is_solid());
            }
            other => panic!("expected InvalidMesh, got {other:?}"),
        }
    }
    let dup = broken::duplicate_vertex();
    assert!(boolean(&dup, &cube, Op::Union, &BooleanOptions::default()).is_err());
    let r = boolean(&dup, &cube, Op::Union, &BooleanOptions { repair: true }).unwrap();
    assert!(r.repairs[0].is_some() && r.repairs[1].is_some());
}

#[test]
fn utm_coordinates_give_the_same_volumes() {
    let a = fixtures::sphere(1);
    let b = moved(&fixtures::cube(), [0.1, 0.2, 0.3], [0.2, 0.1, -0.3]);
    let shift = [700_000.0, 7_500_000.0, 400.0];
    let (ua, ub) = (
        fixtures::translated(&a, shift),
        fixtures::translated(&b, shift),
    );
    for op in [Op::Union, Op::Intersection, Op::Difference] {
        let (v0, v1) = (volume(&run(&a, &b, op)), volume(&run(&ua, &ub, op)));
        assert!(close(v0, v1, 1e-6), "{op:?}: {v0} vs {v1}");
    }
}

#[test]
fn results_do_not_depend_on_thread_count() {
    let a = fixtures::sphere(2);
    let b = moved(&fixtures::torus(10, 7), [0.5, 0.1, 0.2], [0.3, 0.0, 0.1]);
    let go = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| boolean(&a, &b, Op::Difference, &BooleanOptions::default()).unwrap())
    };
    let (r1, r4) = (go(1), go(4));
    assert_eq!(r1.mesh, r4.mesh);
    assert_eq!(r1.face_source, r4.face_source);
}

#[test]
fn face_attributes_follow_their_source() {
    use wireframes_core::Column;
    let a = fixtures::cube()
        .with_face_column("rock", Column::i64(vec![1; 12]))
        .unwrap();
    let b = fixtures::translated(&fixtures::cube(), [0.5, 0.5, 0.5])
        .with_face_column("rock", Column::i64(vec![2; 12]))
        .unwrap();
    let r = boolean(&a, &b, Op::Union, &BooleanOptions::default()).unwrap();
    assert_eq!(r.face_source.len(), r.mesh.n_faces());
    assert_eq!(r.mesh.face_attributes().rows(), r.mesh.n_faces());
}

fn grid() -> impl Strategy<Value = f64> {
    prop_oneof![(-4i32..5).prop_map(|i| i as f64 * 0.25), -1.0..1.0f64]
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

    #[test]
    fn volume_identities_under_rigid_motions(
        angles in proptest::array::uniform3(prop_oneof![Just(0.0), -3.2..3.2f64]),
        d in proptest::array::uniform3(grid()),
        which in 0usize..3,
    ) {
        let a = match which {
            0 => fixtures::cube(),
            1 => fixtures::sphere(1),
            _ => fixtures::box_mesh([0.0; 3], [1.0, 0.5, 0.75]),
        };
        let b = moved(&fixtures::cube(), angles, d);
        let v = [Op::Union, Op::Intersection, Op::Difference, Op::SymmetricDifference]
            .iter()
            .map(|&op| boolean(&a, &b, op, &BooleanOptions::default()).map(|r| r.mesh))
            .collect::<Result<Vec<Mesh>, _>>()
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        for m in &v {
            prop_assert!(closed_manifold(m));
        }
        let (va, vb) = (volume(&a), volume(&b));
        let [vu, vi, vd, vx] = [0, 1, 2, 3].map(|k| volume(&v[k]));
        prop_assert!(close(vu + vi, va + vb, 1e-9), "{} {} {} {}", vu, vi, va, vb);
        prop_assert!(close(vd, va - vi, 1e-9));
        prop_assert!(close(vx, vu - vi, 1e-9));
    }

    #[test]
    fn exactly_touching_and_coplanar_boxes(
        lo in proptest::array::uniform3(-2i32..3),
        size in proptest::array::uniform3(1i32..4),
    ) {
        let a = fixtures::box_mesh([0.0; 3], [2.0, 2.0, 2.0]);
        let lo = lo.map(|x| x as f64 * 0.5);
        let hi = [0, 1, 2].map(|k| lo[k] + size[k] as f64 * 0.5);
        let b = fixtures::box_mesh(lo, hi);
        let ov: f64 = (0..3).map(|k| (hi[k].min(2.0) - lo[k].max(0.0)).max(0.0)).product();
        let (va, vb) = (8.0, volume(&b));
        for (op, want) in [
            (Op::Intersection, ov),
            (Op::Union, va + vb - ov),
            (Op::Difference, va - ov),
        ] {
            let m = boolean(&a, &b, op, &BooleanOptions::default())
                .map_err(|e| TestCaseError::fail(format!("{op:?}: {e}")))?
                .mesh;
            prop_assert!(closed_manifold(&m), "{:?}", op);
            prop_assert!(close(volume(&m), want, 1e-12), "{:?}: {} vs {}", op, volume(&m), want);
        }
    }

    #[test]
    fn boolean_matches_pointwise_location(
        angles in proptest::array::uniform3(-3.2..3.2f64),
        d in proptest::array::uniform3(-0.6..0.6f64),
        seed in 0u64..1000,
    ) {
        let a = fixtures::sphere(1);
        let b = moved(&fixtures::cube(), angles, d);
        let pts = fixtures::random_points(seed, 40, [-1.8; 3], [1.8; 3]);
        let loc = |m: &Mesh, p| locate_point(m.vertices(), m.triangles(), p);
        for (op, f) in [
            (Op::Union, (|x: bool, y: bool| x || y) as fn(bool, bool) -> bool),
            (Op::Intersection, |x, y| x && y),
            (Op::Difference, |x, y| x && !y),
            (Op::SymmetricDifference, |x, y| x != y),
        ] {
            let m = run(&a, &b, op);
            for &p in &pts {
                let (la, lb) = (loc(&a, p), loc(&b, p));
                let lm = if m.is_empty() { 1 } else { loc(&m, p) };
                if la == 0 || lb == 0 || lm == 0 {
                    continue;
                }
                prop_assert_eq!(lm == -1, f(la == -1, lb == -1), "{:?} at {:?}", op, p);
            }
        }
    }

    #[test]
    fn union_and_intersection_commute_in_location(
        angles in proptest::array::uniform3(-3.2..3.2f64),
        d in proptest::array::uniform3(-0.6..0.6f64),
        seed in 0u64..1000,
    ) {
        let a = fixtures::sphere(1);
        let b = moved(&fixtures::cube(), angles, d);
        let pts = fixtures::random_points(seed, 40, [-1.8; 3], [1.8; 3]);
        let loc = |m: &Mesh, p| if m.is_empty() { 1 } else { locate_point(m.vertices(), m.triangles(), p) };
        for op in [Op::Union, Op::Intersection] {
            let (ab, ba) = (run(&a, &b, op), run(&b, &a, op));
            for &p in &pts {
                let (l1, l2) = (loc(&ab, p), loc(&ba, p));
                if l1 == 0 || l2 == 0 {
                    continue;
                }
                prop_assert_eq!(l1, l2, "{:?} at {:?}", op, p);
            }
        }
    }

    #[test]
    fn tetrahedra_sharing_faces_edges_and_vertices(
        p in proptest::array::uniform4(proptest::array::uniform3(-2.0..2.0f64)),
        q in proptest::array::uniform3(-2.0..2.0f64),
        r in proptest::array::uniform3(-2.0..2.0f64),
        share in 0usize..3,
    ) {
        let a = tet([p[0], p[1], p[2], p[3]]);
        let b = match share {
            0 => tet([p[0], p[1], p[2], q]),
            1 => tet([p[0], p[1], q, r]),
            _ => tet([p[0], q, r, p[3]]),
        };
        let (Some(a), Some(b)) = (a, b) else { return Ok(()) };
        let mut vols = vec![];
        for op in [Op::Union, Op::Intersection, Op::Difference, Op::SymmetricDifference] {
            let m = boolean(&a, &b, op, &BooleanOptions::default())
                .map_err(|e| TestCaseError::fail(format!("{op:?}: {e}")))?
                .mesh;
            prop_assert!(closed_manifold(&m), "{:?}", op);
            vols.push(volume(&m));
        }
        let (va, vb) = (volume(&a), volume(&b));
        let [vu, vi, vd, vx] = [vols[0], vols[1], vols[2], vols[3]];
        prop_assert!(close(vu + vi, va + vb, 1e-9), "{} + {} vs {} + {}", vu, vi, va, vb);
        prop_assert!(close(vd, va - vi, 1e-9), "{} vs {} - {}", vd, va, vi);
        prop_assert!(close(vx, vu - vi, 1e-9), "{} vs {} - {}", vx, vu, vi);
    }
}

/// A positively oriented tetrahedron, `None` when flat.
fn tet(p: [[f64; 3]; 4]) -> Option<Mesh> {
    use wireframes_core::predicates::{Sign, side_of_plane};
    let t = match side_of_plane(p[3], p[0], p[1], p[2]) {
        Sign::Zero => return None,
        Sign::Negative => vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]],
        Sign::Positive => vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
    };
    let m = Mesh::new(p.to_vec(), t).ok()?;
    m.validate(0.0).ok()?.summary.is_solid().then_some(m)
}

fn hexahedron(x: [[f64; 2]; 2], y: [f64; 2], z: [f64; 2]) -> Mesh {
    let b = fixtures::box_mesh([0.0; 3], [1.0; 3]);
    let v: Vec<[f64; 3]> = (0..8)
        .map(|i| {
            let (xi, yi, zi) = (i & 1, (i >> 1) & 1, (i >> 2) & 1);
            [x[zi][xi], y[yi], z[zi]]
        })
        .collect();
    Mesh::new(v, b.triangles().to_vec()).unwrap()
}

/// B crosses A's top face in a strip `x in (0.625 + 2^-55, 0.625 + 2^-53)`
/// that holds no `f64`.
fn sliver_crossing_unit_cube() -> (Mesh, Mesh) {
    let u = 2f64.powi(-55);
    let b = hexahedron(
        [[0.5 - 2.0 * u, 0.5], [0.75 + 4.0 * u, 0.75 + 8.0 * u]],
        [0.2, 0.8],
        [0.5, 1.5],
    );
    (fixtures::cube(), b)
}

#[test]
fn sliver_inputs_are_valid_solids() {
    let (a, b) = sliver_crossing_unit_cube();
    assert!(a.validate(0.0).unwrap().summary.is_solid());
    assert!(b.validate(0.0).unwrap().summary.is_solid());
}

#[test]
fn boolean_with_a_patch_thinner_than_f64_spacing_succeeds() {
    let (a, b) = sliver_crossing_unit_cube();
    for op in [Op::Union, Op::Intersection, Op::Difference] {
        let r = boolean(&a, &b, op, &BooleanOptions::default());
        assert!(
            r.is_ok(),
            "{op:?} on valid solids failed: {}",
            r.err().unwrap()
        );
    }
    check_identities(&a, &b, 1e-12);
    check_identities(&b, &a, 1e-12);
}

#[test]
fn symmetric_difference_pieces_touch_only_along_the_curve() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [0.5, 0.25, 0.125]);
    let m = run(&a, &b, Op::SymmetricDifference);
    let s = m.validate(0.0).unwrap().summary;
    assert!(closed_manifold(&m) && s.is_closed, "{s:?}");
    assert_eq!(
        (
            s.non_manifold_edges,
            s.non_manifold_vertices,
            s.inward_shells
        ),
        (0, 0, 0),
        "{s:?}"
    );
    assert_eq!((s.self_intersections, s.coplanar_overlaps), (0, 0), "{s:?}");
    assert_eq!(s.shells, 2, "{s:?}");
    let d = [run(&a, &b, Op::Difference), run(&b, &a, Op::Difference)];
    let curve: std::collections::HashSet<[u64; 3]> = d[0]
        .vertices()
        .iter()
        .filter(|p| d[1].vertices().contains(p))
        .map(|p| p.map(f64::to_bits))
        .collect();
    assert_eq!(s.duplicate_vertices, 0, "{s:?}");
    assert_eq!(s.touching_shell_vertices, curve.len(), "{s:?}");
    assert!(s.touching_shell_edges > 0 && s.is_solid(), "{s:?}");
}

#[test]
fn symmetric_difference_of_overlapping_cubes_is_a_valid_solid() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [0.5, 0.25, 0.125]);
    let m = boolean(&a, &b, Op::SymmetricDifference, &BooleanOptions::default())
        .unwrap()
        .mesh;
    let rep = m.validate(0.0).unwrap();
    assert!(rep.summary.is_solid(), "summary: {:?}", rep.summary);
    let c = fixtures::translated(&a, [-0.3, -0.3, -0.3]);
    let r = boolean(&m, &c, Op::Union, &BooleanOptions::default());
    assert!(r.is_ok(), "reusing the result failed: {}", r.err().unwrap());
}

#[test]
fn an_empty_result_chains_into_further_booleans_and_clips() {
    let a = fixtures::cube();
    let empty = run(
        &a,
        &fixtures::translated(&a, [5.0, 0.0, 0.0]),
        Op::Intersection,
    );
    assert!(empty.is_empty());
    for (x, y, op, v) in [
        (&a, &empty, Op::Union, 1.0),
        (&empty, &a, Op::Union, 1.0),
        (&a, &empty, Op::Difference, 1.0),
        (&a, &empty, Op::SymmetricDifference, 1.0),
        (&a, &empty, Op::Intersection, 0.0),
        (&empty, &a, Op::Difference, 0.0),
        (&empty, &empty, Op::Union, 0.0),
    ] {
        let m = run(x, y, op);
        assert!(validates(&m) && close(volume(&m), v, 1e-12), "{op:?}");
    }
    let ground = fixtures::pit_surface(5, 1);
    let r =
        wireframes_boolean::clip_by_surface(&empty, &ground, &BooleanOptions::default()).unwrap();
    assert!(r.below.is_empty() && r.above.is_empty());
    let e = boolean(
        &fixtures::open_cube(),
        &a,
        Op::Union,
        &BooleanOptions::default(),
    );
    assert!(matches!(e, Err(Error::InvalidMesh { .. })));
}
