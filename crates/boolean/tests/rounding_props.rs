//! Seeded inputs built to hit each rounding mechanism of
//! `docs/design/rounding.md`: walls and caps that agree only to a few ulps,
//! solids rotated by tiny angles or to strike at UTM coordinates, and chains
//! of unions fed with the previous result. Every result must be a solid,
//! satisfy the volume identities against exact rational volumes within the
//! reported `volume_bound`s, and be the same on 1 and 8 threads.

use num_rational::BigRational as Q;
use proptest::prelude::*;
use wireframes_boolean::{BooleanOptions, BooleanResult, Op, boolean};
use wireframes_core::Mesh;

const BOX_FACES: [[u32; 3]; 12] = [
    [0, 2, 1],
    [1, 2, 3],
    [4, 5, 6],
    [5, 7, 6],
    [0, 1, 4],
    [1, 5, 4],
    [2, 6, 3],
    [3, 6, 7],
    [0, 4, 2],
    [2, 4, 6],
    [1, 3, 5],
    [3, 7, 5],
];

/// Corners in the order of `BOX_FACES`: bit 0 is x, bit 1 is y, bit 2 is z.
fn corners(lo: [f64; 3], hi: [f64; 3]) -> Vec<[f64; 3]> {
    (0..8)
        .map(|i| [0, 1, 2].map(|k| if i >> k & 1 == 1 { hi[k] } else { lo[k] }))
        .collect()
}

fn mesh(v: Vec<[f64; 3]>) -> Mesh {
    Mesh::new(v, BOX_FACES.to_vec()).unwrap()
}

/// `p` rotated by `angle` about the vertical through `c`, rounded to `f64`.
fn turn(p: [f64; 3], c: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, k) = angle.sin_cos();
    let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
    [c[0] + k * dx - s * dy, c[1] + s * dx + k * dy, p[2]]
}

/// `x` moved by `n` ulps.
fn nudge(mut x: f64, n: i32) -> f64 {
    for _ in 0..n.abs() {
        x = if n > 0 { x.next_up() } else { x.next_down() };
    }
    x
}

fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}

/// Six times the signed volume of `m`, exactly.
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

fn on_threads(
    n: usize,
    a: &Mesh,
    b: &Mesh,
    op: Op,
) -> std::result::Result<BooleanResult, TestCaseError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .unwrap()
        .install(|| boolean(a, b, op, &BooleanOptions::default()))
        .map_err(|e| TestCaseError::fail(format!("{op:?} on {n} threads: {e}")))
}

fn same(x: &Mesh, y: &Mesh) -> bool {
    let bits =
        |m: &Mesh| -> Vec<u64> { m.vertices().iter().flatten().map(|c| c.to_bits()).collect() };
    bits(x) == bits(y) && x.triangles() == y.triangles()
}

/// Union, intersection and difference of `a` and `b` are solids, identical
/// on 1 and 8 threads, with `vol(a | b) + vol(a & b) = vol(a) + vol(b)` and
/// `vol(a - b) + vol(a & b) = vol(a)` holding exactly within the sum of the
/// results' volume bounds. Returns the union.
fn check(a: &Mesh, b: &Mesh) -> std::result::Result<Mesh, TestCaseError> {
    prop_assert!(solid(a) && solid(b), "inputs");
    let mut out = vec![];
    for op in [Op::Union, Op::Intersection, Op::Difference] {
        let r = on_threads(1, a, b, op)?;
        prop_assert!(solid(&r.mesh), "{op:?}");
        prop_assert!(r.rounding.steps <= r.rounding.step_bound, "{op:?}");
        for k in 0..3 {
            prop_assert!(r.rounding.max_move[k] <= r.rounding.move_bound[k], "{op:?}");
        }
        prop_assert!(
            same(&r.mesh, &on_threads(8, a, b, op)?.mesh),
            "{op:?} threads"
        );
        out.push(r);
    }
    let [u, i, d] =
        [0, 1, 2].map(|k| (volume6(&out[k].mesh), q(6.0 * out[k].rounding.volume_bound)));
    let (va, vb) = (volume6(a), volume6(b));
    let abs = |x: Q| if x < q(0.0) { -x } else { x };
    prop_assert!(
        abs(&u.0 + &i.0 - &va - &vb) <= &u.1 + &i.1,
        "union and intersection"
    );
    prop_assert!(
        abs(&d.0 + &i.0 - &va) <= &d.1 + &i.1,
        "difference and intersection"
    );
    Ok(out.swap_remove(0).mesh)
}

fn utm() -> impl Strategy<Value = [f64; 3]> {
    (0u32..4000, 0u32..4000, 0u32..400).prop_map(|(x, y, z)| {
        [
            700_000.0 + f64::from(x) * 0.37,
            7_500_000.0 + f64::from(y) * 0.29,
            f64::from(z) * 1.3 - 200.0,
        ]
    })
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 12, ..ProptestConfig::default() })]

    /// Two stopes sharing a wall whose vertices on one side lie a few ulps off
    /// it, optionally rotated to strike, at UTM coordinates.
    #[test]
    fn stopes_sharing_a_wall_a_few_ulps_off(
        o in utm(),
        w in 4u32..16,
        dy in -6i32..10,
        dz in -8i32..12,
        nudges in proptest::array::uniform4(-3i32..=3),
        strike in prop_oneof![Just(0.0), 0.0..std::f64::consts::PI],
    ) {
        let (l, h, t) = (20.0, 16.0, 24.0);
        let a: Vec<[f64; 3]> = corners(o, [o[0] + l, o[1] + h, o[2] + t]);
        let lo = [o[0] + l, o[1] + f64::from(dy), o[2] + f64::from(dz)];
        let mut b = corners(lo, [lo[0] + f64::from(w), lo[1] + h, lo[2] + t]);
        for (k, i) in [0usize, 2, 4, 6].into_iter().enumerate() {
            b[i][0] = nudge(b[i][0], nudges[k]);
        }
        let (a, b) = if strike == 0.0 {
            (a, b)
        } else {
            (a.iter().map(|&p| turn(p, o, strike)).collect(), b.iter().map(|&p| turn(p, o, strike)).collect())
        };
        check(&mesh(a), &mesh(b))?;
    }

    /// Prisms whose tops lie on one sloping plane at UTM coordinates, the
    /// smaller one's top corners a few ulps above or below it.
    #[test]
    fn prisms_whose_caps_agree_to_a_few_ulps(
        o in utm(),
        slope in (-0.3f64..0.3, -0.3f64..0.3),
        bumps in proptest::array::uniform4(-2i32..=2),
        shift in (-3.0f64..3.0, -3.0f64..3.0),
    ) {
        let top = |x: f64, y: f64| o[2] + 40.0 + slope.0 * (x - o[0]) + slope.1 * (y - o[1]);
        let prism = |c: [f64; 2], half: f64, bottom: f64| -> Vec<[f64; 3]> {
            (0..8)
                .map(|i| {
                    let x = c[0] + if i & 1 == 1 { half } else { -half };
                    let y = c[1] + if i & 2 == 2 { half } else { -half };
                    [x, y, if i & 4 == 4 { top(x, y) } else { bottom }]
                })
                .collect()
        };
        let a = prism([o[0], o[1]], 10.0, o[2] - 10.0);
        let mut b = prism([o[0] + shift.0, o[1] + shift.1], 6.0, o[2] - 5.0);
        for (k, i) in [4usize, 5, 6, 7].into_iter().enumerate() {
            b[i][2] = nudge(b[i][2], bumps[k]);
        }
        check(&mesh(a), &mesh(b))?;
    }

    /// A box and its copy turned by a tiny angle about a vertical through
    /// its centre, so every face of one lies within ulps of a face of the
    /// other, at UTM coordinates.
    #[test]
    fn boxes_turned_by_tiny_angles(
        o in utm(),
        size in (5.0f64..30.0, 5.0f64..30.0, 5.0f64..30.0),
        angle in prop_oneof![1e-15f64..1e-12, 1e-12f64..1e-9],
        sign in prop_oneof![Just(1.0), Just(-1.0)],
    ) {
        let hi = [o[0] + size.0, o[1] + size.1, o[2] + size.2];
        let a = corners(o, hi);
        let c = [0, 1, 2].map(|k| (o[k] + hi[k]) / 2.0);
        let b: Vec<[f64; 3]> = a.iter().map(|&p| turn(p, c, sign * angle)).collect();
        check(&mesh(a), &mesh(b))?;
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 2, ..ProptestConfig::default() })]

    /// A drive at UTM coordinates turned to strike, united with twelve
    /// crosscuts one at a time, each union fed with the previous result.
    #[test]
    fn chains_of_twelve_unions_stay_solid(
        o in utm(),
        strike in 0.0..std::f64::consts::PI,
        offsets in proptest::collection::vec((0.0f64..1.0, -1.5f64..1.5), 12),
    ) {
        chain(o, strike, &offsets)?;
    }
}

fn chain(
    o: [f64; 3],
    strike: f64,
    offsets: &[(f64, f64)],
) -> std::result::Result<(), TestCaseError> {
    let drive = corners(o, [o[0] + 240.0, o[1] + 5.0, o[2] + 4.0]);
    let mut chain = mesh(drive.iter().map(|&p| turn(p, o, strike)).collect());
    for (k, &(along, dz)) in offsets.iter().enumerate() {
        let x = o[0] + 8.0 + 19.0 * k as f64 + 6.0 * along;
        let lo = [x, o[1] + 2.0, o[2] + 0.5 + dz];
        let cut = corners(lo, [x + 4.0, o[1] + 30.0, lo[2] + 3.0]);
        let cut = mesh(cut.iter().map(|&p| turn(p, o, strike)).collect());
        chain = check(&chain, &cut)?;
    }
    Ok(())
}

/// The chains the property test shrank to while the rounding still refused
/// them: crosscut walls crossing drive walls along vertical lines of points
/// whose x differ by a few ulps.
#[test]
fn chains_the_property_test_found_stay_solid() {
    let found: [([f64; 3], f64, [(f64, f64); 12]); 3] = [
        (
            [700000.0, 7500137.75, 225.10000000000002],
            0.7711623533638733,
            [
                (0.7772062315463661, -0.23551520071828266),
                (0.8868182763684063, -0.4860894529121129),
                (0.47181349132926753, 0.2111520263725513),
                (0.8297191793148249, -0.6973576600905739),
                (0.19283157630112815, 0.15853662503506666),
                (0.7302989146523207, -0.9417129987991567),
                (0.54919682045656, -0.11847482452357522),
                (0.188039017171867, 0.00598783124610597),
                (0.9001580934601976, 0.174069363345005),
                (0.21613949140318706, -0.43766031727896876),
                (0.44977751554131723, 0.7868503533382998),
                (0.014445784486408486, 1.2439203163572596),
            ],
        ),
        (
            [700000.0, 7500000.0, -200.0],
            2.606563183335732,
            [
                (0.5522236288088008, 0.0),
                (0.0, 0.0),
                (0.0, -1.0464448865643472),
                (0.0, 0.0),
                (0.0, 0.0),
                (0.0, 1.0784387140289027),
                (0.3700320023014421, 0.0),
                (0.006937443056138422, 1.3912527735140907),
                (0.8051794571750388, 1.1871808331749372),
                (0.28022479708850717, 0.0),
                (0.0, 0.0),
                (0.0, 0.0),
            ],
        ),
        (
            [700000.0, 7500751.1, -188.3],
            2.606563183335732,
            [
                (0.5522236288088008, 0.4734802105819862),
                (0.38254395752157716, -1.3974438979207637),
                (0.8448702262428947, -1.0464448865643472),
                (0.005894323729215422, -0.6412067102304502),
                (0.12406777733923953, -1.2106956813141987),
                (0.13311270070402323, 1.0784387140289027),
                (0.3700320023014421, 0.1149543156814529),
                (0.006937443056138422, 1.3912527735140907),
                (0.8051794571750388, 1.1871808331749372),
                (0.28022479708850717, -0.20948907870555716),
                (0.3188209833594534, 0.17390279371178335),
                (0.47128310128082523, -1.169940880657193),
            ],
        ),
    ];
    for (o, strike, offsets) in found {
        chain(o, strike, &offsets).unwrap();
    }
}

/// `p` turned to dip `dip` (degrees from horizontal, about the strike axis
/// `x`) and strike `strike` (radians about the vertical), then moved to `o`,
/// each step rounded to `f64`.
fn orient(p: [f64; 3], o: [f64; 3], strike: f64, dip: f64) -> [f64; 3] {
    let (st, ct) = (90.0 - dip).to_radians().sin_cos();
    let (y, z) = (p[1] * ct - p[2] * st, p[1] * st + p[2] * ct);
    let (ss, cs) = strike.sin_cos();
    [
        o[0] + p[0] * cs - y * ss,
        o[1] + p[0] * ss + y * cs,
        o[2] + z,
    ]
}

/// Stopes four along strike in two levels, all with their footwall on the
/// common wall `y = 0`, each with its own height and thickness and each upper
/// stope standing on the top of the one below, oriented and united one at a
/// time; every union is checked to be a solid.
fn stope_chain(
    o: [f64; 3],
    strike: f64,
    dip: f64,
    sizes: &[(f64, f64)],
) -> std::result::Result<Mesh, TestCaseError> {
    let w = 12.0;
    let mut chain: Option<Mesh> = None;
    for (k, &(height, thick)) in sizes.iter().enumerate() {
        let x = w * (k % 4) as f64;
        let z0 = if k >= 4 { sizes[k - 4].0 } else { 0.0 };
        let stope = corners([x, 0.0, z0], [x + w, thick, z0 + height]);
        let stope = mesh(stope.iter().map(|&p| orient(p, o, strike, dip)).collect());
        chain = Some(match chain {
            None => stope,
            Some(c) => {
                let r = on_threads(1, &c, &stope, Op::Union)?;
                prop_assert!(solid(&r.mesh), "union {k}");
                r.mesh
            }
        });
    }
    Ok(chain.unwrap())
}

/// `a - b` and `a & b` are solids within their bounds, the difference is the
/// same on 1 and 4 threads, and `vol(a - b) + vol(a & b) = vol(a)` holds
/// exactly within the sum of their volume bounds.
fn difference_holds(a: &Mesh, b: &Mesh) -> std::result::Result<(), TestCaseError> {
    let d = on_threads(1, a, b, Op::Difference)?;
    let i = on_threads(1, a, b, Op::Intersection)?;
    for (r, op) in [(&d, "difference"), (&i, "intersection")] {
        prop_assert!(solid(&r.mesh), "{op}");
        prop_assert!(r.rounding.steps <= r.rounding.step_bound, "{op}");
        for k in 0..3 {
            prop_assert!(r.rounding.max_move[k] <= r.rounding.move_bound[k], "{op}");
        }
    }
    prop_assert!(
        same(&d.mesh, &on_threads(4, a, b, Op::Difference)?.mesh),
        "difference threads"
    );
    let err = volume6(&d.mesh) + volume6(&i.mesh) - volume6(a);
    let bound = q(6.0 * d.rounding.volume_bound) + q(6.0 * i.rounding.volume_bound);
    let abs = if err < q(0.0) { -err } else { err };
    prop_assert!(abs <= bound, "difference and intersection");
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 6, ..ProptestConfig::default() })]

    /// A chain of eight stopes in two levels against a common wall at a dip
    /// from 45 to 90 degrees, turned to strike at UTM coordinates, minus a
    /// drive along strike on the floor of the lower level or through the
    /// common wall.
    #[test]
    fn stope_chains_minus_a_drive_along_strike(
        o in utm(),
        strike in 0.0..std::f64::consts::PI,
        dip in prop_oneof![Just(90.0), 45.0f64..90.0],
        sizes in proptest::collection::vec((18.0f64..30.0, 4.0f64..8.0), 8),
        through_wall in any::<bool>(),
    ) {
        let a = stope_chain(o, strike, dip, &sizes)?;
        let (lo, hi) = if through_wall {
            ([-6.0, -2.0, 6.0], [54.0, 2.5, 10.0])
        } else {
            ([-6.0, 0.5, -1.5], [54.0, 4.0, 2.5])
        };
        let drive = corners(lo, hi);
        let b = mesh(drive.iter().map(|&p| orient(p, o, strike, dip)).collect());
        difference_holds(&a, &b)?;
    }
}
