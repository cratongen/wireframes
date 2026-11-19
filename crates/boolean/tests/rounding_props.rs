//! Chains of unions at UTM coordinates fed with the previous result. Every
//! result must be a solid, satisfy the volume identities against exact
//! rational volumes within the reported `volume_bound`s, and be the same on
//! 1 and 8 threads.

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
