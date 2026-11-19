mod support;
use support::closed_manifold;

use proptest::prelude::*;
use wireframes_boolean::{BooleanOptions, Op, boolean, split_by_plane};
use wireframes_core::Mesh;
use wireframes_core::fixtures;

fn volume(m: &Mesh) -> f64 {
    let o = m.bounds().map_or([0.0; 3], |b| b.min);
    (0..m.n_faces() as u32)
        .map(|f| {
            let [a, b, d] = m
                .corners(f)
                .map(|p| [p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
            (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0])
                + a[2] * (b[0] * d[1] - b[1] * d[0]))
                / 6.0
        })
        .sum()
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * 1f64.max(a.abs()).max(b.abs())
}

fn run(a: &Mesh, b: &Mesh, op: Op) -> Mesh {
    boolean(a, b, op, &BooleanOptions::default())
        .unwrap_or_else(|e| panic!("{op:?} failed: {e}"))
        .mesh
}

fn ok(m: &Mesh) -> bool {
    m.is_empty() || closed_manifold(m)
}

fn solid(m: &Mesh) -> bool {
    m.is_empty() || m.validate(0.0).unwrap().summary.is_solid()
}

/// The same box with every quad split on the other diagonal.
fn other_diagonals(min: [f64; 3], max: [f64; 3]) -> Mesh {
    let b = fixtures::box_mesh(min, max);
    let t = vec![
        [0, 2, 1],
        [2, 3, 1],
        [4, 5, 6],
        [5, 7, 6],
        [0, 1, 4],
        [1, 5, 4],
        [2, 6, 3],
        [6, 7, 3],
        [0, 4, 2],
        [4, 6, 2],
        [1, 3, 5],
        [3, 7, 5],
    ];
    Mesh::new(b.vertices().to_vec(), t).unwrap()
}

/// Parallelepiped spanned by integer edge vectors from `o`, outward.
fn parallelepiped(o: [f64; 3], u: [f64; 3], v: [f64; 3], w: [f64; 3]) -> Mesh {
    let b = fixtures::cube();
    let p: Vec<[f64; 3]> = b
        .vertices()
        .iter()
        .map(|c| [0, 1, 2].map(|k| o[k] + c[0] * u[k] + c[1] * v[k] + c[2] * w[k]))
        .collect();
    Mesh::new(p, b.triangles().to_vec()).unwrap()
}

fn concat(parts: &[Mesh]) -> Mesh {
    let refs: Vec<&Mesh> = parts.iter().collect();
    Mesh::concat(&refs).unwrap()
}

fn check_all(a: &Mesh, b: &Mesh, u: f64, i: f64, tol: f64) {
    let (va, vb) = (volume(a), volume(b));
    let mu = run(a, b, Op::Union);
    let mi = run(a, b, Op::Intersection);
    let md = run(a, b, Op::Difference);
    let mx = run(a, b, Op::SymmetricDifference);
    for (n, m) in [
        ("union", &mu),
        ("inter", &mi),
        ("diff", &md),
        ("symdiff", &mx),
    ] {
        assert!(ok(m), "{n} not closed manifold");
    }
    assert!(close(volume(&mu), u, tol), "union {} != {u}", volume(&mu));
    assert!(close(volume(&mi), i, tol), "inter {} != {i}", volume(&mi));
    assert!(
        close(volume(&md), va - i, tol),
        "diff {} != {}",
        volume(&md),
        va - i
    );
    assert!(
        close(volume(&mx), u - i, tol),
        "symdiff {} != {}",
        volume(&mx),
        u - i
    );
    let _ = vb;
}

#[test]
fn identical_cubes_with_other_diagonals() {
    let a = fixtures::cube();
    let b = other_diagonals([0.0; 3], [1.0; 3]);
    check_all(&a, &b, 1.0, 1.0, 1e-12);
    assert!(run(&a, &b, Op::Difference).is_empty());
    assert!(run(&a, &b, Op::SymmetricDifference).is_empty());
    assert!(solid(&run(&a, &b, Op::Union)));
}

#[test]
fn half_shift_with_other_diagonals() {
    let a = fixtures::cube();
    let b = other_diagonals([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]);
    check_all(&a, &b, 1.75, 0.25, 1e-12);
}

#[test]
fn filling_the_cavity_gives_the_full_box() {
    let a = fixtures::nested_shells();
    let b = fixtures::box_mesh([1.0; 3], [2.0; 3]);
    let u = run(&a, &b, Op::Union);
    assert!(ok(&u));
    assert!(close(volume(&u), 27.0, 1e-12), "union {}", volume(&u));
    assert_eq!(u.n_faces(), 12, "internal cavity faces left in the union");
    assert!(run(&a, &b, Op::Intersection).is_empty());
    assert!(close(volume(&run(&a, &b, Op::Difference)), 26.0, 1e-12));
    assert!(close(volume(&run(&b, &a, Op::Difference)), 1.0, 1e-12));
}

#[test]
fn box_through_cavity_and_top_face() {
    let a = fixtures::nested_shells();
    let b = fixtures::box_mesh([1.0, 1.0, 1.0], [2.0, 2.0, 3.0]);
    check_all(&a, &b, 27.0, 1.0, 1e-12);
}

#[test]
fn tilted_coplanar_parallelepipeds() {
    let (u, v, w) = ([3.0, 4.0, 0.0], [-4.0, 3.0, 0.0], [0.0, 0.0, 5.0]);
    let a = parallelepiped([0.0; 3], u, v, w);
    let b = parallelepiped([1.5, 2.0, 0.0], u, v, w);
    let b2 = parallelepiped([1.5 - 2.0, 2.0 + 1.5, 2.5], u, v, w);
    check_all(&a, &b, 125.0 * 1.5, 125.0 * 0.5, 1e-12);
    check_all(
        &a,
        &b2,
        125.0 * 2.0 - 125.0 * 0.5 * 0.5 * 0.5,
        125.0 * 0.125,
        1e-12,
    );
}

#[test]
fn utm_sized_coplanar_boxes() {
    let x0 = [712_345.37, 7_481_234.91, 812.13];
    let s = [0.37, 0.41, 0.29];
    let a = fixtures::box_mesh(x0, [0, 1, 2].map(|k| x0[k] + 1.0));
    let lo = [0, 1, 2].map(|k| x0[k] + s[k]);
    let hi = [lo[0] + 1.0, lo[1] + 1.0, x0[2] + 1.0];
    let b = fixtures::box_mesh([lo[0], lo[1], x0[2]], hi);
    let ov: f64 = (0..2).map(|k| (x0[k] + 1.0) - lo[k]).product::<f64>() * 1.0;
    let (va, vb) = (volume(&a), volume(&b));
    check_all(&a, &b, va + vb - ov, ov, 1e-9);
}

#[test]
fn many_cubes_coplanar_with_a_slab() {
    let mut parts = vec![];
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..2 {
                parts.push(fixtures::translated(
                    &fixtures::cube(),
                    [2.0 * i as f64, 2.0 * j as f64, 2.0 * k as f64],
                ));
            }
        }
    }
    let a = concat(&parts);
    let b = fixtures::box_mesh([0.5, 0.0, 1.0], [8.5, 9.0, 2.0]);
    check_all(&a, &b, 50.0 + 72.0, 0.0, 1e-12);
    let b = fixtures::box_mesh([0.5, 0.0, 0.5], [8.5, 9.0, 1.0]);
    check_all(&a, &b, 50.0 + 36.0 - 10.0, 10.0, 1e-12);
}

#[test]
fn cubes_touching_at_a_vertex() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [1.0, 1.0, 1.0]);
    let e = boolean(&a, &b, Op::Union, &BooleanOptions::default()).unwrap_err();
    assert!(e.to_string().contains("not manifold"), "{e}");
    assert!(run(&a, &b, Op::Intersection).is_empty());
    assert!(close(volume(&run(&a, &b, Op::Difference)), 1.0, 1e-12));
    let x = run(&a, &b, Op::SymmetricDifference);
    assert!(ok(&x) && close(volume(&x), 2.0, 1e-12));
}

#[test]
fn wedge_edge_lying_on_a_face() {
    let a = fixtures::cube();
    let b = Mesh::new(
        vec![
            [0.2, 0.5, 1.0],
            [0.8, 0.5, 1.0],
            [0.2, 0.2, 1.5],
            [0.8, 0.2, 1.5],
            [0.2, 0.8, 1.5],
            [0.8, 0.8, 1.5],
        ],
        vec![
            [0, 1, 3],
            [0, 3, 2],
            [0, 4, 5],
            [0, 5, 1],
            [2, 3, 5],
            [2, 5, 4],
            [0, 2, 4],
            [1, 5, 3],
        ],
    )
    .unwrap();
    assert!(b.validate(0.0).unwrap().summary.is_solid());
    assert!(run(&a, &b, Op::Intersection).is_empty());
    let d = run(&a, &b, Op::Difference);
    assert!(ok(&d) && close(volume(&d), 1.0, 1e-12));
    let bd = run(&b, &a, Op::Difference);
    assert!(ok(&bd) && close(volume(&bd), volume(&b), 1e-12));
}

#[test]
fn needle_sliver_tetrahedron_through_a_cube() {
    let h = 1e-12;
    let t = Mesh::new(
        vec![
            [-0.5, 0.3, 0.4],
            [1.5, 0.31, 0.4],
            [0.5, 0.7, 0.4],
            [0.5, 0.5, 0.4 + h],
        ],
        vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
    )
    .unwrap();
    assert!(t.validate(0.0).unwrap().summary.is_solid());
    let a = fixtures::cube();
    let vi = volume(&run(&a, &t, Op::Intersection));
    let vu = volume(&run(&a, &t, Op::Union));
    assert!(close(vu + vi, 1.0 + volume(&t), 1e-12), "{vu} + {vi}");
    assert!(vi > 0.0 && vi < volume(&t));
}

#[test]
fn split_on_a_face_plane() {
    let a = fixtures::cube();
    for (z, neg_v) in [(0.0, 0.0), (1.0, 1.0)] {
        let (neg, pos) = split_by_plane(&a, [0.0, 0.0, z], [0.0, 0.0, 1.0]).unwrap();
        assert!(ok(&neg) && ok(&pos), "z = {z}");
        assert!(
            close(volume(&neg), neg_v, 1e-12),
            "z = {z}: neg {}",
            volume(&neg)
        );
        assert!(
            close(volume(&pos), 1.0 - neg_v, 1e-12),
            "z = {z}: pos {}",
            volume(&pos)
        );
        for m in [&neg, &pos] {
            assert!(
                m.is_empty() || volume(m) > 0.5,
                "z = {z}: a flat zero-volume piece with {} faces",
                m.n_faces()
            );
        }
    }
}

#[test]
fn split_through_a_cube_diagonal_edge_pair() {
    let a = fixtures::cube();
    let (neg, pos) = split_by_plane(&a, [0.0; 3], [1.0, -1.0, 0.0]).unwrap();
    assert!(ok(&neg) && ok(&pos));
    assert!(close(volume(&neg), 0.5, 1e-12) && close(volume(&pos), 0.5, 1e-12));
}

fn half() -> impl Strategy<Value = f64> {
    (-2i32..5).prop_map(|i| i as f64 * 0.5)
}

fn grid_box() -> impl Strategy<Value = ([f64; 3], [f64; 3])> {
    (
        proptest::array::uniform3(half()),
        proptest::array::uniform3(1i32..4),
    )
        .prop_map(|(lo, s)| (lo, [0, 1, 2].map(|k| lo[k] + s[k] as f64 * 0.5)))
}

/// Exact volume of a set of half-cells on the 0.5 grid over [-1, 4)^3.
fn voxels(f: impl Fn([f64; 3]) -> bool) -> f64 {
    let mut n = 0;
    for i in -2..8 {
        for j in -2..8 {
            for k in -2..8 {
                let c = [i, j, k].map(|x| x as f64 * 0.5 + 0.25);
                n += usize::from(f(c));
            }
        }
    }
    n as f64 * 0.125
}

fn inside(b: &([f64; 3], [f64; 3]), p: [f64; 3]) -> bool {
    (0..3).all(|k| b.0[k] < p[k] && p[k] < b.1[k])
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

    #[test]
    fn chained_grid_booleans_match_voxels(
        b0 in grid_box(), b1 in grid_box(), b2 in grid_box(), diag in any::<bool>(),
    ) {
        let m0 = fixtures::box_mesh(b0.0, b0.1);
        let m1 = if diag { other_diagonals(b1.0, b1.1) } else { fixtures::box_mesh(b1.0, b1.1) };
        let m2 = fixtures::box_mesh(b2.0, b2.1);
        let Ok(a) = boolean(&m0, &m1, Op::Union, &BooleanOptions::default()) else {
            return Ok(());
        };
        let a = a.mesh;
        prop_assert!(solid(&a), "union of two boxes is not a valid solid");
        let ina = |p| inside(&b0, p) || inside(&b1, p);
        for (op, f) in [
            (Op::Union, (|x: bool, y: bool| x || y) as fn(bool, bool) -> bool),
            (Op::Intersection, |x, y| x && y),
            (Op::Difference, |x, y| x && !y),
            (Op::SymmetricDifference, |x, y| x != y),
        ] {
            match boolean(&a, &m2, op, &BooleanOptions::default()) {
                Ok(r) => {
                    prop_assert!(ok(&r.mesh), "{:?} not closed manifold", op);
                    let want = voxels(|p| f(ina(p), inside(&b2, p)));
                    prop_assert!(close(volume(&r.mesh), want, 1e-12),
                        "{:?}: {} != {}", op, volume(&r.mesh), want);
                }
                Err(e) => prop_assert!(e.to_string().contains("not manifold"), "{:?}: {}", op, e),
            }
        }
    }

    #[test]
    fn identical_rotated_solids_with_other_diagonals(
        angles in proptest::array::uniform3(-3.2..3.2f64),
    ) {
        let rot = |m: Mesh| {
            let v = m.vertices().iter().map(|p| {
                let (s, c) = (angles.map(f64::sin), angles.map(f64::cos));
                let x = [p[0], c[0] * p[1] - s[0] * p[2], s[0] * p[1] + c[0] * p[2]];
                let y = [c[1] * x[0] + s[1] * x[2], x[1], -s[1] * x[0] + c[1] * x[2]];
                [c[2] * y[0] - s[2] * y[1], s[2] * y[0] + c[2] * y[1], y[2]]
            }).collect();
            m.with_vertices(v).unwrap()
        };
        let a = rot(fixtures::cube());
        let b = rot(other_diagonals([0.0; 3], [1.0; 3]));
        let va = volume(&a);
        for op in [Op::Union, Op::Intersection] {
            let m = boolean(&a, &b, op, &BooleanOptions::default())
                .map_err(|e| TestCaseError::fail(format!("{op:?}: {e}")))?.mesh;
            prop_assert!(ok(&m));
            prop_assert!(close(volume(&m), va, 1e-12), "{:?}: {} vs {}", op, volume(&m), va);
        }
        for op in [Op::Difference, Op::SymmetricDifference] {
            match boolean(&a, &b, op, &BooleanOptions::default()) {
                Ok(r) => {
                    prop_assert!(ok(&r.mesh));
                    prop_assert!(volume(&r.mesh).abs() < 1e-12, "{:?}: {}", op, volume(&r.mesh));
                }
                Err(e) => prop_assert!(e.to_string().contains("not manifold"), "{:?}: {}", op, e),
            }
        }
    }

    #[test]
    fn nested_cavities_and_islands_match_voxels(
        b2 in grid_box(), shift in proptest::array::uniform3(-1i32..2), island in any::<bool>(),
    ) {
        let d = shift.map(|s| s as f64 * 0.5);
        let mut parts = vec![fixtures::nested_shells()];
        if island {
            parts.push(fixtures::box_mesh([1.25; 3], [1.75; 3]));
        }
        let a = fixtures::translated(&concat(&parts), d);
        prop_assert!(solid(&a));
        let ina = |p: [f64; 3]| {
            let q = [0, 1, 2].map(|k| p[k] - d[k]);
            let c = |lo: f64, hi: f64| (0..3).all(|k| lo < q[k] && q[k] < hi);
            (c(0.0, 3.0) && !c(1.0, 2.0)) || (island && c(1.25, 1.75))
        };
        let fine = |f: &dyn Fn([f64; 3]) -> bool| {
            let mut n = 0;
            for i in -8..32 {
                for j in -8..32 {
                    for k in -8..32 {
                        n += usize::from(f([i, j, k].map(|x| x as f64 * 0.125 + 0.0625)));
                    }
                }
            }
            n as f64 / 512.0
        };
        let m2 = fixtures::box_mesh(b2.0, b2.1);
        for (op, f) in [
            (Op::Union, (|x: bool, y: bool| x || y) as fn(bool, bool) -> bool),
            (Op::Intersection, |x, y| x && y),
            (Op::Difference, |x, y| x && !y),
            (Op::SymmetricDifference, |x, y| x != y),
        ] {
            match boolean(&a, &m2, op, &BooleanOptions::default()) {
                Ok(r) => {
                    prop_assert!(ok(&r.mesh), "{:?} not closed manifold", op);
                    let want = fine(&|p| f(ina(p), inside(&b2, p)));
                    prop_assert!(close(volume(&r.mesh), want, 1e-12),
                        "{:?}: {} != {}", op, volume(&r.mesh), want);
                }
                Err(e) => prop_assert!(e.to_string().contains("not manifold"), "{:?}: {}", op, e),
            }
        }
    }
}

mod clipping_degeneracies {
    use super::*;
    use wireframes_boolean::{clip_by_surface, cut_fill, slice_by_plane};

    fn sheet(n: u32, lo: [f64; 2], hi: [f64; 2], z: impl Fn(f64, f64) -> f64) -> Mesh {
        let sx = (hi[0] - lo[0]) / (n - 1) as f64;
        let sy = (hi[1] - lo[1]) / (n - 1) as f64;
        let mut v = vec![];
        for j in 0..n {
            for i in 0..n {
                let (x, y) = (lo[0] + i as f64 * sx, lo[1] + j as f64 * sy);
                v.push([x, y, z(x, y)]);
            }
        }
        let id = |i: u32, j: u32| j * n + i;
        let mut f = vec![];
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        Mesh::new(v, f).unwrap()
    }

    fn clip(solid: &Mesh, s: &Mesh) -> (Mesh, Mesh) {
        let r = clip_by_surface(solid, s, &BooleanOptions::default())
            .unwrap_or_else(|e| panic!("clip failed: {e}"));
        (r.below, r.above)
    }

    #[test]
    fn clip_by_sheet_on_the_bottom_face() {
        let (below, above) = clip(
            &fixtures::cube(),
            &sheet(4, [-1.0; 2], [2.0; 2], |_, _| 0.0),
        );
        assert!(ok(&below) && ok(&above));
        assert!(below.is_empty(), "below has {} faces", below.n_faces());
        assert!(close(volume(&above), 1.0, 1e-12));
    }

    #[test]
    fn clip_by_sheet_folding_off_the_top_face() {
        let s = sheet(
            7,
            [-1.0; 2],
            [2.0; 2],
            |x, _| if x <= 0.5 { 1.0 } else { 1.5 - x },
        );
        let (below, above) = clip(&fixtures::cube(), &s);
        assert!(ok(&below) && ok(&above));
        assert!(
            close(volume(&above), 0.125, 1e-12),
            "above {}",
            volume(&above)
        );
        assert!(
            close(volume(&below), 0.875, 1e-12),
            "below {}",
            volume(&below)
        );
    }

    #[test]
    fn clip_by_sheet_through_two_cube_edges() {
        let s = sheet(7, [-1.0; 2], [2.0; 2], |x, _| x);
        let (below, above) = clip(&fixtures::cube(), &s);
        assert!(ok(&below) && ok(&above));
        assert!(close(volume(&above), 0.5, 1e-12) && close(volume(&below), 0.5, 1e-12));
    }

    #[test]
    fn clip_by_partial_sheet_lying_on_the_top_face() {
        let s = sheet(4, [-1.0, -1.0], [0.5, 2.0], |_, _| 1.0);
        let (below, above) = clip(&fixtures::cube(), &s);
        assert!(ok(&below) && ok(&above));
        assert!(close(volume(&below) + volume(&above), 1.0, 1e-12));
    }

    #[test]
    fn clip_by_sheet_touching_one_vertex() {
        let s = sheet(3, [0.0, 0.0], [2.0, 2.0], |x, y| -(x + y));
        let (below, above) = clip(&fixtures::cube(), &s);
        assert!(ok(&below) && ok(&above));
        assert!(below.is_empty() && close(volume(&above), 1.0, 1e-12));
    }

    #[test]
    fn utm_clip_keeps_volume() {
        let x0 = [712_345.37, 7_481_234.91, 812.13];
        let c = fixtures::box_mesh(x0, [x0[0] + 10.0, x0[1] + 10.0, x0[2] + 10.0]);
        let s = sheet(
            5,
            [x0[0] - 5.0, x0[1] - 5.0],
            [x0[0] + 15.0, x0[1] + 15.0],
            |x, y| x0[2] + 5.0 + 0.1 * (x - x0[0]) - 0.05 * (y - x0[1]),
        );
        let (below, above) = clip(&c, &s);
        assert!(ok(&below) && ok(&above));
        assert!(close(volume(&below) + volume(&above), 1000.0, 1e-9));
        assert!(
            close(volume(&below), 525.0, 1e-9),
            "below {}",
            volume(&below)
        );
    }

    #[test]
    fn slice_on_face_planes_does_not_fail() {
        let c = fixtures::cube();
        for z in [0.0, 1.0] {
            let s = slice_by_plane(&c, [0.0, 0.0, z], [0.0, 0.0, 1.0], &[]);
            assert!(s.is_ok(), "z = {z}: {}", s.err().unwrap());
        }
        let s = slice_by_plane(&c, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], &[]).unwrap();
        let (n, a) = (s.n_parts(), s.area());
        assert!(
            n == 1 && s.closed()[0] && close(a[0], 1.0, 1e-12),
            "{n} parts, areas {a:?}"
        );
        let s = slice_by_plane(&c, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], &[]).unwrap();
        assert_eq!(s.n_parts(), 0);
    }

    #[test]
    fn cut_fill_touching_along_an_edge_has_zero_area() {
        let base = sheet(3, [0.0, 0.0], [2.0, 2.0], |_, _| 0.0);
        let design = sheet(3, [2.0, 0.0], [4.0, 2.0], |_, _| 1.0);
        let r = cut_fill(&base, &design, None).unwrap();
        assert!(r.area == 0.0 && r.cut == 0.0 && r.fill == 0.0, "{r:?}");
    }

    #[test]
    fn utm_cut_fill_matches_a_constant_shift() {
        let (x, y) = (712_345.37, 7_481_234.91);
        let base = sheet(6, [x, y], [x + 100.0, y + 100.0], |a, b| {
            800.0 + 0.01 * (a - x) - 0.02 * (b - y)
        });
        let v: Vec<[f64; 3]> = base
            .vertices()
            .iter()
            .map(|p| [p[0], p[1], p[2] + 2.5])
            .collect();
        let up = base.clone().with_vertices(v).unwrap();
        let r = cut_fill(&base, &up, None).unwrap();
        let area: f64 = (x + 100.0 - x) * (y + 100.0 - y);
        assert!(close(r.area, area, 1e-9), "area {} vs {area}", r.area);
        assert!(
            close(r.fill, 2.5 * area, 1e-9),
            "fill {} vs {}",
            r.fill,
            2.5 * area
        );
        assert!(r.cut == 0.0);
    }
}

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

fn gpt(scale: f64, off: [f64; 3]) -> impl Strategy<Value = [f64; 3]> {
    proptest::array::uniform3(0i32..4)
        .prop_map(move |c| [0, 1, 2].map(|k| off[k] + c[k] as f64 * scale))
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn grid_snapped_tetrahedra_and_boxes(
        p in proptest::array::uniform4(gpt(0.5, [0.0; 3])),
        q in proptest::array::uniform4(gpt(0.5, [0.0; 3])),
        boxed in any::<bool>(),
        utm in any::<bool>(),
    ) {
        let off = if utm { [712_345.0, 7_481_234.0, 812.0] } else { [0.0; 3] };
        let sh = |v: [[f64; 3]; 4]| v.map(|x| [0, 1, 2].map(|k| x[k] + off[k]));
        let Some(a) = tet(sh(p)) else { return Ok(()) };
        let b = if boxed {
            fixtures::box_mesh(off.map(|o| o + 0.5), off.map(|o| o + 1.0))
        } else {
            let Some(b) = tet(sh(q)) else { return Ok(()) };
            b
        };
        let res: Vec<_> = [Op::Union, Op::Intersection, Op::Difference, Op::SymmetricDifference]
            .iter()
            .map(|&op| boolean(&a, &b, op, &BooleanOptions::default()))
            .collect();
        let mut vols = vec![];
        for (k, r) in res.iter().enumerate() {
            match r {
                Ok(r) => {
                    prop_assert!(ok(&r.mesh), "op {} not closed manifold", k);
                    vols.push(Some(volume(&r.mesh)));
                }
                Err(e) => {
                    prop_assert!(e.to_string().contains("not manifold"), "op {}: {}", k, e);
                    vols.push(None);
                }
            }
        }
        let (va, vb) = (volume(&a), volume(&b));
        if let (Some(vu), Some(vi)) = (vols[0], vols[1]) {
            prop_assert!(close(vu + vi, va + vb, 1e-9), "{} + {} vs {} + {}", vu, vi, va, vb);
        }
        if let (Some(vd), Some(vi)) = (vols[2], vols[1]) {
            prop_assert!(close(vd, va - vi, 1e-9), "{} vs {} - {}", vd, va, vi);
        }
        if let (Some(vx), Some(vu), Some(vi)) = (vols[3], vols[0], vols[1]) {
            prop_assert!(close(vx, vu - vi, 1e-9), "{} vs {} - {}", vx, vu, vi);
        }
    }
}

#[test]
fn split_on_cavity_face_planes() {
    let a = fixtures::nested_shells();
    let mut bad = vec![];
    for (z, neg_v) in [(1.0, 9.0), (2.0, 17.0), (1.5, 13.5 - 0.5)] {
        let (neg, pos) = split_by_plane(&a, [0.0, 0.0, z], [0.0, 0.0, 1.0]).unwrap();
        assert!(ok(&neg) && ok(&pos), "z = {z}: not closed manifold");
        assert!(
            close(volume(&neg), neg_v, 1e-12),
            "z = {z}: neg {}",
            volume(&neg)
        );
        assert!(
            close(volume(&pos), 26.0 - neg_v, 1e-12),
            "z = {z}: pos {}",
            volume(&pos)
        );
        for (n, m) in [("neg", &neg), ("pos", &pos)] {
            let s = m.validate(0.0).unwrap().summary;
            if !s.is_solid() {
                bad.push(format!(
                    "z = {z}: {n} has {} shells, {} self-intersections, {} coplanar overlaps",
                    s.shells, s.self_intersections, s.coplanar_overlaps
                ));
            }
        }
    }
    assert!(bad.is_empty(), "{bad:#?}");
}

fn slab_case() -> (Mesh, Mesh) {
    let mut parts = vec![];
    for i in 0..5 {
        for j in 0..5 {
            parts.push(fixtures::translated(
                &other_diagonals([0.0; 3], [1.0; 3]),
                [2.0 * i as f64, 2.0 * j as f64, 0.0],
            ));
        }
    }
    (
        concat(&parts),
        fixtures::box_mesh([0.5, 0.0, 0.5], [8.5, 9.0, 1.0]),
    )
}

#[test]
fn coplanar_heavy_results_do_not_depend_on_thread_count() {
    let (a, b) = slab_case();
    for op in [
        Op::Union,
        Op::Intersection,
        Op::Difference,
        Op::SymmetricDifference,
    ] {
        let go = |n: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| boolean(&a, &b, op, &BooleanOptions::default()).unwrap())
        };
        let (r1, r3) = (go(1), go(3));
        assert_eq!(r1.mesh, r3.mesh, "{op:?}");
        assert_eq!(r1.face_source, r3.face_source, "{op:?}");
    }
}

#[test]
fn face_source_planes_contain_their_output_faces() {
    use wireframes_core::predicates::{Sign, side_of_plane};
    let (a, b) = slab_case();
    for op in [
        Op::Union,
        Op::Intersection,
        Op::Difference,
        Op::SymmetricDifference,
    ] {
        let r = boolean(&a, &b, op, &BooleanOptions::default()).unwrap();
        assert_eq!(r.face_source.len(), r.mesh.n_faces());
        for (f, &(m, s)) in r.face_source.iter().enumerate() {
            let src = if m == 0 { &a } else { &b };
            let c = src.corners(s);
            for p in r.mesh.corners(f as u32) {
                assert_eq!(
                    side_of_plane(p, c[0], c[1], c[2]),
                    Sign::Zero,
                    "{op:?} face {f}"
                );
            }
        }
    }
}

#[test]
fn attributes_on_empty_and_coplanar_results() {
    use wireframes_core::{Column, ColumnData};
    let a = fixtures::cube()
        .with_vertex_column("vid", Column::i64((0..8).collect()))
        .unwrap()
        .with_face_column("rock", Column::i64(vec![1; 12]))
        .unwrap();
    let far = fixtures::translated(&fixtures::cube(), [5.0, 0.0, 0.0]);
    let e = run(&a, &far, Op::Intersection);
    assert!(e.is_empty());
    assert!(
        e.face_attributes().column("rock").is_some(),
        "empty result lost the face schema"
    );
    assert!(
        e.vertex_attributes().column("vid").is_some(),
        "empty result lost the vertex schema"
    );

    let b = fixtures::translated(&fixtures::cube(), [0.5, 0.5, 0.0]);
    let r = boolean(&a, &b, Op::Union, &BooleanOptions::default()).unwrap();
    let vid = r.mesh.vertex_attributes().column("vid").unwrap();
    let ColumnData::Int64(vals) = &vid.data else {
        panic!("vid type")
    };
    for (k, p) in r.mesh.vertices().iter().enumerate() {
        let from_a = a.vertices().iter().position(|q| q == p);
        match from_a {
            Some(i) => {
                assert!(
                    vid.is_valid(k) && vals[k] == i as i64,
                    "vertex {k} at {p:?}"
                )
            }
            None => assert!(
                !vid.is_valid(k),
                "vertex {k} at {p:?} is not from a but has vid"
            ),
        }
    }
    let rock = r.mesh.face_attributes().column("rock").unwrap();
    for (f, &(m, _)) in r.face_source.iter().enumerate() {
        assert_eq!(rock.is_valid(f), m == 0, "face {f}");
    }
}
