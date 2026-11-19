//! Solids made of closed shells that touch along edges or at vertices: valid
//! inputs, valid outputs, and chainable results.

mod support;
use support::closed_manifold;

use wireframes_boolean::{BooleanOptions, Op, boolean, clip_by_surface};
use wireframes_core::fixtures::{self, Rng};
use wireframes_core::{Mesh, MeshSummary};

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

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * 1f64.max(a.abs()).max(b.abs())
}

fn summary(m: &Mesh) -> MeshSummary {
    m.validate(0.0).unwrap().summary
}

/// Runs `op`, and checks the result is empty or a valid closed solid.
fn run(a: &Mesh, b: &Mesh, op: Op) -> Mesh {
    let m = boolean(a, b, op, &BooleanOptions::default())
        .unwrap_or_else(|e| panic!("{op:?} failed: {e}"))
        .mesh;
    let s = summary(&m);
    assert!(
        m.is_empty() || (closed_manifold(&m) && s.is_solid()),
        "{op:?}: {s:?}"
    );
    m
}

const OPS: [Op; 4] = [
    Op::Union,
    Op::Intersection,
    Op::Difference,
    Op::SymmetricDifference,
];

type Inside = Box<dyn Fn([f64; 3]) -> bool>;

fn in_box(lo: [f64; 3], hi: [f64; 3]) -> Inside {
    Box::new(move |p| (0..3).all(|k| lo[k] < p[k] && p[k] < hi[k]))
}

fn combine(op: Op, a: Inside, b: Inside) -> Inside {
    Box::new(move |p| {
        let (x, y) = (a(p), b(p));
        match op {
            Op::Union => x || y,
            Op::Intersection => x && y,
            Op::Difference => x && !y,
            Op::SymmetricDifference => x != y,
        }
    })
}

/// Exact volume of a set bounded by planes on the 0.25 grid, from the
/// centres of the cells of that grid over `[-2, 6]^3`.
fn voxel_volume(inside: &dyn Fn([f64; 3]) -> bool) -> f64 {
    let mut n = 0usize;
    for i in 0..32 {
        for j in 0..32 {
            for k in 0..32 {
                n += usize::from(inside([i, j, k].map(|c| -2.0 + 0.25 * c as f64 + 0.125)));
            }
        }
    }
    n as f64 / 64.0
}

fn check_identities(a: &Mesh, b: &Mesh) -> [Mesh; 4] {
    let [u, i, d, x] = OPS.map(|op| run(a, b, op));
    let (va, vb) = (volume(a), volume(b));
    let [vu, vi, vd, vx] = [&u, &i, &d, &x].map(volume);
    assert!(
        close(vx, va + vb - 2.0 * vi),
        "xor {vx} vs {va} + {vb} - 2 {vi}"
    );
    assert!(
        close(vu + vi, va + vb),
        "union {vu} + inter {vi} vs {va} + {vb}"
    );
    assert!(close(vd, va - vi), "diff {vd} vs {va} - {vi}");
    [u, i, d, x]
}

#[test]
fn edge_contact_pair() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [1.0, 1.0, 0.0]);
    let [u, i, d, x] = check_identities(&a, &b);
    assert!(i.is_empty() && close(volume(&d), 1.0));
    for m in [&u, &x] {
        let s = summary(m);
        assert_eq!(
            (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
            (2, 1, 2)
        );
        assert!(close(volume(m), 2.0));
    }
    let welded = fixtures::edge_touching_cubes(true);
    let c = fixtures::box_mesh([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]);
    let [u, i, ..] = check_identities(&welded, &c);
    assert_eq!(summary(&u).shells, 1);
    assert!(close(volume(&i), 0.5));
}

#[test]
fn vertex_contact_pair() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [1.0, 1.0, 1.0]);
    let [u, i, _, x] = check_identities(&a, &b);
    assert!(i.is_empty());
    for m in [&u, &x] {
        let s = summary(m);
        assert_eq!(
            (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
            (2, 0, 1)
        );
    }
    let t = fixtures::touching_tetrahedra();
    let [u, ..] = check_identities(&t, &fixtures::box_mesh([0.1; 3], [2.0; 3]));
    assert_eq!(summary(&u).shells, 2);
}

#[test]
fn face_contact_pair_merges() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [1.0, 0.0, 0.0]);
    let [u, i, d, x] = check_identities(&a, &b);
    assert!(i.is_empty() && close(volume(&d), 1.0));
    for m in [&u, &x] {
        let s = summary(m);
        assert_eq!((s.shells, s.touching_shell_edges), (1, 0), "{s:?}");
        assert!(close(volume(m), 2.0));
    }
}

#[test]
fn difference_leaving_bodies_joined_at_a_vertex_and_an_edge() {
    let a = fixtures::box_mesh([0.0; 3], [2.0; 3]);
    let checker = Mesh::concat(&[
        &fixtures::box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 2.0]),
        &fixtures::box_mesh([0.0, 1.0, 0.0], [1.0, 2.0, 2.0]),
    ])
    .unwrap();
    assert!(summary(&checker).is_solid());
    let d = run(&a, &checker, Op::Difference);
    let s = summary(&d);
    assert_eq!((s.shells, s.touching_shell_edges), (2, 1), "{s:?}");
    assert!(close(volume(&d), 4.0));
    let corner = Mesh::concat(&[
        &fixtures::box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
        &fixtures::box_mesh([0.0, 1.0, 0.0], [1.0, 2.0, 1.0]),
        &fixtures::box_mesh([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]),
        &fixtures::box_mesh([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
    ])
    .unwrap();
    let d = run(&a, &corner, Op::Difference);
    let s = summary(&d);
    assert_eq!(s.shells, 4, "{s:?}");
    assert!(
        s.touching_shell_edges > 0 && s.touching_shell_vertices > 0,
        "{s:?}"
    );
    let bottom = run(
        &d,
        &fixtures::box_mesh([0.0; 3], [2.0, 2.0, 1.0]),
        Op::Intersection,
    );
    let s = summary(&bottom);
    assert_eq!((s.shells, s.touching_shell_edges), (2, 1), "{s:?}");
    assert!(close(volume(&bottom), 2.0));
}

#[test]
fn seeded_random_pairs_keep_the_volume_identities() {
    let mut rng = Rng::new(20261119);
    let grid = |rng: &mut Rng| -> ([f64; 3], [f64; 3]) {
        let lo = [0, 1, 2].map(|_| rng.below(5) as f64 * 0.5);
        let hi = lo.map(|x| x + (1 + rng.below(4)) as f64 * 0.5);
        (lo, hi)
    };
    let (mut touching, mut nested, mut self_touch) = (0, 0, 0);
    for case in 0..160 {
        let (alo, ahi) = grid(&mut rng);
        let a = fixtures::box_mesh(alo, ahi);
        let mut ia = in_box(alo, ahi);
        let a = match case % 4 {
            0 => a,
            1 => {
                let d = [ahi[0] - alo[0], ahi[1] - alo[1], 0.0];
                let lo2 = [ahi[0], ahi[1], alo[2]];
                let hi2 = [ahi[0] + d[0], ahi[1] + d[1], ahi[2]];
                let prev = ia;
                ia = combine(Op::Union, prev, in_box(lo2, hi2));
                Mesh::concat(&[&a, &fixtures::box_mesh(lo2, hi2)]).unwrap()
            }
            2 => {
                let lo2 = ahi;
                let hi2 = ahi.map(|x| x + 0.5);
                let prev = ia;
                ia = combine(Op::Union, prev, in_box(lo2, hi2));
                Mesh::concat(&[&a, &fixtures::box_mesh(lo2, hi2)]).unwrap()
            }
            _ => {
                let lo2 = [0, 1, 2].map(|k| alo[k] - 0.5);
                let hi2 = [0, 1, 2].map(|k| ahi[k] + 0.5);
                let prev = ia;
                ia = combine(Op::Union, prev, in_box(lo2, hi2));
                fixtures::box_mesh(lo2, hi2)
            }
        };
        assert!(summary(&a).is_solid(), "case {case}");
        let (blo, bhi) = grid(&mut rng);
        let (blo, bhi) = if case % 5 == 0 {
            let lo = alo.map(|x| x + 0.25);
            let hi = [0, 1, 2].map(|k| (ahi[k] - 0.25).max(lo[k] + 0.25));
            (lo, hi)
        } else {
            (blo, bhi)
        };
        let b = fixtures::box_mesh(blo, bhi);
        let ib = in_box(blo, bhi);
        let res: Vec<_> = OPS
            .iter()
            .map(|&op| boolean(&a, &b, op, &BooleanOptions::default()))
            .collect();
        if let Some(Err(e)) = res.iter().find(|r| r.is_err()) {
            assert!(e.to_string().contains("touches itself"), "case {case}: {e}");
            self_touch += 1;
            continue;
        }
        let m: Vec<Mesh> = res.into_iter().map(|r| r.unwrap().mesh).collect();
        for (k, r) in m.iter().enumerate() {
            let s = summary(r);
            assert!(r.is_empty() || s.is_solid(), "case {case} op {k}: {s:?}");
            touching += usize::from(s.touching_shell_edges + s.touching_shell_vertices > 0);
            nested += usize::from(s.shells > 1 && s.touching_shell_vertices == 0);
        }
        let [vu, vi, vd, vx] = [0, 1, 2, 3].map(|k| volume(&m[k]));
        let (va, vb) = (volume(&a), volume(&b));
        assert!(close(vx, va + vb - 2.0 * vi), "case {case}: xor {vx}");
        assert!(close(vu + vi, va + vb), "case {case}: {vu} + {vi}");
        assert!(close(vd, va - vi), "case {case}: diff {vd}");
        let ia: &dyn Fn([f64; 3]) -> bool = &*ia;
        assert!(close(vu, voxel_volume(&|p| ia(p) || ib(p))), "case {case}");
        assert!(close(vx, voxel_volume(&|p| ia(p) != ib(p))), "case {case}");
    }
    println!(
        "{touching} results with touching shells, {nested} with separate shells, \
         {self_touch} of 160 cases refused as a shell touching itself"
    );
    assert!(touching > 0 && nested > 0);
    assert!(self_touch < 16, "{self_touch} self-touching cases");
}

#[test]
fn ten_mixed_booleans_chain_on_their_results() {
    let steps: [(Op, [f64; 3], [f64; 3]); 10] = [
        (Op::SymmetricDifference, [1.0; 3], [3.0; 3]),
        (Op::Union, [3.0, 3.0, 0.0], [4.0, 4.0, 1.0]),
        (Op::Difference, [-1.0, -1.0, -1.0], [0.5, 0.5, 3.0]),
        (Op::Union, [4.0, 2.0, 0.0], [5.0, 3.0, 1.0]),
        (Op::Intersection, [-1.0, -1.0, 0.0], [5.0, 5.0, 2.5]),
        (Op::SymmetricDifference, [1.5, 1.5, 0.5], [2.5, 2.5, 1.5]),
        (Op::Difference, [0.0, 0.0, 2.0], [5.0, 5.0, 2.25]),
        (Op::Union, [-1.0, -1.0, 0.0], [1.0, 1.0, 1.0]),
        (Op::SymmetricDifference, [2.0, 0.0, 0.0], [4.0, 1.0, 1.0]),
        (Op::Difference, [3.5, 2.5, -1.0], [4.5, 3.5, 3.0]),
    ];
    let mut m = fixtures::box_mesh([0.0; 3], [2.0; 3]);
    let mut inside = in_box([0.0; 3], [2.0; 3]);
    let mut contacts = 0;
    for (k, (op, lo, hi)) in steps.into_iter().enumerate() {
        let b = fixtures::box_mesh(lo, hi);
        if let Err(e) = boolean(&m, &b, op, &BooleanOptions::default()) {
            panic!("step {k} {op:?}: {e}");
        }
        let next = run(&m, &b, op);
        inside = combine(op, inside, in_box(lo, hi));
        let want = voxel_volume(&*inside);
        assert!(
            close(volume(&next), want),
            "step {k}: {} vs {want}",
            volume(&next)
        );
        let s = summary(&next);
        contacts += usize::from(s.touching_shell_edges + s.touching_shell_vertices > 0);
        m = next;
    }
    assert!(contacts >= 3, "{contacts} steps gave touching shells");
    let probes = fixtures::random_points(5, 400, [-1.5; 3], [5.5; 3]);
    let index = wireframes_core::topology::locate_point;
    for p in probes {
        let l = index(m.vertices(), m.triangles(), p);
        if l != 0 {
            assert_eq!(l == -1, inside(p), "at {p:?}");
        }
    }
}

#[test]
fn clip_keeps_touching_shells_apart() {
    let solid = fixtures::edge_touching_cubes(false);
    let sheet = Mesh::new(
        vec![
            [-1.0, -1.0, 0.5],
            [3.0, -1.0, 0.5],
            [3.0, 3.0, 0.5],
            [-1.0, 3.0, 0.5],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    let r = clip_by_surface(&solid, &sheet, &BooleanOptions::default()).unwrap();
    for m in [&r.below, &r.above] {
        let s = summary(m);
        assert!(s.is_solid(), "{s:?}");
        assert_eq!((s.shells, s.touching_shell_edges), (2, 1), "{s:?}");
        assert!(close(volume(m), 1.0));
    }
}

#[test]
fn a_shell_touching_itself_is_refused() {
    let l = run(
        &fixtures::box_mesh([0.0; 3], [2.0; 3]),
        &fixtures::box_mesh([1.0; 3], [3.0; 3]),
        Op::Difference,
    );
    let column = fixtures::box_mesh([-1.0; 3], [1.0, 1.0, 3.0]);
    let e = boolean(&l, &column, Op::Difference, &BooleanOptions::default()).unwrap_err();
    assert!(
        e.to_string().contains("touches itself along an edge"),
        "{e}"
    );
}

#[test]
fn repair_option_keeps_touching_inputs_apart() {
    let a = fixtures::edge_touching_cubes(false);
    let b = fixtures::box_mesh([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]);
    let plain = boolean(&a, &b, Op::Union, &BooleanOptions::default()).unwrap();
    let fixed = boolean(&a, &b, Op::Union, &BooleanOptions { repair: true }).unwrap();
    let report = fixed.repairs[0].as_ref().unwrap();
    assert!(report.actions.is_empty(), "{:?}", report.actions);
    assert_eq!(plain.mesh, fixed.mesh);
    let welded = fixtures::edge_touching_cubes(true);
    let r = boolean(&welded, &b, Op::Union, &BooleanOptions { repair: true }).unwrap();
    assert_eq!(r.repairs[0].as_ref().unwrap().actions.len(), 2);
    assert!(close(volume(&r.mesh), 2.0 + 1.0 - 0.5));
}
