//! Solids made of closed shells that touch along edges or at vertices: valid
//! inputs, valid outputs, and chainable results.

mod support;
use support::closed_manifold;

use wireframes_boolean::{
    BooleanOptions, Op, boolean, clip_by_surface, slice_by_plane, split_by_plane,
};
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
fn a_shell_touching_itself_is_a_solid() {
    let l = run(
        &fixtures::box_mesh([0.0; 3], [2.0; 3]),
        &fixtures::box_mesh([1.0; 3], [3.0; 3]),
        Op::Difference,
    );
    let column = fixtures::box_mesh([-1.0; 3], [1.0, 1.0, 3.0]);
    let [_, _, d, _] = check_identities(&l, &column);
    let s = summary(&d);
    assert_eq!((s.shells, s.touching_shell_edges), (1, 1), "{s:?}");
    assert!(close(volume(&d), 5.0));
}

#[test]
fn closing_a_loop_of_cells_in_any_order_gives_the_same_solid() {
    let mut cells = fixtures::EDGE_LOOP.to_vec();
    cells.push([1, 0, 0]);
    let cell = |c: [i32; 3]| fixtures::voxels(&[c]);
    for gap_last in [true, false] {
        let mut order = cells.clone();
        if !gap_last {
            order.swap(4, 5);
        }
        let mut u = cell(order[0]);
        for (k, &c) in order.iter().enumerate().skip(1) {
            u = run(&u, &cell(c), Op::Union);
            assert!(close(volume(&u), (k + 1) as f64), "{gap_last} {k}");
        }
        let s = summary(&u);
        assert_eq!((s.shells, s.touching_shell_edges), (1, 0), "{s:?}");
    }
    let open = fixtures::voxels(&fixtures::EDGE_LOOP);
    let s = summary(&run(&open, &fixtures::cube(), Op::Intersection));
    assert_eq!((s.shells, s.touching_shell_edges), (1, 0));
}

#[test]
fn a_shell_touching_itself_along_an_edge_with_single_umbrella_ends() {
    let mut u = fixtures::voxels(&fixtures::EDGE_LOOP);
    for c in [[0, 0, -1], [1, 0, -1], [1, 1, -1]] {
        u = run(&u, &fixtures::voxels(&[c]), Op::Union);
    }
    assert!(close(volume(&u), 8.0));
    let s = summary(&u);
    assert_eq!(
        (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
        (1, 1, 2),
        "{s:?}"
    );
    let filled = run(&u, &fixtures::voxels(&[[1, 0, 0]]), Op::Union);
    assert_eq!(summary(&filled).touching_shell_edges, 0);
    assert!(close(volume(&filled), 9.0));
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

#[test]
fn slicing_and_splitting_touching_shells() {
    for (name, m) in [
        ("apart", fixtures::edge_touching_cubes(false)),
        ("welded", fixtures::edge_touching_cubes(true)),
    ] {
        let s = slice_by_plane(&m, [0.0, 0.0, 0.5], [0.0, 0.0, 1.0], &[]).unwrap();
        assert_eq!(s.n_parts(), 2, "{name}");
        assert!(s.closed().iter().all(|c| *c), "{name}");
        let areas: f64 = s.area().iter().sum();
        assert!(close(areas, 2.0), "{name}: {areas}");
        let (below, above) = split_by_plane(&m, [0.0, 0.0, 0.5], [0.0, 0.0, 1.0]).unwrap();
        for piece in [&below, &above] {
            let r = summary(piece);
            assert!(r.is_solid(), "{name}: {r:?}");
            assert_eq!(r.shells, 2, "{name}");
            assert!(close(volume(piece), 1.0), "{name}");
        }
        let (left, right) = split_by_plane(&m, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        assert!(
            close(volume(&left), 1.0) && close(volume(&right), 1.0),
            "{name}"
        );
        let diagonal = [1.0, 1.0, 0.0];
        let (a, b) = split_by_plane(&m, [1.0, 1.0, 0.0], diagonal).unwrap();
        assert!(close(volume(&a), 1.0) && close(volume(&b), 1.0), "{name}");
        assert!(summary(&a).is_solid() && summary(&b).is_solid(), "{name}");
    }
}

/// Triangular prism over `z` in `[0, 1]` with one vertical edge on the `z`
/// axis and the others above `a` and `b`, wound outward.
fn wedge(a: [f64; 2], b: [f64; 2]) -> Mesh {
    let v = [[0.0, 0.0], a, b]
        .iter()
        .flat_map(|p| [[p[0], p[1], 0.0], [p[0], p[1], 1.0]])
        .collect();
    let t = vec![
        [0, 4, 2],
        [1, 3, 5],
        [0, 2, 3],
        [0, 3, 1],
        [2, 4, 5],
        [2, 5, 3],
        [4, 0, 1],
        [4, 1, 5],
    ];
    let m = Mesh::new(v, t).unwrap();
    if volume(&m) < 0.0 {
        fixtures::flipped(&m)
    } else {
        m
    }
}

#[test]
fn many_shells_on_one_edge_and_one_vertex() {
    let d = [
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 0.0],
        [-1.0, -1.0],
        [0.0, -1.0],
        [1.0, -1.0],
    ];
    let wedges: Vec<Mesh> = (0..4).map(|k| wedge(d[2 * k], d[2 * k + 1])).collect();
    let fan = Mesh::concat(&wedges.iter().collect::<Vec<_>>()).unwrap();
    let s = summary(&fan);
    assert!(s.is_solid(), "{s:?}");
    assert_eq!(
        (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
        (4, 1, 2)
    );
    let core = fixtures::box_mesh([-0.5, -0.5, 0.25], [0.5, 0.5, 0.75]);
    let i = run(&fan, &core, Op::Intersection);
    assert_eq!(summary(&i).touching_shell_edges, 1);
    assert!(close(volume(&i), 0.25));
    let d = run(&fan, &core, Op::Difference);
    assert!(close(volume(&d), 2.0 - 0.25));
    let u = run(&fan, &core, Op::Union);
    assert!(close(volume(&u), 2.0 + 0.25));
    let s = summary(&u);
    assert_eq!((s.shells, s.touching_shell_edges), (1, 2), "{s:?}");

    let octants = Mesh::concat(&[
        &fixtures::box_mesh([0.0; 3], [1.0; 3]),
        &fixtures::box_mesh([-1.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
        &fixtures::box_mesh([-1.0, 0.0, -1.0], [0.0, 1.0, 0.0]),
        &fixtures::box_mesh([0.0, -1.0, -1.0], [1.0, 0.0, 0.0]),
    ])
    .unwrap();
    let s = summary(&octants);
    assert!(s.is_solid(), "{s:?}");
    assert_eq!(
        (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
        (4, 6, 7)
    );
    let ball = fixtures::box_mesh([-0.5; 3], [0.5; 3]);
    let [i, d, x] = [Op::Intersection, Op::Difference, Op::SymmetricDifference]
        .map(|op| run(&octants, &ball, op));
    assert!(close(volume(&i), 0.5) && close(volume(&d), 3.5) && close(volume(&x), 4.0));
    assert_eq!(summary(&x).shells, 8);
}

#[test]
fn a_cavity_touching_the_outer_shell_at_its_vertices() {
    let m = fixtures::vertex_touching_cavity();
    let want = 4.5 - 1.0 / 6.0;
    assert!(close(volume(&m), want));
    let part = |outer: bool| {
        let keep: Vec<bool> = (0..m.n_faces()).map(|f| (f < 12) == outer).collect();
        m.filter_faces(&keep).unwrap()
    };
    let (outer, inner) = (part(true), fixtures::flipped(&part(false)));
    let d = run(&outer, &inner, Op::Difference);
    let s = summary(&d);
    assert_eq!(
        (s.shells, s.inward_shells, s.touching_shell_vertices),
        (2, 0, 4),
        "{s:?}"
    );
    assert!(close(volume(&d), want));
    let side = fixtures::box_mesh([1.5, -1.0, -1.0], [4.0, 4.0, 4.0]);
    check_identities(&m, &side);
    let cut = fixtures::box_mesh([-1.0, -1.0, 0.5], [4.0, 4.0, 4.0]);
    let [_, i, _, _] = check_identities(&m, &cut);
    assert!(summary(&i).touching_shell_vertices > 0);
    let u = run(&m, &inner, Op::Union);
    assert!(close(volume(&u), 4.5));

    let big = fixtures::box_mesh([0.0; 3], [4.0; 3]);
    let tip = fixtures::tetrahedron();
    let tip = tip
        .clone()
        .with_vertices([[0.0; 3], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]].to_vec())
        .unwrap();
    let tip = if volume(&tip) < 0.0 {
        fixtures::flipped(&tip)
    } else {
        tip
    };
    let d = run(&big, &tip, Op::Difference);
    let s = summary(&d);
    assert_eq!((s.shells, s.touching_shell_vertices), (2, 1), "{s:?}");
    assert!(close(volume(&d), 64.0 - volume(&tip)));
}

#[test]
fn split_through_a_lone_vertex_on_the_plane() {
    let (a, b) = split_by_plane(&fixtures::cube(), [1.0; 3], [1.0; 3]).unwrap();
    assert_eq!((a.n_faces(), b.n_faces()), (12, 0));
    let (a, b) = split_by_plane(&fixtures::vertex_touching_cubes(), [1.0; 3], [1.0; 3]).unwrap();
    assert!(close(volume(&a), 1.0) && close(volume(&b), 1.0));
    for m in [
        fixtures::edge_touching_cubes(false),
        fixtures::edge_touching_cubes(true),
        fixtures::touching_tetrahedra(),
    ] {
        for (o, n) in [
            ([1.0; 3], [1.0; 3]),
            ([0.0; 3], [1.0; 3]),
            ([1.0, 1.0, 1.0], [1.0, -1.0, 0.0]),
        ] {
            let (a, b) = split_by_plane(&m, o, n).unwrap();
            assert!(close(volume(&a) + volume(&b), volume(&m)), "{o:?} {n:?}");
            for piece in [&a, &b] {
                assert!(piece.is_empty() || summary(piece).is_solid(), "{o:?} {n:?}");
            }
        }
    }
}

fn permuted(m: &Mesh, seed: u64) -> Mesh {
    let mut rng = Rng::new(seed);
    let mut order: Vec<u32> = (0..m.n_vertices() as u32).collect();
    rng.shuffle(&mut order);
    let mut slot = vec![0u32; order.len()];
    for (k, &o) in order.iter().enumerate() {
        slot[o as usize] = k as u32;
    }
    let v = order.iter().map(|&o| m.vertices()[o as usize]).collect();
    let mut t: Vec<[u32; 3]> = m
        .triangles()
        .iter()
        .map(|t| {
            let t = t.map(|i| slot[i as usize]);
            let r = rng.below(3);
            [t[r], t[(r + 1) % 3], t[(r + 2) % 3]]
        })
        .collect();
    rng.shuffle(&mut t);
    Mesh::new(v, t).unwrap()
}

/// Faces as coordinate triples starting at the smallest corner, sorted.
fn canonical(m: &Mesh) -> Vec<[[u64; 3]; 3]> {
    let mut out: Vec<[[u64; 3]; 3]> = (0..m.n_faces() as u32)
        .map(|f| {
            let c = m.corners(f).map(|p| p.map(|x| (x + 0.0).to_bits()));
            let r = (0..3).min_by_key(|&i| c[i]).unwrap();
            [c[r], c[(r + 1) % 3], c[(r + 2) % 3]]
        })
        .collect();
    out.sort_unstable();
    out
}

#[test]
fn permuted_inputs_and_thread_counts_give_the_same_results() {
    let cases = [
        (
            fixtures::edge_touching_cubes(false),
            fixtures::box_mesh([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]),
        ),
        (
            fixtures::edge_touching_cubes(true),
            fixtures::box_mesh([0.5, 0.5, -0.5], [1.5, 1.5, 1.5]),
        ),
        (
            fixtures::vertex_touching_cubes(),
            fixtures::box_mesh([0.5; 3], [1.5; 3]),
        ),
        (
            fixtures::vertex_touching_cavity(),
            fixtures::box_mesh([1.5, -1.0, -1.0], [4.0, 4.0, 4.0]),
        ),
        (fixtures::cube(), fixtures::box_mesh([0.5; 3], [1.5; 3])),
    ];
    let opts = BooleanOptions::default();
    let outcome = |a: &Mesh, b: &Mesh, op: Op| {
        boolean(a, b, op, &opts)
            .map(|r| (canonical(&r.mesh), summary(&r.mesh)))
            .map_err(|e| e.to_string())
    };
    let mut refused = 0;
    for (n, (a, b)) in cases.iter().enumerate() {
        for op in OPS {
            let base = outcome(a, b, op);
            match &base {
                Ok((_, s)) => assert!(s.is_solid() || s == &MeshSummary::default(), "{s:?}"),
                Err(e) => {
                    assert!(e.contains("touches itself"), "case {n} {op:?}: {e}");
                    refused += 1;
                }
            }
            for seed in 0..4 {
                let (pa, pb) = (permuted(a, seed), permuted(b, seed + 100));
                assert_eq!(outcome(&pa, &pb, op), base, "case {n} {op:?} seed {seed}");
            }
            let on = |threads: usize| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| boolean(a, b, op, &opts).map(|r| (r.mesh, r.face_source)))
                    .map_err(|e| e.to_string())
            };
            assert_eq!(on(1), on(8), "case {n} {op:?}");
        }
    }
    assert!(refused <= 2, "{refused} refusals");
}

#[test]
fn contacts_at_utm_scale_with_micrometre_offsets() {
    let o = [700000.0, 7500000.0, 100.0];
    let a = fixtures::translated(&fixtures::cube(), o);
    let at = |d: [f64; 3]| fixtures::translated(&fixtures::cube(), [0, 1, 2].map(|k| o[k] + d[k]));
    for (d, edges, points) in [([1.0, 1.0, 0.0], 1, 2), ([1.0, 1.0, 1.0], 0, 1)] {
        let b = at(d);
        let s = summary(&Mesh::concat(&[&a, &b]).unwrap());
        assert!(s.is_solid(), "{s:?}");
        assert_eq!(
            (s.touching_shell_edges, s.touching_shell_vertices),
            (edges, points)
        );
        let [u, ..] = check_identities(&a, &b);
        let s = summary(&u);
        assert_eq!((s.shells, s.touching_shell_edges), (2, edges));
    }
    let gap = at([1.0 + 1e-6, 1.0, 0.0]);
    let s = summary(&Mesh::concat(&[&a, &gap]).unwrap());
    assert!(s.is_solid() && s.touching_shell_vertices == 0, "{s:?}");
    check_identities(&a, &gap);
    let partial = at([1.0, 1.0, 1.0 - 1e-6]);
    let s = summary(&Mesh::concat(&[&a, &partial]).unwrap());
    assert!(!s.is_solid() && s.self_intersections > 0, "{s:?}");
    let [u, ..] = check_identities(&a, &partial);
    let s = summary(&u);
    assert_eq!(
        (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
        (2, 1, 2)
    );
}

#[test]
fn shells_touching_along_a_face_are_a_defect_and_merge_in_a_union() {
    let a = fixtures::cube();
    let b = fixtures::translated(&a, [1.0, 0.0, 0.0]);
    let s = summary(&Mesh::concat(&[&a, &b]).unwrap());
    assert!(!s.is_solid() && s.coplanar_overlaps > 0, "{s:?}");
    let e = boolean(
        &Mesh::concat(&[&a, &b]).unwrap(),
        &a,
        Op::Union,
        &BooleanOptions::default(),
    )
    .unwrap_err();
    assert!(e.to_string().contains("invalid mesh"), "{e}");
}
