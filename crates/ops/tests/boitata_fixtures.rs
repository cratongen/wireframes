//! boitata's mesh validation and repair fixtures (`crates/core/src/mesh.rs` tests),
//! ported with boitata's expected counts. Kinds boitata does not have
//! (`CoplanarOverlap`, `UnusedVertex`) are asserted separately.

use proptest::prelude::*;
use wireframes_core::{Mesh, MeshProblem, MeshProblemKind, MeshSummary, validate};
use wireframes_ops::{HoleFill, RepairOptions, RepairReport, repair, volume};

use MeshProblemKind::*;

fn mesh(vertices: &[[f64; 3]], triangles: &[[u32; 3]]) -> Mesh {
    Mesh::new(vertices.to_vec(), triangles.to_vec()).unwrap()
}

fn cube() -> Mesh {
    let v: Vec<[f64; 3]> = (0..8)
        .map(|i| [(i & 1) as f64, (i >> 1 & 1) as f64, (i >> 2) as f64])
        .collect();
    mesh(
        &v,
        &[
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
        ],
    )
}

fn two_cubes(shift: [f64; 3]) -> Mesh {
    let c = cube();
    let mut v = c.vertices().to_vec();
    v.extend(
        c.vertices()
            .iter()
            .map(|p| [0, 1, 2].map(|k| p[k] + shift[k])),
    );
    let mut t = c.triangles().to_vec();
    t.extend(c.triangles().iter().map(|f| f.map(|i| i + 8)));
    mesh(&v, &t)
}

fn flipped(ts: &[[u32; 3]]) -> Vec<[u32; 3]> {
    ts.iter().map(|&[a, b, c]| [a, c, b]).collect()
}

fn summary(m: &Mesh, tolerance: f64) -> MeshSummary {
    validate(m, tolerance).unwrap().summary
}

fn boitata_kind(k: MeshProblemKind) -> bool {
    !matches!(k, CoplanarOverlap | UnusedVertex)
}

fn kinds(m: &Mesh, tolerance: f64) -> Vec<(MeshProblemKind, Option<u32>, Option<u32>)> {
    validate(m, tolerance)
        .unwrap()
        .problems
        .iter()
        .filter(|p| boitata_kind(p.kind))
        .map(|p| (p.kind, p.face, p.vertex))
        .collect()
}

fn fix(m: &Mesh, weld: f64) -> (Mesh, RepairReport) {
    let opts = RepairOptions {
        weld: Some(weld),
        ..RepairOptions::default()
    };
    repair(m, &opts).unwrap()
}

fn unit_volume(m: &Mesh, tol: f64) -> bool {
    (volume(m).unwrap() - 1.0).abs() < tol
}

#[test]
fn clean_cube_has_no_problems_and_one_shell() {
    let r = validate(&cube(), 1e-9).unwrap();
    assert!(r.problems.is_empty(), "{:?}", r.problems);
    assert_eq!(r.summary.shells, 1);
    assert!(r.summary.is_closed);
}

#[test]
fn one_flipped_triangle_gives_three_inconsistent_edges() {
    let c = cube();
    for (face, wound) in [(2, [4, 6, 5]), (0, [0, 1, 2]), (7, [6, 3, 7])] {
        let mut t = c.triangles().to_vec();
        t[face] = wound;
        let m = mesh(c.vertices(), &t);
        let s = summary(&m, 0.0);
        assert!(s.is_closed);
        assert_eq!((s.inconsistent_edges, s.inward_shells), (3, 0), "{face}");
        let (fixed, _) = fix(&m, 0.0);
        assert!(validate(&fixed, 0.0).unwrap().problems.is_empty());
        assert!(unit_volume(&fixed, 1e-12));
    }
}

#[test]
fn open_square_has_four_boundary_edges() {
    let m = mesh(
        &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
        &[[0, 1, 2], [0, 2, 3]],
    );
    let s = summary(&m, 0.0);
    assert_eq!((s.boundary_edges, s.non_manifold_edges), (4, 0));
    assert!(!s.is_closed);
}

#[test]
fn fin_is_one_non_manifold_edge() {
    let m = mesh(
        &[
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., -1., 0.],
            [0., 0., 1.],
        ],
        &[[0, 1, 2], [0, 1, 3], [0, 1, 4]],
    );
    let s = summary(&m, 0.0);
    assert_eq!(s.non_manifold_edges, 1);
    assert!(!s.is_closed);
}

#[test]
fn repeated_and_collinear_triangles_are_degenerate() {
    let m = mesh(
        &[[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]],
        &[[0, 1, 1], [0, 1, 2]],
    );
    let s = summary(&m, 0.0);
    assert_eq!((s.degenerate_faces, s.boundary_edges), (2, 0));
    assert!(!s.is_closed);
}

#[test]
fn empty_mesh_is_open_and_bad_input_is_rejected() {
    assert!(!summary(&mesh(&[], &[]), 0.0).is_closed);
    assert!(Mesh::new(vec![[0.0; 3]], vec![[0, 0, 9]]).is_err());
    assert!(Mesh::new(vec![[f64::NAN, 0.0, 0.0]], vec![]).is_err());
    assert!(validate(&cube(), -1.0).is_err());
}

#[test]
fn cube_missing_a_triangle_has_three_boundary_edges() {
    let c = cube();
    let mut t = c.triangles().to_vec();
    t.pop();
    let r = validate(&mesh(c.vertices(), &t), 0.0).unwrap();
    assert_eq!(r.summary.boundary_edges, 3);
    assert!(!r.summary.is_closed && r.summary.inward_shells == 0);
    let p = &r.problems[0];
    let face = t[p.face.unwrap() as usize];
    assert_eq!(face[p.edge.unwrap() as usize], p.vertex.unwrap());
}

#[test]
fn inward_cube_is_one_inward_shell_and_repairs_outward() {
    let c = cube();
    let m = mesh(c.vertices(), &flipped(c.triangles()));
    assert_eq!(kinds(&m, 0.0), vec![(InwardShell, Some(0), None)]);
    assert!(validate(&m, 0.0).unwrap().problems.len() == 1);
    assert!(unit_volume(&fix(&m, 0.0).0, 1e-12));
    for weld in [-1.0, f64::NAN] {
        let opts = RepairOptions {
            weld: Some(weld),
            ..RepairOptions::default()
        };
        assert!(repair(&c, &opts).is_err(), "{weld}");
    }
}

#[test]
fn each_kind_with_boitata_faces_and_vertices() {
    let c = cube();
    let mut v = c.vertices().to_vec();
    v.push([1.0 + 1e-7, 1.0, 1.0]);
    let mut t = c.triangles().to_vec();
    t.push(t[5]);
    t.push([0, 0, 1]);
    let m = mesh(&v, &t);
    assert_eq!(
        kinds(&m, 1e-6),
        vec![
            (DegenerateFace, Some(13), None),
            (DuplicateFace, Some(12), None),
            (DuplicateVertex, None, Some(8)),
            (NonManifoldEdge, Some(2), Some(4)),
            (NonManifoldEdge, Some(4), Some(1)),
            (NonManifoldEdge, Some(5), Some(1)),
        ]
    );
    assert!(kinds(&m, 0.0).iter().all(|k| k.0 != DuplicateVertex));
    let s = summary(&m, 1e-6);
    assert_eq!((s.unused_vertices, s.coplanar_overlaps), (1, 0));
}

#[test]
fn tetrahedra_touching_at_a_vertex_have_one_non_manifold_vertex() {
    let tet: [[u32; 3]; 4] = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
    let mut v = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    v.extend([[-1., 0., 0.], [0., -1., 0.], [0., 0., -1.]]);
    let second = tet.map(|t| t.map(|i| if i == 0 { 0 } else { i + 3 }));
    let second = second.map(|[a, b, c]| [a, c, b]);
    let r = validate(&mesh(&v, &[tet, second].concat()), 0.0).unwrap();
    assert_eq!(
        r.problems
            .iter()
            .map(|p| (p.kind, p.vertex))
            .collect::<Vec<_>>(),
        vec![(NonManifoldVertex, Some(0))]
    );
    assert_eq!(r.summary.shells, 2);
    assert!(r.summary.is_closed);
}

fn cavity() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<[u32; 3]>) {
    let c = cube();
    let mut v: Vec<[f64; 3]> = c.vertices().iter().map(|p| p.map(|x| 3.0 * x)).collect();
    v.extend(c.vertices().iter().map(|p| p.map(|x| x + 1.0)));
    let outer = c.triangles().to_vec();
    let inner: Vec<[u32; 3]> = outer.iter().map(|t| t.map(|i| i + 8)).collect();
    (v, outer, inner)
}

#[test]
fn cavity_wound_inward_is_clean_and_outward_is_an_inward_shell() {
    let (v, outer, inner) = cavity();
    let good = mesh(&v, &[outer.clone(), flipped(&inner)].concat());
    assert!(validate(&good, 0.0).unwrap().problems.is_empty());
    let bad = mesh(&v, &[outer, inner].concat());
    assert_eq!(kinds(&bad, 0.0), vec![(InwardShell, Some(12), None)]);
}

#[test]
fn repair_keeps_a_cavity_wound_inward_from_any_start() {
    let (v, outer, inner) = cavity();
    let correct = [outer.clone(), flipped(&inner)].concat();
    let (fixed, report) = fix(&mesh(&v, &correct), 0.0);
    assert_eq!(fixed.triangles(), &correct[..]);
    assert!(report.actions.is_empty());
    assert!((volume(&fixed).unwrap() - 26.0).abs() < 1e-12);
    for start in [
        [outer.clone(), inner.clone()],
        [flipped(&outer), inner.clone()],
        [flipped(&outer), flipped(&inner)],
    ] {
        let (m, _) = fix(&mesh(&v, &start.concat()), 0.0);
        assert!(validate(&m, 0.0).unwrap().problems.is_empty());
        assert!((volume(&m).unwrap() - 26.0).abs() < 1e-12);
    }
}

#[test]
fn missing_side_opens_four_edges_and_fills_with_two_faces() {
    let c = cube();
    let open = mesh(c.vertices(), &c.triangles()[2..]);
    let s = summary(&open, 0.0);
    assert_eq!((s.boundary_edges, s.is_closed), (4, false));
    let opts = RepairOptions {
        fill_holes: Some(HoleFill::Fan),
        ..RepairOptions::default()
    };
    let (filled, report) = repair(&open, &opts).unwrap();
    assert_eq!(filled.n_faces(), 12);
    assert_eq!(report.holes_left, 0);
    assert!(validate(&filled, 0.0).unwrap().problems.is_empty());
    assert!(unit_volume(&filled, 1e-12));
    let (again, report) = repair(&filled, &opts).unwrap();
    assert_eq!(again.triangles(), filled.triangles());
    assert!(report.actions.is_empty());
}

/// A cube as a triangle soup (each triangle its own nearly equal corners), half
/// its faces flipped, with a sliver and a repeated face.
#[test]
fn repair_rebuilds_a_broken_cube() {
    let c = cube();
    let mut vertices = vec![];
    let mut triangles = vec![];
    for (t, tri) in c.triangles().iter().enumerate() {
        let n = vertices.len() as u32;
        vertices.extend(tri.map(|i| c.vertices()[i as usize].map(|x| x + 1e-7 * t as f64)));
        triangles.push(if t % 2 == 0 {
            [n, n + 1, n + 2]
        } else {
            [n, n + 2, n + 1]
        });
    }
    triangles.push(triangles[3]);
    triangles.push([0, 1, 1]);
    let broken = mesh(&vertices, &triangles);
    assert!(!summary(&broken, 0.0).is_closed);

    let (loose, _) = fix(&broken, 1e-9);
    assert_eq!(loose.n_faces(), 12);
    assert!(!summary(&loose, 0.0).is_closed);

    let (fixed, report) = fix(&broken, 1e-5);
    assert_eq!((fixed.n_vertices(), fixed.n_faces()), (8, 12));
    assert!(summary(&fixed, 0.0).is_closed);
    assert!(unit_volume(&fixed, 1e-5));
    let kept: Vec<Option<u32>> = (0..12).map(Some).chain([None, None]).collect();
    assert_eq!(report.face_map, kept);

    let (again, report) = fix(&fixed, 1e-5);
    assert_eq!(again.vertices(), fixed.vertices());
    assert_eq!(again.triangles(), fixed.triangles());
    assert!(report.actions.is_empty());
}

#[test]
fn crossing_triangles_are_one_self_intersection() {
    let v = [
        [0., 0., 0.],
        [2., 0., 0.],
        [0., 2., 0.],
        [0.5, 0.5, -1.],
        [0.5, 0.5, 1.],
        [3., 3., 0.5],
    ];
    let r = validate(&mesh(&v, &[[0, 1, 2], [3, 4, 5]]), 0.0).unwrap();
    assert_eq!(r.summary.self_intersections, 1);
    assert_eq!(
        r.problems.last(),
        Some(&MeshProblem {
            kind: SelfIntersection,
            face: Some(1),
            vertex: None,
            edge: None,
            other: Some(0),
        })
    );
    let touching = summary(&mesh(&v, &[[0, 1, 2], [0, 4, 5]]), 0.0);
    assert_eq!(touching.self_intersections, 0);
}

#[test]
fn overlapping_cubes_cross_and_repair_keeps_the_crossing() {
    let m = two_cubes([0.5, 0.25, 0.25]);
    let s = summary(&m, 0.0);
    assert!(s.self_intersections > 0 && s.is_closed);
    let again = summary(&fix(&m, 0.0).0, 0.0);
    assert_eq!(again.self_intersections, s.self_intersections);
    assert_eq!(
        summary(&two_cubes([2.0, 0.0, 0.0]), 0.0).self_intersections,
        0
    );
}

fn rotation([a, b, c]: [f64; 3]) -> [[f64; 3]; 3] {
    let (sa, ca) = a.sin_cos();
    let (sb, cb) = b.sin_cos();
    let (sc, cc) = c.sin_cos();
    [
        [cc * cb, cc * sb * sa - sc * ca, cc * sb * ca + sc * sa],
        [sc * cb, sc * sb * sa + cc * ca, sc * sb * ca - cc * sa],
        [-sb, cb * sa, cb * ca],
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// A rotated, shifted cube stays clean; two cubes cross only while they
    /// overlap.
    #[test]
    fn crossings_follow_overlap(
        angles in proptest::array::uniform3(-3.0..3.0f64),
        offset in proptest::array::uniform3(-1e3..1e3f64),
        shift in 0.1..1.9f64,
    ) {
        let r = rotation(angles);
        let pair = two_cubes([shift, 0.3, 0.2]);
        let moved: Vec<[f64; 3]> = pair
            .vertices()
            .iter()
            .map(|p| r.map(|row| row[0] * p[0] + row[1] * p[1] + row[2] * p[2]))
            .map(|p| [0, 1, 2].map(|k| p[k] + offset[k]))
            .collect();
        let s = summary(&mesh(&moved, pair.triangles()), 0.0);
        prop_assert_eq!(s.self_intersections > 0, shift < 1.0);
        let one = mesh(&moved, &pair.triangles()[..12]);
        prop_assert_eq!(kinds(&one, 0.0), vec![]);
        prop_assert_eq!(summary(&one, 0.0).unused_vertices, 8);
    }
}
