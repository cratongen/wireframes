use proptest::prelude::*;
use wireframes_boolean::{
    BooleanOptions, clip_by_surface, cut_fill, slice_by_plane, split_by_plane,
};
use wireframes_core::fixtures;
use wireframes_core::{Column, Mesh, Polylines};

mod support;
use support::closed_manifold;

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

fn area(m: &Mesh) -> f64 {
    (0..m.n_faces() as u32)
        .map(|f| {
            let [a, b, c] = m.corners(f);
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
        })
        .sum()
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * 1f64.max(a.abs()).max(b.abs())
}

fn closed(m: &Mesh) -> bool {
    m.is_empty() || closed_manifold(m)
}

/// A flat open grid surface at height `z(x, y)` over `[lo, hi]^2`, normals up.
fn sheet(n: u32, lo: f64, hi: f64, z: impl Fn(f64, f64) -> f64) -> Mesh {
    let step = (hi - lo) / (n - 1) as f64;
    let mut v = vec![];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (lo + i as f64 * step, lo + j as f64 * step);
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

#[test]
fn clip_cube_by_flat_sheet() {
    let cube = fixtures::cube();
    let s = sheet(3, -1.0, 2.0, |_, _| 0.25);
    let r = clip_by_surface(&cube, &s, &BooleanOptions::default()).unwrap();
    assert!(closed(&r.below) && closed(&r.above));
    assert!(close(volume(&r.below), 0.25, 1e-12), "{}", volume(&r.below));
    assert!(close(volume(&r.above), 0.75, 1e-12));
}

#[test]
fn clip_through_a_vertex_plane_and_coplanar_face() {
    let cube = fixtures::cube();
    let s = sheet(4, -1.0, 2.0, |_, _| 1.0);
    let r = clip_by_surface(&cube, &s, &BooleanOptions::default()).unwrap();
    assert!(close(volume(&r.below), 1.0, 1e-12));
    assert!(r.above.is_empty());
}

#[test]
fn clip_surface_that_misses() {
    let cube = fixtures::cube();
    let s = sheet(3, -1.0, 2.0, |_, _| 5.0);
    let r = clip_by_surface(&cube, &s, &BooleanOptions::default()).unwrap();
    assert!(close(volume(&r.below), 1.0, 1e-12) && r.above.is_empty());
    let s = sheet(3, -1.0, 2.0, |_, _| -5.0);
    let r = clip_by_surface(&cube, &s, &BooleanOptions::default()).unwrap();
    assert!(close(volume(&r.above), 1.0, 1e-12) && r.below.is_empty());
}

#[test]
fn inclined_splits_through_every_sphere_vertex_cap_the_section() {
    let m = fixtures::sphere(2);
    for &o in m.vertices() {
        for n in [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 2.0, 3.0],
        ] {
            let (a, b) = split_by_plane(&m, o, n).unwrap_or_else(|e| panic!("{o:?} {n:?}: {e}"));
            assert!(closed(&a) && closed(&b), "{o:?} {n:?}");
            assert!(
                close(volume(&a) + volume(&b), volume(&m), 1e-9),
                "{o:?} {n:?}"
            );
        }
    }
}

#[test]
fn clip_refuses_a_solid_beyond_the_surface_extent() {
    let s = sheet(3, -1.0, 2.0, |_, _| 0.0);
    for z in [50.0, -50.0] {
        let v: Vec<[f64; 3]> = fixtures::cube()
            .vertices()
            .iter()
            .map(|p| [p[0] + 100.0, p[1], p[2] + z])
            .collect();
        let far = fixtures::cube().with_vertices(v).unwrap();
        let e = clip_by_surface(&far, &s, &BooleanOptions::default())
            .err()
            .unwrap();
        assert!(e.to_string().contains("does not reach"), "{e}");
    }
    let v: Vec<[f64; 3]> = fixtures::cube()
        .vertices()
        .iter()
        .map(|p| [p[0] + 0.25, p[1], p[2] + 3.0])
        .collect();
    let r = clip_by_surface(
        &fixtures::cube().with_vertices(v).unwrap(),
        &s,
        &BooleanOptions::default(),
    )
    .unwrap();
    assert!(close(volume(&r.above), 1.0, 1e-12) && r.below.is_empty());
}

#[test]
fn clip_refuses_a_surface_ending_inside() {
    let cube = fixtures::cube();
    let s = sheet(3, 0.25, 0.75, |_, _| 0.5);
    let e = clip_by_surface(&cube, &s, &BooleanOptions::default())
        .err()
        .unwrap();
    assert!(e.to_string().contains("does not cut through"), "{e}");
}

#[test]
fn clip_sphere_by_pit_surface() {
    let pit = fixtures::pit_surface(9, 3);
    let v: Vec<[f64; 3]> = fixtures::sphere(2)
        .vertices()
        .iter()
        .map(|p| [50.0 + 30.0 * p[0], 50.0 + 30.0 * p[1], 85.0 + 30.0 * p[2]])
        .collect();
    let solid = fixtures::sphere(2).with_vertices(v).unwrap();
    let r = clip_by_surface(&solid, &pit, &BooleanOptions::default()).unwrap();
    assert!(closed(&r.below) && closed(&r.above));
    assert!(close(
        volume(&r.below) + volume(&r.above),
        volume(&solid),
        1e-9
    ));
    assert!(volume(&r.below) > 0.0 && volume(&r.above) > 0.0);
}

#[test]
fn split_and_slice_a_box() {
    let b = fixtures::box_mesh([0.0; 3], [2.0, 3.0, 4.0]);
    let (n, p) = split_by_plane(&b, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]).unwrap();
    assert!(n.is_closed() && p.is_closed());
    assert!(close(volume(&n), 6.0, 1e-12) && close(volume(&p), 18.0, 1e-12));
    let s = slice_by_plane(&b, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], &[0.0, 2.0, 10.0]).unwrap();
    assert_eq!(s.len(), 2);
    assert!(s.closed().iter().all(|c| *c));
    for a in s.area() {
        assert!(close(a, 6.0, 1e-12));
    }
    let ring = s.part(0);
    let signed: f64 = (0..ring.len())
        .map(|i| {
            let (a, c) = (ring[i], ring[(i + 1) % ring.len()]);
            a[0] * c[1] - a[1] * c[0]
        })
        .sum();
    assert!(
        signed > 0.0,
        "ring must wind counter-clockwise seen from +z"
    );
    assert_eq!(
        s.attributes().column("offset").unwrap().data,
        Column::f64(vec![0.0, 2.0]).data
    );
}

#[test]
fn split_through_vertices_keeps_volume() {
    let b = fixtures::cube();
    for z in [0.0, 1.0] {
        let (n, p) = split_by_plane(&b, [0.0, 0.0, z], [0.0, 0.0, 1.0]).unwrap();
        assert!(close(volume(&n) + volume(&p), 1.0, 1e-12));
    }
    let (n, p) = split_by_plane(&b, [0.0, 0.0, 0.0], [1.0, 1.0, 0.0]).unwrap();
    assert!(n.is_empty() && close(volume(&p), 1.0, 1e-12));
    let (n, p) = split_by_plane(&b, [0.5, 0.5, 0.5], [1.0, -1.0, 0.0]).unwrap();
    assert!(n.is_closed() && p.is_closed());
    assert!(close(volume(&n), 0.5, 1e-12) && close(volume(&p), 0.5, 1e-12));
}

#[test]
fn slice_splits_lines_where_face_attributes_change() {
    let b = fixtures::cube()
        .with_face_column(
            "zone",
            Column::i64((0..12).map(|f| (f / 4) as i64).collect()),
        )
        .unwrap();
    let s = slice_by_plane(&b, [0.0, 0.0, 0.5], [0.0, 0.0, 1.0], &[]).unwrap();
    assert!(s.len() > 1);
    assert!(s.closed().iter().all(|c| !*c));
    let named = fixtures::cube()
        .with_face_column("offset", Column::f64(vec![0.0; 12]))
        .unwrap();
    assert!(slice_by_plane(&named, [0.0; 3], [0.0, 0.0, 1.0], &[]).is_err());
}

#[test]
fn torus_slice_gives_ring_and_hole() {
    let t = fixtures::torus(24, 12);
    let s = slice_by_plane(&t, [0.0; 3], [0.0, 0.0, 1.0], &[]).unwrap();
    assert_eq!(s.len(), 2);
    assert!(s.closed().iter().all(|c| *c));
    let signed = |r: &[[f64; 3]]| -> f64 {
        (0..r.len())
            .map(|i| {
                let (a, c) = (r[i], r[(i + 1) % r.len()]);
                a[0] * c[1] - a[1] * c[0]
            })
            .sum()
    };
    let (a, b) = (signed(s.part(0)), signed(s.part(1)));
    assert!(a * b < 0.0, "outer ring and hole wind opposite ways");
}

#[test]
fn cut_fill_of_shifted_surfaces() {
    let base = sheet(5, 0.0, 4.0, |x, y| 0.1 * x + 0.05 * y);
    let up = sheet(5, 0.0, 4.0, |x, y| 0.1 * x + 0.05 * y + 2.0);
    let r = cut_fill(&base, &up, None).unwrap();
    assert!(close(r.fill, 32.0, 1e-9) && r.cut == 0.0);
    assert!(close(r.area, 16.0, 1e-12) && close(r.fill_area, 16.0, 1e-12));
    let r2 = cut_fill(&up, &base, None).unwrap();
    assert!(close(r2.cut, r.fill, 1e-12) && r2.fill == 0.0);
    let z = cut_fill(&base, &base, None).unwrap();
    assert!(z.cut == 0.0 && z.fill == 0.0 && close(z.area, 16.0, 1e-12));
}

#[test]
fn cut_fill_crossing_planes_and_boundary() {
    let base = sheet(3, 0.0, 2.0, |_, _| 0.0);
    let tilt = sheet(4, 0.0, 2.0, |x, _| x - 1.0);
    let r = cut_fill(&base, &tilt, None).unwrap();
    assert!(
        close(r.fill, 1.0, 1e-12) && close(r.cut, 1.0, 1e-12),
        "{r:?}"
    );
    assert!(close(r.cut_area, 2.0, 1e-12) && close(r.fill_area, 2.0, 1e-12));
    let square = Polylines::from_parts(
        vec![vec![
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [1.0, 2.0, 0.0],
        ]],
        vec![true],
        None,
    )
    .unwrap();
    let r = cut_fill(&base, &tilt, Some(&square)).unwrap();
    assert!(close(r.fill, 1.0, 1e-12) && r.cut.abs() < 1e-15 && close(r.area, 2.0, 1e-12));
}

#[test]
fn cut_fill_refuses_folded_surfaces() {
    let mut v = sheet(3, 0.0, 2.0, |_, _| 0.0).vertices().to_vec();
    v.push([0.5, 0.5, 1.0]);
    v.push([1.5, 0.5, 1.0]);
    v.push([0.5, 1.5, 1.0]);
    let mut t = sheet(3, 0.0, 2.0, |_, _| 0.0).triangles().to_vec();
    t.push([9, 10, 11]);
    let folded = Mesh::new(v, t).unwrap();
    let flat = sheet(3, 0.0, 2.0, |_, _| 0.0);
    assert!(cut_fill(&folded, &flat, None).is_err());
}

#[test]
fn split_cap_with_rounded_collinear_points_stays_closed() {
    let m = fixtures::torus(10, 6);
    let o = [0.0, 0.981726445642584, 0.7596377723565166];
    let (a, b) = split_by_plane(&m, o, [0.0, 1.0, 0.9190331232657712]).unwrap();
    assert!(closed(&a) && closed(&b));
    assert!(close(volume(&a) + volume(&b), volume(&m), 1e-9));
    let slivers = |x: &Mesh| {
        wireframes_core::validate::degenerate_mask(x)
            .iter()
            .filter(|d| **d)
            .count()
    };
    println!(
        "cap slivers after rounding: {} + {}",
        slivers(&a),
        slivers(&b)
    );
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 32, ..ProptestConfig::default() })]

    #[test]
    fn split_volumes_add_up(
        o in proptest::array::uniform3(prop_oneof![(-2i32..3).prop_map(|x| x as f64 * 0.5), -1.0..1.0f64]),
        n in proptest::array::uniform3(prop_oneof![Just(0.0), Just(1.0), -1.0..1.0f64]),
        which in 0usize..3,
    ) {
        prop_assume!(n.iter().any(|x| *x != 0.0));
        let m = match which { 0 => fixtures::cube(), 1 => fixtures::sphere(1), _ => fixtures::torus(10, 6) };
        let (a, b) = split_by_plane(&m, o, n).unwrap();
        prop_assert!(closed(&a) && closed(&b));
        prop_assert!(close(volume(&a) + volume(&b), volume(&m), 1e-9), "{} + {} vs {}", volume(&a), volume(&b), volume(&m));
        let open = fixtures::open_cube();
        let (a, b) = split_by_plane(&open, o, n).unwrap();
        prop_assert!(close(area(&a) + area(&b), area(&open), 1e-9));
    }

    #[test]
    fn slices_of_closed_meshes_are_rings(z in -1.2..1.2f64, which in 0usize..2) {
        let m = if which == 0 { fixtures::sphere(2) } else { fixtures::torus(12, 8) };
        let s = slice_by_plane(&m, [0.0, 0.0, z], [0.0, 0.0, 1.0], &[]).unwrap();
        prop_assert!(s.closed().iter().all(|c| *c));
    }

    #[test]
    fn clip_volumes_add_up(z in -0.2..1.2f64, tilt in -0.5..0.5f64, angles in proptest::array::uniform3(-3.2..3.2f64)) {
        let base = fixtures::sphere(1);
        let s = sheet(5, -3.0, 3.0, |x, y| z + tilt * x - 0.3 * tilt * y);
        let _ = angles;
        let r = clip_by_surface(&base, &s, &BooleanOptions::default()).unwrap();
        prop_assert!(closed(&r.below) && closed(&r.above));
        prop_assert!(close(volume(&r.below) + volume(&r.above), volume(&base), 1e-9));
    }

    #[test]
    fn cut_fill_swaps_and_shifts(h in 0.0..5.0f64, seed in 0u64..50) {
        let base = fixtures::pit_surface(6, seed);
        let v: Vec<[f64; 3]> = base.vertices().iter().map(|p| [p[0], p[1], p[2] + h]).collect();
        let up = base.clone().with_vertices(v).unwrap();
        let r = cut_fill(&base, &up, None).unwrap();
        prop_assert!(close(r.fill, h * 10_000.0, 1e-9), "{} vs {}", r.fill, h * 1e4);
        let s = cut_fill(&up, &base, None).unwrap();
        prop_assert!(close(s.cut, r.fill, 1e-12) && close(s.fill, r.cut, 1e-12));
        let z = cut_fill(&base, &base, None).unwrap();
        prop_assert!(z.cut == 0.0 && z.fill == 0.0);
    }
}

#[test]
fn clip_and_cut_fill_do_not_depend_on_thread_count() {
    let pit = fixtures::pit_surface(9, 3);
    let v: Vec<[f64; 3]> = fixtures::sphere(2)
        .vertices()
        .iter()
        .map(|p| [50.0 + 30.0 * p[0], 50.0 + 30.0 * p[1], 85.0 + 30.0 * p[2]])
        .collect();
    let solid = fixtures::sphere(2).with_vertices(v).unwrap();
    let flat = sheet(7, 0.0, 100.0, |_, _| 90.0);
    let go = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| {
                let c = clip_by_surface(&solid, &pit, &BooleanOptions::default()).unwrap();
                (c.below, c.above, cut_fill(&pit, &flat, None).unwrap())
            })
    };
    let (a, b) = (go(1), go(4));
    assert_eq!(a.0, b.0);
    assert_eq!(a.1, b.1);
    assert_eq!(a.2, b.2);
}

/// A frustum over a random convex `sides`-gon: bottom at `z = 0`, top scaled
/// by `top` at `z = 1`. Each side is cut into `k` collinear pieces, so each
/// wall is a strip of triangles that is planar only up to the rounding of its
/// corners, and a section crosses it in a long, nearly collinear chain.
fn frustum(seed: u64, sides: usize, k: usize, top: f64) -> Mesh {
    let mut rng = fixtures::Rng::new(seed);
    let mut angles: Vec<f64> = (0..sides)
        .map(|_| rng.uniform(0.0, std::f64::consts::TAU))
        .collect();
    angles.sort_by(f64::total_cmp);
    let corner = |a: f64| [100.0 * a.cos(), 100.0 * a.sin()];
    let ring: Vec<[f64; 2]> = (0..sides)
        .flat_map(|i| {
            let (p, q) = (corner(angles[i]), corner(angles[(i + 1) % sides]));
            (0..k).map(move |j| {
                let s = j as f64 / k as f64;
                [p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])]
            })
        })
        .collect();
    let mut v = vec![];
    for (s, z) in [(1.0, 0.0), (top, 1.0)] {
        v.extend(ring.iter().map(|p| [s * p[0], s * p[1], z]));
    }
    let c = angles.iter().map(|&a| corner(a)).fold([0.0; 2], |s, p| {
        [s[0] + p[0] / sides as f64, s[1] + p[1] / sides as f64]
    });
    v.push([c[0], c[1], 0.0]);
    v.push([top * c[0], top * c[1], 1.0]);
    let n = ring.len() as u32;
    let mut t = vec![];
    for i in 0..n {
        let j = (i + 1) % n;
        t.push([2 * n, j, i]);
        t.push([2 * n + 1, n + i, n + j]);
        t.push([i, j, n + j]);
        t.push([i, n + j, n + i]);
    }
    Mesh::new(v, t).unwrap()
}

#[test]
fn split_caps_of_nearly_collinear_sections_are_valid_solids() {
    for seed in 0..60 {
        let m = frustum(seed, 5, 9, 1.7);
        assert!(m.validate(0.0).unwrap().summary.is_solid(), "seed {seed}");
        for z in [0.1, 0.3, 0.7] {
            let (a, b) = split_by_plane(&m, [0.0, 0.0, z], [0.0, 0.0, 1.0]).unwrap();
            for part in [&a, &b] {
                let s = part.validate(0.0).unwrap().summary;
                assert!(s.is_solid(), "seed {seed}, z {z}: {s:?}");
            }
        }
    }
}

#[test]
fn split_side_test_survives_underflow() {
    let tet = Mesh::new(
        vec![
            [-1e-200, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
    .unwrap();
    assert!(tet.validate(0.0).unwrap().summary.is_solid());
    let (neg, pos) = split_by_plane(&tet, [0.0; 3], [1e-200, 0.0, 0.0]).unwrap();
    assert!(
        neg.n_faces() > 0,
        "vertex at x = -1e-200 lies strictly on the negative side of x = 0"
    );
    let unit = split_by_plane(&tet, [0.0; 3], [1.0, 0.0, 0.0]).unwrap();
    assert_eq!((&neg, &pos), (&unit.0, &unit.1));
    let lines = slice_by_plane(&tet, [0.0; 3], [1e-200, 0.0, 0.0], &[]).unwrap();
    assert_eq!(lines.n_parts(), 1);
}

#[test]
fn split_of_a_sub_underflow_sliver_is_closed() {
    let tet = Mesh::new(
        vec![
            [-1e-200, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
    .unwrap();
    let (neg, pos) = split_by_plane(&tet, [0.0; 3], [1.0, 0.0, 0.0]).unwrap();
    assert!(closed_manifold(&neg) && closed_manifold(&pos));
}

fn height_field_solid(top: &[[f64; 5]; 5]) -> Mesh {
    let mut v = vec![];
    for level in [None, Some(top)] {
        for (i, row) in [0.0f64, 1.0, 2.0, 3.0, 4.0].iter().enumerate() {
            for (j, y) in [0.0f64, 1.0, 2.0, 3.0, 4.0].iter().enumerate() {
                v.push([*row, *y, level.map_or(0.0, |h| h[i][j])]);
            }
        }
    }
    let at = |l: u32, i: u32, j: u32| l * 25 + i * 5 + j;
    let mut t = vec![];
    for i in 0..4 {
        for j in 0..4 {
            let [a, b, c, d] = [
                at(1, i, j),
                at(1, i + 1, j),
                at(1, i + 1, j + 1),
                at(1, i, j + 1),
            ];
            t.push([a, b, c]);
            t.push([a, c, d]);
            let [a, b, c, d] = [
                at(0, i, j),
                at(0, i + 1, j),
                at(0, i + 1, j + 1),
                at(0, i, j + 1),
            ];
            t.push([a, c, b]);
            t.push([a, d, c]);
        }
    }
    let ring: Vec<(u32, u32)> = (0..4)
        .map(|k| (k, 0))
        .chain((0..4).map(|k| (4, k)))
        .chain((0..4).map(|k| (4 - k, 4)))
        .chain((0..4).map(|k| (0, 4 - k)))
        .collect();
    for k in 0..ring.len() {
        let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
        let [a, b] = [at(0, p.0, p.1), at(0, q.0, q.1)];
        let [c, d] = [at(1, q.0, q.1), at(1, p.0, p.1)];
        t.push([a, b, c]);
        t.push([a, c, d]);
    }
    Mesh::new(v, t).unwrap()
}

#[test]
fn split_at_a_level_the_top_touches_from_above_gives_solids() {
    let mut pit = [[2.0; 5]; 5];
    pit[2][2] = 1.0;
    let mut valley = [[2.0; 5]; 5];
    valley[1][2] = 1.0;
    valley[2][2] = 1.0;
    valley[3][2] = 1.0;
    let mut two_pits = [[2.0; 5]; 5];
    two_pits[1][1] = 1.0;
    two_pits[3][3] = 1.0;
    for top in [pit, valley, two_pits] {
        let m = height_field_solid(&top);
        assert!(m.validate(0.0).unwrap().summary.is_solid());
        let (lower, upper) = split_by_plane(&m, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]).unwrap();
        for part in [&lower, &upper] {
            let s = part.validate(0.0).unwrap().summary;
            assert!(s.is_solid(), "{top:?}: {s:?}");
        }
        assert!((volume(&lower) - 16.0).abs() < 1e-12);
        assert!((volume(&lower) + volume(&upper) - volume(&m)).abs() < 1e-12);
        let (upper_flipped, lower_flipped) =
            split_by_plane(&m, [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]).unwrap();
        assert!((volume(&upper_flipped) - volume(&upper)).abs() < 1e-12);
        assert!((volume(&lower_flipped) - 16.0).abs() < 1e-12);
    }
}
