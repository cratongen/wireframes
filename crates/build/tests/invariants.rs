use std::panic::{AssertUnwindSafe, catch_unwind};

use drusa_build::{
    DecimateOptions, DtmOptions, LoftOptions, close_to_elevation, contours, decimate, dtm, loft,
    offset, remesh,
};
use drusa_core::fixtures::{box_mesh, cube, pit_surface, sphere, torus};
use drusa_core::{Error, Mesh, Polylines};
use drusa_query::SurfaceIndex;
use proptest::prelude::*;

#[test]
fn dtm_tiny_min_spacing_does_not_panic() {
    let pts = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
    let r = catch_unwind(|| {
        dtm(
            &pts,
            &DtmOptions {
                min_spacing: 1e-300,
                ..Default::default()
            },
        )
    });
    assert!(r.is_ok(), "dtm panicked for min_spacing = 1e-300");
}

#[test]
fn dtm_crossing_at_a_vertex_of_the_second_breakline_takes_z_of_the_first() {
    let pts = [
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.0, 4.0, 0.0],
        [0.0, 4.0, 0.0],
    ];
    let lines = Polylines::from_parts(
        vec![
            vec![[1.0, 2.0, 10.0], [3.0, 2.0, 20.0]],
            vec![[2.0, 1.0, -5.0], [2.0, 2.0, -5.0], [2.0, 3.0, -5.0]],
        ],
        vec![false, false],
        None,
    )
    .unwrap();
    let m = dtm(
        &pts,
        &DtmOptions {
            breaklines: Some(&lines),
            ..Default::default()
        },
    )
    .unwrap();
    let at = m
        .vertices()
        .iter()
        .find(|p| p[0] == 2.0 && p[1] == 2.0)
        .unwrap();
    assert_eq!(
        at[2], 15.0,
        "crossing of the two breaklines at (2, 2) should take z along the first breakline"
    );
}

#[test]
fn decimate_face_with_repeated_vertex_does_not_panic() {
    let v = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
    ];
    let m = Mesh::new(v, vec![[0, 1, 2], [0, 2, 3], [1, 4, 4]]).unwrap();
    let r = catch_unwind(AssertUnwindSafe(|| {
        decimate(
            &m,
            &DecimateOptions {
                target_faces: Some(1),
                ..Default::default()
            },
        )
    }));
    assert!(
        r.is_ok(),
        "decimate panicked on a face with a repeated vertex"
    );
    assert!(matches!(r.unwrap(), Err(Error::InvalidInput(_))));
    let r = catch_unwind(AssertUnwindSafe(|| remesh(&m, 0.5, 2)));
    assert!(
        r.is_ok(),
        "remesh panicked on a face with a repeated vertex"
    );
    assert!(matches!(r.unwrap(), Err(Error::InvalidInput(_))));
}

#[test]
fn dtm_point_on_a_breakline_keeps_its_z() {
    let pts = [
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.0, 4.0, 0.0],
        [0.0, 4.0, 0.0],
        [2.0, 2.0, 7.0],
    ];
    let lines = Polylines::from_parts(
        vec![vec![[1.0, 2.0, 10.0], [3.0, 2.0, 20.0]]],
        vec![false],
        None,
    )
    .unwrap();
    let m = dtm(
        &pts,
        &DtmOptions {
            breaklines: Some(&lines),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(&m.vertices()[..pts.len()], &pts[..]);
}

#[test]
fn dtm_tiny_min_spacing_on_utm_coordinates_keeps_distinct_points() {
    let pts = [
        [7.0e6, 7.0e6, 0.0],
        [7.0e6 + 1.0, 7.0e6, 1.0],
        [7.0e6, 7.0e6 + 1.0, 2.0],
        [7.0e6, 7.0e6, 3.0],
    ];
    let m = dtm(
        &pts,
        &DtmOptions {
            min_spacing: 1e-13,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(m.vertices(), &pts[..3]);
}

#[test]
fn close_to_elevation_bowtie_is_solid_or_refused() {
    let v = vec![
        [0.0, 0.0, 10.0],
        [1.0, 0.0, 10.0],
        [1.0, 1.0, 10.0],
        [2.0, 1.0, 10.0],
        [2.0, 2.0, 10.0],
    ];
    let m = Mesh::new(v, vec![[0, 1, 2], [2, 3, 4]]).unwrap();
    if let Ok(s) = close_to_elevation(&m, 0.0) {
        let r = s.validate(0.0).unwrap();
        assert!(r.summary.is_solid(), "{:?}", r.summary);
    }
}

#[test]
fn close_to_elevation_volume_equals_cut_fill_against_flat_plane() {
    let top = pit_surface(12, 3);
    let flat = Mesh::new(
        top.vertices().iter().map(|p| [p[0], p[1], 0.0]).collect(),
        top.triangles().to_vec(),
    )
    .unwrap();
    let s = close_to_elevation(&top, 0.0).unwrap();
    let vol = drusa_ops::volume(&s).unwrap();
    let cf = drusa_boolean::cut_fill(&flat, &top, None).unwrap();
    let net = cf.net.abs();
    assert!(
        (vol - net).abs() <= 1e-9 * vol.max(1.0),
        "volume {vol} vs cut_fill net {net} ({cf:?})"
    );
}

#[test]
fn contour_through_a_row_of_vertices_is_one_straight_line_on_the_plane() {
    let n = 5u32;
    let mut v = vec![];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (i as f64, j as f64);
            v.push([x, y, x]);
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
    let m = Mesh::new(v, f).unwrap();
    let p = contours(&m, &[2.0]).unwrap();
    assert_eq!(p.len(), 1, "expected one contour line, got {}", p.len());
    for q in p.part(0) {
        assert_eq!(q[0], 2.0);
        assert_eq!(q[2], 2.0);
    }
}

#[test]
fn loft_squares_with_different_vertex_counts_is_a_box() {
    let a = 3.0;
    let h = 2.0;
    let dense = vec![
        [0.0, 0.0, 0.0],
        [a / 2.0, 0.0, 0.0],
        [a, 0.0, 0.0],
        [a, a / 2.0, 0.0],
        [a, a, 0.0],
        [a / 2.0, a, 0.0],
        [0.0, a, 0.0],
        [0.0, a / 2.0, 0.0],
    ];
    let sparse = vec![[0.0, 0.0, h], [a, 0.0, h], [a, a, h], [0.0, a, h]];
    let p = Polylines::from_parts(vec![dense, sparse], vec![true, true], None).unwrap();
    let m = loft(&p, &LoftOptions::default()).unwrap();
    let got = drusa_ops::volume(&m).unwrap();
    let want = a * a * h;
    assert!((got - want).abs() <= 1e-9 * want, "{got} vs {want}");
}

fn sampled_hausdorff(a: &Mesh, b: &Mesh) -> f64 {
    let samples = |m: &Mesh| -> Vec<[f64; 3]> {
        let mut s = m.vertices().to_vec();
        for t in m.triangles() {
            let c = t.map(|i| m.vertices()[i as usize]);
            s.push([0, 1, 2].map(|k| (c[0][k] + c[1][k] + c[2][k]) / 3.0));
            for e in 0..3 {
                let (p, q) = (c[e], c[(e + 1) % 3]);
                s.push([0, 1, 2].map(|k| (p[k] + q[k]) / 2.0));
            }
        }
        s
    };
    let one =
        |m: &Mesh, idx: &SurfaceIndex| idx.distance(&samples(m)).into_iter().fold(0.0, f64::max);
    one(a, &SurfaceIndex::new(b)).max(one(b, &SurfaceIndex::new(a)))
}

#[test]
fn remesh_keeps_closed_inputs_closed_and_near() {
    let cases: Vec<(&str, Mesh, Vec<f64>)> = vec![
        ("cube", cube(), vec![0.2, 0.35, 0.5, 0.9]),
        (
            "slab",
            box_mesh([0.0; 3], [4.0, 4.0, 0.2]),
            vec![0.3, 0.6, 1.0],
        ),
        ("torus", torus(16, 8), vec![0.3, 0.6]),
        ("sphere", sphere(2), vec![0.15, 0.4]),
    ];
    let mut failures = vec![];
    for (name, m, ls) in cases {
        for l in ls {
            let r = remesh(&m, l, 5).unwrap();
            let rep = r.validate(0.0).unwrap();
            let h = sampled_hausdorff(&m, &r);
            if !r.is_closed() || rep.summary.non_manifold_edges > 0 || h >= l {
                failures.push(format!(
                    "{name} l={l}: closed={} non_manifold_edges={} hausdorff={h}",
                    r.is_closed(),
                    rep.summary.non_manifold_edges
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn decimate_keeps_closed_manifold(which in 0usize..3, frac in 0.05f64..0.9, preserve: bool) {
        let m = match which {
            0 => sphere(2),
            1 => torus(12, 6),
            _ => box_mesh([0.0; 3], [3.0, 1.0, 0.5]),
        };
        let target = ((m.n_faces() as f64) * frac) as usize;
        let d = decimate(&m, &DecimateOptions { target_faces: Some(target), preserve_boundary: preserve, ..Default::default() }).unwrap();
        let s = d.validate(0.0).unwrap().summary;
        prop_assert!(d.is_closed(), "{:?}", s);
        prop_assert_eq!(s.non_manifold_edges, 0);
        prop_assert_eq!(s.non_manifold_vertices, 0);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]
    #[test]
    fn remesh_boxes_stay_closed_and_near(
        sx in 0.1f64..4.0, sy in 0.1f64..4.0, sz in 0.1f64..4.0, l in 0.3f64..1.2,
    ) {
        let m = box_mesh([0.0; 3], [sx, sy, sz]);
        let r = remesh(&m, l, 3).unwrap();
        let s = r.validate(0.0).unwrap().summary;
        prop_assert!(r.is_closed(), "{:?}", s);
        prop_assert_eq!(s.non_manifold_edges, 0);
        let h = sampled_hausdorff(&m, &r);
        prop_assert!(h < l, "hausdorff {} >= edge_length {}", h, l);
    }
}

#[test]
fn remesh_thin_slab_keeps_its_corners_within_edge_length() {
    let m = box_mesh([0.0; 3], [4.0, 4.0, 0.2]);
    let l = 1.0;
    let r = remesh(&m, l, 5).unwrap();
    let in_to_out = SurfaceIndex::new(&r)
        .distance(m.vertices())
        .into_iter()
        .fold(0.0, f64::max);
    let (v0, v1) = (
        drusa_ops::volume(&m).unwrap(),
        drusa_ops::volume(&r).unwrap(),
    );
    assert!(
        in_to_out < l,
        "input corner lies {in_to_out} from the remeshed slab ({} faces, volume {v0} -> {v1})",
        r.n_faces()
    );
}

#[test]
fn results_do_not_depend_on_thread_count() {
    let solid = sphere(2);
    let top = pit_surface(12, 3);
    let run = |n: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| {
                (
                    remesh(&solid, 0.3, 3).unwrap(),
                    decimate(
                        &solid,
                        &DecimateOptions {
                            target_faces: Some(solid.n_faces() / 3),
                            ..Default::default()
                        },
                    )
                    .unwrap(),
                    offset(&solid, 0.05).unwrap(),
                    close_to_elevation(&top, 0.0).unwrap(),
                )
            })
    };
    assert_eq!(run(1), run(4));
}
