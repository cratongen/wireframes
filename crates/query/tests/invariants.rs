use drusa_core::intersect::point_in_triangle;
use drusa_core::predicates::{Sign, orient3d};
use drusa_core::{Mesh, fixtures};
use drusa_query::{Blocks, HitKind, Location, SolidIndex, SurfaceIndex, Trace};
use proptest::prelude::*;

/// Float solid-angle winding number of `m` around `p`.
fn winding(m: &Mesh, p: [f64; 3]) -> f64 {
    (0..m.n_faces() as u32)
        .map(|f| {
            let [a, b, c] = m
                .corners(f)
                .map(|v| [v[0] - p[0], v[1] - p[1], v[2] - p[2]]);
            let n = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            let den = n(a) * n(b) * n(c) + dot(a, b) * n(c) + dot(b, c) * n(a) + dot(c, a) * n(b);
            2.0 * det.atan2(den)
        })
        .sum::<f64>()
        / (4.0 * std::f64::consts::PI)
}

/// Exact winding number of `m` around `p`, which must be off the surface, as the
/// signed count of crossings of the segment from `p` to `far` (exact `orient3d`
/// signs, so a rational computation). `None` when the segment meets an edge, a
/// vertex or a face plane at `far`, or runs in a face plane.
fn exact_winding(m: &Mesh, p: [f64; 3], far: [f64; 3]) -> Option<i32> {
    let mut w = 0;
    for f in 0..m.n_faces() as u32 {
        let [a, b, c] = m.corners(f);
        let (sp, sf) = (
            Sign::of(orient3d(a, b, c, p)),
            Sign::of(orient3d(a, b, c, far)),
        );
        if sf == Sign::Zero {
            return None;
        }
        if sp == Sign::Zero || sp == sf {
            continue;
        }
        let s = [(a, b), (b, c), (c, a)].map(|(u, v)| Sign::of(orient3d(p, far, u, v)));
        if s.contains(&Sign::Positive) && s.contains(&Sign::Negative) {
            continue;
        }
        if s.contains(&Sign::Zero) {
            return None;
        }
        w += if sp == Sign::Positive { 1 } else { -1 };
    }
    Some(w)
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * 1f64.max(a.abs()).max(b.abs())
}

fn volume(m: &Mesh) -> f64 {
    let o = m.bounds().unwrap().min;
    (0..m.n_faces() as u32)
        .map(|f| {
            let [a, b, c] = m
                .corners(f)
                .map(|v| [v[0] - o[0], v[1] - o[1], v[2] - o[2]]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

fn rotation(yaw: f64, pitch: f64, roll: f64) -> [[f64; 3]; 3] {
    let (cy, sy, cp, sp, cr, sr) = (
        yaw.cos(),
        yaw.sin(),
        pitch.cos(),
        pitch.sin(),
        roll.cos(),
        roll.sin(),
    );
    let rz = [[cy, -sy, 0.0], [sy, cy, 0.0], [0.0, 0.0, 1.0]];
    let ry = [[cp, 0.0, sp], [0.0, 1.0, 0.0], [-sp, 0.0, cp]];
    let rx = [[1.0, 0.0, 0.0], [0.0, cr, -sr], [0.0, sr, cr]];
    let mul = |a: [[f64; 3]; 3], b: [[f64; 3]; 3]| {
        [0, 1, 2].map(|i| [0, 1, 2].map(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum::<f64>()))
    };
    mul(mul(rz, ry), rx)
}

/// Block grid in the frame `rotation` around `origin` that covers a ball of radius `reach`.
fn grid(origin: [f64; 3], size: [f64; 3], reach: f64, r: Option<[[f64; 3]; 3]>) -> Blocks {
    let n = |k: usize| (reach / size[k]).ceil() as i64 + 1;
    let rot = r;
    let r = r.unwrap_or([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let mut centroids = vec![];
    for i in -n(0)..=n(0) {
        for j in -n(1)..=n(1) {
            for k in -n(2)..=n(2) {
                let l = [i as f64 * size[0], j as f64 * size[1], k as f64 * size[2]];
                centroids.push(
                    [0, 1, 2].map(|a| origin[a] + (0..3).map(|b| r[a][b] * l[b]).sum::<f64>()),
                );
            }
        }
    }
    Blocks {
        centroids,
        size,
        rotation: rot,
    }
}

fn solids() -> Vec<Mesh> {
    vec![
        fixtures::sphere(2),
        fixtures::torus(12, 8),
        fixtures::tetrahedron(),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn block_fractions_sum_to_solid_volume(
        sx in 0.15f64..0.7, sy in 0.15f64..0.7, sz in 0.15f64..0.7,
        ox in -0.3f64..0.3, oy in -0.3f64..0.3, oz in -0.3f64..0.3,
        yaw in -3.2f64..3.2, pitch in -1.5f64..1.5, roll in -3.2f64..3.2,
        rotated: bool,
    ) {
        for m in solids() {
            let s = SolidIndex::new(&m).unwrap();
            let size = [sx, sy, sz];
            let blocks = grid([ox, oy, oz], size, 2.6, rotated.then(|| rotation(yaw, pitch, roll)));
            let f = s.block_fraction(&blocks).unwrap();
            prop_assert!(f.iter().all(|x| (0.0..=1.0).contains(x)));
            let total: f64 = f.iter().map(|x| x * sx * sy * sz).sum();
            prop_assert!(close(total, volume(&m), 1e-9), "{} vs {}", total, volume(&m));
        }
    }

    #[test]
    fn blocks_inside_and_outside_a_convex_solid(seed in 0u64..1000, h in 0.02f64..0.3) {
        let m = fixtures::sphere(2);
        let s = SolidIndex::new(&m).unwrap();
        let bounds = m.bounds().unwrap();
        let centroids = fixtures::random_points(seed, 40, [-1.5; 3], [1.5; 3]);
        let f = s.block_fraction(&Blocks { centroids: centroids.clone(), size: [h; 3], rotation: None }).unwrap();
        for (c, f) in centroids.iter().zip(&f) {
            let corners: Vec<[f64; 3]> = (0..8)
                .map(|i| [0, 1, 2].map(|k| c[k] + if i >> k & 1 == 1 { h / 2.0 } else { -h / 2.0 }))
                .collect();
            if s.locate_all(&corners).iter().all(|&l| l == Location::Inside) {
                prop_assert!(close(*f, 1.0, 1e-9), "{}", f);
            }
            let disjoint = (0..3).any(|k| c[k] + h / 2.0 < bounds.min[k] || c[k] - h / 2.0 > bounds.max[k]);
            if disjoint {
                prop_assert!(close(*f, 0.0, 1e-9), "{}", f);
            }
        }
    }

    #[test]
    fn trace_inside_length_is_the_chord(
        lo in prop::array::uniform3(-2.0f64..-0.5), hi in prop::array::uniform3(0.5f64..2.0),
        a in prop::array::uniform3(-4.0f64..4.0), b in prop::array::uniform3(-4.0f64..4.0),
        collar in 0.0f64..100.0,
    ) {
        let m = fixtures::box_mesh(lo, hi);
        let s = SolidIndex::new(&m).unwrap();
        let len = (0..3).map(|k| (b[k] - a[k]).powi(2)).sum::<f64>().sqrt();
        prop_assume!(len > 1e-3);
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for k in 0..3 {
            let d = b[k] - a[k];
            if d == 0.0 {
                if a[k] < lo[k] || a[k] > hi[k] { t1 = -1.0; }
                continue;
            }
            let (u, v) = ((lo[k] - a[k]) / d, (hi[k] - a[k]) / d);
            t0 = t0.max(u.min(v));
            t1 = t1.min(u.max(v));
        }
        let chord = if t1 > t0 { (t1 - t0) * len } else { 0.0 };
        let iv = s.trace_intervals(&Trace { xyz: vec![a, b], depth: vec![collar, collar + len] }).unwrap();
        prop_assert_eq!(iv.first().unwrap().from, collar);
        prop_assert_eq!(iv.last().unwrap().to, collar + len);
        for w in iv.windows(2) {
            prop_assert_eq!(w[0].to, w[1].from);
            prop_assert!(w[0].inside != w[1].inside);
            prop_assert!(w[0].to > w[0].from);
        }
        let inside: f64 = iv.iter().filter(|x| x.inside).map(|x| x.to - x.from).sum();
        prop_assert!((inside - chord).abs() <= 1e-9 * (1.0 + collar + len), "{} vs {}", inside, chord);
        prop_assert!(iv.iter().filter(|x| x.inside).count() <= 1);
    }

    #[test]
    fn bent_traces_alternate_and_cover_their_depth(seed in 0u64..1000) {
        let m = fixtures::torus(12, 8);
        let s = SolidIndex::new(&m).unwrap();
        let xyz = fixtures::random_points(seed, 6, [-1.5; 3], [1.5; 3]);
        let mut depth = vec![0.0];
        for w in xyz.windows(2) {
            let d = (0..3).map(|k| (w[1][k] - w[0][k]).powi(2)).sum::<f64>().sqrt();
            depth.push(depth.last().unwrap() + d);
        }
        let iv = s.trace_intervals(&Trace { xyz: xyz.clone(), depth: depth.clone() }).unwrap();
        prop_assert_eq!(iv[0].from, 0.0);
        prop_assert_eq!(iv.last().unwrap().to, *depth.last().unwrap());
        for w in iv.windows(2) {
            prop_assert_eq!(w[0].to, w[1].from);
            prop_assert!(w[0].inside != w[1].inside);
            prop_assert!(w[0].exit_face.is_some() && w[1].entry_face.is_some());
        }
        prop_assert_eq!(iv[0].inside, s.locate(xyz[0]) != Location::Outside);
    }

    #[test]
    fn ray_parity_from_outside_matches_locate(seed in 0u64..1000) {
        let m = fixtures::torus(12, 8);
        let surf = SurfaceIndex::new(&m);
        let s = SolidIndex::new(&m).unwrap();
        let ends = fixtures::random_points(seed, 30, [-1.5; 3], [1.5; 3]);
        let origins = fixtures::random_points(seed + 1, 30, [3.0; 3], [4.0; 3]);
        let hits = surf.segment_hits(&origins, &ends).unwrap();
        let dirs: Vec<[f64; 3]> = origins.iter().zip(&ends).map(|(o, e)| [0, 1, 2].map(|k| e[k] - o[k])).collect();
        let rays = surf.ray_hits(&origins, &dirs).unwrap();
        for ((h, r), e) in hits.iter().zip(&rays).zip(&ends) {
            prop_assume!(h.iter().all(|x| x.kind == HitKind::Interior));
            prop_assert_eq!(h.len() % 2 == 1, s.locate(*e) == Location::Inside);
            prop_assert!(h.windows(2).all(|w| w[0].entering != w[1].entering));
            prop_assert!(r.len() % 2 == 0);
            prop_assert!(r.windows(2).all(|w| w[0].t <= w[1].t));
        }
    }

    #[test]
    fn locate_agrees_with_the_winding_number(seed in 0u64..1000) {
        let m = fixtures::torus(10, 6);
        let s = SolidIndex::new(&m).unwrap();
        for p in fixtures::random_points(seed, 40, [-1.5; 3], [1.5; 3]) {
            prop_assume!(s.signed_distance(&[p])[0].abs() > 1e-6);
            let want = if winding(&m, p).round() != 0.0 { Location::Inside } else { Location::Outside };
            prop_assert_eq!(s.locate(p), want);
        }
    }
}

#[test]
fn block_fraction_is_bitwise_equal_on_one_and_four_threads() {
    let s = SolidIndex::new(&fixtures::torus(16, 10)).unwrap();
    let blocks = grid(
        [0.01, 0.02, 0.03],
        [0.1, 0.12, 0.08],
        1.6,
        Some(rotation(0.3, 0.2, 0.1)),
    );
    let run = |n| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| s.block_fraction(&blocks).unwrap())
    };
    let (a, b) = (run(1), run(4));
    assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
}

#[test]
fn utm_scale_fractions_match_the_origin() {
    let m = fixtures::sphere(2);
    let off = [700_000.0, 7_500_000.0, 1_000.0];
    let far = fixtures::translated(&m, off);
    let bf = grid([0.05, 0.05, 0.05], [0.25; 3], 1.2, None);
    let mut bn = bf.clone();
    bn.centroids = bf
        .centroids
        .iter()
        .map(|c| [0, 1, 2].map(|k| c[k] + off[k]))
        .collect();
    let a = SolidIndex::new(&m).unwrap().block_fraction(&bf).unwrap();
    let z = SolidIndex::new(&far).unwrap().block_fraction(&bn).unwrap();
    let total = |f: &[f64]| f.iter().sum::<f64>() * 0.25f64.powi(3);
    assert!(
        close(total(&a), total(&z), 1e-9),
        "{} vs {}",
        total(&a),
        total(&z)
    );
    assert!(a.iter().zip(&z).all(|(x, y)| (x - y).abs() < 1e-6));
}

fn small_solids() -> Vec<Mesh> {
    vec![
        fixtures::cube(),
        fixtures::tetrahedron(),
        fixtures::sphere(1),
        fixtures::torus(10, 6),
        fixtures::nested_shells(),
    ]
}

/// Points at random heights that share x, y or both with mesh vertices.
fn vertex_plan_points(m: &Mesh, rng: &mut fixtures::Rng, n: usize) -> Vec<[f64; 3]> {
    let verts = m.vertices();
    let mut pts = vec![];
    for _ in 0..n {
        let v = verts[rng.below(verts.len())];
        let w = verts[rng.below(verts.len())];
        let z = rng.uniform(-1.6, 1.6);
        let u = rng.uniform(-1.6, 1.6);
        pts.extend([
            [v[0], v[1], z],
            [v[0], w[1], z],
            [w[0], v[1], z],
            [v[0], u, z],
            [u, v[1], z],
            [v[0], v[1], w[2]],
        ]);
    }
    pts
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn locate_agrees_with_the_exact_winding_number(seed in 0u64..10_000) {
        for m in small_solids() {
            let s = SolidIndex::new(&m).unwrap();
            let mut rng = fixtures::Rng::new(seed);
            let mut pts = vertex_plan_points(&m, &mut rng, 8);
            pts.extend(fixtures::random_points(seed, 16, [-1.6; 3], [1.6; 3]));
            for p in pts {
                if (0..m.n_faces() as u32).any(|f| point_in_triangle(p, m.corners(f))) {
                    prop_assert_eq!(s.locate(p), Location::Boundary);
                    continue;
                }
                let w = (0..32)
                    .find_map(|_| {
                        let d = [0; 3].map(|_| rng.uniform(-1.0, 1.0));
                        exact_winding(&m, p, [0, 1, 2].map(|k| p[k] + 100.0 * d[k]))
                    })
                    .expect("a generic direction exists");
                let want = if w != 0 { Location::Inside } else { Location::Outside };
                prop_assert_eq!(s.locate(p), want, "p = {:?}, winding {}", p, w);
            }
        }
    }

    #[test]
    fn locate_matches_winding_at_vertex_plan_coordinates(seed in 0u64..10_000) {
        for m in [fixtures::torus(10, 6), fixtures::sphere(1), fixtures::nested_shells()] {
            let s = SolidIndex::new(&m).unwrap();
            let mut rng = fixtures::Rng::new(seed);
            let verts = m.vertices();
            for _ in 0..20 {
                let v = verts[rng.below(verts.len())];
                let w = verts[rng.below(verts.len())];
                let z = rng.uniform(-1.6, 1.6);
                for p in [[v[0], v[1], z], [v[0], w[1], z], [w[0], v[1], z]] {
                    let d = s.signed_distance(&[p])[0];
                    if d.abs() < 1e-6 {
                        continue;
                    }
                    let want = if winding(&m, p).round() != 0.0 { Location::Inside } else { Location::Outside };
                    prop_assert_eq!(s.locate(p), want, "p = {:?}", p);
                }
            }
        }
    }

    #[test]
    fn trace_labels_match_locate_on_tiny_pieces(seed in 0u64..10_000, k in 1u32..64) {
        let m = fixtures::cube();
        let s = SolidIndex::new(&m).unwrap();
        let mut rng = fixtures::Rng::new(seed);
        let a = [rng.uniform(0.1, 0.9), rng.uniform(0.1, 0.9), rng.uniform(0.1, 0.9)];
        let axis = rng.below(3);
        let mut b = a;
        b[axis] = if rng.below(2) == 0 { 1.0 + k as f64 * f64::EPSILON } else { -(k as f64) * f64::EPSILON * 0.5 };
        let iv = s.trace_intervals(&Trace { xyz: vec![a, b], depth: vec![0.0, 1.0] }).unwrap();
        prop_assert_eq!(s.locate(b), Location::Outside);
        prop_assert!(iv.first().unwrap().inside);
        prop_assert!(!iv.last().unwrap().inside, "{:?}", iv);
        prop_assert_eq!(iv.len(), 2);
        prop_assert!(iv[0].exit_face.is_some() && iv[1].entry_face.is_some());
    }
}

#[test]
fn trace_leaving_by_one_ulp_has_an_outside_tail() {
    let s = SolidIndex::new(&fixtures::cube()).unwrap();
    let q = [0.5, 0.5, 1.0 + f64::EPSILON];
    assert_eq!(s.locate(q), Location::Outside);
    let iv = s
        .trace_intervals(&Trace {
            xyz: vec![[0.5, 0.5, 0.5], q],
            depth: vec![0.0, 1.0],
        })
        .unwrap();
    assert!(iv.first().unwrap().inside);
    assert!(
        !iv.last().unwrap().inside,
        "the trace ends outside the solid but its last interval is inside: {iv:?}"
    );
}

#[test]
fn closest_tie_on_a_shared_edge_takes_the_smaller_face() {
    let s = SurfaceIndex::new(&fixtures::cube());
    let mut rng = fixtures::Rng::new(7);
    let mut bad = vec![];
    for _ in 0..2000 {
        let x = rng.uniform(0.05, 0.95);
        let p = [x, -1e-20, -1e-20];
        let c = s.closest(&[p])[0];
        if c.face != 1 {
            bad.push((x, c.face));
        }
    }
    assert!(
        bad.is_empty(),
        "{} of 2000 points equidistant from faces 1 and 4 chose another face, first {:?}",
        bad.len(),
        bad.first()
    );
}

#[test]
fn ray_through_two_cube_corners_reports_two_vertex_hits() {
    let s = SurfaceIndex::new(&fixtures::cube());
    let h = s
        .ray_hits(&[[-1.0, -1.0, -1.0]], &[[1.0, 1.0, 1.0]])
        .unwrap();
    let pts: Vec<_> = h[0].iter().map(|x| (x.t, x.point, x.kind)).collect();
    assert_eq!(h[0].len(), 2, "{pts:?}");
    assert_eq!(h[0][0].point, [0.0; 3]);
    assert_eq!(h[0][1].point, [1.0; 3]);
    assert!(h[0].iter().all(|x| x.kind == HitKind::Vertex));
}

#[test]
fn block_equal_to_the_solid_and_its_neighbour() {
    let s = SolidIndex::new(&fixtures::cube()).unwrap();
    let f = s
        .block_fraction(&Blocks {
            centroids: vec![[0.5; 3], [1.5, 0.5, 0.5], [0.5, 0.5, 1.5], [0.5, 0.5, -0.5]],
            size: [1.0; 3],
            rotation: None,
        })
        .unwrap();
    assert_eq!(f, [1.0, 0.0, 0.0, 0.0]);
}

#[test]
fn utm_trace_chord() {
    let off = [700_000.0, 7_500_000.0, 1_000.0];
    let m = fixtures::translated(&fixtures::cube(), off);
    let s = SolidIndex::new(&m).unwrap();
    let a = [off[0] + 0.3, off[1] + 0.4, off[2] + 2.0];
    let b = [off[0] + 0.6, off[1] + 0.7, off[2] - 1.0];
    let len = ((0.3f64).powi(2) * 2.0 + 9.0).sqrt();
    let iv = s
        .trace_intervals(&Trace {
            xyz: vec![a, b],
            depth: vec![0.0, len],
        })
        .unwrap();
    assert_eq!(iv.len(), 3, "{iv:?}");
    let inside = iv[1].to - iv[1].from;
    assert!(
        (inside - len / 3.0).abs() < 1e-9,
        "{inside} vs {}",
        len / 3.0
    );
}

#[test]
fn signed_distance_sign_at_utm_scale() {
    let off = [700_000.0, 7_500_000.0, 1_000.0];
    let m = fixtures::translated(&fixtures::sphere(2), off);
    let s = SolidIndex::new(&m).unwrap();
    let pts = fixtures::random_points(3, 200, [-1.2; 3], [1.2; 3]);
    let far: Vec<_> = pts
        .iter()
        .map(|p| [0, 1, 2].map(|k| p[k] + off[k]))
        .collect();
    let near = SolidIndex::new(&fixtures::sphere(2)).unwrap();
    let a = near.signed_distance(&pts);
    let b = s.signed_distance(&far);
    for (x, y) in a.iter().zip(&b) {
        if x.abs() > 1e-6 {
            assert_eq!(x.signum(), y.signum());
            assert!((x - y).abs() < 1e-6, "{x} vs {y}");
        }
    }
}

#[test]
fn tiny_positive_block_size_gives_a_fraction_in_range() {
    let s = SolidIndex::new(&fixtures::cube()).unwrap();
    let f = s
        .block_fraction(&Blocks {
            centroids: vec![[0.5; 3]],
            size: [1e-110; 3],
            rotation: None,
        })
        .unwrap();
    assert!((0.0..=1.0).contains(&f[0]), "fraction is {}", f[0]);
}

#[test]
fn empty_mesh_is_not_a_solid() {
    assert!(SolidIndex::new(&Mesh::empty()).is_err());
}
