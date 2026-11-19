use proptest::prelude::*;
use wireframes_core::fixtures::{self, Rng};
use wireframes_core::{Aabb, Mesh};
use wireframes_ops::{
    HoleFill, RepairOptions, RepairReport, ThicknessMode, area, box_volume, centroid, components,
    flip, repair, thickness, translate, volume, volume_by_elevation,
};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * 1f64.max(a.abs()).max(b.abs())
}

fn solid(which: usize) -> Mesh {
    match which {
        0 => fixtures::cube(),
        1 => fixtures::sphere(1),
        2 => fixtures::torus(8, 6),
        _ => fixtures::nested_shells(),
    }
}

fn plane(z0: f64, gx: f64, gy: f64) -> Mesh {
    let v = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
        .map(|[x, y]| [x, y, z0 + gx * x + gy * y])
        .to_vec();
    Mesh::new(v, vec![[0, 1, 2], [0, 2, 3]]).unwrap()
}

/// Shuffled faces and vertices, a few duplicated faces, about a third flipped,
/// then exploded into a soup.
fn scrambled(seed: u64, base: &Mesh) -> Mesh {
    let mut rng = Rng::new(seed);
    let mut faces = base.triangles().to_vec();
    for _ in 0..rng.below(4) {
        let f = faces[rng.below(faces.len())];
        faces.push(f);
    }
    rng.shuffle(&mut faces);
    for t in faces.iter_mut() {
        if rng.below(3) == 0 {
            *t = [t[0], t[2], t[1]];
        }
    }
    let mut perm: Vec<u32> = (0..base.n_vertices() as u32).collect();
    rng.shuffle(&mut perm);
    let mut v = vec![[0.0; 3]; base.n_vertices()];
    for (i, &p) in perm.iter().enumerate() {
        v[p as usize] = base.vertices()[i];
    }
    let t = faces.iter().map(|t| t.map(|i| perm[i as usize])).collect();
    fixtures::soup(&Mesh::new(v, t).unwrap())
}

fn check_maps(input: &Mesh, out: &Mesh, r: &RepairReport) -> std::result::Result<(), String> {
    if r.vertex_map.len() != input.n_vertices() || r.face_map.len() != input.n_faces() {
        return Err("map lengths differ from the input".into());
    }
    for (f, m) in r.face_map.iter().enumerate() {
        let Some(g) = m else { continue };
        let mut want =
            input.triangles()[f].map(|i| r.vertex_map[i as usize].expect("corner of a kept face"));
        let mut got = out.triangles()[*g as usize];
        want.sort_unstable();
        got.sort_unstable();
        if want != got {
            return Err(format!("face {f}: {want:?} != {got:?}"));
        }
    }
    Ok(())
}

fn options() -> Vec<RepairOptions> {
    let mut all = vec![RepairOptions::default()];
    for how in [HoleFill::Fan, HoleFill::MinArea] {
        all.push(RepairOptions {
            fill_holes: Some(how),
            ..RepairOptions::default()
        });
    }
    all
}

fn open_surface(which: usize, seed: u64) -> Mesh {
    match which {
        0 => fixtures::open_cube(),
        1 => fixtures::flipped(&fixtures::open_cube()),
        2 => fixtures::pit_surface(6, seed),
        3 => plane(1.0, 0.5, -0.25),
        _ => Mesh::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap(),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn volume_is_translation_invariant_and_flips_sign(
        which in 0usize..4,
        d in proptest::array::uniform3(-1e5..1e5f64),
    ) {
        let m = solid(which);
        let v = volume(&m).unwrap();
        prop_assert!(close(volume(&translate(&m, d).unwrap()).unwrap(), v));
        prop_assert_eq!(volume(&flip(&m)).unwrap(), -v);
    }

    #[test]
    fn box_measures_are_analytic(
        lo in proptest::array::uniform3(-1e3..1e3f64),
        size in proptest::array::uniform3(0.1..50.0f64),
    ) {
        let hi = [lo[0] + size[0], lo[1] + size[1], lo[2] + size[2]];
        let m = fixtures::box_mesh(lo, hi);
        let d = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        prop_assert!(close(volume(&m).unwrap(), d[0] * d[1] * d[2]));
        prop_assert!(close(area(&m), 2.0 * (d[0] * d[1] + d[1] * d[2] + d[0] * d[2])));
        let c = centroid(&m).unwrap();
        for k in 0..3 {
            let want = 0.5 * (lo[k] + hi[k]);
            prop_assert!((c[k] - want).abs() <= 1e-9 * 1f64.max(want.abs()).max(d[k]), "{:?}", c);
        }
    }

    #[test]
    fn elevation_slabs_sum_to_the_volume(which in 0usize..4, cuts in 1usize..8, seed in any::<u64>()) {
        let m = solid(which);
        let b = m.bounds().unwrap();
        let mut rng = Rng::new(seed);
        let mut edges: Vec<f64> = (0..cuts).map(|_| rng.uniform(b.min[2], b.max[2])).collect();
        edges.push(b.min[2] - 1.0);
        edges.push(b.max[2] + 1.0);
        edges.sort_by(f64::total_cmp);
        edges.dedup();
        let slabs = volume_by_elevation(&m, &edges).unwrap();
        prop_assert_eq!(slabs.len(), edges.len() - 1);
        prop_assert!(close(slabs.iter().sum(), volume(&m).unwrap()), "{:?}", slabs);
    }

    #[test]
    fn box_volumes_over_a_covering_grid_sum_to_the_volume(
        which in 0usize..4,
        n in proptest::array::uniform3(1usize..4),
        shift in proptest::array::uniform3(0.0..0.5f64),
    ) {
        let m = solid(which);
        let b = m.bounds().unwrap();
        let lo: [f64; 3] = std::array::from_fn(|k| b.min[k] - shift[k]);
        let step: [f64; 3] = std::array::from_fn(|k| (b.max[k] + 1.0 - lo[k]) / n[k] as f64);
        let faces: Vec<u32> = (0..m.n_faces() as u32).collect();
        let mut total = 0.0;
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let at = [i, j, k];
                    let min: [f64; 3] = std::array::from_fn(|a| lo[a] + at[a] as f64 * step[a]);
                    let max: [f64; 3] = std::array::from_fn(|a| lo[a] + (at[a] + 1) as f64 * step[a]);
                    total += box_volume(&m, &faces, &Aabb { min, max });
                }
            }
        }
        prop_assert!(close(total, volume(&m).unwrap()), "{} {}", total, volume(&m).unwrap());
    }

    #[test]
    fn repair_restores_a_scrambled_solid_and_is_idempotent(which in 0usize..4, seed in any::<u64>()) {
        let base = solid(which);
        let m = scrambled(seed, &base);
        let (out, r) = repair(&m, &RepairOptions::default()).unwrap();
        prop_assert!(r.after.is_solid(), "{:?}", r.after);
        prop_assert!(close(volume(&out).unwrap(), volume(&base).unwrap()));
        check_maps(&m, &out, &r).map_err(TestCaseError::fail)?;
        let (again, r2) = repair(&out, &RepairOptions::default()).unwrap();
        prop_assert!(r2.actions.is_empty(), "{:?}", r2.actions);
        prop_assert_eq!(again, out);
    }

    #[test]
    fn repair_of_open_surfaces_is_idempotent(which in 0usize..5, seed in any::<u64>(), opt in 0usize..3) {
        let opts = &options()[opt];
        let m = open_surface(which, seed);
        let (out, r) = repair(&m, opts).unwrap();
        check_maps(&m, &out, &r).map_err(TestCaseError::fail)?;
        prop_assert_eq!(r.after.degenerate_faces, 0);
        let (again, r2) = repair(&out, opts).unwrap();
        prop_assert!(r2.actions.is_empty(), "{:?}", r2.actions);
        prop_assert_eq!(again, out);
    }

    #[test]
    fn components_partition_faces_and_volume(
        lows in proptest::collection::vec(proptest::array::uniform3(-10.0..10.0f64), 1..5),
        which in 0usize..4,
    ) {
        let mut parts: Vec<Mesh> = lows
            .iter()
            .map(|&lo| fixtures::box_mesh(lo, [lo[0] + 1.0, lo[1] + 2.0, lo[2] + 0.5]))
            .collect();
        parts.push(fixtures::translated(&solid(which), [30.0, 0.0, 0.0]));
        let m = Mesh::concat(&parts.iter().collect::<Vec<_>>()).unwrap();
        let cs = components(&m).unwrap();
        prop_assert_eq!(cs.iter().map(Mesh::n_faces).sum::<usize>(), m.n_faces());
        let sum: f64 = cs.iter().map(|c| volume(c).unwrap()).sum();
        prop_assert!(close(sum, volume(&m).unwrap()), "{} {}", sum, volume(&m).unwrap());
    }

    #[test]
    fn vertical_thickness_of_two_planes_is_analytic(
        a in -50.0..50.0f64, b in -50.0..50.0f64,
        g in proptest::array::uniform4(-1.0..1.0f64),
        p in proptest::array::uniform2(0.0..10.0f64),
    ) {
        let t = thickness(&plane(a, g[0], g[1]), &plane(b, g[2], g[3]), &[p], ThicknessMode::Vertical)
            .unwrap()[0];
        let want = (a + g[0] * p[0] + g[1] * p[1]) - (b + g[2] * p[0] + g[3] * p[1]);
        prop_assert!(close(t, want), "{} {}", t, want);
    }
}

#[test]
fn touching_shells_measure_and_split_into_components() {
    use wireframes_ops::{centroid, components, volume, volume_by_elevation};
    for (name, m, c) in [
        (
            "edge, apart",
            fixtures::edge_touching_cubes(false),
            [1.0, 1.0, 0.5],
        ),
        (
            "edge, welded",
            fixtures::edge_touching_cubes(true),
            [1.0, 1.0, 0.5],
        ),
        ("vertex", fixtures::vertex_touching_cubes(), [1.0, 1.0, 1.0]),
    ] {
        let v = volume(&m).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!((v - 2.0).abs() < 1e-12, "{name}: {v}");
        let g = centroid(&m).unwrap();
        assert!((0..3).all(|k| (g[k] - c[k]).abs() < 1e-12), "{name}: {g:?}");
        let slabs = volume_by_elevation(&m, &[-1.0, 0.5, 1.0, 3.0]).unwrap();
        assert!(
            (slabs.iter().sum::<f64>() - 2.0).abs() < 1e-12,
            "{name}: {slabs:?}"
        );
        let parts = components(&m).unwrap();
        assert_eq!(parts.len(), 2, "{name}");
        for p in &parts {
            assert!((volume(p).unwrap() - 1.0).abs() < 1e-12, "{name}");
            assert!(p.validate(0.0).unwrap().summary.is_solid(), "{name}");
        }
    }
}
