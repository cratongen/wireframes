use proptest::prelude::*;
use wireframes_core::fixtures::{self, broken};
use wireframes_core::{Mesh, MeshProblemKind};

fn rotate(p: [f64; 3], a: [f64; 3]) -> [f64; 3] {
    let (s, c) = (a.map(f64::sin), a.map(f64::cos));
    let x = [p[0], c[0] * p[1] - s[0] * p[2], s[0] * p[1] + c[0] * p[2]];
    let y = [c[1] * x[0] + s[1] * x[2], x[1], -s[1] * x[0] + c[1] * x[2]];
    [c[2] * y[0] - s[2] * y[1], s[2] * y[0] + c[2] * y[1], y[2]]
}

#[test]
fn solids_report_is_solid_and_lose_it_with_one_face() {
    for m in [
        fixtures::cube(),
        fixtures::sphere(3),
        fixtures::torus(16, 9),
    ] {
        let r = m.validate(0.0).unwrap();
        assert!(r.summary.is_solid(), "{:?}", r.summary);
        let keep: Vec<bool> = (0..m.n_faces()).map(|f| f != 5).collect();
        let s = m
            .filter_faces(&keep)
            .unwrap()
            .validate(0.0)
            .unwrap()
            .summary;
        assert_eq!(s.boundary_edges, 3);
    }
    let s = broken::inward_cube().validate(0.0).unwrap().summary;
    assert_eq!(s.inward_shells, 1);
    assert!(s.is_closed && !s.is_solid());
}

proptest! {
    #[test]
    fn mesh_new_rejects_bad_indices_and_coordinates(n in 1usize..20, bad in 0u32..3, value in prop_oneof![Just(f64::NAN), Just(f64::INFINITY), Just(f64::NEG_INFINITY)]) {
        let v = vec![[0.0; 3]; n];
        let t = vec![[0, 0, n as u32 + bad]];
        prop_assert!(Mesh::new(v.clone(), t).is_err());
        let mut w = v.clone();
        w[n - 1][bad as usize] = value;
        prop_assert!(Mesh::new(w, vec![]).is_err());
        prop_assert!(Mesh::new(v, vec![[0, 0, n as u32 - 1]]).is_ok());
    }

    /// A rotated, shifted pair of unit cubes crosses only while they overlap.
    #[test]
    fn crossings_follow_overlap(
        angles in proptest::array::uniform3(-3.0..3.0f64),
        offset in proptest::array::uniform3(-1e3..1e3f64),
        shift in 0.1..1.9f64,
    ) {
        let pair = Mesh::concat(&[&fixtures::cube(), &fixtures::translated(&fixtures::cube(), [shift, 0.3, 0.2])]).unwrap();
        let moved: Vec<[f64; 3]> = pair
            .vertices()
            .iter()
            .map(|&p| rotate(p, angles))
            .map(|p| [0, 1, 2].map(|k| p[k] + offset[k]))
            .collect();
        let pair = pair.with_vertices(moved).unwrap();
        let s = pair.validate(0.0).unwrap().summary;
        prop_assert_eq!(s.self_intersections > 0, shift < 1.0);
        let keep: Vec<bool> = (0..24).map(|f| f < 12).collect();
        let one = pair.filter_faces(&keep).unwrap().validate(0.0).unwrap();
        prop_assert!(one.problems.iter().all(|p| p.kind == MeshProblemKind::UnusedVertex));
    }
}
