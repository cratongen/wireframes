//! Open-pit acceptance cases reduced to their smallest failing inputs.

mod support;
use support::closed_manifold;

use wireframes_boolean::{BooleanOptions, Op, boolean};
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

/// Box with a flat bottom at `bottom` and its top corners on `z = 5 + x / 4 + y / 2`.
fn prism(half: f64, bottom: f64) -> Mesh {
    let mut v = Vec::new();
    for top in [false, true] {
        for y in [-half, half] {
            for x in [-half, half] {
                v.push([x, y, if top { 5.0 + x / 4.0 + y / 2.0 } else { bottom }]);
            }
        }
    }
    Mesh::new(v, BOX_FACES.to_vec()).unwrap()
}

#[test]
#[ignore = "acceptance: rounding folds a sliver of the intersection onto a wall"]
fn booleans_of_prisms_whose_caps_differ_by_one_ulp_are_valid_solids() {
    let a = prism(10.0, -10.0);
    let exact = prism(6.0, -5.0);
    let mut v = exact.vertices().to_vec();
    v[4][2] = v[4][2].next_up();
    let b = Mesh::new(v, exact.triangles().to_vec()).unwrap();
    for m in [&a, &b] {
        assert!(m.validate(0.0).unwrap().summary.is_solid());
    }
    for op in [Op::Union, Op::Intersection, Op::Difference] {
        for (x, y) in [(&a, &b), (&b, &a)] {
            let m = boolean(x, y, op, &BooleanOptions::default()).unwrap().mesh;
            assert!(m.is_empty() || closed_manifold(&m), "{op:?}");
            let s = m.validate(0.0).unwrap().summary;
            assert!(m.is_empty() || s.is_solid(), "{op:?}: {s:?}");
        }
    }
}
