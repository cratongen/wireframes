use wireframes_core::{Aabb, Error, Mesh, Result};

use crate::measure::{check_closed_oriented, sub, sum_by};

/// Sutherland-Hodgman keeping the part where `keep(p[axis])`; `at` is the cut.
fn clip(poly: &[[f64; 3]], axis: usize, at: f64, keep: impl Fn(f64) -> bool) -> Vec<[f64; 3]> {
    let inside = |p: &[f64; 3]| keep(p[axis]);
    let mut out = Vec::with_capacity(poly.len() + 2);
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        if inside(&a) {
            out.push(a);
        }
        if inside(&a) != inside(&b) {
            let t = (at - a[axis]) / (b[axis] - a[axis]);
            let mut p: [f64; 3] = std::array::from_fn(|k| a[k] + t * (b[k] - a[k]));
            p[axis] = at;
            out.push(p);
        }
    }
    out
}

/// `sum over the fan of (signed plan area) * f(mean z)`, twice.
fn fan(poly: &[[f64; 3]], z: impl Fn(&[[f64; 3]; 3]) -> f64) -> f64 {
    (1..poly.len().saturating_sub(1))
        .map(|i| {
            let t = [poly[0], poly[i], poly[i + 1]];
            let a2 = (t[1][0] - t[0][0]) * (t[2][1] - t[0][1])
                - (t[1][1] - t[0][1]) * (t[2][0] - t[0][0]);
            a2 * z(&t)
        })
        .sum()
}

/// Volume of the solid inside an axis-aligned box, by clamped-column integration.
pub fn box_volume(mesh: &Mesh, faces: &[u32], b: &Aabb) -> f64 {
    let d = sub(b.max, b.min);
    let s = sum_by(faces.len(), |k| {
        let mut poly: Vec<[f64; 3]> = mesh.corners(faces[k]).map(|p| sub(p, b.min)).to_vec();
        poly = clip(&poly, 0, 0.0, |x| x >= 0.0);
        poly = clip(&poly, 0, d[0], |x| x <= d[0]);
        poly = clip(&poly, 1, 0.0, |y| y >= 0.0);
        poly = clip(&poly, 1, d[1], |y| y <= d[1]);
        let poly = clip(&poly, 2, 0.0, |z| z >= 0.0);
        let mid = clip(&poly, 2, d[2], |z| z <= d[2]);
        // Strict, so a piece lying in the plane `z = zmax` counts once.
        let top = clip(&poly, 2, d[2], |z| z > d[2]);
        let inside = fan(&mid, |t| (t[0][2] + t[1][2] + t[2][2]) / 3.0);
        let above = fan(&top, |_| d[2]);
        [0.5 * (inside + above)]
    });
    s[0]
}

/// Solid volume between consecutive `edges` (strictly increasing, at least 2).
pub fn volume_by_elevation(mesh: &Mesh, edges: &[f64]) -> Result<Vec<f64>> {
    if edges.len() < 2 {
        return Err(Error::invalid(format!(
            "edges needs at least 2 values, got {}",
            edges.len()
        )));
    }
    if let Some(i) = (0..edges.len()).find(|&i| !edges[i].is_finite()) {
        return Err(Error::invalid(format!(
            "edges must be finite: edge {i} is {}",
            edges[i]
        )));
    }
    if let Some(i) = (1..edges.len()).find(|&i| edges[i] <= edges[i - 1]) {
        return Err(Error::invalid(format!(
            "edges must be strictly increasing: edge {i} is {} after {}",
            edges[i],
            edges[i - 1]
        )));
    }
    check_closed_oriented(mesh, "volume_by_elevation")?;
    let bounds = mesh.bounds().expect("a closed mesh has vertices");
    let faces: Vec<u32> = (0..mesh.n_faces() as u32).collect();
    Ok(edges
        .windows(2)
        .map(|w| {
            let b = Aabb {
                min: [bounds.min[0], bounds.min[1], w[0]],
                max: [bounds.max[0], bounds.max[1], w[1]],
            };
            box_volume(mesh, &faces, &b)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::measure::volume;
    use proptest::prelude::*;
    use wireframes_core::fixtures;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * 1f64.max(a.abs()).max(b.abs())
    }

    fn all(m: &Mesh) -> Vec<u32> {
        (0..m.n_faces() as u32).collect()
    }

    #[test]
    fn box_volume_of_a_cube_part() {
        let m = fixtures::cube();
        let b = Aabb {
            min: [0.5, -1.0, 0.25],
            max: [2.0, 0.5, 0.75],
        };
        assert!(close(box_volume(&m, &all(&m), &b), 0.5 * 0.5 * 0.5));
        let t = fixtures::tetrahedron();
        let b = Aabb {
            min: [-1.0; 3],
            max: [2.0, 2.0, 0.5],
        };
        assert!(close(box_volume(&t, &all(&t), &b), 1.0 / 6.0 - 1.0 / 48.0));
    }

    #[test]
    fn faces_in_the_box_top_count_once() {
        let m = fixtures::cube();
        let b = Aabb {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!(close(box_volume(&m, &all(&m), &b), 1.0));
    }

    #[test]
    fn elevation_slabs_sum_to_the_volume() {
        let s = fixtures::translated(&fixtures::sphere(3), [500000.0, 7500000.0, 300.0]);
        let edges: Vec<f64> = (0..=9).map(|i| 298.5 + 0.3125 * i as f64).collect();
        let slabs = volume_by_elevation(&s, &edges).unwrap();
        let total: f64 = slabs.iter().sum();
        assert!(close(total, volume(&s).unwrap()), "{total}");
        assert!(slabs[0] == 0.0 && slabs[1] > 0.0);
        assert!(volume_by_elevation(&s, &[1.0]).is_err());
        assert!(volume_by_elevation(&s, &[1.0, 1.0]).is_err());
        assert!(volume_by_elevation(&fixtures::open_cube(), &[0.0, 1.0]).is_err());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]
        #[test]
        fn box_grid_sums_to_volume(n in 1usize..5, sub in 0u32..3) {
            let s = fixtures::sphere(sub);
            let b = s.bounds().unwrap();
            let step: [f64; 3] = std::array::from_fn(|k| (b.max[k] - b.min[k]) / n as f64);
            let mut total = 0.0;
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let lo = [
                            b.min[0] + step[0] * i as f64,
                            b.min[1] + step[1] * j as f64,
                            b.min[2] + step[2] * k as f64,
                        ];
                        let hi = [lo[0] + step[0], lo[1] + step[1], lo[2] + step[2]];
                        total += box_volume(&s, &all(&s), &Aabb { min: lo, max: hi });
                    }
                }
            }
            let v = volume(&s).unwrap();
            prop_assert!((total - v).abs() <= 1e-9 * v, "{} {}", total, v);
        }
    }
}
