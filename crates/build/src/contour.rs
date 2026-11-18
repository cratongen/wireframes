use drusa_boolean::slice_by_plane;
use drusa_core::{Error, Mesh, Polylines, Result};

/// `slice_by_plane` through the origin with normal `+z` at offsets `levels`,
/// the `offset` column renamed `elevation`.
pub fn contours(surface: &Mesh, levels: &[f64]) -> Result<Polylines> {
    if surface.face_attributes().column("elevation").is_some() {
        return Err(Error::invalid(
            "surface must not have a face column named \"elevation\"",
        ));
    }
    let lines = slice_by_plane(surface, [0.0; 3], [0.0, 0.0, 1.0], levels)?;
    let attrs = lines.attributes().clone();
    let elevation = attrs
        .column("offset")
        .expect("slice_by_plane writes offset")
        .clone();
    let attrs = attrs
        .without_column("offset")
        .with_column("elevation", elevation)?;
    lines.with_attributes(attrs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use drusa_core::Column;
    use drusa_core::fixtures::{box_mesh, sphere};
    use proptest::prelude::*;

    fn plane(n: u32, a: f64, b: f64, c: f64) -> Mesh {
        let mut v = vec![];
        for j in 0..n {
            for i in 0..n {
                let (x, y) = (i as f64 * 1.37, j as f64 * 0.91);
                v.push([x, y, a * x + b * y + c]);
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

    proptest! {
        #[test]
        fn plane_contours_are_straight_and_on_the_plane(n in 2u32..10, a in -2.0f64..2.0, b in -2.0f64..2.0, c in -5.0f64..5.0, t in 0.05f64..0.95) {
            prop_assume!(a.abs() + b.abs() > 0.1);
            let m = plane(n, a, b, c);
            let zs: Vec<f64> = m.vertices().iter().map(|p| p[2]).collect();
            let (lo, hi) = zs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &z| (l.min(z), h.max(z)));
            let level = lo + t * (hi - lo);
            let p = contours(&m, &[level]).unwrap();
            prop_assert_eq!(p.len(), 1);
            let pts = p.part(0);
            let (s, e) = (pts[0], pts[pts.len() - 1]);
            let scale = 1.0 + (a.abs() + b.abs()) * 15.0 + c.abs();
            for q in pts {
                prop_assert_eq!(q[2], level);
                prop_assert!((a * q[0] + b * q[1] + c - level).abs() <= 1e-9 * scale);
                let cr = (e[0] - s[0]) * (q[1] - s[1]) - (e[1] - s[1]) * (q[0] - s[0]);
                prop_assert!(cr.abs() <= 1e-9 * scale * scale);
            }
        }
    }

    #[test]
    fn solid_slices_are_counter_clockwise_rings() {
        let m = box_mesh([0.0; 3], [2.0, 3.0, 4.0]);
        let p = contours(&m, &[1.0, 2.5]).unwrap();
        assert_eq!(p.len(), 2);
        assert_eq!(p.closed(), &[true, true]);
        for f in 0..2 {
            let r = p.part(f);
            let shoelace: f64 = (0..r.len())
                .map(|i| {
                    let (a, b) = (r[i], r[(i + 1) % r.len()]);
                    a[0] * b[1] - b[0] * a[1]
                })
                .sum::<f64>()
                / 2.0;
            assert!((shoelace - 6.0).abs() < 1e-12, "{shoelace}");
        }
        let e = p.attributes().require("elevation").unwrap();
        assert_eq!(e, &Column::f64(vec![1.0, 2.5]));
        let s = contours(&sphere(2), &[0.3]).unwrap();
        assert_eq!(s.len(), 1);
        assert!(s.closed()[0]);
    }
}
