use drusa_core::predicates::orient2d;
use drusa_core::{Aabb, Bvh, Mesh, Result, units};
use rayon::prelude::*;

use crate::measure::{dot, sub};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThicknessMode {
    Vertical,
    Normal,
}

fn plan_index(surface: &Mesh) -> Bvh {
    let boxes: Vec<Aabb> = surface
        .face_aabbs()
        .into_iter()
        .map(|b| Aabb {
            min: [b.min[0], b.min[1], 0.0],
            max: [b.max[0], b.max[1], 0.0],
        })
        .collect();
    Bvh::new(&boxes)
}

/// z of face `c` above `p` when its closed plan triangle covers `p`.
fn z_on(c: [[f64; 3]; 3], p: [f64; 2]) -> Option<f64> {
    let q = c.map(|v| [v[0], v[1]]);
    let whole = orient2d(q[0], q[1], q[2]);
    if whole == 0.0 {
        return None;
    }
    let w = [
        orient2d(q[1], q[2], p),
        orient2d(q[2], q[0], p),
        orient2d(q[0], q[1], p),
    ];
    if w.iter().any(|&x| x != 0.0 && (x > 0.0) != (whole > 0.0)) {
        return None;
    }
    let (a, b, cc) = (c[0], c[1], c[2]);
    let det = (b[0] - a[0]) * (cc[1] - a[1]) - (b[1] - a[1]) * (cc[0] - a[0]);
    let (dx, dy) = (p[0] - a[0], p[1] - a[1]);
    let s = (dx * (cc[1] - a[1]) - dy * (cc[0] - a[0])) / det;
    let t = ((b[0] - a[0]) * dy - (b[1] - a[1]) * dx) / det;
    Some(a[2] + s * (b[2] - a[2]) + t * (cc[2] - a[2]))
}

fn surface_z_with(surface: &Mesh, index: &Bvh, p: [f64; 2]) -> f64 {
    let probe = Aabb {
        min: [p[0], p[1], 0.0],
        max: [p[0], p[1], 0.0],
    };
    let mut best = f64::NAN;
    index.query(&probe, |f| {
        if let Some(z) = z_on(surface.corners(f), p)
            && (best.is_nan() || z > best)
        {
            best = z;
        }
    });
    best
}

/// Highest surface z at each `(x, y)`; NaN where no face covers it.
pub fn surface_z(surface: &Mesh, xy: &[[f64; 2]]) -> Vec<f64> {
    let index = plan_index(surface);
    xy.par_iter()
        .map(|&p| surface_z_with(surface, &index, p))
        .collect()
}

/// `z - surface_z` per point.
pub fn vertical_distance(surface: &Mesh, points: &[[f64; 3]]) -> Vec<f64> {
    let xy: Vec<[f64; 2]> = points.iter().map(|p| [p[0], p[1]]).collect();
    surface_z(surface, &xy)
        .into_iter()
        .zip(points)
        .map(|(z, p)| p[2] - z)
        .collect()
}

/// Closest point of triangle `t` to `p` (Ericson, Real-Time Collision Detection 5.1.5).
fn closest_on_triangle(p: [f64; 3], [a, b, c]: [[f64; 3]; 3]) -> [f64; 3] {
    let lerp = |u: [f64; 3], v: [f64; 3], t: f64| -> [f64; 3] {
        std::array::from_fn(|k| u[k] + t * (v[k] - u[k]))
    };
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return lerp(a, b, d1 / (d1 - d3));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return lerp(a, c, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return lerp(b, c, (d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = va + vb + vc;
    if denom == 0.0 {
        return a;
    }
    let (v, w) = (vb / denom, vc / denom);
    std::array::from_fn(|k| a[k] + ab[k] * v + ac[k] * w)
}

pub fn thickness(
    top: &Mesh,
    bottom: &Mesh,
    xy: &[[f64; 2]],
    mode: ThicknessMode,
) -> Result<Vec<f64>> {
    units::check_same(
        "top",
        (&top.crs, &top.length_unit),
        (&bottom.crs, &bottom.length_unit),
    )?;
    let zt = surface_z(top, xy);
    let zb = surface_z(bottom, xy);
    let vertical: Vec<f64> = zt.iter().zip(&zb).map(|(a, b)| a - b).collect();
    if mode == ThicknessMode::Vertical {
        return Ok(vertical);
    }
    let o = bottom.bounds().map_or([0.0; 3], |b| b.min);
    let local = bottom
        .clone()
        .with_vertices(bottom.vertices().iter().map(|&p| sub(p, o)).collect())?;
    let index = Bvh::new(&local.face_aabbs());
    Ok((0..xy.len())
        .into_par_iter()
        .map(|i| {
            if vertical[i].is_nan() {
                return f64::NAN;
            }
            let p = sub([xy[i][0], xy[i][1], zt[i]], o);
            let d2 = |f: u32| {
                let q = sub(closest_on_triangle(p, local.corners(f)), p);
                dot(q, q)
            };
            let d = index.nearest(p, d2).map_or(f64::NAN, |(_, d2)| d2.sqrt());
            if vertical[i] < 0.0 { -d } else { d }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use drusa_core::fixtures;
    use proptest::prelude::*;

    fn plane(z0: f64, gx: f64, gy: f64) -> Mesh {
        let v = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
            .map(|[x, y]| [x, y, z0 + gx * x + gy * y])
            .to_vec();
        Mesh::new(v, vec![[0, 1, 2], [0, 2, 3]]).unwrap()
    }

    #[test]
    fn surface_z_takes_the_highest_cover_and_counts_edges() {
        let c = fixtures::cube();
        let z = surface_z(&c, &[[0.5, 0.5], [1.0, 1.0], [0.0, 0.5], [2.0, 2.0]]);
        assert_eq!(&z[..3], &[1.0, 1.0, 1.0]);
        assert!(z[3].is_nan());
        let d = vertical_distance(&c, &[[0.5, 0.5, 3.0]]);
        assert_eq!(d, vec![2.0]);
    }

    #[test]
    fn normal_thickness_of_parallel_tilted_planes() {
        let (top, bottom) = (plane(5.0, 0.0, 1.0), plane(1.0, 0.0, 1.0));
        let t = thickness(&top, &bottom, &[[5.0, 5.0]], ThicknessMode::Normal).unwrap();
        assert!((t[0] - 4.0 / 2f64.sqrt()).abs() < 1e-12, "{t:?}");
        let crossed = thickness(&bottom, &top, &[[5.0, 5.0]], ThicknessMode::Normal).unwrap();
        assert!(crossed[0] < 0.0);
        let utm = plane(1.0, 0.0, 0.0).with_length_unit(Some("m".into()));
        let ft = plane(1.0, 0.0, 0.0).with_length_unit(Some("ft".into()));
        assert!(thickness(&utm, &ft, &[[1.0, 1.0]], ThicknessMode::Vertical).is_err());
    }

    proptest! {
        #[test]
        fn vertical_thickness_of_two_planes_is_analytic(
            a in -50.0..50.0f64, b in -50.0..50.0f64,
            g in proptest::array::uniform4(-1.0..1.0f64),
            p in proptest::array::uniform2(0.0..10.0f64),
        ) {
            let top = plane(a, g[0], g[1]);
            let bottom = plane(b, g[2], g[3]);
            let t = thickness(&top, &bottom, &[p], ThicknessMode::Vertical).unwrap()[0];
            let want = (a + g[0] * p[0] + g[1] * p[1]) - (b + g[2] * p[0] + g[3] * p[1]);
            prop_assert!((t - want).abs() <= 1e-9 * 1f64.max(want.abs()).max(t.abs()), "{} {}", t, want);
        }
    }
}
