use std::collections::HashSet;

use wireframes_boolean::slice_by_plane;
use wireframes_core::{Column, Error, Mesh, Polylines, Result, Table};

fn bits(p: [f64; 3]) -> [u64; 3] {
    p.map(|x| (x + 0.0).to_bits())
}

fn segment(a: [f64; 3], b: [f64; 3]) -> [[u64; 3]; 2] {
    let (a, b) = (bits(a), bits(b));
    if a < b { [a, b] } else { [b, a] }
}

/// The segments of each part, the closing one included for rings.
fn segments(lines: &Polylines, p: usize) -> Vec<[[f64; 3]; 2]> {
    let v = lines.part(p);
    let n = if lines.closed()[p] {
        v.len()
    } else {
        v.len() - 1
    };
    (0..n).map(|i| [v[i], v[(i + 1) % v.len()]]).collect()
}

/// Lines at each level: `slice_by_plane` with normal `+z` (a vertex on the level
/// counts above, so this traces where the surface drops below the level), then
/// each line of `slice_by_plane` with normal `-z` (a vertex on the level counts
/// below) that shares no segment with them, reversed to the `+z` winding: the
/// toe of a flat area at the level that nothing crosses, such as a pit floor.
/// A line that shares a segment is the same contour going round a flat
/// triangle on its other side, and is left out.
pub fn contours(surface: &Mesh, levels: &[f64]) -> Result<Polylines> {
    if surface.face_attributes().column("elevation").is_some() {
        return Err(Error::invalid(
            "surface must not have a face column named \"elevation\"",
        ));
    }
    let levels = if levels.is_empty() {
        &[0.0][..]
    } else {
        levels
    };
    let mut out: Vec<Polylines> = vec![];
    for &level in levels {
        let up = slice_by_plane(surface, [0.0; 3], [0.0, 0.0, 1.0], &[level])?;
        let down = slice_by_plane(surface, [0.0; 3], [0.0, 0.0, -1.0], &[-level])?;
        let seen: HashSet<[[u64; 3]; 2]> = (0..up.n_parts())
            .flat_map(|p| segments(&up, p))
            .map(|[a, b]| segment(a, b))
            .collect();
        let (mut parts, mut closed, mut rows) = (vec![], vec![], vec![]);
        for f in 0..down.len() {
            for p in down.feature_parts(f) {
                let segs: Vec<[[f64; 3]; 2]> = segments(&down, p)
                    .into_iter()
                    .rev()
                    .map(|[a, b]| [b, a])
                    .collect();
                let Some(last) = segs.last() else { continue };
                if segs.iter().all(|s| !seen.contains(&segment(s[0], s[1]))) {
                    let mut line: Vec<[f64; 3]> = segs.iter().map(|s| s[0]).collect();
                    if !down.closed()[p] {
                        line.push(last[1]);
                    }
                    parts.push(line);
                    closed.push(down.closed()[p]);
                    rows.push(f as u32);
                }
            }
        }
        out.push(up);
        if !parts.is_empty() {
            let n = rows.len();
            let attrs = down
                .attributes()
                .take(&rows)
                .without_column("offset")
                .with_column("offset", Column::f64(vec![level; n]))?;
            out.push(Polylines::from_parts(parts, closed, None)?.with_attributes(attrs)?);
        }
    }
    let lines = Polylines::concat(&out.iter().collect::<Vec<_>>())?;
    let mut attrs = Table::new(lines.len());
    for name in lines.attributes().names() {
        let mut c = lines.attributes().column(name).expect("listed").clone();
        if c.valid.as_ref().is_some_and(|v| v.iter().all(|&x| x)) {
            c.valid = None;
        }
        attrs = attrs.with_column(name, c)?;
    }
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
    use proptest::prelude::*;
    use wireframes_core::Column;
    use wireframes_core::fixtures::{box_mesh, sphere};

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

    fn grid(z: impl Fn(f64, f64) -> f64) -> Mesh {
        let mut v = vec![];
        for j in -3..=3 {
            for i in -3..=3 {
                let (x, y) = (i as f64, j as f64);
                v.push([x, y, z(x, y)]);
            }
        }
        let id = |i: u32, j: u32| j * 7 + i;
        let mut f = vec![];
        for j in 0..6 {
            for i in 0..6 {
                f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        Mesh::new(v, f).unwrap()
    }

    fn signed_area(r: &[[f64; 3]]) -> f64 {
        (0..r.len())
            .map(|i| {
                let (a, b) = (r[i], r[(i + 1) % r.len()]);
                a[0] * b[1] - b[0] * a[1]
            })
            .sum::<f64>()
            / 2.0
    }

    #[test]
    fn a_flat_floor_at_the_level_gives_its_toe_once() {
        let pit = grid(|x, y| (x.abs().max(y.abs()) - 1.0).max(0.0));
        let toe = contours(&pit, &[0.0]).unwrap();
        assert_eq!(toe.len(), 1);
        assert!(toe.closed()[0]);
        let just_above = contours(&pit, &[1e-9]).unwrap();
        let (a, b) = (signed_area(toe.part(0)), signed_area(just_above.part(0)));
        assert!((a.abs() - 4.0).abs() < 1e-12 && a * b > 0.0, "{a} {b}");
        let hill = grid(|x, y| -(x.abs().max(y.abs()) - 1.0).max(0.0));
        let crest = contours(&hill, &[0.0]).unwrap();
        assert_eq!(crest.len(), 1);
        assert!((signed_area(crest.part(0)).abs() - 4.0).abs() < 1e-12);
        let bench = grid(|x, y| {
            let r = x.abs().max(y.abs());
            if r <= 1.0 {
                -1.0
            } else if r <= 2.0 {
                0.0
            } else {
                1.0
            }
        });
        let lines = contours(&bench, &[0.0]).unwrap();
        assert_eq!(lines.closed(), &[true]);
        assert!((signed_area(lines.part(0)).abs() - 15.0).abs() < 1e-12);
        let berm = grid(|x, y| {
            let r = x.abs().max(y.abs());
            if r == 0.0 {
                -1.0
            } else if r <= 2.0 {
                0.0
            } else {
                1.0
            }
        });
        let lines = contours(&berm, &[0.0]).unwrap();
        let mut areas: Vec<f64> = (0..lines.len())
            .map(|p| signed_area(lines.part(p)).abs())
            .collect();
        areas.sort_by(f64::total_cmp);
        assert_eq!(lines.closed(), &[true, true]);
        assert_eq!(areas, [3.0, 16.0]);
        let e = lines.attributes().require("elevation").unwrap();
        assert_eq!(e, &Column::f64(vec![0.0; 2]));
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
