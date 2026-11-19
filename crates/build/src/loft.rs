use wireframes_core::{Column, Error, Mesh, Polylines, Result};

use crate::geom::{add, centroid, check_finite, cross, dot, newell, norm, normalize, scale, sub};
use crate::polygon::triangulate_region;

#[derive(Clone, Debug)]
pub struct LoftOptions {
    pub cap_ends: bool,
}

impl Default for LoftOptions {
    fn default() -> LoftOptions {
        LoftOptions { cap_ends: true }
    }
}

fn tri_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    norm(cross(sub(b, a), sub(c, a))) / 2.0
}

/// Area centroid of a planar ring, or its vertex average when the ring has no area.
fn ring_centre(ring: &[[f64; 3]]) -> [f64; 3] {
    let c = centroid(ring);
    let n = newell(ring);
    let (mut sum, mut w) = ([0.0; 3], 0.0);
    for (i, &p) in ring.iter().enumerate() {
        let q = ring[(i + 1) % ring.len()];
        let a = dot(cross(sub(p, c), sub(q, c)), n);
        sum = add(sum, scale(add(add(c, p), q), a / 3.0));
        w += a;
    }
    if w > 0.0 { scale(sum, 1.0 / w) } else { c }
}

/// Places a ring in the frame of a band: centred on its own centre, with the
/// component along the band axis removed, so the band length does not enter the
/// linking cost.
fn overlay(ring: &[[f64; 3]], centre: [f64; 3], axis: [f64; 3]) -> Vec<[f64; 3]> {
    ring.iter()
        .map(|&p| {
            let d = sub(p, centre);
            sub(d, scale(axis, dot(d, axis)))
        })
        .collect()
}

/// Band between rings `a` and `b` (overlay coordinates), both starting at their
/// aligned vertex and running the same way, that minimizes the total overlay
/// triangle area. Returns local corner indices: `(false, i)` is `a[i % m]`,
/// `(true, j)` is `b[j % n]`. A ring of one vertex adds no triangles of its own.
fn link(a: &[[f64; 3]], b: &[[f64; 3]]) -> Vec<[(bool, usize); 3]> {
    let (m, n) = (a.len(), b.len());
    let pa = |i: usize| a[i % m];
    let pb = |j: usize| b[j % n];
    let w = n + 1;
    let mut cost = vec![f64::INFINITY; (m + 1) * w];
    let mut from_a = vec![false; (m + 1) * w];
    cost[0] = 0.0;
    for i in 0..=m {
        for j in 0..=n {
            let here = cost[i * w + j];
            if i < m {
                let t = if m == 1 {
                    0.0
                } else {
                    tri_area(pa(i), pa(i + 1), pb(j))
                };
                let k = (i + 1) * w + j;
                if here + t < cost[k] {
                    cost[k] = here + t;
                    from_a[k] = true;
                }
            }
            if j < n {
                let t = if n == 1 {
                    0.0
                } else {
                    tri_area(pa(i), pb(j + 1), pb(j))
                };
                let k = i * w + j + 1;
                if here + t < cost[k] {
                    cost[k] = here + t;
                    from_a[k] = false;
                }
            }
        }
    }
    let (mut i, mut j) = (m, n);
    let mut out = vec![];
    while i > 0 || j > 0 {
        if from_a[i * w + j] {
            i -= 1;
            if m > 1 {
                out.push([(false, i % m), (false, (i + 1) % m), (true, j % n)]);
            }
        } else {
            j -= 1;
            if n > 1 {
                out.push([(false, i % m), (true, (j + 1) % n), (true, j % n)]);
            }
        }
    }
    out.reverse();
    out
}

pub fn loft(sections: &Polylines, opts: &LoftOptions) -> Result<Mesh> {
    check_finite("sections", sections.vertices())?;
    if sections.len() < 2 {
        return Err(Error::invalid(format!(
            "sections must hold at least 2 features, got {}",
            sections.len()
        )));
    }
    let k = sections.len();
    let mut rings: Vec<Vec<[f64; 3]>> = vec![];
    let mut closed = vec![];
    for f in 0..k {
        let parts = sections.feature_parts(f);
        if parts.len() != 1 {
            return Err(Error::Unsupported(format!(
                "loft links one ring per section; section {f} has {} parts",
                parts.len()
            )));
        }
        let pts = sections.part(parts.start);
        let is_closed = sections.closed()[parts.start];
        if !is_closed {
            return Err(Error::invalid(format!("section {f} must be a closed ring")));
        }
        if is_closed || pts.len() < 3 {
            rings.push(pts.to_vec());
        } else {
            let back = pts[1..pts.len() - 1].iter().rev();
            rings.push(pts.iter().chain(back).copied().collect());
        }
        closed.push(is_closed);
    }
    let centres: Vec<[f64; 3]> = (0..k)
        .map(|s| {
            if closed[s] {
                ring_centre(&rings[s])
            } else {
                centroid(sections.part(sections.feature_parts(s).start))
            }
        })
        .collect();
    for s in 0..k - 1 {
        if centres[s] == centres[s + 1] {
            return Err(Error::invalid(format!(
                "sections {s} and {} have the same centre",
                s + 1
            )));
        }
    }
    for s in 0..k {
        if !closed[s] {
            continue;
        }
        let dir = if s + 1 < k {
            sub(centres[s + 1], centres[s])
        } else {
            sub(centres[s], centres[s - 1])
        };
        let along = dot(newell(&rings[s]), dir);
        if along == 0.0 {
            return Err(Error::invalid(format!(
                "section {s} lies parallel to the direction of the next section"
            )));
        }
        if along < 0.0 {
            rings[s].reverse();
        }
    }
    let mut verts = vec![];
    let mut ids: Vec<Vec<u32>> = vec![];
    for (s, r) in rings.iter().enumerate() {
        let base = verts.len() as u32;
        let n = if closed[s] {
            r.len()
        } else {
            sections.part(sections.feature_parts(s).start).len()
        };
        ids.push(
            (0..r.len())
                .map(|i| base + if i < n { i } else { 2 * (n - 1) - i } as u32)
                .collect(),
        );
        verts.extend(r[..n].iter().copied());
    }
    let mut tris = vec![];
    let mut section = vec![];
    for s in 0..k - 1 {
        let axis = normalize(sub(centres[s + 1], centres[s]));
        let a = overlay(&rings[s], centres[s], axis);
        let b = overlay(&rings[s + 1], centres[s + 1], axis);
        let (mut best, mut bi, mut bj) = (f64::INFINITY, 0, 0);
        for (i, &p) in a.iter().enumerate() {
            for (j, &q) in b.iter().enumerate() {
                let d = sub(p, q);
                let d2 = dot(d, d);
                if d2 < best {
                    (best, bi, bj) = (d2, i, j);
                }
            }
        }
        let (m, n) = (a.len(), b.len());
        let ra: Vec<[f64; 3]> = (0..m).map(|i| a[(i + bi) % m]).collect();
        let rb: Vec<[f64; 3]> = (0..n).map(|j| b[(j + bj) % n]).collect();
        let band: Vec<[u32; 3]> = link(&ra, &rb)
            .into_iter()
            .map(|t| {
                t.map(|(on_b, i)| {
                    if on_b {
                        ids[s + 1][(i + bj) % n]
                    } else {
                        ids[s][(i + bi) % m]
                    }
                })
            })
            .collect();
        section.extend(std::iter::repeat_n(s as i64, band.len()));
        tris.extend(band);
    }
    if opts.cap_ends {
        for (s, label) in [(0, -1i64), (k - 1, -2)] {
            if !closed[s] {
                continue;
            }
            let (cap_v, cap_t) = triangulate_region(std::slice::from_ref(&rings[s]))?;
            if cap_v.len() != rings[s].len() {
                return Err(Error::invalid(format!("section {s} ring crosses itself")));
            }
            for t in cap_t {
                let t = t.map(|i| ids[s][i as usize]);
                tris.push(if s == 0 { [t[0], t[2], t[1]] } else { t });
                section.push(label);
            }
        }
    }
    let mesh = Mesh::new(verts, tris)?.with_face_column("section", Column::i64(section))?;
    if opts.cap_ends {
        let report = mesh.validate(0.0)?;
        if !report.summary.is_solid() {
            return Err(Error::invalid(format!(
                "linked sections do not form a solid: {:?}",
                report.summary.to_pairs()
            )));
        }
    }
    Ok(mesh
        .with_crs(sections.crs.clone())
        .with_length_unit(sections.length_unit.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn square(a: f64, z: f64, ccw: bool) -> Vec<[f64; 3]> {
        let mut r = vec![[0.0, 0.0, z], [a, 0.0, z], [a, a, z], [0.0, a, z]];
        if !ccw {
            r.reverse();
        }
        r
    }

    proptest! {
        #[test]
        fn two_squares_make_a_box(a in 0.5f64..50.0, h in 0.5f64..50.0, ccw0: bool, ccw1: bool) {
            let p = Polylines::from_parts(vec![square(a, 10.0, ccw0), square(a, 10.0 + h, ccw1)], vec![true, true], None).unwrap();
            let m = loft(&p, &LoftOptions::default()).unwrap();
            let want = a * a * h;
            let got = wireframes_ops::volume(&m).unwrap();
            prop_assert!((got - want).abs() <= 1e-9 * want.max(1.0), "{} vs {}", got, want);
        }
    }

    #[test]
    fn rings_of_different_sizes_link_into_a_solid() {
        let circle = |n: usize, r: f64, z: f64| -> Vec<[f64; 3]> {
            (0..n)
                .map(|k| {
                    let t = std::f64::consts::TAU * k as f64 / n as f64;
                    [r * t.cos(), r * t.sin(), z]
                })
                .collect()
        };
        let p = Polylines::from_parts(
            vec![
                circle(7, 3.0, 0.0),
                circle(12, 4.0, 5.0),
                circle(5, 2.0, 9.0),
            ],
            vec![true; 3],
            None,
        )
        .unwrap();
        let m = loft(&p, &LoftOptions::default()).unwrap();
        assert!(wireframes_ops::volume(&m).unwrap() > 0.0);
        let open = loft(&p, &LoftOptions { cap_ends: false }).unwrap();
        assert!(!open.is_closed());
    }

    fn ring(n: usize, r: f64, at: impl Fn(f64, f64) -> [f64; 3]) -> Vec<[f64; 3]> {
        (0..n)
            .map(|k| {
                let t = std::f64::consts::TAU * k as f64 / n as f64;
                at(r * t.cos(), r * t.sin())
            })
            .collect()
    }

    fn polygon_area(n: usize, r: f64) -> f64 {
        0.5 * n as f64 * r * r * (std::f64::consts::TAU / n as f64).sin()
    }

    #[test]
    fn long_band_to_an_oblique_section_is_a_truncated_prism() {
        for n in [4, 8, 11, 16] {
            for (len, slope) in [(5.0, 0.5), (40.0, 0.05), (40.0, 0.125)] {
                let a = ring(n, 2.5, |y, z| [0.0, y, z]);
                let b = ring(n, 2.5, |y, z| [len + slope * y, y, z]);
                let p = Polylines::from_parts(vec![a, b], vec![true; 2], None).unwrap();
                let m = loft(&p, &LoftOptions::default()).unwrap();
                let want = polygon_area(n, 2.5) * len;
                let got = wireframes_ops::volume(&m).unwrap();
                assert!(
                    (got - want).abs() <= 1e-9 * want,
                    "n {n} len {len}: {got} vs {want}"
                );
            }
        }
    }

    #[test]
    fn mitred_rings_at_a_bend_make_the_swept_prisms() {
        let (n, r, len) = (16, 2.5, 40.0);
        for bend in [0.1f64, 0.3, 0.6] {
            let (c, s, h) = (bend.cos(), bend.sin(), (bend / 2.0).tan());
            let a = ring(n, r, |y, z| [0.0, y, z]);
            let b = ring(n, r, |y, z| [len - h * y, y, z]);
            let d = ring(n, r, |y, z| [len + c * len - s * y, s * len + c * y, z]);
            let p = Polylines::from_parts(vec![a, b, d], vec![true; 3], None).unwrap();
            let m = loft(&p, &LoftOptions::default()).unwrap();
            let want = polygon_area(n, r) * 2.0 * len;
            let got = wireframes_ops::volume(&m).unwrap();
            assert!(
                (got - want).abs() <= 1e-9 * want,
                "bend {bend}: {got} vs {want}"
            );
        }
    }

    #[test]
    fn multipart_section_is_unsupported() {
        let p = Polylines::from_parts(
            vec![
                square(1.0, 0.0, true),
                square(1.0, 1.0, true),
                square(1.0, 2.0, true),
            ],
            vec![true; 3],
            Some(vec![0, 1, 1]),
        )
        .unwrap();
        assert!(matches!(
            loft(&p, &LoftOptions::default()),
            Err(Error::Unsupported(_))
        ));
    }
}
