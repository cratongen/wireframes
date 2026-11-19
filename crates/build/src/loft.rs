use wireframes_core::{Column, Error, Mesh, Polylines, Result};

use crate::geom::{centroid, check_finite, cross, dot, newell, norm, sub};
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

/// Minimum-area band between rings `a` and `b`, both starting at their aligned
/// vertex and running the same way.
fn link(verts: &[[f64; 3]], a: &[u32], b: &[u32]) -> Vec<[u32; 3]> {
    let (m, n) = (a.len(), b.len());
    let pa = |i: usize| verts[a[i % m] as usize];
    let pb = |j: usize| verts[b[j % n] as usize];
    let w = n + 1;
    let mut cost = vec![f64::INFINITY; (m + 1) * w];
    let mut from_a = vec![false; (m + 1) * w];
    cost[0] = 0.0;
    for i in 0..=m {
        for j in 0..=n {
            let here = cost[i * w + j];
            if i < m {
                let c = here + tri_area(pa(i), pa(i + 1), pb(j));
                let k = (i + 1) * w + j;
                if c < cost[k] {
                    cost[k] = c;
                    from_a[k] = true;
                }
            }
            if j < n {
                let c = here + tri_area(pa(i), pb(j + 1), pb(j));
                let k = i * w + j + 1;
                if c < cost[k] {
                    cost[k] = c;
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
            out.push([a[i % m], a[(i + 1) % m], b[j % n]]);
        } else {
            j -= 1;
            out.push([a[i % m], b[(j + 1) % n], b[j % n]]);
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
    let mut rings: Vec<Vec<[f64; 3]>> = vec![];
    for f in 0..sections.len() {
        let parts = sections.feature_parts(f);
        if parts.len() != 1 {
            return Err(Error::Unsupported(format!(
                "loft links one ring per section; section {f} has {} parts",
                parts.len()
            )));
        }
        if !sections.closed()[parts.start] {
            return Err(Error::invalid(format!("section {f} must be a closed ring")));
        }
        rings.push(sections.part(parts.start).to_vec());
    }
    let centres: Vec<[f64; 3]> = rings.iter().map(|r| centroid(r)).collect();
    let k = rings.len();
    for s in 0..k {
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
    for r in &rings {
        let base = verts.len() as u32;
        ids.push((0..r.len() as u32).map(|i| base + i).collect());
        verts.extend(r.iter().copied());
    }
    let mut tris = vec![];
    let mut section = vec![];
    for s in 0..k - 1 {
        let (a, b) = (&ids[s], &ids[s + 1]);
        let (mut best, mut bi, mut bj) = (f64::INFINITY, 0, 0);
        for i in 0..a.len() {
            for j in 0..b.len() {
                let d = sub(verts[a[i] as usize], verts[b[j] as usize]);
                let d2 = dot(d, d);
                if d2 < best {
                    (best, bi, bj) = (d2, i, j);
                }
            }
        }
        let ra: Vec<u32> = (0..a.len()).map(|i| a[(i + bi) % a.len()]).collect();
        let rb: Vec<u32> = (0..b.len()).map(|j| b[(j + bj) % b.len()]).collect();
        let band = link(&verts, &ra, &rb);
        section.extend(std::iter::repeat_n(s as i64, band.len()));
        tris.extend(band);
    }
    if opts.cap_ends {
        for (s, label) in [(0, -1i64), (k - 1, -2)] {
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
