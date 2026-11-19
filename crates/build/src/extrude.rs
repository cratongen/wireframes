use std::collections::HashMap;

use wireframes_core::{Column, Error, Mesh, Polylines, Result};

use crate::geom::{add, check_finite, dist, dot, ekey, newell};
use crate::polygon::triangulate_region;

/// Directed edges used by exactly one triangle of `tris`.
pub(crate) fn boundary_half_edges(tris: &[[u32; 3]]) -> Vec<[u32; 2]> {
    let mut count: HashMap<(u32, u32), u32> = HashMap::new();
    for t in tris {
        for e in 0..3 {
            *count.entry(ekey(t[e], t[(e + 1) % 3])).or_default() += 1;
        }
    }
    let mut out = vec![];
    for t in tris {
        for e in 0..3 {
            let (a, b) = (t[e], t[(e + 1) % 3]);
            if count[&ekey(a, b)] == 1 {
                out.push([a, b]);
            }
        }
    }
    out
}

pub fn extrude(polygons: &Polylines, vector: [f64; 3]) -> Result<Mesh> {
    check_finite("polygons", polygons.vertices())?;
    check_finite("vector", &[vector])?;
    let mut verts: Vec<[f64; 3]> = vec![];
    let mut tris: Vec<[u32; 3]> = vec![];
    let mut feature = vec![];
    for f in 0..polygons.len() {
        let rings: Vec<Vec<[f64; 3]>> = polygons
            .feature_parts(f)
            .filter(|&p| polygons.closed()[p])
            .map(|p| polygons.part(p).to_vec())
            .collect();
        if rings.is_empty() {
            continue;
        }
        let (cap_v, mut cap_t) = triangulate_region(&rings)?;
        let n: [f64; 3] = cap_t
            .iter()
            .map(|t| newell(&t.map(|i| cap_v[i as usize])))
            .fold([0.0; 3], add);
        let along = dot(n, vector);
        if along == 0.0 {
            return Err(Error::invalid(format!(
                "vector must not be parallel to the plane of feature {f}"
            )));
        }
        if along > 0.0 {
            for t in &mut cap_t {
                t.swap(1, 2);
            }
        }
        let base = verts.len() as u32;
        let m = cap_v.len() as u32;
        let start = tris.len();
        verts.extend(cap_v.iter().copied());
        verts.extend(cap_v.iter().map(|&p| add(p, vector)));
        tris.extend(cap_t.iter().map(|t| t.map(|i| i + base)));
        tris.extend(
            cap_t
                .iter()
                .map(|t| [t[0] + base + m, t[2] + base + m, t[1] + base + m]),
        );
        for [a, b] in boundary_half_edges(&cap_t) {
            let (a, b) = (a + base, b + base);
            let (at, bt) = (a + m, b + m);
            let v = |i: u32| verts[i as usize];
            if dist(v(b), v(at)) <= dist(v(a), v(bt)) {
                tris.push([b, a, at]);
                tris.push([b, at, bt]);
            } else {
                tris.push([b, a, bt]);
                tris.push([a, at, bt]);
            }
        }
        feature.extend(std::iter::repeat_n(f as i64, tris.len() - start));
    }
    let mesh = Mesh::new(verts, tris)?.with_face_column("feature", Column::i64(feature))?;
    Ok(mesh
        .with_crs(polygons.crs.clone())
        .with_length_unit(polygons.length_unit.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use wireframes_core::fixtures::Rng;

    fn ring(seed: u64, n: usize) -> Vec<[f64; 3]> {
        let mut rng = Rng::new(seed);
        (0..n)
            .map(|k| {
                let a = std::f64::consts::TAU * k as f64 / n as f64;
                let r = rng.uniform(5.0, 10.0);
                [1000.0 + r * a.cos(), 2000.0 + r * a.sin(), 300.0]
            })
            .collect()
    }

    proptest! {
        #[test]
        fn prism_volume_is_area_times_height(seed in 0u64..10_000, n in 3usize..20, h in -50.0f64..50.0, cw: bool) {
            prop_assume!(h.abs() > 1e-3);
            let mut r = ring(seed, n);
            if cw { r.reverse(); }
            let hole: Vec<[f64; 3]> = vec![[999.0, 1999.0, 300.0], [1001.0, 1999.0, 300.0], [1000.0, 2001.0, 300.0]];
            let p = Polylines::from_parts(vec![r, hole], vec![true, true], Some(vec![0, 0])).unwrap();
            let m = extrude(&p, [0.3, -0.2, h]).unwrap();
            let report = m.validate(0.0).unwrap();
            prop_assert!(report.summary.is_solid(), "{:?}", report.summary);
            let want = p.area()[0] * h.abs();
            let got = wireframes_ops::volume(&m).unwrap();
            prop_assert!((got - want).abs() <= 1e-9 * want.max(1.0), "{} vs {}", got, want);
        }
    }

    #[test]
    fn parallel_vector_errors() {
        let p = Polylines::from_parts(vec![ring(1, 5)], vec![true], None).unwrap();
        assert!(extrude(&p, [1.0, 0.0, 0.0]).is_err());
    }
}
