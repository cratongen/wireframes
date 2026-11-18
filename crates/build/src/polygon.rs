use drusa_core::{Column, Error, Mesh, Polylines, Result};

use crate::geom::{check_finite, newell, triangulate_rings};

/// Triangles of one planar region bounded by `rings` (even-odd), wound along
/// the Newell normal of the ring with the largest vector area.
pub(crate) fn triangulate_region(
    rings: &[Vec<[f64; 3]>],
) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
    let normals: Vec<[f64; 3]> = rings.iter().map(|r| newell(r)).collect();
    let size = |n: [f64; 3]| n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    let outer = (0..rings.len())
        .max_by(|&a, &b| {
            size(normals[a])
                .total_cmp(&size(normals[b]))
                .then(b.cmp(&a))
        })
        .ok_or_else(|| Error::invalid("polygon has no closed ring"))?;
    let n = normals[outer];
    let drop = (0..3)
        .max_by(|&a, &b| n[a].abs().total_cmp(&n[b].abs()).then(b.cmp(&a)))
        .unwrap();
    if n[drop] == 0.0 {
        return Err(Error::invalid("polygon ring has zero area"));
    }
    let (verts, mut tris) = triangulate_rings(rings, drop)?;
    if n[drop] < 0.0 {
        for t in &mut tris {
            t.swap(1, 2);
        }
    }
    Ok((verts, tris))
}

pub fn triangulate_polygons(polygons: &Polylines) -> Result<Mesh> {
    check_finite("polygons", polygons.vertices())?;
    let mut verts = vec![];
    let mut tris = vec![];
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
        let (v, t) = triangulate_region(&rings)?;
        let base = verts.len() as u32;
        tris.extend(t.iter().map(|t| t.map(|i| i + base)));
        feature.extend(std::iter::repeat_n(f as i64, t.len()));
        verts.extend(v);
    }
    let mesh = Mesh::new(verts, tris)?.with_face_column("feature", Column::i64(feature))?;
    Ok(mesh
        .with_crs(polygons.crs.clone())
        .with_length_unit(polygons.length_unit.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{dot, face_normal};
    use drusa_core::fixtures::Rng;
    use proptest::prelude::*;

    fn star(rng: &mut Rng, c: [f64; 2], r: f64, n: usize, ccw: bool) -> Vec<[f64; 3]> {
        let mut ring: Vec<[f64; 3]> = (0..n)
            .map(|k| {
                let a = std::f64::consts::TAU * k as f64 / n as f64;
                let rr = r * rng.uniform(0.5, 1.0);
                [c[0] + rr * a.cos(), c[1] + rr * a.sin(), 0.0]
            })
            .collect();
        if !ccw {
            ring.reverse();
        }
        ring
    }

    fn plan_area(m: &Mesh) -> f64 {
        (0..m.n_faces() as u32)
            .map(|f| face_normal(m.corners(f))[2] / 2.0)
            .sum()
    }

    proptest! {
        #[test]
        fn area_matches_polylines_area(seed in 0u64..10_000, n in 3usize..30, holes in 0usize..3, ccw: bool) {
            let mut rng = Rng::new(seed);
            let mut parts = vec![star(&mut rng, [0.0, 0.0], 100.0, n, ccw)];
            for h in 0..holes {
                let c = [-20.0 + 40.0 * h as f64, 0.0];
                parts.push(star(&mut rng, c, 10.0, 3 + h, !ccw));
            }
            let k = parts.len();
            let p = Polylines::from_parts(parts, vec![true; k], Some(vec![0; k])).unwrap();
            let m = triangulate_polygons(&p).unwrap();
            let want = p.area()[0];
            let got = plan_area(&m);
            let sign = if ccw { 1.0 } else { -1.0 };
            prop_assert!((got * sign - want).abs() <= 1e-9 * want.max(1.0), "{} vs {}", got, want);
            prop_assert_eq!(m.n_vertices(), p.n_vertices());
        }
    }

    #[test]
    fn vertical_polygon_follows_newell_normal() {
        let ring = vec![
            [0.0, 5.0, 0.0],
            [2.0, 5.0, 0.0],
            [2.0, 5.0, 3.0],
            [1.0, 5.0, 4.0],
            [0.0, 5.0, 3.0],
        ];
        let n = newell(&ring);
        let p = Polylines::from_parts(vec![ring], vec![true], None).unwrap();
        let m = triangulate_polygons(&p).unwrap();
        assert_eq!(m.n_faces(), 3);
        for f in 0..3 {
            assert!(dot(face_normal(m.corners(f)), n) > 0.0);
        }
        let area: f64 = (0..3)
            .map(|f| crate::geom::norm(face_normal(m.corners(f))) / 2.0)
            .sum();
        assert!((area - 7.0).abs() < 1e-12);
    }
}
