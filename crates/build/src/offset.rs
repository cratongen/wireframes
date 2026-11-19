use wireframes_core::{Error, Mesh, Result};

use crate::geom::{add, cross, dot, norm, normalize, scale, sub};

/// Angle-weighted vertex normals (unit length, zero where undefined).
pub(crate) fn vertex_normals(mesh: &Mesh) -> Vec<[f64; 3]> {
    let v = mesh.vertices();
    let mut n = vec![[0.0; 3]; v.len()];
    for t in mesh.triangles() {
        let c = t.map(|i| v[i as usize]);
        let fnorm = normalize(cross(sub(c[1], c[0]), sub(c[2], c[0])));
        for k in 0..3 {
            let e1 = sub(c[(k + 1) % 3], c[k]);
            let e2 = sub(c[(k + 2) % 3], c[k]);
            let (l1, l2) = (norm(e1), norm(e2));
            if l1 == 0.0 || l2 == 0.0 {
                continue;
            }
            let angle = (dot(e1, e2) / (l1 * l2)).clamp(-1.0, 1.0).acos();
            n[t[k] as usize] = add(n[t[k] as usize], scale(fnorm, angle));
        }
    }
    n.into_iter().map(normalize).collect()
}

pub fn offset(mesh: &Mesh, distance: f64) -> Result<Mesh> {
    if !distance.is_finite() {
        return Err(Error::invalid(format!(
            "distance must be finite, got {distance}"
        )));
    }
    let normals = vertex_normals(mesh);
    let v: Vec<[f64; 3]> = mesh
        .vertices()
        .iter()
        .zip(&normals)
        .map(|(&p, &n)| add(p, scale(n, distance)))
        .collect();
    let out = mesh.clone().with_vertices(v)?;
    let report = out.validate(0.0)?;
    let s = &report.summary;
    if s.self_intersections > 0 || s.coplanar_overlaps > 0 || s.inward_shells > 0 {
        return Err(Error::InvalidMesh {
            context: "offset result self-intersects".into(),
            report: Box::new(report),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wireframes_core::fixtures::{cube, sphere};

    #[test]
    fn sphere_grows_and_shrinks() {
        let m = sphere(3);
        let v0 = wireframes_ops::volume(&m).unwrap();
        let big = offset(&m, 0.1).unwrap();
        let small = offset(&m, -0.1).unwrap();
        let r = |d: f64| (1.0f64 + d).powi(3);
        assert!((wireframes_ops::volume(&big).unwrap() / v0 / r(0.1) - 1.0).abs() < 0.01);
        assert!((wireframes_ops::volume(&small).unwrap() / v0 / r(-0.1) - 1.0).abs() < 0.01);
        assert!(big.validate(0.0).unwrap().summary.is_solid());
    }

    #[test]
    fn collapsing_offset_is_refused() {
        assert!(matches!(
            offset(&cube(), -2.0),
            Err(Error::InvalidMesh { .. })
        ));
    }
}
