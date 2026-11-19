use wireframes_core::{EdgeMap, Error, Mesh, Result};

/// Vertex neighbours from unique edges, and whether each vertex is on a boundary
/// or non-manifold edge.
pub(crate) fn neighbours(mesh: &Mesh) -> (Vec<Vec<u32>>, Vec<bool>) {
    let map = EdgeMap::new(mesh, None);
    let mut nb = vec![vec![]; mesh.n_vertices()];
    let mut fixed = vec![false; mesh.n_vertices()];
    for (e, &[a, b]) in map.edges.iter().enumerate() {
        nb[a as usize].push(b);
        nb[b as usize].push(a);
        if map.faces_of_edge(e as u32).len() != 2 {
            fixed[a as usize] = true;
            fixed[b as usize] = true;
        }
    }
    (nb, fixed)
}

fn laplacian_step(v: &[[f64; 3]], nb: &[Vec<u32>], fixed: &[bool], w: f64) -> Vec<[f64; 3]> {
    v.iter()
        .enumerate()
        .map(|(i, &p)| {
            if fixed[i] || nb[i].is_empty() {
                return p;
            }
            let k = nb[i].len() as f64;
            let mut m = [0.0; 3];
            for &j in &nb[i] {
                for c in 0..3 {
                    m[c] += v[j as usize][c] - p[c];
                }
            }
            [0, 1, 2].map(|c| p[c] + w * m[c] / k)
        })
        .collect()
}

/// With `preserve_volume`, the mesh must be a closed outward solid; the smoothed
/// vertices are scaled about their volume centroid so the volume matches the
/// input's to rounding.
pub fn smooth_taubin(
    mesh: &Mesh,
    iterations: usize,
    lambda: f64,
    mu: f64,
    preserve_volume: bool,
) -> Result<Mesh> {
    if !(lambda.is_finite() && lambda > 0.0) {
        return Err(Error::invalid(format!("lambda must be > 0, got {lambda}")));
    }
    if !(mu.is_finite() && mu < -lambda) {
        return Err(Error::invalid(format!("mu must be < -lambda, got {mu}")));
    }
    let v0 = if preserve_volume {
        let v0 = wireframes_ops::volume(mesh)?;
        if !(v0.is_finite() && v0 > 0.0) {
            return Err(Error::invalid(format!(
                "preserve_volume needs a positive volume, got {v0}"
            )));
        }
        Some(v0)
    } else {
        None
    };
    let (nb, fixed) = neighbours(mesh);
    let mut v = mesh.vertices().to_vec();
    for _ in 0..iterations {
        v = laplacian_step(&v, &nb, &fixed, lambda);
        v = laplacian_step(&v, &nb, &fixed, mu);
    }
    let out = mesh.clone().with_vertices(v)?;
    let Some(v0) = v0 else { return Ok(out) };
    let v1 = wireframes_ops::volume(&out)?;
    if !(v1.is_finite() && v1 > 0.0) {
        return Err(Error::invalid(format!(
            "smoothing left a volume of {v1}, cannot rescale to {v0}"
        )));
    }
    let c = wireframes_ops::centroid(&out)?;
    let s = (v0 / v1).cbrt();
    let scaled = out
        .vertices()
        .iter()
        .map(|p| [0, 1, 2].map(|k| c[k] + s * (p[k] - c[k])))
        .collect();
    out.with_vertices(scaled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use wireframes_core::fixtures::{flipped, pit_surface, sphere};

    proptest! {
        #[test]
        fn boundary_fixed_and_topology_kept(n in 3u32..12, seed in 0u64..1000, it in 0usize..10) {
            let m = pit_surface(n, seed);
            let s = smooth_taubin(&m, it, 0.5, -0.53, false).unwrap();
            let (_, fixed) = neighbours(&m);
            prop_assert_eq!(s.triangles(), m.triangles());
            for i in 0..m.n_vertices() {
                if fixed[i] {
                    prop_assert_eq!(s.vertices()[i], m.vertices()[i]);
                }
            }
        }
    }

    #[test]
    fn closed_sphere_stays_a_solid_and_keeps_its_size() {
        let m = sphere(2);
        let s = smooth_taubin(&m, 10, 0.5, -0.53, false).unwrap();
        assert!(s.validate(0.0).unwrap().summary.is_solid());
        let v0 = wireframes_ops::volume(&m).unwrap();
        let v1 = wireframes_ops::volume(&s).unwrap();
        assert!((v1 / v0 - 1.0).abs() < 0.05, "{v0} -> {v1}");
    }

    #[test]
    fn preserve_volume_keeps_the_volume_to_rounding() {
        let m = sphere(2);
        let v0 = wireframes_ops::volume(&m).unwrap();
        let drift =
            wireframes_ops::volume(&smooth_taubin(&m, 10, 0.5, -0.53, false).unwrap()).unwrap();
        assert!((drift / v0 - 1.0).abs() > 1e-6, "{v0} -> {drift}");
        let s = smooth_taubin(&m, 10, 0.5, -0.53, true).unwrap();
        assert!(s.validate(0.0).unwrap().summary.is_solid());
        assert_eq!(s.triangles(), m.triangles());
        let v1 = wireframes_ops::volume(&s).unwrap();
        assert!((v1 / v0 - 1.0).abs() < 1e-12, "{v0} -> {v1}");
    }

    #[test]
    fn preserve_volume_refuses_open_and_inward_meshes() {
        assert!(smooth_taubin(&pit_surface(5, 1), 1, 0.5, -0.53, true).is_err());
        assert!(smooth_taubin(&flipped(&sphere(1)), 1, 0.5, -0.53, true).is_err());
    }

    #[test]
    fn bad_weights_error() {
        let m = sphere(0);
        assert!(smooth_taubin(&m, 1, 0.0, -0.5, false).is_err());
        assert!(smooth_taubin(&m, 1, 0.5, -0.4, false).is_err());
    }
}
