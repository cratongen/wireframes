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

pub fn smooth_taubin(mesh: &Mesh, iterations: usize, lambda: f64, mu: f64) -> Result<Mesh> {
    if !(lambda.is_finite() && lambda > 0.0) {
        return Err(Error::invalid(format!("lambda must be > 0, got {lambda}")));
    }
    if !(mu.is_finite() && mu < -lambda) {
        return Err(Error::invalid(format!("mu must be < -lambda, got {mu}")));
    }
    let (nb, fixed) = neighbours(mesh);
    let mut v = mesh.vertices().to_vec();
    for _ in 0..iterations {
        v = laplacian_step(&v, &nb, &fixed, lambda);
        v = laplacian_step(&v, &nb, &fixed, mu);
    }
    mesh.clone().with_vertices(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use wireframes_core::fixtures::{pit_surface, sphere};

    proptest! {
        #[test]
        fn boundary_fixed_and_topology_kept(n in 3u32..12, seed in 0u64..1000, it in 0usize..10) {
            let m = pit_surface(n, seed);
            let s = smooth_taubin(&m, it, 0.5, -0.53).unwrap();
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
        let s = smooth_taubin(&m, 10, 0.5, -0.53).unwrap();
        assert!(s.validate(0.0).unwrap().summary.is_solid());
        let v0 = wireframes_ops::volume(&m).unwrap();
        let v1 = wireframes_ops::volume(&s).unwrap();
        assert!((v1 / v0 - 1.0).abs() < 0.05, "{v0} -> {v1}");
    }

    #[test]
    fn bad_weights_error() {
        let m = sphere(0);
        assert!(smooth_taubin(&m, 1, 0.0, -0.5).is_err());
        assert!(smooth_taubin(&m, 1, 0.5, -0.4).is_err());
    }
}
