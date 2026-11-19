use wireframes_core::predicates::orient2d;
use wireframes_core::{Error, HalfEdges, Mesh, Result};

use crate::extrude::boundary_half_edges;
use crate::geom::triangulate_rings;

pub fn close_to_elevation(surface: &Mesh, z: f64) -> Result<Mesh> {
    if !z.is_finite() {
        return Err(Error::invalid(format!("z must be finite, got {z}")));
    }
    if surface.is_empty() {
        return Err(Error::invalid("surface must have faces"));
    }
    let v = surface.vertices();
    let mut up = None;
    for (f, t) in surface.triangles().iter().enumerate() {
        let c = t.map(|i| [v[i as usize][0], v[i as usize][1]]);
        let o = orient2d(c[0], c[1], c[2]);
        let s = o > 0.0;
        if o == 0.0 || up.is_some_and(|u| u != s) {
            return Err(Error::invalid(format!(
                "surface must be a height field: face {f} is vertical or folds back in plan"
            )));
        }
        up = Some(s);
    }
    let up = up.unwrap_or(true);
    let loops = HalfEdges::new(surface)?.boundary_loops(surface);
    if loops.len() != 1 {
        return Err(Error::invalid(format!(
            "surface must have exactly one boundary loop, found {}",
            loops.len()
        )));
    }
    let ring = &loops[0];
    let above = v[ring[0] as usize][2] > z;
    let used: Vec<u32> = surface.triangles().iter().flatten().copied().collect();
    if let Some(&bad) = used
        .iter()
        .find(|&&i| (v[i as usize][2] > z) != above || v[i as usize][2] == z)
    {
        return Err(Error::invalid(format!(
            "surface must lie entirely above or entirely below z = {z}: vertex {bad} is at {}",
            v[bad as usize][2]
        )));
    }
    let mut tris: Vec<[u32; 3]> = surface.triangles().to_vec();
    if up != above {
        for t in &mut tris {
            t.swap(1, 2);
        }
    }
    let mut verts = v.to_vec();
    let mut base_of = vec![u32::MAX; v.len()];
    let mut base_ring = vec![];
    for &i in ring {
        base_of[i as usize] = verts.len() as u32;
        let p = v[i as usize];
        verts.push([p[0], p[1], z]);
        base_ring.push([p[0], p[1], z]);
    }
    let n_top = tris.len();
    let rim = boundary_half_edges(&tris);
    for &[a, b] in &rim {
        let (a0, b0) = (base_of[a as usize], base_of[b as usize]);
        tris.push([b, a, a0]);
        tris.push([b, a0, b0]);
    }
    let (base_v, mut base_t) = triangulate_rings(&[base_ring.clone()], 2)?;
    if base_v.len() != base_ring.len() {
        return Err(Error::invalid(
            "surface boundary must not cross itself in plan",
        ));
    }
    let [a, b] = rim[0];
    let want = [base_of[b as usize], base_of[a as usize]];
    let local = |i: u32| base_of[ring[i as usize] as usize];
    let has = base_t
        .iter()
        .any(|t| (0..3).any(|e| [local(t[e]), local(t[(e + 1) % 3])] == want));
    if !has {
        for t in &mut base_t {
            t.swap(1, 2);
        }
    }
    tris.extend(base_t.iter().map(|t| t.map(local)));
    let face = surface.face_attributes().take_or_null(
        &(0..tris.len())
            .map(|f| (f < n_top).then_some(f as u32))
            .collect::<Vec<_>>(),
    );
    let vertex = surface.vertex_attributes().take_or_null(
        &(0..verts.len())
            .map(|i| (i < v.len()).then_some(i as u32))
            .collect::<Vec<_>>(),
    );
    let mesh = Mesh::new(verts, tris)?.with_attributes(vertex, face)?;
    Ok(mesh
        .with_crs(surface.crs.clone())
        .with_length_unit(surface.length_unit.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn plane_grid(n: u32, a: f64, b: f64, c: f64) -> Mesh {
        let mut v = vec![];
        for j in 0..n {
            for i in 0..n {
                let (x, y) = (i as f64, j as f64);
                v.push([x, y, a * x + b * y + c]);
            }
        }
        let mut f = vec![];
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let id = |i: u32, j: u32| j * n + i;
                f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        Mesh::new(v, f).unwrap()
    }

    proptest! {
        #[test]
        fn volume_is_the_integral_of_height(n in 2u32..8, a in -1.0f64..1.0, b in -1.0f64..1.0, c in 20.0f64..30.0, below: bool) {
            let m = plane_grid(n, a, b, c);
            let z = if below { 100.0 } else { 0.0 };
            let s = close_to_elevation(&m, z).unwrap();
            let report = s.validate(0.0).unwrap();
            prop_assert!(report.summary.is_solid(), "{:?}", report.summary);
            let l = (n - 1) as f64;
            let mean = a * l / 2.0 + b * l / 2.0 + c;
            let want = l * l * (mean - z).abs();
            let got = wireframes_ops::volume(&s).unwrap();
            prop_assert!((got - want).abs() <= 1e-9 * want.max(1.0), "{} vs {}", got, want);
        }
    }

    #[test]
    fn pit_surface_closes_to_a_solid() {
        let m = wireframes_core::fixtures::pit_surface(12, 3);
        let s = close_to_elevation(&m, 0.0).unwrap();
        assert!(s.validate(0.0).unwrap().summary.is_solid());
        assert!(
            close_to_elevation(&m, 90.0).is_err(),
            "the pit floor dips below 90"
        );
        assert!(
            close_to_elevation(&m, 110.0)
                .unwrap()
                .validate(0.0)
                .unwrap()
                .summary
                .is_solid()
        );
    }
}
