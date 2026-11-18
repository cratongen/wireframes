use drusa_core::{Error, Mesh, Result, units};

fn with_geometry(mesh: &Mesh, vertices: Vec<[f64; 3]>, flip: bool) -> Result<Mesh> {
    let triangles = if flip {
        mesh.triangles()
            .iter()
            .map(|&[a, b, c]| [a, c, b])
            .collect()
    } else {
        mesh.triangles().to_vec()
    };
    Ok(Mesh::new(vertices, triangles)?
        .with_attributes(
            mesh.vertex_attributes().clone(),
            mesh.face_attributes().clone(),
        )?
        .with_crs(mesh.crs.clone())
        .with_length_unit(mesh.length_unit.clone()))
}

/// Affine `m`, row-major, applied to column vectors; flips faces when `det < 0`.
pub fn transform(mesh: &Mesh, m: [[f64; 4]; 4]) -> Result<Mesh> {
    if m.iter().flatten().any(|x| !x.is_finite()) {
        return Err(Error::invalid("m must be finite"));
    }
    if m[3] != [0.0, 0.0, 0.0, 1.0] {
        return Err(Error::invalid(format!(
            "m must be affine: last row is {:?}, expected [0, 0, 0, 1]",
            m[3]
        )));
    }
    let a = |i: usize, j: usize| m[i][j];
    let det = a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
        - a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0))
        + a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0));
    let vertices = mesh
        .vertices()
        .iter()
        .map(|p| {
            std::array::from_fn(|i| m[i][0] * p[0] + m[i][1] * p[1] + m[i][2] * p[2] + m[i][3])
        })
        .collect();
    with_geometry(mesh, vertices, det < 0.0)
}

/// `InvalidInput` when `d` or a translated coordinate is not finite.
pub fn translate(mesh: &Mesh, d: [f64; 3]) -> Result<Mesh> {
    if d.iter().any(|x| !x.is_finite()) {
        return Err(Error::invalid(format!("d must be finite: {d:?}")));
    }
    let vertices = mesh
        .vertices()
        .iter()
        .map(|p| [p[0] + d[0], p[1] + d[1], p[2] + d[2]])
        .collect();
    mesh.clone().with_vertices(vertices)
}

/// `InvalidInput` when the mesh declares no unit.
pub fn to_length_unit(mesh: &Mesh, unit: &str) -> Result<Mesh> {
    let from = mesh
        .length_unit
        .as_deref()
        .ok_or_else(|| Error::invalid("mesh declares no length unit to convert from"))?;
    let k = units::length_factor(from)? / units::length_factor(unit)?;
    let vertices = mesh.vertices().iter().map(|p| p.map(|x| x * k)).collect();
    Ok(mesh
        .clone()
        .with_vertices(vertices)?
        .with_length_unit(Some(unit.to_string())))
}

pub fn flip(mesh: &Mesh) -> Mesh {
    with_geometry(mesh, mesh.vertices().to_vec(), true).expect("same vertices and faces")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::measure::volume;
    use drusa_core::fixtures;

    #[test]
    fn mirror_keeps_volume_positive() {
        let m = fixtures::box_mesh([0.0; 3], [1.0, 2.0, 3.0]).with_crs(Some("EPSG:31983".into()));
        let mirror = [
            [-1.0, 0.0, 0.0, 10.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 2.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let t = transform(&m, mirror).unwrap();
        assert_eq!(volume(&t).unwrap(), 12.0);
        assert_eq!(t.crs, m.crs);
        assert_eq!(t.vertices()[1][0], 9.0);
        let mut bad = mirror;
        bad[3][0] = 1.0;
        assert!(transform(&m, bad).is_err());
    }

    #[test]
    fn translate_flip_and_units() {
        let m = fixtures::cube();
        assert_eq!(volume(&flip(&m)).unwrap(), -1.0);
        assert_eq!(
            translate(&m, [1.0, 2.0, 3.0]).unwrap().vertices()[0],
            [1.0, 2.0, 3.0]
        );
        assert!(translate(&m, [f64::NAN, 0.0, 0.0]).is_err());
        assert!(
            translate(
                &m.clone().with_vertices(vec![[f64::MAX; 3]; 8]).unwrap(),
                [f64::MAX; 3]
            )
            .is_err()
        );
        assert!(to_length_unit(&m, "ft").is_err());
        let ft = to_length_unit(&m.with_length_unit(Some("m".into())), "ft").unwrap();
        assert_eq!(ft.length_unit.as_deref(), Some("ft"));
        assert!(
            (ft.vertices()[6].iter().fold(0f64, |a, &x| a.max(x)) - 1.0 / 0.3048).abs() < 1e-12
        );
    }
}
