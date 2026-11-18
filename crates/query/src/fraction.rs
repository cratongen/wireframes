use drusa_core::{Aabb, Bvh, Error, Mesh, Result};
use drusa_ops::box_volume;
use rayon::prelude::*;

use crate::SurfaceIndex;
use crate::distance::check_finite;

#[derive(Clone, Debug, PartialEq)]
pub struct Blocks {
    pub centroids: Vec<[f64; 3]>,
    /// Positive and finite.
    pub size: [f64; 3],
    /// `rotation[i][k]` is row `i`, column `k`. Columns are the block axes in
    /// world coordinates; orthonormal with determinant +1 (checked to 1e-12).
    pub rotation: Option<[[f64; 3]; 3]>,
}

fn check(blocks: &Blocks) -> Result<()> {
    check_finite("centroids", &blocks.centroids)?;
    if blocks.size.iter().any(|s| !(s.is_finite() && *s > 0.0)) {
        return Err(Error::invalid(format!(
            "size must be positive and finite: size is {:?}",
            blocks.size
        )));
    }
    if let Some(r) = blocks.rotation {
        let col = |k: usize| [r[0][k], r[1][k], r[2][k]];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let orthonormal = (0..3).all(|i| {
            (0..3).all(|j| (dot(col(i), col(j)) - if i == j { 1.0 } else { 0.0 }).abs() <= 1e-12)
        });
        let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
        if !orthonormal || (det - 1.0).abs() > 1e-12 || r.iter().flatten().any(|x| !x.is_finite()) {
            return Err(Error::invalid(format!(
                "rotation must be orthonormal with determinant +1: rotation is {r:?}"
            )));
        }
    }
    Ok(())
}

pub(crate) fn fractions(surface: &SurfaceIndex, blocks: &Blocks) -> Result<Vec<f64>> {
    check(blocks)?;
    let Some(r) = blocks.rotation else {
        return Ok(axis_fractions(
            &surface.mesh,
            &surface.bvh,
            &blocks.centroids,
            blocks.size,
        ));
    };
    let Some(&c0) = blocks.centroids.first() else {
        return Ok(vec![]);
    };
    let local = |v: &[f64; 3]| {
        let d = [v[0] - c0[0], v[1] - c0[1], v[2] - c0[2]];
        [0, 1, 2].map(|k| r[0][k] * d[0] + r[1][k] * d[1] + r[2][k] * d[2])
    };
    let mesh = surface
        .mesh
        .clone()
        .with_vertices(surface.mesh.vertices().iter().map(local).collect())?;
    let bvh = Bvh::new(&mesh.face_aabbs());
    let centroids: Vec<[f64; 3]> = blocks.centroids.iter().map(local).collect();
    Ok(axis_fractions(&mesh, &bvh, &centroids, blocks.size))
}

fn axis_fractions(mesh: &Mesh, bvh: &Bvh, centroids: &[[f64; 3]], size: [f64; 3]) -> Vec<f64> {
    let volume = size[0] * size[1] * size[2];
    centroids
        .par_iter()
        .map(|c| {
            let b = Aabb {
                min: [0, 1, 2].map(|k| c[k] - 0.5 * size[k]),
                max: [0, 1, 2].map(|k| c[k] + 0.5 * size[k]),
            };
            let column = Aabb {
                min: b.min,
                max: [b.max[0], b.max[1], f64::INFINITY],
            };
            let mut faces = vec![];
            bvh.query(&column, |f| faces.push(f));
            (box_volume(mesh, &faces, &b) / volume).clamp(0.0, 1.0)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SolidIndex;
    use drusa_core::fixtures;

    #[test]
    fn unit_cube_fractions() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let blocks = Blocks {
            centroids: vec![
                [0.5, 0.5, 0.5],
                [0.5, 0.5, 1.0],
                [1.0, 1.0, 1.0],
                [3.0, 0.5, 0.5],
                [0.5, 0.5, 0.25],
                [0.5, 0.5, -0.5],
            ],
            size: [1.0, 1.0, 1.0],
            rotation: None,
        };
        let f = s.block_fraction(&blocks).unwrap();
        assert_eq!(f, [1.0, 0.5, 0.125, 0.0, 0.75, 0.0]);
        let tiny = Blocks {
            centroids: vec![[0.5; 3]],
            size: [0.1; 3],
            rotation: None,
        };
        assert_eq!(s.block_fraction(&tiny).unwrap(), [1.0]);
    }

    #[test]
    fn rotation_about_z_by_45_degrees() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let h = 0.5f64.sqrt();
        let blocks = Blocks {
            centroids: vec![[0.5, 0.5, 0.5]],
            size: [1.0, 1.0, 1.0],
            rotation: Some([[h, -h, 0.0], [h, h, 0.0], [0.0, 0.0, 1.0]]),
        };
        let f = s.block_fraction(&blocks).unwrap()[0];
        let want = 2.0 * 2f64.sqrt() - 2.0;
        assert!((f - want).abs() < 1e-12, "{f} vs {want}");
    }

    #[test]
    fn bad_blocks_are_refused() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let mut b = Blocks {
            centroids: vec![[0.0; 3]],
            size: [1.0, 0.0, 1.0],
            rotation: None,
        };
        assert!(s.block_fraction(&b).is_err());
        b.size = [1.0; 3];
        b.rotation = Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]);
        assert!(s.block_fraction(&b).is_err());
        b.rotation = Some([[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(s.block_fraction(&b).is_err());
        b.rotation = None;
        b.centroids = vec![[f64::NAN; 3]];
        assert!(s.block_fraction(&b).is_err());
    }
}
