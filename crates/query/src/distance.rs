use drusa_core::{Mesh, Result};

use crate::Hit;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Closest {
    pub point: [f64; 3],
    pub face: u32,
    pub distance: f64,
}

pub struct SurfaceIndex {}

impl SurfaceIndex {
    /// Any mesh.
    pub fn new(mesh: &Mesh) -> SurfaceIndex {
        todo!()
    }

    pub fn closest(&self, points: &[[f64; 3]]) -> Vec<Closest> {
        todo!()
    }

    pub fn distance(&self, points: &[[f64; 3]]) -> Vec<f64> {
        todo!()
    }

    pub fn ray_hits(&self, origins: &[[f64; 3]], dirs: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>> {
        todo!()
    }

    pub fn segment_hits(&self, a: &[[f64; 3]], b: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>> {
        todo!()
    }
}
