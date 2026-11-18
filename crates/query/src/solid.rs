use drusa_core::{Mesh, Result};

use crate::{Blocks, Interval, Trace};

/// Sign agrees with signed distance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i8)]
pub enum Location {
    Inside = -1,
    Boundary = 0,
    Outside = 1,
}

pub struct SolidIndex {}

impl SolidIndex {
    /// `InvalidMesh { context: "solid" }` unless `validate(0.0)` says it is a solid.
    pub fn new(solid: &Mesh) -> Result<SolidIndex> {
        todo!()
    }

    pub fn mesh(&self) -> &Mesh {
        todo!()
    }

    pub fn locate(&self, p: [f64; 3]) -> Location {
        todo!()
    }

    pub fn locate_all(&self, points: &[[f64; 3]]) -> Vec<Location> {
        todo!()
    }

    pub fn contains_all(&self, points: &[[f64; 3]], boundary_inside: bool) -> Vec<bool> {
        todo!()
    }

    /// Negative inside, 0 on the boundary.
    pub fn signed_distance(&self, points: &[[f64; 3]]) -> Vec<f64> {
        todo!()
    }

    pub fn block_fraction(&self, blocks: &Blocks) -> Result<Vec<f64>> {
        todo!()
    }

    pub fn trace_intervals(&self, trace: &Trace) -> Result<Vec<Interval>> {
        todo!()
    }
}
