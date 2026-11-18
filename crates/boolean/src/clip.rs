use drusa_core::{Mesh, Result};
use drusa_ops::RepairReport;

use crate::BooleanOptions;

pub struct ClipResult {
    pub below: Mesh,
    pub above: Mesh,
    pub repairs: [Option<RepairReport>; 2],
}

/// Splits a solid by an open surface; `above` is on the side the surface normals point to.
pub fn clip_by_surface(solid: &Mesh, surface: &Mesh, opts: &BooleanOptions) -> Result<ClipResult> {
    todo!()
}
