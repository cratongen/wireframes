use drusa_core::{Aabb, Mesh, Result};

/// Volume of the solid inside an axis-aligned box, by clamped-column integration.
pub fn box_volume(mesh: &Mesh, faces: &[u32], b: &Aabb) -> f64 {
    todo!()
}

/// Solid volume between consecutive `edges` (strictly increasing, at least 2).
pub fn volume_by_elevation(mesh: &Mesh, edges: &[f64]) -> Result<Vec<f64>> {
    todo!()
}
