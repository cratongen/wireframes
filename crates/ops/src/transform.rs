use drusa_core::{Mesh, Result};

/// Affine `m`, row-major, applied to column vectors; flips faces when `det < 0`.
pub fn transform(mesh: &Mesh, m: [[f64; 4]; 4]) -> Result<Mesh> {
    todo!()
}

pub fn translate(mesh: &Mesh, d: [f64; 3]) -> Mesh {
    todo!()
}

/// `InvalidInput` when the mesh declares no unit.
pub fn to_length_unit(mesh: &Mesh, unit: &str) -> Result<Mesh> {
    todo!()
}

pub fn flip(mesh: &Mesh) -> Mesh {
    todo!()
}
