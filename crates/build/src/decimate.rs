use drusa_core::{Mesh, Result};

#[derive(Clone, Debug, Default)]
pub struct DecimateOptions {
    pub target_faces: Option<usize>,
    pub max_error: Option<f64>,
    pub preserve_boundary: bool,
}

pub fn decimate(mesh: &Mesh, opts: &DecimateOptions) -> Result<Mesh> {
    todo!()
}
