use drusa_core::{Mesh, Result};

/// Edge-connected component label per face, numbered by first face.
pub fn component_labels(mesh: &Mesh) -> Vec<u32> {
    todo!()
}

/// One mesh per label, faces and vertices in input order.
pub fn components(mesh: &Mesh) -> Result<Vec<Mesh>> {
    todo!()
}
