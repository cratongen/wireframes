use drusa_core::{Mesh, Polylines, Result};

/// `(negative, positive)`; points on the plane go to positive.
pub fn split_by_plane(mesh: &Mesh, origin: [f64; 3], normal: [f64; 3]) -> Result<(Mesh, Mesh)> {
    todo!()
}

/// Section lines at each offset along `normal`; empty `offsets` means `[0.0]`.
pub fn slice_by_plane(
    mesh: &Mesh,
    origin: [f64; 3],
    normal: [f64; 3],
    offsets: &[f64],
) -> Result<Polylines> {
    todo!()
}
