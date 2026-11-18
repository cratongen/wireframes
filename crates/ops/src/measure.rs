use drusa_core::{Mesh, Result};

pub fn area(mesh: &Mesh) -> f64 {
    todo!()
}

pub fn face_areas(mesh: &Mesh) -> Vec<f64> {
    todo!()
}

/// Unit normals; `[0, 0, 0]` for degenerate faces.
pub fn face_normals(mesh: &Mesh) -> Vec<[f64; 3]> {
    todo!()
}

pub fn face_centroids(mesh: &Mesh) -> Vec<[f64; 3]> {
    todo!()
}

/// Angle-weighted unit normals; `[0, 0, 0]` for isolated vertices.
pub fn vertex_normals(mesh: &Mesh) -> Vec<[f64; 3]> {
    todo!()
}

/// Per face `(dip, dip direction)` in degrees: dip 0 to 90, direction 0 to 360
/// clockwise from `+y`; NaN direction for horizontal faces.
pub fn dip(mesh: &Mesh) -> Vec<(f64, f64)> {
    todo!()
}

/// Signed volume; `InvalidMesh` unless closed and consistently wound.
pub fn volume(mesh: &Mesh) -> Result<f64> {
    todo!()
}

/// Volume centroid, same precondition as `volume`; `InvalidInput` when the volume is 0.
pub fn centroid(mesh: &Mesh) -> Result<[f64; 3]> {
    todo!()
}

/// `None` when the area is 0.
pub fn area_centroid(mesh: &Mesh) -> Option<[f64; 3]> {
    todo!()
}
