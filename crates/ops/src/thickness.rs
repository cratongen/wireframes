use drusa_core::{Mesh, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThicknessMode {
    Vertical,
    Normal,
}

/// Highest surface z at each `(x, y)`; NaN where no face covers it.
pub fn surface_z(surface: &Mesh, xy: &[[f64; 2]]) -> Vec<f64> {
    todo!()
}

/// `z - surface_z` per point.
pub fn vertical_distance(surface: &Mesh, points: &[[f64; 3]]) -> Vec<f64> {
    todo!()
}

pub fn thickness(
    top: &Mesh,
    bottom: &Mesh,
    xy: &[[f64; 2]],
    mode: ThicknessMode,
) -> Result<Vec<f64>> {
    todo!()
}
