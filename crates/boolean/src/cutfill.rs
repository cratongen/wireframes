use drusa_core::{Mesh, Polylines, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CutFill {
    pub cut: f64,
    pub fill: f64,
    pub net: f64,
    pub cut_area: f64,
    pub fill_area: f64,
    pub area: f64,
}

pub fn cut_fill(base: &Mesh, design: &Mesh, boundary: Option<&Polylines>) -> Result<CutFill> {
    todo!()
}
