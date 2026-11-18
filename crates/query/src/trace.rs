/// Desurveyed trace vertices and their measured depths, non-decreasing.
#[derive(Clone, Debug, PartialEq)]
pub struct Trace {
    pub xyz: Vec<[f64; 3]>,
    pub depth: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    pub from: f64,
    pub to: f64,
    pub inside: bool,
    pub entry_face: Option<u32>,
    pub exit_face: Option<u32>,
}
